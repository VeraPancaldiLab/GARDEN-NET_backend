import gzip
import json
import os
import re
import shelve
import shutil
import subprocess
import tempfile
import traceback

from celery import Celery
from celery import states
from celery.exceptions import Ignore
from flask import abort
from flask import Flask
from flask import jsonify
from flask import request
from flask import url_for
from flask_cors import CORS
from werkzeug.utils import secure_filename

UPLOAD_FOLDER = "/tmp/flask_uploads"
ALLOWED_EXTENSIONS = set(["txt", "txt.gz", "bed", "bed.gz"])

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})
# Add redis broker to Flask app
app.config["CELERY_BROKER_URL"] = "redis://redis:6379/0"
app.config["CELERY_RESULT_BACKEND"] = "redis://redis:6379/0"
# Initialize celery distributed task queue
celery = Celery(app.name, broker=app.config["CELERY_BROKER_URL"])
celery.conf.update(app.config)

# Upload folder setting
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER


# SANITIZE_PATTERN = re.compile('[^-a-zA-Z0-9:]')


@app.route("/", defaults={"search": ""})
@app.route("/")
def main():
    search = request.args.get("search")
    organism = request.args.get("organism")
    cell_type = request.args.get("cell_type")

    # Open or create a simple cache
    shelve_cache = shelve.open(".shelve_cache")

    # Valid URLs:
    #   '127.0.0.1:5000/'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&features'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&search=Y_581553'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&search=Y_581553&features'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&search=Hoxa1'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&search=Hoxa1&features'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&search=6:52155590-52158317'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&search=6:52155590-52158317&nearest'
    #   '127.0.0.1:5000/?organism=Mus_musculus&cell_type=Embryonic_stem_cells&search=6:52155590-52158317&expand=20000'

    # Generate the keys for the cache

    key = "|".join([search, organism, cell_type])

    if key not in shelve_cache:
        cmd_list = ["./search_query.R"]

        if search:
            # sanitized_search = SANITIZE_PATTERN.sub('', search.split()[0])
            sanitized_search = search.split()[0]
            cmd_list.append("--search=" + "'" + sanitized_search + "'")

        if organism:
            cmd_list.append("--organism=" + organism)

        if cell_type:
            cmd_list.append("--cell_type=" + cell_type)

        other_cmd = []
        other_cmd.append(
            r"sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]*\.\?[[:digit:]]\+\)\"/\1/'"
        )
        other_cmd.append("./layout_enricher/layout_enricher")
        other_cmd.append("jq --compact-output .")

        all_cmds = " ".join(cmd_list) + " | " + " | ".join(other_cmd)

        print(all_cmds)
        try:
            output = subprocess.check_output(all_cmds, shell=True, encoding="UTF-8")
            shelve_cache[key] = output
        except:
            # 408 Request Timeout
            return abort(408)
    else:
        output = shelve_cache[key]

    shelve_cache.close()

    if len(output) == 3:
        return abort(404)

    return output


@app.route("/upload_features", methods=["POST"])
def upload_features():
    print("upload_features")
    # http://flask.pocoo.org/docs/1.0/patterns/fileuploads/
    organism = request.args.get("organism")
    cell_type = request.args.get("cell_type")
    print(request.files["features"])
    features_file_object = request.files["features"]
    features_filename = secure_filename(features_file_object.filename)
    tmp_dir = tempfile.mkdtemp()
    features_path = os.path.join(tmp_dir, features_filename)
    features_file_object.save(features_path)
    # https://stackoverflow.com/a/28305785
    # uncompressed = gzip.decompress(features_file_object.read())
    cat_command = "cat"
    if features_filename.endswith(".gz"):
        cat_command = "zcat"
    headers_number = int(
        subprocess.check_output(
            " ".join(
                [
                    cat_command,
                    features_path,
                    "|",
                    "head -n1",
                    "|",
                    "sed 's/[^\t]//g'",
                    "|",
                    "awk '{print length + 1}'",
                ]
            ),
            shell=True,
        ).strip()
    )

    features_file_type = "unknown"
    if headers_number == 2:
        features_file_type = "features_on_nodes"
    if headers_number == 4:
        try:
            last_column = float(
                subprocess.check_output(
                    " ".join(
                        [
                            cat_command,
                            features_path,
                            "|",
                            "head -n1",
                            "|",
                            "awk '{print $NF}'",
                        ]
                    ),
                    shell=True,
                ).strip()
            )
            features_file_type = "bed3"
        except:
            features_file_type = "chromhmm"

    elif headers_number == 6:
        features_file_type = "bed6"
    elif headers_number == 9 or headers_number == 10:
        features_file_type = "macs2"

    # print("Header: Number of columns = " + str(headers_number))
    task = processing_features.apply_async(
        args=(tmp_dir, organism, cell_type, features_path, features_file_type)
    )
    return (
        jsonify({}),
        202,
        {
            "Access-Control-Expose-Headers": "Location",
            "Location": url_for("features_task", task_id=task.id),
        },
    )


@celery.task(bind=True)
def processing_features(
    self, tmp_dir, organism, cell_type, features_file, features_file_type
):

    if features_file_type == "unknown":
        self.update_state(
            state="FAILURE",
            # https://www.distributedpython.com/2018/09/28/celery-task-states/
            meta={
                "percentage": 1,
                "total": 1,
                "message": "Error unknown feature file format",
                "exc_message": "Error unknown feature file format",
                "exc_type": "str",
            },
        )
    fifo_file = os.path.join(tmp_dir, "fifo")
    os.mkfifo(fifo_file)

    R_process = subprocess.Popen(
        " ".join(
            [
                "Rscript merge_features.R",
                "--fifo_file",
                fifo_file,
                "--organism",
                organism,
                "--cell_type",
                cell_type,
                "--features_file",
                features_file,
                "--features_file_type",
                features_file_type,
            ]
        ),
        shell=True,
        encoding="UTF-8",
    )

    fifo_data = ""
    while fifo_data != "QUIT":
        with open(fifo_file, "r") as fifo:
            while True:
                fifo_data = fifo.read().strip()
                if not fifo_data or fifo_data == "QUIT":
                    break
                counter = re.search(r":\s*(\d+)", fifo_data)[1]

                r_progress_info = fifo_data.split(":")[0]
                self.update_state(
                    state="PROGRESS",
                    meta={
                        "percentage": counter,
                        "total": 100,
                        "message": r_progress_info + "...",
                    },
                )

    return_code = R_process.returncode

    if return_code != 0 and return_code is not None:
        self.update_state(
            state="FAILURE",
            # https://www.distributedpython.com/2018/09/28/celery-task-states/
            meta={
                "percentage": 1,
                "total": 1,
                "message": "Error proccesing the feature file",
                "exc_message": "Error proccesing the feature file",
                "exc_type": "str",
            },
        )
        # ignore the task so no other state is recorded
        raise Ignore()

    features = None
    features_metadata = None
    try:
        with open(os.path.join(tmp_dir, "features.json"), "r") as f:
            features = json.loads(f.read())
        with open(os.path.join(tmp_dir, "features_metadata.json"), "r") as f:
            features_metadata = json.loads(f.read())

        os.remove(fifo_file)
        shutil.rmtree(tmp_dir)
    except:
        self.update_state(
            state="FAILURE",
            # https://www.distributedpython.com/2018/09/28/celery-task-states/
            meta={
                "percentage": 1,
                "total": 1,
                "message": "Error proccesing the feature file",
                "exc_message": "Error proccesing the feature file",
                "exc_type": "str",
            },
        )
        # ignore the task so no other state is recorded
        raise Ignore()

    return {
        "percentage": 100,
        "total": 100,
        "message": "Features file processed successfully",
        "result": {"features": features, "features_metadata": features_metadata},
    }


@app.route("/status/<task_id>")
def features_task(task_id):
    task = processing_features.AsyncResult(task_id)
    if task.state == "PENDING":
        # job did not start yet
        response = {
            "state": task.state,
            "percentage": 0,
            "total": 1,
            "message": "Uploading features file...",
        }
    elif task.state != "FAILURE":
        response = {
            "state": task.state,
            "percentage": task.info.get("percentage", 0),
            "total": task.info.get("total", 1),
            "message": task.info.get("message", ""),
        }
        if "result" in task.info:
            response["result"] = task.info["result"]
    else:
        # something went wrong in the background job
        # celery.backends.base.str -> str
        message = (
            str(task.info)
            .replace("'", "")
            .replace(", ", "")
            .replace("(", "")
            .replace(")", "")
        )
        response = {
            "state": task.state,
            "percentage": 1,
            "total": 1,
            "message": message,  # this is the exception raised
        }
    return jsonify(response)


if __name__ == "__main__":
    app.run(host="0.0.0.0")
