from subprocess import check_output, CalledProcessError
from flask import Flask, request, abort
from flask_cors import CORS
import shelve
import re

app = Flask(__name__)
CORS(app)

sanitize_pattern = re.compile('\W')

@app.route("/", defaults={'search': ''})
@app.route("/")
def main():
    expand   = request.args.get('expand')
    features = request.args.get('features')
    nearest  = request.args.get('nearest')
    search   = request.args.get('search')

    # Open or create a simple cache
    shelve_cache = shelve.open('.shelve_cache')

    # Valid URLs:
    #   '127.0.0.1:5000/'
    #   '127.0.0.1:5000/?features'
    #   '127.0.0.1:5000/?search=Y_581553'
    #   '127.0.0.1:5000/?search=Y_581553&features'
    #   '127.0.0.1:5000/?search=Hoxa1'
    #   '127.0.0.1:5000/?search=Hoxa1&features'
    #   '127.0.0.1:5000/?search=6:52155590-52158317'
    #   '127.0.0.1:5000/?search=6:52155590-52158317&nearest'
    #   '127.0.0.1:5000/?search=6:52155590-52158317&expand=20000'

    # Generate the keys for the cache
    nearest_key  = nearest if nearest is not None else ''
    features_key = features if features is not None else ''
    expand_key   = expand if expand is not None else ''
    key = '|'.join([search, features_key, nearest_key, expand_key])

    if  key not in shelve_cache:
        cmd_list = ["./network_generator.R ~/R_DATA/ChAs/PCHiC_interaction_map.txt"]

        if features is not None:
            cmd_list.append("--features")
            cmd_list.append("~/R_DATA/ChAs/Features_mESC.txt")

        if search:
            sanitized_search = sanitize_pattern.sub('', search)
            cmd_list.append('--search')
            cmd_list.append("'" + sanitized_search + "'")

        if nearest is not None:
            cmd_list.append('--nearest')

        if expand:
            cmd_list.append('--expand')
            cmd_list.append(expand)

        other_cmd = []
        other_cmd.append(r"sed -e 's/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/'")
        other_cmd.append('./layout_enricher/layout_enricher')
        other_cmd.append('jq --compact-output .')

        all_cmds = " ".join(cmd_list) + " | " + " | ".join(other_cmd)

        print(all_cmds)
        try:
            output = check_output(all_cmds, shell=True, encoding='UTF-8')
            shelve_cache[key] = output
        except:
            return abort(404)
    else:
        output = shelve_cache[key]


    shelve_cache.close()

    if len(output) == 3:
        return abort(404)

    return output
