from subprocess import check_output, CalledProcessError
from flask import Flask, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

@app.route("/", defaults={'search': ''})
@app.route("/")
def main():
    expand   = request.args.get('expand')
    features = request.args.get('features')
    nearest  = request.args.get('nearest')
    search   = request.args.get('search')
    print(nearest)
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

    cmd_list = ["./network_generator.R ~/R_DATA/ChAs/PCHiC_interaction_map.txt"]

    if features is not None:
        cmd_list.append("--features")
        cmd_list.append("~/R_DATA/ChAs/Features_mESC.txt")

    if search:
        cmd_list.append('--search')
        cmd_list.append(search)

    if nearest is not None:
        cmd_list.append('--nearest')

    if expand:
        cmd_list.append('--expand')
        cmd_list.append(expand)


    print(" ".join(cmd_list))
    output = check_output(" ".join(cmd_list), shell=True, encoding='UTF-8')
    return output
