# [GARDEN-NET](https://github.com/VeraPancaldiLab/GARDEN-NET) Utils

Set of utilities and helpers for generating json networks in *cytoscape.js json* format from R dataframes to TSV format in [GARDEN-NET](https://github.com/VeraPancaldiLab/GARDEN-NET)

## network_generator.R

This tool converts R dataframes in TSV of connected chromosome fragments and from optionally with a features dataframe file to *cytoscape json* format

### Dependencies
- [optparse](https://github.com/trevorld/r-optparse)
- [igraph](https://igraph.org/r/)
- [rjson](https://github.com/alexcb/rjson/)
- [tidyverse](https://www.tidyverse.org/)

### Usage
```
Usage: ./network_generator.R [options]
Separated values file to cytoscape json mapper

Options:
        -h, --help
                Show this help message and exit

        --PCHiC=PCHIC
                Separated values file PCHiC as input file

        --wt_threshold=WT_THRESHOLD
                The minimun value for considering the edge [default: 5]

        --features=FEATURES
                Separated values file of features as input file

        --search=SEARCH
                Search node by name or fragment position in the graph to generate a neighborhood subgraph

        --chromosome=CHROMOSOME
                Filter by chromosome

        --no-features-binarization
                Features will be binarized by default

        --nearest
                Search the nearest range

        --expand=EXPAND
                Number of bases to expand the search by range
        --pipeline=FOLDER
                Run the pipeline mode:
                                It generates the folder structure datasets/Organism/Cell_type

                                It takes the PCHiC file Organism_Cell_type.PCHiC and the features file Organism_Cell_type.features if exists

                                The metadata files: datasets/Organism/Cell_type/{search.Rdata, suggestions.json, features.json}

                                The chromosomes: datasets/Organism/Cell_type/chromosomes/chrNN.json (according to the organisms Homo_sapiens and Mus_musculus)

        --organism=ORGANISM
                Select an organism, only for searcher_query.R

        --cell_type=CELL_TYPE
                Select a cell_type, only for searcher_query.R

        --only_pp_interactions
                Use only promoter-promoter interactions from the network
```
## layout_enricher

This tool enrichs R *json cytoscape* networks with the *Cose* positions layout using Cytoscape.js and NodeJS

### Dependencies
- [yarn](https://yarnpkg.com/en)
    - [commander](https://github.com/tj/commander.js)
    - [cytoscape.js](http://js.cytoscape.org/)

### Usage
`yarn install`
```
Usage: layout_enricher json_cytoscape_file
  or   layout_enricher < json_cytoscape_file
  or   | layout_enricher

Options:
  -h, --help  output usage information
```

## layout_api_enricher

This tool enrichs R *json cytoscape* networks with the *Prefuse Force Directed* positions layout using CyRest and the network is generated in json or an image for each feature

### Dependencies
- [Cytoscape](http://cytoscape.org/)

### Usage
```
usage: layout_api_enricher [-h] [-p PORT] [-f {json,png,pdf,svg,all_images}]
                           [-n NAME] [-u URL] [-av API_VERSION] [-d DIRECTORY]
                           [-l LAYOUT]
                           [cytoscape_json_file]

Separated values file to cytoscape json mapper

positional arguments:
  cytoscape_json_file   Cytoscape Json format

optional arguments:
  -h, --help            show this help message and exit
  -p PORT, --port PORT  Cytoscape REST port
  -f {json,png,pdf,svg,all_images}, --format {json,png,pdf,svg,all_images}
                        The default output is 'json'
  -n NAME, --name NAME  Filename
  -u URL, --url URL     Cytoscape REST url
  -av API_VERSION, --api_version API_VERSION
                        Cytoscape REST version
  -d DIRECTORY, --directory DIRECTORY
                        Default directory
  -l LAYOUT, --layout LAYOUT
                        Cytoscape layout
```

## Usage
### Dependencies
  - [Cytoscape](https://cytoscape.org/)
  - [GNU parallel](https://www.gnu.org/software/parallel)
  - [Wget](https://www.gnu.org/software/wget/)
  - [jq](https://stedolan.github.io/jq)

### Usage - generate JSONs
#### Nodejs headless mode (very slow)
`parallel --eta ./network_generator.R "--PCHiC PCHiC_interaction_map.txt --chromosome {} --features Features_mESC.txt | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_enricher/layout_enricher | jq --monochrome-output --compact-output . > chromosomes/chr{}.json" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT`
#### CyRest not headless mode (very fast)
`cytoscape -R 1234 &`

`parallel --eta ./network_generator.R "--PCHiC PCHiC_interaction_map.txt --chromosome {} --features Features_mESC.txt | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_api_enricher | jq --monochrome-output --compact-output .elements > chromosomes/chr{}.json" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y`

### Usage - Generate search
`./network_generator.R --PCHiC PCHiC_interaction_map.txt --features Features_mESC.txt --search 'Hoxa1' | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_enricher/layout_enricher | jq --monochrome-output --compact-output .`

### Usage - Search query
`./search_query.R --PCHiC PCHiC_interaction_map.txt --features Features_mESC.txt --search 'Hoxa1' | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_enricher/layout_enricher | jq --monochrome-output --compact-output .`

### Usage - generate images
`cytoscape -R 1234 &`

`parallel --eta ./network_generator.R "--PCHiC PCHiC_interaction_map.txt --chromosome {} --features Features_mESC.txt | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_api_enricher -f png -d chromosomes/chr{}" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y`

## Pipeline
`Usage: ./pipeline.sh [-m] [-p] -i input_folder [-o output_folder]`

(-m is the flag for generating only metadata)

All files in input_folder need to have the next format:
- Organism-Cell_type.tsv
- Organism-Cell_type.features

Only Mus_musculus and Homo_sapiens are implemented by know


## backend.<span/>py
Backend for the [network_generator.R](network_generator.R) script

### Dependencies
  - [Flask](http://flask.pocoo.org/)
  - [Flask-CORS](https://flask-cors.readthedocs.io/)
  - [gunicorn](https://github.com/benoitc/gunicorn/)
  - [celery](http://www.celeryproject.org/)

### Usage
gunicorn --bind 0.0.0.0:5000 wsgi:app
celery -A backend.celery worker -l info

## Docker deployment
### Build image

docker build -t GARDEN-NET_utils

### Run container
#### wsgi<span/>.py development
`docker run --rm --interactive --tty --user "$(id -u):$(id -g)" --volume "$(pwd):/cytoscape_utils" --workdir /GARDEN-NET_utils GARDEN-NET_utils sh -c "gunicorn --workers 9 --bind unix:backend.sock wsgi:app"`

#### wsgi<span/>.py deployment
docker-compose up -d

##### Note about service management policy
The responsible of serve and restart the container is docker thanks to the restart=always policy instead of using the classic systemd unit file

