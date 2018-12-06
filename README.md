# Cytoscape Utils

Set of utilities and helpers for generating json networks in *cytoscape.js json* format from R dataframes in TSV format

## network_generator.R

This tool converts R dataframes in TSV of connected chromosome fragments and from optionally with a features dataframe file to *cytoscape json* format

### Dependencies
- [argparse](https://github.com/trevorld/argparse)
- [igraph](https://igraph.org/r/)
- [rjson](https://github.com/alexcb/rjson/)
- [tidyverse](https://www.tidyverse.org/)

### Usage
```
usage: ./network_generator.R [-h] [--wt_threshold WT_THRESHOLD]
                             [--features FEATURES] [--search SEARCH]
                             [--chromosome CHROMOSOME]
                             [--no-features-binarization] [--nearest]
                             [--expand EXPAND]
                             PCHiC file

Separated values file to cytoscape json mapper

positional arguments:
  PCHiC file            Separated values file PCHiC as input file

optional arguments:
  -h, --help            show this help message and exit
  --wt_threshold WT_THRESHOLD
                        The minimun value for considering the edge
  --features FEATURES   Separated values file of features as input file

  --search SEARCH       Search node by name or fragment position in the graph
                        to generate a neighborhood subgraph
  --chromosome CHROMOSOME
                        Filter by chromosome
  --no-features-binarization
                        Features will be binarized by default
  --nearest             Search the nearest range
  --expand EXPAND       Number of bases to expand the search by range

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

## Pipeline
### Dependencies
  - [Cytoscape](https://cytoscape.org/)
  - [GNU parallel](https://www.gnu.org/software/parallel)
  - [Wget](https://www.gnu.org/software/wget/)
  - [jq](https://stedolan.github.io/jq)

### Usage - generate JSONs
#### Nodejs headless mode (very slow)
`parallel --eta ./network_generator.R "PCHiC_interaction_map.txt --chromosome {} --features Features_mESC.txt | sed -e 's/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_enricher/layout_enricher | jq --compact-output . > chromosomes/chr{}.json" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y`
#### CyRest not headless mode (very fast)
`parallel --eta -j 1 ./network_generator.R "PCHiC_interaction_map.txt --chromosome {} --features Features_mESC.txt | sed -e 's/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_api_enricher | jq --compact-output .elements > chromosomes/chr{}.json" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y`
### Usage - generate images
`cytoscape -R 1234 &`

`parallel --eta -j 1 ./network_generator.R "PCHiC_interaction_map.txt --chromosome {} --features Features_mESC.txt | sed -e 's/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_api_enricher -f png -d chromosomes/chr{}" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y`

## backend.py
Backend for the [network_generator.R](network_generator.R) script

### Dependencies
  - [Flask](http://flask.pocoo.org/)

### Usage
FLASK_APP=backend.py flask run
