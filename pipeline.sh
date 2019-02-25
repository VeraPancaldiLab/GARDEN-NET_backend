#!/usr/bin/env bash

# Exit inmediately if any command has non-zero exit status
# Never use variables unset variables
# Exit inmediately if any command in a pipe chain has non-zero exit status
set -euo pipefail

# One input parameter is required
if [ "$1" == "" ]; then
  echo "An input folder is required: $0 input_datasets"
  exit 1
fi

# The input parameter has to be a folder
if [ ! -d "$1" ]; then
  echo 'A folder is required as input parameter'
  exit  1
fi

# Define the number of chromosomes for each organism
declare -A chromosomes_initial_number
chromosomes_initial_number['Homo_sapiens']=22
chromosomes_initial_number['Mus_musculus']=19
# Proccess each file
# Extract the realpath of the folder is needed to remove the last separator
for file in $(realpath "$1"/*); do
  case $file in
    # Ignore features files
    *.features) true;;
    *)
      # Remove from the end all characters until the last dot
      filename=$(basename "${file%%.*}")
      # Remove from the end all characters until the last dash
      organism=${filename%%-*}
      # Remove from the beggining all characters until the last dash
      cell_type=${filename##*-}
      # Default output_folder . if there is not a second input parameter
      output_folder=$(realpath ${2:-GARDEN-NET_DATA})
      mkdir -p "$output_folder/$organism/$cell_type/chromosomes"
      echo "$file:"
      # Generate chromosomes sequence
      chromosomes_seq_string="$(seq --separator ' ' 1 ${chromosomes_initial_number[$organism]}) X Y"
      # Size of the chromosomes sequence
      chromosomes_seq_size="$(wc -w <<< $chromosomes_seq_string)"
      printf "\t%s - %s\n" "$organism" "$cell_type"
      printf "\t%s chromosomes: %s\n" "$chromosomes_seq_size" "$chromosomes_seq_string"
      # Obtein features file using which has the same filename but with .features extension
      features="${file%.*}.features"
      features_parameter=""
      if [ -f "$features" ]; then
        printf "\tFeatures file: %s\n" "$features"
        features_parameter="--features $features"
      else
        printf "\tNo features file found"
      fi
      echo
      parallel --eta ./network_generator.R "--PCHiC $file $features_parameter --chromosome {} --pipeline $output_folder | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]\+\)\"/\1/' | ./layout_api_enricher | jq --monochrome-output --compact-output .elements > $output_folder/$organism/$cell_type/chromosomes/chr{}.json" ::: $chromosomes_seq_string
  esac
done
