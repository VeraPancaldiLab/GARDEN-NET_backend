#!/usr/bin/env bash

# Exit inmediately if any command has non-zero exit status
# Never use variables unset variables
# Exit inmediately if any command in a pipe chain has non-zero exit status
set -euo pipefail

input=""
output=""
only_metadata=false
only_pp_interactions=false

usage() { echo "Usage: $0 [-m] [-p] -i input_folder [-o output_folder]" 1>&2; exit 1; }

while getopts ":i:o:mp" opts; do
  case "$opts" in
    i) input=${OPTARG};;
    o) output=${OPTARG};;
    m) only_metadata=true;; # only generate metadata
    p) only_pp_interactions=true;; # only generate P-P networks
    *) usage;;
  esac
done
shift $((OPTIND-1))

# One input parameter is required
if [ -z "$input" ] || [ ! -d "$input" ]; then
  usage
  exit 1
fi

# Define the number of chromosomes for each organism
declare -A chromosomes_initial_number
chromosomes_initial_number['Homo_sapiens']=22
chromosomes_initial_number['Mus_musculus']=19
# Proccess each file
# Extract the realpath of the folder is needed to remove the last separator
for file in $(realpath "$input"/*); do
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
      output_folder=$(realpath "${output:-GARDEN-NET_DATA}")
      mkdir -p "$output_folder/$organism/$cell_type/"{chromosomes,metadata}
      echo "${file##*/}:"
      # Generate chromosomes sequence
      chromosomes_seq_string="$(seq --separator ' ' 1 ${chromosomes_initial_number[$organism]}) X Y PP"
      # Size of the chromosomes sequence
      chromosomes_seq_size=$(wc -w <<< "$chromosomes_seq_string")
      printf "\t%s - %s\n" "$organism" "$cell_type"
      printf "\t%s chromosomes: %s\n" "$chromosomes_seq_size" "$chromosomes_seq_string"
      # Obtein features file using which has the same filename but with .features extension
      features="${file%.*}.features"
      features_parameter=""
      if [ -f "$features" ]; then
        printf "\tFeatures file: %s\n" "${features##*/}"
        features_parameter="--features $features"
      else
        printf "\tNo features file found"
      fi
      echo

      only_pp_interactions_parameter=""
      if $only_pp_interactions; then
        only_pp_interactions_parameter="--only_pp_interactions"
      fi

      if ! $only_metadata; then
        # $chromosomes_seq_string has to be really splited by spaces in words so disable linter here
        # shellcheck disable=SC2086
        parallel --eta ./network_generator.R "--PCHiC $file $features_parameter --chromosome {} --pipeline $output_folder $only_pp_interactions_parameter | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]*\.\?[[:digit:]]\+\)\"/\1/' | ./layout_api_enricher | jq --monochrome-output --compact-output .elements > $output_folder/$organism/$cell_type/chromosomes/chr{}.json" ::: $chromosomes_seq_string
        # Always verify at the end all chromosomes are well generated
        echo "Verifying if all chromosomes are well generated..."
        ./chromosomes_positions_checker.sh "$output_folder"
      else
        rmdir "$output_folder/$organism/$cell_type/chromosomes" 2> /dev/null
        # $chromosomes_seq_string has to be really splited by spaces in words so disable linter here
        # shellcheck disable=SC2086
        parallel --eta ./network_generator.R "--PCHiC $file $features_parameter --chromosome {} --pipeline $output_folder $only_pp_interactions_parameter | sed -e '/chr/! s/\"[[:space:]]*\([[:digit:]]*\.\?[[:digit:]]\+\)\"/\1/' > /dev/null" ::: $chromosomes_seq_string
      fi
  esac
done

