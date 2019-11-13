#!/usr/bin/env bash

# Exit inmediately if any command has non-zero exit status
# Never use variables unset variables
# Exit inmediately if any command in a pipe chain has non-zero exit status
set -euo pipefail

input=""
output=""
only_metadata=false
only_broken_chromosomes=false

usage() {
  echo "Usage: $0 [-m] [-c] -i input_folder [-o output_folder]" 1>&2
  exit 1
}

while getopts "i:o:mc" opts; do
  case "$opts" in
  i) input=${OPTARG} ;;
  o) output=${OPTARG} ;;
  m) only_metadata=true ;; # only generate metadata
  c) only_broken_chromosomes=true ;; # only generate broken chromosomes
  *) usage ;;
  esac
done
shift $((OPTIND - 1))

# One input parameter is required
if [ -z "$input" ] || [ ! -d "$input" ]; then
  usage
  exit 1
fi

# Define the number of chromosomes for each organism
declare -A chromosomes_initial_number
chromosomes_initial_number['Homo_sapiens']=22
chromosomes_initial_number['Mus_musculus']=19
if [[ $only_broken_chromosomes == true ]]; then
  # Only chromosmes 1 and 2 are broken
  chromosomes_initial_number['Homo_sapiens']=2
elif [[ $only_metadata == true ]]; then
  # Only metadata is generated when chromosome 1 is required
  chromosomes_initial_number['Homo_sapiens']=1
  chromosomes_initial_number['Mus_musculus']=1
fi

# Proccess each file
# Extract the realpath of the folder is needed to remove the last separator
for file in $(realpath "$input"/*); do
  case $file in
  # Ignore features files
  *.features) true ;;
  *)
    if [[ $only_broken_chromosomes == true && ! $file =~ FoeT|naCD4|nCD8|tB|tCD4 ]]; then
      # Only generate broken chromosomes for FoeT|naCD4|nCD8|tB|tCD4 cell types
      continue
    fi
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
    if [[ $only_broken_chromosomes == true ]]; then
      chromosomes_seq_string="$(seq --separator ' ' 1 ${chromosomes_initial_number[$organism]})"
    elif [[ $only_metadata == true ]]; then
      chromosomes_seq_string="1"
    fi

    # Size of the chromosomes sequence
    chromosomes_seq_size=$(wc -w <<<"$chromosomes_seq_string")
    printf "\t%s - %s\n" "$organism" "$cell_type"
    if [[ $only_metadata == true ]]; then
      printf "\tMetadata only\n"
    elif [[ $only_broken_chromosomes == true ]]; then
      printf "\tBroken chromosomes: %s\n" "$chromosomes_seq_string"
    else
      printf "\t%s chromosomes: %s\n" "$chromosomes_seq_size" "$chromosomes_seq_string"
    fi

    # Obtain features file using which has the same filename but with .features extension
    features="${file%.*}.features"
    features_parameter=""
    if [ -f "$features" ]; then
      printf "\tFeatures file: %s\n" "${features##*/}"
      features_parameter="--features $features"
    else
      printf "\tNo features file found"
    fi
    echo

    intronic_parameter=""
    bait_names_parameter=""
    if [ "$organism" = 'Homo_sapiens' ]; then
      intronic_parameter="--intronic_regions ./intronic_regions.tsv"
      bait_names_parameter="--bait_names ./HindIII_annotation_ens37.txt"
    fi

    if [[ $only_broken_chromosomes == true ]]; then
      rmdir "$output_folder/$organism/$cell_type/metadata" 2>/dev/null
      parallel --eta ./network_generator.R "--PCHiC $file $features_parameter --chromosome {} $intronic_parameter --alias ./alias_databases/${organism}.tsv ${bait_names_parameter} | sed -e '/chr/! s/\"[[:space:]]*\(-\?[[:digit:]]*\.\?[[:digit:]]\+\)\"/\1/' > $output_folder/$organism/$cell_type/chromosomes/chr{}.json" ::: $chromosomes_seq_string
    elif [[ $only_metadata == true ]]; then
      rmdir "$output_folder/$organism/$cell_type/chromosomes" 2>/dev/null
      # $chromosomes_seq_string has to be really splited by spaces in words so disable linter here
      # shellcheck disable=SC2086
      parallel --eta ./network_generator.R "--PCHiC $file $features_parameter --chromosome {} --pipeline $output_folder $intronic_parameter --alias ./alias_databases/${organism}.tsv ${bait_names_parameter} | sed -e '/chr/! s/\"[[:space:]]*\(-\?[[:digit:]]*\.\?[[:digit:]]\+\)\"/\1/' > /dev/null" ::: 1
    else
      # $chromosomes_seq_string has to be really splited by spaces in words so disable linter here
      # shellcheck disable=SC2086
      parallel --eta ./network_generator.R "--PCHiC $file $features_parameter --chromosome {} --pipeline $output_folder $intronic_parameter --alias ./alias_databases/${organism}.tsv ${bait_names_parameter} | sed -e '/chr/! s/\"[[:space:]]*\(-\?[[:digit:]]*\.\?[[:digit:]]\+\)\"/\1/' | ./layout_api_enricher | jq --monochrome-output --compact-output .elements > $output_folder/$organism/$cell_type/chromosomes/chr{}.json" ::: $chromosomes_seq_string
    fi
    ;;
  esac
done

# Always verify at the end all chromosomes are well generated
if [[ $only_metadata == false && $only_broken_chromosomes == false ]]; then
  echo "Verifying if all chromosomes are well generated..."
  ./chromosomes_positions_checker.sh "$output_folder"
fi
