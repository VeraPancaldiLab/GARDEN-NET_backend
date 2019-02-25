#!/usr/bin/env bash

if [ ! -d "$1" ]; then
  echo "Usage: $0 input_data_root_folder"
fi

# Let globbing to be recursive
shopt -s globstar
for chromosome_file in "$1"/**/chr*.json; do
  case $chromosome_file in
    # Remove chromosomes list json files from the processing
    *chromosomes.json) true ;;
    *)
      x=$(jq '.nodes[0].position.x' < "$chromosome_file")
      y=$(jq '.nodes[0].position.y' < "$chromosome_file")
      bad_positions=$(bc <<< "$x == $y && $x == 0")
      if [ "$bad_positions" -eq 1 ]; then
        echo "There are not generated positions in chromosome file: $chromosome_file"
      fi
      ;;
  esac
done
