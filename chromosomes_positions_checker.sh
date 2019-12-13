#!/usr/bin/env bash

if [ ! -d "$1" ]; then
  echo "Usage: $0 input_data_root_folder"
fi

# Let globbing to be recursive
shopt -s globstar
for chromosome_file in "$(realpath "$1")"/**/chromosomes/chr*.json; do
  case $chromosome_file in
  *)
    x=$(jq '.nodes[0].position.x' <"$chromosome_file")
    y=$(jq '.nodes[0].position.y' <"$chromosome_file")
    bad_positions=$(bc 2>/dev/null <<<"$x == $y && $x == 0")
    start=$(jq '.nodes[0].data.start' <"$chromosome_file")
    end=$(jq '.nodes[0].data.end' <"$chromosome_file")
    bad_chr_positions=$(bc 2>/dev/null <<<"$start == $end && $start == 1")
    if [ "$bad_positions" == "" ]; then
      echo -e "$chromosome_file: \e[31m\e[1mis empty!\e[0m"
    elif [ "$bad_positions" -eq 1 ]; then
      echo -e "$chromosome_file: \e[31m\e[1mnot generated positions!\e[0m"
    elif [[ "$bad_positions" == "" || "$bad_chr_positions" -eq 1 ]]; then
      echo -e "$chromosome_file: \e[31m\e[1mexponential number as chromosome positions!\e[0m"
    fi
    ;;
  esac
done
