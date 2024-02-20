#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1

# Create an array to store line numbers to be filtered out
declare -a excluded_lines

# Process the input file and store line numbers to be excluded in the array
awk -F '' '{
    count_0 = 0
    count_1 = 0
    for (i = 1; i <= NF; i++) {
        if ($i == "0") {
            count_0++
        } else if ($i == "1") {
            count_1++
        }
    }
    total_fields = count_0 + count_1
    if (count_0 <= (total_fields * 0.05) || count_0 >= (total_fields * 0.95)) {
        excluded_lines[NR] = 1
    }
} END {
    # Print the line numbers to be excluded
    for (line_number in excluded_lines) {
        print line_number"d"
    }
}' "$input_file"

