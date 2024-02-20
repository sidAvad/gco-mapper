#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1

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
    if (count_0 >= (total_fields * 0.05) && count_0 <= (total_fields * 0.95)) {
        print
    }
}' "$input_file"

