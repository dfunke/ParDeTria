#!/usr/bin/env bash

filename="$1"
runNumber=$2

while read -r run
do
    bin/Benchmark --run-string "$run" --run-number $runNumber --gen-only
done < "$filename"