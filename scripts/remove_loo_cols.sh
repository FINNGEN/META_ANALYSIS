#!/bin/bash

# Remove leave one out columns from summary stats
# Usage: bash remove_leave_one_out_cols.sh gs://path/to/summary/stats/file
filename=$(basename $1)
gsdir=$(dirname $1)
rmv_cols=$(gsutil cat $1 | zcat -f | head -1 | tr "\t" "\n" | grep -n "^leave_" | cut -f 1 -d ":" | tr "\n" "," | sed 's/,$//')
if [ -z "$rmv_cols" ]
then
    echo "No leave one out columns found in $1"
    exit 1
fi
gcloud storage cat $1 | zcat -f | cut --complement -f $(echo $rmv_cols) | bgzip > $filename
tabix -s 1 -b 2 -e 2 $filename
gcloud storage cp $filename $filename.tbi $gsdir
rm $filename $filename.tbi