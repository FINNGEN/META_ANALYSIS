#!/bin/bash

# Take cromwell output files from the meta-analysis workflow and copy them to bucket with appropriate directory structure

# Usage: bash copy_cromwell_outputs.sh <cromwell_meta_output_prefix> <bucket>
# Example: bash copy_cromwell_outputs.sh /cromwell/output/path/6223ac01-a676-488c-8321-5d1588f9315c.meta_analysis. gs://my-bucket

# Check for required input
if [ $# -ne 2 ]; then
    echo "Usage: bash copy_cromwell_outputs.sh <cromwell_meta_output_prefix> <bucket>"
    exit 1
fi

# Get cromwell outputs and bucket
cromwell_outputs=$1
bucket=$2

# Check for trailing slashes
bucket=$(echo $bucket | sed 's/\/$//')

# Input files
metas=$cromwell_outputs"metas_with_rsids"
metas_filt=$cromwell_outputs"filtered_metas"
pdfs=$cromwell_outputs"pdfs"
pngs=$cromwell_outputs"pngs"
qc=$cromwell_outputs"qc"
qc_hits=$cromwell_outputs"qc_hits"
qc_xlsx=$cromwell_outputs"qc_xlsx"
gathered_qc=$cromwell_outputs"gathered_qc"

if [ -f $metas ]; then
    sed 's/$/*/' $metas | gcloud storage cp -I $bucket"/summary_stats/"
else
    echo "No metas_with_rsids file found at '$metas'"
fi

if [ -f $metas_filt ]; then
    sed 's/$/*/' $metas_filt | gcloud storage cp -I $bucket"/summary_stats_filtered/"
else
    echo "No filtered_metas file found at '$metas_filt'"
fi

if [ -f $pdfs ]; then
    cat $pdfs | gcloud storage cp -I $bucket"/plots/"
else
    echo "No pdfs file found at '$pdfs'"
fi

if [ -f $pngs ]; then
    cat $pngs | gcloud storage cp -I $bucket"/plots/"
else
    echo "No pngs file found at '$pngs'"
fi

if [ -f $qc ]; then
    cat $qc | gcloud storage cp -I $bucket"/qc/"
else
    echo "No qc file found at '$qc'"
fi

if [ -f $qc_hits ]; then
    cat $qc_hits | gcloud storage cp -I $bucket"/qc/"
else
    echo "No qc_hits file found at '$qc_hits'"
fi

if [ -f $qc_xlsx ]; then
    cat $qc_xlsx | gcloud storage cp -I $bucket"/qc/"
else
    echo "No qc_xlsx file found at '$qc_xlsx'"
fi

if [ -f $gathered_qc ]; then
    cat $gathered_qc | gcloud storage cp -I $bucket"/qc/"
else
    echo "No gathered_qc file found at '$gathered_qc'"
fi
