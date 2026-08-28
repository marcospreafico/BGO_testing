#!/bin/bash

BATCH=$1
RUNLIST="run_list.dat"

if [ -z "$BATCH" ]; then
    echo "Usage: $0 <batch>"
    exit 1
fi

if [ ! -f "$RUNLIST" ]; then
    echo "Error: $RUNLIST not found"
    exit 1
fi

RUNS=$(awk -v batch="$BATCH" '
    $1 == batch {
        found = 1
        next
    }

    found && $1 ~ /^[0-9]+$/ && length($1) <= 2 {
        exit
    }

    found {
        print $1
    }
' "$RUNLIST")

if [ -z "$RUNS" ]; then
    echo "No runs found for batch $BATCH"
    exit 1
fi

echo "Batch $BATCH"
echo "Runs:"
echo "$RUNS"

for RUN in $RUNS; do
    CALIB_FILE="./calib/phe_${RUN}.dat"

    # Check that the calibration file exists
    if [ ! -f "$CALIB_FILE" ]; then
        echo "ERROR: calibration file $CALIB_FILE not found. Skipping run $RUN."
        continue
    fi
    echo "Running scripts for run $RUN"
    root -l -q -b "preprocess.C($RUN, $BATCH)"
    root -l -q -b "event_builder.C($RUN, $BATCH)"

    EVT_FILE="./out/${RUN}_evt.root"

    if [ -f "$EVT_FILE" ]; then
        ROOT_FILES+=("$EVT_FILE")
    else
        echo "WARNING: $EVT_FILE not found."
    fi

done

# Merge all produced ROOT files
if [ ${#ROOT_FILES[@]} -gt 0 ]; then
    echo "Merging files into batch_${BATCH}_evt.root"
    rm "./out/batch${BATCH}.root"
    /opt/homebrew/Cellar/root/6.36.04/bin/hadd -o "./out/batch${BATCH}.root" "${ROOT_FILES[@]}"
else
    echo "ERROR: no ROOT files found to merge."
    exit 1
fi

root -l -q -b "ana_fit.C($BATCH)"

open "pdf/batch${BATCH}_fit.pdf"

