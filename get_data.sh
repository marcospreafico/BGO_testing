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

echo
echo "Fetching data for batch $BATCH"
echo "Run list:"
echo "$RUNS"
echo 

ALL_GOOD=true

for RUN in $RUNS; do
    CALIB_FILE="./calib/phe_${RUN}.dat"

    if [ ! -f "$CALIB_FILE" ]; then
        printf "ERROR: calibration file %s not found.\nPlease run the calibration for run %s first.\n\n" \
            "$CALIB_FILE" "$RUN"
        ALL_GOOD=false
    fi

    DATA_FILE="./data/${RUN}.root"

    if [ ! -f "$DATA_FILE" ]; then
        echo "ERROR: data file $DATA_FILE missing."
        echo "Please transfer it from the DAQ machine using:"
        echo "scp streamdaq@193.206.147.141:/home/streamdaq/Documents/BGO_test/DAQ/${RUN}/RAW/SDataR_${RUN}.root ./data/${RUN}.root"
        echo
        ALL_GOOD=false
    fi
done

if [ "$ALL_GOOD" = true ]; then
    echo "All good! All calibration and data files are available."
    echo "You can run the data analysis now!"
else
    echo "Some files are missing. Please fix the errors above before running the analysis."
fi
echo 
