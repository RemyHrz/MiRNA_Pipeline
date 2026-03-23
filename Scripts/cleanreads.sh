#!/bin/bash

INPUT_DIR=$1
OUTPUT_DIR=$2
QUALITY_THRESHOLD=$3

MIN_SIZE=$4
MAX_SIZE=$5

TRIM_FRONT=$6
TRIM_TAIL=$7

LOGS_DIR=$8

ADAPTER=$9

THREADS=${10}

for file in "$INPUT_DIR"/*.fastq*; do
    if [[ -f "$file" ]]; then
        filename=$(basename "$file")
        filename=$(echo "$filename" | cut -f 1 -d '.')
		outname="${filename}.${QUALITY_THRESHOLD}.test"
        OUTPUT_FILE="${OUTPUT_DIR}/${outname}.fastq"
        REPORT_HTML="${LOGS_DIR}/${outname}_fastp.html"
        REPORT_HTML2="${LOGS_DIR}/${outname}_fastp2.html"
        REPORT_JSON="/dev/null"

        fastp -i "$file" -o "$OUTPUT_FILE" \
              --average_qual $QUALITY_THRESHOLD \
              -a "$ADAPTER" -w "$THREADS" \
              --html "$REPORT_HTML" \
              --json "$REPORT_JSON"
    fi
done
