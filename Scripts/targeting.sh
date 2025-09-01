#!/bin/bash

MIRNAs=$1
TARGETS=$2
JOBS=$3
JOBS=7

cat $MIRNAs | xargs -P $JOBS -n 2 bash -c 'targetfinder.pl -s "$0" -d "'"$TARGETS"'" -q "$1" -t 2 -p table 2>> /dev/null' | grep -v 'No results for'

