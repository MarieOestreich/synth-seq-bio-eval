#!/bin/bash 

PATH_TO_REAL='path-to-real-data-csv'
PATH_TO_SYNTH='path-to-synthetic-data-csv'
EXPERIMENT_NAME='my-comparison'

python -u bio-eval.py -rd=$PATH_TO_REAL -sd=$PATH_TO_SYNTH -d=$EXPERIMENT_NAME