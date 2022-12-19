#!/bin/bash 

PATH_TO_REAL=
PATH_TO_SYNTH=

python -u bio-eval.py -rd=$PATH_TO_REAL -sd=$PATH_TO_SYNTH > report.txt