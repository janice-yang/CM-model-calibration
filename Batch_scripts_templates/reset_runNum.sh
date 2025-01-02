#!/bin/bash

## Load Package
module load matlab/R2019a

## Reset runNum to 0 and save in GA folder
matlab -nodisplay -nosplash -r "clear; close all; runNum = -1; save ../GA/runNum runNum; quit"
