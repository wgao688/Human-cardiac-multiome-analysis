#!/bin/bash

# gzip the STAR output files, since this is the format the SoupX and scanpy expect
nohup gzip ../raw_data/*/STAR/Solo.out/Gene/filtered/* &
nohup gzip ../raw_data/*/STAR/Solo.out/GeneFull/filtered/* &
nohup gzip ../raw_data/*/STAR/Solo.out/Gene/raw/* &
nohup gzip ../raw_data/*/STAR/Solo.out/GeneFull/raw/* &
