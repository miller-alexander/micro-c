# generate_inputs.py
## Overview
Python script for generating inputs.json
## Usage
  ```bash
  $ python generate_inputs.py [options]
  ```
  > Note: json, pandas, and getopt packages required.

## Options
```
--complexity                        #Optional. Default: False. Determines whether lib_complexity task will be called
--r1=/path/to/read_1.fastq          #Optional. For running individual samples
--r2=/path/to/read_2.fastq          #Optional. For running individual samples
--clean                             #Optional. Runs an intermediate cleanup step post-processing
--norun                             #Optional. If flagged, no commands will be executed. Input JSON will be printed to terminal.
-i /directory/containing/fastqs     #Required. Directory containing fastq files
-s /path/to/Samplemap.csv           #Optional. For running all samples within a Samplemap.csv
-r /path/to/hg38_reference.tar.gz   #Required.
-d docker-image                     #Required. Default: atex91/micro-c (recommended)
-j jobgroup                         #Required.
-w /path/to/wdl                     #Required.
```
