# generate_inputs.py
## Overview
Python script for generating inputs.json
## Usage
  ```bash
  $ python generate_inputs.py [options]
  ```
  > Note: json, pandas, and getopt packages required.
**Options**
```
--complexity #Required. Determines whether lib_complexity task will be called
--r1=/path/to/read_1.fastq #Optional. For running individual samples
--r2=/path/to/read_2.fastq #Optional. For running individual samples
--bulk #Optional. For running all samples in a directory. -i required
-i /directory/containing/fastqs #Optional. For use with --bulk
-s /path/to/Samplemap.csv #Optional. For running all samples within a Samplemap.csv
-r /path/to/hg38_reference.tar.gz #Required.
-d docker-image (default: atex91/micro-c) #Required.
-j jobgroup #Required.
```
