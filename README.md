# micro-c

## Overview
WDL implementation of [Dovetail micro-c](https://micro-c.readthedocs.io/en/latest/index.html) analysis pipeline.

## Installation
1. Git clone the WDL 
  ```bash
  $ git clone https://github.com/miller-alexander/micro-c
  ```
2. Install [Cromwell](https://github.com/broadinstitute/cromwell)
> Follow this [guide](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/)
3. [Configure](https://cromwell.readthedocs.io/en/stable/Configuring/) Cromwell for use with your compute framework
> **Note:** Configuration requires runtime attributes Int cpu, Int memory_gb, String jobGroup, and String docker.

## Usage
1. Generate the input json
> See the launcher.py readme for information on how to generate inputs.
2. Launch the WDL with Cromwell
  ```bash
  $ java -Dconfig.file=/path/to/backend_configuration -jar /cromwell/cromwell.jar run -t wdl -m outfile.json -i inputs.json micro-c.wdl
  ```
