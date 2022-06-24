#!/usr/local/bin/python

import sys, os, getopt, json, pandas as pd

def arg_parse():
	r1 = None #read 1 filepath. For specifying individual samples
	r2 = None #read 2 filepath. For specifying individual samples
	ref = "/storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/all_sequences.tar.gz" #reference tarball
	samplemap = None #samplemap.csv path.
	inputs = None #directory containing fastq pairs. 
	bulk = False #Optional. Use --bulk to run all samples in a directory when no Samplemap.csv is present
	complexity = False #Optional. Use --complexity to run the library complexity estimation task.
	docker = "atex91/micro-c" #sets the virtual environment to run the majority of the workflow in.
	jobgroup = "/spencerlab/micro-c" #sets the jobgroup parameter for LSF backends.

	argv = sys.argv[1:]

	try:
		opts, args = getopt.getopt(argv, "s:i:r:d:j:", ["bulk", "r1=", "r2=", "complexity"])
	except getopt.GetoptError as err:
		print(err)
		opts = []

	for opt, arg in opts:
		if opt in ["--bulk"]:
			bulk = True
		if opt in ["--r1"]:
			r1 = arg
		if opt in ["--r2"]:
			r2 = arg
		if opt in ["--complexity"]:
			complexity = True
		if opt in ["-s"]:
			samplemap = arg
		if opt in ["-i"]:
			inputs = arg
		if opt in ["-r"]:
			ref = arg
		if opt in ["-d"]:
			docker = arg
		if opt in ["-j"]:
			jobgroup = arg

	return r1, r2, samplemap, inputs, ref, bulk, complexity, docker, jobgroup

def parse_samplemap(samplemap, inputs):
	#Parses Samplemap.csv information into the appropriate format for use with the WDL
	samplemap = pd.read_csv(samplemap)
	ones = smap[smap['File Name'].str.contains('R1')]['File Name']
	twos = smap[smap['File Name'].str.contains('R2')]['File Name']
	fastq_pairs = [{"read_1": f"{inputs}/{read_1}", "read_2": f"{inputs}/{read_2}"} for read_1, read_2 in zip(ones, twos)]
	sample_names = list(smap['Sample Name'].unique())
	return fastq_pairs, sample_names

def parse_directory(inputs):
	#Parses the data directory and generates a list of fastq pairs for use with the WDL. Also generates sample names for gathering output
	files = [file for file in os.listdir(inputs) if '.fastq' in file]
	ones = [file for file in files if '_R1_' in file]
	twos = [file for file in files if '_R2_' in file]
	fastq_pairs = [{"read_1": f"{inputs}/{read_1}", "read_2": f"{inputs}/{read_2}"} for read_1, read_2 in zip(sorted(ones), sorted(twos))]
	sample_names = [fname.split('_R1_')[0] for fname in sorted(ones)]
	return fastq_pairs, sample_names

def main():
	r1, r2, samplemap, inputs, ref, bulk, complexity, docker, jobgroup = arg_parse()
	
	if samplemap:
		fastq_pairs, sample_names = parse_samplemap(samplemap, inputs)
	elif bulk:
		fastq_pairs, sample_names = parse_directory(inputs)
	else:
		fastq_pairs = [{"read_1": r1, "read_2": r2}]

	input_json = {}
	input_json["micro_c.fastq"] = fastq_pairs
	input_json["micro_c.reference_index"] = ref
	input_json["micro_c.run_lib_complexity"] = complexity
	input_json["micro_c.docker"] = docker
	input_json["micro_c.jobGroup"] = jobgroup

	with open("inputs.json", "w") as outfile:
		json.dump(input_json, outfile)

if __name__ == "__main__":
	main()
