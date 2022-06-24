#!/usr/local/bin/python

import sys, os, getopt, json, pandas as pd

def arg_parse():
	r1 = None #read 1 filepath. Only use if not using --samplemap
	r2 = None #read 2 filepath. Only use if not using --samplemap
	sample_names = None
	ref = "/storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/all_sequences.tar.gz" #reference tarball
	samplemap = False #Optional. Use --samplemap to automate launch for several fastqs from a samplemap.csv. Requires -m
	map_path = None #samplemap.csv path. Required if using --samplemap
	data_dir = None #directory containing fastq pairs. 
	wdl_path = None #directory containing the .wdl, inputs.json, lsf backend configuration
	outdir = None #directory to output files to
	clean = False #Optional. Use --clean to delete cromwell executions folder for each sample ran
	bulk = False #Optional. Use -b to run all samples in a directory when no Samplemap.csv is present
	complexity = False #Optional. Use --complexity to run the library complexity estimation task.

	argv = sys.argv[1:]

	try:
		opts, args = getopt.getopt(argv, "w:o:m:d:r:n:b:", ["r1=", "r2=", "samplemap", "clean", "complexity"])
	except getopt.GetoptError as err:
		print(err)
		opts = []

	for opt, arg in opts:
		if opt in ["--r1"]:
			r1 = arg
		if opt in ["--r2"]:
			r2 = arg
		if opt in ["--samplemap"]:
			samplemap = True
		if opt in ["--clean"]:
			clean = True
		if opt in ["--complexity"]:
			complexity = True
		if opt in ["-w"]:
			wdl_path = arg
		if opt in ["-o"]:
			outdir = arg
		if opt in ["-m"]:
			map_path = arg
		if opt in ["-d"]:
			data_dir = arg
		if opt in ["-r"]:
			ref = arg
		if opt in ["-n"]:
			sample_names = [arg]
		if opt in ["-b"]:
			bulk = True

	return r1, r2, samplemap, map_path, wdl_path, outdir, clean, data_dir, ref, sample_names, bulk, complexity

def parse_samplemap(map_path, data_dir):
	#Parses Samplemap.csv information into the appropriate format for use with the WDL
	smap = pd.read_csv(map_path)
	ones = smap[smap['File Name'].str.contains('R1')]['File Name']
	twos = smap[smap['File Name'].str.contains('R2')]['File Name']
	fastq_pairs = [{"read_1": f"{data_dir}/{read_1}", "read_2": f"{data_dir}/{read_2}"} for read_1, read_2 in zip(ones, twos)]
	sample_names = list(smap['Sample Name'].unique())
	return fastq_pairs, sample_names

def parse_directory(data_dir):
	#Parses the data directory and generates a list of fastq pairs for use with the WDL. Also generates sample names for gathering output
	files = [file for file in os.listdir(data_dir) if '.fastq' in file]
	ones = [file for file in files if '_R1_' in file]
	twos = [file for file in files if '_R2_' in file]
	fastq_pairs = [{"read_1": f"{data_dir}/{read_1}", "read_2": f"{data_dir}/{read_2}"} for read_1, read_2 in zip(sorted(ones), sorted(twos))]
	sample_names = [fname.split('_R1_')[0] for fname in sorted(ones)]
	return fastq_pairs, sample_names

def get_output(outdir, sample_names):
	#Retrieves the output specified by the wdl workflow.
	with open(f"{outdir}/outfile.json", "r") as infile:
		my_json = json.load(infile)
	for key in my_json['outputs']:
		if type(my_json['outputs'][key]) is list:
			for path, sample in zip(my_json['outputs'][key], sample_names):
				if os.path.exists(f"{outdir}/{sample}"):
					os.system(f"mv {path} {outdir}/{sample}")
				else:
					os.mkdir(f"{outdir}/{sample}")
					os.system(f"mv {path} {outdir}/{sample}")
		else:
			for sample in sample_names:
				if os.path.exists(f"{outdir}/{sample}"):
					os.system(f"cp {my_json['outputs'][key]} {outdir}/{sample}")
				else:
					os.mkdir(f"{outdir}/{sample}")
					os.system(f"cp {my_json['outputs'][key]} {outdir}/{sample}")
			os.system(f"rm {my_json['outputs'][key]}")

def cleanup(outdir):
	#Removes the cromwell-executions directory for the sample.
	with open(f"{outdir}/outfile.json", "r") as infile:
		my_json = json.load(infile)
	executions = my_json['workflowRoot']
	os.system(f"rm -r {executions}")

def main():
	r1, r2, samplemap, map_path, wdl_path, outdir, clean, data_dir, ref, sample_names, bulk, complexity = arg_parse()

	with open(f"{wdl_path}/inputs.json") as infile:
		input_json = json.load(infile)
	
	if samplemap:
		fastq_pairs, sample_names = parse_samplemap(map_path, data_dir)
	elif bulk:
		fastq_pairs, sample_names = parse_directory(data_dir)
	else:
		fastq_pairs = [{"read_1": r1, "read_2": r2}]

	input_json["micro_c.fastq"] = fastq_pairs
	input_json["micro_c.reference_index"] = ref
	input_json["micro_c.run_lib_complexity"] = complexity

	with open(f"{wdl_path}/inputs.json", "w") as outfile:
		json.dump(input_json, outfile)

	os.system(f"LSF_DOCKER_VOLUMES=\"/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer $HOME:$HOME\" bsub -oo {outdir}/out.log -eo {outdir}/err.log -K -q dspencer -G compute-oncology -g /spencerlab/micro-c -R \"span[hosts=1]\" -a \"docker(henrycwong/cromwell)\" /usr/bin/java -Dconfig.file={wdl_path}/application.conf -jar /cromwell/cromwell.jar run -t wdl -m {outdir}/outfile.json -i {wdl_path}/inputs.json {wdl_path}/micro_c.wdl")

	get_output(outdir, sample_names)

	if clean:
		cleanup(outdir)

if __name__ == "__main__":
	main()