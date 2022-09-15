#!/usr/local/bin/python

import sys, os, getopt, json, pandas as pd

def arg_parse():
	r1 = None #read 1 filepath. For specifying individual samples
	r2 = None #read 2 filepath. For specifying individual samples
	sample_names = None
	ref = "/storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/all_sequences.tar.gz" #reference tarball
	samplemap = None #samplemap.csv path.
	inputs = None #directory containing fastq pairs. 
	wdl = "/storage1/fs1/dspencer/Active/spencerlab/alex/micro_c/pipeline/micro_c.wdl" #path to the .wdl
	outdir = None #directory to output files to
	complexity = False #Optional. Use --complexity to run the library complexity estimation task.
	docker = "atex91/micro-c" #sets the virtual environment to run the majority of the workflow in.
	jobgroup = "/spencerlab/micro-c-tasks" #sets the jobGroup parameter for LSF backends.
	config = "/storage1/fs1/dspencer/Active/spencerlab/alex/micro_c/pipeline/application.conf" #directory to LSF backend config file
	norun = False #Boolean. If True, no command will be run, but the input JSON will be printed
	clean = False #Boolean. If True, a final cleanup step will be run to remove intermediate files

	argv = sys.argv[1:]

	try:
		opts, args = getopt.getopt(argv, "w:o:s:i:r:n:d:j:c:", ["clean", "r1=", "r2=", "complexity", "norun"])
	except getopt.GetoptError as err:
		print(err)
		opts = []

	for opt, arg in opts:
		if opt in ["--clean"]:
			clean = True
		if opt in ["--r1"]:
			r1 = arg
		if opt in ["--r2"]:
			r2 = arg
		if opt in ["--complexity"]:
			complexity = True
		if opt in ["--norun"]:
			norun = True
		if opt in ["-w"]:
			wdl = arg
		if opt in ["-o"]:
			outdir = arg
		if opt in ["-s"]:
			samplemap = arg
		if opt in ["-i"]:
			inputs = arg
		if opt in ["-r"]:
			ref = arg
		if opt in ["-n"]:
			sample_names = [arg]
		if opt in ["-d"]:
			docker = arg
		if opt in ["-j"]:
			jobgroup = arg
		if opt in ["-c"]:
			config = arg

	return r1, r2, samplemap, wdl, outdir, inputs, ref, sample_names, complexity, docker, jobgroup, config, norun, clean

def parse_samplemap(samplemap, inputs):
	samplemap = pd.read_csv(samplemap)
	samples = []
	shortnames = list(set([name.split("_rep")[0] for name in samplemap['Sample Name'].unique()]))
	names = samplemap['Sample Name'].unique()
	for shortname in shortnames:
		df = samplemap[samplemap['Sample Name'].str.contains(shortname)]
		fastq_replicates = []
		for name in names:
			subdf = df[df['Sample Name'].str.contains(name)]
			if len(subdf) > 0:
				sample_name = subdf['Sample Name'].unique()[0]
				read_1 = list(subdf[subdf['File Name'].str.contains('R1')]['File Name'])
				read_2 = list(subdf[subdf['File Name'].str.contains('R2')]['File Name'])
				new_dict = {'name': sample_name, 'read_1': [f"{inputs}/{read}" for read in read_1], 'read_2': [f"{inputs}/{read}" for read in read_2]}
				fastq_replicates.append(new_dict)
		samples.append(fastq_replicates)
	return samples

def main():
	r1, r2, samplemap, wdl, outdir, inputs, ref, sample_names, complexity, docker, jobgroup, config, norun, clean = arg_parse()
	
	if samplemap:
		samples = parse_samplemap(samplemap, inputs)
	else:
		samples = [{"name": sample_names[0], "read_1": [r1], "read_2": [r2]}]


	for sample in samples:
		sample_id = sample[0]['name'].split("_rep")[0]
		input_json = {}
		input_json["micro_c.fastq"] = sample
		input_json["micro_c.reference_index"] = ref
		input_json["micro_c.run_lib_complexity"] = complexity
		input_json["micro_c.run_cleanup"] = clean
		input_json["micro_c.outdir"] = f"{outdir}/{sample_id}"
		input_json["micro_c.docker"] = docker
		input_json["micro_c.jobGroup"] = jobgroup
		
		if norun:
			print(json.dumps(input_json, indent=4, sort_keys=True))
		else:
			os.system(f"mkdir -p {outdir}/{sample_id}")

			with open(f"{outdir}/{sample_id}/inputs.json", "w") as outfile:
				json.dump(input_json, outfile)

			os.system(f"LSF_DOCKER_VOLUMES=\"/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer $HOME:$HOME\" bsub -oo {outdir}/{sample_id}/out.log -eo {outdir}/{sample_id}/err.log -q dspencer -G compute-oncology -g /spencerlab/micro-c -R \"span[hosts=1]\" -a \"docker(henrycwong/cromwell)\" /usr/bin/java -Dconfig.file={config} -jar /cromwell/cromwell.jar run -t wdl -m {outdir}/{sample_id}/outfile.json -i {outdir}/{sample_id}/inputs.json {wdl}")

if __name__ == "__main__":
	main()
