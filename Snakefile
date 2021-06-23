## Configuration file
import os
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("".join(["Make sure there is a config.yaml file in ", os.getcwd(), 
			" or specify one with the --configfile commandline parameter."]))

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['annotation', 'organism', 'build', 'release', 'txome', 'genome', 'gtf', 'STARindex', 'readlength', 'fldMean', 'fldSD', 'metatxt', 'ncores', 'FASTQ', 'fqext1', 'fqext2', 'fqsuffix', 'output', 'useCondaR', 'Rbin', 'run_trimming', 'run_STAR']
for k in configvars:
	if k not in config:
		config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
	if str is None:
		str = ''
	return str

config['gtf'] = sanitizefile(config['gtf'])
config['genome'] = sanitizefile(config['genome'])
config['STARindex'] = sanitizefile(config['STARindex'])
config['metatxt'] = sanitizefile(config['metatxt'])
config['salmonindex'] = sanitizefile(config['salmonindex'])

## Read metadata
if not os.path.isfile(config["metatxt"]):
	sys.exit("".join(["Metadata file ", config["metatxt"], " does not exist."]))

import pandas as pd
samples = pd.read_csv(config["metatxt"], sep='\t')

if not set(['names','type']).issubset(samples.columns):
	sys.exit("".join(["Make sure 'names' and 'type' are columns in ", config["metatxt"]]))


## Sanitize provided input and output directories
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

outputdir = getpath(config["output"])
FASTQdir = getpath(config["FASTQ"])

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

## Define the R binary
Rbin = config["Rbin"]

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		os.path.join(outputdir, "MultiQC", "multiqc_report.html")

rule setup:
	input:
		os.path.join(outputdir, "Rout", "pkginstall_state.txt"),
		os.path.join(outputdir, "Rout", "softwareversions.done")

## Install R packages
rule pkginstall:
	input:
		script = "scripts/install_pkgs.R"
	output:
	  	os.path.join(outputdir, "Rout", "pkginstall_state.txt")
	params:
		flag = config["annotation"],
		ncores = config["ncores"],
		organism = config["organism"],
		Rbin = Rbin
	priority:
		50
	conda:
		Renv
	log:
		os.path.join(outputdir, "Rout", "install_pkgs.Rout")
	benchmark:
	  	os.path.join(outputdir, "benchmarks", "install_pkgs.txt")
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args outtxt='{output}' ncores='{params.ncores}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext1"]), "_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext2"]), "_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext1"]), "_val_1_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext2"]), "_val_2_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "{sample}_trimmed_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist())

## STAR alignment
rule runstar:
	input:
		expand(os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam.bai"), sample = samples.names.values.tolist()),
		expand(os.path.join(outputdir, "STARbigwig", "{sample}_Aligned.sortedByCoord.out.bw"), sample = samples.names.values.tolist())

rule run_removedup:
	input:
		expand(os.path.join(outputdir, "BAM_deduplicated", "{sample}_deduplicated.bam.bai"), sample = samples.names.values.tolist())

rule run_filter_second_read:
	input:
		expand(os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.bam.bai"), sample = samples.names.values.tolist())

rule run_clipper:
	input:
		expand(os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.ucsc.bam.bai"), sample = samples.names.values.tolist()),
		expand(os.path.join(outputdir, config["clipper"], "{sample}_deduplicated.r2.ucsc.clipper_peaks.bed"), sample = samples.names.values.tolist())

rule run_norm_bigwig:
	input:
		expand(os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.ucsc.replicatesNorm.bw"), sample = samples.names.values.tolist())

rule run_bigwig_second_read:
	input:
		expand(os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.bw"), sample = samples.names.values.tolist())

## Salmon quantification
rule runsalmonquant:
	input:
		expand(os.path.join(outputdir, "salmon", "{sample}", "quant.sf"), sample = samples.names.values.tolist()),
		os.path.join(outputdir, "outputR", "tximeta_se.rds")


## List all the packages that were used by the R analyses
rule listpackages:
	log:
		os.path.join(outputdir, "Rout", "list_packages.Rout")
	params:
		Routdir = os.path.join(outputdir, "Rout"),
		outtxt = os.path.join(outputdir, "R_package_versions.txt"),
		script = "scripts/list_packages.R",
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {params.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	output:
		touch(os.path.join(outputdir, "Rout", "softwareversions.done"))
	log:
		os.path.join(outputdir, "logs", "softversions.log")
	conda:
		"envs/environment.yaml"
	shell:
		"echo -n 'ARMOR version ' && cat version; "
		"salmon --version; trim_galore --version; "
		"echo -n 'cutadapt ' && cutadapt --version; "
		"fastqc --version; STAR --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##
## Generate STAR index
rule starindex:
	input:
		genome = config["genome"],
		gtf = config["gtf"]
	output:
		os.path.join(config["STARindex"], "SA"),
		os.path.join(config["STARindex"], "chrNameLength.txt")
	log:
		os.path.join(outputdir, "logs", "STAR_index.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_index.txt")
	params:
		STARindex = lambda wildcards, output: os.path.dirname(output[0]),   ## dirname of first output
		readlength = config["readlength"],
		starextraparams = config["additional_star_index"]
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.STARindex} "
		"--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength} "
		"{params.starextraparams}"

## Generate Salmon index from merged cDNA and ncRNA files
rule salmonindex:
	input:
		txome = config["txome"]
	output:
		os.path.join(config["salmonindex"], "versionInfo.json")
	log:
		os.path.join(outputdir, "logs", "salmon_index.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "salmon_index.txt")
	params:
		salmonoutdir = lambda wildcards, output: os.path.dirname(output[0]),   ## dirname of first output
		anno = config["annotation"],
		salmonextraparams = config["additional_salmon_index"]
	conda:
		"envs/environment.yaml"
	shell:
	  """
	  if [ {params.anno} == "Gencode" ]; then
      echo 'Salmon version:\n' > {log}; salmon --version >> {log};
  	  salmon index -t {input.txome} -i {params.salmonoutdir} --gencode {params.salmonextraparams}
    else
  	  echo 'Salmon version:\n' > {log}; salmon --version >> {log};
      salmon index -t {input.txome} -i {params.salmonoutdir} {params.salmonextraparams}
    fi
    """

## Generate linkedtxome mapping
rule linkedtxome:
	input:
		txome = config["txome"],
		gtf = config["gtf"],
		salmonidx = os.path.join(config["salmonindex"], "versionInfo.json"),
		script = "scripts/generate_linkedtxome.R",
		install = os.path.join(outputdir, "Rout", "pkginstall_state.txt")
	log:
		os.path.join(outputdir, "Rout", "generate_linkedtxome.Rout")
	benchmark:
		os.path.join(outputdir, "benchmarks", "generate_linkedtxome.txt")
	output:
		"".join([config["salmonindex"], ".json"])
	params:
		flag = config["annotation"],
		organism = config["organism"],
		release = str(config["release"]),
		build = config["build"],
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args transcriptfasta='{input.txome}' salmonidx='{input.salmonidx}' gtf='{input.gtf}' annotation='{params.flag}' organism='{params.organism}' release='{params.release}' build='{params.build}' output='{output}'" {input.script} {log}'''


## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip")
	params:
		FastQC = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "fastqc_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "fastqc_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## FastQC, trimmed reads
rule fastqctrimmed:
	input:
		fastq = os.path.join(outputdir, "FASTQtrimmed", "{sample}.fq.gz")
	output:
		os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip")
	params:
		FastQC = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "fastqc_trimmed_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "fastqc_trimmed_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"



# The config.yaml files determines which steps should be performed
def multiqc_input(wildcards):
	input = []
	input.extend(expand(os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist()))
	input.extend(expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext1"]), "_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
	input.extend(expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext2"]), "_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
	if config["run_trimming"]:
		input.extend(expand(os.path.join(outputdir, "FASTQtrimmed", "{sample}_trimmed.fq.gz"), sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
		# input.extend(expand(os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1_val_1.fq.gz"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
		# input.extend(expand(os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2_val_2.fq.gz"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FastQC", "{sample}_trimmed_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext1"]), "_val_1_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext2"]), "_val_2_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
		# input.extend(expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext1"]), "_val_1_val_1_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
		# input.extend(expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext2"]), "_val_2_val_2_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()))
	if config["run_STAR"]:
		input.extend(expand(os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam.bai"), sample = samples.names.values.tolist()))
	return input

## Determine the input directories for MultiQC depending on the config file
def multiqc_params(wildcards):
	param = [os.path.join(outputdir, "FastQC")]
	if config["run_trimming"]:
		param.append(os.path.join(outputdir, "FASTQtrimmed"))
	if config["run_STAR"]:
		param.append(os.path.join(outputdir, "STAR"))
	return param

## MultiQC
rule multiqc:
	input:
		multiqc_input
	output:
		os.path.join(outputdir, "MultiQC", "multiqc_report.html")
	params:
		inputdirs = multiqc_params,
		MultiQCdir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "multiqc.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "multiqc.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "FASTQtrimmed", "{sample}_trimmed.fq.gz")
	params:
		FASTQtrimmeddir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "trimgalore_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"

rule trimgalorePE:
	input:
		fastq1 = os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])),
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"]))
	params:
		FASTQtrimmeddir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "trimgalore_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2}"

## second round of adapter trimming to remove duplicate adapter ligation events
## we specify the TruSeq adapters taken from https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
# --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
rule trimgalorePE_2nd_round:
	input:
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"]))
	output:
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1_val_1.fq.gz"])),
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2_val_2.fq.gz"]))
	params:
		FASTQtrimmeddir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "trimgalore_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2} --illumina"


## ------------------------------------------------------------------------------------ ##
## Salmon abundance estimation
## ------------------------------------------------------------------------------------ ##
# Estimate abundances with Salmon
rule salmonSE:
	input:
		index = os.path.join(config["salmonindex"], "versionInfo.json"),
		fastq = os.path.join(outputdir, "FASTQtrimmed", "{sample}_trimmed.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "salmon", "{sample}", "quant.sf")
	log:
		os.path.join(outputdir, "logs", "salmon_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "salmon_{sample}.txt")
	threads:
		config["ncores"]
	params:
		salmonindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		salmondir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		salmonextraparams = config["additional_salmon_quant"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -r {input.fastq} "
		"-o {params.salmondir}/{wildcards.sample} -p {threads} {params.salmonextraparams}"

rule salmonPE:
	input:
		index = os.path.join(config["salmonindex"], "versionInfo.json"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "salmon", "{sample}", "quant.sf")
	log:
		os.path.join(outputdir, "logs", "salmon_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "salmon_{sample}.txt")
	threads:
		config["ncores"]
	params:
		salmonindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		salmondir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		salmonextraparams = config["additional_salmon_quant"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -1 {input.fastq1} -2 {input.fastq2} "
		"-o {params.salmondir}/{wildcards.sample} -p {threads} {params.salmonextraparams}"


## ------------------------------------------------------------------------------------ ##
## Transcript quantification
## ------------------------------------------------------------------------------------ ##
## tximeta
rule tximeta:
	input:
	    os.path.join(outputdir, "Rout", "pkginstall_state.txt"),
		expand(os.path.join(outputdir, "salmon", "{sample}", "quant.sf"), sample = samples.names.values.tolist()),
		metatxt = config["metatxt"],
		salmonidx = os.path.join(config["salmonindex"], "versionInfo.json"),
		json = "".join([config["salmonindex"], ".json"]),
		script = "scripts/run_tximeta.R"
	output:
		os.path.join(outputdir, "outputR", "tximeta_se.rds")
	log:
		os.path.join(outputdir, "Rout", "tximeta_se.Rout")
	benchmark:
		os.path.join(outputdir, "benchmarks", "tximeta_se.txt")
	params:
		salmondir = lambda wildcards, input: os.path.dirname(os.path.dirname(input[1])),   ## dirname of second output
		flag = config["annotation"],
		organism = config["organism"],
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args salmondir='{params.salmondir}' json='{input.json}' metafile='{input.metatxt}' outrds='{output}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''


## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
rule starSE:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq = os.path.join(outputdir, "FASTQtrimmed", "{sample}_trimmed.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs", "STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_{sample}.txt")
	params:
		STARindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		STARdir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		starextraparams = config["additional_star_align"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"

rule starPE:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs", "STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_{sample}.txt")
	params:
		STARindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		STARdir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		starextraparams = config["additional_star_align"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"

rule starPE_2nd_round:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1_val_1.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2_val_2.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "STAR_round2", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs", "STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_{sample}.txt")
	params:
		STARindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		STARdir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		starextraparams = config["additional_star_align_2nd_round"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"


## Index bam files
rule bamindex:
	input:
		bam = os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	output:
		os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam.bai")
	log:
		os.path.join(outputdir, "logs", "samtools_index_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "samtools_index_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam"),
		chrl = os.path.join(config["STARindex"], "chrNameLength.txt")
	output:
		os.path.join(outputdir, "STARbigwig", "{sample}_Aligned.sortedByCoord.out.bw")
	params:
		STARbigwigdir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "bigwig_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "bigwig_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | LC_COLLATE=C sort -k1,1 -k2,2n > "
		"{params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## Normalize the two replicates to the mean library size 
## the scaling factors were computed with R (clipper_analysis.Rmd)
def get_sf(wildcards):
    return config["sf"][wildcards.sample] 

rule normalize_replicates_bigwig:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated", "{sample}_deduplicated.r2.ucsc.bam"),
		chrl = "reference/hg38.chrom.sizes_GRCh38.txt"
	output:
		os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.ucsc.replicatesNorm.bw")
	params:
		bigwigdir = os.path.join(outputdir, "bigwig"),
		bedgraph = "{sample}_deduplicated.r2.ucsc.bedGraph",
		sf = get_sf
	log:
		"logs/bigwig_filtered_{sample}.log"
	benchmark:
		"benchmarks/bigwig_{sample}.txt"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -scale {params.sf} -ibam {input.bam} -bg > {params.bigwigdir}/{params.bedgraph}; "
		"LC_ALL='C' sort -k1,1 -k2,2n {params.bigwigdir}/{params.bedgraph} > {params.bigwigdir}/{params.bedgraph}_sorted; "
		"bedGraphToBigWig {params.bigwigdir}/{params.bedgraph}_sorted {input.chrl} {output}; "
		"rm -f {params.bigwigdir}/{params.bedgraph}"
		"rm -f {params.bigwigdir}/{params.bedgraph}_sorted"

## we computed the library sizes with 
## samtools view -c {input.bam}
## and computed the mean library size for each replicate pair 

## ------------------------------------------------------------------------------------ ##
## PCR duplicate removal
## ------------------------------------------------------------------------------------ ##

rule remove_PCR_duplicates:
	input:
		os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam")
	output:
		os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.bam")
	params:
		metrics_file = os.path.join(outputdir, "BAM_deduplicated/{sample}_marked_dup_metrics.txt")
	shell:
		"picard MarkDuplicates "
		"I={input} "
		"O={output} "
		"M={params.metrics_file} "
		"REMOVE_DUPLICATES=true "


## Index bam files
rule dedupidx:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated/{sample}.bam")
	output:
		os.path.join(outputdir, "BAM_deduplicated/{sample}.bam.bai")
	shell:
		"samtools index {input.bam}"



## ------------------------------------------------------------------------------------ ##
## Filter second read in pair
## ------------------------------------------------------------------------------------ ##
## the read-pair orientation is F2R1, the forward read is the second in the pair
rule filter_second_read:
	input:
		os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.bam")
	output:
		os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.bam")
	shell:
		"samtools view -f 128 -b -o {output} {input}"

## Convert BAM files to bigWig
rule bigwig_second_read:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated", "{sample}_deduplicated.r2.bam"),
		chrl = os.path.join(config["STARindex"], "chrNameLength.txt")
	output:
		os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.bw")
	params:
		bigwigdir = os.path.join(outputdir, "bigwig"),
		bedgraph = "{sample}_deduplicated.r2.bedGraph"
	log:
		os.path.join(outputdir, "logs", "bigwig_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "bigwig_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg > {params.bigwigdir}/{params.bedgraph}; "
		"LC_ALL='C' sort -k1,1 -k2,2n {params.bigwigdir}/{params.bedgraph} > {params.bigwigdir}/{params.bedgraph}_sorted; "
		"bedGraphToBigWig {params.bigwigdir}/{params.bedgraph}_sorted {input.chrl} {output}; "
		"rm -f {params.bigwigdir}/{params.bedgraph};6M1_200219_L1-2_deduplicated.r2.bam "
		"rm -f {params.bigwigdir}/{params.bedgraph}_sorted"


## ----------------------------------------------------------------------------
## Peak calling with CLIPper
## ----------------------------------------------------------------------------

## convert the Ensembl BAM file to UCSC
rule convert_BAM_UCSC:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.bam")
	output: 
		os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.ucsc.bam")
	shell:
		"samtools view -H {input.bam} | perl -lpe 's/SN:([0-9]+|[XY]|MT)\\b/SN:chr$1/' > {input.bam}_header.ucsc.sam; "
		"sed -i 's/SN:chrMT/SN:chrM/g' {input.bam}_header.ucsc.sam; "
		"samtools reheader {input.bam}_header.ucsc.sam {input.bam} > {output}"

## We use GRCh38_fixed which is the internal CLIPper annotation from which we had removed genome patches that were not present in our samples
rule clipper:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.ucsc.bam")
	output:
		peaks = os.path.join(outputdir, config["clipper"], "{sample}_deduplicated.r2.ucsc.clipper_peaks.bed")
	threads:
		config["ncores"]
	shell:
		"set +eu;"
		"source activate clipper3; "
		"clipper -b {input.bam} -s GRCh38_fixed --processors={threads} -o {output.peaks};"
		"source deactivate;"
		"set -eu"

## We are using `set +eu` to fix a conda error as suggested here: https://github.com/conda/conda/issues/8186 



## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
