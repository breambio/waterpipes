import pandas as pd
#import yaml
#from pathlib import Path
#config = yaml.safe_load(Path("config/config.yaml").read_text())
samples = pd.read_table(config["SAMPLES"])
units = pd.read_table(config["UNITS"])

units["Raw"] = units["Name"] + "_" + units["Unit"].astype(str)

def get_lib(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_units(wildcards):
	return units.loc[units["Name"] == wildcards.raw, "Raw"].unique()

def get_fq1(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Fastq1"].unique()[0]

def get_fq2(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Fastq2"].unique()[0]

def get_contol(wildcards):
	 return samples.loc[samples["Name"] == wildcards.raw, "Control"].unique()[0]

# >>> `map.smk` functions >>>
def get_fqs(wildcards):
	source = wildcards.raw.find("SRR") != -1
	lib = get_lib(wildcards)
	if source:
		srr = get_fq1(wildcards)
		if lib == "Single":
			return f"sra-data/{srr}_1.fastq.gz"
		elif lib == "Paired":
			return f"sra-data/{srr}_1.fastq.gz", f"sra-data/{srr}_2.fastq.gz"
	elif source:
		fq1 = get_fq1(wildcards)
		fq2 = get_fq2(wildcards)
		if lib == "Single":
			return fq1
		elif lib == "Paired":
			return fq1, fq2

def get_filter_p(wildcards):
    lib = get_lib(wildcards)
    if lib == "Single":
        return "-F 3852"
    elif lib == "Paired":
        return "-F 3852 -f 2"

def get_reps(wildcards):
    reps = get_units(wildcards)
    return expand("results_{{ref}}/mapping/{rep}.filtered.bam", rep=reps)
# <<< `map.smk` functions <<<

# >>> `peak.smk` functions >>>
def get_macs_p(wildcards):
	lib = samples.loc[samples["Name"] == wildcards.raw, "Library"].unique()[0]
	input_c = get_contol(wildcards)
	if input_c == "-":
		param = ""
	else:
		param = f"-c results_{wildcards.ref}/mapping/{input_c}.final.bam"
	if lib == "Single":
		return param + " -f BAM"
	elif lib == "Paired":
		return param + " -f BAMPE"
# >>> `peak.smk` functions

ref = config["REF"]["NAME"]
bwNorm = config["OUTPUT"]["BW_NORMALIZATIONS"]

outputs = []

if config["OUTPUT"]["RUN"]["QC"]:
    outputs += ["qc/multiqc_report.html"]

if config["OUTPUT"]["RUN"]["PEAKS"]:
    outputs += [
        f"results_{ref}/peaks/{raw}_peaks.narrowPeak"
        for raw, c in zip(samples["Name"],samples["Control"])
        if c != "-"
    ]

if config["OUTPUT"]["RUN"]["BWS"]:
    outputs += [
        f"results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
        for raw in samples["Name"]
        for norm in bwNorm
    ]
