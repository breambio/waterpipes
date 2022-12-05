import pandas as pd
#config = yaml.safe_load(Path('config/config.yaml').read_text())
samples = pd.read_table(config["SAMPLES"])
units = pd.read_table(config["UNITS"])

units['Raw'] = units['Name'] + '_' + units['Unit'].astype(str)

def get_lib(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_source(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Source"].unique()[0]

def get_units(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Unit"].unique()

def get_fq1(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Fastq1"].unique()[0]

def get_fq2(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Fastq2"].unique()[0]

def get_contol(wildcards):
	raw_inp = samples.loc[samples["Name"] == wildcards.raw, "Control"].unique()[0]

# >>> `map.smk` functions >>>
def get_fqs(wildcards):
	source = get_source(wildcards)
	lib = get_lib(wildcards)
	if source == 'SRA':
		srr = get_fq1(wildcards)
		if lib == 'Single':
			return f"sra-data/{srr}_1.fastq.gz"
		elif lib == 'Paired':
			return f"sra-data/{srr}_1.fastq.gz", f"sra-data/{srr}_2.fastq.gz"
	elif source == 'Path':
		fq1 = get_fq1(wildcards)
		fq2 = get_fq2(wildcards)
		if lib == 'Single':
			return fq1
		elif lib == 'Paired':
			return fq1, fq2

def get_filters(wildcards):
    lib = get_library(wildcards)
    if lib == "Single":
        return "-F 3852"
    elif lib == "Paired":
        return "-F 3852 -f 2"

def get_reps(wildcards):
    reps = get_units(wildcards)
    return expand(f"results/mapping/{wildcards.raw}_{{rep}}.filtered.bam", rep=reps)
# <<< `map.smk` functions <<<

# >>> `peak.smk` functions >>>
def get_macs_p(wildcards):
    lib = samples.loc[samples["Name"] == wildcards.raw, "Library"].unique()[0]
    if lib == "Single":
        return "-f BAM"
    elif lib == "Paired":
        return "-f BAMPE"

def get_macs_i(wildcards):
	input_c = get_contol(wildcards)
	if input_c == '-':
		return "results_{ref}/mapping/{raw}.final.bam"
	else:
		return f"results_{{ref}}/mapping/{{raw}}.final.bam -c results_{{ref}}/mapping/{c}.final.bam"
# >>> `peak.smk` functions

ref = config["REF"]["NAME"]
bwNorm = config["OUTPUT"]["BW_NORMALIZATIONS"]

outputs = []

if config["OUTPUT"]["RUN"]["QC"]:
    outputs += ["qc/multiqc_report.html"]

if config["OUTPUT"]["RUN"]["PEAKS"]:
    outputs += [
        f"results_{ref}/peaks/{raw}_peaks.narrowPeak"
        for raw, c in zip(samples["Name"],samples["InputControl"])
        if c != '-'
    ]

if config["OUTPUT"]["RUN"]["BWS"]:
    outputs += [
        f"results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
        for raw in samples["Name"]
        for norm in bwNorm
    ]