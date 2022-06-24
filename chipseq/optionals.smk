
i = 0
################################################################################
## FASTQC
################################################################################


rule Fastqc:
    input:
        "links/{raw}.{ref}.{run}.fastq.gz"
    output:
        "results/qc/{raw}.{ref}.{run}_fastqc.zip"
    shell:
        """
        fastqc {input} -o results/qc
        """


def getMultiqc(wildcards):
    out = []
    for raw in sampleDF["Raw"]:
        fq = f"results/qc/{raw}"
        ref = getRef(raw)
        fq += f".{ref}"
        lib = getLib(raw)
        if lib == "Single":
            fq += f".1_fastqc.zip"
            out.append(fq)
        elif lib == "Paired":
            fq1 = fq + f".1_fastqc.zip"
            out.append(fq1)
            fq2 = fq + f".2_fastqc.zip"
            out.append(fq2)
    return expand(out)


rule Multiqc:
    input:
        getMultiqc
    output:
        "results/qc/multiqc_report.html"
    shell:
        """
        cd results/qc/ && multiqc .
        """



################################################################################
## MERGE
################################################################################



def getRepsToMerge(wildcards):
    raw = sampleDF.loc[
    (sampleDF["Cell"] == wildcards.cell) &
    (sampleDF["Target"] == wildcards.target) &
    (sampleDF["Lab"] == wildcards.lab),
    "Raw"].to_list()
    return expand("results/mapping/{raw}.{{ref}}.filtered.bam", raw=raw)




rule BamMerge:
    input:
        getRepsToMerge
    output:
        "results/mapping/{cell}.{lab}.{target}.merges.{ref}.bam"
    threads:
        32
    shell:
        """
        samtools merge -@ {threads} {output} {input}
        """



################################################################################
## DEEPTOOLS BIGWIGS
################################################################################

def getInput(wildcards):
    raw = wildcards.raw.rsplit(".", 1)[0]
    inp = sampleDF.loc[(sampleDF[on].str.find(raw) != -1), "Control"].unique()[0]
    return f"results/mapping/{inp}.{{ref}}.filtered.bam"


def getInputBai(wildcards):
    raw = wildcards.raw.rsplit(".", 1)[0]
    inp = sampleDF.loc[(sampleDF[on].str.find(raw) != -1), "Control"].unique()[0]
    return f"results/mapping/{inp}.{{ref}}.bam.bai"




def getParamsBamCom(wildcards):
    raw = wildcards.raw.rsplit(".", 1)[0]
    lib = getLib(raw)
    if lib == "Single":
        return " --extendReads 150 "
    elif lib == "Paired":
        return "--extendReads "



rule BamCom:
    input:
        bam="results/mapping/{raw}.{ref}.bam",
        bai="results/mapping/{raw}.{ref}.bam.bai",
        inp=getInput,
        bii=getInputBai
    output:
        bw="results/bigwig/{raw}.{ref}.{type}.bamcom.bigWig"
    params:
        getParamsBamCom
    threads:
        32
    shell:
        """
        bamCompare -b1 {input.bam} -b2 {input.inp} -o {output.bw} \
        -p {threads} \
        --scaleFactorsMethod None --normalizeUsing {wildcards.type} \
        {params}
        """

rule BamCov:
    input:
        bam="results/mapping/{raw}.{ref}.bam",
    output:
        bw="results/bigwig/{raw}.{ref}.{type}.bamcov.bigWig"
    params:
        getParamsBamCom
    threads:
        32
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} \
        -p {threads} \
        --normalizeUsing {wildcards.type} \
        {params}
        """




################################################################################
## MACS
################################################################################


def getParamsMACS(wildcards):
    raw = wildcards.raw.rsplit(".", 1)[0]
    lib = getLib(raw)
    if lib == "Single":
        return "-f BAM"
    elif lib == "Paired":
        return "-f BAMPE"



rule Macs:
    input:
        bam="results/mapping/{raw}.{ref}.filtered.bam",
        inp=getInput
    output:
        "results/peaks/{raw}.{ref}_peaks.narrowPeak"
    threads:
        16
    params:
        getParamsMACS
    shell:
        """
        macs3 callpeak \
          -t {input.bam} \
          -c {input.inp} \
          -g hs -n results/peaks/{wildcards.raw}.{wildcards.ref} -B -q 0.00001 {params}
        """



################################################################################
## ROSE
################################################################################

rule npk2gff:
    input:
        "results/peaks/{raw}.{ref}_peaks.narrowPeak"
    output:
        "results/gff/{raw}.{ref}.gff"
    shell:
        """
        awk -F'\t' 'BEGIN {{OFS = FS}} {{if ($9 > 20) print $1,$4,"",$2,$3,"",".","",$4}}' {input} > {output}
        """




rule Rose:
    input:
        bam="results/mapping/{raw}.{ref}.filtered.bam",
        inp=getInput,
        gff="results/gff/{raw}.{ref}.gff"
    output:
        "results/SE/{raw}.{ref}_SuperEnhancers.table.txt"
    threads:
        4
    shell:
        """
        ROSE_main.py -g {wildcards.ref} \
        -i {input.gff} -r {input.bam} -c {input.inp} \
        -o {wildcards.raw}.{wildcards.ref}
        """

