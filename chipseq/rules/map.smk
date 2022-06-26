def getFq(wildcards):
    trim = ""
    if config["CUT_ADAPTERS"]:
        trim = ".trimmed"

    lib = sampleDF.loc[sampleDF["Raw"].str.find(wildcards.raw) != -1, "Library"].unique()[0]
    if lib == "Single":
        return f"links/{{raw}}_1{trim}.fastq.gz"
    elif lib == "Paired":
        return f"links/{{raw}}_1{trim}.fastq.gz", f"links/{{raw}}_2{trim}.fastq.gz"

rule bwa_mem:
    input:
        getFq
    output:
        "results_{ref}/mapping/{raw}.raw.bam"
    params:
        idx = config["REF"]["BWA_IDX"]
    threads:
        32
    shell:
        """
        bwa mem -t {threads} \
            {params.idx} \
            {input} | \
            samtools view -bS - > {output}
        """

rule BamProcess:
    input:
        "results_{ref}/mapping/{raw}.raw.bam"
    output:
        "results_{ref}/mapping/{raw}.coorsorted.bam"
    threads:
        32
    params:
        config["BAMPROCESS_PARAMS"]#{config[PARAMS][BamProcess]}" # Note that you can add parameters as "-q 30 -F 1804"
    shell:
        """
        cat <(samtools view -H {input}) <(samtools view {params} {input}) | \
        samtools fixmate -m -@ {threads} - - | \
        samtools sort -@ {threads} -m 10G - | \
        samtools markdup -@ {threads} - {output}
        """

def getFilterParams(wildcards):
    lib = sampleDF.loc[sampleDF["Raw"].str.find(wildcards.raw) != -1, "Library"].unique()[0]
    if lib == "Single":
        return "-F 3852"
    elif lib == "Paired":
        return "-F 3852 -f 2"

rule Filter:
    input:
        "results_{ref}/mapping/{raw}.coorsorted.bam",
    output:
        "results_{ref}/mapping/{raw}.filtered.bam"
    threads:
        32
    params:
        config["REF"]["FA"],
        getFilterParams
    shell:
        """
        samtools view {input} | egrep -v "chrM" | \
        samtools view -b -@ {threads} -T {params} > {output}


        samtools index -@ {threads} {output}
        """

