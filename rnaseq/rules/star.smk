def getFq(wildcards):
    trim = ""
    if config["CUT_ADAPTERS"]:
        trim = ".trimmed"

    lib = sampleDF.loc[sampleDF["Raw"].str.find(wildcards.raw) != -1, "Library"].unique()[0]
    if lib == "Single":
        return f"links/{{raw}}_1{trim}.fastq.gz"
    elif lib == "Paired":
        return f"links/{{raw}}_1{trim}.fastq.gz", f"links/{{raw}}_2{trim}.fastq.gz"

rule star:
    input:
        getFq
    output:
        "results_{ref}/star/{raw}Aligned.sortedByCoord.out.bam",
        "results_{ref}/star/{raw}Aligned.toTranscriptome.out.bam"
    params:
        staridx = config["REF"]["STAR_IDX"]
    threads:
        8
    shell:
        """
        STAR --genomeDir {params.staridx} \
        --runThreadN {threads} \
        --readFilesIn {input} \
        --outFileNamePrefix results_{wildcards.ref}/star/{wildcards.raw} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --readFilesCommand zcat\
        --quantMode TranscriptomeSAM

        samtools index -@{threads} {output}
        """

rule featureCounts:
    input:
        expand("results_{{ref}}/star/{raw}Aligned.sortedByCoord.out.bam", raw = sampleDF["Raw"].tolist())
    output:
        "results_{ref}/countTable.featureCounts.tsv"
    params:
        gtf = config["REF"]["GTF"],
        fea = config["OUTPUT"]["FEATURECOUNTS_FEATURE"],
        att = config["OUTPUT"]["FEATURECOUNTS_ATTRIBUTE"]
    threads:
        64
    shell:
        """
        featureCounts -T {threads} -t {param.fea} -g {param.att} -a {param.gtf} -o featureCounts.tsv {input}

        cut -f1,7- featureCounts.tsv | tail -n+2 > {output}

        rm featureCounts.tsv
        """

rule Salmon:
    input:
        "results_{ref}/star/{raw}Aligned.toTranscriptome.out.bam"
    output:
        "results_{ref}/salmon/{raw}/quant.sf"
    params:
        attotfa = config["REF"]["TRANSCRIPTS_FA"]
    threads:
        4
    shell:
        """
        salmon quant -t {params.attotfa} -l A -a {input} -o {output}
        """

# TODO: Segmentation fault
rule SalmonALL:
    input:
        expand("results_{{ref}}/star/{raw}Aligned.toTranscriptome.out.bam", raw = sampleDF["Raw"].tolist())
    output:
        "results_{ref}/salmon.all/quant.sf"
    params:
        attotfa = config["REF"]["TRANSCRIPTS_FA"]
    threads:
        4
    shell:
        """
        salmon quant -t {params.attotfa} -l A -a {input} -o {output}
        """

