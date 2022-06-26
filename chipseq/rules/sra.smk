rule SRAprefetch:
        output:
                "rawData/{SRA}/{SRA}.sra"
        shell:
                """
                prefetch -O rawData {wildcards.SRA}
                """

rule ParallelFastqDump:
    input:
        "rawData/{srr}/{srr}.sra"
    output:
        R1="rawData/{srr}_1.fastq.gz",
        R2="rawData/{srr}_2.fastq.gz"
    threads:
        64
    run:
        lib = sampleDF.loc[sampleDF["RawSample"] == wildcards.srr, "Library"].unique()[0]
        if lib == "Single":
            shell("""
            parallel-fastq-dump -t {threads} --split-files --gzip -s {input} -O rawData
            touch {output.R2}
            """)
        elif lib == "Paired":
            shell("""
            parallel-fastq-dump -t {threads} --split-files --gzip -s {input} -O rawData
            """)

