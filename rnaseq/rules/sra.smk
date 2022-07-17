rule SRAprefetch:
        output:
                "rawData/{srr}/{srr}.sra"
        shell:
                """
                prefetch -O rawData {wildcards.srr}
                """

rule ParallelFastqDump:
    input:
        "rawData/{srr}/{srr}.sra"
    output:
        r1="rawData/{srr}_1.fastq.gz",
        r2="rawData/{srr}_2.fastq.gz"
    threads:
        64
    run:
        lib = sampleDF.loc[sampleDF["RawSample"] == wildcards.srr, "Library"].unique()[0]
        if lib == "Single":
            shell("""
            parallel-fastq-dump -t {threads} --split-files --gzip -s {input} -O rawData
            touch {output.r2}
            """)
        elif lib == "Paired":
            shell("""
            parallel-fastq-dump -t {threads} --split-files --gzip -s {input} -O rawData
            """)

