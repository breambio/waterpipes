
rule ParallelFastqDump:
    output:
        r1 = "sra-data/{srr}_1.fastq.gz",
        r2 = "sra-data/{srr}_2.fastq.gz"
    threads:
        64
    run:
        lib = units.loc[units["Fastq1"] == wildcards.srr, "Library"].unique()[0]
        if lib == "Single":
            shell("""
            parallel-fastq-dump -t {threads} --split-files --gzip -s {wildcards.srr} -O sra-data
            touch {output.r2}
            """)
        elif lib == "Paired":
            shell("""
            parallel-fastq-dump -t {threads} --split-files --gzip -s {wildcards.srr} -O sra-data
            """)
