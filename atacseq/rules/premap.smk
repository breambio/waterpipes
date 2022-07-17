# TODO: fastq vs fq accession
def getLink(wildcards):
    r = wildcards.run
    link = sampleDF.loc[sampleDF["Raw"] == wildcards.raw, "Path"].to_list()[0]
    return f"{link}_{r}.fastq.gz"

rule cutadapt:
    input:
        r1="rawData/{raw}_1.fastq.gz",
        r2="rawData/{raw}_2.fastq.gz"
    output:
        r1="rawData/{raw}_1.trimmed.fastq.gz",
        r2="rawData/{raw}_2.trimmed.fastq.gz"
    params:
        fwd = config["ADAPTER_FWD"],
        rev = config["ADAPTER_REV"]
    threads:
        32
    run:
        lib = sampleDF.loc[sampleDF["RawSample"] == wildcards.raw, "Library"].unique()[0]
        if lib == "Single":
            shell("""
            cutadapt -a {params.fwd} -o {output.r1} -j {threads} {input}
            touch {output.R2}
            """)
        elif lib == "Paired":
            shell("""
            cutadapt -a {params.fwd} -A {params.rev} -o {output.r1} -p {output.r2} -j {threads} {input}
            """)

rule Links:
    input:
        getLink
    output:
        "links/{raw}_{run}.fastq.gz"
    shell:
        """
        ln -s {input} {output}
        """
