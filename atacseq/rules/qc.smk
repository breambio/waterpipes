





rule Fastqc:
    input:
        "links/{raw}_{run}.fastq.gz"
    output:
        "qc/{raw}.{run}_fastqc.zip"
    shell:
        """
        fastqc {input} -o qc
        """




def getMultiqc(wildcards):
    out = []
    for raw in sampleDF["Raw"]:
        fq = f"qc/{raw}"
        lib = sampleDF.loc[sampleDF["Raw"].str.find(raw) != -1, "Library"].unique()[0]
        if lib == "Single":
            fq += f".1_fastqc.zip"
            out.append(fq)
        elif lib == "Paired":
            fq1 = fq + f"_1_fastqc.zip"
            out.append(fq1)
            fq2 = fq + f"_2_fastqc.zip"
            out.append(fq2)
    return expand(out)


rule Multiqc:
    input:
        getMultiqc
    output:
        "qc/multiqc_report.html"
    shell:
        """
        cd qc/ && multiqc .
        """

