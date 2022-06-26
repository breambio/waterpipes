def getMACSparam(wildcards):
    lib = sampleDF.loc[sampleDF["Raw"].str.find(wildcards.raw) != -1, "Library"].unique()[0]
    if lib == "Single":
        return "-f BAM"
    elif lib == "Paired":
        return "-f BAMPE"

def getMACSinput(wildcards):
    inp = ["results_{ref}/mapping/{raw}.filtered.bam"]
    raw_inp = sampleDF.loc[(sampleDF["Raw"].str.find(wildcards.raw) != -1), "Control"].unique()[0]
    inp.append(f"results_{{ref}}/mapping/{raw_inp}.filtered.bam")
    return inp

rule Macs:
    input:
        getMACSinput
    output:
        "results_{ref}/peaks/{raw}_peaks.narrowPeak"
    threads:
        16
    params:
        getMACSparam
    shell:
        """
        macs3 callpeak \
          -t {input[0]} \
          -c {input[1]} \
          -g hs -n results_{wildcards.ref}/peaks/{wildcards.raw} -B -q 0.00001 {params}
        """

