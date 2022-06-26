


def getMACSinput(wildcards):
    inp = ["results_{ref}/mapping/{raw}.filtered.bam"]
    if wildcards.type in ("tVSc"):
        raw_inp = sampleDF.loc[(sampleDF["Raw"].str.find(wildcards.raw) != -1), "Control"].unique()[0]
        inp.append(f"results_{{ref}}/mapping/{raw_inp}.filtered.bam")
    return inp


rule Macs:
        input:
            getMACSinput
        output:
            "results_{ref}/peaks/{raw}.{type}_peaks.narrowPeak",
            temp("results_{ref}/peaks/{raw}.{type}_treat_pileup.bdg")
        threads:
            16
        params:
            config["MACS3_PARAMS"]
        run:
            if wildcards.type == "nomodel":
                shell("""
                macs3 callpeak -t {input} \
                -n results_{wildcards.ref}/peaks/{wildcards.raw}.nomodel \
                {params} -g hs -q 0.01 --nomodel --shift -75 --extsize 150 \
                --keep-dup all -B --SPMR
                """)
            elif wildcards.type == "callsummit":
                shell("""
                macs3 callpeak -t {input} \
                -n results_{wildcards.ref}/peaks/{wildcards.raw}.callsummit \
                {params} -g hs -q 0.01 --call-summits -B
                """)
            elif wildcards.type == "tVSc":
                shell("""
                macs3 callpeak \
                  -t {input[0]} \
                  -c {input[1]} \
                  -g hs -n results_{wildcards.ref}/peaks/{wildcards.raw}.tVSc -B -q 0.01 {params}
                """)
            elif wildcards.type == "genrich":
                shell("""
                samtools view -h {input} \
                | Genrich -t - -o {output} -j
                """)



# TODO: genomeAnnotations issue
# TODO: pipe sort (clip -> sort)
rule MACSbw:
    input:
        "results_{ref}/peaks/{raw}.{type}_treat_pileup.bdg"
    output:
        bg=temp("results_{ref}/bw/{raw}.{type}.bg")
        bw="results_{ref}/bw/{raw}.{type}.bw"
    params:
        chrSizes=config[wildcards.ref]["CHROM_SIZES"]
    shell:
        """
        bedtools slop -i {input} -g {params.chrSizes} -b 0 \
        | bedClip stdin {params.chrSizes} stdout \
        | sort -k1,1 -k2,2n > {output.bg}

        bedGraphToBigWig {output.bg} {params.chrSizes} {output.bw}
        """

