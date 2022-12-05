rule GenomeCov:
    input:
        bam = "results_{ref}/mapping/{raw}.final.bam",
    output:
        bg = temp("results_{ref}/bigwig/{raw}.genomecov.{norm}.bg"),
        bw = "results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
    params:
        chrSizes = config["REF"]["CHROM_SIZES"]
    threads:
        16
    run:
        if wildcards.norm == "RPM":
            shell("""
            N=`samtools view -@{threads} {input.bam} | wc -l`

            bedtools genomecov \
            -scale `bc -l <<< 1000000/$N` \
            -bg -ibam {input.bam} | sort -k1,1 -k2,2n --parallel={threads} > {output.bg}

            bedGraphToBigWig {output.bg} {params.chrSizes} {output.bw}
            """)
        elif wildcards.norm == "rawcount":
            shell("""
            bedtools genomecov \
            -bg -ibam {input.bam} | sort -k1,1 -k2,2n --parallel={threads} > {output.bg}

            bedGraphToBigWig {output.bg} {params.chrSizes} {output.bw}
            """)
