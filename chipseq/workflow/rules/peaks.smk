rule Macs:
    input:
        get_macs_i
    output:
        "results_{ref}/peaks/{raw}_peaks.narrowPeak"
    threads:
        16
    params:
        get_macs_p
    shell:
        """
        macs3 callpeak \
          -t {input} \
          -g hs -n results_{wildcards.ref}/peaks/{wildcards.raw} -B -q 0.00001 {params}
        """
