rule macs:
    input:
        get_macs_i
    output:
        "results_{ref}/peaks/{raw}_{q}_peaks.narrowPeak"
    threads:
        16
    params:
        get_macs_p
    shell:
        """
        macs3 callpeak \
          {params}  \
          -g hs -n results_{wildcards.ref}/peaks/{wildcards.raw}_{wildcards.q} -B -q {wildcards.q}
        """
