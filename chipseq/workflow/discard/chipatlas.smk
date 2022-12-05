rule ChipAtlas:
    output:
        "results_{ref}/cabigwig/{srx}.bw"
    threads:
        16
    shell:
        """
        url=http://dbarchive.biosciencedbc.jp/kyushu-u/{wildcards.ref}/eachData/bw/{wildcards.srx}.bw
        if [[ $(wget $url -O-) ]] 2>/dev/null; then
            /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
            $url \
            -o {output}
        else
            touch {output}
        fi
        """

rule ChipAtlasBed:
    output:
        "results_{ref}/cabeds/{srx}.{threshold}.bed"
    threads:
        16
    shell:
        """
        url=http://dbarchive.biosciencedbc.jp/kyushu-u/{wildcards.ref}/eachData/bed{wildcards.threshold}/{wildcards.srx}.{wildcards.threshold}.bed
        if [[ $(wget $url -O-) ]] 2>/dev/null;
            then
                /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
                $url \
                -o {output}
            else
                touch {output}
        fi
        """