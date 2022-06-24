

import pandas as pd



sampleDF = pd.read_table("samples.tsv")

# TODO: is this necessary?
exclude = []

# TODO: merge issue could be solved like this.
on = "Raw"
genomecov = [
    f"results/bigwig/{raw}.{ref}.genomecov.{type}.bigWig"
    for raw in sampleDF.loc[ ~sampleDF[on].isin(exclude),on]
    for ref in ["hg38"]
    for type in ["rawcount"]
]

macs = [
    f"results/peaks/{raw}.{ref}_peaks.narrowPeak"
    for raw in sampleDF.loc[ ~sampleDF[on].isin(exclude),on]
    for ref in ["hg38"]
    if raw.find("input") == -1
]


qc = ["results/qc/multiqc_report.html"]


desiredOutputList = genomecov + macs 

