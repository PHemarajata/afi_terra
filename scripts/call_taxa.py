import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sample")
parser.add_argument("--metrics")
parser.add_argument("--ntc")
parser.add_argument("--out")

args = parser.parse_args()

metrics = pd.read_csv(args.metrics, sep="\t")
ntc = pd.read_csv(args.ntc, sep="\t")

results = []

for _, row in metrics.iterrows():

    reads = row["mapped_reads"]
    breadth = row["max_breadth"]

    ntc_reads = ntc.loc[
        ntc["genus"] == row["genus"],
        "mapped_reads"
    ].values[0]

    if reads >= 1000 and breadth >= 0.25 and reads >= 10 * ntc_reads:
        call = "Confirmed"

    elif reads >= 50 and breadth >= 0.20 and reads > ntc_reads:
        call = "Probable"

    elif reads >= 50 and reads <= ntc_reads:
        call = "Not_confirmed"

    else:
        call = "Negative"

    results.append({
        "sample": args.sample,
        "genus": row["genus"],
        "reads": reads,
        "breadth": breadth,
        "ntc_reads": ntc_reads,
        "call": call
    })

pd.DataFrame(results).to_csv(args.out, sep="\t", index=False)
