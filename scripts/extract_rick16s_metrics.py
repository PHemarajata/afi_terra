import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--bam")
parser.add_argument("--panel")
parser.add_argument("--out")

args = parser.parse_args()

data = [
    {"genus": "Orientia", "mapped_reads": 100, "max_breadth": 0.25},
    {"genus": "Rickettsia", "mapped_reads": 20, "max_breadth": 0.05}
]

pd.DataFrame(data).to_csv(args.out, sep="\t", index=False)
