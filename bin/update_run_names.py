#!/usr/bin/env python3
from recmast import update_run_name
import sys
import re

run_name = sys.argv[1]
aln_paths = sys.argv[2]

"""
Convert aln_paths to list
"""
aln_paths = re.sub("\[|\]", "", aln_paths)
aln_paths = aln_paths.split(",")

"""
Remove trailing/leading spaces from elements
"""
aln_paths = [a.strip() for a in aln_paths]

"""
Sort aln_paths so they correspond to the parts
"""
aln_paths = sorted(aln_paths, key = lambda x: x.rsplit('/', 1)[-1])

"""
For each part, generate new names for the next run.
Match with the corresponding part name and save to dict, or throw error
"""
rs = run_name.split("_")
assert len(rs[1:]) == len(aln_paths)

NewRuns = {}
for p, aln in zip(rs[1:], aln_paths):
    new_name = update_run_name(rs, p)
    NewRuns[new_name] = aln

with open("NewRuns.tsv", "w") as f:
    for run, aln in NewRuns.items():
        f.write(f"{run}\t{aln}\n")
