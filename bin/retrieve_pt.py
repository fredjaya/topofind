#!/usr/bin/env python3
import json
import sys
import ast 

with open(sys.argv[1]) as j:
    PartitionedTrees = json.load(j)

run_name = sys.argv[2]
parts = run_name.split("_")[1:]
trees = [PartitionedTrees[k] for k in parts if k in PartitionedTrees]
trees = [t for t in trees if t != "None"]

if len(parts) != len(trees):
    print("No new partitions.")
    sys.exit(1)

else:
    with open(f"{run_name}.tre", "w") as f:
        for t in trees:
            f.write(f"{str(t)}\n")
