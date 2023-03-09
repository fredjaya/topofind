#!/usr/bin/env python3
import sys
import re
from collections import OrderedDict
import json

"""
Read args
"""
part_names = sys.argv[1].split("_")[2:]
trees = sorted(sys.argv[2:])

"""
Read each tree file and save to list
"""
nwk = []
for t in trees:
    with open(t, "r") as f:
        """
        Remove branch lengths
        """
        f = f.read().strip()
        f = re.sub(":0\.\d+", "", f)
        nwk.append(f)

"""
Add partition name and corresponding tree to dict
"""
PartitionedTrees = OrderedDict()
PartitionedTrees = {k:v for k,v in zip(part_names, nwk)}

"""
Save to JSON
"""
with open("PartitionedTrees.json", "w") as j:
    json.dump(PartitionedTrees, j)

