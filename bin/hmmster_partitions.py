#!/usr/bin/env python
import re
import sys

partitions = dict()

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("["):
            line = re.sub(",", "-", line.strip())
            line = re.sub("\[|\]", "", line)
            v,k = line.split("\t")
            if k in partitions:
                partitions[k].append(v)
            else:
                partitions[k] = [v]
            
sorted_parts = {}
for key in sorted(partitions):
    sorted_parts[key] = partitions[key]
    
as_txt = ""
for k, v in sorted_parts.items():
    as_txt += "class_" + k + " = " + " ".join(v) + "\n"

with open("hmmster.partitions_amas", "w") as f:
    f.write(as_txt)
