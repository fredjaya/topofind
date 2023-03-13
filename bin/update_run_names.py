#!/usr/bin/env python3
import sys
import ast
import re

def recurse_trees(run_names):
    try:
        new_list=[]
        for old_run in run_names:
            old_parts = old_run.split("_")
            del old_parts[0]
            for t0 in old_parts:
                n = old_parts.copy()
                n.remove(t0)
                t1 = f"{t0}A"
                t2 = f"{t0}B"
                n += [t1, t2]
                new_list.append(n)
        return new_list
    except AttributeError:
        '''When tree_list==None, you cant .copy()'''
        pass

def form_names(new_list):
    formed = []
    for l in new_list:
        name = str(len(l)) + "_" + "_".join(l)
        formed.append(name)
    return formed

run_names = sys.argv[1]
if run_names == "null":
    print("2_A_B") 
else:
    run_names = ast.literal_eval(run_names)
    new_list = recurse_trees(run_names)
    print(form_names(new_list))
