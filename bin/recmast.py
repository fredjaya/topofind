import sys
import ast
import re

def splitted_run_name(run_name):
    return run_name.split("_")

def update_run_name(splitted_run_name, p_to_split):
    """
    Output the new run_name when partitioning by p_to_split 
    e.g. "2_A_B", "A" -> "3_B_AA_AB"
    """
    new = ""
    n = int(splitted_run_name[0])
    n += 1
    p = splitted_run_name[1:]
    p = [i for i in p if i != p_to_split]
    p = "_".join(p)
    p1 = f"{p_to_split}A"
    p2 = f"{p_to_split}B"
    return f"{n}_{p}_{p1}_{p2}"

def form_names(new_list):
    formed = []
    for l in new_list:
        name = str(len(l)) + "_" + "_".join(l)
        formed.append(name)
    return formed

# TODO: move to test
assert update_run_name([2, 'A', 'B'], "B") == "3_A_BA_BB"
assert update_run_name([3, 'A', 'BA', 'BB'], "A") == "4_BA_BB_AA_AB"
