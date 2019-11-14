import os
import glob
import pandas
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import stats

code_to_letter = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "ASX": "B",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLX": "Z",
        "GLY": "G",
        "HIS": "H",
        "HID": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        }


f = glob.glob("../data/fragment_compound_pairs/*")

#values should be ({overall aa counts}, key appearance counts)
fragment_to_aa = {}
compound_to_aa = {}

overall_matrix = []

fragment_to_protein = {}
for pair_dir in f:

    stem = pair_dir.split("/")[-1]

    s = stem.split("_")
    protein_name = s[0]
    fragment_name = s[1]
    ligand_name = s[2]

    if fragment_name not in fragment_to_protein:
        fragment_to_protein[fragment_name] = set()

    fragment_to_protein[fragment_name].add(protein_name)


    compound_filename = pair_dir + "/" + stem + "_substr_interact.txt"
    frag_filename = pair_dir + "/" + "_".join(stem.split("_")[:2]) + "_frag_interact.txt"

    align_filename = pair_dir + "/" + stem + "_align_bindsite.fasta"
    bindsite_filename = pair_dir + "/" + "_".join(stem.split("_")[:2]) + "_bindsite.fasta"


    compound_file = open(compound_filename, "r")
    frag_file = open(frag_filename, "r")

    align_file = open(align_filename, "r")
    bindsite_file = open(bindsite_filename, "r")

    for is_fragment, interact_file in [(True, frag_file),(False, compound_file)]:
        for line in interact_file:
            if "(" in line:
                s = line.split()
                interaction_type = s[0]
                count = int(s[3])
                id_start = line.find(":") + 1
                code_start = line.find("(")+1
                res_id = line[id_start:code_start - 1]
                code = line[code_start:code_start + 3]
                letter = code_to_letter[code]
                overall_matrix.append((is_fragment, fragment_name, protein_name, ligand_name, interaction_type, code, res_id))

df = pandas.DataFrame(overall_matrix, columns = ["is_fragment", "fragment_name", "protein_name", "ligand_name", "interaction_type", "amino_acid", "residue_id"])


fragment_names = df.fragment_name.unique()
interaction_types = df.interaction_type.unique()

empty_aa_counts = {}
for fragment_name in fragment_names:
    empty_aa_counts[fragment_name] = {}
    for aa_name in code_to_letter.keys():
        empty_aa_counts[fragment_name][aa_name] = 0



aa_names = list(code_to_letter.keys())
x_vals = np.array(range(len(aa_names)))
width = 0.4

normalize = True

stem = "../plots/"

if not os.path.isdir(stem):
    print("'plots' directory does not exist in above directory, exiting...")

f = plt.figure(figsize = (16,9))

x = []
y = []

for fragment_name in fragment_names:
    ligand_count = len(df[(df.fragment_name == fragment_name)].ligand_name.unique())
    x.append(fragment_name)
    y.append(ligand_count)
    continue


    if not os.path.isdir(stem + fragment_name):
        os.mkdir(stem + fragment_name)
    for interaction_type in interaction_types:


        plt.clf()

        output_filename = stem + f"{fragment_name}/{fragment_name}_{interaction_type}_hist.png"

        fragment_aa_counts = copy.deepcopy(empty_aa_counts)
        compound_aa_counts = copy.deepcopy(empty_aa_counts)

        frag = df[(df.fragment_name == fragment_name) & (df.interaction_type == interaction_type) & (df.is_fragment == True)]
        frag = frag.sort_values(["is_fragment", "protein_name", "ligand_name", "residue_id"])
        compound = df[(df.fragment_name == fragment_name) & (df.interaction_type == interaction_type) & (df.is_fragment == False)]

        frag_count = len(frag.protein_name.unique())
        compound_count = len(compound.ligand_name.unique())

        if frag_count == 0 and compound_count == 0:
            continue

        print(output_filename)

        for val in frag.amino_acid:
            fragment_aa_counts[fragment_name][val] += 1
        for val in compound.amino_acid:
            compound_aa_counts[fragment_name][val] += 1

        frag_counts = np.array(list(fragment_aa_counts[fragment_name].values()))
        compound_counts = np.array(list(compound_aa_counts[fragment_name].values()))

        if(normalize):
            if(sum(frag_counts) != 0):
                frag_counts = frag_counts / sum(frag_counts)
            if(sum(compound_counts) != 0):
                compound_counts = compound_counts / sum(compound_counts)

        #kld = stats.entropy(frag_counts, compound_counts)

        plt.bar(x_vals - width/2, frag_counts, width, label = f"Fragment (N = {frag_count})")
        plt.bar(x_vals + width/2, compound_counts, width, label = f"Compound (N = {compound_count})")
        plt.xticks(x_vals, aa_names)
        plt.title(f"{fragment_name}, {interaction_type}")
        plt.ylabel("Frequency")
        plt.legend()
        plt.savefig(output_filename, bbox_inches = "tight")

x = np.array(x)
y = np.array(y)
i = np.argsort(y)
x = x[i]
y = y[i]
for i in range(len(x)):
    print(f"{x[i]},{y[i]}")
