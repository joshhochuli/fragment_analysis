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


def get_sequence_from_file(filename):

    f = open(filename, "r")
    lines = f.readlines()
    sequence = lines[1].strip()
    return sequence

#reads through all files and makes one dataframe
#can save/load if parsing is slow
def get_dataframe():

    matrix = []

    f = glob.glob("../data/fragment_compound_pairs/*")

    for pair_dir in f: #read all files into one matrix

        stem = pair_dir.split("/")[-1]

        s = stem.split("_")
        protein_name = s[0]
        fragment_name = s[1]
        ligand_name = s[2]

        compound_filename = pair_dir + "/" + stem + "_substr_interact.txt"
        frag_filename = pair_dir + "/" + "_".join(stem.split("_")[:2]) + "_frag_interact.txt"

        frag_file = open(frag_filename, "r")
        compound_file = open(compound_filename, "r")

        compound_sequence_filename = pair_dir + "/" + stem + "_align_bindsite.fasta"
        fragment_sequence_filename = pair_dir + "/" + "_".join(stem.split("_")[:2]) + "_bindsite.fasta"

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

                    #lot of rereading the same fragment file, can optimize if needed
                    if is_fragment:
                        sequence = get_sequence_from_file(fragment_sequence_filename)
                    else:
                        sequence = get_sequence_from_file(compound_sequence_filename)

                    matrix.append((is_fragment, fragment_name,
                        protein_name, ligand_name, interaction_type, code,
                        res_id, sequence))

    df = pandas.DataFrame(matrix, columns = ["is_fragment",
        "fragment_name", "protein_name", "ligand_name", "interaction_type",
        "amino_acid", "residue_id", "sequence"])

    return df

#returns dataframe of (fragment_name, number of occurrences, number of associated ligands)
#sorted by number of associated ligands
def get_fragment_counts(df):

    fragment_names = df.fragment_name.unique()
    fragment_counts = []
    ligand_counts = []

    for fragment_name in fragment_names:
        fragment_count = len(df[(df.is_fragment == True) & (df.fragment_name ==
            fragment_name)].protein_name.unique())
        ligand_count = len(df[(df.is_fragment == False) & (df.fragment_name ==
            fragment_name)].ligand_name.unique())
        fragment_counts.append(fragment_count)
        ligand_counts.append(ligand_count)

    fragment_counts = np.array(fragment_counts)
    ligand_counts = np.array(ligand_counts)

    i = list(reversed(np.argsort(ligand_counts)))

    fragment_names = fragment_names[i]
    fragment_counts = fragment_counts[i]
    ligand_counts = ligand_counts[i]

    return pandas.DataFrame(list(zip(fragment_names, fragment_counts, ligand_counts)),
            columns = ["fragment_name", "fragment_count", "ligand_count"])

#uses interactions from file
#one plot for each interaction type
def make_interaction_histograms(df):

    fragment_names = df.fragment_name.unique()
    interaction_types = df.interaction_type.unique()

    #copied for each fragment
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

    for fragment_name in fragment_names:

        if not os.path.isdir(stem + fragment_name):
            os.mkdir(stem + fragment_name)

        for interaction_type in interaction_types:

            plt.clf()

            output_filename = stem + f"{fragment_name}/{fragment_name}_{interaction_type}_hist.png"

            fragment_aa_counts = copy.deepcopy(empty_aa_counts)
            compound_aa_counts = copy.deepcopy(empty_aa_counts)

            frag = df[
                    (df.fragment_name == fragment_name) & 
                    (df.interaction_type == interaction_type) & 
                    (df.is_fragment == True)]

            compound = df[
                    (df.fragment_name == fragment_name) & 
                    (df.interaction_type == interaction_type) & 
                    (df.is_fragment == False)]

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

            plt.bar(x_vals - width/2, frag_counts, width, label = f"Fragment (N = {frag_count})")
            plt.bar(x_vals + width/2, compound_counts, width, label = f"Compound (N = {compound_count})")
            plt.xticks(x_vals, aa_names)
            plt.title(f"{fragment_name}, {interaction_type}")
            plt.ylabel("Frequency")
            plt.legend()
            plt.savefig(output_filename, bbox_inches = "tight")


def main():

    df = get_dataframe()
    make_interaction_histograms(df)

if __name__ == "__main__":
    main()
