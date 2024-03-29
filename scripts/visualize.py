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

letter_to_code = dict([[v,k] for k,v in code_to_letter.items()])
letter_to_code["H"] = "HIS"

code_to_group = {
        "ARG":"Positively Charged",
        "HIS":"Positively Charged",
        "LYS":"Positively Charged",
        "ASP":"Negatively Charged",
        "GLU":"Negatively Charged",
        "SER":"Polar",
        "THR":"Polar",
        "ASN":"Polar",
        "GLN":"Polar",
        "ALA":"Hydrophobic",
        "ILE":"Hydrophobic",
        "LEU":"Hydrophobic",
        "MET":"Hydrophobic",
        "PHE":"Hydrophobic",
        "TRP":"Hydrophobic",
        "TYR":"Hydrophobic",
        "VAL":"Hydrophobic",
        "CYS":"CYS",
        "GLY":"GLY",
        "PRO":"PRO",
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

def make_sequence_histograms(df, normalize = True,
        count_multiple_occurrences = True, grouped = True):

    fragment_names = df.fragment_name.unique()

    empty_counts = {}
    if grouped:
        for group_name in code_to_group.values():
            empty_counts[group_name] = 0
    else:
        for aa_name in code_to_letter.keys():
            empty_counts[aa_name] = 0

    aa_names = list(code_to_letter.keys())
    group_names = ["Positively Charged", "Negatively Charged", "Polar",
    "Hydrophobic", "CYS", "GLY", "PRO"]

    #for plots
    if grouped:
        x_vals = np.array(range(len(group_names)))
    else:
        x_vals = np.array(range(len(aa_names)))
    width = 0.4

    stem = "../plots/"

    if not os.path.isdir(stem):
        print("'plots' directory does not exist in above directory, exiting...")

    if grouped:
        stem = stem + "grouped_sequence_histograms/"
    else:
        stem = stem + "sequence_histograms/"

    if not os.path.isdir(stem):
        os.mkdir(stem)

    f = plt.figure(figsize = (16,9))

    for fragment_name in fragment_names:

        plt.clf()

        if grouped:
            output_filename = stem + f"{fragment_name}_grouped_sequence_hist.png"
        else:
            output_filename = stem + f"{fragment_name}_sequence_hist.png"

        fragment_counts = copy.deepcopy(empty_counts)
        ligand_counts = copy.deepcopy(empty_counts)

        frag_sequences = df[
                (df.fragment_name == fragment_name) &
                (df.is_fragment == True)].sequence.unique()

        ligand_sequences = df[
                (df.fragment_name == fragment_name) &
                (df.is_fragment == False)].sequence.unique()

        frag_count = len(frag_sequences)
        ligand_count = len(ligand_sequences)

        if frag_count == 0 and ligand_count == 0:
            continue

        print(output_filename)

        for sequence in frag_sequences:
            if not count_multiple_occurrences:
                sequence = set(sequence)
            for aa in sequence:
                code = letter_to_code[aa]
                if grouped:
                    code = code_to_group[code]
                fragment_counts[code] += 1
        for sequence in ligand_sequences:
            if not count_multiple_occurrences:
                sequence = set(sequence)
            for aa in sequence:
                code = letter_to_code[aa]
                if grouped:
                    code = code_to_group[code]
                ligand_counts[code] += 1

        total_frag_count = sum([len(seq) for seq in frag_sequences])
        total_ligand_count = sum([len(seq) for seq in ligand_sequences])

        if(normalize):

            if(total_frag_count != 0):
                fragment_counts = {x: fragment_counts[x] / total_frag_count for x in
                        fragment_counts}

            if(total_ligand_count != 0):
                ligand_counts = {x: ligand_counts[x] / total_ligand_count for x in
                        ligand_counts}

        frag_counts_list = []
        ligand_counts_list = []

        if grouped:
            #be sure that group order matches count order
            for group in group_names:
                frag_counts_list.append(fragment_counts[group])
                ligand_counts_list.append(ligand_counts[group])

        else:
            #be sure that amino acid order matches count order
            for aa in aa_names:
                frag_counts_list.append(fragment_counts[aa])
                ligand_counts_list.append(ligand_counts[aa])

        plt.bar(x_vals - width/2, frag_counts_list, width, label = f"Fragment (N = {frag_count})", color = "black")
        plt.bar(x_vals + width/2, ligand_counts_list, width, label = f"Ligand (N = {ligand_count})", color = "gray")
        if grouped:
            plt.xticks(x_vals, group_names)
        else:
            plt.xticks(x_vals, aa_names)
        plt.title(f"{fragment_name}")
        plt.ylabel("Frequency")
        plt.legend()
        plt.savefig(output_filename, bbox_inches = "tight")


#uses interactions from file
#one plot for each interaction type
def make_interaction_histograms(df, normalize = True):

    fragment_names = df.fragment_name.unique()
    interaction_types = df.interaction_type.unique()

    #copied for each fragment
    empty_aa_counts = {}
    for aa_name in code_to_letter.keys():
        empty_aa_counts[aa_name] = 0

    aa_names = list(code_to_letter.keys())

    #for plots
    x_vals = np.array(range(len(aa_names)))
    width = 0.4

    stem = "../plots/"

    if not os.path.isdir(stem):
        print("'plots' directory does not exist in above directory, exiting...")

    stem = stem + "interaction_histograms/"

    if not os.path.isdir(stem):
        os.mkdir(stem)

    f = plt.figure(figsize = (16,9))

    for fragment_name in fragment_names:

        if not os.path.isdir(stem + fragment_name):
            os.mkdir(stem + fragment_name)

        for interaction_type in interaction_types:

            plt.clf()

            output_filename = stem + f"{fragment_name}/{fragment_name}_{interaction_type}_hist.png"

            fragment_aa_counts = copy.deepcopy(empty_aa_counts)
            ligand_aa_counts = copy.deepcopy(empty_aa_counts)

            frag = df[
                    (df.fragment_name == fragment_name) &
                    (df.interaction_type == interaction_type) &
                    (df.is_fragment == True)]

            ligand = df[
                    (df.fragment_name == fragment_name) &
                    (df.interaction_type == interaction_type) &
                    (df.is_fragment == False)]

            frag_count = len(frag.protein_name.unique())
            ligand_count = len(ligand.ligand_name.unique())

            if frag_count == 0 and ligand_count == 0:
                continue

            print(output_filename)

            for aa in frag.amino_acid:
                fragment_aa_counts[aa] += 1
            for aa in ligand.amino_acid:
                ligand_aa_counts[aa] += 1

            total_frag_count = len(frag.amino_acid)
            total_ligand_count = len(ligand.amino_acid)

            if(normalize):

                if(total_frag_count != 0):
                    fragment_aa_counts = {x: fragment_aa_counts[x] / total_frag_count for x in
                            fragment_aa_counts}

                if(total_ligand_count != 0):
                    ligand_aa_counts = {x: ligand_aa_counts[x] / total_ligand_count for x in
                            ligand_aa_counts}


            frag_counts_list = []
            ligand_counts_list = []

            #be sure that amino acid order matches count order
            for aa in aa_names:
                frag_counts_list.append(fragment_aa_counts[aa])
                ligand_counts_list.append(ligand_aa_counts[aa])

            plt.bar(x_vals - width/2, frag_counts_list, width, label = f"Fragment (N = {frag_count})", color = "black")
            plt.bar(x_vals + width/2, ligand_counts_list, width, label = f"Ligand (N = {ligand_count})", color = "gray")
            plt.xticks(x_vals, aa_names)
            plt.title(f"{fragment_name}, {interaction_type}")
            plt.ylabel("Frequency")
            plt.legend()
            plt.savefig(output_filename, bbox_inches = "tight")


def main():

    df = get_dataframe()
    make_sequence_histograms(df, grouped = False)

if __name__ == "__main__":
    main()
