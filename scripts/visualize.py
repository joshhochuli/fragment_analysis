import glob


f = glob.glob("../data/fragment_compound_pairs/*")

for pair_dir in f:

    stem = pair_dir.split("/")[-1]
    print(stem)

    substr_filename = pair_dir + "/" + stem + "_substr_interact.txt"
    frag_filename = pair_dir + "/" + "_".join(stem.split("_")[:2]) + "_frag_interact.txt"

    align_filename = pair_dir + "/" + stem + "_align_bindsite.fasta"
    bindsite_filename = pair_dir + "/" + "_".join(stem.split("_")[:2]) + "_bindsite.fasta"

    substr_file = open(substr_filename, "r")
    frag_file = open(frag_filename, "r")

    align_file = open(align_filename, "r")
    bindsite_file = open(bindsite_filename, "r")

    letter_to_code = {
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

    for seq_file, interact_file in [(bindsite_file, frag_file),(align_file, substr_file)]:
        seq = seq_file.readlines()[1].strip()
        print("\nSequence: " + seq)

        print("Interaction, Res_id, Res, Count")
        for line in interact_file:
            if "(" in line:
                s = line.split()
                interaction_type = s[0]
                count = int(s[3])
                id_start = line.find(":") + 1
                code_start = line.find("(")+1
                res_id = line[id_start:code_start - 1]
                code = line[code_start:code_start + 3]
                print(f"{interaction_type}, {res_id}, {code}, {count}")
                letter = letter_to_code[code]
                if letter not in seq:
                    print(pair_dir)
                    print(seq)
                    print(letter)
                    exit()

    break


