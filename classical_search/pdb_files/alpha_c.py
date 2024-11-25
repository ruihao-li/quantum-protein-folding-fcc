from Bio.PDB import *
import re

pdb_id = "6VHB"
parser = PDBParser()
io = PDBIO()
structure = parser.get_structure(pdb_id, pdb_id + ".pdb")
header = parser.get_header()
aa_seq = header["compound"]["1"]["molecule"]
# Amino acid three letter to one letter conversion
aa_dict = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
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

# Parse the structure to get the CA atoms and write them to an xyz file
fw = open(pdb_id + "_ca" + ".xyz", "w")
# Write the number of amino acids to the first line of the file
fw.write(str(len(aa_seq)) + "\n")
# Leave an empty line for comments
fw.write("\n")
# Write the one-letter amino acid ID and the x, y, z coordinates to the file line by line
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "CA":
                    x, y, z = atom.get_coord()
                    fw.write(
                        aa_dict[residue.get_resname()]
                        + " "
                        + str(x)
                        + " "
                        + str(y)
                        + " "
                        + str(z)
                        + "\n"
                    )


# for chain in structure.get_chains():
#     io.set_structure(chain)
#     io.save(chain.get_id() + ".pdb")


# def CA(pdbFileName):
#     fw = open("CA.xyz", "w")
#     fr = open(pdbFileName, "r")

#     for record in fr:
#         if re.search(r"^ATOM\s+\d+\s+CA\s+", record):
#             fw.write(record)

#     fw.close()
#     fr.close()


# CA("A.pdb")
