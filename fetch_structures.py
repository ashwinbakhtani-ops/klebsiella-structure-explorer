from Bio.PDB import PDBList


STRUCTURES = {
    "NDM-1": "3SPU",
    "KPC-2": "2OV5",
    "OXA-48": "3HBR",
    "CTX-M-15": "4HBT",
}

pdbl = PDBList(verbose=False)
for name, pid in STRUCTURES.items():
    path = pdbl.retrieve_pdb_file(pid, pdir="structures", file_format="pdb")
    print(name, pid, "->", path)
