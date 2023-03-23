import os
import glob

import urllib.request

import spreadsheetor

def get_names_smiless():
    spreadsheet = spreadsheetor.get_spreadsheet()
    smiless_with_names_range = "A:B"
    names_smiless_value_range = spreadsheetor.get_values(spreadsheetor.SPREADSHEET_ID, smiless_with_names_range)
    names_smiless = names_smiless_value_range['values']
    return names_smiless

def get_rcsb_pdb(id):
    rcsb_structure_url = 'https://www.rcsb.org/structure/'
    urllib.request.urlretrieve(rcsb_structure_url + id)

def dock(ligand_smiles, receptor_path):
    pass

if __name__ == "__main__":
    names_smiless = get_names_smiless()
    smiless_index = 1
    rcsb_apo_pdb_ids = {'AFABP': '3RZY', 'EFABP': '4LKP'}
    for receptor, pdb_id in rcsb_apo_pdb_ids.items():
        for smiles in names_smiless[smiless_index]:
            dock(smiles, pdb_id)
