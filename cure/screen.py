import os
import shutil

import chem
import rcsb_pdb
import spreadsheetor

rcsb_apo_pdb_ids = {'AFABP': '3RZY', 'EFABP': '4LKP'}


def candidates():
    spreadsheet = spreadsheetor.get_spreadsheet()
    smiless_with_names_range = "A:B"
    names_smiless_value_range = spreadsheetor.get_values(spreadsheetor.SPREADSHEET_ID, smiless_with_names_range)
    names_smiless = names_smiless_value_range['values']
    names_smiless = list(zip(*names_smiless))
    data_start = 4
    return names_smiless[data_start:]


def write_candidates(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

    for candidate in candidates():
        chem.sdf(*candidate, path=path)


def write_targets(path):
    destination_path = os.path.join(path, '')
    if not os.path.exists(destination_path):
        os.makedirs(destination_path)

    for receptor in ['AFABP', 'EFABP']:
        rcsb_pdb.write_pdb(rcsb_pdb.apo_id(receptor), path=path)


def dock_candidates_to_targets():
    for receptor, pdb_id in rcsb_apo_pdb_ids.items():
        os.mkdir()
        for i, (name, smiles) in enumerate(names_smiless):
            gnina.dock(sdf, pdb_id)


if __name__ == "__main__":
    try:
        # write_candidates(os.path.join(os.getcwd(), "candidates/"))
        write_targets(os.path.join(os.getcwd(), "../targets/"))
        # dock_candidates_to_targets()
    except Exception as e:
        print(e)
