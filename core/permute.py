import itertools
import os
import glob
from ase import Atoms
# from dscribe.descriptors import CoulombMatrix
from sklearn.preprocessing import StandardScaler
import tqdm
import numpy as np
import random

import core.tools as tools
import core.representatives as rep
import core.connectivity as net
os.environ["OMP_NUM_THREADS"] = "1"
def get_templates(params: dict) -> list:
    natoms = params["NUMELEM1"]+params["NUMELEM2"]
    combinations = itertools.combinations(range(natoms), params["NUMELEM1"])
    return list(combinations)

# def gen_permutations(params: dict, files: str, combinations: list) -> None:
#     """
#     Generate all permutations of the alloys
#     """

#     X = []
#     E = []

#     for file in files:
#         natoms, atomtypes, coords, energy = tools.xyzRead(file)
#         X.append([atomtypes, coords])
#         E.append(energy)

#     mols = [Atoms(positions=coordinates, symbols=symbols) for (symbols, coordinates) in X]
#     cm = CoulombMatrix(n_atoms_max=natoms, permutation="eigenspectrum")
#     coulomb = cm.create(mols)

#     # natoms = params["NUMELEM1"]+params["NUMELEM2"]
#     tmpfolder = params["MOD2"]["TMP_FOLDER"]
#     base_name = os.path.basename(file)
#     base_name = os.path.splitext(base_name)[0]

#     os.system(f"rm -rf {tmpfolder}/unfiltered/*")
#     os.system(f"rm -rf {tmpfolder}/filtered/*")

#     # combinations = itertools.combinations(range(natoms), params["NUMELEM1"])

#     os.system(f"rm -rf {tmpfolder}/unfiltered/*")
#     for index, combination in enumerate(combinations):
#         outfile = tmpfolder+"/unfiltered/"+base_name+"_"+str(index)+".xyz"
#         os.system(f"cp {file} {outfile}")
#         tools.replace_atoms_xyz(params, outfile, combination)
    
#     net.integrity_test(tmpfolder) # unfiltered -> filtered

#     rep.get_representatives(params["MOD2"])
#     os.system(f"mv {tmpfolder}/filtered/*.xyz {tmpfolder}/filtered_all/")


# VERSION 02
# def gen_permutations(params: dict, files: str) -> None:
#     """
#     Generate all permutations of the alloys
#     """
#     atom1 = params["ELEM1"]
#     atom2 = params["ELEM2"]
#     natoms = params["NUMELEM1"]+params["NUMELEM2"]
#     outfolder = params["MOD2"]["TMP_FOLDER"]+"/unfiltered/"
#     tmpfolder = params["MOD2"]["TMP_FOLDER"]

#     combinations = list(itertools.combinations(range(natoms), params["NUMELEM1"]))

#     natoms, atomtypes, coords, energy = tools.xyzRead(files[0])
#     list_atoms = []
#     for index, combination in enumerate(combinations):
#         for id in range(len(atomtypes)):
#             if id in combination:
#                 atomtypes[id] = atom1
#             else:
#                 atomtypes[id] = atom2
#         # print(index, atomtypes)
#         list_atoms.append(atomtypes[:])

#     for file in files:
#         pfiles = []
#         X = []
#         list_coords = []
#         coulomb = []
#         base_name = os.path.basename(file)
#         base_name = os.path.splitext(base_name)[0]
#         natoms, atomtypes, coords, energy = tools.xyzRead(file)
#         for index, combination in enumerate(combinations):
#             for id in range(len(atomtypes)):
#                 if id in combination:
#                     atomtypes[id] = atom1
#                 else:
#                     atomtypes[id] = atom2
#             pfile = base_name+"_P"+str(index)
#             pfiles.append(pfile[:])
#             list_coords.append(coords[:])
#             coulomb.append(tools.eigenCoulomb(natoms, list_atoms[index], coords))

#         coulomb = np.array(coulomb)
#         coulomb = StandardScaler().fit_transform(coulomb)
#         sel_samples = rep.get_representatives(params["MOD2"], coulomb, pfiles)

#         for idx in sel_samples:
#             # print(idx, list_atoms[idx], list_coords[idx][0])
#             tools.generateXYZ(list_atoms[idx], list_coords[idx], 0.0, pfiles[idx], outfolder)

#     # net.integrity_test(tmpfolder)


# # VERSION 03
# def gen_permutations(params: dict, files: str) -> None:
#     """
#     Generate all permutations of the alloys
#     """
#     atom1 = params["ELEM1"]
#     atom2 = params["ELEM2"]
#     natoms = params["NUMELEM1"]+params["NUMELEM2"]
#     outfolder = params["MOD2"]["TMP_FOLDER"]+"/unfiltered/"
#     maxgen_frame = params["MOD2"]["MAX_GEN_PER_FRAME"]
#     pre_selection = params["MOD2"]["PRE_SELECTION"]
#         # "RUN_MOD_ZERO": false,

#     combinations = list(itertools.combinations(range(natoms), params["NUMELEM1"]))

#     for file in tqdm.tqdm(files):
#         pfiles = []
#         X = []
#         list_coords = []
#         list_atoms = []
#         coulomb = []
#         base_name = os.path.basename(file)
#         base_name = os.path.splitext(base_name)[0]
#         natoms, atomtypes, coords, energy = tools.xyzRead(file)
#         random.shuffle(combinations)



#         # if pre_selection:
#         #     prot_atoms = []
#         #     for index, combination in enumerate(combinations):
#         #         for id in range(len(atomtypes)):
#         #             if id in combination:
#         #                 atomtypes[id] = atom1
#         #             else:
#         #                 atomtypes[id] = atom2
#         #         # print(index, atomtypes)
#         #         prot_atoms.append(atomtypes[:])


#         for index, combination in enumerate(combinations):
#             if index < maxgen_frame:
#                 for id in range(len(atomtypes)):
#                     if id in combination:
#                         atomtypes[id] = atom1
#                     else:
#                         atomtypes[id] = atom2
#                 pfile = base_name+"_P"+str(index)
#                 if pre_selection:
#                     pfiles.append(pfile[:])
#                     list_coords.append(coords[:])
#                     list_atoms.append()
#                     coulomb.append(tools.eigenCoulomb(natoms, list_atoms[index], coords))
#                 else:
#                     tools.generateXYZ(atomtypes, coords, 0.0, pfile, outfolder)

#         if pre_selection:
#             coulomb = np.array(coulomb)
#             coulomb = StandardScaler().fit_transform(coulomb)
#             sel_samples = rep.get_representatives(params["MOD2"], coulomb, pfiles)

#             for idx in sel_samples:
#                 # print(idx, list_atoms[idx], list_coords[idx][0])
#                 tools.generateXYZ(list_atoms[idx], list_coords[idx], 0.0, pfiles[idx], outfolder)

#     # net.integrity_test(tmpfolder)

# VERSION 04
def gen_permutations(params: dict, files: list) -> None:
    """
    Generate all permutations of the alloys
    """
    atom1 = params["ELEM1"]
    atom2 = params["ELEM2"]
    natoms = params["NUMELEM1"] + params["NUMELEM2"]
    outfolder = params["MOD2"]["TMP_FOLDER"] + "/unfiltered"
    maxgen_frame = params["MOD2"]["MAX_GEN_PER_FRAME"]
    pre_selection = params["MOD2"]["PRE_SELECTION"]

    combinations = list(itertools.combinations(range(natoms), params["NUMELEM1"]))
    ids_combinations = list(range(len(combinations)))

    for file in tqdm.tqdm(files, desc='Processing Frames', position=0):
        pfiles = []
        X = []
        list_coords = []
        list_atoms = []
        coulomb = []
        base_name = os.path.basename(file)
        base_name = os.path.splitext(base_name)[0]
        natoms, atomtypes, coords, energy = tools.xyzRead(file)
        random.shuffle(ids_combinations)

        for index, id_combination in enumerate(ids_combinations):
            if index < maxgen_frame:
                for id_atom in range(len(atomtypes)):
                    if id_atom in combinations[id_combination]:
                        atomtypes[id_atom] = atom1
                    else:
                        atomtypes[id_atom] = atom2
                pfile = base_name + "_P" + str(index)
                if pre_selection:
                    pfiles.append(pfile)
                    list_coords.append(coords.copy())  # Используем копию
                    list_atoms.append(atomtypes.copy())  # Используем копию
                    coulomb.append(tools.eigenCoulomb(natoms, list_atoms[index], coords))
                else:
                    tools.generateXYZ(atomtypes, coords, 0.0, pfile, outfolder)

        if pre_selection:
            coulomb = np.array(coulomb)
            coulomb = StandardScaler().fit_transform(coulomb)
            energies = tools.getEnergies(files)

            # 🛠 Исправленный вызов get_representatives
            sel_samples = rep.get_representatives(params["MOD2"], coulomb, energies, pfiles)

            # 🚀 Исправление: проверяем, является ли sel_samples кортежем
            if isinstance(sel_samples, tuple) and len(sel_samples) > 0:
                sel_samples = sel_samples[0]  # Берем только индексы
            elif isinstance(sel_samples, np.ndarray):
                sel_samples = sel_samples.tolist()
            elif not isinstance(sel_samples, list):
                raise TypeError(f"Unexpected type for sel_samples: {type(sel_samples)}")

            # 🛠 Дополнительная проверка индексов
            valid_indices = [int(idx) for idx in sel_samples if 0 <= int(idx) < len(pfiles)]

            # 📌 Отладочный вывод
            print(f"Selected valid indices: {valid_indices}")

            # 🔥 Сохраняем только валидные структуры
            for idx in valid_indices:
                tools.generateXYZ(list_atoms[idx], list_coords[idx], 0.0, pfiles[idx], outfolder)