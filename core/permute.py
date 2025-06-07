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

def load_atomic_numbers(filename="mol.txt"):
    """Загружает атомные номера из файла mol.txt"""
    atomic_numbers = {}
    with open(filename, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                element, number = parts[0], int(parts[1])
                atomic_numbers[element] = number
    return atomic_numbers

def gen_permutations(params: dict, files: list) -> None:
    """Generate all permutations of the alloys"""
    atom1 = params["ELEM1"]
    atom2 = params["ELEM2"]
    natoms = params["NUMELEM1"] + params["NUMELEM2"]
    outfolder = params["MOD2"]["TMP_FOLDER"] + "/unfiltered"
    maxgen_frame = params["MOD2"]["MAX_GEN_PER_FRAME"]
    pre_selection = params["MOD2"]["PRE_SELECTION"]

    # Загружаем атомные номера
    atomic_numbers = load_atomic_numbers()

    combinations = list(itertools.combinations(range(natoms), params["NUMELEM1"]))
    ids_combinations = list(range(len(combinations)))

    for file in tqdm.tqdm(files, desc="Processing Frames", position=0):
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

                # Заменяем атомные символы на номера
                atomtypes_numeric = [atomic_numbers[atom] for atom in atomtypes]

                pfile = base_name + "_P" + str(index)
                if pre_selection:
                    pfiles.append(pfile)
                    list_coords.append(coords.copy())
                    list_atoms.append(atomtypes_numeric)  # Сохраняем как номера
                    coulomb.append(tools.eigenCoulomb(natoms, atomtypes, coords))
                else:
                    tools.generateXYZ(atomtypes, coords, 0.0, pfile, outfolder)

        if pre_selection:
            coulomb = np.array(coulomb)
            coulomb = StandardScaler().fit_transform(coulomb)
            energies = tools.getEnergies(files)

            sel_samples = rep.get_representatives(params["MOD2"], coulomb, energies, pfiles)

            if isinstance(sel_samples, tuple) and len(sel_samples) > 0:
                sel_samples = sel_samples[0]
            elif isinstance(sel_samples, np.ndarray):
                sel_samples = sel_samples.tolist()
            elif not isinstance(sel_samples, list):
                raise TypeError(f"Unexpected type for sel_samples: {type(sel_samples)}")

            valid_indices = [int(idx) for idx in sel_samples if 0 <= int(idx) < len(pfiles)]

            for idx in valid_indices:
                atomic_symbols = {28: "Ni", 31: "Ga"}  # Или загрузи из mol.txt, но в обратном порядке
                list_atoms_symbols = [[atomic_symbols[num] for num in atoms] for atoms in list_atoms]
                tools.generateXYZ(list_atoms_symbols[idx], list_coords[idx], 0.0, pfiles[idx], outfolder)




