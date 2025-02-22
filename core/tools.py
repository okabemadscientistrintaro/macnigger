from sklearn.preprocessing import StandardScaler
import numpy as np
import numpy.typing as npt

import tqdm
import glob
import math
import ast
import os

def xyzRead(fname: str):
    with open(fname, "r") as fin:
        line1 = fin.readline().split()
        natoms = int(line1[0])
        comments = fin.readline().strip()
        
        # Попытка считать энергию
        energy = None
        if comments:
            try:
                energy = float(comments.split()[-1])  # Делаем float, а не int
            except ValueError:
                energy = None

        coords = np.zeros([natoms, 3], dtype="float64")
        atomtypes = []
        
        for x in coords:
            line = fin.readline().split()
            atomtypes.append(line[0])
            x[:] = list(map(float, line[1:4]))
    
    return natoms, atomtypes, coords, energy

def getEnergies(files):
    return [xyzRead(file)[3] for file in files]  # Читаем 4-й элемент (energy)


# def setAtomNull(params: dict, folder: str):
#     files=glob.glob(folder+"/*.xyz")
#     for f in files:
#         with open(f, 'r') as file:
#             lines = file.readlines()
#         file.close()
#         lines = [w.replace(params["ELEM1"], "At") for w in lines]
#         lines = [w.replace(params["ELEM2"], "At") for w in lines]
#         newfile = open(f, 'w')
#         for l in lines:
#             newfile.write(l) 
#         newfile.close()

def replaceAtomSymbols(params: dict, fname: str, permut: tuple):
    with open(fname, 'r') as file:
        line = file.readline()
        modified = [line]
        line = file.readline()
        modified.append(line)
        lines = file.readlines()
        for at, lin in enumerate(lines):
            if at in permut:
                modified.append(lin.replace('At', params["ELEM1"]))
            else:
                modified.append(lin.replace('At', params["ELEM2"]))

    with open(fname, 'w') as file:
        file.writelines(modified)
    file.close()


def getCharge(element):
    f = open("mol.txt")
    atomicnum = [line.split()[1] for line in f if line.split()[0] == element]
    f.close()
    return int(atomicnum[0])

def coulombMatrix(natoms, atomtypes, coords):
    i=0 ; j=0    
    colM = np.zeros((natoms,natoms))
    chargearray = np.zeros((natoms,1))
    charge = [getCharge(symbol)  for symbol in atomtypes]
    for i in range(0,natoms):
        colM[i,i]=0.5*charge[i]**2.4   # Diagonal term described by Potential energy of isolated atom
        for j in range(i+1,natoms):
            dist= np.linalg.norm(coords[i,:] - coords[j,:])   
            colM[j,i] = charge[i]*charge[j]/dist   #Pair-wise repulsion 
            colM[i,j] = colM[j,i]
    return colM
 
def eigenCoulomb(natoms, atomtypes, coords):
    sCoulomb = coulombMatrix(natoms, atomtypes, coords)
    sCoulomb = sCoulomb.astype(int)
    eigValues = -np.sort(-np.linalg.eigvals(sCoulomb))
    if np.any(np.iscomplex(eigValues)) == False:
        return eigValues#[0:natoms]
    else:
        # print('\nWARNING: complex (coulomb matrix) engenvalues for ' + fname + ' employing np.linalg.eigvals function.\n') 
        eigValues = -np.sort(-np.linalg.eigvalsh(sCoulomb))
        if np.any(np.iscomplex(eigValues)):
            # print('\nWARNING: complex (coulomb matrix) engenvalues for ' + fname + ' employing np.linalg.eigvalsh function.\n WARNING: only the real part will be returned.\n')
            return eigValues.real#[0:natoms]
        else :
            return eigValues#[0:atoms]
        
def generateXYZ(list_atoms, list_coords, energy, pfile, outfolder):
    """
    Generate an XYZ file from atomic data and save it to the specified folder.

    Parameters:
        list_atoms (list): List of atomic symbols (e.g., ['Fe', 'S', 'Fe']).
        list_coords (list): List of atomic coordinates (e.g., [[x1, y1, z1], [x2, y2, z2], ...]).
        energy (float): Energy value to include in the file.
        pfile (str): Base name for the output file (e.g., "structure_1").
        outfolder (str): Directory where the file will be saved.
    """
    # Ensure the output folder exists
    os.makedirs(outfolder, exist_ok=True)

    # Construct the full file path using os.path.join
    fname = os.path.join(outfolder, f"{pfile}.xyz")

    # Write the XYZ file
    with open(fname, 'w') as file:
        # Line 1: Number of atoms
        file.write(f"{len(list_atoms)}\n")

        # Line 2: Comment line (energy)
        file.write(f"Energy = {energy}\n")

        # Subsequent lines: Atomic symbols and coordinates
        for idx in range(len(list_atoms)):
            file.write(f"{list_atoms[idx]} {' '.join(map(str, list_coords[idx]))}\n")

    # print(f"File saved: {fname}")

def getCoulombEig(files):
    coulomb = []
    energies = []
    
    # Сначала собираем все собственные значения и энергии
    for file in files:
        natoms, atomtypes, coords, energy = xyzRead(file)
        eig_values = eigenCoulomb(natoms, atomtypes, coords)
        coulomb.append(eig_values)
        energies.append(energy)
    
    # Находим максимальную длину собственных значений
    max_length = max(len(eig) for eig in coulomb)
    
    # Дополняем собственные значения нулями до одинаковой длины
    coulomb_padded = [np.pad(eig, (0, max_length - len(eig)), mode='constant') for eig in coulomb]
    
    # Преобразуем в массив NumPy
    coulomb = np.array(coulomb_padded)
    energies = np.array(energies)
    
    # Масштабируем данные (если это необходимо)
    coulomb = StandardScaler().fit_transform(coulomb)
    
    return coulomb, energies
