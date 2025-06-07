#!/usr/bin/python

#   Copyright 2023 Raphael Bühler, Max Schütz, Karla F. Andriani,
#   João Paulo A. de Mendonça, Vivianne K. Ocampo-Restrepo, Marcos G. Quiles 
#   Christian Gemel, Juarez L. F. Da Silva, Roland A. Fischer
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import os
# import time
import sys
import random
# import matplotlib.cm as cm
import json
import glob
import numpy as np
import tqdm
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import shutil
# Local packages
# import core.abcluster as abcluster
import core.connectivity as net
import core.representatives as rep
import core.permute as perm
import core.tools as tools
import core.frames as gen
import core.complexes as complexes
from tqdm import tqdm


def start_message() -> bool:
	"""
	Open info about the program parameters
	"""

	print("""\t~?JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ?!     
	P&#P55555555555555555555555555555555555B&#:    
	P&P                                    ?&B:    
	P&P            ^YPPPPPPPP5!            J&B:    
	P&P           7#&5??????J#&J           J&B:    
	P&P          J&#7        ~B&5:         J&B:    
	P&P   ^!~~~~5&B~          :P&G!~~~~^   J&B:    
	P&P   JBBBBB##J            !#&BBBBBP.  J&B:    
	P&P   ......7##J          7##Y......   J&B:    
	P&P          ~B&5:      .J##?          J&B:    
	P&P           :P&BGGGGGGG&#!           J&B:    
	P&P            .!7777777?B#5.          J&B:    
	P&P                      ^G&G^         J&B:    
	P&P                       .?!.         ?&B:    
	P&B?????????????^ .^^^^. .?????????????P&#:    
	7PPPPPPPPPPPPPPG! .~~~~: :PPPPPPPPPPPPPPPY.   
	""")

	print("\n	QTNano Cluster Assembler \n\n")


	if len(sys.argv) != 2:
		print("\tUsage: \n")
		print(f"\t\t$ python {sys.argv[0]} parameters.json\n")
		print(f"\t\t\tparameters.json - a json file that contains all the simulation parameters")
		return False	
	return True

def load_parameters() -> dict:
	"""
	Read the simulation parameters from a json file
	"""
	try:
		with open(sys.argv[1]) as f:
			j = json.load(f)
	except FileNotFoundError:
		print("[ERROR] input json file not found")
		exit()
	except json.decoder.JSONDecodeError:
		print("[ERROR] input json format is invalid")
		exit()
	except Exception:
		print("[ERROR] an unexpected error happened trying to read the input")
		exit()

	return j




def M0_Selection(params: dict, adhoc: bool) -> None:
	print("\n\t\tModule 0: Representative Sample Selection:")

	if adhoc:
		inputfolder = params["INPUT_FOLDER"]
		outfolder = params["OUTPUT_FOLDER"]
	else:
		inputfolder = params["TMP_FOLDER"] + "/filtered"
		outfolder = params["TMP_FOLDER"] + "/selected"

	if os.path.exists(outfolder):
		shutil.rmtree(outfolder)  # Удаляет папку и ее содержимое
		os.makedirs(outfolder, exist_ok=True)  # Создает новую пустую папку
		print(f"\n\t\t\tLoading files from: {inputfolder}")

	files = sorted(glob.glob(f"{inputfolder}/*.xyz"))
	if not files:
		print(f"\n\t\t\tFail - Folder {inputfolder} is empty\n\n")
		exit()

	print(f"\n\t\t\tSample Selection - Files Loaded {len(files)}")
	if len(files) > 1:
		coulomb, energies = tools.getCoulombEig(files)
		sel_samples = rep.get_representatives(params, coulomb, energies, files)
	else:
		sel_samples = [0]

	# Проверка и исправление sel_samples
	if isinstance(sel_samples, (list, np.ndarray)):
		sel_samples = [int(sample[0]) if isinstance(sample, (list, np.ndarray)) else int(sample) for sample in sel_samples]

	# Копирование выбранных файлов
	for id in sel_samples:
		infile = files[id]
		outfile = f"{outfolder}/{os.path.basename(infile)}"
		shutil.copy(infile, outfile)

	inFiles = len(glob.glob(f"{outfolder}/*.xyz"))
	print(f"\n\t\t\tTotal of Selected Samples: {inFiles}")

	print(f"\n\t\t\tSelected Files are available at {outfolder}")

	print("\n\t\tEnd of module 0 - Representative Selection")


def M1_frame_family(params: dict) -> None:
	print("\n\nModule 1: Frame Family:")

	tmpfolder = params["MOD1"]["TMP_FOLDER"]
	outfolder = params["MOD1"]["OUTPUT_FOLDER"]
	threshold = params["MOD1"]["INTEGRITY"]

	# Очистка папок
	shutil.rmtree(tmpfolder, ignore_errors=True)
	shutil.rmtree(outfolder, ignore_errors=True)
	os.makedirs(tmpfolder + "/unfiltered", exist_ok=True)
	os.makedirs(tmpfolder + "/filtered", exist_ok=True)
	os.makedirs(tmpfolder + "/selected", exist_ok=True)
	os.makedirs(outfolder, exist_ok=True)

	print(f"\t\tGenerating frames in: {tmpfolder}/unfiltered")

	# Параметры генерации
	num_structures = params["MOD1"]["NUMGEN"]  # Кол-во кластеров
	atom = params["ELEM1"]  # Тип атома (Fe, S и т. д.)
	n_atoms = params["NUMELEM1"] + params["NUMELEM2"]  # Число атомов в кластере
	shape = params["MOD1"]["SHAPE"]  # "CUBE" или "SPHERE"
	factor = params["MOD1"]["FACTOR"]  # Коэффициент размера
	gamma = params["MOD1"]["GAMMA"]  # Допуск по расстояниям

	# Генерация структур
	gen.genSamples(num_structures, shape, atom, n_atoms, factor, gamma, tmpfolder+"/unfiltered")

	print(f"\t\tRunning integrity test on generated frames...")

	net.integrity_test(tmpfolder = tmpfolder, threshold = threshold) 

	print(f"\t\tFrames are available at {outfolder}")
print("\nEnd of module 01 - Frame Family")
	
	
def M2_core_family(params: dict) -> None:
	print("\n\nModule 2: Core Family:")
	inputfolder = params["MOD2"]["INPUT_FOLDER"]
	tmpfolder = params["MOD2"]["TMP_FOLDER"]
	outfolder = params["MOD2"]["OUTPUT_FOLDER"]
	threshold = params["MOD2"]["INTEGRITY"]
	shutil.rmtree(tmpfolder, ignore_errors=True)
	shutil.rmtree(outfolder, ignore_errors=True)
	os.makedirs(tmpfolder + "/unfiltered", exist_ok=True)
	os.makedirs(tmpfolder + "/filtered", exist_ok=True)
	os.makedirs(tmpfolder + "/selected", exist_ok=True)
	os.makedirs(outfolder, exist_ok=True)
	print(f"\t\tGenerating Cores from: {inputfolder}")
	perm.gen_permutations(params, glob.glob(f"{inputfolder}/*.xyz"))
	net.integrity_test(tmpfolder, threshold)
	print(f"\t\tCores are available at {outfolder}")
	print("\nEnd of module 02 - Core Family")

def M3_add_ligants(params: dict) -> None:
	print("\n\nModule 3: Complexes Generation:")

	# Получение параметров
	coresfolder = params["MOD3"]["CORES_FOLDER"]
	ligandsfolder = params["MOD3"]["LIGANDS_FOLDER"]
	tmpfolder = params["MOD3"]["TMP_FOLDER"]
	outfolder = params["MOD3"]["OUTPUT_FOLDER"]
	threshold = params["MOD3"]["INTEGRITY"]
	lig_distribution = params["MOD3"]["LIGANDS_DISTRIBUTION"]
	lig_orientation = params["MOD3"]["LIGANDS_ORIENTATION"]
	numsim = params["MOD3"]["N_SAMPLES"]
	ncores = params["MOD3"]["N_CORES"]
	ligdist = params["MOD3"]["LIGANDS_DISTANCE"]
	deformation = params["MOD3"]["DEFORMATION"]

	# Безопасное создание директорий
	unfiltered_path = os.path.abspath(os.path.join(tmpfolder, "unfiltered"))  # Абсолютный путь
	filtered_path = os.path.abspath(os.path.join(tmpfolder, "filtered"))
	selected_path = os.path.abspath(os.path.join(tmpfolder, "selected"))
	os.makedirs(unfiltered_path, exist_ok=True)
	os.makedirs(filtered_path, exist_ok=True)
	os.makedirs(selected_path, exist_ok=True)
	os.makedirs(outfolder, exist_ok=True)

	# Загрузка ядер
	cores = sorted(glob.glob(os.path.join(coresfolder, "*.xyz")))
	if len(cores) < ncores:
		ncores = len(cores)
	cores = cores[:ncores]

	# Загрузка лигандов
	ligands = sorted(glob.glob(os.path.join(ligandsfolder, "*.xyz")))

	# Проверка доступности файлов
	if not cores:
		print(f"\n\t\tFail - No cores found in {coresfolder}\n")
		return
	if not ligands:
		print(f"\n\t\tFail - No ligands found in {ligandsfolder}\n")
		return

	print(f"\n\t\tLoaded {len(cores)} cores and {len(ligands)} ligands.")

	# Основной цикл обработки
	for core in tqdm(cores, desc='Processing Cores'):
		try:
			base_name = os.path.splitext(os.path.basename(core))[0]
			_, core_atomtypes, core_coords, _ = tools.xyzRead(core)

			for sim in range(numsim):
				complex_atoms = core_atomtypes.copy()
				complex_coords = core_coords.copy()

				# Генерация позиций лигандов
				for lig_idx, numlig in enumerate(lig_distribution):
					if numlig <= 0:
						continue

					lig_path = ligands[lig_idx % len(ligands)]
					_, lig_atomtypes, lig_coords, _ = tools.xyzRead(lig_path)
					
					# Позиционирование нескольких лигандов
					for _ in range(numlig):
						shift = np.random.rand(3) * ligdist
						rotated_ligand = complexes.rotate_atoms(lig_coords)
						placed_ligand = rotated_ligand + shift

						complex_atoms.extend(lig_atomtypes)
						complex_coords = np.concatenate((complex_coords, placed_ligand))

				# Корректный путь без двойного расширения
				pfile = os.path.join(unfiltered_path, f"{base_name}_C{sim}.xyz")
				tools.generateXYZ(complex_atoms, complex_coords, 0.0, pfile, unfiltered_path)

		except Exception as e:
			print(f"Error processing {core}: {str(e)}")
			continue

	# Проверка целостности
	generated_files = glob.glob(os.path.join(unfiltered_path, "*.xyz"))
	print(f"\t\t\tTotal Generated Complexes: {len(generated_files)}")

	print("\n\t\tTesting the integrity of the complexes")
	net.integrity_test_complexes(tmpfolder, threshold)
	inFiles = glob.glob(os.path.join(filtered_path, "*.xyz"))
	print(f"\t\t\tTotal of Complexes After The Integrity Test: {len(inFiles)}")
	runZero= params["MOD3"]["RUN_MOD_ZERO"]
	if runZero:
		print("\n\t\tSelecting representative complexes (k-means)")
		M0_Selection(params["MOD3"], False)
		# files=sorted(glob.glob(f"{tmpfolder}/filtered/*.xyz"))
		# coulomb = tools.getCoulombEig(files)
		# sel_samples = rep.get_representatives(params["MOD3"], coulomb, files) # /filtered to /selected
		filfold = tmpfolder+"/selected"
		# for id in sel_samples:
		# 	infile = files[id]
		# 	outfile = f"{filfold}/{os.path.basename(infile)}"
		# 	os.system(f"cp {infile} {outfile}") #/selected to /output_folder(params)
		# inFiles = len(glob.glob(f"{tmpfolder}/selected/*.xyz"))
		# print(f"\t\t\tTotal of Selected Complexes: {inFiles}")
	else:
		filfold = tmpfolder+"/filtered"
		for file in glob.glob(os.path.join(filfold, "*.xyz")):
			shutil.copy(file, outfolder)
	
	print(f"\n\t\tComplexes are available at {outfolder}")
	print("\nEnd of module 03 - Complexes Generation")

def main() -> None:
	"""
	Main routines
	"""
	if start_message() == False:
		exit("---DONE---")
	
	params = load_parameters()

	# HPC = params["RUNHPC"]

	if 0 in params["MODULES"]:
		M0_Selection(params["MOD0"], True)

	if 1 in params["MODULES"]:
		M1_frame_family(params)

	if 2 in params["MODULES"]:
		M2_core_family(params)

	if 3 in params["MODULES"]:
		M3_add_ligants(params)

	print("\n\nE N D !")


if __name__ == "__main__":
	main()

