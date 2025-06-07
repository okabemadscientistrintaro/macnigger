import ase.neighborlist as ase_n
# from ase import NeighborList, neighbor_list
# from ase import geometry
# from ase.build import molecule
from ase.io import read
from scipy import sparse
import os
import glob
import shutil
import tqdm
import numpy as np
import tqdm
def integrity_test( tmpfolder: str, threshold: float):
        """
        Perform an integrity test on the generated structures.

        Parameters:
            tmpfolder (str): Path to the temporary folder containing the unfiltered files.
            threshold (float): Integrity threshold for filtering structures.
        """
        # Get all .xyz files in the unfiltered folder
        files = glob.glob(os.path.join(tmpfolder, "./unfiltered", "*.xyz"))
        folder_out = os.path.join(tmpfolder, "./filtered")

        # Create the output folder if it doesn't exist
        if not os.path.exists(folder_out):
            print(f"Creating folder: {folder_out}")
            os.makedirs(folder_out, exist_ok=True)
        print(f"DEBUG: Найдено {len(files)} файлов в {tmpfolder}/unfiltered")
        # Process each file
        for file in tqdm.tqdm(files, desc="Running integrity test"):
            try:
                mol = read(file)  # Read the molecule from the XYZ file
                cutOff = ase_n.natural_cutoffs(mol)  # Get natural cutoffs for each atom
                rmax = np.max(cutOff)  # Maximum cutoff distance
                i, j, d = ase_n.neighbor_list('ijd', mol, rmax * 2.0)  # Get neighbor list

                # Check integrity based on the threshold
                file_ok = True
                for index_i, index_j, dist in zip(i, j, d):
                    if dist * threshold < (cutOff[index_i] + cutOff[index_j]):
                        # print(f"DEBUG: Файл {file}  НЕ прошел тест")
                        file_ok = False
                        break

                # Save the file if it passes the integrity test
                if file_ok:
                    # print(f"DEBUG: Файл {file} прошел тест")
                    shutil.copy(file, os.path.join(folder_out, os.path.basename(file)))
            except Exception as e:
                print(f"Error processing file {file}: {e}")

                
                
# def integrity_test_complexes(folder: str, threshold: float):
#     # Исправление пути: ищем файлы в подпапке unfiltered
#     input_folder = os.path.join(folder, "unfiltered")
#     files = glob.glob(os.path.join(input_folder, "*.xyz"))
    
#     # Создание выходной папки
#     output_folder = os.path.join(folder, "filtered")
#     os.makedirs(output_folder, exist_ok=True)
    
#     # Проверка наличия файлов
#     if not files:
#         print(f"ERROR: No files found in {input_folder}")
#         return

#     print(f"Processing {len(files)} complexes...")
    
#     for file in tqdm.tqdm(files, desc="Integrity Check"):
#         try:
#             mol = read(file)
#             cutOff = ase_n.natural_cutoffs(mol)
#             rmax = np.max(cutOff)
            
#             # Получение списка соседей
#             i, j, d = ase_n.neighbor_list('ijd', mol, rmax * threshold)  # Исправлено!     
#             file_ok = True
#             for idx_i, idx_j, dist in zip(i, j, d):
#                 # Исправленное условие проверки
#                 if dist < (cutOff[idx_i] + cutOff[idx_j]) * threshold:
#                     file_ok = False
#                     break
                    
#             if file_ok:
#                 output_path = os.path.join(output_folder, os.path.basename(file))
#                 shutil.copy(file, output_path)
                
#         except Exception as e:
#             print(f"Error processing {file}: {str(e)}")
    
#     print(f"Saved {len(os.listdir(output_folder))} valid complexes")


# def integrity_test_complexes(folder: str, threshold: float):
#     input_dir = os.path.join(folder, "unfiltered")
#     output_dir = os.path.join(folder, "filtered")
    
#     # Проверка существования папок
#     if not os.path.exists(input_dir):
#         raise FileNotFoundError(f"Directory {input_dir} does not exist!")
#     os.makedirs(output_dir, exist_ok=True)
    
#     files = glob.glob(os.path.join(input_dir, "*.xyz"))
#     if not files:
#         print(f"No files found in {input_dir}")
#         return

#     passed = 0
#     for file in tqdm.tqdm(files, desc="Integrity Check"):
#         try:
#             mol = read(file)
#             cutOff = ase_n.natural_cutoffs(mol)
#             rmax = np.max(cutOff)
#             i, j, d = ase_n.neighbor_list('ijd', mol, rmax * threshold)
            
#             min_dist = np.min(d) if len(d) > 0 else 0
#             sum_radii = (cutOff[i] + cutOff[j]).min() if len(i) > 0 else 0
#             threshold_val = sum_radii * threshold
            
#             print(f"\nFile: {os.path.basename(file)}")
#             print(f"Min distance: {min_dist:.2f} Å")
#             print(f"Threshold: {threshold_val:.2f} Å")
            
#             file_ok = min_dist >= threshold_val
#             if file_ok:
#                 shutil.copy(file, os.path.join(output_dir, os.path.basename(file)))
#                 passed += 1
                
#         except Exception as e:
#             print(f"Error processing {file}: {str(e)}")
    
#     print(f"\nPassed {passed}/{len(files)} complexes")


def integrity_test_complexes(folder: str, threshold: float):
    input_dir = os.path.join(folder, "unfiltered")
    output_dir = os.path.join(folder, "filtered")
    
    os.makedirs(output_dir, exist_ok=True)
    
    files = glob.glob(os.path.join(input_dir, "*.xyz"))
    if not files:
        print(f"No files found in {input_dir}")
        return

    
    for file in tqdm.tqdm(files, desc="Integrity Check"):
        try:
            mol = read(file)
            cutOff = ase_n.natural_cutoffs(mol)
            i, j, d = ase_n.neighbor_list('ijd', mol, cutoff=np.max(cutOff)*2*threshold)
            
            # Преобразуем индексы в целые числа
            i = i.astype(int)
            j = j.astype(int)
            
            file_ok = True
            for idx_i, idx_j, dist in zip(i, j, d):
                if idx_i >= len(cutOff) or idx_j >= len(cutOff):
                    continue
                if dist < (cutOff[idx_i] + cutOff[idx_j]) * (1 / threshold):
                    file_ok = False
                    break
            if file_ok:
                shutil.copy(file, os.path.join(output_dir, os.path.basename(file)))
                passed += 1
                
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
    
    
