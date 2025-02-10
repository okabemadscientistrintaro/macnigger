$ conda env create -f environment.yml

$ conda activate py396
$ python main.py parameters.json



#### Main Settings
| Parameter   | Type           | Description                                                                 | Default Value |
|-------------|----------------|-----------------------------------------------------------------------------|---------------|
| `MODULES`   | list of integers | Select which module (or modules) will be run                                | `[1, 2, 3]`   |
| `ELEM1`     | string         | The first species of the clusters                                           | `"Cu"`        |
| `ELEM2`     | string         | The second species of the clusters                                          | `"Zn"`        |
| `NUMELEM1`  | integer        | Number of atoms of the first species                                        | `3`           |
| `NUMELEM2`  | integer        | Number of atoms of the second species                                       | `4`           |

---

#### Common Parameters (Modules 0, 1, 2, and/or 3)
| Parameter       | Type           | Description                                                                 | Default Value       |
|-----------------|----------------|-----------------------------------------------------------------------------|---------------------|
| `KMEANS`        | list of integers | Define the number of seeds (k) for k-means. A single value sets the exact number; a range (e.g., `[2, 50, 1]`) triggers Silhouette analysis to find the optimal k. | `[10]` or `[2, 50, 1]` |
| `RUN_KMEANS`    | integer        | Number of k-means runs                                                      | `10`                |
| `MAXSAMPLES`    | integer        | Number of samples selected per cluster                                      | `1`                 |
| `INPUT_FOLDER`  | string         | Input data folder                                                           | `./input_data`      |
| `OUTPUT_FOLDER` | string         | Output data folder                                                          | `./output_data`     |
| `TMP_FOLDER`    | string         | Temporary data folder for intermediate steps                                | `./var_data`        |
| `RUN_MOD_ZERO`  | boolean        | Whether to run Module 0 (k-means clustering for representative structures)  | `false`             |

---

#### Module 1 - Frame Family
| Parameter        | Type    | Description                                                                 | Default Value |
|------------------|---------|-----------------------------------------------------------------------------|---------------|
| `INTEGRITY`      | float   | Tolerance for covalent radius to form a structure (e.g., 1 = covalent radius limit) | `1.5`         |
| `NUMGEN`         | integer | Number of generated frames                                                  | `100`         |
| `RADIUS_FACTOR`  | float   | Controls the size of the box/sphere                                         | `0.5`         |
| `GAMMA`          | float   | Tolerance for species covalent radius                                       | `0.2`         |

---

#### Module 2 - Core Family
| Parameter              | Type    | Description                                                                 | Default Value |
|------------------------|---------|-----------------------------------------------------------------------------|---------------|
| `MAX_GEN_PER_FRAME`    | integer | Number of cores to generate for each input frame                            | `100`         |
| `INTEGRITY`            | float   | Tolerance for covalent radius to form a structure                           | `1.5`         |

---

#### Module 3 - Complexes Generation
| Parameter               | Type           | Description                                                                 | Default Value       |
|-------------------------|----------------|-----------------------------------------------------------------------------|---------------------|
| `DEFORMATION`           | boolean        | If `true`, sites are adjusted to the surface of the core                    | `true`              |
| `INTEGRITY`             | float          | Tolerance for covalent radius to form complexes                             | `1.25`              |
| `LIGANDS_DISTRIBUTION`  | list of integers | Defines ligands and their quantities. The first number is ligand 1, second is ligand 2, etc. | `[3, 2]`            |
| `LIGANDS_ORIENTATION`   | list of lists  | Ligand orientation: `[-1]` for random, `[0, 1]` for diatomic molecules, `[1, 2, 3]` for specific atoms | `[[-1], [0, 1, 2]]` |
| `N_SAMPLES`             | integer        | Number of complexes generated for each input core                           | `100`               |
| `N_CORES`               | integer        | Number of cores used in the simulation (limited by files in `CORES_FOLDER`) | `3`                 |
| `LIGANDS_DISTANCE`      | float          | Distance between the ligand and the core                                    | `2.25`              |
| `CORES_FOLDER`          | string         | Folder containing core (XYZ) files                                          | `./cores`           |
| `LIGANDS_FOLDER`        | string         | Folder containing ligand (XYZ) files                                        | `./ligands`         |




Modules
Module 0 - MOD0
This module is a k-means clustering tool to select the representatives unary cores, binary cores, and complexes. This module uses the eigenvalues of the Coulomb matrix of the provided structures to select a set of representative ones.

Module 1 - MOD1.
In this module, unary cores are generated, processed, analyzed, and clustered. We included a series of user options that allow combining external sourcing and internally generating the unary cores as needed by the case of study. The cores are checked according to a integrity filter (covalent radius) and clustered by k-means if the RUN_MOD_ZERO flag is activated.

Module 2 - MOD2.
Here, permutations (i.e., all possible or user-specified) are performed in the atomic species of each of the selected unary cores. The binary cores are subsequently checked for integrity and clustered by k-means if the RUN_MOD_ZERO flag is activated.

Module 3 - MOD3.
Module responsable to distribute the ligands around the selected binary cores. The ligands XYZ are imported from ligands/ folder. Several ligands might be used in the process. Their orientation regarding the core is controled by the user, which can select a random or specific orientation (provided by a set of base atoms). An overlap filter iteratively excludes structures with overlapping atoms (integrity). The default threshold is set to t=1.25, and can be manually tuned. The filtered complexes are again clustered by k-means if the RUN_MOD_ZERO flag is activated. Finally, the selected structures are the representatives within the living library family of structures, and can be externally optimized using methods such as Density Functional Theory.
