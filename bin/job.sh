#export I_MPI_FABRICS=shm
#module load gcc/5.5
#source /opt/intel/parallel_studio_xe_2020.0.088/psxevars.sh

export OMP_NUM_THREADS=1
PATH_ABCLUSTER/isomer FNAME.inp > FNAME.out

# module purge
