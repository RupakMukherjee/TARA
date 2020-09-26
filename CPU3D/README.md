mhd_3d.f95: 3D MHD solver using PS scheme + time solver AB + openMP parallel

mhd_3d_mpi.f95: 3D MHD solver using PS scheme + time solver AB + MPI parallel (Old version, contains bugs)

mhd_3d_mpi_alloc.f95: 3D MHD solver using PS scheme + time solver AB + MPI parallel WITH I/O parallelization

mhd_3d_mpi_new.f95: 3D MHD solver using PS scheme + time solver AB + MPI parallel without I/O parallelization

jobfile_mpi.sh: Jobfile to run the code in PBS based clusters

PSMHD3: A folder containing the same 3D MHD solver + with a MAKEFILE and subroutines splitted in different files (only a few subroutines).

PSMHD3_FINAL: A folder containing the same 3D MHD solver + with a MAKEFILE and subroutines splitted in different files (heavily splitted, but contains a bug).
