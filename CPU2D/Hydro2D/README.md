All the codes in this folder are serial codes. There are no  parallel codes in this file.


vorticity_wo_tracer.f95: incompressible single fluid with PS scheme (serial code + time solver AB / RK-4 + NO tracer particles) :: This is the most fundamental code!

Two_demensional_Navier_Stokes_Vorticity_StreamFunction.f95: What is the difference between vorticity_wo_tracer.f95 code?

2d_two_fluid_incomp.f95: incompressible two fluid problem with PS scheme (serial code + time solver AB)

hybrid_fluid.f95: incompressible single fluid with PS + WENO scheme (serial code + time solver TVD-RK = Total Variation Diminished RK)

incomp_hydro_2d_mpi.f95: incompressible single fluid with PS scheme (MPI code without I/O parallelization)

incomp_hydro_2d_mpi_alloc.f95: incompressible single fluid with PS scheme (MPI code WITH I/O parallelization)

vorticity.f95: incompressible single fluid with PS scheme (serial code + time solver AB + tracer particles)

compressible_hydro2d.f95: compressible single fluid with PS scheme (openmp code + time solver AB / Taylor/ Predoctor-Corrector / RK-4)
