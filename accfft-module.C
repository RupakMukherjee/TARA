#include <accfft_gpu.h>
#include <cuda_runtime_api.h>
#include <mpi.h>
#define CHECK_CUDA(call)                                                      \
{                                                                              \
     cudaError_t err;                                                          \
     if ( (err = (call)) != CUDA_SUCCESS)                                      \
     {   fprintf(stderr, "Got CUDA error %d at %s:%d\n", err, __FILE__,        \
                 __LINE__);                                                     \
         exit(1);                                                               \
     }                                                                          \
}
 
// Declared extern "C" to disable C++ name mangling

extern "C" void setGPUDevice(int gpuid)
{
    CHECK_CUDA(cudaSetDevice(gpuid));
}

extern "C" void accfftCreateComm(int nprocs, void **comm)
{
    int c_dims[2];
    MPI_Comm *communicator;
    MPI_Comm **comm_ptr;
//Creating a slab or 1D decomposition

    c_dims[0] = nprocs;
    c_dims[1] = 1;
    communicator = (MPI_Comm *) malloc (sizeof(MPI_Comm));
    accfft_create_comm(MPI_COMM_WORLD, c_dims, communicator);
    comm_ptr = (MPI_Comm **)comm;
    *comm_ptr = communicator;
}

extern "C" void accfftDestroyComm(void **comm)
{
    MPI_Comm **comm_ptr;

    comm_ptr = (MPI_Comm **)comm;
   
    MPI_Comm_free(*comm_ptr);
    free(*comm_ptr);
    *comm_ptr = 0;
}
    
extern "C" void accfftLocalSizeDFTR2CGPU(void *n, void *isize, void *istart, void *osize, void *ostart, void *comm, void *allocsize)
{
   int alloc_max = 0;
   int *allocsize_ptr = (int *)allocsize;
   MPI_Comm *comm_ptr = (MPI_Comm *)comm;
   
   alloc_max = accfft_local_size_dft_r2c_gpu((int *)n, (int *)isize, (int *)istart, (int *)osize, (int *)ostart,*comm_ptr); 
   *allocsize_ptr = alloc_max;
} 

extern "C" void accfftCreatePlan3DR2CGPU(void *n, void *idata, void *odata, void *comm, void **plan)
{
    accfft_plan_gpu **plan_ptr;
    accfft_plan_gpu *plan_gpu;
    MPI_Comm *comm_ptr = (MPI_Comm *)comm;

   

    plan_ptr = (accfft_plan_gpu **)plan;
  
    plan_gpu = accfft_plan_dft_3d_r2c_gpu((int *)n, (double *)idata, (double *)odata, *comm_ptr, ACCFFT_MEASURE);

   *plan_ptr = plan_gpu;
} 

extern "C" void accfftDestroyPlan3DR2CGPU(void *plan)
{
    accfft_plan_gpu *plan_ptr = (accfft_plan_gpu *)plan;
    
    accfft_destroy_plan_gpu(plan_ptr);
    accfft_cleanup_gpu();
}

extern "C" void accfftExecuteR2CGPU(void *plan, void *idata, void *odata, void *timer)
{
    accfft_execute_r2c_gpu((accfft_plan_gpu *)plan, (double *)idata, (Complex *)odata, (double *)timer);
    cudaDeviceSynchronize();
}

extern "C" void accfftExecuteC2RGPU(void *plan, void *idata, void *odata, void *timer)
{
    accfft_execute_c2r_gpu((accfft_plan_gpu *)plan, (Complex *)idata, (double *)odata, (double *)timer);
    cudaDeviceSynchronize();
}

