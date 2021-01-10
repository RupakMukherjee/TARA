#include <cufft.h>
#include "cufftUtils.h"
 
// Declared extern "C" to disable C++ name mangling

extern "C" void createCUFFTPlan2D(void *plan, int nx, int ny, int planType, void *stream)
{
    cufftHandle *cPlan;

    cPlan = (cufftHandle *)plan;
    CHECK_CUFFT(cufftPlan2d(cPlan, nx, ny, (cufftType)planType));
    CHECK_CUFFT(cufftSetStream(*cPlan,  (cudaStream_t)stream));
}

extern "C" void executeCUFFT2D(void *plan, void *iData, void *oData, int planType)
{
    cufftHandle *cPlan;

    cPlan = ((cufftHandle *)plan);
    switch (planType)
    {
        case CUFFT_D2Z: CHECK_CUFFT(cufftExecD2Z(*cPlan, (cufftDoubleReal *)iData, (cufftDoubleComplex *)oData));
                        break;
        case CUFFT_Z2D: CHECK_CUFFT(cufftExecZ2D(*cPlan, (cufftDoubleComplex *)iData, (cufftDoubleReal *)oData));
                        break;
    }
}

extern "C" void destroyCUFFTPlan2D(void *plan)
{
    cufftHandle *cPlan;

    cPlan = (cufftHandle *)plan;
    CHECK_CUFFT(cufftDestroy(*cPlan));
}  
