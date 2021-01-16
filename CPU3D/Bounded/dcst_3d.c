#include <fftw3.h>

int dcst3dx_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=*n0-1;i0>0;i0--) {
        for(int i1=0;i1<*n1;i1++) {
            for(int i2=0;i2<*n2;i2++) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[((i0-1)*(*n1)+i1)*(*n2)+i2] * (i0*0.5) / dn;
            }
        }
    }
    for(int i1=0;i1<*n1;i1++) {
        for(int i2=0;i2<*n2;i2++) {
            out3cs[i1*(*n2)+i2] = 0.0;
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_REDFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dy_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=*n1-1;i1>0;i1--) {
            for(int i2=0;i2<*n2;i2++) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[(i0*(*n1)+i1-1)*(*n2)+i2] * (i1*0.5) / dn;
            }
        }
    }
    for(int i0=0;i0<*n0;i0++) {
        for(int i2=0;i2<*n2;i2++) {
            out3cs[i0*(*n1)*(*n2)+i2] = 0.0;
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_RODFT00, FFTW_REDFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dz_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            for(int i2=*n2-1;i2>0;i2--) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[(i0*(*n1)+i1-1)*(*n2)+i2-1] * (i2*0.5) / dn;
            }
        }
    }
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            out3cs[(i0*(*n1)+i1)*(*n2)] = 0.0;
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dxy_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=*n0-1;i0>0;i0--) {
        for(int i1=*n1-1;i1>0;i1--) {
            for(int i2=0;i2<*n2;i2++) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[((i0-1)*(*n1)+i1-1)*(*n2)+i2] * (i0*i1*0.25) / dn;
            }
        }
    }
    for(int i1=0;i1<*n1;i1++) {
        for(int i2=0;i2<*n2;i2++) {
            out3cs[i1*(*n2)+i2] = 0.0;
        }
    }
    for(int i0=0;i0<*n0;i0++) {
        for(int i2=0;i2<*n2;i2++) {
            out3cs[i0*(*n1)*(*n2)+i2] = 0.0;
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_REDFT00, FFTW_REDFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dyz_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=*n1-1;i1>0;i1--) {
            for(int i2=*n2-1;i2>0;i2--) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[(i0*(*n1)+i1-1)*(*n2)+i2-1] * (i1*i2*0.25) / dn;
            }
        }
    }
    for(int i0=0;i0<*n0;i0++) {
        for(int i2=0;i2<*n2;i2++) {
            out3cs[i0*(*n1)*(*n2)+i2] = 0.0;
        }
    }
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            out3cs[(i0*(*n1)+i1)*(*n2)] = 0.0;
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_RODFT00, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dxz_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=*n0-1;i0>0;i0--) {
        for(int i1=0;i1<*n1;i1++) {
            for(int i2=*n2-1;i2>0;i2--) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[((i0-1)*(*n1)+i1)*(*n2)+i2-1] * (i0*i2*0.25) / dn;
            }
        }
    }
    for(int i1=0;i1<*n1;i1++) {
        for(int i2=0;i2<*n2;i2++) {
            out3cs[i1*(*n2)+i2] = 0.0;
        }
    }
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            out3cs[(i0*(*n1)+i1)*(*n2)] = 0.0;
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_REDFT00, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dxx_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            for(int i2=0;i2<*n2;i2++) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[(i0*(*n1)+i1)*(*n2)+i2] * (-(i0+1)*(i0+1)*0.25) / dn;
            }
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dyy_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            for(int i2=0;i2<*n2;i2++) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[(i0*(*n1)+i1)*(*n2)+i2] * (-(i1+1)*(i1+1)*0.25) / dn;
            }
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dzz_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            for(int i2=0;i2<*n2;i2++) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[(i0*(*n1)+i1)*(*n2)+i2] * (-(i2+1)*(i2+1)*0.25) / dn;
            }
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}

int dcst3dpoi_(int *n0, int *n1, int *n2, double *in3cs)
{
    double *out3cs;
    out3cs = (double*) fftw_malloc(sizeof(double) * (*n0) * (*n1) * (*n2));
    fftw_plan p3cs;
    p3cs = fftw_plan_r2r_3d(*n0, *n1, *n2, in3cs, out3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3cs);
    fftw_destroy_plan(p3cs);
    
    double dn = 8.0 * (*n0) * (*n1) * (*n2);
    for(int i0=0;i0<*n0;i0++) {
        for(int i1=0;i1<*n1;i1++) {
            for(int i2=0;i2<*n2;i2++) {
                out3cs[(i0*(*n1)+i1)*(*n2)+i2] = out3cs[(i0*(*n1)+i1)*(*n2)+i2] / ( (-(i0+1)*(i0+1)*0.25) + (-(i1+1)*(i1+1)*0.25) + (-(i2+1)*(i2+1)*0.25) ) / dn;
            }
        }
    }
    
    fftw_plan p3csi;
    p3csi = fftw_plan_r2r_3d(*n0, *n1, *n2, out3cs, in3cs, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
    
    fftw_execute(p3csi);
    fftw_destroy_plan(p3csi);
    fftw_free(out3cs);
    
    return 0;
}
