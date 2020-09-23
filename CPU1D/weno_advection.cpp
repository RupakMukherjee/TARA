/*
** 1-D FLUID solver using WENO(Weighted Essentially Non-Oscillatory) algorithm**
@file                weno_advection.cpp
@author              Sayan Adhikari <sayan.adhikari@fys.uio.no>
@source              Rupak Mukherjee <rupakmukherjee06@gmail.com>
@date                23.09.2020
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <stdlib.h>

using namespace std;


/* Define universal constants */
const double pi = 3.14159265359;

/* Define Simulation Parameters */
int N = 128;

/* Simulation Parameters */
double L = 2.0*pi;
double dx = L/(double)(N);

double time_min = 0.0;
double time_max = 1.5;
double dt = 0.0010;
/***************************/

void Write_data(FILE *file1, double x[], double u_bar[]);
void weno_reconstruction(double u_bar[], double RHS[]);

int main()
{


  /* House Keeping */
  system("rm -rf output");
  /*create an output folder*/
  system("mkdir output");

  int i;

  double gamma_1, gamma_2, gamma_3;
  double epslion;
  double gamma_1R, gamma_2R, gamma_3R;
  // double time, time_min, time_max, dt;


  double x[N];   //index = 0-(N-1)
  double u[N+4];  //index = 0-(N+3)
  double u_bar[N+4];
  double u_bar_dum[N+4];
  double u_1[N+4];
  double u_2[N+4];
  double uhp[N+1]; //index = 0-N
  // double uhR[N+1];

  double uhp_1[N];
  double uhp_2[N];
  double uhp_3[N];
  double beta_1[N];
  double beta_2[N];
  double beta_3[N];

  double w1_bar[N];
  double w2_bar[N];
  double w3_bar[N];
  double w_d[N];
  double w1[N];
  double w2[N];
  double w3[N];

  // double uhR_1[N];
  // double uhR_2[N];
  // double uhR_3[N];
  // double beta_1R[N];
  // double beta_2R[N];
  // double beta_3R[N];

  // double w1R_bar[N];
  // double w2R_bar[N];
  // double w3R_bar[N];
  // double wR_d[N];
  // double w1R[N];
  // double w2R[N];
  // double w3R[N];
  double RHS[N];

  /* Data write Initialize*/
  FILE *file1;
  char NAME[50];
/**************************/

  int Nt = round((time_max - time_min)/dt);
  /***************************/

  for (int i = 0; i < N; i++) {
  x[i] = (double)i*dx;
  u[i+2] = sin(x[i]);
  }

// Periodic boundary
  u[0] = u[N-1];
  u[1] = u[N];
  u[N+2] = u[3];
  u[N+3] = u[4];

  gamma_1 = 1.0/16.0;
  gamma_2 = 5.0/8.0;
  gamma_3 = 5.0/16.0;

  epslion = pow((1.0/10.0),6);


  for (int i = 2; i <= N+1; i++) {
    uhp_1[i-2] = (3.0/8.0)*u[i-2] - (5.0/4.0)*u[i-1] + (15.0/8.0)*u[i];
    uhp_2[i-2] = (-1.0/8.0)*u[i-1] + (3.0/4.0)*u[i] + (3.0/8.0)*u[i+1];
    uhp_3[i-2] = (3.0/8.0)*u[i] + (3.0/4.0)*u[i+1] - (1.0/8.0)*u[i+2];

    beta_1[i-2] = (1.0/3.0)*(4.0*u[i-2]*u[i-2] - 19.0*u[i-2]*u[i-1] + 25.0*u[i-1]*u[i-1]
                             + 11.0*u[i-2]*u[i] - 31.0*u[i-1]*u[i] + 10.0*u[i]*u[i]);
    beta_2[i-2] = (1.0/3.0)*(4.0*u[i-1]*u[i-1] - 13.0*u[i-1]*u[i] + 13.0*u[i]*u[i]
                             + 5.0*u[i-1]*u[i+1] - 13.0*u[i]*u[i+1] + 4.0*u[i+1]*u[i+1]);
    beta_3[i-2] = (1.0/3.0)*(10.0*u[i]*u[i] - 31.0*u[i]*u[i+1] + 25.0*u[i+1]*u[i+1]
                             + 11.0*u[i]*u[i+2] - 19.0*u[i+1]*u[i+2] + 4.0*u[i+2]*u[i+2]);

    w1_bar[i-2] = (gamma_1)/pow((epslion+beta_1[i-2]),2);
    w2_bar[i-2] = (gamma_2)/pow((epslion+beta_2[i-2]),2);
    w3_bar[i-2] = (gamma_3)/pow((epslion+beta_3[i-2]),2);

    w_d[i-2] = w1_bar[i-2] + w2_bar[i-2] + w3_bar[i-2];

    w1[i-2] = w1_bar[i-2]/w_d[i-2];
    w2[i-2] = w2_bar[i-2]/w_d[i-2];
    w3[i-2] = w3_bar[i-2]/w_d[i-2];

    uhp[i-1] = w1[i-2]*uhp_1[i-2] + w2[i-2]*uhp_2[i-2] + w3[i-2]*uhp_3[i-2];
  }

  uhp[0] = uhp[N];

  /*************************************************************************************/

  for (int i = 0; i < N; i++) {
    x[i] = (double)i*dx;
    u_bar[i+2] = sin(x[i]);
  }

  // sprintf(NAME,"output/i%d.dat",ts);
  file1 = fopen("output/time00.dat","w");
  Write_data(file1,x,u_bar);
  fclose(file1);

  for (int time = 0; time < Nt; time++) {

    weno_reconstruction(u_bar,RHS);

    for (int i = 2; i <= N+1; i++) {
      u_bar_dum[i] = u_bar[i];
      u_1[i] = u_bar[i] + dt*RHS[i-2];
    }

    weno_reconstruction(u_1,RHS);

    for (int i = 2; i <= N+1; i++) {
      u_2[i] = 3.0*u_bar_dum[i]/4.0 + u_1[i]/4.0 + dt * RHS[i-2]/4.0;
    }

    weno_reconstruction(u_2,RHS);

    for (int i = 2; i <= N+1; i++) {
      u_bar[i] = u_bar_dum[i]/3.0 + 2.0*u_2[i]/3.0 + dt * 2.0 * RHS[i-2]/3.0;
    }
    /*    DATA Print          */
    printf("Time = %d\n",time);
    sprintf(NAME,"output/time%d.dat",time);
    file1 = fopen(NAME,"w");
    Write_data(file1,x,u_bar);
    fclose(file1);
  }
}




void weno_reconstruction(double u_bar[], double RHS[]) {

  int i;
  // double gamma_1, gamma_2, gamma_3;
  double epslion;
  double gamma_1R, gamma_2R, gamma_3R;
  // double time, time_min, time_max, dt;


  // double x[N];
  // double u[N+4];
  // double u_bar[N+4];
  // double u_bar_dum[N+4];
  // double u_1[N+4];
  // double u_2[N+4];
  // double uhp[N+1];
  double uhR[N+1];

  // double uhp_1[N];
  // double uhp_2[N];
  // double uhp_3[N];
  // double beta_1[N];
  // double beta_2[N];
  // double beta_3[N];

  // double w1_bar[N];
  // double w2_bar[N];
  // double w3_bar[N];
  // double w_d[N];
  // double w1[N];
  // double w2[N];
  // double w3[N];

  double uhR_1[N];
  double uhR_2[N];
  double uhR_3[N];
  double beta_1R[N];
  double beta_2R[N];
  double beta_3R[N];

  double w1R_bar[N];
  double w2R_bar[N];
  double w3R_bar[N];
  double wR_d[N];
  double w1R[N];
  double w2R[N];
  double w3R[N];
//  double RHS[N];



  gamma_1R = 1.0/10.0;
  gamma_2R = 3.0/5.0;
  gamma_3R = 3.0/10.0;

  epslion = pow((1.0/10.0),6);

/* Periodic Boundary */
  u_bar[0] = u_bar[N-1];
  u_bar[1] = u_bar[N];
  u_bar[N+2] = u_bar[3];
  u_bar[N+3] = u_bar[4];



  for (int i = 2; i <= N+1; i++) {

    uhR_1[i-2] = (1.0/3.0)*u_bar[i-2] - (7.0/6.0)*u_bar[i-1] + (11.0/6.0)*u_bar[i];
    uhR_2[i-2] = (-1.0/6.0)*u_bar[i-1] + (5.0/6.0)*u_bar[i] + (1.0/3.0)*u_bar[i+1];
    uhR_3[i-2] = (1.0/3.0)*u_bar[i] + (5.0/6.0)*u_bar[i+1] - (1.0/6.0)*u_bar[i+2];

    beta_1R[i-2] = (13.0/12.0)*pow((u_bar[i-2] - 2.0*u_bar[i-1] + u_bar[i]),2)
                +(1.0/4.0)*pow((u_bar[i-2] - 4.0*u_bar[i-1] + 3.0*u_bar[i]),2);
    beta_2R[i-2] = (13.0/12.0)*pow((u_bar[i-1] - 2.0*u_bar[i] + u_bar[i+1]),2)
                +(1.0/4.0)*pow((u_bar[i-1] - u_bar[i+1] ),2);
    beta_3R[i-2] = (13.0/12.0)*pow((u_bar[i] - 2.0*u_bar[i+1] + u_bar[i+2]),2)
                +(1.0/4.0)*pow((3.0*u_bar[i] - 4.0*u_bar[i+1] + u_bar[i+2]),2);

    w1R_bar[i-2] = (gamma_1R)/pow((epslion+beta_1R[i-2]),2);
    w2R_bar[i-2] = (gamma_2R)/pow((epslion+beta_2R[i-2]),2);
    w3R_bar[i-2] = (gamma_3R)/pow((epslion+beta_3R[i-2]),2);

    wR_d[i-2] = w1R_bar[i-2] + w2R_bar[i-2] + w3R_bar[i-2];

    w1R[i-2] = w1R_bar[i-2]/wR_d[i-2];
    w2R[i-2] = w2R_bar[i-2]/wR_d[i-2];
    w3R[i-2] = w3R_bar[i-2]/wR_d[i-2];

    uhR[i-1] = w1R[i-2]*uhR_1[i-2] + w2R[i-2]*uhR_2[i-2] + w3R[i-2]*uhR_3[i-2];

  }

  uhR[0] = uhR[N];

  for (int i = 0; i < N; i++) {
    RHS[i] = - (uhR[i+1] - uhR[i]) / dx;
  }

}




void Write_data(FILE *file1, double x[], double u_bar[])
{
  for (int i = 0; i < N; i++) {
    fprintf(file1,"%g \t %g\n",x[i], u_bar[i+2]);
  }
  fflush(file1);
}
