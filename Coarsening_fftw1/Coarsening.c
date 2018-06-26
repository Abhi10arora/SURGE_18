#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main()
{
  FILE *f, *fp, *fp1, *f1;
  f = fopen("initial", "w");
  // f1 = fopen("oneD", "w");
  char* file[10] = {"0","1","2","3","4","5","6","7","8","9"};
  // char* file_initial[3] = {"a","b","c"};

  int nx = 128, ny = 128, x, y, t, j = 0, no_etas = 3, total_time = 10000, print_time = 1000, update_print_time = 1000;
  double dt = 0.1, dx = 1.0, dy = 1.0, K_conc = 3.5, K_eta = 0.625, M = 1.0, L = 1.0, w = 0.1, kx, ky, temp = 0;
  fftw_plan plan_f_conc, plan_f_g_conc, plan_b_conc, plan_b_g_conc, plan_f_eta1, plan_f_g_eta1, plan_b_eta1, plan_b_g_eta1, plan_f_eta2, plan_f_g_eta2, plan_b_eta2, plan_b_g_eta2;
  fftw_complex *conc, *g_conc, *eta1, *g_eta1, *eta2, *g_eta2;
  double rad_1 = 20, rad_2 = 10;
  double g_alpha_a = 0, g_alpha_b = 1.0, g_alpha_c = -2.0, g_alpha_d = 1.0, g_beta_a = 0, g_beta_b = 1.0, g_beta_c = 0.0, g_beta_d = 0.0;

  eta1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  g_eta1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  eta2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  g_eta2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  conc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  g_conc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  double* kpow2 = (double*)malloc(sizeof(double)*(nx*ny));
  double* kpow4 = (double*)malloc(sizeof(double)*(nx*ny));

  plan_f_conc = fftw_plan_dft_2d(nx ,ny , conc, conc, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_f_g_conc = fftw_plan_dft_2d(nx, ny, g_conc, g_conc, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_b_conc = fftw_plan_dft_2d(nx, ny, conc, conc, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_b_g_conc = fftw_plan_dft_2d(nx, ny, g_conc, g_conc, FFTW_BACKWARD, FFTW_ESTIMATE);

  plan_f_eta1 = fftw_plan_dft_2d(nx ,ny , eta1, eta1, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_f_g_eta1 = fftw_plan_dft_2d(nx, ny, g_eta1, g_eta1, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_b_eta1 = fftw_plan_dft_2d(nx, ny, eta1, eta1, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_b_g_eta1 = fftw_plan_dft_2d(nx, ny, g_eta1, g_eta1, FFTW_BACKWARD, FFTW_ESTIMATE);

  plan_f_eta2 = fftw_plan_dft_2d(nx ,ny , eta2, eta2, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_f_g_eta2 = fftw_plan_dft_2d(nx, ny, g_eta2, g_eta2, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_b_eta2 = fftw_plan_dft_2d(nx, ny, eta2, eta2, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_b_g_eta2 = fftw_plan_dft_2d(nx, ny, g_eta2, g_eta2, FFTW_BACKWARD, FFTW_ESTIMATE);


  for(x = 0; x<nx; x++)
  {
    if(x<nx/2)
    kx = 2.0*M_PI*(double)x/(double)(nx*dx);
    else
    kx = 2.0*M_PI*(double)(x-nx)/(double)(nx*dx);

    for(y = 0; y<ny; y++)
    {
      int i = x*nx + y;
      if(pow(x-nx*0.25,2) + pow(y-ny/2,2) <= pow(rad_1,2))
      {
        conc[i] = 1.0;
        eta1[i] = 1.0;
      }
      else if(pow(x-nx*0.75,2) + pow(y-ny/2,2) <= pow(rad_2,2))
      {
        conc[i] = 1.0;
        eta2[i] = 1.0;
      }
      else
      {
        conc[i] = 0.05;
        eta1[i] = 0.0;
        eta2[i] = 0.0;
      }
      fprintf(f, "%d\t%d\t%lf\t%lf\t%lf\n", x, y, creal(conc[i]), creal(eta1[i]), creal(eta2[i]));

      if(y<ny/2)
      ky = 2.0*M_PI*(double)y/(double)(ny*dx);
      else
      ky = 2.0*M_PI*(double)(y-ny)/(double)(ny*dy);

      kpow2[i] = ((kx*kx)+(ky*ky));
      kpow4[i] = kpow2[i]*kpow2[i];
    }
    fprintf(f, "\n");
  }

  for(t = 1; t<=total_time; t++)
  {
    if(t == print_time)
    {
      fp = fopen(file[j], "w");
    }
    for(x = 0; x<nx; x++)
    {
      for(y = 0; y<ny; y++)
      {
        double sum_gconc = 0, sum_geta = 0;
        int i = x*nx + y;
        double h_eta1 = 30*pow(creal(eta1[i]),2)*(pow(creal(eta1[i]),2) - 2*creal(eta1[i]) +1);
        double h_eta2 = 30*pow(creal(eta2[i]),2)*(pow(creal(eta2[i]),2) - 2*creal(eta2[i]) +1);
        sum_gconc = pow(creal(eta1[i]),3)*(6*pow(creal(eta1[i]),2) - 15*creal(eta1[i]) + 10)
                  + pow(creal(eta2[i]),3)*(6*pow(creal(eta2[i]),2) - 15*creal(eta2[i]) + 10);
        sum_geta = creal(eta1[i]) + creal(eta2[i]);

        g_conc[i] = sum_gconc*(3*g_alpha_a*pow(creal(conc[i]), 2) + 2*g_alpha_b*creal(conc[i]) + g_alpha_c) + (1-sum_gconc)*(3*g_beta_a*pow(creal(conc[i]), 2) + 2*g_beta_b*creal(conc[i]) + g_beta_c);
        g_eta1[i] = h_eta1*(g_alpha_a*pow(creal(conc[i]),3) + g_alpha_b*pow(creal(conc[i]),2) + g_alpha_c*creal(conc[i]) + g_alpha_d)
                    - h_eta1*(g_beta_a*pow(creal(conc[i]),3) + g_beta_b*pow(creal(conc[i]),2) + g_beta_c*creal(conc[i]) + g_beta_d)
                    + w*(sum_geta - creal(eta1[i]));
        g_eta2[i] = h_eta2*(g_alpha_a*pow(creal(conc[i]),3) + g_alpha_b*pow(creal(conc[i]),2) + g_alpha_c*creal(conc[i]) + g_alpha_d)
                    - h_eta2*(g_beta_a*pow(creal(conc[i]),3) + g_beta_b*pow(creal(conc[i]),2) + g_beta_c*creal(conc[i]) + g_beta_d)
                    + w*(sum_geta - creal(eta2[i]));
      }
    }

    fftw_execute(plan_f_conc);
    fftw_execute(plan_f_g_conc);
    fftw_execute(plan_f_eta1);
    fftw_execute(plan_f_g_eta1);
    fftw_execute(plan_f_eta2);
    fftw_execute(plan_f_g_eta2);

    for(x = 0; x<nx; x++)
    {
      for(y = 0; y<ny; y++)
      {
        int i = x*nx + y;
        conc[i] = (conc[i] - M*dt*kpow2[i]*g_conc[i])/(1 + 2*M*K_conc*dt*kpow4[i]);
        eta1[i] = (eta1[i] - L*dt*g_eta1[i])/(1 + 2*K_eta*L*dt*kpow2[i]);
        eta2[i] = (eta2[i] - L*dt*g_eta2[i])/(1 + 2*K_eta*L*dt*kpow2[i]);
      }
    }

    fftw_execute(plan_b_conc);
    fftw_execute(plan_b_g_conc);
    fftw_execute(plan_b_eta1);
    fftw_execute(plan_b_g_eta1);
    fftw_execute(plan_b_eta2);
    fftw_execute(plan_b_g_eta2);

    for(x = 0; x<nx; x++)
    {
      for(y = 0; y<ny; y++)
      {
        int i = x*nx + y;
        conc[i] = creal(conc[i])/((double)(nx*ny));
        eta1[i] = creal(eta1[i])/((double)(nx*ny));
        eta2[i] = creal(eta2[i])/((double)(nx*ny));
        if(t == print_time)
        {
          fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\n", x, y, creal(conc[i]), creal(eta1[i]), creal(eta2[i]));
        }
      }
      if(t == print_time)
      fprintf(fp, "\n");
    }
    if(t == print_time)
    {
      print_time = print_time + update_print_time;
      j++;
      fclose(fp);
    }
  }

  fclose(f);
  // fclose(f1);
  free(kpow2);
  free(kpow4);
  fftw_destroy_plan(plan_f_conc);
  fftw_destroy_plan(plan_f_g_conc);
  fftw_destroy_plan(plan_b_conc);
  fftw_destroy_plan(plan_b_g_conc);
  fftw_destroy_plan(plan_f_eta1);
  fftw_destroy_plan(plan_f_g_eta1);
  fftw_destroy_plan(plan_b_eta1);
  fftw_destroy_plan(plan_b_g_eta1);
  fftw_destroy_plan(plan_f_eta2);
  fftw_destroy_plan(plan_f_g_eta2);
  fftw_destroy_plan(plan_b_eta2);
  fftw_destroy_plan(plan_b_g_eta2);
  fftw_free(conc);
  fftw_free(g_conc);
  fftw_free(eta1);
  fftw_free(g_eta1);
  fftw_free(eta2);
  fftw_free(g_eta2);

  return 0;
}
