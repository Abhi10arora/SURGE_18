#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main()
{
  FILE *f, *fp, *fp1, *f1;
  f = fopen("initial", "r");
  f1 = fopen("test_ini", "w");
  char* file[10] = {"0","1","2","3","4","5","6","7","8","9"};
  char* file_initial[3] = {"a","b","c"};

  int nx = 1024, ny = 1024, x, y, t, j, no_etas = 3, total_time = 50000, print_time = 5000, update_print_time = 5000;
  double dt = 0.1, dx = 1.0, dy = 1.0, K_conc = 3.5, K_eta = 0.625, M = 1.0, L = 1.0, w = 0.1, kx, ky, temp = 0;
  fftw_plan plan_f_conc, plan_f_g_conc, plan_b_conc, plan_b_g_conc, plan_f_eta[no_etas], plan_f_g_eta[no_etas], plan_b_eta[no_etas], plan_b_g_eta[no_etas];
  fftw_complex *conc, *g_conc;
  fftw_complex **eta, **g_eta;
  double g_alpha_a, g_alpha_b, g_alpha_c, g_alpha_d, g_beta_a, g_beta_b, g_beta_c, g_beta_d;

  eta = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*)*no_etas);
  g_eta = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*)*no_etas);
  for(j = 0; j<no_etas; j++)
  {
    eta[j] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
    g_eta[j] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  }
  conc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  g_conc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
  double* kpow2 = (double*)malloc(sizeof(double)*(nx*ny));
  double* kpow4 = (double*)malloc(sizeof(double)*(nx*ny));

  plan_f_conc = fftw_plan_dft_2d(nx ,ny , conc, conc, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_f_g_conc = fftw_plan_dft_2d(nx, ny, g_conc, g_conc, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_b_conc = fftw_plan_dft_2d(nx, ny, conc, conc, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_b_g_conc = fftw_plan_dft_2d(nx, ny, g_conc, g_conc, FFTW_BACKWARD, FFTW_ESTIMATE);
  for(j = 0; j<no_etas; j++)
  {
    plan_f_eta[j] = fftw_plan_dft_2d(nx ,ny , eta[j], eta[j], FFTW_FORWARD, FFTW_ESTIMATE);
    plan_f_g_eta[j] = fftw_plan_dft_2d(nx, ny, g_eta[j], g_eta[j], FFTW_FORWARD, FFTW_ESTIMATE);
    plan_b_eta[j] = fftw_plan_dft_2d(nx, ny, eta[j], eta[j], FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_b_g_eta[j] = fftw_plan_dft_2d(nx, ny, g_eta[j], g_eta[j], FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  for(x = 0; x<nx; x++)
  {
    if(x<nx/2)
    kx = 2.0*M_PI*(double)x/(double)(nx*dx);
    else
    kx = 2.0*M_PI*(double)(x-nx)/(double)(nx*dx);

    for(y = 0; y<ny; y++)
    {
      int i = x*nx + y;
      fscanf(f, "%lf", &temp);
      conc[i] = temp;
      // fprintf(f1, "%lf\n", creal(conc[i]));

      if(y<ny/2)
      ky = 2.0*M_PI*(double)y/(double)(ny*dx);
      else
      ky = 2.0*M_PI*(double)(y-ny)/(double)(ny*dy);

      kpow2[i] = ((kx*kx)+(ky*ky));
      kpow4[i] = kpow2[i]*kpow2[i];
    }
    // fprintf(f1, "\n");
  }

  for(j = 0; j<no_etas; j++)
  {
    fp = fopen(file_initial[j], "r");
    // fp1 = fopen(file[j], "w");
    for(x = 0; x<nx; x++)
    {
      for(y = 0; y<ny; y++)
      {
        int i = x*nx + y;
        fscanf(fp, "%lf", &temp);
        eta[j][i] = temp;
        // fprintf(fp1, "%lf\n", creal(eta[j][i]));
      }
      // fprintf(fp1, "\n");
    }
    fclose(fp);
    // fclose(fp1);
  }

  for(t = 1; t<=total_time; t++)
  {
    for(x = 0; x<nx; x++)
    {
      for(y = 0; y<ny; y++)
      {
        double sum = 0;
        int i = x*nx + y;
        for(j = 0; j<no_etas; j++)
        {
          sum = sum + pow(creal(eta[j][i]),3)*(6*pow(creal(eta[j][i]),2) - 15*creal(eta[j][i]) + 10);
        }
        g_conc[i] = sum*(3*g_alpha_a*pow(creal(conc[i]), 2) + 2*g_alpha_b*creal(conc[i]) + g_alpha_c) + (1-sum)*(3*g_beta_a*pow(creal(conc[i]), 2) + 2*g_beta_b*creal(conc[i]) + g_beta_c);
      }
    }
    for(j = 0; j<no_etas; j++)
    {
      for(x = 0; x<nx; x++)
      {
        for(y = 0; y<ny; y++)
        {
          int i = x*nx + y;
          double h_eta = 30*pow(creal(eta[j][i]),2)*(pow(creal(eta[j][i]),2) - 2*creal(eta[j][i]) +1);
          double sum = creal(eta[0][i]) + creal(eta[1][i]) + creal(eta[2][i]);

          g_eta[j][i] = h_eta*(g_alpha_a*pow(creal(conc[i]),3) + g_alpha_b*pow(creal(conc[i]),2) + g_alpha_c*creal(conc[i]) + g_alpha_d)
                        - h_eta*(g_beta_a*pow(creal(conc[i]),3) + g_beta_b*pow(creal(conc[i]),2) + g_beta_c*creal(conc[i]) + g_beta_d)
                        + w*(sum - creal(eta[j][i]));
        }
      }
    }

    fftw_execute(plan_f_conc);
    fftw_execute(plan_f_g_conc);
    for(j = 0; j<no_etas; j++)
    {
      fftw_execute(plan_f_eta[j]);
      fftw_execute(plan_f_g_eta[j]);
    }

    for(x = 0; x<nx; x++)
    {
      for(y = 0; y<ny; y++)
      {
        int i = x*nx + y;
        conc[i] = (conc[i] - M*dt*kpow2[i]*g_conc[i])/(1 + 2*M*K_conc*dt*kpow4[i]);
      }
    }
    for(j = 0; j<no_etas; j++)
    {
      for(x = 0; x<nx; x++)
      {
        for(y = 0; y<ny; y++)
        {
          int i = x*nx + y;
          eta[j][i] = (eta[j][i] - L*dt*g_eta[j][i])/(1 + 2*K_eta*L*dt*kpow2[i]);
        }
      }
    }

    fftw_execute(plan_b_conc);
    fftw_execute(plan_b_g_conc);
    for(j = 0; j<no_etas; j++)
    {
      fftw_execute(plan_b_eta[j]);
      fftw_execute(plan_b_g_eta[j]);
    }

    for(x = 0; x<nx; x++)
    {
      for(y = 0; y<ny; y++)
      {
        int i = x*nx + y;
        conc[i] = creal(conc[i])/((double)(nx*ny));
      //   if(t == print_time)
      //   {
      //     fprintf(fp, "%lf\n", creal(conc[i]));
      //   }
      }
      // if(t == print_time)
      // fprintf(fp, "\n");
    }
    for(j = 0; j<no_etas; j++)
    {
      for(x = 0; x<nx; x++)
      {
        for(y = 0; y<ny; y++)
        {
          int i = x*nx + y;
          eta[j][i] = creal(eta[j][i])/((double)(nx*ny));
        }
      }
    }
  }

  fclose(f);
  fclose(f1);
  free(kpow2);
  free(kpow4);
  fftw_destroy_plan(plan_f_conc);
  fftw_destroy_plan(plan_f_g_conc);
  fftw_destroy_plan(plan_b_conc);
  fftw_destroy_plan(plan_b_g_conc);
  for(j = 0; j<no_etas; j++)
  {
    fftw_destroy_plan(plan_f_eta[j]);
    fftw_destroy_plan(plan_f_g_eta[j]);
    fftw_destroy_plan(plan_b_eta[j]);
    fftw_destroy_plan(plan_b_g_eta[j]);
    fftw_free(eta[j]);
    fftw_free(g_eta[j]);
  }
  fftw_free(conc);
  fftw_free(g_conc);
  fftw_free(eta);
  fftw_free(g_eta);

  return 0;
}
