#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"

int main()
{
  FILE *f, *fp;
  f = fopen("initial", "w");
  char* file[10] = {"0","1","2","3","4","5","6","7","8","9"};

  int t = 0, x, y, j = 0, nx = 512, ny = 512, total_time = 80000, print_time = 40000, update_print_time = 40000;
  float dt = 0.0025, M = 1.0, K = 1.0, dx = 1.0, dy = 1.0, A = 1.0, temp = 0, total_noise = 0;
  float conc_initial = 0.5;

  float c_x_y_; // conc[x-1][y-1]
  float c_xpy_; // conc[x+1][y-1]
  float c_x_yp; // conc[x-1][y+1]
  float c_y_1; // conc[x][y-1]
  float c_y_2; // conc[x][y-2]
  long int *seed;
  long int a;
  seed = &a;
  *seed = -2;

  float** conc = (float**)malloc((nx+2)*sizeof(float*));
  float* c_x_1 = (float*)malloc((ny+2)*sizeof(float)); // array of conc[x-1][y]
  float* c_x_2 = (float*)malloc((ny+2)*sizeof(float)); // array of conc[x-2][y]
  for(x = 0; x<nx+2; x++)
  {
    conc[x] = (float*)malloc((ny+2)*sizeof(float));
    for(y = 0; y<ny+2; y++)
    {
      conc[x][y] = conc_initial + 0.02*ran2(seed) - 0.01;
    }
  }
  for(x = 0; x < nx; x++)
  {
    for(y = 0; y<ny; y++)
    {
      total_noise = total_noise + conc[x][y] - conc_initial;
    }
  }
  for(x = 0; x < nx; x++)
  {
    for(y = 0; y<ny; y++)
    {
      conc[x][y] = conc[x][y] - total_noise/(nx*ny);
         fprintf(f, "%f\n", conc[x][y]);
    }
       fprintf(f, "\n");
  }

  for(t = 0; t<=total_time; t++)
  {
    if(t == print_time)
    {
      fp = fopen(file[j], "w");
    }

    c_x_y_ = conc[nx-1][ny-1]; // when x == 0 && y == 0

    for(x = 0; x<nx+2; x++)
    {
      conc[x][ny] = conc[x][0];
      conc[x][ny+1] = conc[x][1];
    }
    for(y = 0; y<ny+2; y++)
    {
      conc[nx][y] = conc[0][y];
      conc[nx+1][y] = conc[1][y];
      c_x_1[y] = conc[nx-1][y];
      c_x_2[y] = conc[nx-2][y];
    }
    for(x = 0; x<nx; x++)
    {
      c_y_1 = conc[x][ny-1];
      c_y_2 = conc[x][ny-2];
      // c_x_1[ny] = c_x_1[0];
      // c_x_1[ny+1] = c_x_1[1];
      for(y = 0; y<ny; y++)
      {
        temp = conc[x][y];
        if(y == 0)
        c_xpy_ = conc[x+1][ny-1];
        else
        c_xpy_ = conc[x+1][y-1];

        if(x != 0 && y == ny-1)
        c_x_yp = conc[x-1][y+1];
        else
        c_x_yp = c_x_1[y+1];

        conc[x][y] = conc[x][y] + (M*dt*A*12*(2*conc[x][y]-1)*pow(conc[x+1][y]-c_x_1[y],2))/(4*dx*dx) + (M*dt*A*12*(2*conc[x][y]-1)*pow(conc[x][y+1]-c_y_1,2))/(4*dy*dy)
        + (M*dt*A*(2+12*pow(conc[x][y],2)-12*conc[x][y])*(conc[x+1][y]-2*conc[x][y]+c_x_1[y]))/(dx*dx) + (M*dt*A*(2+12*pow(conc[x][y],2)-12*conc[x][y])*(conc[x][y+1]-2*conc[x][y]+c_y_1))/(dy*dy)
        - (4*M*dt*K*(conc[x+1][y+1]-2*conc[x+1][y]+c_xpy_-2*conc[x][y+1]+4*conc[x][y]-2*c_y_1+c_x_yp-2*c_x_1[y]+c_x_y_))/(dx*dx*dy*dy)
        - (2*M*dt*K*(conc[x+2][y]-4*conc[x+1][y]+6*conc[x][y]-4*c_x_1[y]+c_x_2[y]))/pow(dx,4) - (2*M*dt*K*(conc[x][y+2]-4*conc[x][y+1]+6*conc[x][y]-4*c_y_1+c_y_2))/pow(dy,4);

        c_x_2[y] = c_x_1[y];
        c_x_y_ = c_x_1[y];
        c_x_1[y] = temp;
        c_y_2 = c_y_1;
        c_y_1 = temp;

        if(t == print_time)
        {
          fprintf(fp, "%f\n", conc[x][y]);
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
  // fclose(fp);
  for(x = 0; x<nx+2; x++)
  {
    free(conc[x]);
  }
  free(conc);
  free(c_x_2);
  free(c_x_1);
  return 0;
}
