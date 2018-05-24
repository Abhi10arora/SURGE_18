#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float radius(float c1, float x1, float c2, float x2)
{
  float rad = ((0.5-c1)*(x1-x2))/(c1-c2) + x1;
  return rad;
}

int main()
{
  FILE *f, *fp;
  f = fopen("conc", "w");
  fp = fopen("area", "w");

  int t, x, y, nx = 1024, ny = 1024, total_time = 4000000, print_time = 40000, update_print_time = 40000;
  float dt = 0.0025, M = 1.0, K = 1.0, dx = 1.0, dy = 1.0, A = 1.0, temp_1 = 0, temp_2 = 0, temp = 0;
  float rad_0 = 10, c11 = 0, c12 = 1, c21 = 0, c22 = 1, x11 = 0, x12 = 0, x21 = 0, x22 = 0, c_x_y_;

  fprintf(fp, "%f\t%f\n", 0.0,pow(rad_0,2));

  float** conc = (float**)malloc((nx+2)*sizeof(float*));
  float* c_xini_1 = (float*)malloc((ny+2)*sizeof(float));
  float* c_xini_2 = (float*)malloc((ny+2)*sizeof(float));
  float* c_yini_1 = (float*)malloc((nx+2)*sizeof(float));
  float* c_yini_2 = (float*)malloc((nx+2)*sizeof(float));
  for(x = 0; x<nx+2; x++)
  {
    conc[x] = (float*)malloc((ny+2)*sizeof(float));
    for(y = 0; y<ny+2; y++)
    {
      if((pow(y-ny/2,2) + pow(x-nx/2,2)) <= pow(rad_0,2))
      conc[x][y] = 1.0;
      else
      conc[x][y] = 0.1;
    }
  }

  for(t = 0; t<=total_time; t++)
  {
    c_x_y_ = 0;
    for(y = 0; y<ny+2; y++)
    {
      conc[nx][y] = conc[0][y];
      conc[nx+1][y] = conc[1][y];
      c_xini_1[y] = conc[nx-1][y];
      c_xini_2[y] = conc[nx-2][y];
    }
    for(x = 0; x<nx+2; x++)
    {
      conc[x][ny] = conc[x][0];
      conc[x][ny+1] = conc[x][1];
      c_yini_1[x] = conc[x][ny-1];
      c_yini_2[x] = conc[x][ny-2];
    }
    for(x = 0; x<nx; x++)
    {
      temp_1 = conc[x][ny-1];
      temp_2 = conc[x][ny-2];
      c_xini_1[ny] = c_xini_1[0];
      c_xini_1[ny+1] = c_xini_1[1];
      for(y = 0; y<ny; y++)
      {
        temp = conc[x][y];
        if(y == 0)
        {
          conc[x][y] = conc[x][y] + (M*dt*A*12*(2*conc[x][y]-1)*pow(conc[x+1][y]-c_xini_1[y],2))/(4*dx*dx) + (M*dt*A*12*(2*conc[x][y]-1)*pow(conc[x][y+1]-temp_1,2))/(4*dy*dy)
          + (M*dt*A*(2+12*pow(conc[x][y],2)-12*conc[x][y])*(conc[x+1][y]-2*conc[x][y]+c_xini_1[y]))/(dx*dx) + (M*dt*A*(2+12*pow(conc[x][y],2)-12*conc[x][y])*(conc[x][y+1]-2*conc[x][y]+temp_1))/(dy*dy)
          - (4*M*dt*K*(conc[x+1][y+1]-2*conc[x+1][y]+c_yini_1[x+1]-2*conc[x][y+1]+4*conc[x][y]-2*temp_1+c_xini_1[y+1]-2*c_xini_1[y]+c_x_y_))/(dx*dx*dy*dy)
          - (2*M*dt*K*(conc[x+2][y]-4*conc[x+1][y]+6*conc[x][y]-4*c_xini_1[y]+c_xini_2[y]))/pow(dx,4) - (2*M*dt*K*(conc[x][y+2]-4*conc[x][y+1]+6*conc[x][y]-4*temp_1+temp_2))/pow(dy,4);
        }
        else
        {
          conc[x][y] = conc[x][y] + (M*dt*A*12*(2*conc[x][y]-1)*pow(conc[x+1][y]-c_xini_1[y],2))/(4*dx*dx) + (M*dt*A*12*(2*conc[x][y]-1)*pow(conc[x][y+1]-temp_1,2))/(4*dy*dy)
          + (M*dt*A*(2+12*pow(conc[x][y],2)-12*conc[x][y])*(conc[x+1][y]-2*conc[x][y]+c_xini_1[y]))/(dx*dx) + (M*dt*A*(2+12*pow(conc[x][y],2)-12*conc[x][y])*(conc[x][y+1]-2*conc[x][y]+temp_1))/(dy*dy)
          - (4*M*dt*K*(conc[x+1][y+1]-2*conc[x+1][y]+conc[x+1][y-1]-2*conc[x][y+1]+4*conc[x][y]-2*temp_1+c_xini_1[y+1]-2*c_xini_1[y]+c_x_y_))/(dx*dx*dy*dy)
          - (2*M*dt*K*(conc[x+2][y]-4*conc[x+1][y]+6*conc[x][y]-4*c_xini_1[y]+c_xini_2[y]))/pow(dx,4) - (2*M*dt*K*(conc[x][y+2]-4*conc[x][y+1]+6*conc[x][y]-4*temp_1+temp_2))/pow(dy,4);
        }
        c_xini_2[y] = c_xini_1[y];
        c_x_y_ = c_xini_1[y];
        c_xini_1[y] = temp;
        temp_2 = temp_1;
        temp_1 = temp;

        if(t == print_time)
        {
          if(conc[x][y] > 0.1 && conc[x][y] < 0.5 && x < nx/2 && y == ny/2)
          {
            float min = 0.5 - conc[x][y];
            if(0.5-c11 > min)
            c11 = conc[x][y];
            x11 = x;
          }
          if(conc[x][y] > 0.5 && conc[x][y] < 1.0 && x < nx/2 && y == ny/2)
          {
            float min = conc[x][y] - 0.5;
            if(c12-0.5 > min)
            c12 = conc[x][y];
            x12 =  x;
          }
          if(conc[x][y] > 0.1 && conc[x][y] < 0.5 && x > nx/2 && y == ny/2)
          {
            float min = 0.5 - conc[x][y];
            if(0.5-c21 > min)
            c21 = conc[x][y];
            x21 = x;
          }
          if(conc[x][y] > 0.5 && conc[x][y] < 1.0 && x > nx/2 && y == ny/2)
          {
            float min = conc[x][y] - 0.5;
            if(c22-0.5 > min)
            c22 = conc[x][y];
            x22 = x;
          }
        }
        if(t == total_time)
        {
          fprintf(f, "%f\n", conc[x][y]);
        }
      }
      fprintf(f, "\n");
    }
    if(t == print_time)
    {
      float rad1 = radius(c11, x11, c12, x12);
      float rad2 = radius(c21, x21, c22, x22);
      float radius = (rad2 - rad1)/2;
      fprintf(fp, "%f\t%f\n", dt*t, pow(radius,2));
      print_time = print_time + update_print_time;
      c11 = 0; c12 = 1; c21 = 0; c22 = 1; x11 = 0; x12 = 0; x21 = 0; x22 = 0;
    }
  }
  fclose(f);
  fclose(fp);
  for(x = 0; x<nx+2; x++)
  {
    free(conc[x]);
  }
  free(conc);
  free(c_xini_2);
  free(c_xini_1);
  free(c_yini_1);
  free(c_yini_2);
  return 0;
}
