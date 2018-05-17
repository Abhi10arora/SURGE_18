#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
  FILE *f;
  char* file[10] = {"0","1","2","3","4","5","6","7","8","9"};

  int t, x, j = 0, nx = 100, total_time = 400000, print_time = 40000, update_print_time = 40000;
  float dt = 0.0025, M = 1.0, K = 1.0, dx = 1.0, dy = 1.0, A = 1.0, temp_1 = 0, temp_2 = 0, temp = 0;

  float* conc = (float*)malloc((nx+2)*sizeof(float));
  for(x = 0; x<nx; x++)
  {
    if(x>nx*0.25 && x<nx*0.75)
    conc[x] = 1;
    else
    conc[x] = 0;
  }

  for(t = 0; t<=total_time; t++)
  {
    if(t == print_time)
    {
      f = fopen(file[j], "w");
    }
    conc[nx] = conc[0];
    conc[nx+1] = conc[1];
    temp_1 = conc[nx-1];
    temp_2 = conc[nx-2];
    for(x = 0; x<nx; x++)
    {
      temp = conc[x];
      conc[x] = conc[x] + (M*dt*A*12*(2*conc[x]-1)*pow(conc[x+1]-temp_1,2))/(4*dx*dx) + (M*dt*A*(2+12*pow(conc[x],2)-12*conc[x])*(conc[x+1]-2*conc[x]+temp_1))/(dx*dx) - (2*M*dt*K*(conc[x+2]-4*conc[x+1]+6*conc[x]-4*temp_1+temp_2))/pow(dx,4);
      temp_2 = temp_1;
      temp_1 = temp;
      if(t == print_time)
      {
        fprintf(f,"%f\n", conc[x]);
      }
    }
    if(t == print_time)
    {
      print_time = print_time + update_print_time;
      j++;
      fclose(f);
    }
  }
  free(conc);
  return 0;
}
