#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"

int main()
{
  FILE *f, *fp, *f1;
  char* file[3] = {"a","b","c"};
  f = fopen("initial","w");
  // fp = fopen("tries","w");
  // f1 = fopen("centres","w");
  int nx = 1024, ny = 1024, x, y, i, j, total_centre = 2000, no_tries = 20000, no_centre = 0, counter = 0, region, rad_0 = 36, rad_sep = 400, count = 0;
  int  x_rand, y_rand, x_temp1, y_temp1, no_etas = 3, eta_n;
  float vol = 0, max_vol = 6, matrix_conc = 0.0, particle_conc = 1.0, eta_conc = 1.0;
  long int *seed;
  long int a;
  seed = &a;
  *seed = -2;

  float*** eta = (float***)malloc(no_etas*sizeof(float**));
  int* centre_x = (int*)malloc(total_centre*sizeof(int));
  int* centre_y = (int*)malloc(total_centre*sizeof(int));
  float** conc = (float**)malloc(nx*sizeof(float*));
  for(x = 0; x<nx; x++)
  {
    conc[x] = (float*)malloc(ny*sizeof(float));
    for(y = 0; y<ny; y++)
    {
      conc[x][y] = matrix_conc;
    }
  }
  for(i = 0; i<no_etas; i++)
  {
    eta[i] = (float**)malloc(nx*sizeof(float*));
    for(x = 0; x<nx; x++)
    {
      eta[i][x] = (float*)malloc(ny*sizeof(float));
      for(y = 0; y<ny; y++)
      {
        eta[i][x][y] = 0.0;
      }
    }
  }

  for(i = 0; i < no_tries; i++)
  {
    eta_n = floor(ran2(seed)*3);
    x_rand = round(ran2(seed)*1023);
    y_rand = round(ran2(seed)*1023);
    if(x_rand >= 13 && x_rand <= 1010 && y_rand >= 13 && y_rand <= 1010)
    {
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 5;
    }
    else if(x_rand < 13 && y_rand >= 13 && y_rand <= 1010)
    {
      x_temp1 = x_rand + nx;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_temp1-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 4;
    }
    else if(x_rand > 1010 && y_rand >= 13 && y_rand <= 1010)
    {
      x_temp1 = x_rand - nx;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_temp1-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 6;
    }
    else if(y_rand < 13 && x_rand >= 13 && x_rand <= 1010)
    {
      y_temp1 = y_rand + ny;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_rand-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 2;
    }
    else if(y_rand > 1010 && x_rand >= 13 && x_rand <= 1010)
    {
      y_temp1 = y_rand - ny;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_rand-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 8;
    }
    else if(x_rand < 13 && y_rand < 13)
    {
      x_temp1 = x_rand + nx;
      y_temp1 = y_rand + ny;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_rand-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep ||
            (pow(x_temp1-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_temp1-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 1;
    }
    else if(x_rand < 13 && y_rand > 1010)
    {
      x_temp1 = x_rand + nx;
      y_temp1 = y_rand - ny;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_rand-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep ||
            (pow(x_temp1-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_temp1-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 7;
    }
    else if(x_rand > 1010 && y_rand < 13)
    {
      x_temp1 = x_rand - nx;
      y_temp1 = y_rand + ny;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_rand-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep ||
            (pow(x_temp1-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_temp1-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 3;
    }
    else if(x_rand > 1010 && y_rand > 1010)
    {
      x_temp1 = x_rand - nx;
      y_temp1 = y_rand - ny;
      for(j = 0; j < no_centre; j++)
      {
        if((pow(x_rand-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_rand-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep ||
            (pow(x_temp1-centre_x[j], 2) + pow(y_rand-centre_y[j], 2)) < rad_sep || (pow(x_temp1-centre_x[j], 2) + pow(y_temp1-centre_y[j], 2)) < rad_sep)
        {
          counter = 1;
          break;
        }
      }
      region = 9;
    }
    if(counter == 1)
    {
      counter = 0;
      continue;
    }
    else
    {
      centre_x[no_centre] = x_rand;
      centre_y[no_centre] = y_rand;
      // fprintf(f1, "%d\t%d\n", x_rand, y_rand);
      no_centre++;
      switch (region)
      {
        case 1:
        {
          x_temp1 = x_rand + nx;
          y_temp1 = y_rand + ny;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_rand-x, 2) + pow(y_temp1-y, 2)) <= rad_0 ||
                  (pow(x_temp1-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_temp1-x, 2) + pow(y_temp1-y, 2)) <= rad_0)
                  {
                    conc[x][y] = particle_conc;
                    eta[eta_n][x][y] = eta_conc;
                  }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 2:
        {
          y_temp1 = y_rand + ny;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_rand-x, 2) + pow(y_temp1-y, 2)) <= rad_0)
              {
                conc[x][y] = particle_conc;
                eta[eta_n][x][y] = eta_conc;
              }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 3:
        {
          x_temp1 = x_rand - nx;
          y_temp1 = y_rand + ny;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_rand-x, 2) + pow(y_temp1-y, 2)) <= rad_0 ||
                  (pow(x_temp1-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_temp1-x, 2) + pow(y_temp1-y, 2)) <= rad_0)
                  {
                    conc[x][y] = particle_conc;
                    eta[eta_n][x][y] = eta_conc;
                  }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 4:
        {
          x_temp1 = x_rand + nx;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_temp1-x, 2) + pow(y_rand-y, 2)) <= rad_0)
              {
                conc[x][y] = particle_conc;
                eta[eta_n][x][y] = eta_conc;
              }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 5:
        {
             // printf("%d\n", 1);
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0)
              {
                conc[x][y] = particle_conc;
                eta[eta_n][x][y] = eta_conc;
              }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 6:
        {
          x_temp1 = x_rand - nx;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_temp1-x, 2) + pow(y_rand-y, 2)) <= rad_0)
              {
                conc[x][y] = particle_conc;
                eta[eta_n][x][y] = eta_conc;
              }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 7:
        {
          x_temp1 = x_rand + nx;
          y_temp1 = y_rand - ny;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_rand-x, 2) + pow(y_temp1-y, 2)) <= rad_0 ||
                  (pow(x_temp1-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_temp1-x, 2) + pow(y_temp1-y, 2)) <= rad_0)
                  {
                    conc[x][y] = particle_conc;
                    eta[eta_n][x][y] = eta_conc;
                  }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 8:
        {
          y_temp1 = y_rand - ny;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_rand-x, 2) + pow(y_temp1-y, 2)) <= rad_0)
              {
                conc[x][y] = particle_conc;
                eta[eta_n][x][y] = eta_conc;
              }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
        case 9:
        {
          x_temp1 = x_rand - nx;
          y_temp1 = y_rand - ny;
          for(x = 0; x < nx; x++)
          {
            for(y = 0; y < ny; y++)
            {
              if((pow(x_rand-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_rand-x, 2) + pow(y_temp1-y, 2)) <= rad_0 ||
                  (pow(x_temp1-x, 2) + pow(y_rand-y, 2)) <= rad_0 || (pow(x_temp1-x, 2) + pow(y_temp1-y, 2)) <= rad_0)
                  {
                    conc[x][y] = particle_conc;
                    eta[eta_n][x][y] = eta_conc;
                  }

              if(conc[x][y] == particle_conc) count++;
            }
          }
        }break;
      }
      vol = ((float)count/(nx*ny))*100;
      count = 0;
      if(vol >= max_vol)
      break;
    }
  }
  for(x = 0; x < nx; x++)
  {
    for(y = 0; y < ny; y++)
    {
      fprintf(f, "%f\n", conc[x][y]);
    }
    fprintf(f, "\n");
  }
  for(i = 0; i<no_etas; i++)
  {
    fp = fopen(file[i], "w");
    for(x = 0; x < nx; x++)
    {
      for(y = 0; y < ny; y++)
      {
        fprintf(fp, "%f\n", eta[i][x][y]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  printf("%f\n", vol);
  printf("%d\n", no_centre);
  for(x = 0; x<nx; x++)
  {
    free(conc[x]);
  }
  for(i = 0; i<no_etas; i++)
  {

    for(x = 0; x<nx; x++)
    {
      free(eta[i][x]);
    }
    free(eta[i]);
  }
  free(eta);
  free(conc);
  free(centre_x);
  free(centre_y);
  fclose(f);
  // fclose(f1);
  // fclose(fp);
  return 0;
}
