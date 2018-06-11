#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>


void main()
{
 int ts,Nx,Ny,total_timesteps,timebreak,processors;
 double deltax,deltay,deltat;
 int LoopCount;
 FILE *fpr1,*fpr2,*fpr3,*fpw1,*fpw2,*fpw3;
double *AspectRatio, *Scalar;

 char tmp1[20];
 fpr1 = fopen("./Input/SystemSize.dat", "r");
 fscanf(fpr1,"%s%d %s%d %s%lf %s%lf ", &tmp1[20],&Nx,&tmp1[20],&Ny,&tmp1[20],&deltax,&tmp1[20],&deltay);
 fclose(fpr1);

 fpr2 = fopen("./Input/TimeInfo.dat", "r");
 fscanf(fpr2,"%s%d %s%d %s%d %s%lf ",&tmp1[20],&processors,&tmp1[20],&total_timesteps, &tmp1[20],&timebreak, &tmp1[20],&deltat);
 fclose(fpr2); 

 int i,j,ij,ijk;
 int *label,*Area;
 Scalar      = (double*)malloc(sizeof(double)*(Nx*Ny));
 AspectRatio = (double*)malloc(sizeof(double)*ts);
 label       = (int*)malloc(sizeof(int)*(Nx*Ny));
 Area       = (int*)malloc(sizeof(int)*(Nx*Ny));

double HalfMaxima,m1,m2,x1,x2,y1,y2,Rx,Ry;
int A,B,C,D,ijA,ijB,ijC,ijD,ijright,ijup; 
HalfMaxima = 0.1838;
char FileName[100], FileName2[100];


int count = 0;
int top,bottom,right,left,x,y,xy;
int min,temp,value;
double radius,AvgRadius,VolFrac;


sprintf(FileName2,"./Output/Comp/StatsVsTime.dat");
fpw3 = fopen(FileName2,"w");

// ---------- Time loop starts ------------------

for(ts=0;ts<=total_timesteps;ts=ts+(timebreak*10))
  {

//printf("step 1 cleared! \n");
//ts = 10000;
   sprintf(FileName2,"./Output/Comp/PptLabel%d.dat",ts);
   fpw1 = fopen(FileName2,"w");
   sprintf(FileName2,"./Output/Comp/PptSize%d.dat",ts);
   fpw2 = fopen(FileName2,"w");

   sprintf(FileName,"./Output/Comp/Comp%d.dat",ts);  // printf("%d %s\n",ts,FileName);
   if( (fpr1 = fopen(FileName,"r")) == NULL)
    {
     printf("Unable to open the corresponding file Comp%d.dat\n",ts);
     printf("Exiting\n");
     fclose(fpw1); exit(0);
    }
   else
    {
     fpr1 = fopen(FileName,"r");
    }
    
   for(i=0;i<Nx;i++)
     {
      for(j=0;j<Ny;j++)
        {
         ij = i*Ny+j;
         Scalar[ij] = 0.0;
         fscanf(fpr1,"%d %d %lf",&i,&j,&Scalar[ij] );
        }
     }
   fclose(fpr1);
//printf("step 2 cleared! \n");
//--------- initializing labels ---------------
   for(i=0;i<Nx;i++)
     {
      for(j=0;j<Ny;j++)
        {
         ij = i*Ny+j;  
         label[ij] = 0;
        }
     }
//printf("step 3 cleared! \n"); 
//------ allot labels for the first time -----------
   for(i=0;i<Nx;i++)
     {
      for(j=0;j<Ny;j++)
        {
         ij     =             i*Ny+          j;  
         bottom =             i*Ny+(j-1+Ny)%Ny;
         left   = ((i-1+Nx)%Nx)*Ny+          j;
         if(Scalar[ij]>HalfMaxima)
          {
           if( (label[bottom] == 0) && (label[left] == 0) )
            {
             count ++;
             label[ij] = count;
            }
           else if( (label[bottom] != 0) && (label[left] == 0) )
            {
             label[ij] = label[bottom];
            }
           else if( (label[bottom] == 0) && (label[left] != 0) )
            {
             label[ij] = label[left];
            }
           else if(label[bottom] > label[left])
            {
             label[ij] = label[left];
            }
           else
            {
             label[ij] = label[bottom];
            }
          }
        }
     }
//printf("step 4 cleared! \n");
//-------- refine labels --------
   for(i=0;i<Nx;i++)
     {
      for(j=0;j<Ny;j++)
        {
         ij     =             i*Ny+          j;  //printf("%d ",ij);
         top    =             i*Ny+(j+1+Ny)%Ny;
         right  = ((i+1+Nx)%Nx)*Ny+          j;
         bottom =             i*Ny+(j-1+Ny)%Ny;
         left   = ((i-1+Nx)%Nx)*Ny+          j;
         if(label[ij] !=0)
          {
           min = 0;
           value = label[ij]; min = label[ij];
           if( (label[top] < min) && (label[top] != 0) ) 
            {
             min = label[top];
            }
           if( (label[right] < min) && (label[right] != 0) )
            {
             min = label[right];
            }
           if( (label[bottom] < min) && (label[bottom] != 0) )
            {
             min = label[bottom];
            }
           if( (label[left] < min) && (label[left] != 0) )
            {
             min = label[left];
            }    
          // printf("%d %d %d %d %d %d %d\n",label[left],label[right],label[top],label[bottom],label[ij],value,min);
/*       
           if( ( (label[top] != 0) && (label[top] < label[ij]) && (label[right] == 0) ) ||( (label[right] != 0) && (label[top] != 0) && (label[right] > label[top]) ) )
            {
             min = label[top];
            }
           else if( ( (label[right] != 0) && (label[right] < label[ij]) && (label[top] == 0) ) || ( (label[right] != 0) && (label[top] != 0) && (label[right] < label[top]) ) )
            {
             min = label[right];
            }
*/
if(min<label[ij])
 {
           for(x=0;x<Nx;x++)
             {
              for(y=0;y<Ny;y++)
                {
                 xy = x*Nx+y;
                 if(label[xy] == value)
                  {
                   label[xy] = min;
                  }
                }
             }
 }
          }
       } // j loop
     } //i loop
//printf("step 5 cleared! \n");

//----- calculate area and keep count of number of precipitates
   for(i=0;i<=count;i++)
     {
      Area[i] = 0;             
     }

   for(i=0;i<Nx;i++)
     {
      for(j=0;j<Ny;j++)
        {
         ij = i*Ny+j;
         fprintf(fpw1,"%d %d %d \n",i,j,label[ij]); 
        }
       fprintf(fpw1,"\n");
     }

//printf("step 6 cleared! \n");
for(x=1;x<=count;x++)
  {
   for(i=0;i<Nx;i++)
     {
      for(j=0;j<Ny;j++)
        {
         ij = i*Ny+j;
         if(label[ij] == x)
          {
           Area[x]++;
          }
        //printf("%d ",ij); //exit(0);
        }
     }
 }
  //printf("step 7 cleared! \n");
temp = 0;
AvgRadius = 0.0;
VolFrac = 0.0;
   for(i=1;i<=count;i++)
     {
      if(Area[i] != 0)
       {
        temp++;
        radius = 0.0; radius = sqrt(Area[i]/M_PI);
        fprintf(fpw2,"%d %d %d %lf\n",temp,i,Area[i],radius);
        AvgRadius += radius;
        VolFrac += Area[i];
       }             
     }
AvgRadius = AvgRadius/temp;
VolFrac = VolFrac/(Nx*Ny);
fprintf(fpw3,"%d %d %lf %lf %lf\n",ts,temp,VolFrac,AvgRadius,(AvgRadius*AvgRadius*AvgRadius));
printf("time %d ppt count = %d \n",ts,temp);
fclose(fpw1);
fclose(fpw2);
} // t loop
//printf("step 8 cleared! \n");
fclose(fpw3);

free(AspectRatio);
free(Scalar);
free(label);
free(Area);


} // main loop

