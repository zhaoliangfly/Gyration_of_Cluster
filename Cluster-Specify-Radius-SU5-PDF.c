#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h"
#include "math.h"
#include "time.h"
#include "Cluster_Parameters_SU5_OH.h"
#include "Cluster_Functions_SU5_OH.h"

#define Specify_Cluster_Max_Number_Search_CM   0
#define Max_Cluster_Number      1
#define All_Cluster_Number      0

#define Specify_CM_Max_Cluster_Number_Gyration_Radius  1
#define predefine_max_cluster_number  15



#define TIME 0
#define TIME_FI  200000
#define Deltar   0.001
#define Dotsnumber  2200 
#define Pdf_Cutoff   2.2 
#define Cluster_Atoms_Distance   0.60
#define Cluster_CM_SU5_OH_Distance  0.9
#define Boundary_Margin 0.7

void SEAL(int i);
int jacobi(int n, double a[],double u[], double eps);

int natoms,step,natom,read_return;
float time_1,prec;
//float prec;
matrix box;
rvec *x;
XDRFILE *xtc;
 
int cluster;
int i,j;
          struct data
          {
                 int label;
                 int location_x;
                 int location_y;
                 int location_z;
          };
         struct data array[Molnumber_SU5_OH];   /** Notice in this array, the array[0] is also used to store the data! **/

int main ()
{
   long t1;
   t1=clock();
   
   char path_xtc[]="F:/GMXSimulationData/SuperComputation6202WatersRealSU5/42/traj.xtc";

    xtc=xdrfile_open (path_xtc,"r");
    read_xtc_natoms (path_xtc,&natoms);
    x = calloc(natoms, sizeof (x[0]));


  int time_initial=0;  /* The starting point of our investigation depends on the real time of the .xtc file! **/
  int time_statistic=TIME;
  double pdistribution[2201];
    for(i=0;i<=Dotsnumber;i++)
          pdistribution[i]=0.00;
  
  double All_Cluster[Molnumber_SU5_OH+1];
    for(i=0;i<=Molnumber_SU5_OH;i++)
         All_Cluster[i]=0.00;


while (1)
{
  read_return=read_xtc (xtc,natoms,&step,&time_1,box,x,&prec);
            if (read_return!=0) break;
   

  
#if 0
    for (natom=1;natom<=natoms;natom++)
    {
    printf ("%d %f %d %f %f %f\n",step,time,natom,x[natom-1][0],x[natom-1][1],x[natom-1][2]);
    }
#endif  
  
   

   
   if(time_statistic==time_initial)  
   {
         j=-1;
         for(i=0;i<=Molnumber_SU5_OH*Atomnumber_SU5_OH-1;i+=1)
         {
              if(i%Atomnumber_SU5_OH==0)
                   j++;
              MOL_SU5_OH_48[j].MOL_SU5_OH_1[i%Atomnumber_SU5_OH].x=x[i][0];
              MOL_SU5_OH_48[j].MOL_SU5_OH_1[i%Atomnumber_SU5_OH].y=x[i][1];
              MOL_SU5_OH_48[j].MOL_SU5_OH_1[i%Atomnumber_SU5_OH].z=x[i][2];
         }

         j=-1;
         for(i=Molnumber_SU5_OH*Atomnumber_SU5_OH;i<=Molnumber_SU5_OH*Atomnumber_SU5_OH+3*Molnumber_SOL-1;i+=1)
         {
              if((i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3==0)
                   j++;
              MOL_SOL_3559[j].MOL_SOL_1[(i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3].x=x[i][0];
              MOL_SOL_3559[j].MOL_SOL_1[(i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3].y=x[i][1];
              MOL_SOL_3559[j].MOL_SOL_1[(i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3].z=x[i][2];
          }



         for(i=0;i<=Molnumber_SU5_OH-1;i++)
         {
               array[i].label=0;
               array[i].location_x=0;
               array[i].location_y=0;
               array[i].location_z=0;
         }
         
     /** For the convienience of Computation  ***
      ** we define the CM matrix explicitly   ***
      ** We must to ensure the GMX will set   ***
      ** the molecules in the same way that   ***
      ** keeping the molecules with an entity ***
      ** !                                    ***
      **/
      
         Position CM_Of_SU5_OH[Molnumber_SU5_OH];
           for(i=0;i<=Molnumber_SU5_OH-1;i++)
           {
                 CM_Of_SU5_OH[i].x=CM_SU5_OH_X(&(MOL_SU5_OH_48[i].MOL_SU5_OH_1[0]));
                 CM_Of_SU5_OH[i].y=CM_SU5_OH_Y(&(MOL_SU5_OH_48[i].MOL_SU5_OH_1[0]));
                 CM_Of_SU5_OH[i].z=CM_SU5_OH_Z(&(MOL_SU5_OH_48[i].MOL_SU5_OH_1[0]));
           }

           
         
         Position CM_Of_SOL[Molnumber_SOL];
           for(i=0;i<=Molnumber_SOL-1;i++)
           {
               CM_Of_SOL[i].x=CM_SOL_X(&(MOL_SOL_3559[i].MOL_SOL_1[0]));
               CM_Of_SOL[i].y=CM_SOL_Y(&(MOL_SOL_3559[i].MOL_SOL_1[0]));
               CM_Of_SOL[i].z=CM_SOL_Z(&(MOL_SOL_3559[i].MOL_SOL_1[0]));
            }
         
        /** In order to compute the radius of the Max_Cluster_Number,   **
         ** we must first to translate all the CM of the molecules into **
         ** the central box!                                            **
         ** The symbol "CM" refers to the center of molecules consisting  **
         ** of the atoms with same mass rather the weighted atoms!      **
         **/ 

         for(i=0;i<=Molnumber_SU5_OH-1;i++)
         {
               if(CM_Of_SU5_OH[i].x>box[0][0])
               {
                     CM_Of_SU5_OH[i].x -= box[0][0];
                     for(j=0;j<=Atomnumber_SU5_OH-1;j++)
                            MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x -= box[0][0];
               }
               if(CM_Of_SU5_OH[i].x<box[0][0])
               {
                     CM_Of_SU5_OH[i].x += box[0][0];
                     for(j=0;j<=Atomnumber_SU5_OH-1;j++)
                            MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x += box[0][0];
               }

               if(CM_Of_SU5_OH[i].y>box[1][1])
               {
                     CM_Of_SU5_OH[i].y -= box[1][1];
                     for(j=0;j<=Atomnumber_SU5_OH-1;j++)
                            MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y -= box[1][1];
               }
               if(CM_Of_SU5_OH[i].y<box[1][1])
               {
                     CM_Of_SU5_OH[i].y += box[1][1];
                     for(j=0;j<=Atomnumber_SU5_OH-1;j++)
                            MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y += box[1][1];
               }

               if(CM_Of_SU5_OH[i].z>box[2][2])
               {
                     CM_Of_SU5_OH[i].z -= box[2][2];
                     for(j=0;j<=Atomnumber_SU5_OH-1;j++)
                            MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z -= box[2][2];
               }
               if(CM_Of_SU5_OH[i].z<box[2][2])
               {
                     CM_Of_SU5_OH[i].z += box[2][2];
                     for(j=0;j<=Atomnumber_SU5_OH-1;j++)
                            MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z += box[2][2];
               }
         }

        /** Finishing the translation!  **
         **/ 

             cluster=0;
              for(j=0;j<=Molnumber_SU5_OH-1;j++)
              {
                   if(array[j].label==0)   /***To prevent the repeated labeling******/
                        {
                           cluster++;
                           SEAL(j); 
                         }  
              }

    /** End of the Cluster Labeling **/
      
    #if 0
  
      for(i=0;i<=Molnumber_SU5_OH-1;i++)
               printf("    1SU5  A   %5d%8.3f%8.3f%8.3f\n",i+1,CM_Of_SU5_OH[i].x,CM_Of_SU5_OH[i].y,CM_Of_SU5_OH[i].z);

      printf("   %8.3f   %8.3f   %8.3f\n", box[0][0], box[1][1], box[2][2]);

    #endif
  

#if  Max_Cluster_Number

    int Recording_Cluster[Molnumber_SU5_OH+1];
    int Recording_Label;    /** This number can tell us which label corresponds to the max Cluster Number! **/
    for(i=0;i<=Molnumber_SU5_OH;i++)
        Recording_Cluster[i]=0;

    for(i=0;i<=Molnumber_SU5_OH-1;i++)
        Recording_Cluster[array[i].label]++;
    int max=1;
    for(i=0;i<=Molnumber_SU5_OH;i++)
    {
          if(Recording_Cluster[i]>=max)
          {
                max=Recording_Cluster[i];
                Recording_Label=i;
           }
     }

     //printf("%d %d\n",time_statistic, max);
#endif

#if Specify_CM_Max_Cluster_Number_Gyration_Radius
 //if(max==predefine_max_cluster_number)
 //{
    Position Center_Of_Cluster,Center_Of_Tails,Center_Of_Heads;
    Center_Of_Cluster.x=0.00;     Center_Of_Cluster.y=0.00;    Center_Of_Cluster.z=0.00;
    Center_Of_Tails.x=0.00;    Center_Of_Tails.y=0.00;   Center_Of_Tails.z=0.00;
    Center_Of_Heads.x=0.00;    Center_Of_Heads.y=0.00;   Center_Of_Heads.z=0.00;
  
    Position CG_Of_SU5_OH[Molnumber_SU5_OH], CG_Of_Tails[Molnumber_SU5_OH], CG_Of_Heads[Molnumber_SU5_OH];
    for(i=0;i<=Molnumber_SU5_OH-1;i++)
    {
          CG_Of_SU5_OH[i].x=0.00;   CG_Of_SU5_OH[i].y=0.00;   CG_Of_SU5_OH[i].z=0.00;
          CG_Of_Tails[i].x=0.00;   CG_Of_Tails[i].y=0.00;   CG_Of_Tails[i].z=0.00;
          CG_Of_Heads[i].x=0.00;   CG_Of_Heads[i].y=0.00;   CG_Of_Heads[i].z=0.00;
    }

    

    /**  We must first to translate all the relative **
     **  molecules into a unified form.              **
     **/

     for(i=0;i<=Molnumber_SU5_OH-1;i++)
     {
        if(array[i].label==Recording_Label)
        {
         
          for(j=0;j<=Atomnumber_SU5_OH-1;j++)
          {
                 MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x += array[i].location_x*box[0][0];
                 MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y += array[i].location_y*box[1][1];
                 MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z += array[i].location_z*box[2][2];

                 Center_Of_Cluster.x += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
                 Center_Of_Cluster.y += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
                 Center_Of_Cluster.z += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;

                 CG_Of_SU5_OH[i].x += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
                 CG_Of_SU5_OH[i].y += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
                 CG_Of_SU5_OH[i].z += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;

            /** We can also calculate the CG_Of_Tails[] or **
             ** CG_Of_Heads[]!                             **
             **/
             
              if(j<=3)
              {
                  Center_Of_Tails.x += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
                  Center_Of_Tails.y += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
                  Center_Of_Tails.z += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;

                 CG_Of_Tails[i].x += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
                 CG_Of_Tails[i].y += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
                 CG_Of_Tails[i].z += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;
              }

              if((j>=4)&&(j<=6))
              {   
                  Center_Of_Heads.x += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
                  Center_Of_Heads.y += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
                  Center_Of_Heads.z += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;

                 CG_Of_Heads[i].x += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
                 CG_Of_Heads[i].y += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
                 CG_Of_Heads[i].z += MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;
              }    
                         
          }
                 CG_Of_SU5_OH[i].x /= Atomnumber_SU5_OH;
                 CG_Of_SU5_OH[i].y /= Atomnumber_SU5_OH;
                 CG_Of_SU5_OH[i].z /= Atomnumber_SU5_OH;
                 CG_Of_Tails[i].x /= 4;
                 CG_Of_Tails[i].y /= 4;
                 CG_Of_Tails[i].z /= 4;
                 CG_Of_Heads[i].x /= 3;
                 CG_Of_Heads[i].y /= 3;
                 CG_Of_Heads[i].z /= 3;
                 
        }  
      }


     /** End of the easy translation!             **
      **/

    /** We can now calculate the CM of the Cluster! **
     **/

        Center_Of_Cluster.x /= max*Atomnumber_SU5_OH;
        Center_Of_Cluster.y /= max*Atomnumber_SU5_OH;
        Center_Of_Cluster.z /= max*Atomnumber_SU5_OH;
        Center_Of_Tails.x /= max*4;
        Center_Of_Tails.y /= max*4;
        Center_Of_Tails.z /= max*4;
        Center_Of_Heads.x /= max*3;
        Center_Of_Heads.y /= max*3;
        Center_Of_Heads.z /= max*3;

   /** In order to compute the CG of Tails and **
    ** Heads, we newly create the 2 array!     **
    ** These two arrays are the copies of the  **
    ** MOL_SU5_OH_48[].MOL_SU5_OH_1[]...       **
    **/
    double Tails[Molnumber_SU5_OH][Atomnumber_SU5_OH][3], Heads[Molnumber_SU5_OH][Atomnumber_SU5_OH][3];
    for(i=0;i<=Molnumber_SU5_OH-1;i++)
    {
        for(j=0;j<=Atomnumber_SU5_OH-1;j++)
        {
             Tails[i][j][0]=MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
             Tails[i][j][1]=MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
             Tails[i][j][2]=MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;
             Heads[i][j][0]=MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x;
             Heads[i][j][1]=MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y;
             Heads[i][j][2]=MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z;
        }
     }         
    
       for(i=0;i<=Molnumber_SU5_OH-1;i++)
       {
        if(array[i].label==Recording_Label)
        {
            /** The CG of the molecules, Tails and **
             ** Heads have been subjected to the   **
             ** translations!                      **
             **/
                 CG_Of_SU5_OH[i].x -= Center_Of_Cluster.x;
                 CG_Of_SU5_OH[i].y -= Center_Of_Cluster.y;
                 CG_Of_SU5_OH[i].z -= Center_Of_Cluster.z;
                 CG_Of_Tails[i].x -= Center_Of_Tails.x;
                 CG_Of_Tails[i].y -= Center_Of_Tails.y;
                 CG_Of_Tails[i].z -= Center_Of_Tails.z;
                 CG_Of_Heads[i].x -= Center_Of_Heads.x;
                 CG_Of_Heads[i].y -= Center_Of_Heads.y;
                 CG_Of_Heads[i].z -= Center_Of_Heads.z;
                 
          for(j=0;j<=Atomnumber_SU5_OH-1;j++)
          {
                 MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].x -= Center_Of_Cluster.x;
                 MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].y -= Center_Of_Cluster.y;
                 MOL_SU5_OH_48[i].MOL_SU5_OH_1[j].z -= Center_Of_Cluster.z;
                 Tails[i][j][0] -= Center_Of_Tails.x;
                 Tails[i][j][1] -= Center_Of_Tails.y;
                 Tails[i][j][2] -= Center_Of_Tails.z;
                 Heads[i][j][0] -= Center_Of_Heads.x;
                 Heads[i][j][1] -= Center_Of_Heads.y;
                 Heads[i][j][2] -= Center_Of_Heads.z;
          }
        }
      }

   #if 0
      printf("This is a gro file!\n");
      printf("%d\n",max);
     j=0;
     for(i=0;i<=Molnumber_SU5_OH-1;i++)
     {
           if(array[i].label==Recording_Label)
           {  
               j++;
               printf("    1SU5  A   %5d%8.3f%8.3f%8.3f\n",j,CG_Of_SU5_OH[i].x,CG_Of_SU5_OH[i].y,CG_Of_SU5_OH[i].z);
           }
      }

      printf("   %8.3f   %8.3f   %8.3f\n", box[0][0], box[1][1], box[2][2]);
  
    #endif

   #if 0
        double TEMP_CENTER[3];
            for(j=0;j<=2;j++)
               TEMP_CENTER[j]=0.00;


                  TEMP_CENTER[0]=Center_Of_Cluster.x;
                  TEMP_CENTER[1]=Center_Of_Cluster.y;
                  TEMP_CENTER[2]=Center_Of_Cluster.z;
   #endif

    //printf("%d %lf %lf %lf\n",time_statistic, Center_Of_Cluster.x, Center_Of_Cluster.y, Center_Of_Cluster.z);
  
   double Gyration[3][3],Gyration_Tails[3][3],Gyration_Heads[3][3];
   for(i=0;i<=2;i++)
      for(j=0;j<=2;j++)
      {
           Gyration[i][j]=0.00;
           Gyration_Tails[i][j]=0.00;
           Gyration_Heads[i][j]=0.00;
      }

  /** For the convienience of computation, we **
   ** store the labels.x, .y, .z into the     **
   ** matrix!                                 **
   **/
    
   int var_1,var_2;
   double TEMP_MOL[Molnumber_SU5_OH][Atomnumber_SU5_OH][3];
   for(var_1=0;var_1<=Molnumber_SU5_OH-1;var_1++)
   {
      for(var_2=0;var_2<=Atomnumber_SU5_OH-1;var_2++)
      {
            TEMP_MOL[var_1][var_2][0]=MOL_SU5_OH_48[var_1].MOL_SU5_OH_1[var_2].x;
            TEMP_MOL[var_1][var_2][1]=MOL_SU5_OH_48[var_1].MOL_SU5_OH_1[var_2].y;
            TEMP_MOL[var_1][var_2][2]=MOL_SU5_OH_48[var_1].MOL_SU5_OH_1[var_2].z;
      } 
   }

   for(i=0;i<=2;i++)
   {
      for(j=0;j<=2;j++)
      {
           for(var_1=0;var_1<=Molnumber_SU5_OH-1;var_1++)
           {
               if(array[var_1].label==Recording_Label)
               {
                  /** Calculation of the whole molecule ! **
                   **/
                    for(var_2=0;var_2<=Atomnumber_SU5_OH-1;var_2++)
                    {
                         Gyration[i][j] += (TEMP_MOL[var_1][var_2][i])*(TEMP_MOL[var_1][var_2][j]);
                       /** Tails are from 0 to 3 ! **
                        **/
                       if(var_2<=3)    
                         Gyration_Tails[i][j] += Tails[var_1][var_2][i]*Tails[var_1][var_2][j];
                         
                       /** Heads are from 4 to 6 ! **
                        **/
                       if(var_2>=4) 
                         Gyration_Heads[i][j] += Heads[var_1][var_2][i]*Heads[var_1][var_2][j];
                    }      
   
               }
           }
                         Gyration[i][j] /= max*Atomnumber_SU5_OH;
                         Gyration_Tails[i][j] /= max*4.0;
                         Gyration_Heads[i][j] /= max*3.0;
                         
       }
    }

    /** The codes below are cited from a professional **
     ** programming book! It is designed specially to **
     ** digonalizing a matrix of which for our aim is **
     ** named Gyration[i][j]!                         **
     **/
  
    double vv[3][3],vv_tail[3][3], vv_head[3][3];
    double eps=0.000001;
    int calnumber,calnumber_Tails,calnumber_Heads;
    calnumber=jacobi(3,&Gyration[0][0],&vv[0][0],eps);
    calnumber_Tails=jacobi(3,&Gyration_Tails[0][0],&vv_tail[0][0],eps);
    calnumber_Heads=jacobi(3,&Gyration_Heads[0][0],&vv_head[0][0],eps);
    //printf("%d %d %d\n", calnumber, calnumber_Tails, calnumber_Heads);

    if((calnumber<10000)&&(calnumber_Tails<10000)&&(calnumber_Heads<10000))
    {      
       int s1,s2;
       double temp=0.00;
       for(s1=0;s1<=1;s1++)
       {
            if(Gyration[s1+1][s1+1]<Gyration[s1][s1])
            {
                  temp=Gyration[s1][s1];
                  Gyration[s1][s1]=Gyration[s1+1][s1+1];
                  Gyration[s1+1][s1+1]=temp;
            }
        }

        if(Gyration[0][0]>Gyration[1][1])
        {
             temp=Gyration[0][0];
             Gyration[0][0]=Gyration[1][1];
             Gyration[1][1]=temp;
        }

        double asphericity=0.00,acylindricity=0.00,eccentricity=0.00,alpha=0.00,belta=0.00,kai=0.00;
        asphericity=Gyration[2][2]-0.5*(Gyration[1][1]+Gyration[0][0]);
        acylindricity=Gyration[1][1]-Gyration[0][0];
        eccentricity=1-(sqrt(Gyration[0][0])/sqrt((Gyration[0][0]+Gyration[1][1]+Gyration[2][2])*0.3333));
        kai=(asphericity*asphericity+0.75*acylindricity*acylindricity)/pow((Gyration[0][0]+Gyration[1][1]+Gyration[2][2]),2);
        //alpha=Gyration[1][1]/Gyration[2][2];
        //belta=Gyration[0][0]/Gyration[2][2];
        //printf("%d %lf %lf %lf %lf %lf %lf\n",time_statistic, sqrt(Gyration[0][0]+Gyration[1][1]+Gyration[2][2]),alpha,belta,eccentricity,asphericity, acylindricity);
        //printf("%d %d %lf %lf %lf %lf\n",time_statistic, max, sqrt(Gyration[0][0]+Gyration[1][1]+Gyration[2][2]),eccentricity,asphericity, acylindricity);   

       for(s1=0;s1<=1;s1++)
       {
            if(Gyration_Tails[s1+1][s1+1]<Gyration_Tails[s1][s1])
            {
                  temp=Gyration_Tails[s1][s1];
                  Gyration_Tails[s1][s1]=Gyration_Tails[s1+1][s1+1];
                  Gyration_Tails[s1+1][s1+1]=temp;
            }
        }

        if(Gyration_Tails[0][0]>Gyration_Tails[1][1])
        {
             temp=Gyration_Tails[0][0];
             Gyration_Tails[0][0]=Gyration_Tails[1][1];
             Gyration_Tails[1][1]=temp;
        }

       for(s1=0;s1<=1;s1++)
       {
            if(Gyration_Heads[s1+1][s1+1]<Gyration_Heads[s1][s1])
            {
                  temp=Gyration_Heads[s1][s1];
                  Gyration_Heads[s1][s1]=Gyration_Heads[s1+1][s1+1];
                  Gyration_Heads[s1+1][s1+1]=temp;
            }
        }

        if(Gyration_Heads[0][0]>Gyration_Heads[1][1])
        {
             temp=Gyration_Heads[0][0];
             Gyration_Heads[0][0]=Gyration_Heads[1][1];
             Gyration_Heads[1][1]=temp;
        }
        printf("%d %d %lf %lf %lf %lf %lf\n",time_statistic, max, sqrt(Gyration[0][0]+Gyration[1][1]+Gyration[2][2]),sqrt(Gyration_Tails[0][0]+Gyration_Tails[1][1]+Gyration_Tails[2][2]), sqrt(Gyration_Heads[0][0]+Gyration_Heads[1][1]+Gyration_Heads[2][2]),kai, asphericity);
        
     }     

 //}
 /** End of the predefine_max_cluster_number 
   **/    
#endif 

#if  All_Cluster_Number

    int Recording_Cluster[Molnumber_SU5_OH+1]; /** The label begins from 1 ! **/
    int Recording_Label;    /** This number can tell us which label corresponds to the max Cluster Number! **/
    for(i=0;i<=Molnumber_SU5_OH;i++)
        Recording_Cluster[i]=0;

    for(i=0;i<=Molnumber_SU5_OH-1;i++)
        Recording_Cluster[array[i].label]++;

    for(i=0;i<=Molnumber_SU5_OH;i++)
    { 
    	if(Recording_Cluster[i]!=0)	
        All_Cluster[Recording_Cluster[i]]++;
    }    
      

        
#endif






   
      

           /** Here is the end of the labeling! **/

           time_statistic++;
          if(time_statistic==TIME_FI+1)
               break;

    }  /** End of the time_statistic ***/

    time_initial++;

}   /** End of the loops ***/

#if All_Cluster_Number
      for(i=0;i<=Molnumber_SU5_OH;i++)
      {
      	 if(All_Cluster[i]!=0)
           printf("%d %f\n",i,All_Cluster[i]/((double)(TIME_FI-TIME)));
      }     
           
#endif           
  
    xdrfile_close (xtc);
    free(x);

    //getchar();

    long t2;
    t2=clock();
    printf("%d\n",t2-t1);

        return 0;


}   /** End of the main() ***/

     /****The sub-code below is the coral codes for labeling********/    
         void SEAL( int s)    
         {
             array[s].label=cluster;
             int neighbour[Molnumber_SU5_OH];
             for(i=0;i<=Molnumber_SU5_OH-1;i++)
                neighbour[i]=-1;
             int u=0;

             int t_1;
             Position TEMP_1[Atomnumber_SU5_OH];
             for(t_1=0;t_1<=Atomnumber_SU5_OH-1;t_1++)
             {
                      TEMP_1[t_1].x=MOL_SU5_OH_48[s].MOL_SU5_OH_1[t_1].x+array[s].location_x*box[0][0];
                      TEMP_1[t_1].y=MOL_SU5_OH_48[s].MOL_SU5_OH_1[t_1].y+array[s].location_y*box[1][1];
                      TEMP_1[t_1].z=MOL_SU5_OH_48[s].MOL_SU5_OH_1[t_1].z+array[s].location_z*box[2][2];
             }

             for(i=0;i<=Molnumber_SU5_OH-1;i++)    /****We try to judge whether  the molecules can be relative to objetive atom*****/ 
             {
                   if(i==s)
                       continue;
                   else
                   {
                       int px;
                       int py;
                       int pz;
                       int PX;
                       int PY;
                       int PZ;
                       double Dis=20.0;

                       int v1,v2;
                   for(v1=0;v1<=3;v1++)
                   {
                     for(v2=0;v2<=3;v2++)
                     {

                       for(px=-1;px<=1;px++)    /****We must take the peoredical boxes into considerations******/
                       {
                          for(py=-1;py<=1;py++)
                          {
                             for(pz=-1;pz<=1;pz++)
                             {

                                   if(distance(TEMP_1[v1].x,TEMP_1[v1].y,TEMP_1[v1].z,MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].x+px*box[0][0],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].y+py*box[1][1],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].z+pz*box[2][2])<=Dis)
                                   {
                                           Dis=distance(MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].x,MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].y,MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].z,MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].x+px*box[0][0],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].y+py*box[1][1],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].z+pz*box[2][2]);
                                           PX=px;
                                           PY=py;
                                           PZ=pz;
                                   }

                               }
                            }
                         }
                        }
                       }
       
                       
                                         
                                   if((Dis<=Cluster_Atoms_Distance)&&(array[i].label==0))
                   /****The condition at the last can ensure the finity of loops. And it is very important!**************/                  
                                   {
                                          array[i].location_x=PX;
                                          array[i].location_y=PY;
                                          array[i].location_z=PZ;
                                          neighbour[u]=i;
                                          
                                          array[i].label=cluster;
                                          u++;
                                          
                                   }
           
                        }
                     } 

         
               int k;
               
               for(k=0;k<=Molnumber_SU5_OH-1;k++)
               {      
                      if(array[neighbour[k]].label==-1)
                               break;
                      if(array[neighbour[k]].label==cluster)
                               SEAL(neighbour[k]);
               }
                
          }
