#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "Cluster_Parameters_SU5_OH.h"
#include "Cluster_Functions_SU5_OH.h"

double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double s;
	s=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	return s;
}

double  Distance(Position *Atom_1, Position *Atom_2)
{
     double distance;
     distance=sqrt((Atom_1->x-Atom_2->x)*(Atom_1->x-Atom_2->x)+(Atom_1->y-Atom_2->y)*(Atom_1->y-Atom_2->y)+(Atom_1->z-Atom_2->z)*(Atom_1->z-Atom_2->z));
     return distance;
}


double LJ_2Atoms_SU5_OH_SU5_OH_Cutoff(Position *A, Position *B, int label_A, int label_B)
{
    double energy_1;
    double r_1;
    double C6_A;
    double C12_A;
    double C6_B;
    double C12_B;
    if(label_A==0||label_A==1||label_A==2||label_A==3)
    {
      C6_A=C6_CAA;
      C12_A=C12_CAA;
    }
    
    if(label_A==4)
    {
       C6_A=C6_CBB;
       C12_A=C12_CBB;
    } 
    
    if(label_A==5)
    {
       C6_A=C6_OAA;
       C12_A=C12_OAA;
    }
    
    if(label_A==6)
    {
       C6_A=C6_HAA;
       C12_A=C12_HAA;
     }
     
    if(label_B==0||label_B==1||label_B==2||label_B==3)
    {
      C6_B=C6_CAA;
      C12_B=C12_CAA;
    }
    
    if(label_B==4)
    {
       C6_B=C6_CBB;
       C12_B=C12_CBB;
    } 
    
    if(label_B==5)
    {
       C6_B=C6_OAA;
       C12_B=C12_OAA;
    }
    
    if(label_B==6)
    {
       C6_B=C6_HAA;
       C12_B=C12_HAA;
     }
     
     r_1=Distance(A,B);
   if(r_1<=LJ_Cutoff)
   {
     energy_1=sqrt(C12_A*C12_B)/pow(r_1,12)-sqrt(C6_A*C6_B)/pow(r_1,6);
     return energy_1;
   }
   else
     return 0.000;
}     

double LJ_2Atoms_SU5_OH_SU5_OH(Position *A, Position *B, int label_A, int label_B)
{
    double energy_1;
    double r_1;
    double C6_A;
    double C12_A;
    double C6_B;
    double C12_B;
    if(label_A==0||label_A==1||label_A==2||label_A==3)
    {
      C6_A=C6_CAA;
      C12_A=C12_CAA;
    }
    
    if(label_A==4)
    {
       C6_A=C6_CBB;
       C12_A=C12_CBB;
    } 
    
    if(label_A==5)
    {
       C6_A=C6_OAA;
       C12_A=C12_OAA;
    }
    
    if(label_A==6)
    {
       C6_A=C6_HAA;
       C12_A=C12_HAA;
     }
     
    if(label_B==0||label_B==1||label_B==2||label_B==3)
    {
      C6_B=C6_CAA;
      C12_B=C12_CAA;
    }
    
    if(label_B==4)
    {
       C6_B=C6_CBB;
       C12_B=C12_CBB;
    } 
    
    if(label_B==5)
    {
       C6_B=C6_OAA;
       C12_B=C12_OAA;
    }
    
    if(label_B==6)
    {
       C6_B=C6_HAA;
       C12_B=C12_HAA;
     }
     
     r_1=Distance(A,B);

     energy_1=sqrt(C12_A*C12_B)/pow(r_1,12)-sqrt(C6_A*C6_B)/pow(r_1,6);
     return energy_1;

}    
     
double LJ_SU5_OH_SU5_OH_Cutoff(Position *C, Position *D)
{
     double energy_2=0;
     int label_C;
     int label_D;
     for(label_C=0;label_C<=Atomnumber_SU5_OH-1;label_C+=1)
     {
          for(label_D=0;label_D<=Atomnumber_SU5_OH-1;label_D+=1)
               energy_2=energy_2+LJ_2Atoms_SU5_OH_SU5_OH_Cutoff(C+label_C, D+label_D, label_C, label_D);
      }
      
      return energy_2;
} 

double LJ_SU5_OH_SU5_OH(Position *C, Position *D)
{
     double energy_2=0;
     int label_C;
     int label_D;
     for(label_C=0;label_C<=Atomnumber_SU5_OH-1;label_C+=1)
     {
          for(label_D=0;label_D<=Atomnumber_SU5_OH-1;label_D+=1)
               energy_2=energy_2+LJ_2Atoms_SU5_OH_SU5_OH(C+label_C, D+label_D, label_C, label_D);
      }
      
      return energy_2;
}   

double Coul_2Atoms_SU5_OH_SU5_OH_Cutoff(Position *E, Position *F, int label_E, int label_F)
{
    double energy_3;
    double q_E;
    double q_F;
    double r_3;

    if(label_E==0||label_E==1||label_E==2||label_E==3)
       q_E=q_CAA;
    if(label_E==4)
       q_E=q_CBB;
    if(label_E==5)
       q_E=q_OAA;
    if(label_E==6)
       q_E=q_HAA;
       
    if(label_F==0||label_F==1||label_F==2||label_F==3)
       q_F=q_CAA;
    if(label_F==4)
       q_F=q_CBB;
    if(label_F==5)
       q_F=q_OAA;
    if(label_F==6)
       q_F=q_HAA;
    
    r_3=Distance(E,F); 
   if(r_3<=Coul_Cutoff)
   {
     energy_3=138.935485*(q_E*q_F)/r_3;
     return energy_3;
    }
   else
     return 0.000; 
} 

double Coul_2Atoms_SU5_OH_SU5_OH(Position *E, Position *F, int label_E, int label_F)
{
    double energy_3;
    double q_E;
    double q_F;
    double r_3;

    if(label_E==0||label_E==1||label_E==2||label_E==3)
       q_E=q_CAA;
    if(label_E==4)
       q_E=q_CBB;
    if(label_E==5)
       q_E=q_OAA;
    if(label_E==6)
       q_E=q_HAA;
       
    if(label_F==0||label_F==1||label_F==2||label_F==3)
       q_F=q_CAA;
    if(label_F==4)
       q_F=q_CBB;
    if(label_F==5)
       q_F=q_OAA;
    if(label_F==6)
       q_F=q_HAA;
    
    r_3=Distance(E,F); 

     energy_3=138.935485*(q_E*q_F)/r_3;
     return energy_3;
}

double Coul_SU5_OH_SU5_OH_Cutoff(Position *G, Position *H)
{
     double energy_4=0;
     int label_G;
     int label_H;
     for(label_G=0;label_G<=Atomnumber_SU5_OH-1;label_G+=1)
     {
          for(label_H=0;label_H<=Atomnumber_SU5_OH-1;label_H+=1)
               energy_4=energy_4+Coul_2Atoms_SU5_OH_SU5_OH_Cutoff(G+label_G, H+label_H, label_G, label_H);
      }
      
      return energy_4;
}                                       
     
double LJ_2Atoms_SOL_SOL_Cutoff(Position *I, Position *J, int label_I, int label_J)
{
    double energy_5;
    double r_5;
    double C6_I;
    double C12_I;
    double C6_J;
    double C12_J;
    if(label_I==0)
    {
      C6_I=C6_OW;
      C12_I=C12_OW;
    }
    
    if(label_I==1||label_I==2)
    {
       C6_I=C6_HW;
       C12_I=C12_HW;
    } 
    
    if(label_J==0)
    {
      C6_J=C6_OW;
      C12_J=C12_OW;
    }
    
    if(label_J==1||label_J==2)
    {
       C6_J=C6_HW;
       C12_J=C12_HW;
    } 
     
     r_5=Distance(I,J);
     if(r_5<=LJ_Cutoff)
     {
          energy_5=sqrt(C12_I*C12_J)/pow(r_5,12)-sqrt(C6_I*C6_J)/pow(r_5,6);
          return energy_5;
     }
     else
          return 0.000;
} 

double LJ_SOL_SOL_Cutoff(Position *K, Position *L)
{
     double energy_6=0;
     int label_K;
     int label_L;
     for(label_K=0;label_K<=Atomnumber_SOL-1;label_K+=1)
     {
          for(label_L=0;label_L<=Atomnumber_SOL-1;label_L+=1)
               energy_6=energy_6+LJ_2Atoms_SOL_SOL_Cutoff(K+label_K, L+label_L, label_K, label_L);
      }
      
      return energy_6;
}

double Coul_2Atoms_SOL_SOL_Cutoff(Position *M, Position *N, int label_M, int label_N)
{
    double energy_7;
    double q_M;
    double q_N;
    double r_7;
    if(label_M==0)
       q_M=q_OW;
    if(label_M==1||label_M==2)
       q_M=q_HW;
    if(label_N==0)
       q_N=q_OW;
    if(label_N==1||label_N==2)
       q_N=q_HW;
    
    r_7=Distance(M,N);
    if(r_7<=Coul_Cutoff)
    {   
         energy_7=138.935485*(q_M*q_N)/r_7;
         return energy_7;
    }
    else
         return 0.000;
}

double Coul_SOL_SOL_Cutoff(Position *O, Position *P)
{
     double energy_8=0;
     int label_O;
     int label_P;
     for(label_O=0;label_O<=Atomnumber_SOL-1;label_O+=1)
     {
          for(label_P=0;label_P<=Atomnumber_SOL-1;label_P+=1)
               energy_8=energy_8+Coul_2Atoms_SOL_SOL_Cutoff(O+label_O, P+label_P, label_O, label_P);
      }
      return energy_8;
}                       

double LJ_2Atoms_SU5_OH_SOL_Cutoff(Position *Q, Position *R, int label_Q, int label_R)
{
    double energy_9;
    double r_9;
    double C6_Q;
    double C12_Q;
    double C6_R;
    double C12_R;
    if(label_Q==0||label_Q==1||label_Q==2||label_Q==3)
    {
      C6_Q=C6_CAA;
      C12_Q=C12_CAA;
    }
    
    if(label_Q==4)
    {
       C6_Q=C6_CBB;
       C12_Q=C12_CBB;
    }
    
    if(label_Q==5)
    {
       C6_Q=C6_OAA;
       C12_Q=C12_OAA;
     }

    if(label_Q==6)
    {
       C6_Q=C6_HAA;
       C12_Q=C12_HAA;
     }
     
    if(label_R==0)
    {
      C6_R=C6_OW;
      C12_R=C12_OW;
    }
    
    if(label_R==1||label_R==2)
    {
       C6_R=C6_HW;
       C12_R=C12_HW;
    } 
     
     r_9=Distance(Q,R);
     if(r_9<=LJ_Cutoff)
     {
           energy_9=sqrt(C12_Q*C12_R)/pow(r_9,12)-sqrt(C6_Q*C6_R)/pow(r_9,6);
           return energy_9;
     }
     else
           return 0.000;
}

double LJ_2Atoms_SU5_OH_SOL(Position *Q, Position *R, int label_Q, int label_R)
{
    double energy_9;
    double r_9;
    double C6_Q;
    double C12_Q;
    double C6_R;
    double C12_R;
    if(label_Q==0||label_Q==1||label_Q==2||label_Q==3)
    {
      C6_Q=C6_CAA;
      C12_Q=C12_CAA;
    }
    
    if(label_Q==4)
    {
       C6_Q=C6_CBB;
       C12_Q=C12_CBB;
    }
    
    if(label_Q==5)
    {
       C6_Q=C6_OAA;
       C12_Q=C12_OAA;
     }

    if(label_Q==6)
    {
       C6_Q=C6_HAA;
       C12_Q=C12_HAA;
     }
     
    if(label_R==0)
    {
      C6_R=C6_OW;
      C12_R=C12_OW;
    }
    
    if(label_R==1||label_R==2)
    {
       C6_R=C6_HW;
       C12_R=C12_HW;
    } 
     
     r_9=Distance(Q,R);
     energy_9=sqrt(C12_Q*C12_R)/pow(r_9,12)-sqrt(C6_Q*C6_R)/pow(r_9,6);
           
     return energy_9;

}

double LJ_SU5_OH_SOL_Cutoff(Position *S, Position *T)
{
     double energy_10=0;
     int label_S;
     int label_T;
     for(label_S=0;label_S<=Atomnumber_SU5_OH-1;label_S+=1)
     {
          for(label_T=0;label_T<=Atomnumber_SOL-1;label_T+=1)
               energy_10=energy_10+LJ_2Atoms_SU5_OH_SOL_Cutoff(S+label_S, T+label_T, label_S, label_T);
      }
      
      return energy_10;
}  

double LJ_SU5_OH_SOL(Position *S, Position *T)
{
     double energy_10=0;
     int label_S;
     int label_T;
     for(label_S=0;label_S<=Atomnumber_SU5_OH-1;label_S+=1)
     {
          for(label_T=0;label_T<=Atomnumber_SOL-1;label_T+=1)
               energy_10=energy_10+LJ_2Atoms_SU5_OH_SOL(S+label_S, T+label_T, label_S, label_T);
      }
      
      return energy_10;
}  
double Coul_2Atoms_SU5_OH_SOL_Cutoff(Position *U, Position *V, int label_U, int label_V)
{
    double energy_11;
    double q_U;
    double q_V;
    double r_11;

    if(label_U==0||label_U==1||label_U==2||label_U==3)
       q_U=q_CAA;
    if(label_U==4)
       q_U=q_CBB;
    if(label_U==5)
       q_U=q_OAA;
    if(label_U==6)
       q_U=q_HAA;

       
    if(label_V==0)
       q_V=q_OW;
    if(label_V==1||label_V==2)
       q_V=q_HW;
  
    r_11=Distance(U,V);
    if(r_11<=Coul_Cutoff) 
    {  
         energy_11=138.935485*(q_U*q_V)/r_11;
         return energy_11;
    }
    else
         return  0.000;
}

double Coul_2Atoms_SU5_OH_SOL(Position *U, Position *V, int label_U, int label_V)
{
    double energy_11;
    double q_U;
    double q_V;
    double r_11;

    if(label_U==0||label_U==1||label_U==2||label_U==3)
       q_U=q_CAA;
    if(label_U==4)
       q_U=q_CBB;
    if(label_U==5)
       q_U=q_OAA;
    if(label_U==6)
       q_U=q_HAA;

       
    if(label_V==0)
       q_V=q_OW;
    if(label_V==1||label_V==2)
       q_V=q_HW;
  
    r_11=Distance(U,V); 
    energy_11=138.935485*(q_U*q_V)/r_11;
    return energy_11;
}

double Coul_SU5_OH_SOL_Cutoff(Position *W, Position *X)
{
     double energy_12=0;
     int label_W;
     int label_X;
     for(label_W=0;label_W<=Atomnumber_SU5_OH-1;label_W+=1)
     {
          for(label_X=0;label_X<=Atomnumber_SOL-1;label_X+=1)
               energy_12=energy_12+Coul_2Atoms_SU5_OH_SOL_Cutoff(W+label_W, X+label_X, label_W, label_X);
      }
      
      return energy_12;
}

double CM_SU5_OH_X(Position *Y)
{
    double x1=0;
    int    ii;
    double m;
    double m_Total=0;
    for(ii=0;ii<=Atomnumber_SU5_OH-1;ii+=1)
    {
         if(ii==0||ii==1||ii==2||ii==3)
            m=m_CAA;
         if(ii==4)
            m=m_CBB;
         if(ii==5)
            m=m_OAA;
         if(ii==6)
            m=m_HAA;
         x1=x1+((Y+ii)->x)*m;
         m_Total=m_Total+m;
    }
    
    x1=x1/m_Total;
    return x1;
}  

double CQ_SU5_OH_Chargegroup5_X(Position *Y)
{
    double x1=0;
    int    ii;
    for(ii=4;ii<=6;ii+=1)
           x1=x1+((Y+ii)->x);

    x1=x1/3;
    return x1;
}

double CM_SU5_OH_Chargegroup5_X(Position *Y)
{
    double x1=0;
    int    ii;
    double m;
    double m_Total=0;
    for(ii=4;ii<=Atomnumber_SU5_OH-1;ii+=1)
    {
         if(ii==4)
            m=m_CBB;
         if(ii==5)
            m=m_OAA;
         if(ii==6)
            m=m_HAA;
         x1=x1+((Y+ii)->x)*m;
         m_Total=m_Total+m;
    }
    
    x1=x1/m_Total;
    return x1;
}  
    
double CM_SU5_OH_Y(Position *Z)
{
   double y=0;
   int jj;
   double m1;
   double m1_Total=0;
   for(jj=0;jj<=Atomnumber_SU5_OH-1;jj+=1)
   {
         if(jj==0||jj==1||jj==2||jj==3)
            m1=m_CAA;
         if(jj==4)
            m1=m_CBB;
         if(jj==5)
            m1=m_OAA;
         if(jj==6)
            m1=m_HAA;
         y=y+((Z+jj)->y)*m1;
         m1_Total=m1_Total+m1; 
   }
   
   y=y/m1_Total;
   return y;
}  

double CQ_SU5_OH_Chargegroup5_Y(Position *Z)
{
    double y=0;
    int    jj;
    for(jj=4;jj<=6;jj+=1)
           y=y+((Z+jj)->y);

    y=y/3;
    return y;
}

double CM_SU5_OH_Chargegroup5_Y(Position *Z)
{
   double y=0;
   int jj;
   double m1;
   double m1_Total=0;
   for(jj=4;jj<=Atomnumber_SU5_OH-1;jj+=1)
   {
         if(jj==4)
            m1=m_CBB;
         if(jj==5)
            m1=m_OAA;
         if(jj==6)
            m1=m_HAA;
         y=y+((Z+jj)->y)*m1;
         m1_Total=m1_Total+m1; 
   }
   
   y=y/m1_Total;
   return y;
}  
    
   
double CM_SU5_OH_Z(Position *X1)
{
   double z=0;
   int kk;
   double m2;
   double m2_Total=0;
   for(kk=0;kk<=Atomnumber_SU5_OH-1;kk+=1)
   {
         if(kk==0||kk==1||kk==2||kk==3)
            m2=m_CAA;
         if(kk==4)
            m2=m_CBB;
         if(kk==5)
            m2=m_OAA;
         if(kk==6)
            m2=m_HAA;
         z=z+((X1+kk)->z)*m2;
         m2_Total=m2_Total+m2; 
   }
   
   z=z/m2_Total;
   return z;
}   

double CQ_SU5_OH_Chargegroup5_Z(Position *Z2)
{
    double z=0;
    int    kk;
    for(kk=4;kk<=6;kk+=1)
           z=z+((Z2+kk)->z);

    z=z/3;
    return z;
}

double CM_SU5_OH_Chargegroup5_Z(Position *X1)
{
   double z=0;
   int kk;
   double m2;
   double m2_Total=0;
   for(kk=4;kk<=Atomnumber_SU5_OH-1;kk+=1)
   {
         if(kk==4)
            m2=m_CBB;
         if(kk==5)
            m2=m_OAA;
         if(kk==6)
            m2=m_HAA;
         z=z+((X1+kk)->z)*m2;
         m2_Total=m2_Total+m2; 
   }
   
   z=z/m2_Total;
   return z;
}   

double CM_SOL_X(Position *X2)
{
    double xx=0;
    int uu;
    double M;
    double M_Total=0;
    for(uu=0;uu<=Atomnumber_SOL-1;uu+=1)
    {
        if(uu==0)
          M=m_OW;
        if(uu==1||uu==2)
          M=m_HW;
        xx=xx+M*((X2+uu)->x);
        M_Total=M_Total+M;
    }
    
    xx=xx/M_Total;
    return xx;
}

double CQ_SOL_Chargegroup1_X(Position *X2)
{
    double xx=0;
    int uu;
    for(uu=0;uu<=Atomnumber_SOL-1;uu+=1)
        xx=xx+((X2+uu)->x);

    xx=xx/3;
    return xx;
}

double CM_SOL_Y(Position *X3)
{
    double yy=0;
    int vv;
    double M1;
    double M1_Total;
    for(vv=0;vv<=Atomnumber_SOL-1;vv+=1)
    {
        if(vv==0)
          M1=m_OW;
        if(vv==1||vv==2)
          M1=m_HW;
        yy=yy+M1*((X3+vv)->y);
        M1_Total=M1_Total+M1;
    }
    
    yy=yy/M1_Total;
    return yy;
}

double CQ_SOL_Chargegroup1_Y(Position *X3)
{
    double yy=0;
    int vv;
    for(vv=0;vv<=Atomnumber_SOL-1;vv+=1)
        yy=yy+((X3+vv)->y);
    
    yy=yy/3;
    return yy;
}            

double CM_SOL_Z(Position *X4)
{
    double zz=0;
    int ww;
    double M2;
    double M2_Total;
    for(ww=0;ww<=Atomnumber_SOL-1;ww+=1)
    {
        if(ww==0)
          M2=m_OW;
        if(ww==1||ww==2)
          M2=m_HW;
        zz=zz+M2*((X4+ww)->z);
        M2_Total=M2_Total+M2;
    }
    
    zz=zz/M2_Total;
    return zz;
} 
 
double CQ_SOL_Chargegroup1_Z(Position *X4)
{
    double zz=0;
    int ww;
    for(ww=0;ww<=Atomnumber_SOL-1;ww+=1)
        zz=zz+((X4+ww)->z);
    
    zz=zz/3;
    return zz;
} 

double Error(double r)
{ 
     double error_value;
     double delta_r=0.0001;
     int step;
     for(step=0;step*delta_r<=r;step++)
            error_value=error_value+pow(2.71823,-1*Alpha*step*delta_r*step*delta_r)*delta_r;
 
     error_value=error_value*2*sqrt(Alpha/Pi);
     return error_value;
}
   
double CM_SU5_OH_Tailgroup_X(Position *X1)
{
   double z=0;
   int kk;
   double m2;
   double m2_Total=0;
   for(kk=0;kk<=3;kk++)
   {
        if(kk==0||kk==1||kk==2||kk==3)
            m2=m_CAA;

         z=z+((X1+kk)->z)*m2;
         m2_Total=m2_Total+m2;
   }

   z=z/m2_Total;
   return z; 
}
     
double CM_SU5_OH_Tailgroup_Y(Position *Z)
{
   double y=0;
   int jj;
   double m1;
   double m1_Total=0;
   for(jj=0;jj<=3;jj+=1)
   {
         if(jj==0||jj==1||jj==2||jj==3)
            m1=m_CAA;

         y=y+((Z+jj)->y)*m1;
         m1_Total=m1_Total+m1; 
   }
   
   y=y/m1_Total;
   return y;
}  
 
double CM_SU5_OH_Tailgroup_Z(Position *X1)
{
   double z=0;
   int kk;
   double m2;
   double m2_Total=0;
   for(kk=0;kk<=3;kk+=1)
   {
         if(kk==0||kk==1||kk==2||kk==3)
            m2=m_CAA;

         z=z+((X1+kk)->z)*m2;
         m2_Total=m2_Total+m2; 
   }
   
   z=z/m2_Total;
   return z;
}                   
