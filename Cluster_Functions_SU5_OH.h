double distance(double x1, double y1, double z1, double x2, double y2, double z2);
double Distance(Position *Atom_1, Position *Atom_2);
double Error(double r);

double LJ_2Atoms_SU5_OH_SU5_OH_Cutoff(Position *A, Position *B, int label_A, int label_B);
double LJ_2Atoms_SU5_OH_SU5_OH(Position *A, Position *B, int label_A, int label_B);
double LJ_SU5_OH_SU5_OH_Cutoff(Position *C, Position *D);
double LJ_SU5_OH_SU5_OH(Position *C, Position *D);
double Coul_2Atoms_SU5_OH_SU5_OH_Cutoff(Position *E, Position *F, int label_E, int label_F);
double Coul_2Atoms_SU5_OH_SU5_OH(Position *E, Position *F, int label_E, int label_F);
double Coul_SU5_OH_SU5_OH_Cutoff(Position *G, Position *H);
double LJ_2Atoms_SOL_SOL_Cutoff(Position *I, Position *J, int label_I, int label_J);
double LJ_2Atoms_SOL_SOL(Position *I, Position *J, int label_I, int label_J);
double LJ_SOL_SOL_Cutoff(Position *K, Position *L);
double Coul_2Atoms_SOL_SOL_Cutoff(Position *M, Position *N, int label_M, int label_N);
double Coul_2Atoms_SOL_SOL(Position *M, Position *N, int label_M, int label_N);
double Coul_SOL_SOL_Cutoff(Position *O, Position *P);
double LJ_2Atoms_SU5_OH_SOL_Cutoff(Position *Q, Position *R, int label_Q, int label_R);
double LJ_2Atoms_SU5_OH_SOL(Position *Q, Position *R, int label_Q, int label_R);
double LJ_SU5_OH_SOL_Cutoff(Position *S, Position *T);
double LJ_SU5_OH_SOL(Position *S, Position *T);
double Coul_2Atoms_SU5_OH_SOL_Cutoff(Position *U, Position *V, int label_U, int label_V);
double Coul_2Atoms_SU5_OH_SOL(Position *U, Position *V, int label_U, int label_V);
double Coul_SU5_OH_SOL_Cutoff(Position *W, Position *X);

double CM_SU5_OH_X(Position *Y);
double CM_SU5_OH_Y(Position *Z);
double CM_SU5_OH_Z(Position *X1);
double CM_SOL_X(Position *X2);
double CM_SOL_Y(Position *X3);
double CM_SOL_Z(Position *X4);

double CQ_SU5_OH_Chargegroup5_X(Position *Y);
double CQ_SU5_OH_Chargegroup5_Y(Position *Z);
double CQ_SU5_OH_Chargegroup5_Z(Position *X1);
double CQ_SOL_Chargegroup1_X(Position *X2);
double CQ_SOL_Chargegroup1_Y(Position *X3);
double CQ_SOL_Chargegroup1_Z(Position *X4);

double CM_SU5_OH_Tailgroup_X(Position *X1);
double CM_SU5_OH_Tailgroup_Y(Position *Z);
double CM_SU5_OH_Tailgroup_Z(Position *X1);

double CM_SU5_OH_Chargegroup5_X(Position *Y);
double CM_SU5_OH_Chargegroup5_Y(Position *Z);
double CM_SU5_OH_Chargegroup5_Z(Position *X1);
