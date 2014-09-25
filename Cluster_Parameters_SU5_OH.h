#define Alpha   0.34
#define Pi 3.1415926
#define Atomnumber_SU5_OH   7
#define Chargegroupnumber_SU5_OH  4   
#define Molnumber_SU5_OH    42
#define Atomnumber_SOL   3
#define Chargegroupnumber_SOL  1
#define Molnumber_SOL    6202
#define LJ_Cutoff    1.2
#define Coul_Cutoff  1.2

#define m_CAA      14.0270
#define m_CBB      14.0270 
#define m_OAA      15.9994
#define m_HAA      1.0080 

#define m_OW      15.9994
#define m_HW      1.0080

#define q_CAA        0.000
#define q_CBB        +0.265
#define q_OAA        -0.683
#define q_HAA        +0.418

#define q_OW        -0.8476
#define q_HW         0.4238

#define C6_CAA       0.90975E-02
#define C12_CAA      0.35333E-04
#define C6_CBB       0.90975E-02
#define C12_CBB      0.35333E-04
#define C6_OAA       0.90975E-02
#define C12_OAA      0.35333E-04
#define C6_HAA        0.0000
#define C12_HAA       0.0000

#define C6_OW        0.26171E-02
#define C12_OW       0.26331E-05 
#define C6_HW        0.0000
#define C12_HW       0.0000
struct data_1
{
    double x;
    double y;
    double z;
    };
    
typedef struct data_1 Position;

struct data_2
{
     Position MOL_SU5_OH_1[Atomnumber_SU5_OH];
     };
     
struct data_4
{
     Position Chargegroup_SU5_OH_1[Chargegroupnumber_SU5_OH];
     };

struct data_2  MOL_SU5_OH_48[Molnumber_SU5_OH];
struct data_4  Chargegroup_SU5_OH_18[Molnumber_SU5_OH];

struct data_3
{
    Position MOL_SOL_1[Atomnumber_SOL];
    };

struct data_5
{
    Position Chargegroup_SOL_1;
};

struct data_3 MOL_SOL_3559[Molnumber_SOL];
struct data_5 Chargegroup_SOL_2881[Molnumber_SOL];


