#include "Bwd_l5p77_400016R15.h"

namespace S400016R15
{
float txfit_l5p77                             (float *x,int m){
    int ncoeff=  2;
    float avdat= -0.7499914E-03;
    float xmin[10]={
        -0.70706E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61191E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[  3]={
        -0.71681151E-02, 0.11636502E+00,
              0.      };
    int ientry=0;

    if (ientry==0){
        ientry=1;
        for(int i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
//  set up monomials   functions
    float x11 = x1;

//                 function

    float v_txfit_l5p77                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        ;

    return v_txfit_l5p77                             ;
}
float delta_l5p77                             (float *x,int m){
    int ncoeff= 27;
    float avdat=  0.1578709E-02;
    float xmin[10]={
        -0.70706E+00,-0.25517E-01,-0.65639E-01,-0.48941E-01,-0.14998E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61191E+00, 0.19623E-01, 0.51919E-01, 0.45656E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 28]={
        -0.49783653E-02, 0.48054159E-01, 0.42634504E-02, 0.25760327E-02,
         0.41639758E-02,-0.16310059E-02,-0.53034062E-02, 0.54592703E-03,
         0.34373938E-03, 0.42329361E-02, 0.52952487E-03, 0.16619867E-02,
        -0.48438041E-03, 0.34888133E-02,-0.46465420E-02, 0.45760351E-03,
        -0.24997706E-02,-0.22411499E-03,-0.34702087E-02,-0.10637053E-02,
         0.61008893E-03, 0.25165523E-02, 0.45019976E-03, 0.21699520E-02,
         0.86164189E-03,-0.58194902E-03,-0.86804741E-03,
              0.      };
    int ientry=0;

    if (ientry==0){
        ientry=1;
        for(int i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x13 = x12*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x51 = x5;

//                 function

    float v_delta_l5p77                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]*x12
        +coeff[  3]                *x51
        +coeff[  4]*x11*x21
        +coeff[  5]    *x22
        +coeff[  6]    *x21*x31
        +coeff[  7]    *x21
    ;
    v_delta_l5p77                             =v_delta_l5p77
        +coeff[  8]        *x31
        +coeff[  9]    *x21    *x41
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]*x11    *x31
        +coeff[ 13]        *x32
        +coeff[ 14]        *x31*x41
        +coeff[ 15]*x11*x22
        +coeff[ 16]    *x22*x31
    ;
    v_delta_l5p77                             =v_delta_l5p77
        +coeff[ 17]*x11    *x32
        +coeff[ 18]    *x21*x32
        +coeff[ 19]            *x41*x51
        +coeff[ 20]        *x32    *x51
        +coeff[ 21]            *x42
        +coeff[ 22]*x13
        +coeff[ 23]    *x22    *x41
        +coeff[ 24]    *x24
        +coeff[ 25]*x12*x21*x31
    ;
    v_delta_l5p77                             =v_delta_l5p77
        +coeff[ 26]*x11*x22*x31
        ;

    return v_delta_l5p77                             ;
}
float theta_l5p77                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.7331615E-03;
    float xmin[10]={
        -0.70706E+00,-0.25517E-01,-0.65639E-01,-0.48941E-01,-0.14998E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61191E+00, 0.19623E-01, 0.51919E-01, 0.45656E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.81565762E-02,-0.56474269E-02,-0.52074187E-01, 0.82109859E-02,
        -0.69648209E-02,-0.28562781E-02, 0.62991390E-02,-0.49963966E-02,
        -0.14625980E-02,-0.31614598E-01, 0.26740136E-01, 0.53838962E-02,
        -0.65623461E-02, 0.65444941E-02, 0.12026409E-01,-0.90852519E-02,
         0.16636739E-02, 0.13395566E-02, 0.36237102E-02,-0.70804032E-03,
        -0.17937729E-01,-0.41030951E-01, 0.11248547E-01, 0.67316987E-01,
        -0.32225769E-01, 0.18819688E-01,-0.68452102E-02,-0.11676223E-02,
        -0.14954571E-02,-0.48536980E-02, 0.76865479E-02, 0.31956437E-02,
        -0.69178189E-02, 0.14133289E-02,-0.62527051E-02,-0.10778537E-01,
         0.76834126E-02,-0.61449909E-03, 0.97391888E-03, 0.12626122E-01,
        -0.78572221E-02,-0.12231137E-02,-0.18905034E-01,-0.84468322E-02,
        -0.67249439E-02, 0.22040002E-01, 0.74832812E-02, 0.12277166E-01,
        -0.61917151E-02,-0.10023297E-01,
              0.      };
    int ientry=0;

    if (ientry==0){
        ientry=1;
        for(int i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x51 = x5;

//                 function

    float v_theta_l5p77                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]            *x41
        +coeff[  5]*x12
        +coeff[  6]*x11*x21
        +coeff[  7]    *x22
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[  8]*x11    *x31
        +coeff[  9]    *x21*x31
        +coeff[ 10]    *x21    *x41
        +coeff[ 11]    *x23
        +coeff[ 12]*x11*x21    *x41
        +coeff[ 13]                *x51
        +coeff[ 14]        *x31    *x51
        +coeff[ 15]            *x41*x51
        +coeff[ 16]*x11        *x41*x51
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 17]*x11        *x41
        +coeff[ 18]*x11*x22
        +coeff[ 19]*x12    *x31
        +coeff[ 20]    *x22*x31
        +coeff[ 21]    *x21*x32
        +coeff[ 22]    *x22    *x41
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]    *x21    *x42
        +coeff[ 25]    *x23*x31
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 26]*x11*x21*x32
        +coeff[ 27]        *x33*x41
        +coeff[ 28]*x11            *x51
        +coeff[ 29]    *x22        *x51
        +coeff[ 30]        *x32    *x51
        +coeff[ 31]        *x32
        +coeff[ 32]*x11    *x32
        +coeff[ 33]    *x21        *x51
        +coeff[ 34]        *x31*x41*x51
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 35]    *x22*x31    *x51
        +coeff[ 36]    *x23*x31    *x51
        +coeff[ 37]*x12*x21
        +coeff[ 38]        *x33
        +coeff[ 39]*x11    *x31*x41
        +coeff[ 40]*x11        *x42
        +coeff[ 41]*x11*x23
        +coeff[ 42]    *x22*x32
        +coeff[ 43]    *x21*x33
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 44]    *x23    *x41
        +coeff[ 45]    *x22*x31*x41
        +coeff[ 46]    *x21*x32*x41
        +coeff[ 47]*x11*x21    *x42
        +coeff[ 48]*x12*x21*x32
        +coeff[ 49]*x11*x22*x32
        ;

    return v_theta_l5p77                             ;
}
float phi_l5p77                               (float *x,int m){
    int ncoeff= 70;
    float avdat= -0.2398741E-02;
    float xmin[10]={
        -0.70706E+00,-0.25517E-01,-0.65639E-01,-0.48941E-01,-0.14998E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61191E+00, 0.19623E-01, 0.51919E-01, 0.45656E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 71]={
        -0.13976179E-02,-0.31211702E-01,-0.22858491E-02, 0.53233933E-02,
         0.16578278E-01, 0.16679488E-01,-0.18125884E-01,-0.10006619E-01,
         0.77650188E-02,-0.18048480E-02, 0.69258837E-02,-0.13907677E-02,
         0.46054861E-02, 0.65421709E-02,-0.26116692E-02, 0.32321739E-02,
         0.71469401E-02,-0.68336227E-02,-0.59263315E-02,-0.30937768E-02,
        -0.43955143E-02,-0.25933281E-01, 0.20497985E-01, 0.37296664E-02,
        -0.65897743E-03, 0.44823988E-03, 0.87754271E-03, 0.15141899E-01,
        -0.12467002E-01,-0.18372780E-02,-0.42931255E-03, 0.19603594E-03,
        -0.11385423E-02, 0.49205222E-02,-0.60501052E-02,-0.25983933E-01,
        -0.52419980E-02,-0.89130588E-02,-0.18064268E-01,-0.21098186E-02,
         0.20989105E-02,-0.12161921E-01, 0.11759900E-01, 0.19019049E-02,
         0.52825594E-02, 0.15789657E-02, 0.34372604E-02,-0.31847162E-02,
         0.24019531E-02, 0.37572403E-02,-0.10327641E-02, 0.42424139E-03,
        -0.11881156E-02, 0.24350598E-01,-0.87954830E-02,-0.27574881E-02,
         0.11496405E-02,-0.32424409E-03, 0.33459605E-02, 0.37657410E-01,
        -0.18636493E-01,-0.14948267E-01, 0.24083437E-01,-0.12118699E-01,
        -0.42583109E-02, 0.27185397E-02,-0.14521902E-02,-0.35441627E-02,
         0.57855956E-02,-0.19248476E-02,
              0.      };
    int ientry=0;

    if (ientry==0){
        ientry=1;
        for(int i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x13 = x12*x1;
    float x14 = x13*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x35 = x34*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;

//                 function

    float v_phi_l5p77                               =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]*x11
        +coeff[  4]*x11    *x31
        +coeff[  5]    *x21*x31
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x22
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[  8]    *x21*x32
        +coeff[  9]*x12    *x31
        +coeff[ 10]*x12        *x41
        +coeff[ 11]    *x23
        +coeff[ 12]*x11*x22*x31
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]        *x32    *x51
        +coeff[ 15]    *x21
        +coeff[ 16]        *x32
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 17]        *x31*x41
        +coeff[ 18]    *x21    *x41
        +coeff[ 19]*x11*x21
        +coeff[ 20]*x11*x21*x31
        +coeff[ 21]    *x22*x31
        +coeff[ 22]    *x22    *x41
        +coeff[ 23]*x11*x22
        +coeff[ 24]                *x51
        +coeff[ 25]*x14*x21
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 26]*x11            *x51
        +coeff[ 27]    *x21*x31    *x51
        +coeff[ 28]    *x21    *x41*x51
        +coeff[ 29]*x11*x21        *x51
        +coeff[ 30]        *x33    *x51
        +coeff[ 31]    *x24*x34
        +coeff[ 32]*x11*x21*x33    *x51
        +coeff[ 33]*x11    *x32
        +coeff[ 34]    *x21*x31*x41
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 35]    *x22*x32
        +coeff[ 36]    *x23*x31
        +coeff[ 37]*x11*x22    *x41
        +coeff[ 38]        *x35
        +coeff[ 39]        *x31    *x51
        +coeff[ 40]            *x41*x51
        +coeff[ 41]*x13    *x33
        +coeff[ 42]*x13    *x32*x41
        +coeff[ 43]        *x31*x41*x51
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 44]    *x21*x32    *x51
        +coeff[ 45]    *x22*x31    *x51
        +coeff[ 46]*x11*x21    *x41*x51
        +coeff[ 47]    *x23        *x51
        +coeff[ 48]    *x21    *x42
        +coeff[ 49]*x11*x21    *x41
        +coeff[ 50]*x12*x21
        +coeff[ 51]*x11        *x43
        +coeff[ 52]*x11*x21*x32
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 53]    *x22*x31*x41
        +coeff[ 54]    *x22    *x42
        +coeff[ 55]*x12*x21*x31
        +coeff[ 56]*x13*x21
        +coeff[ 57]*x12*x22
        +coeff[ 58]    *x24
        +coeff[ 59]        *x34*x41
        +coeff[ 60]        *x33*x42
        +coeff[ 61]*x11    *x34
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 62]*x11    *x33*x41
        +coeff[ 63]*x11    *x32*x42
        +coeff[ 64]    *x22*x33
        +coeff[ 65]*x12*x22*x31
        +coeff[ 66]    *x24*x31
        +coeff[ 67]*x11*x23    *x41
        +coeff[ 68]    *x24    *x41
        +coeff[ 69]*x11    *x35
        ;

    return v_phi_l5p77                               ;
}
float y00_l5p77                               (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.8507284E-03;
    float xmin[10]={
        -0.70706E+00,-0.25517E-01,-0.65639E-01,-0.48941E-01,-0.14998E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61191E+00, 0.19623E-01, 0.51919E-01, 0.45656E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
         0.11265546E-01,-0.50351784E-01, 0.89693256E-01,-0.57443939E-02,
        -0.10224757E-02,-0.12878326E-01, 0.96639236E-02,-0.48283283E-01,
        -0.58032971E-01, 0.13336676E-01, 0.20299055E-01, 0.31205101E-03,
         0.41882363E-02, 0.12323117E-01, 0.33031099E-01, 0.65613881E-01,
        -0.17241608E-01,-0.15005548E-01,-0.69857193E-02, 0.23756901E-01,
         0.10319058E-02, 0.13112995E-01, 0.54768495E-01,-0.10671430E-01,
        -0.17679501E-01,-0.58352537E-02, 0.16759699E-01,-0.16028875E-01,
         0.23450082E-01, 0.74566957E-02, 0.54125892E-04,-0.84056929E-02,
         0.62331962E-02,-0.13558227E-01, 0.10974362E-01,-0.23789873E-01,
         0.39941496E-02,-0.10024508E-01,-0.42256131E-02, 0.36216790E-02,
        -0.79320371E-01,-0.21203481E-01,-0.28698711E-01, 0.13872656E-02,
        -0.13463191E-02, 0.24311421E-02, 0.98399585E-03, 0.49600075E-02,
        -0.12430261E-01, 0.83605889E-02, 0.15125021E-01, 0.13346126E-02,
        -0.11234370E-02, 0.91716666E-02,-0.73113129E-02,-0.62612658E-02,
         0.57180561E-02,-0.29618612E-02, 0.14841467E-01, 0.75547458E-02,
              0.      };
    int ientry=0;

    if (ientry==0){
        ientry=1;
        for(int i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x13 = x12*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;

//                 function

    float v_y00_l5p77                               =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]*x11
        +coeff[  4]    *x21
        +coeff[  5]        *x32
        +coeff[  6]        *x31*x41
        +coeff[  7]*x11    *x31
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[  8]    *x21*x31
        +coeff[  9]*x11        *x41
        +coeff[ 10]    *x21    *x41
        +coeff[ 11]*x12
        +coeff[ 12]*x11*x21
        +coeff[ 13]    *x22
        +coeff[ 14]        *x33
        +coeff[ 15]        *x31*x42
        +coeff[ 16]            *x43
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 17]    *x21*x32
        +coeff[ 18]*x11    *x31*x41
        +coeff[ 19]    *x21*x31*x41
        +coeff[ 20]*x12    *x31
        +coeff[ 21]*x11*x21*x31
        +coeff[ 22]    *x22*x31
        +coeff[ 23]*x12        *x41
        +coeff[ 24]*x11*x21    *x41
        +coeff[ 25]*x11*x22
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 26]*x11    *x33
        +coeff[ 27]*x11    *x32*x41
        +coeff[ 28]    *x22*x32
        +coeff[ 29]*x11*x21*x31*x41
        +coeff[ 30]        *x31    *x51
        +coeff[ 31]    *x21        *x51
        +coeff[ 32]        *x32    *x51
        +coeff[ 33]        *x31*x41*x51
        +coeff[ 34]            *x42*x51
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 35]    *x21*x31    *x51
        +coeff[ 36]*x11*x21        *x51
        +coeff[ 37]    *x21*x32    *x51
        +coeff[ 38]*x11        *x42*x51
        +coeff[ 39]            *x42
        +coeff[ 40]        *x32*x41
        +coeff[ 41]    *x21    *x42
        +coeff[ 42]    *x22    *x41
        +coeff[ 43]                *x51
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 44]*x11    *x31*x42
        +coeff[ 45]*x11*x21    *x42
        +coeff[ 46]*x13    *x31
        +coeff[ 47]*x12*x21*x31
        +coeff[ 48]*x11*x22*x31
        +coeff[ 49]    *x23*x31
        +coeff[ 50]*x11*x22    *x41
        +coeff[ 51]            *x41*x51
        +coeff[ 52]*x11            *x51
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 53]    *x22*x33
        +coeff[ 54]*x11*x22*x32
        +coeff[ 55]*x12*x22    *x41
        +coeff[ 56]*x11*x23    *x41
        +coeff[ 57]*x11    *x31    *x51
        +coeff[ 58]    *x21    *x41*x51
        +coeff[ 59]*x13    *x31*x42
        ;

    return v_y00_l5p77                               ;
}
}
