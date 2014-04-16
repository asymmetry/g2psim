#include "Bwd_l5p77_484816R15.h"

namespace S484816R15
{
  float txfit_l5p77                             (float *x,int m){
    int ncoeff=  2;
    float avdat=  0.2301719E-02;
    float xmin[10]={
        -0.64213E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[  3]={
        -0.61417315E-02, 0.11132351E+00,
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
    int ncoeff= 24;
    float avdat=  0.1806338E-02;
    float xmin[10]={
        -0.70676E+00,-0.14199E+00,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.11708E+00, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 25]={
        -0.55805203E-02, 0.46362121E-01,-0.38440626E-01, 0.25946326E-02,
         0.70649199E-01, 0.17127871E-02,-0.25737366E-01, 0.32485768E-02,
        -0.17152056E-01, 0.19063091E-01,-0.46577128E-02, 0.13716113E-02,
         0.76822750E-03,-0.27712427E-01, 0.17454170E-01,-0.19477982E-01,
         0.59087313E-03,-0.70724229E-03, 0.24144603E-01, 0.27827884E-02,
         0.41169554E-03,-0.18073943E-02, 0.19482767E-02,-0.81350130E-03,
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
        +coeff[  5]    *x21
        +coeff[  6]    *x22
        +coeff[  7]        *x32
    ;
    v_delta_l5p77                             =v_delta_l5p77
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x21    *x41
        +coeff[ 10]        *x31*x41
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]            *x41
        +coeff[ 13]    *x21*x31
        +coeff[ 14]*x11    *x31*x41
        +coeff[ 15]    *x21*x31*x41
        +coeff[ 16]            *x42*x51
    ;
    v_delta_l5p77                             =v_delta_l5p77
        +coeff[ 17]        *x31
        +coeff[ 18]*x11    *x31
        +coeff[ 19]            *x42
        +coeff[ 20]*x13
        +coeff[ 21]*x11            *x51
        +coeff[ 22]    *x21        *x51
        +coeff[ 23]            *x41*x51
        ;

    return v_delta_l5p77                             ;
}
float theta_l5p77                             (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.1079361E-02;
    float xmin[10]={
        -0.70676E+00,-0.14199E+00,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.11708E+00, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.89760758E-02, 0.26518983E+00,-0.30147970E+00,-0.16516908E-02,
        -0.68220288E-01, 0.10328466E+00,-0.33962116E-01, 0.64165043E-02,
         0.11395771E-02, 0.16839723E+00,-0.18981226E+00, 0.41584983E-02,
        -0.15703995E-01, 0.21142804E-02, 0.10030282E-01, 0.38751900E+00,
        -0.59100832E-02, 0.65394957E-02,-0.15180048E+00, 0.17056018E+00,
        -0.15870733E-01, 0.72081965E+00,-0.38013077E+00, 0.17758415E-02,
         0.14910127E-01,-0.83619561E-02, 0.27232796E-02,-0.42041591E+00,
        -0.34775296E+00, 0.36046317E+00,-0.46538478E-02, 0.12247270E-02,
         0.75528747E-02, 0.11336331E-01,-0.33976433E+00, 0.31582156E+00,
        -0.61238980E+00, 0.29014546E+00,-0.42595020E-02,-0.47356826E-02,
         0.49268809E-03, 0.37244367E-02, 0.31088439E-02, 0.77589885E-02,
         0.17871469E-01,-0.21302927E-01, 0.57030655E-02,-0.54498017E-02,
         0.15321963E-01,-0.14660956E-01,
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
    float x43 = x42*x4;
    float x51 = x5;

//                 function

    float v_theta_l5p77                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]            *x41
        +coeff[  4]*x12
        +coeff[  5]*x11*x21
        +coeff[  6]    *x22
        +coeff[  7]                *x51
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[  8]        *x31
        +coeff[  9]*x11    *x31
        +coeff[ 10]    *x21*x31
        +coeff[ 11]            *x43
        +coeff[ 12]*x11    *x32*x41
        +coeff[ 13]*x11            *x51
        +coeff[ 14]        *x31    *x51
        +coeff[ 15]*x11    *x33*x41
        +coeff[ 16]        *x33*x42
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 17]        *x32
        +coeff[ 18]*x11        *x41
        +coeff[ 19]    *x21    *x41
        +coeff[ 20]        *x31*x41
        +coeff[ 21]*x11*x21*x31
        +coeff[ 22]    *x22*x31
        +coeff[ 23]        *x33
        +coeff[ 24]    *x21*x32*x41
        +coeff[ 25]            *x41*x51
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 26]*x11*x23    *x41
        +coeff[ 27]    *x21*x33*x41
        +coeff[ 28]*x11    *x32*x42
        +coeff[ 29]    *x21*x32*x42
        +coeff[ 30]*x11    *x31    *x51
        +coeff[ 31]        *x32    *x51
        +coeff[ 32]        *x33*x43
        +coeff[ 33]            *x42
        +coeff[ 34]*x12    *x31
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 35]*x12        *x41
        +coeff[ 36]*x11*x21    *x41
        +coeff[ 37]    *x22    *x41
        +coeff[ 38]    *x21        *x51
        +coeff[ 39]*x11*x23*x31
        +coeff[ 40]*x12            *x51
        +coeff[ 41]    *x21*x31    *x51
        +coeff[ 42]*x11        *x41*x51
        +coeff[ 43]*x12    *x33*x41
    ;
    v_theta_l5p77                             =v_theta_l5p77
        +coeff[ 44]*x12*x22    *x42
        +coeff[ 45]*x11*x23    *x42
        +coeff[ 46]*x11        *x42*x51
        +coeff[ 47]    *x21    *x42*x51
        +coeff[ 48]        *x33*x41*x51
        +coeff[ 49]        *x32*x42*x51
        ;

    return v_theta_l5p77                             ;
}
float phi_l5p77                               (float *x,int m){
    int ncoeff= 70;
    float avdat= -0.2389198E-02;
    float xmin[10]={
        -0.70676E+00,-0.14199E+00,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.11708E+00, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 71]={
         0.73294697E-03,-0.29323922E-01,-0.38734819E-02, 0.30638243E-02,
         0.25506027E-02, 0.46345540E-02,-0.43044104E-02,-0.26002793E-01,
         0.46644673E-01,-0.24319105E-01,-0.55266094E+00, 0.34358278E+00,
        -0.11951858E-03,-0.25375397E-02, 0.82476856E-02,-0.27947050E+00,
        -0.21764202E-01, 0.26257660E-01,-0.70737326E+00,-0.74286926E+00,
         0.41042957E+00,-0.72959776E-03,-0.11094283E-01,-0.33788979E-02,
         0.38308281E-01,-0.46228209E+00, 0.10341206E+01,-0.57900506E+00,
        -0.42749097E-03,-0.19630446E-03,-0.34679431E-01, 0.57933331E-01,
         0.68035103E-01, 0.26186211E-02,-0.14962874E-01, 0.17969711E-01,
         0.12962069E-01,-0.13791754E-01,-0.41386494E-02,-0.37103381E-01,
        -0.77449745E-02, 0.62984377E+00,-0.35454360E+00, 0.12502778E+01,
         0.30348673E-02, 0.16150456E-01, 0.39287176E-03,-0.33003844E-02,
         0.79820509E-03,-0.13175854E-02, 0.14894949E-01,-0.33962376E-01,
        -0.77330723E-03, 0.39675817E-01,-0.12546374E+00,-0.70146341E-02,
        -0.35471821E-03, 0.94335701E-03, 0.35570469E-01, 0.15731642E-01,
        -0.10844129E-01, 0.29652534E-01,-0.34474470E-01, 0.17555293E-01,
         0.99828164E-03, 0.11051828E+00,-0.23970830E+00, 0.12991732E+00,
        -0.19850605E-02,-0.15788807E-01,
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
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;

//                 function

    float v_phi_l5p77                               =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]*x11
        +coeff[  4]    *x21
        +coeff[  5]        *x32
        +coeff[  6]        *x31*x41
        +coeff[  7]*x11    *x31
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[  8]    *x21*x31
        +coeff[  9]*x11        *x41
        +coeff[ 10]*x12    *x31
        +coeff[ 11]*x12        *x41
        +coeff[ 12]                *x51
        +coeff[ 13]            *x42
        +coeff[ 14]    *x21    *x41
        +coeff[ 15]*x12
        +coeff[ 16]*x11        *x42
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]    *x22*x31
        +coeff[ 19]*x11*x21    *x41
        +coeff[ 20]    *x22    *x41
        +coeff[ 21]    *x24
        +coeff[ 22]*x11    *x34
        +coeff[ 23]*x11    *x33*x41
        +coeff[ 24]    *x21        *x51
        +coeff[ 25]*x12    *x34
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 26]*x11*x21*x34
        +coeff[ 27]    *x22*x34
        +coeff[ 28]        *x32    *x51
        +coeff[ 29]        *x31*x41*x51
        +coeff[ 30]*x11    *x31    *x51
        +coeff[ 31]*x14*x22
        +coeff[ 32]*x12*x24
        +coeff[ 33]        *x35*x43
        +coeff[ 34]*x11    *x32    *x51
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 35]    *x21*x32    *x51
        +coeff[ 36]*x11*x22        *x51
        +coeff[ 37]    *x23        *x51
        +coeff[ 38]        *x35*x44
        +coeff[ 39]*x11*x21    *x42*x51
        +coeff[ 40]        *x35    *x51
        +coeff[ 41]*x11*x21
        +coeff[ 42]    *x22
        +coeff[ 43]*x11*x21*x31
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 44]        *x33*x41
        +coeff[ 45]*x12    *x32
        +coeff[ 46]*x13    *x31
        +coeff[ 47]    *x23*x31
        +coeff[ 48]        *x35
        +coeff[ 49]            *x45
        +coeff[ 50]    *x21*x34
        +coeff[ 51]*x11            *x51
        +coeff[ 52]            *x42*x51
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 53]    *x21*x31    *x51
        +coeff[ 54]*x13*x23
        +coeff[ 55]            *x43*x51
        +coeff[ 56]*x14*x23
        +coeff[ 57]*x11    *x32*x41*x51
        +coeff[ 58]*x12        *x42*x51
        +coeff[ 59]        *x34*x41*x51
        +coeff[ 60]        *x33
        +coeff[ 61]        *x32*x41
    ;
    v_phi_l5p77                               =v_phi_l5p77
        +coeff[ 62]        *x31*x42
        +coeff[ 63]            *x43
        +coeff[ 64]*x11    *x32
        +coeff[ 65]*x12*x21
        +coeff[ 66]*x11*x22
        +coeff[ 67]    *x23
        +coeff[ 68]*x11    *x33
        +coeff[ 69]*x11*x21*x32
        ;

    return v_phi_l5p77                               ;
}
float y00_l5p77                               (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.8433543E-03;
    float xmin[10]={
        -0.70676E+00,-0.14199E+00,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.11708E+00, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
         0.51898803E-02,-0.45000613E-01, 0.83123080E-01,-0.21205880E-01,
         0.17449414E-01,-0.55657234E-02, 0.33911907E-02, 0.26452346E-02,
         0.21202976E+00,-0.28974342E+00,-0.78234538E-01, 0.10419550E+00,
         0.33739161E+00,-0.76296061E+00, 0.43020010E+00, 0.27601756E-01,
        -0.58417872E-01, 0.38904984E-01,-0.78678560E-02, 0.10476364E+01,
        -0.23654768E+01, 0.13305575E+01,-0.40879053E+00, 0.94235259E+00,
        -0.55490744E+00, 0.33707710E-03, 0.74281573E-01,-0.12171090E+00,
         0.47966007E-01,-0.81393969E+00,-0.40344306E-03, 0.45571040E-01,
        -0.51744234E-01, 0.46187114E-01,-0.52700669E-01,-0.11114685E-01,
         0.15646171E-01,-0.33334427E+01, 0.13017499E+01, 0.23931463E-02,
        -0.11326126E-01, 0.40250045E-03,-0.51122308E-02, 0.25377404E-02,
        -0.29752667E-02, 0.28543179E+01,-0.73232176E-03,-0.30405615E-01,
        -0.48195686E-01, 0.13214016E+00,-0.37833545E-01, 0.10928995E-02,
         0.19268262E-02, 0.29751265E-02,-0.18696904E-01, 0.33975657E-01,
        -0.15377354E-01,-0.19122839E-01, 0.47707442E-01,-0.36763873E-01,
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
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
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
        +coeff[  7]            *x42
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[  8]*x11    *x31
        +coeff[  9]    *x21*x31
        +coeff[ 10]*x11        *x41
        +coeff[ 11]    *x21    *x41
        +coeff[ 12]*x12
        +coeff[ 13]*x11*x21
        +coeff[ 14]    *x22
        +coeff[ 15]        *x33
        +coeff[ 16]        *x32*x41
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 17]        *x31*x42
        +coeff[ 18]            *x43
        +coeff[ 19]*x12    *x31
        +coeff[ 20]*x11*x21*x31
        +coeff[ 21]    *x22*x31
        +coeff[ 22]*x12        *x41
        +coeff[ 23]*x11*x21    *x41
        +coeff[ 24]    *x22    *x41
        +coeff[ 25]                *x51
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 26]    *x21*x32*x41
        +coeff[ 27]    *x21*x31*x42
        +coeff[ 28]    *x21    *x43
        +coeff[ 29]*x13    *x31
        +coeff[ 30]        *x31    *x51
        +coeff[ 31]*x11            *x51
        +coeff[ 32]    *x21        *x51
        +coeff[ 33]*x11    *x31    *x51
        +coeff[ 34]    *x21*x31    *x51
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 35]*x12    *x32
        +coeff[ 36]*x12    *x31*x41
        +coeff[ 37]*x11*x22*x31
        +coeff[ 38]    *x23*x31
        +coeff[ 39]*x13*x21
        +coeff[ 40]*x11*x22*x31*x41
        +coeff[ 41]        *x31*x41*x51
        +coeff[ 42]            *x42*x51
        +coeff[ 43]*x11        *x42*x51
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 44]        *x34
        +coeff[ 45]*x12*x21*x31
        +coeff[ 46]            *x41*x51
        +coeff[ 47]    *x22*x33
        +coeff[ 48]*x12    *x32*x41
        +coeff[ 49]*x11*x21*x32*x41
        +coeff[ 50]*x12    *x31*x42
        +coeff[ 51]        *x32*x44
        +coeff[ 52]        *x32    *x51
    ;
    v_y00_l5p77                               =v_y00_l5p77
        +coeff[ 53]*x13*x21*x31
        +coeff[ 54]*x12            *x51
        +coeff[ 55]*x11*x21        *x51
        +coeff[ 56]    *x22        *x51
        +coeff[ 57]*x13    *x33
        +coeff[ 58]*x13    *x32*x41
        +coeff[ 59]*x12*x21*x32*x41
        ;

    return v_y00_l5p77                               ;
}
}
