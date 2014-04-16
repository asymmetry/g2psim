#include "Fwd_l5p77_484816R15.h"

namespace S484816R15
{
float x_l5p77_den                             (float *x,int m){
    int ncoeff= 60;
    float avdat= -0.5101204E+01;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
        -0.81959576E-03, 0.14616318E-01,-0.98075286E-01,-0.53911959E-02,
        -0.35031866E-01, 0.77880039E-04, 0.61130745E-03,-0.23459137E-04,
         0.43797037E-02, 0.30268552E-02,-0.29669735E-02,-0.79280958E-02,
         0.48737573E-02,-0.12855270E-01,-0.87778335E-02, 0.14148776E-01,
         0.15808228E-04,-0.63506799E-03, 0.11159170E-02,-0.12405650E-02,
        -0.12748021E-03,-0.15586975E-02,-0.64287642E-02, 0.95892251E-02,
         0.13710288E-01, 0.44949524E-03,-0.10495341E-02, 0.16915092E-02,
         0.83082626E-02,-0.37247715E-04, 0.12938236E-03, 0.21778829E-03,
        -0.68257464E-03, 0.43610873E-03,-0.12840992E-02, 0.17239571E-02,
         0.24081802E-02,-0.17574205E-02, 0.17097516E-02, 0.28689683E-02,
        -0.37466320E-02, 0.41881847E-03,-0.12611572E-03,-0.42280144E-03,
         0.33193556E-03, 0.25664747E-03, 0.79829112E-03,-0.52885065E-03,
        -0.41877739E-02, 0.10065008E-03,-0.11722506E-03,-0.10694713E-03,
        -0.16093919E-03, 0.20047436E-03, 0.30508311E-03,-0.52931212E-03,
        -0.13133122E-02, 0.13737759E-03,-0.63821906E-03,-0.56684104E-03,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_x_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]    *x21        *x51
        +coeff[  4]    *x21    *x41
        +coeff[  5]*x13        *x41
        +coeff[  6]*x12*x21*x31
        +coeff[  7]    *x21*x31    *x52
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[  8]    *x23*x31
        +coeff[  9]            *x41
        +coeff[ 10]*x11    *x31
        +coeff[ 11]*x11        *x41
        +coeff[ 12]    *x22
        +coeff[ 13]    *x21*x31
        +coeff[ 14]    *x23
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x13*x22
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 17]    *x23*x31*x41
        +coeff[ 18]        *x31
        +coeff[ 19]*x11            *x51
        +coeff[ 20]*x13
        +coeff[ 21]*x12*x21
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]    *x23    *x41
        +coeff[ 25]*x11*x21
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 26]            *x42
        +coeff[ 27]    *x21*x32
        +coeff[ 28]*x11*x22    *x41
        +coeff[ 29]*x13        *x42
        +coeff[ 30]                *x51
        +coeff[ 31]*x12
        +coeff[ 32]        *x31*x41
        +coeff[ 33]    *x21        *x52
        +coeff[ 34]*x11*x21    *x41
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]*x11        *x42
        +coeff[ 37]    *x22    *x41
        +coeff[ 38]*x12*x21    *x41
        +coeff[ 39]*x11*x22*x31
        +coeff[ 40]    *x21    *x43
        +coeff[ 41]    *x23*x31*x42
        +coeff[ 42]        *x32
        +coeff[ 43]*x11*x21*x31
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 44]*x11    *x32
        +coeff[ 45]    *x22        *x51
        +coeff[ 46]    *x21    *x41*x51
        +coeff[ 47]    *x22*x31
        +coeff[ 48]    *x21*x31*x42
        +coeff[ 49]*x12*x21*x31    *x51
        +coeff[ 50]*x12*x21*x32*x41
        +coeff[ 51]    *x23*x32*x41
        +coeff[ 52]*x12        *x41
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 53]*x11        *x41*x51
        +coeff[ 54]    *x21*x31    *x51
        +coeff[ 55]    *x21    *x42*x51
        +coeff[ 56]    *x21*x32*x41
        +coeff[ 57]*x13    *x31    *x51
        +coeff[ 58]    *x23    *x41*x51
        +coeff[ 59]*x12*x21*x31*x41*x51
        ;

    return v_x_l5p77_den                             ;
}
float t_l5p77_den                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.3631719E+01;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
        -0.43407299E-01,-0.11986437E+00, 0.60004592E+00,-0.52316073E-01,
         0.13283540E+00, 0.23689860E+00, 0.35526983E-01, 0.84690526E-01,
         0.11431852E+00,-0.82707144E-02,-0.19466067E-01, 0.10390907E-04,
        -0.40151677E-02,-0.20679504E-01, 0.61280057E-02, 0.88862635E-01,
         0.55775896E-01, 0.39016146E-01, 0.90872636E-02, 0.29496202E-01,
         0.37783079E-01,-0.91575772E-01,-0.78446150E-03,-0.78376746E-02,
         0.20080270E-01,-0.79515977E-02,-0.58639422E-01,-0.63559286E-01,
         0.14300658E-02,-0.67807171E-02, 0.39902152E-02, 0.10101167E-01,
         0.12350409E-01,-0.90109063E-02, 0.15868146E-01,-0.16610034E-01,
         0.10909669E-01,-0.12265643E-01,-0.58163695E-01,-0.37707224E-01,
        -0.27788588E-03, 0.39537679E-02,-0.10168036E-01,-0.40800693E-02,
         0.43022251E-02,-0.32907110E-02,-0.15767174E-01,-0.17550176E-01,
         0.12533536E-01, 0.20041151E-01,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_t_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]*x11*x21
        +coeff[  4]    *x22
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x21        *x51
        +coeff[  7]    *x23
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[  8]    *x22    *x41
        +coeff[  9]*x12*x21*x31
        +coeff[ 10]    *x23*x31
        +coeff[ 11]*x13        *x41
        +coeff[ 12]*x13*x22*x31
        +coeff[ 13]            *x41
        +coeff[ 14]                *x51
        +coeff[ 15]    *x21*x31
        +coeff[ 16]*x11        *x41
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 17]            *x42
        +coeff[ 18]*x11            *x51
        +coeff[ 19]*x11*x22
        +coeff[ 20]    *x22*x31
        +coeff[ 21]    *x21    *x42
        +coeff[ 22]*x12*x21*x31*x41
        +coeff[ 23]        *x31
        +coeff[ 24]*x11    *x31
        +coeff[ 25]            *x41*x51
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 26]    *x21*x31*x41
        +coeff[ 27]    *x23    *x41
        +coeff[ 28]*x12*x23
        +coeff[ 29]*x13*x22    *x41
        +coeff[ 30]*x12
        +coeff[ 31]        *x31*x41
        +coeff[ 32]*x12*x21
        +coeff[ 33]    *x21*x32
        +coeff[ 34]*x11*x21    *x41
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 35]*x11        *x42
        +coeff[ 36]    *x22        *x51
        +coeff[ 37]*x12*x21    *x41
        +coeff[ 38]*x11*x22    *x41
        +coeff[ 39]    *x22    *x42
        +coeff[ 40]*x13    *x31*x41
        +coeff[ 41]*x11*x21*x31
        +coeff[ 42]*x11    *x31*x41
        +coeff[ 43]    *x21    *x41*x51
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 44]            *x42*x51
        +coeff[ 45]    *x21        *x52
        +coeff[ 46]*x11*x22*x31
        +coeff[ 47]    *x22*x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x12*x21*x31*x42
        ;

    return v_t_l5p77_den                             ;
}
float y_l5p77_den                             (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1376437E-01;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.10560826E-01, 0.13908137E+00,-0.51848944E-02, 0.13416602E-01,
        -0.21691246E-01, 0.75193397E-02, 0.30727601E-01, 0.16652556E-01,
         0.13431000E-01,-0.21796111E-01,-0.11705220E-02,-0.12248415E-01,
        -0.89207133E-02, 0.21676805E-02,-0.78807473E-02,-0.94183721E-02,
         0.22132862E-02,-0.16512698E-02, 0.15709029E-02,-0.11353276E-02,
        -0.36742168E-02, 0.21085655E-02,-0.59659878E-03, 0.18513274E-03,
         0.19752274E-02,-0.13278343E-02, 0.10743908E-01, 0.82834531E-02,
         0.26762311E-03, 0.11033304E-02,-0.12029695E-02, 0.76899258E-03,
         0.47655329E-02,-0.52289017E-02, 0.36043744E-02,-0.42709871E-02,
        -0.13780061E-02,-0.51026436E-03,-0.39583287E-03,-0.21474050E-03,
         0.15394699E-02, 0.20489910E-03,-0.84689265E-03, 0.67307026E-03,
         0.13779759E-03, 0.12051559E-02,-0.46736651E-03, 0.56629663E-03,
         0.20413914E-03, 0.22659956E-03, 0.16989287E-03,-0.94518582E-04,
         0.84624902E-04, 0.11183446E-03, 0.59249013E-03,-0.38260716E-03,
        -0.28353180E-02, 0.35530704E-03,-0.53520501E-03,-0.21010765E-03,
         0.99375215E-03, 0.19037200E-02, 0.18513374E-03,-0.32838227E-03,
        -0.41745941E-04,-0.12507378E-03, 0.10070985E-03,-0.77944016E-04,
        -0.26275555E-03,-0.10035514E-03, 0.62249077E-04, 0.26787040E-03,
        -0.37359309E-03,-0.40231788E-03,-0.20860706E-03,-0.10055188E-03,
        -0.20010686E-02,-0.33252864E-03, 0.86953858E-03,-0.56402688E-03,
        -0.24150852E-03,-0.17097710E-03,-0.13705708E-03, 0.73474913E-03,
         0.58138056E-03,-0.33180724E-03,-0.31729511E-03,-0.33272477E-02,
        -0.22818735E-02,-0.50161918E-04, 0.18187467E-03,-0.66261142E-04,
         0.35807543E-03, 0.55606666E-04,-0.14320156E-03, 0.87851600E-04,
         0.66420872E-03,
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
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x22
        +coeff[  7]            *x41*x51
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[  8]*x11*x21
        +coeff[  9]    *x22    *x41
        +coeff[ 10]*x11
        +coeff[ 11]            *x42
        +coeff[ 12]        *x31*x41
        +coeff[ 13]        *x31    *x51
        +coeff[ 14]    *x22*x31
        +coeff[ 15]*x11*x21    *x41
        +coeff[ 16]    *x21*x31
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 17]        *x32
        +coeff[ 18]*x12
        +coeff[ 19]                *x52
        +coeff[ 20]*x11*x21*x31
        +coeff[ 21]    *x22        *x51
        +coeff[ 22]    *x21        *x51
        +coeff[ 23]            *x43
        +coeff[ 24]    *x23
        +coeff[ 25]            *x41*x52
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 26]    *x22    *x42
        +coeff[ 27]    *x22*x31*x41
        +coeff[ 28]            *x42*x51
        +coeff[ 29]*x11*x22
        +coeff[ 30]*x12        *x41
        +coeff[ 31]*x11*x21        *x51
        +coeff[ 32]*x11*x21    *x42
        +coeff[ 33]    *x24
        +coeff[ 34]*x11*x21*x31*x41
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 35]*x11*x23
        +coeff[ 36]*x12*x22
        +coeff[ 37]    *x21    *x41*x51
        +coeff[ 38]*x12    *x31
        +coeff[ 39]        *x31    *x52
        +coeff[ 40]    *x22*x32
        +coeff[ 41]                *x53
        +coeff[ 42]    *x22    *x41*x51
        +coeff[ 43]*x11*x21*x32
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 44]*x11    *x31
        +coeff[ 45]        *x31*x42
        +coeff[ 46]    *x21*x31*x41
        +coeff[ 47]        *x32*x41
        +coeff[ 48]*x11        *x41*x51
        +coeff[ 49]*x12*x21
        +coeff[ 50]    *x21        *x52
        +coeff[ 51]*x11        *x43
        +coeff[ 52]*x12            *x51
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 53]    *x21    *x44
        +coeff[ 54]*x12        *x42
        +coeff[ 55]    *x22*x31    *x51
        +coeff[ 56]    *x22    *x43
        +coeff[ 57]*x12    *x31*x41
        +coeff[ 58]*x11*x21    *x41*x51
        +coeff[ 59]*x11*x21*x31    *x51
        +coeff[ 60]    *x22    *x42*x51
        +coeff[ 61]*x11*x23    *x41
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 62]*x11        *x41
        +coeff[ 63]    *x21    *x42
        +coeff[ 64]*x11            *x51
        +coeff[ 65]    *x21*x32
        +coeff[ 66]        *x33
        +coeff[ 67]*x11    *x31*x41
        +coeff[ 68]    *x21    *x43
        +coeff[ 69]    *x21*x31*x42
        +coeff[ 70]*x11    *x31    *x51
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 71]            *x43*x51
        +coeff[ 72]    *x21    *x42*x51
        +coeff[ 73]            *x45
        +coeff[ 74]    *x21*x31*x41*x51
        +coeff[ 75]    *x23        *x51
        +coeff[ 76]    *x22*x31*x42
        +coeff[ 77]    *x22        *x52
        +coeff[ 78]    *x24    *x41
        +coeff[ 79]    *x22*x32*x41
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 80]*x11*x21*x31*x42
        +coeff[ 81]*x13*x21
        +coeff[ 82]*x11*x21        *x52
        +coeff[ 83]*x11*x23*x31
        +coeff[ 84]*x12*x22    *x41
        +coeff[ 85]*x11*x24    *x41
        +coeff[ 86]*x11*x22    *x44
        +coeff[ 87]    *x24    *x44
        +coeff[ 88]    *x24*x31*x43
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 89]        *x32    *x51
        +coeff[ 90]        *x31*x44
        +coeff[ 91]*x11        *x42*x51
        +coeff[ 92]        *x32*x43
        +coeff[ 93]*x11*x22        *x51
        +coeff[ 94]        *x31*x43*x51
        +coeff[ 95]            *x41*x53
        +coeff[ 96]    *x24*x31
        ;

    return v_y_l5p77_den                             ;
}
float p_l5p77_den                             (float *x,int m){
    int ncoeff= 70;
    float avdat= -0.3827681E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 71]={
         0.58807656E-02, 0.11644445E-02,-0.19197173E-01,-0.78174062E-01,
         0.46046120E-02,-0.66275080E-02,-0.45138337E-02,-0.26437987E-01,
         0.38391920E-02, 0.25411530E-02, 0.17036267E-01,-0.31332595E-02,
        -0.22210025E-02, 0.45270566E-02,-0.63175336E-02, 0.53854579E-04,
         0.25571776E-02,-0.64527891E-02, 0.45086895E-02, 0.68167220E-02,
         0.26594482E-02,-0.35743276E-02, 0.70858496E-03,-0.27509364E-02,
        -0.54175587E-03, 0.14282792E-02,-0.10575869E-02, 0.14126075E-02,
        -0.18266928E-02,-0.11553977E-01, 0.30745294E-04,-0.36014864E-03,
        -0.13498738E-02,-0.45381009E-03,-0.92078251E-03,-0.56274799E-02,
         0.90351253E-03,-0.43216007E-03,-0.27145648E-02, 0.48746078E-04,
        -0.56821934E-03, 0.19627299E-03, 0.42125303E-03, 0.34949949E-03,
        -0.34331020E-02, 0.44930371E-03,-0.95206685E-03, 0.29303227E-02,
        -0.10077470E-02,-0.15851908E-02,-0.14819382E-02,-0.74724900E-03,
        -0.70803246E-03,-0.27503995E-02, 0.42978927E-03, 0.23518944E-03,
        -0.22076585E-03,-0.16282115E-02,-0.92764363E-04,-0.16602644E-03,
         0.32190792E-03,-0.22031290E-03, 0.16357728E-03, 0.16495051E-02,
         0.12551008E-03,-0.29729647E-03,-0.19277162E-02, 0.17213543E-02,
         0.10545523E-02,-0.31388173E-03,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21*x31
        +coeff[  7]    *x21    *x41
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[  8]            *x42
        +coeff[  9]    *x21        *x51
        +coeff[ 10]            *x41*x51
        +coeff[ 11]*x11*x21
        +coeff[ 12]                *x52
        +coeff[ 13]*x11        *x41
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]        *x33    *x51
        +coeff[ 16]        *x31    *x51
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]    *x22        *x51
        +coeff[ 19]*x11*x21    *x41
        +coeff[ 20]        *x31*x41
        +coeff[ 21]    *x23
        +coeff[ 22]*x11    *x31
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]*x11            *x51
        +coeff[ 25]    *x21    *x41*x51
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 26]*x11*x22
        +coeff[ 27]*x11*x21*x31
        +coeff[ 28]*x11        *x42
        +coeff[ 29]    *x22    *x42
        +coeff[ 30]            *x45
        +coeff[ 31]*x12
        +coeff[ 32]            *x41*x52
        +coeff[ 33]    *x22*x32
        +coeff[ 34]*x11    *x31*x41
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 35]    *x22*x31*x41
        +coeff[ 36]*x11*x21        *x51
        +coeff[ 37]*x11        *x41*x51
        +coeff[ 38]*x11*x21    *x42
        +coeff[ 39]        *x31*x44
        +coeff[ 40]        *x31*x43*x51
        +coeff[ 41]*x11
        +coeff[ 42]        *x32
        +coeff[ 43]    *x22*x31
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 44]            *x43
        +coeff[ 45]    *x21*x31    *x51
        +coeff[ 46]    *x24
        +coeff[ 47]    *x21    *x43
        +coeff[ 48]    *x22    *x41*x51
        +coeff[ 49]*x11*x22    *x41
        +coeff[ 50]*x11*x21*x31*x41
        +coeff[ 51]*x11*x21    *x41*x51
        +coeff[ 52]*x12*x21    *x41
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 53]*x11*x23    *x41
        +coeff[ 54]    *x23*x31*x42
        +coeff[ 55]*x12            *x53
        +coeff[ 56]    *x21*x32
        +coeff[ 57]        *x31*x42
        +coeff[ 58]        *x32    *x51
        +coeff[ 59]    *x21        *x52
        +coeff[ 60]    *x23    *x41
        +coeff[ 61]        *x31    *x52
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 62]                *x53
        +coeff[ 63]    *x21*x31*x42
        +coeff[ 64]*x11            *x52
        +coeff[ 65]    *x22        *x52
        +coeff[ 66]    *x24    *x41
        +coeff[ 67]    *x23    *x42
        +coeff[ 68]    *x22    *x42*x51
        +coeff[ 69]        *x34*x41
        ;

    return v_p_l5p77_den                             ;
}
float l_l5p77_den                             (float *x,int m){
    int ncoeff= 67;
    float avdat= -0.6224317E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 68]={
         0.46872459E-02, 0.18921329E+00,-0.28190628E-01,-0.15882835E-01,
         0.68502791E-01, 0.10354555E-01, 0.15272221E-01,-0.87174373E-02,
         0.30470602E-03, 0.24914751E-01,-0.62741116E-02, 0.57436265E-02,
         0.17474033E-01,-0.26311491E-01, 0.12468064E-01, 0.14335813E-02,
         0.20085818E-02, 0.23628864E-02,-0.36281324E-02,-0.18904585E-01,
         0.30252619E-02,-0.26920205E-01,-0.16310599E-01, 0.10717124E-02,
        -0.33054508E-02,-0.46356018E-02,-0.55758846E-02, 0.69865373E-05,
        -0.23960799E-03,-0.31088435E-03,-0.65991428E-03,-0.87311305E-03,
        -0.33460113E-02, 0.73680207E-02,-0.32594688E-02, 0.87439606E-03,
         0.14256408E-03,-0.72910730E-03, 0.27957026E-03,-0.68325264E-03,
        -0.70240931E-03,-0.13299590E-02, 0.58968458E-03,-0.64532016E-03,
         0.13699167E-02,-0.12537552E-02, 0.21618085E-02, 0.77536628E-02,
        -0.11742075E-02, 0.14497850E-03,-0.40496129E-03, 0.71940722E-03,
        -0.11600467E-03,-0.71296201E-03, 0.83114888E-03,-0.39269417E-03,
         0.35975574E-03, 0.20291936E-03, 0.10715252E-02, 0.25271955E-02,
         0.12179526E-02,-0.53553039E-03, 0.12821495E-02,-0.22490486E-02,
        -0.24001782E-03, 0.99723320E-03, 0.22847746E-02,
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
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_l_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]*x11
        +coeff[  3]    *x22
        +coeff[  4]    *x21    *x41
        +coeff[  5]    *x21        *x51
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x23*x31
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[  8]    *x21*x33
        +coeff[  9]    *x21*x31
        +coeff[ 10]            *x42
        +coeff[ 11]*x11    *x31
        +coeff[ 12]    *x23
        +coeff[ 13]    *x21    *x42
        +coeff[ 14]*x11*x22
        +coeff[ 15]    *x23*x31*x41
        +coeff[ 16]            *x41*x51
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 17]*x11            *x51
        +coeff[ 18]    *x22    *x41
        +coeff[ 19]    *x21*x31*x41
        +coeff[ 20]*x12*x21
        +coeff[ 21]    *x23    *x41
        +coeff[ 22]*x11*x22    *x41
        +coeff[ 23]*x11*x21
        +coeff[ 24]    *x21*x32
        +coeff[ 25]*x11        *x42
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 26]*x11*x22*x31
        +coeff[ 27]*x11    *x33*x41
        +coeff[ 28]        *x31
        +coeff[ 29]                *x51
        +coeff[ 30]*x12
        +coeff[ 31]    *x21        *x52
        +coeff[ 32]*x11    *x31*x41
        +coeff[ 33]    *x21    *x43
        +coeff[ 34]*x12*x21    *x41
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 35]    *x23*x31*x42
        +coeff[ 36]        *x32
        +coeff[ 37]        *x31*x41
        +coeff[ 38]        *x31    *x51
        +coeff[ 39]    *x22*x31
        +coeff[ 40]    *x22        *x51
        +coeff[ 41]    *x21    *x41*x51
        +coeff[ 42]            *x42*x51
        +coeff[ 43]*x11    *x32
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 44]*x11*x21    *x41
        +coeff[ 45]    *x24
        +coeff[ 46]    *x22    *x42
        +coeff[ 47]    *x21*x31*x42
        +coeff[ 48]*x12*x21*x31
        +coeff[ 49]    *x23*x31    *x51
        +coeff[ 50]*x11*x23*x31
        +coeff[ 51]    *x23*x32*x41
        +coeff[ 52]                *x52
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 53]    *x21*x31    *x51
        +coeff[ 54]*x11*x21*x31
        +coeff[ 55]*x11        *x41*x51
        +coeff[ 56]*x12        *x41
        +coeff[ 57]*x13
        +coeff[ 58]    *x22*x31*x41
        +coeff[ 59]    *x21*x32*x41
        +coeff[ 60]    *x21    *x42*x51
        +coeff[ 61]*x11*x23
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 62]    *x24    *x41
        +coeff[ 63]    *x21    *x44
        +coeff[ 64]*x11    *x31    *x53
        +coeff[ 65]    *x23*x31*x41*x51
        +coeff[ 66]    *x23    *x43*x51
        ;

    return v_l_l5p77_den                             ;
}
float x_l5p77_dex                             (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.8277598E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
        -0.16539904E-02,-0.51949944E-01, 0.13304844E+00, 0.26510850E+00,
         0.34274541E-01, 0.10771694E+00,-0.91737820E-04,-0.24439530E-02,
        -0.85843104E-03,-0.14473050E-01,-0.88240001E-02, 0.90227034E-02,
         0.24168666E-01, 0.40417608E-01, 0.27572906E-01,-0.45233313E-01,
         0.57089020E-03,-0.85752872E-04,-0.32340235E-02,-0.34940941E-02,
        -0.56396304E-02,-0.21698342E-02, 0.46627200E-02, 0.18500468E-01,
         0.11312559E-01,-0.28026065E-01,-0.41185569E-01, 0.54270062E-02,
         0.43398314E-02,-0.47124336E-02,-0.25727561E-01, 0.22508464E-03,
         0.69474679E-03, 0.12551142E-02, 0.35800149E-02,-0.51066875E-02,
        -0.69898884E-02, 0.13759040E-02,-0.57043075E-02,-0.87979967E-02,
         0.58218115E-03, 0.78067940E-03, 0.95581706E-03,-0.89717482E-03,
        -0.93666825E-03, 0.20170223E-02, 0.69342307E-02,-0.25026514E-02,
         0.55922917E-02, 0.42985403E-02, 0.68716295E-02, 0.49708541E-02,
        -0.24165375E-03, 0.24351590E-03, 0.46054964E-03, 0.11432369E-02,
        -0.78090600E-03, 0.85848785E-03, 0.11498345E-02,-0.95927465E-03,
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
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_x_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]    *x21
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x21    *x41
        +coeff[  6]*x13        *x41
        +coeff[  7]*x12*x21*x31
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[  8]    *x21*x31    *x52
        +coeff[  9]    *x23*x31
        +coeff[ 10]            *x41
        +coeff[ 11]*x11    *x31
        +coeff[ 12]*x11        *x41
        +coeff[ 13]    *x21*x31
        +coeff[ 14]    *x23
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x13*x22
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 17]    *x23*x31*x41
        +coeff[ 18]        *x31
        +coeff[ 19]                *x52
        +coeff[ 20]*x11*x21
        +coeff[ 21]            *x42
        +coeff[ 22]*x12*x21
        +coeff[ 23]*x11*x22
        +coeff[ 24]    *x22    *x41
        +coeff[ 25]    *x21*x31*x41
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 26]    *x23    *x41
        +coeff[ 27]    *x21    *x41*x51
        +coeff[ 28]    *x22*x31
        +coeff[ 29]    *x21*x32
        +coeff[ 30]*x11*x22    *x41
        +coeff[ 31]*x13        *x42
        +coeff[ 32]*x11            *x51
        +coeff[ 33]            *x41*x51
        +coeff[ 34]*x11*x21    *x41
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]*x11        *x42
        +coeff[ 37]    *x21*x31    *x51
        +coeff[ 38]*x12*x21    *x41
        +coeff[ 39]*x11*x22*x31
        +coeff[ 40]*x13*x21*x31
        +coeff[ 41]    *x22
        +coeff[ 42]        *x31*x41
        +coeff[ 43]    *x21        *x52
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 44]*x11    *x32
        +coeff[ 45]    *x23        *x51
        +coeff[ 46]    *x21    *x43
        +coeff[ 47]    *x23    *x42
        +coeff[ 48]*x12*x21*x31*x42
        +coeff[ 49]    *x21*x31*x42*x52
        +coeff[ 50]    *x23*x31*x42
        +coeff[ 51]*x14*x21*x32*x41
        +coeff[ 52]*x12
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 53]        *x32
        +coeff[ 54]*x11*x21        *x51
        +coeff[ 55]*x11*x21*x31
        +coeff[ 56]            *x42*x51
        +coeff[ 57]*x11*x22        *x51
        +coeff[ 58]*x11*x23
        +coeff[ 59]    *x22    *x41*x51
        ;

    return v_x_l5p77_dex                             ;
}
float t_l5p77_dex                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.5904654E+00;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.36752112E-07, 0.12164116E-01,-0.82909934E-01, 0.23661680E-02,
         0.30996365E-01,-0.30015010E-01,-0.64140265E-02, 0.47539760E-03,
         0.37444334E-02,-0.26743834E-04,-0.23946529E-02,-0.11197585E-01,
        -0.64525837E-02,-0.76424237E-02, 0.12294588E-01, 0.40806137E-03,
         0.97796042E-03, 0.20403126E-02,-0.17597103E-02,-0.77117042E-03,
        -0.13428333E-02,-0.53343442E-02,-0.35741548E-02, 0.69696899E-02,
        -0.88616624E-04, 0.11605892E-01, 0.12983896E-02,-0.12956207E-02,
         0.11165227E-02, 0.18330553E-03, 0.68765241E-02,-0.29092353E-04,
        -0.13246357E-03,-0.66145294E-03, 0.17246031E-02, 0.58553036E-03,
        -0.45535158E-03,-0.10876654E-02, 0.39935362E-03, 0.52175281E-03,
         0.21900057E-02, 0.14969263E-02,-0.67498295E-05, 0.18584823E-02,
        -0.41122164E-03,-0.16958050E-03, 0.36454134E-03,-0.23998028E-03,
         0.22396354E-03,-0.11125259E-02,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_t_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x21        *x51
        +coeff[  7]*x12*x21*x31
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[  8]    *x23*x31
        +coeff[  9]*x13        *x41
        +coeff[ 10]*x11    *x31
        +coeff[ 11]    *x21*x31
        +coeff[ 12]*x11        *x41
        +coeff[ 13]    *x23
        +coeff[ 14]    *x21    *x42
        +coeff[ 15]*x13*x22
        +coeff[ 16]        *x31
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 17]    *x22
        +coeff[ 18]            *x42
        +coeff[ 19]*x11            *x51
        +coeff[ 20]                *x52
        +coeff[ 21]*x11*x22
        +coeff[ 22]    *x22    *x41
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]*x13*x21
        +coeff[ 25]    *x23    *x41
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 26]*x11*x21
        +coeff[ 27]*x12*x21
        +coeff[ 28]    *x21*x32
        +coeff[ 29]        *x31*x41*x51
        +coeff[ 30]*x11*x22    *x41
        +coeff[ 31]*x11*x23    *x41
        +coeff[ 32]        *x31    *x51
        +coeff[ 33]    *x22*x31
        +coeff[ 34]*x11        *x42
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 35]    *x22        *x51
        +coeff[ 36]    *x21    *x41*x51
        +coeff[ 37]            *x42*x51
        +coeff[ 38]    *x21        *x52
        +coeff[ 39]            *x41*x52
        +coeff[ 40]*x11*x22*x31
        +coeff[ 41]*x12*x21    *x41
        +coeff[ 42]            *x43*x51
        +coeff[ 43]*x13    *x31*x41
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 44]*x11*x21    *x43
        +coeff[ 45]        *x31*x41
        +coeff[ 46]            *x41*x51
        +coeff[ 47]*x11*x21*x31
        +coeff[ 48]*x11    *x32
        +coeff[ 49]*x11*x21    *x41
        ;

    return v_t_l5p77_dex                             ;
}
float y_l5p77_dex                             (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1534991E-01;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.62709488E-02, 0.76020487E-01,-0.67677116E-02,-0.14770331E-01,
        -0.26289424E-01,-0.12165310E-01,-0.44100504E-01,-0.90986686E-02,
         0.37428480E-01,-0.43318984E-02, 0.84040277E-02, 0.55203311E-01,
         0.15371344E-01, 0.53282049E-02, 0.78260479E-02,-0.14293845E-01,
        -0.56517674E-02,-0.40239174E-01,-0.48015369E-02,-0.10633140E-01,
         0.10535141E-01,-0.12414184E-02,-0.17367488E-02,-0.53540883E-02,
         0.16492626E-02,-0.13527067E-02,-0.66725286E-02,-0.64228880E-02,
        -0.42649936E-02, 0.29171067E-02, 0.11885333E-02,-0.49987053E-02,
        -0.19103317E-02,-0.31179807E-02,-0.95467735E-02,-0.45566331E-02,
        -0.54315040E-02,-0.11531046E-02,-0.53710176E-03,-0.25511959E-02,
        -0.14830417E-02, 0.68600085E-02,-0.15822313E-02,-0.21957958E-02,
         0.50642039E-02,-0.69096580E-03, 0.76267205E-03,-0.29695868E-02,
        -0.10264616E-02,-0.20660926E-02,-0.14164847E-02, 0.44068604E-03,
         0.26201655E-03, 0.50189364E-03,-0.49082830E-03, 0.53514290E-03,
         0.32709502E-02, 0.30788337E-02, 0.99324761E-03, 0.12836908E-02,
         0.45486912E-03, 0.29374133E-02,-0.14731056E-02, 0.64594834E-02,
        -0.12638023E-02,-0.88943302E-03,-0.91926742E-03,-0.35268455E-02,
         0.42881328E-02, 0.44958820E-02, 0.34968802E-03,-0.38637492E-03,
        -0.23112196E-03,-0.28112289E-03, 0.23085192E-03, 0.10754794E-02,
         0.12474755E-02, 0.12279433E-02, 0.10760819E-02, 0.86160272E-03,
        -0.18707296E-02, 0.78211329E-03, 0.15672154E-03,-0.10377495E-02,
         0.83952682E-03, 0.60038135E-03,-0.14677554E-02, 0.47283288E-03,
         0.52632991E-03, 0.11642950E-02, 0.26446269E-02,-0.43672373E-03,
        -0.13550030E-02, 0.21159581E-02, 0.18599241E-02, 0.10336188E-03,
         0.88820067E-04,
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
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]            *x42
        +coeff[  6]    *x21    *x41
        +coeff[  7]        *x31*x41
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[  8]    *x22
        +coeff[  9]    *x21*x31
        +coeff[ 10]*x11        *x41
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]        *x31    *x51
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]                *x52
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]*x11        *x42
        +coeff[ 19]    *x22*x31
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]*x11
        +coeff[ 22]        *x32
        +coeff[ 23]            *x43
        +coeff[ 24]*x12
        +coeff[ 25]*x11            *x51
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 26]    *x21*x31*x41
        +coeff[ 27]    *x23
        +coeff[ 28]            *x41*x52
        +coeff[ 29]*x11*x21        *x51
        +coeff[ 30]*x11    *x31
        +coeff[ 31]*x11*x21    *x41
        +coeff[ 32]    *x21    *x41*x51
        +coeff[ 33]*x11*x21*x31
        +coeff[ 34]    *x24
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 35]    *x22    *x41*x51
        +coeff[ 36]*x11*x23
        +coeff[ 37]        *x31*x42
        +coeff[ 38]            *x42*x51
        +coeff[ 39]*x11    *x31*x41
        +coeff[ 40]*x11*x22
        +coeff[ 41]    *x21    *x43
        +coeff[ 42]*x12        *x41
        +coeff[ 43]    *x22    *x42
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 44]    *x21*x31*x42
        +coeff[ 45]        *x31    *x52
        +coeff[ 46]                *x53
        +coeff[ 47]*x11*x22    *x41
        +coeff[ 48]*x12*x21    *x41
        +coeff[ 49]*x11*x21    *x41*x51
        +coeff[ 50]        *x31*x41*x51
        +coeff[ 51]    *x21*x31    *x51
        +coeff[ 52]*x11        *x41*x51
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 53]*x12*x21
        +coeff[ 54]*x12    *x31
        +coeff[ 55]    *x21        *x52
        +coeff[ 56]    *x22*x31*x41
        +coeff[ 57]*x11*x21    *x42
        +coeff[ 58]    *x23*x31
        +coeff[ 59]    *x22*x32
        +coeff[ 60]*x12            *x51
        +coeff[ 61]*x11*x21*x31*x41
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 62]    *x22*x31    *x51
        +coeff[ 63]    *x23    *x42
        +coeff[ 64]*x12*x22
        +coeff[ 65]*x11*x21*x31    *x51
        +coeff[ 66]    *x22        *x52
        +coeff[ 67]    *x24    *x41
        +coeff[ 68]*x11*x22    *x42
        +coeff[ 69]    *x22    *x42*x51
        +coeff[ 70]        *x32*x41
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 71]    *x21*x32
        +coeff[ 72]*x11    *x32
        +coeff[ 73]        *x32    *x51
        +coeff[ 74]*x11    *x31    *x51
        +coeff[ 75]*x11        *x43
        +coeff[ 76]            *x43*x51
        +coeff[ 77]    *x23    *x41
        +coeff[ 78]    *x21*x32*x41
        +coeff[ 79]*x11    *x31*x42
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 80]    *x21    *x42*x51
        +coeff[ 81]        *x31*x42*x51
        +coeff[ 82]*x11            *x52
        +coeff[ 83]    *x21*x31*x41*x51
        +coeff[ 84]*x12        *x42
        +coeff[ 85]*x11*x21*x32
        +coeff[ 86]    *x22    *x43
        +coeff[ 87]*x12    *x31*x41
        +coeff[ 88]    *x21    *x41*x52
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 89]*x11*x21    *x43
        +coeff[ 90]    *x23*x31*x41
        +coeff[ 91]*x11*x21        *x52
        +coeff[ 92]*x11*x23    *x41
        +coeff[ 93]*x11*x22*x31*x41
        +coeff[ 94]    *x22*x31*x41*x51
        +coeff[ 95]        *x33
        +coeff[ 96]*x13
        ;

    return v_y_l5p77_dex                             ;
}
float p_l5p77_dex                             (float *x,int m){
    int ncoeff= 90;
    float avdat=  0.1597260E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.21376721E-03,-0.88761683E-03,-0.59889657E-02,-0.32703893E-02,
        -0.31230147E-02, 0.43510804E-02,-0.99780422E-03,-0.10603494E-01,
        -0.93305937E-03,-0.11843165E-02, 0.12085094E-02, 0.15223607E-02,
         0.12192381E-01, 0.17734681E-02,-0.15378136E-02, 0.17763462E-02,
        -0.74469163E-02,-0.27533120E-02, 0.29898391E-02,-0.14911565E-02,
        -0.54484105E-03,-0.40182094E-04,-0.10054874E-02, 0.70514536E-03,
         0.17376764E-03, 0.25958123E-03, 0.22104379E-04,-0.12394720E-03,
         0.21916740E-03,-0.14743642E-02,-0.92975155E-03,-0.12714370E-02,
         0.19107462E-03,-0.55240106E-03,-0.17490018E-02,-0.14896087E-02,
        -0.17955904E-03,-0.18632027E-02,-0.33724998E-03,-0.43848285E-03,
        -0.29004796E-03,-0.16079142E-02, 0.32197681E-03,-0.26938517E-03,
         0.10219772E-02,-0.46075793E-03, 0.20977945E-03,-0.68374793E-03,
        -0.21894311E-03,-0.83206582E-03,-0.55362576E-04, 0.44954999E-03,
        -0.14977117E-04, 0.28908550E-03, 0.17079449E-03, 0.14184955E-03,
         0.15161523E-03, 0.13410375E-02,-0.36381034E-03,-0.27919863E-03,
        -0.45114162E-03,-0.20201041E-03, 0.26371866E-03, 0.20022362E-04,
         0.18121893E-03,-0.76306482E-04,-0.14775210E-03,-0.45637797E-04,
        -0.29219719E-03, 0.81499969E-03, 0.40750179E-04,-0.46298673E-03,
        -0.51772680E-04,-0.11382720E-03, 0.13318447E-02, 0.12784991E-03,
         0.16495833E-03, 0.20429394E-03,-0.19209793E-03, 0.38852877E-03,
        -0.86569664E-04,-0.89764908E-04,-0.31095603E-03, 0.28283274E-03,
        -0.20953032E-03, 0.62773551E-03, 0.15828678E-03, 0.33194377E-03,
         0.22255664E-03,-0.31275427E-03,
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
    float x24 = x23*x2;
    float x25 = x24*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21*x31
        +coeff[  7]    *x21    *x41
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[  8]        *x31*x41
        +coeff[  9]            *x42
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]            *x41*x51
        +coeff[ 13]*x11*x21
        +coeff[ 14]                *x52
        +coeff[ 15]*x11        *x41
        +coeff[ 16]    *x22    *x41
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]    *x22        *x51
        +coeff[ 19]    *x21    *x41*x51
        +coeff[ 20]*x11*x22
        +coeff[ 21]*x11    *x32
        +coeff[ 22]*x11        *x42
        +coeff[ 23]*x11*x21        *x51
        +coeff[ 24]    *x25
        +coeff[ 25]    *x24*x31
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 26]*x11            *x53
        +coeff[ 27]*x11
        +coeff[ 28]*x11    *x31
        +coeff[ 29]    *x22*x31
        +coeff[ 30]    *x21*x31*x41
        +coeff[ 31]            *x43
        +coeff[ 32]*x12
        +coeff[ 33]            *x41*x52
        +coeff[ 34]    *x22    *x42
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 35]    *x22    *x41*x51
        +coeff[ 36]        *x32
        +coeff[ 37]    *x23
        +coeff[ 38]*x11            *x51
        +coeff[ 39]        *x31*x42
        +coeff[ 40]        *x31*x41*x51
        +coeff[ 41]    *x24
        +coeff[ 42]    *x21        *x52
        +coeff[ 43]*x11*x21*x31
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 44]    *x23    *x41
        +coeff[ 45]*x11    *x31*x41
        +coeff[ 46]*x11        *x41*x51
        +coeff[ 47]*x11*x23
        +coeff[ 48]*x12        *x41
        +coeff[ 49]*x11*x22    *x41
        +coeff[ 50]    *x21*x33    *x51
        +coeff[ 51]    *x25    *x41*x51
        +coeff[ 52]*x12            *x53
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 53]    *x23*x31
        +coeff[ 54]*x11*x21    *x41
        +coeff[ 55]    *x22*x32
        +coeff[ 56]                *x53
        +coeff[ 57]    *x21    *x43
        +coeff[ 58]    *x22*x31    *x51
        +coeff[ 59]    *x24        *x51
        +coeff[ 60]*x11*x21    *x41*x51
        +coeff[ 61]*x12*x21    *x41
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 62]    *x23*x31*x42
        +coeff[ 63]*x11*x23*x31    *x51
        +coeff[ 64]    *x21*x31    *x51
        +coeff[ 65]        *x32    *x51
        +coeff[ 66]            *x42*x51
        +coeff[ 67]        *x31    *x52
        +coeff[ 68]    *x23        *x51
        +coeff[ 69]    *x21*x31*x42
        +coeff[ 70]*x11    *x31    *x51
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 71]    *x21    *x42*x51
        +coeff[ 72]*x12    *x31
        +coeff[ 73]    *x21    *x41*x52
        +coeff[ 74]    *x23    *x42
        +coeff[ 75]*x12            *x51
        +coeff[ 76]*x11    *x31*x42
        +coeff[ 77]*x11        *x43
        +coeff[ 78]*x11*x21*x31    *x51
        +coeff[ 79]    *x22    *x42*x51
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 80]*x12*x22
        +coeff[ 81]    *x21*x32*x41*x51
        +coeff[ 82]*x11*x23    *x41
        +coeff[ 83]    *x22    *x41*x52
        +coeff[ 84]    *x24*x31*x41
        +coeff[ 85]*x11*x22    *x42
        +coeff[ 86]*x11*x21*x32*x41
        +coeff[ 87]    *x23*x32*x41
        +coeff[ 88]*x11*x21    *x43
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 89]    *x23*x31*x41*x51
        ;

    return v_p_l5p77_dex                             ;
}
float l_l5p77_dex                             (float *x,int m){
    int ncoeff= 67;
    float avdat= -0.5406602E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 68]={
         0.64887367E-02,-0.38065177E+00,-0.98172158E-01, 0.67930199E-01,
        -0.10367237E-01,-0.14510611E+00,-0.34144595E-01,-0.33121336E-01,
         0.11630415E-02, 0.17565282E-01,-0.68865472E-03, 0.18207202E-01,
        -0.52902669E-01, 0.11228444E-01,-0.12448854E-01,-0.37470281E-01,
         0.62446009E-01,-0.25589779E-01,-0.24631924E-02, 0.63313884E-02,
        -0.31941582E-02,-0.27846921E-01, 0.40703721E-01,-0.66303080E-02,
         0.55587307E-01, 0.29429523E-03, 0.35044014E-01,-0.92851678E-02,
         0.17606631E-02,-0.85612917E-02, 0.70579085E-02,-0.62046307E-02,
         0.97562531E-02, 0.14149994E-03, 0.11826725E-01, 0.60940736E-04,
        -0.44782273E-02, 0.48164758E-03,-0.26351076E-02,-0.22035337E-02,
         0.71602450E-02, 0.72939279E-02, 0.73525482E-02, 0.19751461E-02,
         0.13149999E-02, 0.14414343E-02,-0.38911358E-02,-0.66476536E-03,
        -0.17033571E-01,-0.15030618E-01, 0.14258667E-02,-0.31729804E-02,
         0.27090104E-02, 0.36539557E-05, 0.15376926E-03, 0.12541312E-02,
        -0.65344438E-03, 0.55333332E-03,-0.58152468E-03,-0.46034515E-03,
        -0.61262180E-02, 0.49716826E-02,-0.11954006E-02, 0.32780068E-02,
         0.97128347E-03,-0.48288852E-02, 0.37035823E-02,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_l_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]                *x51
        +coeff[  3]*x11
        +coeff[  4]    *x22
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x21        *x51
        +coeff[  7]*x11        *x41
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[  8]            *x43
        +coeff[  9]    *x23*x31
        +coeff[ 10]    *x21*x33
        +coeff[ 11]            *x41
        +coeff[ 12]    *x21*x31
        +coeff[ 13]*x11*x21
        +coeff[ 14]*x11    *x31
        +coeff[ 15]    *x23
        +coeff[ 16]    *x21    *x42
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 17]*x11*x22
        +coeff[ 18]    *x23*x31*x41
        +coeff[ 19]        *x31
        +coeff[ 20]*x11            *x51
        +coeff[ 21]    *x22    *x41
        +coeff[ 22]    *x21*x31*x41
        +coeff[ 23]*x12*x21
        +coeff[ 24]    *x23    *x41
        +coeff[ 25]            *x44
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 26]*x11*x22    *x41
        +coeff[ 27]            *x42
        +coeff[ 28]                *x52
        +coeff[ 29]    *x22*x31
        +coeff[ 30]    *x21*x32
        +coeff[ 31]*x11*x21    *x41
        +coeff[ 32]*x11        *x42
        +coeff[ 33]        *x33*x41
        +coeff[ 34]*x11*x22*x31
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 35]*x11    *x33*x41
        +coeff[ 36]        *x31*x41
        +coeff[ 37]            *x41*x51
        +coeff[ 38]    *x21    *x41*x51
        +coeff[ 39]*x11*x21*x31
        +coeff[ 40]*x11    *x31*x41
        +coeff[ 41]*x12*x21    *x41
        +coeff[ 42]    *x23    *x44
        +coeff[ 43]            *x42*x51
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 44]    *x21        *x52
        +coeff[ 45]*x11    *x32
        +coeff[ 46]    *x24
        +coeff[ 47]        *x34
        +coeff[ 48]    *x21*x31*x42
        +coeff[ 49]    *x21    *x43
        +coeff[ 50]    *x22    *x41*x51
        +coeff[ 51]*x11*x23
        +coeff[ 52]*x12*x21*x31
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 53]    *x23*x32*x41
        +coeff[ 54]        *x31    *x51
        +coeff[ 55]        *x31*x42
        +coeff[ 56]    *x22        *x51
        +coeff[ 57]        *x31*x41*x51
        +coeff[ 58]            *x41*x52
        +coeff[ 59]*x13
        +coeff[ 60]    *x21*x32*x41
        +coeff[ 61]    *x22    *x42
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 62]    *x23        *x51
        +coeff[ 63]    *x24    *x41
        +coeff[ 64]    *x22*x32*x41
        +coeff[ 65]    *x23    *x42
        +coeff[ 66]    *x24*x31*x41
        ;

    return v_l_l5p77_dex                             ;
}
float x_l5p77_fp                              (float *x,int m){
    int ncoeff= 49;
    float avdat=  0.1037642E-01;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 50]={
         0.43306914E-02,-0.37461292E-01, 0.65247756E+00,-0.53094141E-01,
         0.43631431E-01, 0.10533310E-01, 0.34232263E-01, 0.13113662E-01,
        -0.13049352E-01,-0.28302468E-01,-0.31553069E-02, 0.71598189E-02,
        -0.35838955E-02, 0.55334871E-02, 0.13583186E-01,-0.63093682E-03,
        -0.92198449E-03,-0.84599864E-03,-0.23410739E-02, 0.23430891E-02,
         0.23224298E-02,-0.36624696E-02, 0.83011631E-02,-0.18580198E-01,
         0.57500147E-05, 0.64136641E-03, 0.13909411E-02, 0.46201651E-02,
         0.29571308E-02, 0.49823741E-02,-0.10104933E-01,-0.80614034E-02,
        -0.60537928E-02,-0.41016904E-02, 0.64769480E-02, 0.82695931E-02,
        -0.50676841E-03, 0.90390758E-03, 0.91962179E-03,-0.10226306E-02,
         0.95622899E-03,-0.16535607E-02, 0.91654388E-03, 0.95703517E-03,
         0.29396894E-02,-0.18439010E-02,-0.37418739E-02, 0.26034771E-02,
        -0.20723392E-02,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_x_l5p77_fp                              =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]                *x52
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x21*x31
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[  8]            *x42
        +coeff[  9]    *x21    *x42
        +coeff[ 10]            *x41
        +coeff[ 11]*x11        *x41
        +coeff[ 12]        *x31*x41
        +coeff[ 13]                *x53
        +coeff[ 14]    *x21    *x41*x51
        +coeff[ 15]    *x23*x31*x41
        +coeff[ 16]    *x21
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 17]        *x31
        +coeff[ 18]*x11            *x51
        +coeff[ 19]*x11    *x31
        +coeff[ 20]            *x41*x51
        +coeff[ 21]    *x21        *x52
        +coeff[ 22]    *x23
        +coeff[ 23]    *x23    *x41
        +coeff[ 24]*x13*x22
        +coeff[ 25]        *x31    *x51
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 26]*x12*x21
        +coeff[ 27]*x11*x22
        +coeff[ 28]    *x21*x31    *x51
        +coeff[ 29]    *x22    *x41
        +coeff[ 30]    *x21*x31*x41
        +coeff[ 31]*x11*x22    *x41
        +coeff[ 32]*x11*x21    *x42
        +coeff[ 33]    *x21    *x42*x51
        +coeff[ 34]    *x22    *x42
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 35]*x11*x23    *x42
        +coeff[ 36]*x12
        +coeff[ 37]*x11            *x52
        +coeff[ 38]*x11*x21    *x41
        +coeff[ 39]*x11    *x31*x41
        +coeff[ 40]            *x42*x51
        +coeff[ 41]    *x21*x32
        +coeff[ 42]            *x43
        +coeff[ 43]*x11*x23
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 44]    *x23        *x51
        +coeff[ 45]    *x21*x31*x41*x51
        +coeff[ 46]    *x23*x31
        +coeff[ 47]*x13*x22        *x51
        +coeff[ 48]*x13*x22*x31
        ;

    return v_x_l5p77_fp                              ;
}
float t_l5p77_fp                              (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1062570E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.11077350E-02,-0.33733402E-02,-0.18544450E-01, 0.17614244E-04,
         0.53349595E-04, 0.11301004E+00, 0.53365729E-02,-0.10616464E-01,
         0.12934288E-02,-0.23490486E-02, 0.15437859E-03,-0.76126872E-03,
        -0.77129283E-03,-0.22160634E-02, 0.19795976E-02,-0.79794251E-03,
         0.11125805E-02,-0.42822622E-03, 0.36771514E-03,-0.72328359E-04,
        -0.18052294E-03,-0.26631175E-03, 0.14183560E-03,-0.42838341E-03,
         0.37646177E-03, 0.36626146E-03, 0.29599539E-03, 0.83670061E-03,
        -0.77796273E-03, 0.54008586E-04, 0.13073875E-02,-0.92502400E-04,
        -0.24945231E-03, 0.34539666E-04,-0.27379484E-03,-0.20809581E-03,
         0.31029331E-03, 0.15211779E-03,-0.31875214E-03, 0.51623362E-03,
         0.11671873E-03,-0.56759396E-03,-0.65548613E-03, 0.55199582E-03,
        -0.73636409E-04, 0.68878988E-04, 0.15802577E-03,-0.38055772E-04,
         0.10173125E-03, 0.98714750E-04,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_l5p77_fp                              =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]            *x41
        +coeff[  5]                *x51
        +coeff[  6]    *x21        *x51
        +coeff[  7]                *x52
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[  8]    *x22
        +coeff[  9]            *x42
        +coeff[ 10]*x11*x21
        +coeff[ 11]    *x21    *x41
        +coeff[ 12]        *x31*x41
        +coeff[ 13]    *x21    *x42
        +coeff[ 14]    *x21    *x41*x51
        +coeff[ 15]    *x21        *x52
        +coeff[ 16]                *x53
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 17]*x11            *x51
        +coeff[ 18]            *x41*x51
        +coeff[ 19]*x12
        +coeff[ 20]    *x21*x31
        +coeff[ 21]*x11        *x41
        +coeff[ 22]        *x31    *x51
        +coeff[ 23]*x11*x22
        +coeff[ 24]*x11        *x42
        +coeff[ 25]    *x21*x31    *x51
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 26]            *x42*x51
        +coeff[ 27]*x11*x23
        +coeff[ 28]    *x23    *x41
        +coeff[ 29]    *x22*x31*x41
        +coeff[ 30]    *x22    *x42
        +coeff[ 31]*x11    *x31
        +coeff[ 32]    *x23
        +coeff[ 33]*x11*x21*x31
        +coeff[ 34]    *x22*x31
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 35]    *x22    *x41
        +coeff[ 36]            *x43
        +coeff[ 37]*x11            *x52
        +coeff[ 38]    *x21    *x43
        +coeff[ 39]*x11*x22        *x51
        +coeff[ 40]    *x23        *x51
        +coeff[ 41]    *x21    *x42*x51
        +coeff[ 42]    *x23*x31*x41*x51
        +coeff[ 43]    *x23        *x53
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 44]*x12*x21
        +coeff[ 45]*x11    *x31*x41
        +coeff[ 46]        *x31*x42
        +coeff[ 47]        *x32    *x51
        +coeff[ 48]*x11        *x41*x51
        +coeff[ 49]        *x31*x41*x51
        ;

    return v_t_l5p77_fp                              ;
}
float y_l5p77_fp                              (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.7994905E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
         0.52612768E-02,-0.57620980E-01, 0.33640799E-02, 0.54813368E-03,
         0.12045373E-01, 0.76661306E-02,-0.15869256E-01,-0.21559799E-02,
        -0.39944821E-02,-0.64756460E-02,-0.26825587E-02,-0.46120640E-02,
        -0.28282397E-02, 0.99318596E-02,-0.57701692E-02, 0.24748833E-04,
         0.66861110E-02, 0.53515318E-02, 0.41033053E-02, 0.29814264E-02,
         0.12191844E-03, 0.44704559E-02, 0.16130925E-03, 0.78950386E-03,
         0.31576527E-03,-0.26144295E-02, 0.77399111E-03,-0.77533396E-03,
         0.37326076E-03,-0.15962466E-02, 0.18648774E-02, 0.43839901E-02,
         0.14665494E-02, 0.35643280E-02, 0.95052342E-03, 0.90258906E-03,
         0.14716295E-02, 0.27762835E-02, 0.34358655E-03,-0.43472103E-02,
        -0.36774177E-03, 0.87315863E-03, 0.89557172E-03,-0.10848453E-02,
        -0.38841788E-02,-0.26848079E-02, 0.22077179E-02,-0.49715969E-02,
         0.19980804E-02,-0.14557256E-02,-0.10368233E-02,-0.11291369E-02,
         0.42034776E-03,-0.23467009E-03,-0.12318480E-02,-0.20351436E-02,
        -0.35835567E-03, 0.63619530E-03,-0.13908915E-02,-0.53971203E-03,
         0.81699312E-03, 0.10928685E-02, 0.57058013E-03, 0.30318982E-03,
        -0.75380172E-04,-0.56452909E-03,-0.74342033E-03,-0.19152529E-03,
         0.21839396E-02, 0.75005635E-03, 0.52176462E-03,-0.14441682E-03,
        -0.12458458E-02,-0.93929307E-03, 0.17294186E-03,-0.17310370E-03,
         0.67673338E-03,-0.51215968E-04,-0.45142876E-03,-0.53760281E-03,
         0.69358974E-03,-0.31132935E-03,-0.34575808E-03, 0.37477849E-03,
        -0.17718188E-03,-0.14595651E-02, 0.86034561E-03,-0.63177536E-03,
        -0.75414102E-03, 0.34884442E-03,-0.22573124E-03,-0.65047829E-03,
        -0.53947623E-03, 0.50399051E-03, 0.71309641E-03,-0.11509071E-02,
        -0.92779614E-04,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_fp                              =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]*x11
        +coeff[  4]                *x51
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x22
        +coeff[  7]*x11        *x41
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[  8]            *x41*x51
        +coeff[  9]*x11*x21
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]                *x52
        +coeff[ 13]    *x22    *x41
        +coeff[ 14]    *x21    *x41*x51
        +coeff[ 15]            *x44
        +coeff[ 16]            *x41*x52
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 17]            *x42
        +coeff[ 18]        *x31*x41
        +coeff[ 19]    *x21    *x42
        +coeff[ 20]        *x31*x42
        +coeff[ 21]*x11*x21    *x41
        +coeff[ 22]        *x31*x42*x52
        +coeff[ 23]    *x22    *x41*x53
        +coeff[ 24]    *x24*x31    *x52
        +coeff[ 25]        *x31
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 26]        *x32
        +coeff[ 27]*x12
        +coeff[ 28]*x11            *x51
        +coeff[ 29]            *x42*x51
        +coeff[ 30]    *x23
        +coeff[ 31]    *x22*x31
        +coeff[ 32]*x11*x21*x31
        +coeff[ 33]    *x22        *x51
        +coeff[ 34]*x11        *x41*x51
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 35]*x11*x21        *x51
        +coeff[ 36]    *x21        *x52
        +coeff[ 37]        *x31    *x52
        +coeff[ 38]            *x45
        +coeff[ 39]            *x41*x53
        +coeff[ 40]    *x21*x31
        +coeff[ 41]    *x21*x31*x41
        +coeff[ 42]*x11        *x42
        +coeff[ 43]        *x31*x41*x51
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 44]    *x22    *x42
        +coeff[ 45]    *x22*x31*x41
        +coeff[ 46]    *x24
        +coeff[ 47]    *x22    *x41*x51
        +coeff[ 48]*x11*x23
        +coeff[ 49]    *x23        *x51
        +coeff[ 50]    *x22*x31    *x51
        +coeff[ 51]    *x21    *x44*x51
        +coeff[ 52]*x12        *x41
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 53]*x12*x21
        +coeff[ 54]            *x43*x51
        +coeff[ 55]*x11*x21    *x42
        +coeff[ 56]*x11            *x52
        +coeff[ 57]                *x53
        +coeff[ 58]*x11*x21*x31*x41
        +coeff[ 59]*x11        *x42*x51
        +coeff[ 60]            *x42*x52
        +coeff[ 61]    *x21    *x41*x52
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 62]*x12*x22
        +coeff[ 63]    *x21*x31    *x52
        +coeff[ 64]*x11    *x31    *x52
        +coeff[ 65]    *x21        *x53
        +coeff[ 66]        *x31    *x53
        +coeff[ 67]*x11    *x33*x41
        +coeff[ 68]    *x22    *x41*x52
        +coeff[ 69]            *x43
        +coeff[ 70]*x11    *x31*x41
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 71]        *x32    *x51
        +coeff[ 72]    *x21    *x43
        +coeff[ 73]    *x21*x31*x42
        +coeff[ 74]*x12    *x31
        +coeff[ 75]*x11    *x31    *x51
        +coeff[ 76]    *x23    *x41
        +coeff[ 77]*x11    *x31*x42
        +coeff[ 78]    *x23*x31
        +coeff[ 79]    *x22*x32
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 80]*x11*x22    *x41
        +coeff[ 81]*x11*x21*x32
        +coeff[ 82]        *x32*x43
        +coeff[ 83]*x12*x21    *x41
        +coeff[ 84]*x11*x21*x31    *x51
        +coeff[ 85]    *x22    *x42*x51
        +coeff[ 86]    *x23    *x41*x51
        +coeff[ 87]    *x24        *x51
        +coeff[ 88]*x11*x22    *x41*x51
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 89]*x11        *x42*x52
        +coeff[ 90]*x11*x21    *x41*x52
        +coeff[ 91]    *x21    *x41*x53
        +coeff[ 92]    *x22        *x53
        +coeff[ 93]*x11        *x41*x53
        +coeff[ 94]    *x24    *x43
        +coeff[ 95]*x11*x23    *x43
        +coeff[ 96]*x11    *x31
        ;

    return v_y_l5p77_fp                              ;
}
float p_l5p77_fp                              (float *x,int m){
    int ncoeff= 90;
    float avdat= -0.6457705E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.25734296E-02, 0.32411383E-02, 0.60911034E-02,-0.30777883E-01,
         0.10644548E-01,-0.14892758E-01, 0.13716331E-02, 0.17594201E-01,
         0.36614942E-02,-0.23709750E-02,-0.38738709E-02,-0.20642057E-01,
        -0.61980207E-02, 0.14962205E-02,-0.35300148E-02, 0.15989168E-01,
         0.52628675E-02,-0.32384943E-02, 0.59830816E-03, 0.17583127E-02,
         0.41104437E-03,-0.17174665E-03,-0.23711554E-03,-0.22580120E-03,
        -0.99513323E-04, 0.29054083E-03, 0.51706872E-03, 0.49399687E-02,
         0.43355818E-02, 0.55438938E-03, 0.22564046E-02,-0.70602150E-03,
         0.28910018E-02, 0.74571947E-03,-0.47682339E-03, 0.20309547E-02,
         0.35047298E-02, 0.11528452E-02, 0.22490008E-02, 0.87272003E-03,
        -0.80529059E-03, 0.19641961E-02, 0.13964284E-02, 0.28606895E-02,
        -0.11645516E-03, 0.53290260E-03, 0.10502943E-03,-0.24808542E-03,
         0.80540805E-03,-0.22020822E-02, 0.57125848E-03,-0.81395579E-03,
         0.46219447E-03,-0.30892456E-03, 0.12505148E-03,-0.49934589E-03,
        -0.90310158E-03,-0.15855072E-02,-0.19793648E-03,-0.15011493E-02,
        -0.11804224E-03, 0.63604955E-03,-0.11224736E-03, 0.17360311E-03,
        -0.94497664E-03,-0.19829329E-02,-0.15271318E-03, 0.18633236E-03,
         0.57227793E-03, 0.14452003E-03,-0.12087314E-02, 0.48471385E-03,
         0.34719027E-03,-0.11468756E-02,-0.54965715E-03, 0.31892494E-04,
        -0.15535601E-02,-0.67248126E-03, 0.10828866E-03,-0.29323765E-03,
         0.83991799E-04,-0.74399740E-03,-0.20709293E-03,-0.10470693E-03,
         0.29516927E-03, 0.10062858E-03,-0.19615554E-03, 0.25176935E-03,
         0.71487582E-03, 0.20488727E-03,
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
    float x24 = x23*x2;
    float x25 = x24*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_fp                              =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21*x31
        +coeff[  7]    *x21    *x41
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[  8]        *x31*x41
        +coeff[  9]    *x21        *x51
        +coeff[ 10]        *x31    *x51
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]                *x52
        +coeff[ 14]*x11        *x41
        +coeff[ 15]    *x22    *x41
        +coeff[ 16]    *x21    *x42
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 17]    *x22        *x51
        +coeff[ 18]*x11*x22
        +coeff[ 19]*x11        *x42
        +coeff[ 20]    *x22    *x42
        +coeff[ 21]        *x32*x42
        +coeff[ 22]            *x44
        +coeff[ 23]    *x25
        +coeff[ 24]    *x24*x31
        +coeff[ 25]        *x32*x44
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 26]*x11
        +coeff[ 27]            *x42
        +coeff[ 28]    *x22*x31
        +coeff[ 29]*x11            *x51
        +coeff[ 30]            *x43
        +coeff[ 31]*x12
        +coeff[ 32]            *x41*x52
        +coeff[ 33]        *x32
        +coeff[ 34]*x11    *x31
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 35]    *x21*x31*x41
        +coeff[ 36]    *x24
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]*x11*x21    *x41
        +coeff[ 39]        *x31    *x52
        +coeff[ 40]*x11*x21        *x51
        +coeff[ 41]*x11*x23
        +coeff[ 42]*x11*x22    *x41
        +coeff[ 43]    *x23
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 44]        *x32*x41
        +coeff[ 45]        *x31*x42
        +coeff[ 46]        *x31*x41*x51
        +coeff[ 47]            *x42*x51
        +coeff[ 48]*x11    *x31*x41
        +coeff[ 49]    *x21    *x43
        +coeff[ 50]*x12        *x41
        +coeff[ 51]            *x41*x53
        +coeff[ 52]*x12*x21    *x41
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 53]    *x21    *x41*x51
        +coeff[ 54]    *x21        *x52
        +coeff[ 55]    *x23*x31
        +coeff[ 56]    *x23    *x41
        +coeff[ 57]    *x22*x31*x41
        +coeff[ 58]    *x21*x32*x41
        +coeff[ 59]    *x21*x31*x42
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]    *x22    *x41*x51
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 62]*x11            *x52
        +coeff[ 63]*x12    *x31
        +coeff[ 64]*x11*x21    *x42
        +coeff[ 65]    *x23    *x42
        +coeff[ 66]*x12            *x51
        +coeff[ 67]    *x23*x31    *x51
        +coeff[ 68]*x11*x21    *x41*x51
        +coeff[ 69]*x11    *x31*x41*x51
        +coeff[ 70]    *x22    *x42*x51
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 71]*x12*x22
        +coeff[ 72]    *x24*x32
        +coeff[ 73]*x11*x22    *x42
        +coeff[ 74]*x11*x21    *x43
        +coeff[ 75]*x12*x23
        +coeff[ 76]*x11*x23*x31*x41
        +coeff[ 77]            *x43*x53
        +coeff[ 78]    *x21*x32
        +coeff[ 79]    *x21*x31    *x51
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 80]*x11    *x32
        +coeff[ 81]    *x22*x32
        +coeff[ 82]    *x23        *x51
        +coeff[ 83]                *x53
        +coeff[ 84]    *x22*x31    *x51
        +coeff[ 85]*x11        *x41*x51
        +coeff[ 86]*x12*x21
        +coeff[ 87]    *x22        *x52
        +coeff[ 88]    *x24    *x41
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 89]    *x21*x31    *x52
        ;

    return v_p_l5p77_fp                              ;
}
float l_l5p77_fp                              (float *x,int m){
    int ncoeff= 82;
    float avdat= -0.1011234E-01;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 83]={
         0.12967309E-01,-0.24596201E+00, 0.14018445E-01,-0.31420160E-01,
         0.41501168E-01,-0.14260489E-01,-0.87614104E-01,-0.26619842E-01,
         0.96398061E-02,-0.20599537E-01, 0.10297282E-01,-0.44491267E-03,
         0.53376313E-02,-0.32726616E-01,-0.12085279E-01,-0.79004960E-02,
        -0.22222411E-01,-0.28312158E-01, 0.36639787E-01,-0.15956016E-01,
        -0.13235122E-02,-0.31007160E-02,-0.18545851E-02,-0.85584857E-02,
         0.24096835E-01, 0.49483934E-02,-0.40051197E-02, 0.33379640E-01,
         0.21799650E-01,-0.37693859E-02, 0.17785559E-02, 0.42977934E-02,
        -0.64118840E-02, 0.58440422E-02, 0.11097166E-03, 0.33454859E-03,
        -0.62351080E-03,-0.24643447E-02, 0.45210510E-02,-0.45893723E-02,
         0.71560610E-02, 0.43278639E-02, 0.13837379E-02, 0.16375623E-02,
         0.12353054E-02,-0.13710224E-02, 0.92269119E-03, 0.11916523E-02,
         0.69854427E-02,-0.10842112E-01,-0.89773117E-02, 0.10642965E-02,
        -0.39104391E-02,-0.11174220E-02, 0.16225617E-02, 0.71113650E-02,
         0.27819281E-02, 0.41882161E-03,-0.11369223E-02, 0.20611615E-03,
         0.47109381E-03,-0.24518318E-03, 0.12383243E-02, 0.16801758E-02,
        -0.76396868E-03, 0.49226370E-03,-0.30085209E-03, 0.29678186E-03,
         0.31723554E-03,-0.31795976E-03, 0.42904890E-02,-0.43207975E-02,
        -0.13085817E-02,-0.73265738E-03, 0.35302114E-03,-0.84127689E-03,
         0.42955152E-03,-0.91743184E-03, 0.20087708E-02, 0.15999562E-02,
         0.26319672E-02, 0.16398882E-02,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;
    float x54 = x53*x5;

//                 function

    float v_l_l5p77_fp                              =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]*x11
        +coeff[  5]    *x22
        +coeff[  6]    *x21    *x41
        +coeff[  7]                *x52
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[  8]*x11*x21
        +coeff[  9]*x11        *x41
        +coeff[ 10]    *x23*x31
        +coeff[ 11]    *x21*x33
        +coeff[ 12]        *x31
        +coeff[ 13]    *x21*x31
        +coeff[ 14]            *x42
        +coeff[ 15]*x11    *x31
        +coeff[ 16]    *x23
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]*x11*x22
        +coeff[ 20]    *x23*x31*x41
        +coeff[ 21]    *x21        *x51
        +coeff[ 22]*x11            *x51
        +coeff[ 23]    *x22*x31
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]                *x53
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 26]*x12*x21
        +coeff[ 27]    *x23    *x41
        +coeff[ 28]*x11*x22    *x41
        +coeff[ 29]        *x31*x41
        +coeff[ 30]            *x41*x51
        +coeff[ 31]    *x21*x32
        +coeff[ 32]*x11*x21    *x41
        +coeff[ 33]*x11        *x42
        +coeff[ 34]*x11    *x33*x41
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 35]*x11*x24*x31
        +coeff[ 36]        *x32
        +coeff[ 37]    *x21        *x52
        +coeff[ 38]*x11    *x31*x41
        +coeff[ 39]    *x24
        +coeff[ 40]*x11*x22*x31
        +coeff[ 41]*x12*x21    *x41
        +coeff[ 42]    *x21*x31    *x51
        +coeff[ 43]    *x21    *x41*x51
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 44]        *x31*x41*x51
        +coeff[ 45]*x11*x21*x31
        +coeff[ 46]*x11    *x32
        +coeff[ 47]*x11        *x41*x51
        +coeff[ 48]    *x22    *x42
        +coeff[ 49]    *x21*x31*x42
        +coeff[ 50]    *x21    *x43
        +coeff[ 51]    *x21    *x42*x51
        +coeff[ 52]*x11*x23
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 53]*x12*x22
        +coeff[ 54]*x12*x21*x31
        +coeff[ 55]    *x24    *x41
        +coeff[ 56]    *x21    *x44
        +coeff[ 57]            *x44*x51
        +coeff[ 58]    *x24*x31*x41
        +coeff[ 59]    *x23*x32*x41
        +coeff[ 60]    *x21*x34*x41
        +coeff[ 61]        *x31    *x51
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 62]        *x31*x42
        +coeff[ 63]            *x43
        +coeff[ 64]    *x22        *x51
        +coeff[ 65]            *x41*x52
        +coeff[ 66]*x11*x21        *x51
        +coeff[ 67]*x11    *x31    *x51
        +coeff[ 68]*x11            *x52
        +coeff[ 69]*x13
        +coeff[ 70]    *x22*x31*x41
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 71]    *x21*x32*x41
        +coeff[ 72]    *x21    *x41*x52
        +coeff[ 73]            *x42*x52
        +coeff[ 74]    *x21        *x53
        +coeff[ 75]                *x54
        +coeff[ 76]*x11*x22        *x51
        +coeff[ 77]*x11        *x42*x51
        +coeff[ 78]    *x24*x31
        +coeff[ 79]*x11*x24
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 80]*x11*x23    *x41
        +coeff[ 81]*x11*x23    *x42
        ;

    return v_l_l5p77_fp                              ;
}
float x_l5p77_q1en                            (float *x,int m){
    int ncoeff= 45;
    float avdat=  0.9247321E-03;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25772E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 46]={
        -0.10315748E-02, 0.96696965E-01, 0.19487011E-02, 0.76378630E-02,
         0.37627149E-04, 0.70906531E-05, 0.17782291E-02,-0.73287993E-04,
        -0.14848435E-04,-0.72262774E-04, 0.12862576E-01, 0.30075454E-02,
        -0.72480959E-03, 0.74155198E-03, 0.23941428E-02,-0.27864827E-02,
        -0.61561586E-04,-0.27187905E-03, 0.31865327E-03,-0.64240926E-03,
        -0.20493537E-02, 0.31518142E-04,-0.20902485E-02,-0.98779004E-04,
        -0.32928088E-03, 0.22349374E-03, 0.40226191E-03,-0.60385867E-03,
        -0.21952430E-03,-0.37413673E-03, 0.56038320E-03,-0.14184947E-02,
        -0.94859145E-03, 0.18256123E-04,-0.84142150E-04,-0.32507349E-04,
         0.78213670E-04, 0.15425986E-03, 0.21133533E-03,-0.40329382E-03,
         0.15858028E-03,-0.26234999E-03,-0.47105015E-03, 0.55953005E-03,
         0.87792438E-03,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_x_l5p77_q1en                            =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]*x11        *x41
        +coeff[  3]    *x21    *x41
        +coeff[  4]*x13
        +coeff[  5]*x11            *x52
        +coeff[  6]*x11*x22
        +coeff[  7]*x11    *x32
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[  8]*x13            *x52
        +coeff[  9]*x13*x22
        +coeff[ 10]*x11
        +coeff[ 11]    *x21*x31
        +coeff[ 12]            *x41
        +coeff[ 13]*x11    *x31
        +coeff[ 14]    *x23
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x12*x21*x31*x41
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[ 17]        *x31
        +coeff[ 18]    *x21        *x51
        +coeff[ 19]    *x22
        +coeff[ 20]    *x21*x31*x41
        +coeff[ 21]*x13*x21
        +coeff[ 22]    *x23    *x41
        +coeff[ 23]*x11*x22    *x43
        +coeff[ 24]*x11*x21
        +coeff[ 25]            *x42
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[ 26]*x12*x21
        +coeff[ 27]*x11        *x42
        +coeff[ 28]    *x21    *x41*x51
        +coeff[ 29]    *x21*x32
        +coeff[ 30]    *x22    *x41
        +coeff[ 31]*x11*x22    *x41
        +coeff[ 32]    *x23*x31
        +coeff[ 33]        *x31*x43
        +coeff[ 34]*x13    *x31*x41
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[ 35]                *x51
        +coeff[ 36]*x11            *x51
        +coeff[ 37]        *x31*x41
        +coeff[ 38]*x11*x21    *x41
        +coeff[ 39]*x11    *x31*x41
        +coeff[ 40]    *x22*x31
        +coeff[ 41]*x12*x21    *x41
        +coeff[ 42]*x11*x22*x31
        +coeff[ 43]    *x21    *x43
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[ 44]    *x23*x31*x42
        ;

    return v_x_l5p77_q1en                            ;
}
float t_l5p77_q1en                            (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.7094675E-03;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25772E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
        -0.21852220E-03,-0.23623419E-02, 0.36474217E-01, 0.99846041E-02,
        -0.17395499E-03,-0.11166484E-02,-0.27803428E-04,-0.85469667E-03,
         0.84557547E-03, 0.36334158E-02, 0.23019994E-02,-0.31357704E-03,
        -0.70729438E-03, 0.36669630E-03, 0.26941334E-02,-0.31430349E-02,
        -0.34460365E-02,-0.41249023E-04,-0.76629649E-05,-0.71206974E-04,
        -0.33394288E-03, 0.47634672E-04, 0.45219352E-03, 0.19445403E-02,
         0.77563256E-03,-0.23383251E-02,-0.63587626E-03,-0.20820936E-02,
        -0.45202058E-04, 0.22448457E-03, 0.26096342E-03, 0.32021268E-03,
        -0.46754497E-03,-0.36238594E-03,-0.71400072E-03,-0.47839279E-03,
         0.60385319E-04, 0.12551575E-02, 0.12128768E-03, 0.41078663E-03,
        -0.32157350E-04, 0.13845872E-03, 0.81699167E-04, 0.10726824E-03,
        -0.80251943E-04,-0.41752323E-03,-0.12757660E-03, 0.10066279E-02,
         0.13499477E-03, 0.64346468E-03,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;

//                 function

    float v_t_l5p77_q1en                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]    *x21    *x41
        +coeff[  4]*x12*x21*x31
        +coeff[  5]    *x23*x31
        +coeff[  6]*x13        *x41
        +coeff[  7]            *x41
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[  8]*x11    *x31
        +coeff[  9]    *x21*x31
        +coeff[ 10]*x11        *x41
        +coeff[ 11]        *x31
        +coeff[ 12]    *x22
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]    *x23
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]    *x23    *x41
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 17]*x13*x22
        +coeff[ 18]*x12*x21*x31*x41
        +coeff[ 19]    *x23*x31*x41
        +coeff[ 20]*x11*x21
        +coeff[ 21]*x13
        +coeff[ 22]*x12*x21
        +coeff[ 23]*x11*x22
        +coeff[ 24]    *x22    *x41
        +coeff[ 25]    *x21*x31*x41
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 26]*x11        *x42
        +coeff[ 27]*x11*x22    *x41
        +coeff[ 28]    *x23*x32
        +coeff[ 29]            *x42
        +coeff[ 30]    *x22*x31
        +coeff[ 31]*x11*x21    *x41
        +coeff[ 32]*x11    *x31*x41
        +coeff[ 33]    *x21    *x41*x51
        +coeff[ 34]*x11*x22*x31
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 35]*x12*x21    *x41
        +coeff[ 36]*x12    *x31*x41
        +coeff[ 37]    *x21    *x43
        +coeff[ 38]*x12*x21*x31*x42
        +coeff[ 39]    *x23*x31*x42
        +coeff[ 40]                *x51
        +coeff[ 41]        *x31*x41
        +coeff[ 42]*x11            *x51
        +coeff[ 43]*x11*x21*x31
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 44]*x11    *x32
        +coeff[ 45]    *x21*x32
        +coeff[ 46]    *x21*x31    *x51
        +coeff[ 47]    *x21*x31*x42
        +coeff[ 48]*x12*x21*x32*x41
        +coeff[ 49]    *x23*x32*x41
        ;

    return v_t_l5p77_q1en                            ;
}
float y_l5p77_q1en                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.6351584E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25772E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.76303845E-02, 0.70959516E-01, 0.17093671E-01,-0.41945037E-02,
         0.65124105E-02, 0.30298170E-02,-0.11493454E-02,-0.49026869E-02,
        -0.27331761E-02,-0.20429711E-02,-0.18616440E-02,-0.24288183E-02,
        -0.28066771E-03, 0.87917043E-03, 0.38175055E-03,-0.97860466E-03,
         0.33508902E-03,-0.39119524E-03, 0.21107090E-03, 0.28755059E-02,
         0.19603190E-02, 0.55942328E-05,-0.19438843E-03,-0.39417503E-03,
        -0.30477403E-03, 0.32553132E-03, 0.23790900E-03,-0.26846130E-03,
        -0.30237716E-03, 0.12946153E-02,-0.88134978E-03, 0.38061815E-03,
         0.90861128E-03,-0.84069552E-03,-0.38898797E-05, 0.21366021E-03,
         0.46164954E-04,-0.89942921E-04, 0.49254304E-03, 0.53687498E-03,
         0.20862829E-03,-0.10948780E-03,-0.12482876E-03,-0.31995261E-03,
         0.86162261E-04,-0.84338433E-04, 0.97699725E-04,-0.57945399E-04,
        -0.66088716E-04, 0.85286585E-04, 0.57972673E-04,-0.16746671E-03,
        -0.12811457E-03, 0.16282087E-03, 0.11973150E-03, 0.18084611E-03,
         0.87420172E-04, 0.61798393E-03, 0.61257702E-03, 0.11309263E-04,
         0.31280666E-04,-0.21527008E-04,-0.12276181E-04, 0.15641182E-04,
        -0.54526165E-04,-0.17307668E-04, 0.72611991E-04,-0.52854052E-03,
         0.79215148E-04,-0.72237704E-03,-0.27227183E-03,-0.40888826E-04,
         0.24248517E-03, 0.19423531E-03,-0.15638424E-03,-0.31067751E-03,
        -0.13184229E-03,-0.27361390E-03,-0.27038099E-03,-0.10939185E-03,
         0.15168807E-03, 0.14128155E-03, 0.77412798E-04, 0.85273445E-04,
        -0.11657823E-04,-0.31560357E-04, 0.87390612E-06,-0.23662751E-04,
         0.31287727E-04, 0.18496299E-03,-0.46922531E-04, 0.36649762E-04,
        -0.53634831E-04,-0.64310196E-04,-0.11252577E-03, 0.55374174E-04,
        -0.40869947E-03,
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
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_q1en                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]        *x31
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]*x11*x21
        +coeff[  6]    *x21
        +coeff[  7]    *x22    *x41
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[  8]            *x42
        +coeff[  9]        *x31*x41
        +coeff[ 10]    *x22*x31
        +coeff[ 11]*x11*x21    *x41
        +coeff[ 12]*x11
        +coeff[ 13]    *x21    *x41
        +coeff[ 14]*x12
        +coeff[ 15]*x11*x21*x31
        +coeff[ 16]    *x21*x31
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 17]        *x32
        +coeff[ 18]                *x52
        +coeff[ 19]    *x22    *x42
        +coeff[ 20]    *x22*x31*x41
        +coeff[ 21]*x13        *x41
        +coeff[ 22]            *x41*x51
        +coeff[ 23]    *x21    *x42
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]    *x23
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 26]*x11*x22
        +coeff[ 27]    *x22        *x51
        +coeff[ 28]*x12        *x41
        +coeff[ 29]*x11*x21    *x42
        +coeff[ 30]    *x24
        +coeff[ 31]    *x22*x32
        +coeff[ 32]*x11*x21*x31*x41
        +coeff[ 33]*x11*x23
        +coeff[ 34]*x13    *x31
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 35]*x11        *x41
        +coeff[ 36]    *x21        *x51
        +coeff[ 37]        *x31    *x51
        +coeff[ 38]            *x43
        +coeff[ 39]        *x31*x42
        +coeff[ 40]        *x32*x41
        +coeff[ 41]*x12    *x31
        +coeff[ 42]*x11*x21        *x51
        +coeff[ 43]*x12*x22
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 44]*x11    *x31
        +coeff[ 45]*x11        *x42
        +coeff[ 46]            *x42*x51
        +coeff[ 47]    *x21*x32
        +coeff[ 48]*x11    *x31*x41
        +coeff[ 49]        *x31*x41*x51
        +coeff[ 50]*x12*x21
        +coeff[ 51]    *x23    *x41
        +coeff[ 52]*x11*x22    *x41
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 53]    *x22    *x41*x51
        +coeff[ 54]*x12        *x42
        +coeff[ 55]*x11*x21*x32
        +coeff[ 56]*x12    *x31*x41
        +coeff[ 57]    *x24    *x41
        +coeff[ 58]*x11*x23    *x41
        +coeff[ 59]*x11            *x51
        +coeff[ 60]        *x33
        +coeff[ 61]    *x21    *x41*x51
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 62]*x11    *x32
        +coeff[ 63]        *x32    *x51
        +coeff[ 64]    *x23*x31
        +coeff[ 65]*x12            *x51
        +coeff[ 66]    *x22*x31    *x51
        +coeff[ 67]    *x22    *x43
        +coeff[ 68]*x11*x21    *x41*x51
        +coeff[ 69]    *x22*x31*x42
        +coeff[ 70]    *x22*x32*x41
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 71]*x13*x21
        +coeff[ 72]*x11*x23*x31
        +coeff[ 73]*x12*x22    *x41
        +coeff[ 74]*x11*x21*x31*x43
        +coeff[ 75]*x11*x23    *x42
        +coeff[ 76]*x11*x21*x32*x42
        +coeff[ 77]*x11*x21    *x45
        +coeff[ 78]*x11*x21*x31*x44
        +coeff[ 79]    *x22    *x44*x51
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 80]    *x24*x31*x42
        +coeff[ 81]    *x24    *x44
        +coeff[ 82]*x12*x24*x31
        +coeff[ 83]    *x21*x31*x42
        +coeff[ 84]                *x53
        +coeff[ 85]*x11*x22*x31
        +coeff[ 86]    *x21*x31*x43
        +coeff[ 87]*x12*x21    *x41
        +coeff[ 88]*x11*x21*x31    *x51
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 89]    *x24*x31
        +coeff[ 90]    *x22*x33
        +coeff[ 91]    *x21    *x45
        +coeff[ 92]*x11*x21*x32*x41
        +coeff[ 93]    *x22*x31*x41*x51
        +coeff[ 94]    *x22    *x44
        +coeff[ 95]    *x21*x32*x43
        +coeff[ 96]    *x24    *x42
        ;

    return v_y_l5p77_q1en                            ;
}
float p_l5p77_q1en                            (float *x,int m){
    int ncoeff= 51;
    float avdat=  0.5847551E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25772E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 52]={
        -0.44612591E-02,-0.14153930E-02, 0.23939665E-02, 0.35504423E-01,
        -0.62013390E-02, 0.86985286E-02,-0.33833012E-02, 0.37311856E-02,
        -0.52377414E-02,-0.33952983E-03,-0.24103001E-02,-0.19754607E-02,
        -0.24613829E-02, 0.87307335E-03,-0.42224475E-03, 0.31692837E-03,
         0.45609515E-03,-0.94859506E-03, 0.34289968E-02, 0.21641233E-03,
         0.32890221E-03,-0.18364046E-03, 0.18169245E-03,-0.37084395E-03,
         0.13875937E-02, 0.13990713E-03, 0.40516601E-03,-0.48800246E-03,
         0.24402094E-03,-0.13869136E-02, 0.21509952E-02,-0.15005608E-03,
        -0.11133562E-02,-0.26134733E-03, 0.88748609E-03,-0.52778363E-04,
        -0.32677560E-03,-0.29098379E-04, 0.11162267E-03, 0.55249591E-04,
        -0.64535401E-04, 0.61491271E-04,-0.28620142E-03, 0.49954955E-03,
         0.53723715E-03, 0.60906797E-04, 0.28181271E-03,-0.82507700E-04,
         0.14735512E-03, 0.75435943E-04, 0.14849620E-03,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_p_l5p77_q1en                            =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]            *x42
        +coeff[  7]*x11*x21
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[  8]    *x22    *x41
        +coeff[  9]*x11
        +coeff[ 10]        *x31*x41
        +coeff[ 11]    *x22*x31
        +coeff[ 12]*x11*x21    *x41
        +coeff[ 13]    *x21    *x41
        +coeff[ 14]        *x32
        +coeff[ 15]                *x52
        +coeff[ 16]*x12
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 17]*x11*x21*x31
        +coeff[ 18]    *x22    *x42
        +coeff[ 19]    *x24*x31*x41
        +coeff[ 20]    *x21*x31
        +coeff[ 21]            *x41*x51
        +coeff[ 22]*x11        *x41
        +coeff[ 23]    *x22        *x51
        +coeff[ 24]*x11*x21    *x42
        +coeff[ 25]*x11*x23*x31*x41
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 26]    *x23
        +coeff[ 27]    *x21    *x42
        +coeff[ 28]*x11*x22
        +coeff[ 29]    *x24
        +coeff[ 30]    *x22*x31*x41
        +coeff[ 31]*x11*x21        *x51
        +coeff[ 32]*x11*x23
        +coeff[ 33]*x12        *x41
        +coeff[ 34]*x11*x21*x31*x41
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 35]    *x23*x31*x41
        +coeff[ 36]*x12*x22
        +coeff[ 37]            *x45
        +coeff[ 38]    *x24*x32
        +coeff[ 39]    *x21        *x51
        +coeff[ 40]        *x31    *x51
        +coeff[ 41]*x11    *x31
        +coeff[ 42]    *x21*x31*x41
        +coeff[ 43]        *x31*x42
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 44]            *x43
        +coeff[ 45]            *x42*x51
        +coeff[ 46]    *x22*x32
        +coeff[ 47]*x12    *x31
        +coeff[ 48]*x11*x21*x32
        +coeff[ 49]    *x22*x32*x41
        +coeff[ 50]        *x34*x41
        ;

    return v_p_l5p77_q1en                            ;
}
float l_l5p77_q1en                            (float *x,int m){
    int ncoeff= 49;
    float avdat= -0.9177745E-03;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25772E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 50]={
         0.26036895E-03, 0.18447302E-02, 0.57269596E-02,-0.11647895E-03,
        -0.18957240E-02,-0.14118374E-03,-0.12122416E-02,-0.68314030E-03,
         0.13591786E-03, 0.18183996E-03,-0.25370048E-03, 0.97302531E-04,
        -0.17047727E-03,-0.29524143E-04,-0.55305016E-04, 0.31545360E-03,
         0.12635181E-04,-0.10283304E-04, 0.20534806E-04, 0.41600171E-04,
        -0.39679504E-04, 0.86952277E-04, 0.10569551E-03,-0.46362024E-05,
        -0.64720902E-04, 0.93915587E-05,-0.14694689E-04, 0.39090369E-05,
        -0.44764315E-05, 0.12825153E-04,-0.20515134E-04, 0.69949194E-04,
         0.22667042E-04,-0.12830175E-04,-0.14112372E-03,-0.38614464E-04,
         0.18084278E-03,-0.10931842E-03, 0.52784126E-04,-0.28884151E-04,
         0.12835687E-03,-0.87133434E-04, 0.90881898E-04, 0.16078362E-04,
         0.12119412E-04, 0.65800468E-05,-0.71966874E-05, 0.21245029E-04,
        -0.19008054E-04,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_l_l5p77_q1en                            =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]        *x31*x41
        +coeff[  6]            *x42
        +coeff[  7]    *x22    *x41
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[  8]            *x41*x51
        +coeff[  9]*x11*x21
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]    *x21    *x41
        +coeff[ 12]    *x22*x31
        +coeff[ 13]    *x21
        +coeff[ 14]*x11*x21*x31
        +coeff[ 15]    *x22    *x42
        +coeff[ 16]    *x24*x31*x41
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 17]*x11
        +coeff[ 18]    *x21*x31
        +coeff[ 19]    *x23
        +coeff[ 20]    *x21    *x42
        +coeff[ 21]            *x43
        +coeff[ 22]*x11*x21    *x42
        +coeff[ 23]*x11        *x43
        +coeff[ 24]    *x22*x31*x42
        +coeff[ 25]*x11*x23*x31*x41
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 26]        *x32
        +coeff[ 27]        *x31    *x51
        +coeff[ 28]                *x52
        +coeff[ 29]*x12
        +coeff[ 30]    *x21*x31*x41
        +coeff[ 31]        *x31*x42
        +coeff[ 32]*x11*x22
        +coeff[ 33]*x12        *x41
        +coeff[ 34]    *x24
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 35]    *x23    *x41
        +coeff[ 36]    *x22*x31*x41
        +coeff[ 37]*x11*x23
        +coeff[ 38]*x11*x21*x31*x41
        +coeff[ 39]*x12*x22
        +coeff[ 40]    *x24    *x41
        +coeff[ 41]    *x22    *x43
        +coeff[ 42]*x11*x23    *x41
        +coeff[ 43]*x11        *x41
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 44]        *x32*x41
        +coeff[ 45]    *x22        *x51
        +coeff[ 46]            *x41*x52
        +coeff[ 47]    *x22*x32
        +coeff[ 48]*x11*x22    *x41
        ;

    return v_l_l5p77_q1en                            ;
}
float x_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 52;
    float avdat=  0.1030081E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 53]={
        -0.33856413E-03, 0.43503325E-02, 0.11492133E+00, 0.29791938E-02,
         0.20091338E-01,-0.16988031E-04, 0.46920585E-02, 0.73638731E-02,
        -0.17507772E-02, 0.17822186E-02, 0.53844550E-02,-0.70402212E-02,
        -0.27623530E-04, 0.31997974E-05,-0.66244189E-03,-0.14634107E-02,
         0.67007568E-04, 0.95662603E-03, 0.39676963E-02,-0.51887766E-02,
         0.10958033E-04,-0.72087720E-02, 0.39663399E-03,-0.71869150E-03,
        -0.96356252E-03, 0.16284148E-02,-0.44298405E-02,-0.22325530E-02,
         0.99564881E-04,-0.59727990E-05,-0.73963136E-04, 0.52065583E-03,
         0.69493585E-03,-0.10679403E-02,-0.15123088E-02, 0.55216486E-03,
        -0.10882330E-02,-0.70801766E-05,-0.15331858E-02,-0.14868940E-04,
        -0.58186450E-03,-0.85591426E-04, 0.72952462E-04, 0.40847101E-03,
         0.23730814E-03,-0.21423367E-03,-0.34025966E-03, 0.21683532E-02,
         0.21207125E-02,-0.21677099E-03, 0.34719074E-03, 0.10737858E-02,
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
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_x_l5p77_q1ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]    *x21        *x51
        +coeff[  4]    *x21    *x41
        +coeff[  5]    *x21*x31    *x52
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21*x31
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[  8]            *x41
        +coeff[  9]*x11    *x31
        +coeff[ 10]    *x23
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]*x13*x22
        +coeff[ 13]    *x23*x31*x41
        +coeff[ 14]        *x31
        +coeff[ 15]    *x22
        +coeff[ 16]*x13
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 17]*x12*x21
        +coeff[ 18]*x11*x22
        +coeff[ 19]    *x21*x31*x41
        +coeff[ 20]*x13*x21
        +coeff[ 21]    *x23    *x41
        +coeff[ 22]*x11            *x51
        +coeff[ 23]*x11*x21
        +coeff[ 24]    *x21*x32
        +coeff[ 25]    *x22    *x41
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 26]*x11*x22    *x41
        +coeff[ 27]    *x23*x31
        +coeff[ 28]*x13        *x42
        +coeff[ 29]*x11    *x33*x41
        +coeff[ 30]                *x51
        +coeff[ 31]            *x42
        +coeff[ 32]*x11*x21    *x41
        +coeff[ 33]*x11    *x31*x41
        +coeff[ 34]*x11        *x42
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 35]    *x22*x31
        +coeff[ 36]*x12*x21    *x41
        +coeff[ 37]*x12    *x31*x41
        +coeff[ 38]*x11*x22*x31
        +coeff[ 39]        *x33*x41
        +coeff[ 40]    *x21    *x41*x53
        +coeff[ 41]*x12
        +coeff[ 42]        *x32
        +coeff[ 43]        *x31*x41
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 44]*x11*x21*x31
        +coeff[ 45]*x11    *x32
        +coeff[ 46]*x12*x21*x31
        +coeff[ 47]    *x21*x31*x42
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]    *x21*x31    *x53
        +coeff[ 50]*x12*x21*x32*x41
        +coeff[ 51]    *x23*x32*x41
        ;

    return v_x_l5p77_q1ex                            ;
}
float t_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.2191480E-03;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.19173509E-03,-0.63087870E-02,-0.60998476E-02, 0.49153129E-02,
         0.22201587E-02, 0.86188029E-05,-0.82357525E-04,-0.73989801E-03,
        -0.16440164E-04,-0.93911140E-05,-0.45035495E-05,-0.39451956E-03,
        -0.30978629E-03, 0.18311139E-02, 0.10543996E-02, 0.23005916E-03,
         0.10422883E-02,-0.17651013E-05,-0.18773357E-02,-0.11178919E-04,
        -0.19727470E-02,-0.29017065E-04, 0.82760762E-05, 0.66583772E-04,
        -0.13651642E-03,-0.14105179E-03, 0.38561440E-03, 0.23403096E-04,
         0.19029644E-03, 0.81410445E-03, 0.43412193E-03,-0.11562252E-02,
        -0.12040631E-02,-0.29096625E-04,-0.74095010E-05,-0.16381953E-04,
         0.10179790E-03, 0.12057424E-03,-0.17755132E-03, 0.14771552E-03,
        -0.20180014E-03,-0.30275289E-03, 0.15026855E-03,-0.96407435E-04,
        -0.39650366E-03,-0.22646521E-03, 0.21162061E-04, 0.37637432E-03,
         0.53382741E-03, 0.65025342E-04,
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
    float x52 = x51*x5;

//                 function

    float v_t_l5p77_q1ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x21        *x51
        +coeff[  5]*x12        *x41
        +coeff[  6]*x12*x21*x31
        +coeff[  7]    *x23*x31
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[  8]    *x21*x33
        +coeff[  9]*x13        *x41
        +coeff[ 10]*x13*x22*x31
        +coeff[ 11]            *x41
        +coeff[ 12]    *x22
        +coeff[ 13]    *x21*x31
        +coeff[ 14]*x11        *x41
        +coeff[ 15]*x11            *x51
        +coeff[ 16]    *x23
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 17]        *x33
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]*x13    *x31
        +coeff[ 20]    *x23    *x41
        +coeff[ 21]*x13*x22
        +coeff[ 22]*x12*x21*x31*x41
        +coeff[ 23]    *x23*x31*x41
        +coeff[ 24]        *x31
        +coeff[ 25]*x11*x21
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 26]*x11    *x31
        +coeff[ 27]*x13
        +coeff[ 28]*x12*x21
        +coeff[ 29]*x11*x22
        +coeff[ 30]    *x22    *x41
        +coeff[ 31]    *x21*x31*x41
        +coeff[ 32]*x11*x22    *x41
        +coeff[ 33]*x13        *x42
        +coeff[ 34]    *x21*x32    *x52
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 35]                *x51
        +coeff[ 36]            *x42
        +coeff[ 37]    *x22*x31
        +coeff[ 38]    *x21*x32
        +coeff[ 39]*x11*x21    *x41
        +coeff[ 40]*x11    *x31*x41
        +coeff[ 41]*x11        *x42
        +coeff[ 42]    *x21    *x41*x51
        +coeff[ 43]    *x21        *x52
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 44]*x11*x22*x31
        +coeff[ 45]*x12*x21    *x41
        +coeff[ 46]*x12    *x31*x41
        +coeff[ 47]    *x21    *x43
        +coeff[ 48]    *x23*x31*x42
        +coeff[ 49]        *x31*x41
        ;

    return v_t_l5p77_q1ex                            ;
}
float y_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1199349E-01;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.12116483E-01, 0.17086735E+00, 0.29141398E-01,-0.18324599E-01,
         0.26657723E-01, 0.11845203E-01,-0.16037821E-02,-0.44924407E-02,
        -0.10762360E-01,-0.79489481E-02,-0.18423496E-01,-0.90165213E-02,
        -0.10532123E-02, 0.31484666E-02,-0.27426304E-02, 0.14233699E-02,
         0.11072924E-02,-0.70225829E-02,-0.35166501E-02, 0.12372417E-02,
        -0.14791369E-02,-0.72438503E-03, 0.11830587E-01, 0.83027026E-02,
         0.76903601E-03, 0.13150476E-02, 0.90231607E-03,-0.12948896E-02,
        -0.10988686E-02, 0.50632781E-02,-0.39477190E-02, 0.34891688E-02,
        -0.34398641E-02,-0.28817769E-05, 0.29827325E-03, 0.22317436E-03,
         0.20166531E-02, 0.19891914E-02,-0.11339659E-02, 0.72260742E-03,
        -0.38850080E-03,-0.56718854E-03, 0.13954008E-02,-0.12950744E-02,
        -0.34999117E-03, 0.48263196E-03,-0.20798581E-03,-0.25018057E-03,
         0.39778874E-03, 0.20729996E-03, 0.61664934E-03, 0.70217432E-03,
         0.49355098E-04, 0.87608169E-04,-0.10585998E-03,-0.44432152E-04,
         0.13685727E-03,-0.66755019E-03,-0.74218573E-04,-0.44655288E-03,
         0.54081820E-03, 0.27323927E-03, 0.34962015E-03, 0.28969671E-03,
         0.21785772E-02,-0.18439612E-03, 0.23901451E-02, 0.90783346E-03,
         0.81204707E-04, 0.85195870E-03,-0.16424428E-03, 0.26973258E-03,
         0.34967252E-04,-0.20963191E-03,-0.68196256E-04,-0.11658897E-03,
        -0.14159172E-02, 0.13181606E-03,-0.54365862E-03,-0.36774919E-03,
         0.79747871E-03,-0.25361613E-03,-0.19591952E-03,-0.10801711E-02,
        -0.18748422E-02, 0.25777187E-03,-0.35751832E-02,-0.43584909E-02,
        -0.23048541E-03,-0.30282215E-03,-0.13238007E-02,-0.13962633E-02,
        -0.42178563E-03, 0.63591672E-03,-0.82226179E-04,-0.34249260E-04,
         0.76162905E-04,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_q1ex                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]        *x31
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]*x11*x21
        +coeff[  6]    *x21    *x42
        +coeff[  7]    *x21
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[  8]            *x42
        +coeff[  9]        *x31*x41
        +coeff[ 10]    *x22    *x41
        +coeff[ 11]*x11*x21    *x41
        +coeff[ 12]*x11
        +coeff[ 13]    *x21    *x41
        +coeff[ 14]            *x41*x51
        +coeff[ 15]*x12
        +coeff[ 16]                *x52
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 17]    *x22*x31
        +coeff[ 18]*x11*x21*x31
        +coeff[ 19]    *x21*x31
        +coeff[ 20]        *x32
        +coeff[ 21]        *x31    *x51
        +coeff[ 22]    *x22    *x42
        +coeff[ 23]    *x22*x31*x41
        +coeff[ 24]*x11        *x41
        +coeff[ 25]    *x23
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 26]*x11*x22
        +coeff[ 27]    *x22        *x51
        +coeff[ 28]*x12        *x41
        +coeff[ 29]*x11*x21    *x42
        +coeff[ 30]    *x24
        +coeff[ 31]*x11*x21*x31*x41
        +coeff[ 32]*x11*x23
        +coeff[ 33]*x11*x23        *x51
        +coeff[ 34]*x11    *x31
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 35]    *x21        *x51
        +coeff[ 36]            *x43
        +coeff[ 37]        *x31*x42
        +coeff[ 38]    *x21*x31*x41
        +coeff[ 39]        *x32*x41
        +coeff[ 40]*x12    *x31
        +coeff[ 41]*x11*x21        *x51
        +coeff[ 42]    *x22*x32
        +coeff[ 43]*x12*x22
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 44]*x11        *x42
        +coeff[ 45]            *x42*x51
        +coeff[ 46]    *x21*x32
        +coeff[ 47]*x11    *x31*x41
        +coeff[ 48]        *x31*x41*x51
        +coeff[ 49]*x12*x21
        +coeff[ 50]    *x22    *x41*x51
        +coeff[ 51]*x11*x21*x32
        +coeff[ 52]*x11            *x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 53]        *x33
        +coeff[ 54]    *x21    *x41*x51
        +coeff[ 55]*x11    *x32
        +coeff[ 56]            *x41*x52
        +coeff[ 57]    *x23    *x41
        +coeff[ 58]*x12            *x51
        +coeff[ 59]*x11*x22    *x41
        +coeff[ 60]*x12        *x42
        +coeff[ 61]    *x22*x31    *x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 62]*x12    *x31*x41
        +coeff[ 63]*x11*x21    *x41*x51
        +coeff[ 64]    *x24    *x41
        +coeff[ 65]*x13*x21
        +coeff[ 66]*x11*x23    *x41
        +coeff[ 67]*x11*x23*x31
        +coeff[ 68]        *x34    *x51
        +coeff[ 69]*x12*x22    *x41
        +coeff[ 70]    *x22    *x42*x53
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 71]*x12*x24*x31
        +coeff[ 72]        *x31    *x52
        +coeff[ 73]    *x23*x31
        +coeff[ 74]                *x53
        +coeff[ 75]*x11*x22*x31
        +coeff[ 76]    *x22*x31*x42
        +coeff[ 77]*x11*x21*x31    *x51
        +coeff[ 78]    *x22*x32*x41
        +coeff[ 79]    *x22    *x42*x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 80]    *x24*x31
        +coeff[ 81]    *x22*x31*x41*x51
        +coeff[ 82]*x11*x21    *x42*x51
        +coeff[ 83]*x11*x23    *x42
        +coeff[ 84]    *x22    *x45
        +coeff[ 85]*x12*x22*x32
        +coeff[ 86]    *x24    *x44
        +coeff[ 87]    *x24*x31*x43
        +coeff[ 88]*x11*x22*x32*x41*x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 89]*x13*x21*x31*x42
        +coeff[ 90]    *x24*x32*x42
        +coeff[ 91]*x11*x23*x31*x43
        +coeff[ 92]*x11*x21*x34*x42
        +coeff[ 93]    *x22*x32*x45
        +coeff[ 94]            *x44
        +coeff[ 95]    *x21*x31    *x51
        +coeff[ 96]    *x21    *x43
        ;

    return v_y_l5p77_q1ex                            ;
}
float p_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 46;
    float avdat=  0.6160242E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 47]={
        -0.54003615E-02, 0.98044639E-02, 0.72380595E-01,-0.94280513E-02,
         0.13594504E-01,-0.22021141E-02, 0.60144593E-02, 0.60799276E-03,
        -0.22299474E-02,-0.52890470E-02,-0.86255549E-02,-0.40193922E-02,
         0.36054496E-02, 0.38810274E-04,-0.54733176E-03, 0.14081155E-02,
        -0.68606331E-03,-0.38561518E-02,-0.47574149E-03, 0.66860771E-03,
        -0.32883077E-02, 0.73636771E-03,-0.15383859E-02, 0.54120482E-02,
         0.52280223E-03,-0.82632998E-03, 0.30644526E-03,-0.67110523E-03,
         0.41855083E-03,-0.21115267E-02, 0.50317333E-03,-0.34453499E-03,
        -0.17730450E-02,-0.44917219E-03,-0.80176024E-03, 0.20707180E-02,
         0.16762057E-03,-0.54207171E-03, 0.79924241E-04, 0.11261086E-03,
         0.11166939E-03, 0.85992413E-03, 0.86012407E-03,-0.13997426E-03,
         0.13041710E-02, 0.28082583E-03,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_q1ex                            =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]            *x41*x51
        +coeff[  6]*x11*x21
        +coeff[  7]    *x23
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[  8]    *x21
        +coeff[  9]            *x42
        +coeff[ 10]    *x22    *x41
        +coeff[ 11]*x11*x21    *x41
        +coeff[ 12]    *x22*x31*x41
        +coeff[ 13]        *x33*x41
        +coeff[ 14]*x11
        +coeff[ 15]    *x21    *x41
        +coeff[ 16]        *x32
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[ 17]        *x31*x41
        +coeff[ 18]        *x31    *x51
        +coeff[ 19]                *x52
        +coeff[ 20]    *x22*x31
        +coeff[ 21]*x12
        +coeff[ 22]*x11*x21*x31
        +coeff[ 23]    *x22    *x42
        +coeff[ 24]    *x21*x31
        +coeff[ 25]    *x22        *x51
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[ 26]*x11        *x41
        +coeff[ 27]    *x21    *x42
        +coeff[ 28]*x11*x22
        +coeff[ 29]    *x24
        +coeff[ 30]    *x22*x32
        +coeff[ 31]*x11*x21        *x51
        +coeff[ 32]*x11*x23
        +coeff[ 33]*x12        *x41
        +coeff[ 34]    *x23*x31*x41
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[ 35]*x11*x21    *x42
        +coeff[ 36]    *x21        *x53
        +coeff[ 37]*x12*x22
        +coeff[ 38]            *x45
        +coeff[ 39]*x11*x23*x31*x41
        +coeff[ 40]*x11    *x31
        +coeff[ 41]        *x31*x42
        +coeff[ 42]            *x43
        +coeff[ 43]*x12    *x31
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[ 44]*x11*x21*x31*x41
        +coeff[ 45]        *x34*x41
        ;

    return v_p_l5p77_q1ex                            ;
}
float l_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 90;
    float avdat= -0.1734428E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.10988429E-02, 0.56253476E-02,-0.22384019E-02, 0.58745274E-04,
        -0.64355605E-04, 0.31490577E-03,-0.89249259E-03,-0.42425916E-02,
         0.97112235E-03, 0.21568265E-03, 0.67590481E-05,-0.42837235E-03,
         0.58917449E-05,-0.21090026E-02, 0.45059231E-03, 0.55045116E-03,
        -0.14233701E-04,-0.14510991E-03,-0.77069009E-03,-0.51025268E-06,
         0.82401010E-06, 0.16815179E-03,-0.22237041E-05,-0.28081588E-05,
        -0.45028728E-05, 0.50612484E-05, 0.24074229E-05, 0.26299590E-05,
        -0.39406877E-05,-0.69176003E-05, 0.18397116E-02, 0.10520311E-03,
         0.91670981E-05, 0.11121243E-02,-0.39563423E-04,-0.65125008E-04,
         0.58188194E-04, 0.16491862E-03, 0.60375902E-03, 0.44305826E-03,
        -0.51236552E-05,-0.10621751E-03,-0.28170642E-04, 0.86592117E-04,
         0.10118484E-03,-0.83787825E-04, 0.90173460E-04,-0.73166600E-04,
        -0.29501511E-03,-0.22335621E-03, 0.21414357E-03, 0.54228038E-03,
         0.38194709E-03,-0.67550054E-05,-0.12127724E-04,-0.73058487E-04,
        -0.14832414E-03, 0.49620150E-04,-0.27384673E-04, 0.11643446E-04,
         0.67937552E-04,-0.15926638E-03,-0.81521670E-04,-0.57474652E-04,
        -0.55606110E-03, 0.12937858E-03,-0.19347257E-03,-0.62089998E-05,
         0.10440312E-03, 0.12616974E-04,-0.38695256E-04,-0.98495060E-04,
        -0.24955074E-04,-0.21869166E-04, 0.36137109E-04,-0.31229368E-03,
        -0.15231170E-03, 0.10618510E-03,-0.11357298E-03,-0.20254913E-03,
        -0.42647684E-05,-0.12342843E-04, 0.71846716E-05,-0.29055553E-04,
         0.27539732E-04,-0.99636636E-04,-0.18071396E-04, 0.18839493E-04,
        -0.13967471E-04,-0.11690571E-04,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_l_l5p77_q1ex                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x22
        +coeff[  3]    *x21*x31
        +coeff[  4]        *x32
        +coeff[  5]    *x21    *x41
        +coeff[  6]        *x31*x41
        +coeff[  7]            *x42
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[  8]            *x41*x51
        +coeff[  9]*x11*x21
        +coeff[ 10]*x11    *x31
        +coeff[ 11]    *x22*x31
        +coeff[ 12]        *x33
        +coeff[ 13]    *x22    *x41
        +coeff[ 14]        *x31*x42
        +coeff[ 15]            *x43
        +coeff[ 16]        *x31    *x52
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 17]*x11*x21*x31
        +coeff[ 18]*x11*x21    *x41
        +coeff[ 19]*x12    *x31
        +coeff[ 20]        *x32    *x52
        +coeff[ 21]    *x24*x31
        +coeff[ 22]    *x22*x33
        +coeff[ 23]        *x34*x41
        +coeff[ 24]        *x33*x42
        +coeff[ 25]    *x22*x31    *x52
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 26]        *x33    *x52
        +coeff[ 27]*x12    *x33
        +coeff[ 28]    *x24*x33
        +coeff[ 29]    *x22*x33    *x52
        +coeff[ 30]        *x31
        +coeff[ 31]        *x31    *x51
        +coeff[ 32]                *x53
        +coeff[ 33]    *x22    *x42
        +coeff[ 34]    *x21
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 35]                *x52
        +coeff[ 36]*x11        *x41
        +coeff[ 37]    *x22        *x51
        +coeff[ 38]    *x22*x31*x41
        +coeff[ 39]*x11*x21    *x42
        +coeff[ 40]    *x21    *x44
        +coeff[ 41]                *x51
        +coeff[ 42]    *x21        *x51
        +coeff[ 43]    *x23
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 44]            *x42*x51
        +coeff[ 45]            *x41*x52
        +coeff[ 46]*x11*x21        *x51
        +coeff[ 47]*x12        *x41
        +coeff[ 48]    *x24
        +coeff[ 49]*x11*x23
        +coeff[ 50]*x11*x21*x31*x41
        +coeff[ 51]    *x24    *x41
        +coeff[ 52]*x11*x23    *x41
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 53]*x11            *x51
        +coeff[ 54]*x12
        +coeff[ 55]    *x21*x31*x41
        +coeff[ 56]    *x21    *x42
        +coeff[ 57]*x11*x22
        +coeff[ 58]*x11        *x42
        +coeff[ 59]*x12            *x51
        +coeff[ 60]    *x22*x32
        +coeff[ 61]    *x23    *x41
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 62]*x11*x22    *x41
        +coeff[ 63]*x12*x22
        +coeff[ 64]    *x22    *x43
        +coeff[ 65]*x11*x23*x31
        +coeff[ 66]    *x24*x31*x42
        +coeff[ 67]*x11
        +coeff[ 68]        *x32*x41
        +coeff[ 69]        *x31*x41*x51
        +coeff[ 70]    *x23*x31
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 71]            *x44
        +coeff[ 72]    *x22*x31    *x51
        +coeff[ 73]    *x22        *x52
        +coeff[ 74]*x12        *x42
        +coeff[ 75]    *x22*x31*x42
        +coeff[ 76]*x11*x21    *x43
        +coeff[ 77]*x12*x22    *x41
        +coeff[ 78]    *x24*x32*x41
        +coeff[ 79]*x11*x23*x31*x42
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 80]        *x32    *x51
        +coeff[ 81]*x11    *x31*x41
        +coeff[ 82]*x12*x21
        +coeff[ 83]        *x32*x42
        +coeff[ 84]    *x21    *x43
        +coeff[ 85]        *x31*x43
        +coeff[ 86]*x11*x22*x31
        +coeff[ 87]*x11*x21*x32
        +coeff[ 88]*x11*x21*x31    *x51
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 89]*x11*x21    *x41*x51
        ;

    return v_l_l5p77_q1ex                            ;
}
float x_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.1319563E+01;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
         0.34062919E-03,-0.23021912E-01, 0.18848316E+00, 0.10512160E-01,
         0.62990844E-01,-0.17079886E-03,-0.11258640E-02,-0.17204927E-03,
        -0.33806890E-03,-0.53733587E-02, 0.52980636E-02, 0.14263028E-01,
         0.23173526E-01, 0.15442866E-01,-0.24356155E-01, 0.16049891E-03,
         0.18666520E-02,-0.19800647E-02, 0.20899375E-02,-0.43782652E-02,
         0.11341004E-01,-0.17199323E-01, 0.85337888E-04,-0.23865797E-01,
        -0.21368079E-02, 0.28192503E-02,-0.30030175E-02, 0.54104105E-02,
        -0.14740811E-01,-0.76095252E-02, 0.47422101E-03,-0.26400242E-03,
         0.14497383E-02,-0.77383069E-03, 0.23157443E-02,-0.31511770E-02,
        -0.47379956E-02, 0.18068737E-02,-0.31571628E-02, 0.10303967E-04,
        -0.49816398E-02,-0.86855325E-04, 0.11400218E-02, 0.79746393E-03,
        -0.59395714E-03,-0.26663151E-03, 0.61718078E-03, 0.49703852E-02,
         0.39604111E-02, 0.35451469E-02,-0.56246059E-04, 0.71972725E-03,
        -0.38592680E-03, 0.10736925E-03, 0.17733639E-03,-0.13104896E-03,
         0.19269776E-03,-0.48166697E-03,-0.72620198E-03, 0.39501321E-02,
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
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_x_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]    *x21        *x51
        +coeff[  4]    *x21    *x41
        +coeff[  5]*x13        *x41
        +coeff[  6]*x12*x21*x31
        +coeff[  7]    *x21*x31    *x52
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[  8]    *x21*x33
        +coeff[  9]            *x41
        +coeff[ 10]*x11    *x31
        +coeff[ 11]*x11        *x41
        +coeff[ 12]    *x21*x31
        +coeff[ 13]    *x23
        +coeff[ 14]    *x21    *x42
        +coeff[ 15]*x13*x22
        +coeff[ 16]    *x21*x33*x41
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 17]        *x31
        +coeff[ 18]*x11            *x51
        +coeff[ 19]    *x22
        +coeff[ 20]*x11*x22
        +coeff[ 21]    *x21*x31*x41
        +coeff[ 22]*x13*x21
        +coeff[ 23]    *x23    *x41
        +coeff[ 24]*x11*x21
        +coeff[ 25]*x12*x21
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 26]    *x21*x32
        +coeff[ 27]    *x22    *x41
        +coeff[ 28]*x11*x22    *x41
        +coeff[ 29]    *x23*x31
        +coeff[ 30]*x13        *x42
        +coeff[ 31]                *x51
        +coeff[ 32]            *x42
        +coeff[ 33]    *x21        *x52
        +coeff[ 34]*x11*x21    *x41
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]*x11        *x42
        +coeff[ 37]    *x22*x31
        +coeff[ 38]*x12*x21    *x41
        +coeff[ 39]*x12    *x31*x41
        +coeff[ 40]*x11*x22*x31
        +coeff[ 41]        *x33*x41
        +coeff[ 42]        *x31*x41
        +coeff[ 43]*x11*x21*x31
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 44]*x11    *x32
        +coeff[ 45]*x14
        +coeff[ 46]    *x21    *x42*x51
        +coeff[ 47]    *x21    *x43
        +coeff[ 48]    *x21*x32*x42
        +coeff[ 49]    *x21*x31*x43
        +coeff[ 50]*x12*x21*x31*x42
        +coeff[ 51]    *x21*x31*x42*x52
        +coeff[ 52]    *x23*x31*x42
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 53]            *x41*x51
        +coeff[ 54]        *x32
        +coeff[ 55]*x11            *x52
        +coeff[ 56]*x12        *x41
        +coeff[ 57]    *x21*x31    *x51
        +coeff[ 58]    *x21    *x41*x51
        +coeff[ 59]    *x21*x31*x42
        ;

    return v_x_l5p77_q2ex                            ;
}
float t_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.5766578E+00;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.63234655E-03,-0.11202579E-01, 0.55861399E-01,-0.30467319E-02,
         0.22987105E-01, 0.30485024E-02,-0.93670169E-04,-0.68493740E-03,
        -0.29334591E-02,-0.55467510E-04,-0.30587727E-03,-0.63191185E-03,
        -0.18282143E-02, 0.86160414E-02, 0.50515146E-02, 0.55339122E-02,
        -0.91661150E-02,-0.67231573E-04,-0.40832259E-04,-0.43304724E-04,
         0.18821404E-02, 0.86258043E-03, 0.58149712E-04, 0.40850998E-02,
        -0.56832572E-02,-0.90597039E-02,-0.31391170E-03,-0.14863664E-03,
         0.95849816E-03,-0.89790171E-03,-0.14499655E-02,-0.50950469E-02,
        -0.15793089E-03,-0.89145738E-04,-0.14370163E-03,-0.23740334E-03,
         0.42420934E-03, 0.57620159E-03, 0.73895033E-03, 0.90090994E-03,
        -0.94869698E-03,-0.33393005E-03, 0.67639071E-05,-0.16114464E-02,
        -0.78536541E-03, 0.11978844E-02, 0.16132119E-02, 0.65187211E-04,
         0.70877191E-04, 0.15715024E-03,
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
    float x52 = x51*x5;

//                 function

    float v_t_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]    *x22
        +coeff[  4]    *x21    *x41
        +coeff[  5]    *x21        *x51
        +coeff[  6]            *x43
        +coeff[  7]*x12*x21*x31
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[  8]    *x23*x31
        +coeff[  9]*x13        *x41
        +coeff[ 10]*x13*x22*x31
        +coeff[ 11]        *x31
        +coeff[ 12]            *x41
        +coeff[ 13]    *x21*x31
        +coeff[ 14]*x11        *x41
        +coeff[ 15]    *x23
        +coeff[ 16]    *x21    *x42
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 17]*x11    *x33
        +coeff[ 18]*x13*x22
        +coeff[ 19]*x12*x21*x31*x41
        +coeff[ 20]*x11    *x31
        +coeff[ 21]*x11            *x51
        +coeff[ 22]*x13
        +coeff[ 23]*x11*x22
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]    *x23    *x41
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 26]*x13*x22    *x41
        +coeff[ 27]                *x51
        +coeff[ 28]*x12*x21
        +coeff[ 29]    *x21*x32
        +coeff[ 30]*x11        *x42
        +coeff[ 31]*x11*x22    *x41
        +coeff[ 32]*x11*x21    *x42
        +coeff[ 33]*x13    *x31*x41
        +coeff[ 34]*x12
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 35]*x11*x21
        +coeff[ 36]        *x31*x41
        +coeff[ 37]            *x42
        +coeff[ 38]*x11*x21    *x41
        +coeff[ 39]    *x22    *x41
        +coeff[ 40]*x11    *x31*x41
        +coeff[ 41]    *x21        *x52
        +coeff[ 42]*x13    *x31
        +coeff[ 43]*x11*x22*x31
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 44]*x12*x21    *x41
        +coeff[ 45]    *x21    *x43
        +coeff[ 46]*x12*x21*x31*x42
        +coeff[ 47]        *x32
        +coeff[ 48]            *x41*x51
        +coeff[ 49]*x11*x21*x31
        ;

    return v_t_l5p77_q2ex                            ;
}
float y_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1736706E-01;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.15126196E-01, 0.20832644E+00,-0.64423284E-02, 0.28345061E-01,
        -0.27041193E-01, 0.38644888E-01, 0.16990565E-01,-0.15502153E-01,
        -0.11408136E-01, 0.54653101E-02,-0.26655482E-01,-0.99575175E-02,
        -0.12599463E-01,-0.15030716E-02, 0.44823596E-02, 0.20266508E-02,
        -0.49470915E-02, 0.47276976E-06, 0.17345562E-02, 0.16920801E-01,
         0.12513955E-01,-0.21072689E-02, 0.10939785E-02, 0.43393424E-03,
         0.34230639E-03,-0.22932400E-02, 0.44027815E-03,-0.16004379E-02,
         0.18966967E-02, 0.13106553E-02,-0.15380180E-02, 0.72679068E-02,
        -0.58134119E-02, 0.53841439E-02,-0.50289342E-02, 0.23547665E-02,
         0.61598769E-03,-0.50821115E-03,-0.54609781E-03, 0.19243517E-02,
         0.15733061E-03, 0.10695230E-02,-0.18563939E-02,-0.90160093E-03,
        -0.49029820E-03,-0.29101805E-03,-0.35259957E-03, 0.27307525E-03,
         0.30664194E-03, 0.74395898E-03, 0.90721273E-03, 0.51018019E-03,
         0.57301906E-04, 0.77217887E-03, 0.11172565E-03,-0.14990107E-03,
        -0.15956006E-03,-0.86102472E-03,-0.63963863E-03, 0.15099294E-03,
         0.24144682E-02,-0.23935165E-03, 0.30405307E-02, 0.12796574E-02,
         0.10042959E-02, 0.47313044E-03, 0.14496668E-02,-0.60496015E-04,
        -0.52947802E-04, 0.20869798E-04,-0.16080536E-03,-0.26832633E-02,
        -0.28247382E-02,-0.92628085E-04,-0.10021243E-02, 0.94708841E-03,
        -0.37186390E-02,-0.29666156E-02,-0.24231621E-02,-0.13922520E-03,
         0.65077496E-04,-0.20078991E-02,-0.15974883E-03, 0.34031901E-03,
        -0.33814111E-03, 0.12490350E-03, 0.11715180E-04,-0.21900583E-03,
        -0.13066421E-03, 0.99418692E-04, 0.72818861E-04,-0.23875944E-03,
        -0.24768831E-04, 0.66614310E-04,-0.93106166E-04, 0.22451920E-03,
        -0.70438211E-04,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]*x11*x21
        +coeff[  7]            *x42
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[  8]        *x31*x41
        +coeff[  9]            *x41*x51
        +coeff[ 10]    *x22    *x41
        +coeff[ 11]    *x22*x31
        +coeff[ 12]*x11*x21    *x41
        +coeff[ 13]*x11
        +coeff[ 14]    *x21    *x41
        +coeff[ 15]*x12
        +coeff[ 16]*x11*x21*x31
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 17]        *x34
        +coeff[ 18]    *x21*x31
        +coeff[ 19]    *x22    *x42
        +coeff[ 20]    *x22*x31*x41
        +coeff[ 21]        *x32
        +coeff[ 22]*x11        *x41
        +coeff[ 23]*x11    *x31
        +coeff[ 24]        *x31    *x51
        +coeff[ 25]    *x21    *x42
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 26]                *x52
        +coeff[ 27]    *x21*x31*x41
        +coeff[ 28]    *x23
        +coeff[ 29]*x11*x22
        +coeff[ 30]*x12        *x41
        +coeff[ 31]*x11*x21    *x42
        +coeff[ 32]    *x24
        +coeff[ 33]*x11*x21*x31*x41
        +coeff[ 34]*x11*x23
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 35]        *x31*x42
        +coeff[ 36]            *x42*x51
        +coeff[ 37]            *x41*x52
        +coeff[ 38]*x12    *x31
        +coeff[ 39]    *x22*x32
        +coeff[ 40]            *x45
        +coeff[ 41]        *x32*x43
        +coeff[ 42]*x12*x22
        +coeff[ 43]        *x32*x45
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 44]*x11        *x42
        +coeff[ 45]    *x21*x32
        +coeff[ 46]*x11    *x31*x41
        +coeff[ 47]        *x31*x41*x51
        +coeff[ 48]*x12*x21
        +coeff[ 49]*x12        *x42
        +coeff[ 50]*x11*x21*x32
        +coeff[ 51]*x12    *x31*x41
        +coeff[ 52]    *x21        *x51
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 53]        *x32*x41
        +coeff[ 54]        *x33
        +coeff[ 55]    *x22        *x51
        +coeff[ 56]*x11*x21        *x51
        +coeff[ 57]    *x23    *x41
        +coeff[ 58]*x11*x22    *x41
        +coeff[ 59]        *x31*x44
        +coeff[ 60]    *x24    *x41
        +coeff[ 61]*x13*x21
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 62]*x11*x23    *x41
        +coeff[ 63]*x11*x23*x31
        +coeff[ 64]*x12*x22    *x41
        +coeff[ 65]*x12*x24*x31
        +coeff[ 66]            *x43
        +coeff[ 67]*x11    *x32
        +coeff[ 68]        *x31    *x52
        +coeff[ 69]                *x53
        +coeff[ 70]*x11*x22*x31
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 71]    *x22    *x43
        +coeff[ 72]    *x22*x31*x42
        +coeff[ 73]    *x22        *x52
        +coeff[ 74]    *x22*x32*x41
        +coeff[ 75]    *x24*x31
        +coeff[ 76]    *x24    *x42
        +coeff[ 77]    *x24*x31*x41
        +coeff[ 78]*x11*x23    *x42
        +coeff[ 79]    *x23*x33
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 80]        *x32    *x53
        +coeff[ 81]*x11*x23*x31*x41
        +coeff[ 82]*x11*x22*x33
        +coeff[ 83]*x12*x22*x32
        +coeff[ 84]*x11*x21*x33*x42
        +coeff[ 85]*x12*x21*x31*x44
        +coeff[ 86]*x11            *x51
        +coeff[ 87]            *x44
        +coeff[ 88]        *x31*x43
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 89]    *x21*x31*x42
        +coeff[ 90]            *x43*x51
        +coeff[ 91]    *x23*x31
        +coeff[ 92]*x12            *x51
        +coeff[ 93]    *x22*x31    *x51
        +coeff[ 94]*x12*x21    *x41
        +coeff[ 95]    *x23    *x42
        +coeff[ 96]*x12*x21*x31
        ;

    return v_y_l5p77_q2ex                            ;
}
float p_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 53;
    float avdat= -0.1712184E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 54]={
         0.22002454E-02,-0.65767500E-02,-0.30710399E-01, 0.24132314E-02,
        -0.35889824E-02, 0.14257259E-02, 0.49618641E-02,-0.16892413E-02,
        -0.99062592E-04, 0.62987965E-03, 0.10304960E-02, 0.78351970E-03,
        -0.71579841E-03, 0.19094155E-02, 0.10208240E-02, 0.12315810E-02,
         0.16914439E-03, 0.31068039E-03, 0.87979488E-03, 0.15428435E-03,
        -0.22656280E-03,-0.23111630E-03,-0.20263444E-03, 0.46150171E-03,
        -0.35126085E-03,-0.18487768E-02, 0.40116708E-03,-0.14202378E-03,
        -0.18129284E-03,-0.61425795E-04, 0.45872488E-03,-0.71155839E-03,
        -0.13475586E-03,-0.46059137E-03, 0.47601949E-03, 0.98932280E-04,
         0.77500830E-04,-0.61513722E-03,-0.10457665E-03, 0.38566656E-03,
        -0.11053132E-03, 0.22877216E-03,-0.40922649E-03,-0.83586609E-04,
        -0.10351271E-03,-0.70320784E-04,-0.96399098E-03, 0.78537232E-04,
        -0.14244817E-03,-0.34602600E-03,-0.18956109E-03, 0.12583853E-03,
        -0.88525150E-04,
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
    float x24 = x23*x2;
    float x25 = x24*x2;
    float x26 = x25*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]            *x42
        +coeff[  6]            *x41*x51
        +coeff[  7]*x11*x21
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[  8]    *x23
        +coeff[  9]    *x21
        +coeff[ 10]        *x31*x41
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]                *x52
        +coeff[ 13]    *x22    *x41
        +coeff[ 14]    *x22        *x51
        +coeff[ 15]*x11*x21    *x41
        +coeff[ 16]*x11
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 17]    *x21    *x41
        +coeff[ 18]    *x22*x31
        +coeff[ 19]        *x32
        +coeff[ 20]    *x21        *x51
        +coeff[ 21]*x11        *x41
        +coeff[ 22]*x12
        +coeff[ 23]*x11*x21*x31
        +coeff[ 24]            *x41*x52
        +coeff[ 25]    *x22    *x42
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 26]*x11*x21        *x51
        +coeff[ 27]    *x22    *x43
        +coeff[ 28]    *x24*x31*x41
        +coeff[ 29]*x11    *x31
        +coeff[ 30]    *x21    *x42
        +coeff[ 31]            *x43
        +coeff[ 32]        *x31*x41*x51
        +coeff[ 33]    *x22    *x41*x51
        +coeff[ 34]*x11*x23
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 35]*x12        *x41
        +coeff[ 36]    *x23*x31*x41
        +coeff[ 37]*x11*x21    *x42
        +coeff[ 38]    *x22*x31*x42
        +coeff[ 39]    *x26
        +coeff[ 40]*x11*x23*x31*x41
        +coeff[ 41]    *x21*x31*x41
        +coeff[ 42]        *x31*x42
        +coeff[ 43]            *x42*x51
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 44]*x11*x22
        +coeff[ 45]        *x31    *x52
        +coeff[ 46]    *x22*x31*x41
        +coeff[ 47]                *x53
        +coeff[ 48]    *x22*x31    *x51
        +coeff[ 49]*x11*x21*x31*x41
        +coeff[ 50]*x11*x21    *x41*x51
        +coeff[ 51]*x12*x22
        +coeff[ 52]        *x34*x41
        ;

    return v_p_l5p77_q2ex                            ;
}
float l_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 90;
    float avdat= -0.2904264E-02;
    float xmin[10]={
        -0.14996E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.46084E-01, 0.14970E-01, 0.25393E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.23306834E-02, 0.18473481E-02, 0.55993362E-02,-0.42647952E-02,
         0.60337433E-03,-0.19430276E-02,-0.76624537E-02, 0.22427004E-03,
         0.18082279E-02, 0.10496762E-02,-0.12120018E-02,-0.49900133E-02,
        -0.12796808E-02, 0.22871459E-02,-0.13412113E-03,-0.11941694E-03,
        -0.13018276E-03, 0.96631714E-03, 0.40376847E-03,-0.77435031E-03,
         0.13386451E-02,-0.51709299E-03, 0.77697769E-05,-0.96901073E-04,
         0.19165689E-04, 0.11823583E-03,-0.14568059E-03,-0.49714989E-04,
         0.66630084E-04, 0.22397947E-03, 0.70643478E-05,-0.28050540E-03,
         0.77536615E-03, 0.11465104E-03,-0.18964037E-03,-0.20253929E-03,
         0.15180127E-03, 0.68995880E-03, 0.13662823E-02,-0.15606138E-04,
        -0.13844589E-03, 0.18274314E-03, 0.69555637E-04, 0.85494539E-04,
        -0.50112172E-04, 0.33375443E-04, 0.15644683E-03,-0.32859910E-03,
         0.33011916E-03, 0.40247085E-03,-0.94007183E-03, 0.78733696E-03,
        -0.92349372E-04,-0.29978839E-04, 0.19237988E-04,-0.35778710E-04,
         0.24424586E-04,-0.87405693E-04,-0.55579796E-04,-0.79645077E-04,
        -0.13604912E-03,-0.89510919E-04,-0.43247519E-05,-0.67925727E-03,
         0.21757541E-03,-0.27051702E-03,-0.12933179E-03, 0.89527803E-05,
         0.10653551E-03,-0.85266169E-04,-0.75386022E-04,-0.31770243E-04,
        -0.30521282E-04,-0.29207789E-04,-0.34920027E-04,-0.28636405E-04,
         0.42862983E-04,-0.13267179E-03,-0.13231569E-03, 0.14066981E-03,
         0.26196578E-04,-0.18861210E-03,-0.63364052E-04, 0.49613733E-04,
        -0.25386305E-04,-0.26279851E-04, 0.80769851E-05,-0.71562913E-05,
         0.54620646E-05, 0.64387270E-04,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;
    float x54 = x53*x5;

//                 function

    float v_l_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]    *x22
        +coeff[  4]    *x21    *x41
        +coeff[  5]        *x31*x41
        +coeff[  6]            *x42
        +coeff[  7]        *x31    *x51
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[  8]            *x41*x51
        +coeff[  9]*x11*x21
        +coeff[ 10]    *x22*x31
        +coeff[ 11]    *x22    *x41
        +coeff[ 12]*x11*x21    *x41
        +coeff[ 13]    *x22    *x42
        +coeff[ 14]    *x21
        +coeff[ 15]                *x52
        +coeff[ 16]*x12
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 17]            *x43
        +coeff[ 18]            *x42*x51
        +coeff[ 19]    *x24
        +coeff[ 20]    *x22*x31*x41
        +coeff[ 21]*x11*x23
        +coeff[ 22]        *x33*x42
        +coeff[ 23]                *x51
        +coeff[ 24]*x11
        +coeff[ 25]    *x21*x31
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 26]        *x32
        +coeff[ 27]    *x21        *x51
        +coeff[ 28]*x11        *x41
        +coeff[ 29]    *x23
        +coeff[ 30]        *x33
        +coeff[ 31]    *x21    *x42
        +coeff[ 32]        *x31*x42
        +coeff[ 33]        *x31*x41*x51
        +coeff[ 34]            *x41*x52
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 35]*x11*x21*x31
        +coeff[ 36]*x11*x21        *x51
        +coeff[ 37]*x11*x21    *x42
        +coeff[ 38]    *x24    *x41
        +coeff[ 39]*x11            *x51
        +coeff[ 40]    *x21*x31*x41
        +coeff[ 41]        *x32*x41
        +coeff[ 42]    *x22        *x51
        +coeff[ 43]*x11*x22
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 44]*x12        *x41
        +coeff[ 45]*x12            *x51
        +coeff[ 46]    *x22*x32
        +coeff[ 47]    *x23    *x41
        +coeff[ 48]*x11*x21*x31*x41
        +coeff[ 49]    *x24*x31
        +coeff[ 50]    *x22    *x43
        +coeff[ 51]*x11*x23    *x41
        +coeff[ 52]    *x24*x31*x42
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 53]        *x31    *x52
        +coeff[ 54]                *x53
        +coeff[ 55]*x11        *x42
        +coeff[ 56]*x12    *x31
        +coeff[ 57]    *x23*x31
        +coeff[ 58]    *x22*x31    *x51
        +coeff[ 59]            *x43*x51
        +coeff[ 60]*x11*x22    *x41
        +coeff[ 61]*x12*x22
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 62]*x13    *x31
        +coeff[ 63]    *x22*x31*x42
        +coeff[ 64]*x11*x23*x31
        +coeff[ 65]*x11*x21    *x43
        +coeff[ 66]*x11*x23*x31*x42
        +coeff[ 67]    *x21*x31    *x51
        +coeff[ 68]    *x21    *x43
        +coeff[ 69]            *x44
        +coeff[ 70]        *x31*x42*x51
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 71]    *x22        *x52
        +coeff[ 72]            *x42*x52
        +coeff[ 73]*x11*x22*x31
        +coeff[ 74]*x11*x21*x31    *x51
        +coeff[ 75]*x11*x21    *x41*x51
        +coeff[ 76]*x12        *x42
        +coeff[ 77]    *x22*x32*x41
        +coeff[ 78]*x11*x21*x31*x42
        +coeff[ 79]*x12*x22    *x41
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 80]    *x23*x31*x42
        +coeff[ 81]    *x22*x31*x43
        +coeff[ 82]    *x22*x32*x41*x51
        +coeff[ 83]*x11*x23*x32
        +coeff[ 84]*x11*x21        *x54
        +coeff[ 85]    *x23*x34
        +coeff[ 86]    *x21        *x52
        +coeff[ 87]*x11    *x31*x41
        +coeff[ 88]*x11        *x41*x51
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 89]    *x21*x31*x42
        ;

    return v_l_l5p77_q2ex                            ;
}
float x_l5p77_q3en                            (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.6418221E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.45081E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
        -0.24700824E-02,-0.35758495E-01, 0.13826573E+00, 0.16514866E+00,
        -0.57779360E-02, 0.26701506E-01, 0.68964630E-01,-0.66488348E-04,
        -0.12644759E-02, 0.80956648E-04, 0.31146253E-03,-0.65638572E-02,
        -0.60261148E-02, 0.59409463E-02, 0.15962906E-01, 0.73446832E-02,
         0.25031475E-01, 0.17760701E-01,-0.31773336E-01, 0.11618328E-02,
        -0.23717629E-02,-0.24837812E-02, 0.11755564E-01, 0.11612300E-01,
        -0.19822408E-01,-0.26349774E-01, 0.31223930E-02, 0.50109951E-02,
         0.43313960E-02,-0.33811429E-02,-0.16824959E-01,-0.82319519E-02,
         0.96271106E-03, 0.26993568E-02,-0.45582359E-02, 0.14980320E-02,
         0.12109120E-02,-0.35249765E-02,-0.56076925E-02, 0.20313571E-03,
        -0.23591719E-03, 0.55005500E-03,-0.66259108E-03,-0.32739197E-02,
        -0.60441205E-03, 0.15786620E-02, 0.74872375E-02, 0.64492030E-02,
         0.17416790E-03, 0.13562129E-03, 0.37563793E-03, 0.39626390E-03,
        -0.39985409E-03, 0.26272528E-03, 0.10937292E-02, 0.63835050E-03,
         0.10012973E-02, 0.26186940E-03,-0.53187821E-03, 0.26132320E-02,
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
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_x_l5p77_q3en                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]    *x21
        +coeff[  4]*x11*x21
        +coeff[  5]    *x21        *x51
        +coeff[  6]    *x21    *x41
        +coeff[  7]*x13        *x41
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[  8]*x12*x21*x31
        +coeff[  9]    *x21*x31    *x52
        +coeff[ 10]    *x21*x33
        +coeff[ 11]            *x41
        +coeff[ 12]                *x52
        +coeff[ 13]*x11    *x31
        +coeff[ 14]*x11        *x41
        +coeff[ 15]    *x22
        +coeff[ 16]    *x21*x31
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 17]    *x23
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]    *x23*x31*x41
        +coeff[ 20]        *x31
        +coeff[ 21]            *x42
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x22    *x41
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]    *x23    *x41
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 26]*x12*x21
        +coeff[ 27]    *x21    *x41*x51
        +coeff[ 28]    *x22*x31
        +coeff[ 29]    *x21*x32
        +coeff[ 30]*x11*x22    *x41
        +coeff[ 31]    *x23*x31
        +coeff[ 32]            *x41*x51
        +coeff[ 33]*x11*x21    *x41
        +coeff[ 34]*x11        *x42
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 35]    *x22        *x51
        +coeff[ 36]    *x21*x31    *x51
        +coeff[ 37]*x12*x21    *x41
        +coeff[ 38]*x11*x22*x31
        +coeff[ 39]*x13*x21*x31
        +coeff[ 40]*x13    *x31*x41
        +coeff[ 41]        *x31*x41
        +coeff[ 42]*x11    *x32
        +coeff[ 43]*x11    *x31*x41
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 44]            *x42*x51
        +coeff[ 45]    *x23        *x51
        +coeff[ 46]    *x21*x31*x42
        +coeff[ 47]    *x21    *x43
        +coeff[ 48]*x12*x21*x32*x41
        +coeff[ 49]        *x32
        +coeff[ 50]                *x53
        +coeff[ 51]*x11        *x41*x51
        +coeff[ 52]    *x21        *x52
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 53]            *x41*x52
        +coeff[ 54]*x11*x21*x31
        +coeff[ 55]*x11*x22        *x51
        +coeff[ 56]*x11*x23
        +coeff[ 57]*x13            *x52
        +coeff[ 58]    *x22*x31*x41
        +coeff[ 59]    *x21*x32*x41
        ;

    return v_x_l5p77_q3en                            ;
}
float t_l5p77_q3en                            (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.2362935E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.45081E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.12230999E-02, 0.91067282E-02,-0.61261922E-01, 0.21553643E-02,
         0.22829112E-01, 0.19521090E-02,-0.21964027E-01,-0.39225258E-02,
         0.45515152E-03, 0.26138106E-02, 0.11176683E-03,-0.15484167E-02,
        -0.18672834E-02,-0.81293881E-02,-0.50908467E-02,-0.14217411E-02,
        -0.55307648E-02, 0.80991220E-02,-0.11815460E-03, 0.25654439E-03,
         0.74827322E-03,-0.59690879E-03,-0.10023933E-02,-0.37603187E-02,
        -0.37509636E-02, 0.53211395E-02, 0.76916264E-02,-0.14732672E-02,
         0.86694380E-03, 0.50945478E-02,-0.12697290E-03,-0.58456510E-03,
        -0.57406822E-03,-0.83034486E-03, 0.11053775E-02, 0.16442481E-02,
         0.39362800E-03, 0.17163526E-02, 0.11865050E-02,-0.85160637E-03,
         0.72986313E-05, 0.90122521E-04,-0.34260118E-03, 0.21512150E-03,
        -0.45768636E-04, 0.34670165E-03,-0.21497806E-03, 0.16055314E-03,
        -0.12080681E-03, 0.56673115E-03,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_l5p77_q3en                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]*x11*x21
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x21        *x51
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[  8]*x12*x21*x31
        +coeff[  9]    *x23*x31
        +coeff[ 10]*x13        *x41
        +coeff[ 11]    *x22
        +coeff[ 12]*x11    *x31
        +coeff[ 13]    *x21*x31
        +coeff[ 14]*x11        *x41
        +coeff[ 15]                *x52
        +coeff[ 16]    *x23
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]*x13*x22
        +coeff[ 19]    *x23*x31*x41
        +coeff[ 20]        *x31
        +coeff[ 21]*x11            *x51
        +coeff[ 22]*x12*x21
        +coeff[ 23]*x11*x22
        +coeff[ 24]    *x22    *x41
        +coeff[ 25]    *x21*x31*x41
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 26]    *x23    *x41
        +coeff[ 27]    *x22*x31
        +coeff[ 28]    *x21*x32
        +coeff[ 29]*x11*x22    *x41
        +coeff[ 30]*x13        *x42
        +coeff[ 31]        *x31*x41
        +coeff[ 32]            *x42
        +coeff[ 33]*x11*x21    *x41
        +coeff[ 34]*x11    *x31*x41
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 35]*x11        *x42
        +coeff[ 36]    *x21        *x52
        +coeff[ 37]*x11*x22*x31
        +coeff[ 38]*x12*x21    *x41
        +coeff[ 39]    *x21    *x43
        +coeff[ 40]*x11*x23*x31
        +coeff[ 41]        *x31    *x51
        +coeff[ 42]*x11*x21*x31
        +coeff[ 43]*x11    *x32
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 44]        *x31*x41*x51
        +coeff[ 45]            *x42*x51
        +coeff[ 46]            *x41*x52
        +coeff[ 47]                *x53
        +coeff[ 48]*x12    *x32
        +coeff[ 49]    *x22    *x42
        ;

    return v_t_l5p77_q3en                            ;
}
float y_l5p77_q3en                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1749039E-01;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.45081E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.63453768E-02, 0.76745503E-01,-0.90235546E-02,-0.20516876E-01,
        -0.30301195E-01,-0.13529117E-01,-0.53768326E-01,-0.10207305E-01,
         0.42131603E-01,-0.45224125E-02, 0.10308760E-01, 0.68421677E-01,
         0.17372131E-01, 0.65107918E-02, 0.98821567E-02,-0.16647203E-01,
        -0.70749191E-02,-0.46321023E-01,-0.56016101E-02,-0.76031177E-02,
        -0.11873765E-01, 0.12987575E-01,-0.14396918E-02,-0.19888324E-02,
        -0.67635351E-02, 0.18721714E-02,-0.18570682E-02,-0.69626211E-02,
        -0.34794044E-02,-0.57142093E-02, 0.35461232E-02, 0.13641753E-02,
        -0.58347783E-02,-0.19176210E-02,-0.35556993E-02,-0.10720186E-01,
        -0.39769234E-02,-0.48962031E-02,-0.59834239E-02,-0.17566563E-02,
        -0.25684625E-03,-0.25441356E-02, 0.83784899E-02,-0.17122980E-02,
        -0.34017947E-02, 0.57119527E-02, 0.10252097E-02,-0.88453671E-03,
         0.56290429E-03, 0.10513006E-02,-0.13570490E-02,-0.24736621E-02,
         0.91692165E-03,-0.11745356E-02, 0.71254192E-03,-0.35700656E-03,
         0.62502874E-03,-0.54028042E-03, 0.32336749E-02, 0.37079796E-02,
         0.33982960E-02, 0.15060988E-02, 0.16109351E-02, 0.33709500E-02,
        -0.17145651E-02, 0.78183878E-02,-0.13762846E-02,-0.10564011E-02,
        -0.11660439E-02,-0.30594340E-02, 0.15304087E-02, 0.37391356E-02,
         0.41121636E-02, 0.45235839E-03, 0.94771414E-03, 0.57595753E-03,
         0.22699629E-03,-0.28970937E-03, 0.28521343E-03, 0.10524103E-02,
         0.10316517E-02, 0.10663021E-02, 0.72210608E-03,-0.27645475E-02,
         0.63886278E-03,-0.16009271E-02, 0.78524044E-03, 0.68777846E-03,
        -0.24048469E-03,-0.19252710E-02, 0.38434006E-03, 0.59907930E-03,
        -0.28657666E-03, 0.36991623E-03, 0.14206953E-02,-0.57226117E-03,
        -0.54822973E-03,
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
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_q3en                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]            *x42
        +coeff[  6]    *x21    *x41
        +coeff[  7]        *x31*x41
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[  8]    *x22
        +coeff[  9]    *x21*x31
        +coeff[ 10]*x11        *x41
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]        *x31    *x51
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]                *x52
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]*x11        *x42
        +coeff[ 19]    *x23
        +coeff[ 20]    *x22*x31
        +coeff[ 21]    *x22        *x51
        +coeff[ 22]*x11
        +coeff[ 23]        *x32
        +coeff[ 24]            *x43
        +coeff[ 25]*x12
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 26]*x11            *x51
        +coeff[ 27]    *x21*x31*x41
        +coeff[ 28]    *x21    *x41*x51
        +coeff[ 29]            *x41*x52
        +coeff[ 30]*x11*x21        *x51
        +coeff[ 31]*x11    *x31
        +coeff[ 32]*x11*x21    *x41
        +coeff[ 33]*x11*x22
        +coeff[ 34]*x11*x21*x31
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 35]    *x24
        +coeff[ 36]*x11*x22    *x41
        +coeff[ 37]    *x22    *x41*x51
        +coeff[ 38]*x11*x23
        +coeff[ 39]        *x31*x42
        +coeff[ 40]            *x42*x51
        +coeff[ 41]*x11    *x31*x41
        +coeff[ 42]    *x21    *x43
        +coeff[ 43]*x12        *x41
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 44]    *x22    *x42
        +coeff[ 45]    *x21*x31*x42
        +coeff[ 46]    *x21        *x52
        +coeff[ 47]        *x31    *x52
        +coeff[ 48]*x12            *x51
        +coeff[ 49]                *x53
        +coeff[ 50]*x12*x21    *x41
        +coeff[ 51]*x11*x21    *x41*x51
        +coeff[ 52]*x13        *x41*x51
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 53]        *x31*x41*x51
        +coeff[ 54]    *x21*x31    *x51
        +coeff[ 55]        *x32    *x51
        +coeff[ 56]*x12*x21
        +coeff[ 57]*x12    *x31
        +coeff[ 58]    *x23    *x41
        +coeff[ 59]    *x22*x31*x41
        +coeff[ 60]*x11*x21    *x42
        +coeff[ 61]    *x23*x31
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 62]    *x22*x32
        +coeff[ 63]*x11*x21*x31*x41
        +coeff[ 64]    *x22*x31    *x51
        +coeff[ 65]    *x23    *x42
        +coeff[ 66]*x12*x22
        +coeff[ 67]*x11*x21*x31    *x51
        +coeff[ 68]    *x22        *x52
        +coeff[ 69]    *x24    *x41
        +coeff[ 70]    *x23*x31*x41
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 71]*x11*x22    *x42
        +coeff[ 72]    *x22    *x42*x51
        +coeff[ 73]    *x21    *x44*x51
        +coeff[ 74]    *x21*x31*x43*x51
        +coeff[ 75]            *x43*x53
        +coeff[ 76]        *x32*x41
        +coeff[ 77]*x11    *x32
        +coeff[ 78]*x11    *x31    *x51
        +coeff[ 79]*x11        *x43
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 80]            *x43*x51
        +coeff[ 81]    *x21*x32*x41
        +coeff[ 82]*x11    *x31*x42
        +coeff[ 83]    *x21    *x42*x51
        +coeff[ 84]        *x31*x42*x51
        +coeff[ 85]    *x21*x31*x41*x51
        +coeff[ 86]*x12        *x42
        +coeff[ 87]*x11*x21*x32
        +coeff[ 88]    *x23        *x51
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 89]    *x22    *x43
        +coeff[ 90]*x12    *x31*x41
        +coeff[ 91]    *x21    *x41*x52
        +coeff[ 92]*x12*x21*x31
        +coeff[ 93]*x11*x22        *x51
        +coeff[ 94]*x11*x21    *x43
        +coeff[ 95]    *x23*x32
        +coeff[ 96]*x11*x21        *x52
        ;

    return v_y_l5p77_q3en                            ;
}
float p_l5p77_q3en                            (float *x,int m){
    int ncoeff= 90;
    float avdat=  0.2092765E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.45081E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
        -0.57652127E-04,-0.11287207E-02,-0.59440755E-02, 0.10323259E-02,
        -0.40268791E-02, 0.55196160E-02,-0.12180732E-02,-0.10881913E-01,
        -0.12163619E-02,-0.15501649E-02, 0.99416194E-03, 0.17017678E-02,
         0.13772140E-01, 0.23088390E-02,-0.17342857E-02,-0.16844011E-02,
         0.18696013E-02,-0.83976267E-02,-0.31997368E-02, 0.31601847E-02,
        -0.14460272E-02,-0.97234308E-03,-0.44951477E-03,-0.46570178E-04,
        -0.11208013E-02, 0.78238867E-03, 0.22557403E-03,-0.20617964E-03,
        -0.23874053E-03, 0.33874283E-03,-0.14773656E-02,-0.30494275E-03,
         0.25294544E-03,-0.37649256E-03,-0.55702287E-03,-0.21247317E-02,
        -0.14491341E-02,-0.84984247E-04,-0.15790163E-02,-0.41280067E-03,
        -0.36101980E-03,-0.18338655E-02, 0.25710493E-03, 0.70850871E-03,
        -0.58716588E-03, 0.20471173E-03,-0.48897602E-03,-0.87629468E-03,
        -0.27896170E-03,-0.69904956E-03,-0.57030102E-03, 0.13935160E-03,
         0.14370799E-02,-0.96385146E-03, 0.24142163E-03, 0.15925942E-02,
         0.13401617E-03,-0.21128615E-03,-0.17943207E-03, 0.44529632E-03,
         0.10405375E-02, 0.42840015E-03,-0.12523070E-03,-0.89831163E-04,
        -0.18388417E-03,-0.47635778E-04, 0.16057907E-03,-0.12401365E-03,
         0.88975689E-03, 0.74770192E-04,-0.25536693E-03, 0.34414006E-04,
        -0.65820677E-04,-0.89465793E-04, 0.55619533E-03, 0.12537453E-03,
         0.30061326E-03, 0.38825991E-03, 0.49031113E-03,-0.12932959E-03,
         0.10045538E-03, 0.52232289E-03,-0.25812356E-03, 0.39130781E-03,
         0.25611889E-03,-0.29738725E-03, 0.17153154E-03, 0.26389590E-03,
         0.20067014E-03,-0.47911750E-03,
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
    float x24 = x23*x2;
    float x25 = x24*x2;
    float x26 = x25*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_q3en                            =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21*x31
        +coeff[  7]    *x21    *x41
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[  8]        *x31*x41
        +coeff[  9]            *x42
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]            *x41*x51
        +coeff[ 13]*x11*x21
        +coeff[ 14]                *x52
        +coeff[ 15]    *x22*x31
        +coeff[ 16]*x11        *x41
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]    *x22        *x51
        +coeff[ 20]            *x43
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]*x11*x22
        +coeff[ 23]*x11    *x32
        +coeff[ 24]*x11        *x42
        +coeff[ 25]*x11*x21        *x51
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 26]    *x25
        +coeff[ 27]*x11
        +coeff[ 28]        *x32
        +coeff[ 29]*x11    *x31
        +coeff[ 30]    *x21*x31*x41
        +coeff[ 31]*x11            *x51
        +coeff[ 32]*x12
        +coeff[ 33]*x11*x21*x31
        +coeff[ 34]            *x41*x52
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 35]    *x22    *x42
        +coeff[ 36]    *x22    *x41*x51
        +coeff[ 37]    *x26
        +coeff[ 38]    *x23
        +coeff[ 39]        *x31*x42
        +coeff[ 40]        *x31*x41*x51
        +coeff[ 41]    *x24
        +coeff[ 42]    *x21        *x52
        +coeff[ 43]    *x23    *x41
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 44]*x11    *x31*x41
        +coeff[ 45]                *x53
        +coeff[ 46]    *x22*x31    *x51
        +coeff[ 47]*x11*x23
        +coeff[ 48]*x12        *x41
        +coeff[ 49]*x11*x22    *x41
        +coeff[ 50]*x11*x21    *x41*x51
        +coeff[ 51]    *x21*x31    *x51
        +coeff[ 52]    *x21    *x43
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 53]    *x24    *x41
        +coeff[ 54]*x11*x21*x31*x41
        +coeff[ 55]    *x23    *x42
        +coeff[ 56]*x12            *x51
        +coeff[ 57]*x11*x21*x31    *x51
        +coeff[ 58]*x12*x21    *x41
        +coeff[ 59]    *x22    *x41*x52
        +coeff[ 60]*x11*x22    *x42
        +coeff[ 61]    *x23*x31*x42
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 62]        *x32    *x53
        +coeff[ 63]    *x21*x32
        +coeff[ 64]            *x42*x51
        +coeff[ 65]        *x31    *x52
        +coeff[ 66]    *x22*x32
        +coeff[ 67]    *x23        *x51
        +coeff[ 68]    *x21*x31*x42
        +coeff[ 69]*x11        *x41*x51
        +coeff[ 70]    *x21    *x42*x51
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 71]*x12*x21
        +coeff[ 72]*x12    *x31
        +coeff[ 73]*x11*x22*x31
        +coeff[ 74]    *x23*x31*x41
        +coeff[ 75]*x11*x21    *x42
        +coeff[ 76]*x11        *x43
        +coeff[ 77]    *x23    *x41*x51
        +coeff[ 78]    *x22    *x42*x51
        +coeff[ 79]*x12*x22
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 80]    *x25*x31
        +coeff[ 81]*x11*x22*x31*x41
        +coeff[ 82]    *x24*x31*x41
        +coeff[ 83]    *x23*x32*x41
        +coeff[ 84]*x11*x21    *x43
        +coeff[ 85]    *x23*x31*x41*x51
        +coeff[ 86]*x12*x21    *x42
        +coeff[ 87]*x11*x22*x31*x42
        +coeff[ 88]*x11    *x31*x44
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 89]*x11*x25    *x41
        ;

    return v_p_l5p77_q3en                            ;
}
float l_l5p77_q3en                            (float *x,int m){
    int ncoeff= 71;
    float avdat= -0.4175426E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32620E-01,-0.49941E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.45081E-01, 0.14970E-01, 0.25342E-01, 0.49990E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 72]={
         0.74865641E-02,-0.24521141E+00,-0.31385407E-01, 0.41853048E-01,
        -0.10910466E-01,-0.90028211E-01,-0.15446720E-01,-0.21047676E-01,
         0.15051760E-02, 0.10567617E-01,-0.53240912E-03, 0.51016905E-02,
         0.14568469E-01,-0.32866105E-01, 0.92079658E-02,-0.79061612E-02,
        -0.22437869E-01,-0.23726957E-01, 0.38276702E-01,-0.16531397E-01,
         0.16865115E-03,-0.14618038E-02,-0.10837511E-01,-0.29789147E-02,
        -0.70916880E-02, 0.25087101E-01,-0.41224249E-02, 0.33357549E-01,
         0.21918034E-01,-0.45979284E-02, 0.11156984E-02, 0.44019865E-02,
        -0.49679675E-02, 0.62054428E-02, 0.28046910E-03,-0.19016974E-03,
        -0.16191558E-02, 0.44574835E-02, 0.73915031E-02, 0.45454935E-02,
         0.48850514E-02,-0.57023182E-03,-0.30575381E-03,-0.10175398E-02,
         0.13549023E-02, 0.84000180E-03, 0.93917298E-03,-0.39551808E-02,
         0.52709561E-02,-0.11262064E-01,-0.99290516E-02,-0.26410695E-02,
         0.16889244E-02, 0.71216928E-03,-0.85048063E-03, 0.38589264E-04,
        -0.17686297E-03, 0.15228303E-02, 0.70613774E-03, 0.53324044E-03,
        -0.31681516E-03, 0.47260287E-03,-0.31860141E-03, 0.35181190E-02,
        -0.40199957E-02, 0.67023397E-03, 0.26754415E-03, 0.33945332E-02,
         0.53527433E-03,-0.35981452E-02, 0.12837733E-02,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_l_l5p77_q3en                            =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]                *x51
        +coeff[  3]*x11
        +coeff[  4]    *x22
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x21        *x51
        +coeff[  7]*x11        *x41
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[  8]            *x43
        +coeff[  9]    *x23*x31
        +coeff[ 10]    *x21*x33
        +coeff[ 11]        *x31
        +coeff[ 12]            *x41
        +coeff[ 13]    *x21*x31
        +coeff[ 14]*x11*x21
        +coeff[ 15]*x11    *x31
        +coeff[ 16]    *x23
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]*x11*x22
        +coeff[ 20]            *x44
        +coeff[ 21]    *x23*x31*x41
        +coeff[ 22]            *x42
        +coeff[ 23]*x11            *x51
        +coeff[ 24]    *x22*x31
        +coeff[ 25]    *x21*x31*x41
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 26]*x12*x21
        +coeff[ 27]    *x23    *x41
        +coeff[ 28]*x11*x22    *x41
        +coeff[ 29]        *x31*x41
        +coeff[ 30]            *x41*x51
        +coeff[ 31]    *x21*x32
        +coeff[ 32]*x11*x21    *x41
        +coeff[ 33]*x11        *x42
        +coeff[ 34]*x11    *x33*x41
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 35]*x11*x24*x31
        +coeff[ 36]*x11*x21*x31
        +coeff[ 37]*x11    *x31*x41
        +coeff[ 38]*x11*x22*x31
        +coeff[ 39]*x12*x21    *x41
        +coeff[ 40]    *x23    *x44
        +coeff[ 41]        *x32
        +coeff[ 42]                *x52
        +coeff[ 43]    *x22        *x51
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 44]            *x42*x51
        +coeff[ 45]    *x21        *x52
        +coeff[ 46]*x11    *x32
        +coeff[ 47]    *x24
        +coeff[ 48]    *x22    *x42
        +coeff[ 49]    *x21*x31*x42
        +coeff[ 50]    *x21    *x43
        +coeff[ 51]*x11*x23
        +coeff[ 52]*x12*x21*x31
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 53]        *x33*x41*x51
        +coeff[ 54]    *x24*x31*x41
        +coeff[ 55]    *x23*x32*x41
        +coeff[ 56]*x12
        +coeff[ 57]        *x31*x42
        +coeff[ 58]    *x21*x31    *x51
        +coeff[ 59]    *x21    *x41*x51
        +coeff[ 60]            *x41*x52
        +coeff[ 61]*x11        *x41*x51
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 62]*x13
        +coeff[ 63]    *x22*x31*x41
        +coeff[ 64]    *x21*x32*x41
        +coeff[ 65]    *x22    *x41*x51
        +coeff[ 66]*x12    *x31    *x51
        +coeff[ 67]    *x24    *x41
        +coeff[ 68]        *x34*x41
        +coeff[ 69]    *x23    *x42
        +coeff[ 70]*x11*x24
        ;

    return v_l_l5p77_q3en                            ;
}
float x_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 56;
    float avdat=  0.7319451E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 57]={
         0.10439183E-02,-0.27561534E-01, 0.32711279E+00, 0.52344508E-01,
        -0.22388972E-01, 0.28292978E-01, 0.36700785E-01, 0.72228219E-02,
         0.13365578E-01,-0.35148771E-02, 0.81691006E-02,-0.63538733E-02,
         0.91950642E-02,-0.22498796E-01,-0.12046620E-02, 0.30162190E-02,
        -0.14073421E-02, 0.55933627E-02, 0.76646856E-02, 0.60086949E-02,
        -0.10767043E-01,-0.17213091E-01,-0.15088307E-02,-0.29952340E-02,
         0.13123141E-02, 0.26512244E-02, 0.16173535E-02, 0.17497826E-02,
         0.26413652E-02,-0.92483005E-02,-0.46082973E-03, 0.16325900E-02,
         0.80655265E-03, 0.18045115E-02,-0.17634347E-02, 0.20321819E-02,
        -0.16248162E-02,-0.42555630E-02, 0.24127178E-02, 0.40181377E-03,
        -0.21785892E-03, 0.25004064E-03,-0.37196145E-03, 0.67090249E-03,
        -0.80455234E-03, 0.73374732E-03,-0.14913535E-02,-0.16000012E-02,
         0.51724666E-03,-0.64147222E-04,-0.12513391E-02, 0.16301681E-02,
        -0.28892606E-02, 0.35914148E-02, 0.27192484E-02, 0.24287154E-02,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_x_l5p77_q3ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]    *x21
        +coeff[  4]                *x52
        +coeff[  5]    *x21        *x51
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x22
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[  8]    *x21*x31
        +coeff[  9]            *x41
        +coeff[ 10]*x11        *x41
        +coeff[ 11]            *x42
        +coeff[ 12]    *x23
        +coeff[ 13]    *x21    *x42
        +coeff[ 14]        *x31
        +coeff[ 15]*x11    *x31
        +coeff[ 16]        *x31*x41
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 17]*x11*x22
        +coeff[ 18]    *x21    *x41*x51
        +coeff[ 19]    *x22    *x41
        +coeff[ 20]    *x21*x31*x41
        +coeff[ 21]    *x23    *x41
        +coeff[ 22]*x11            *x51
        +coeff[ 23]*x11*x21
        +coeff[ 24]            *x41*x51
        +coeff[ 25]                *x53
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 26]*x12*x21
        +coeff[ 27]    *x21*x31    *x51
        +coeff[ 28]*x11*x23
        +coeff[ 29]*x11*x22    *x41
        +coeff[ 30]    *x23*x32    *x52
        +coeff[ 31]*x11*x21    *x41
        +coeff[ 32]    *x22        *x51
        +coeff[ 33]    *x22*x31
        +coeff[ 34]    *x21*x32
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 35]    *x23        *x51
        +coeff[ 36]    *x21    *x42*x51
        +coeff[ 37]    *x23*x31
        +coeff[ 38]    *x22    *x42
        +coeff[ 39]*x13*x22*x31
        +coeff[ 40]*x12
        +coeff[ 41]        *x31    *x51
        +coeff[ 42]*x11*x21        *x51
        +coeff[ 43]*x11        *x41*x51
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 44]    *x21        *x52
        +coeff[ 45]*x11*x21*x31
        +coeff[ 46]*x11    *x31*x41
        +coeff[ 47]*x11        *x42
        +coeff[ 48]            *x43
        +coeff[ 49]*x13    *x31
        +coeff[ 50]*x12*x21    *x41
        +coeff[ 51]*x11*x22        *x51
        +coeff[ 52]*x11*x22*x31
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 53]    *x21*x31*x42
        +coeff[ 54]    *x21    *x43
        +coeff[ 55]    *x23*x32*x41
        ;

    return v_x_l5p77_q3ex                            ;
}
float t_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1062619E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.11076281E-02,-0.33732629E-02,-0.18544443E-01, 0.17689139E-04,
         0.53132735E-04, 0.11300905E+00, 0.53365822E-02,-0.10616443E-01,
         0.12934103E-02,-0.23486356E-02, 0.15443529E-03,-0.76080963E-03,
        -0.77146955E-03,-0.22164071E-02, 0.19790994E-02,-0.79686492E-03,
         0.11157951E-02,-0.42814895E-03, 0.36759514E-03,-0.72377596E-04,
        -0.18061096E-03,-0.26637354E-03, 0.14153871E-03,-0.42832809E-03,
         0.37643430E-03, 0.36649834E-03, 0.29694886E-03, 0.83649915E-03,
        -0.77822956E-03, 0.53860633E-04, 0.13068011E-02,-0.92486713E-04,
        -0.24946514E-03, 0.34415469E-04,-0.27364687E-03,-0.20774089E-03,
         0.31041173E-03, 0.15146896E-03,-0.31936832E-03, 0.51576545E-03,
         0.11680031E-03,-0.56778878E-03,-0.65380242E-03, 0.55188366E-03,
        -0.73630501E-04, 0.68980276E-04, 0.15793226E-03,-0.37967733E-04,
         0.10160460E-03, 0.98176468E-04,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_l5p77_q3ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]            *x41
        +coeff[  5]                *x51
        +coeff[  6]    *x21        *x51
        +coeff[  7]                *x52
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[  8]    *x22
        +coeff[  9]            *x42
        +coeff[ 10]*x11*x21
        +coeff[ 11]    *x21    *x41
        +coeff[ 12]        *x31*x41
        +coeff[ 13]    *x21    *x42
        +coeff[ 14]    *x21    *x41*x51
        +coeff[ 15]    *x21        *x52
        +coeff[ 16]                *x53
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 17]*x11            *x51
        +coeff[ 18]            *x41*x51
        +coeff[ 19]*x12
        +coeff[ 20]    *x21*x31
        +coeff[ 21]*x11        *x41
        +coeff[ 22]        *x31    *x51
        +coeff[ 23]*x11*x22
        +coeff[ 24]*x11        *x42
        +coeff[ 25]    *x21*x31    *x51
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 26]            *x42*x51
        +coeff[ 27]*x11*x23
        +coeff[ 28]    *x23    *x41
        +coeff[ 29]    *x22*x31*x41
        +coeff[ 30]    *x22    *x42
        +coeff[ 31]*x11    *x31
        +coeff[ 32]    *x23
        +coeff[ 33]*x11*x21*x31
        +coeff[ 34]    *x22*x31
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 35]    *x22    *x41
        +coeff[ 36]            *x43
        +coeff[ 37]*x11            *x52
        +coeff[ 38]    *x21    *x43
        +coeff[ 39]*x11*x22        *x51
        +coeff[ 40]    *x23        *x51
        +coeff[ 41]    *x21    *x42*x51
        +coeff[ 42]    *x23*x31*x41*x51
        +coeff[ 43]    *x23        *x53
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 44]*x12*x21
        +coeff[ 45]*x11    *x31*x41
        +coeff[ 46]        *x31*x42
        +coeff[ 47]        *x32    *x51
        +coeff[ 48]*x11        *x41*x51
        +coeff[ 49]        *x31*x41*x51
        ;

    return v_t_l5p77_q3ex                            ;
}
float y_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1058558E-01;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.22179575E-02, 0.30679686E-01,-0.61627096E-02,-0.20193111E-01,
        -0.95450354E-03,-0.18435284E-01,-0.83072186E-02,-0.42888533E-01,
        -0.62276549E-02, 0.27265819E-01,-0.44361525E-02, 0.77226786E-02,
         0.54597128E-01, 0.11403156E-01, 0.11477981E-02, 0.38964257E-02,
         0.62486329E-02,-0.11936365E-01,-0.13387307E-02,-0.69533517E-02,
        -0.34398612E-01,-0.55984054E-02,-0.78904703E-02, 0.12249307E-01,
         0.32222418E-02,-0.12169634E-02,-0.53675198E-02, 0.11461754E-02,
        -0.49003251E-02,-0.43293345E-02,-0.50430344E-02,-0.57921247E-02,
        -0.83772438E-02,-0.61638271E-02,-0.15636451E-02,-0.22577422E-02,
        -0.19449932E-02,-0.17002567E-02,-0.21911792E-02,-0.10933450E-02,
         0.28938407E-02,-0.30626529E-02,-0.41992767E-02,-0.19467649E-02,
         0.10838457E-02,-0.15528521E-02, 0.57471171E-03, 0.60760686E-02,
        -0.12052527E-02, 0.73381886E-03, 0.39888173E-02, 0.13675777E-02,
         0.11469406E-02, 0.53033972E-03, 0.78023202E-03,-0.85925480E-03,
        -0.22022398E-02, 0.51019182E-02,-0.20654474E-02,-0.35703639E-03,
        -0.34428609E-03, 0.33726142E-03,-0.26572737E-03, 0.15950733E-03,
         0.29174107E-03, 0.14874148E-02, 0.89725322E-03, 0.17171473E-02,
        -0.10860784E-02, 0.48243409E-03, 0.53771952E-03,-0.23856859E-03,
        -0.66809164E-03,-0.89019228E-03,-0.35958330E-02,-0.95348846E-03,
         0.31420346E-02, 0.37964934E-02,-0.18613223E-02,-0.26117082E-03,
        -0.23108204E-02, 0.16111142E-02,-0.11572446E-02, 0.75071666E-03,
         0.96173689E-03, 0.10593235E-02, 0.64477127E-03, 0.81696472E-03,
         0.41957755E-03,-0.46277660E-03,-0.96211123E-03,-0.49800478E-03,
        -0.23981890E-03,-0.44968061E-03,-0.23795356E-03, 0.11370225E-02,
         0.14593868E-02,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p77_q3ex                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]*x11
        +coeff[  5]                *x51
        +coeff[  6]            *x42
        +coeff[  7]    *x21    *x41
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[  8]        *x31*x41
        +coeff[  9]    *x22
        +coeff[ 10]    *x21*x31
        +coeff[ 11]*x11        *x41
        +coeff[ 12]            *x41*x51
        +coeff[ 13]*x11*x21
        +coeff[ 14]*x11    *x31
        +coeff[ 15]    *x21        *x51
        +coeff[ 16]        *x31    *x51
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]*x11            *x51
        +coeff[ 19]                *x52
        +coeff[ 20]    *x22    *x41
        +coeff[ 21]    *x23
        +coeff[ 22]    *x22*x31
        +coeff[ 23]    *x22        *x51
        +coeff[ 24]*x11*x21        *x51
        +coeff[ 25]        *x32
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 26]            *x43
        +coeff[ 27]*x12
        +coeff[ 28]    *x21*x31*x41
        +coeff[ 29]*x11        *x42
        +coeff[ 30]    *x21    *x41*x51
        +coeff[ 31]    *x22    *x42
        +coeff[ 32]    *x24
        +coeff[ 33]    *x22    *x41*x51
        +coeff[ 34]        *x31*x42
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]        *x31*x41*x51
        +coeff[ 37]*x11*x22
        +coeff[ 38]*x11*x21*x31
        +coeff[ 39]            *x41*x52
        +coeff[ 40]    *x23    *x41
        +coeff[ 41]*x11*x22    *x41
        +coeff[ 42]*x11*x23
        +coeff[ 43]    *x22*x31    *x51
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 44]*x11*x21    *x43
        +coeff[ 45]*x11*x21    *x41
        +coeff[ 46]    *x21*x31    *x51
        +coeff[ 47]    *x21    *x43
        +coeff[ 48]*x12        *x41
        +coeff[ 49]*x11        *x41*x51
        +coeff[ 50]    *x21*x31*x42
        +coeff[ 51]    *x21        *x52
        +coeff[ 52]    *x23*x31
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 53]*x12            *x51
        +coeff[ 54]                *x53
        +coeff[ 55]*x12*x21    *x41
        +coeff[ 56]*x11*x21    *x41*x51
        +coeff[ 57]    *x23    *x42
        +coeff[ 58]    *x23    *x42*x51
        +coeff[ 59]    *x21*x32
        +coeff[ 60]        *x32    *x51
        +coeff[ 61]*x12*x21
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 62]*x12    *x31
        +coeff[ 63]*x11    *x31    *x51
        +coeff[ 64]        *x31    *x52
        +coeff[ 65]*x11*x21    *x42
        +coeff[ 66]    *x22*x32
        +coeff[ 67]*x11*x21*x31*x41
        +coeff[ 68]    *x21*x31*x41*x51
        +coeff[ 69]*x12        *x42
        +coeff[ 70]    *x21    *x41*x52
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 71]            *x44*x51
        +coeff[ 72]*x12*x22
        +coeff[ 73]*x11*x21*x31    *x51
        +coeff[ 74]    *x24    *x41
        +coeff[ 75]            *x41*x53
        +coeff[ 76]*x11*x22    *x42
        +coeff[ 77]    *x22    *x42*x51
        +coeff[ 78]*x11*x23    *x41
        +coeff[ 79]*x11    *x34
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 80]    *x21    *x44*x51
        +coeff[ 81]    *x23*x32*x41
        +coeff[ 82]            *x42*x51
        +coeff[ 83]*x11        *x43
        +coeff[ 84]            *x43*x51
        +coeff[ 85]    *x22*x31*x41
        +coeff[ 86]*x11    *x31*x42
        +coeff[ 87]        *x31*x42*x51
        +coeff[ 88]*x11*x21*x32
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 89]    *x23        *x51
        +coeff[ 90]    *x22    *x43
        +coeff[ 91]    *x22        *x52
        +coeff[ 92]*x12        *x41*x51
        +coeff[ 93]*x11*x21        *x52
        +coeff[ 94]        *x31    *x53
        +coeff[ 95]*x11*x22*x31*x41
        +coeff[ 96]    *x22*x31*x41*x51
        ;

    return v_y_l5p77_q3ex                            ;
}
float p_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 90;
    float avdat= -0.6457438E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.25730138E-02, 0.32415991E-02, 0.60916599E-02,-0.30777715E-01,
         0.10645136E-01,-0.14892458E-01, 0.13714773E-02, 0.17593469E-01,
         0.36614637E-02,-0.23707608E-02,-0.38732507E-02,-0.20645456E-01,
        -0.61977007E-02, 0.14971299E-02,-0.35303810E-02, 0.15988667E-01,
         0.52631153E-02,-0.32376049E-02, 0.59840636E-03, 0.17585995E-02,
         0.41208844E-03,-0.17084695E-03,-0.23653502E-03,-0.22428195E-03,
        -0.99699944E-04, 0.28946073E-03, 0.51708549E-03, 0.49395226E-02,
         0.43350607E-02, 0.55417040E-03, 0.22564705E-02,-0.70590456E-03,
         0.28923708E-02, 0.74568909E-03,-0.47667092E-03, 0.20311850E-02,
         0.35040362E-02, 0.11528519E-02, 0.22485389E-02, 0.86866145E-03,
        -0.80523774E-03, 0.19639281E-02, 0.13967928E-02, 0.28594614E-02,
        -0.11667373E-03, 0.53261861E-03, 0.10454656E-03,-0.24772240E-03,
         0.80552074E-03,-0.22031036E-02, 0.57110889E-03,-0.79803192E-03,
         0.46215364E-03,-0.30708118E-03, 0.12361330E-03,-0.49875164E-03,
        -0.90202794E-03,-0.15844725E-02,-0.19764972E-03,-0.15012909E-02,
        -0.11753599E-03, 0.63503836E-03,-0.11226493E-03, 0.17362258E-03,
        -0.94413833E-03,-0.19837855E-02,-0.15234466E-03, 0.18428157E-03,
         0.57225529E-03, 0.14407429E-03,-0.12133894E-02, 0.48473550E-03,
         0.34731519E-03,-0.11481058E-02,-0.54794003E-03, 0.31638869E-04,
        -0.15531606E-02,-0.67631115E-03, 0.10818742E-03,-0.29264635E-03,
         0.83974919E-04,-0.74409711E-03,-0.20696620E-03,-0.10879805E-03,
         0.29472035E-03, 0.99772129E-04,-0.19626722E-03, 0.25435979E-03,
         0.71431481E-03, 0.20328704E-03,
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
    float x24 = x23*x2;
    float x25 = x24*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_q3ex                            =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21*x31
        +coeff[  7]    *x21    *x41
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[  8]        *x31*x41
        +coeff[  9]    *x21        *x51
        +coeff[ 10]        *x31    *x51
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]                *x52
        +coeff[ 14]*x11        *x41
        +coeff[ 15]    *x22    *x41
        +coeff[ 16]    *x21    *x42
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 17]    *x22        *x51
        +coeff[ 18]*x11*x22
        +coeff[ 19]*x11        *x42
        +coeff[ 20]    *x22    *x42
        +coeff[ 21]        *x32*x42
        +coeff[ 22]            *x44
        +coeff[ 23]    *x25
        +coeff[ 24]    *x24*x31
        +coeff[ 25]        *x32*x44
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 26]*x11
        +coeff[ 27]            *x42
        +coeff[ 28]    *x22*x31
        +coeff[ 29]*x11            *x51
        +coeff[ 30]            *x43
        +coeff[ 31]*x12
        +coeff[ 32]            *x41*x52
        +coeff[ 33]        *x32
        +coeff[ 34]*x11    *x31
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 35]    *x21*x31*x41
        +coeff[ 36]    *x24
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]*x11*x21    *x41
        +coeff[ 39]        *x31    *x52
        +coeff[ 40]*x11*x21        *x51
        +coeff[ 41]*x11*x23
        +coeff[ 42]*x11*x22    *x41
        +coeff[ 43]    *x23
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 44]        *x32*x41
        +coeff[ 45]        *x31*x42
        +coeff[ 46]        *x31*x41*x51
        +coeff[ 47]            *x42*x51
        +coeff[ 48]*x11    *x31*x41
        +coeff[ 49]    *x21    *x43
        +coeff[ 50]*x12        *x41
        +coeff[ 51]            *x41*x53
        +coeff[ 52]*x12*x21    *x41
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 53]    *x21    *x41*x51
        +coeff[ 54]    *x21        *x52
        +coeff[ 55]    *x23*x31
        +coeff[ 56]    *x23    *x41
        +coeff[ 57]    *x22*x31*x41
        +coeff[ 58]    *x21*x32*x41
        +coeff[ 59]    *x21*x31*x42
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]    *x22    *x41*x51
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 62]*x11            *x52
        +coeff[ 63]*x12    *x31
        +coeff[ 64]*x11*x21    *x42
        +coeff[ 65]    *x23    *x42
        +coeff[ 66]*x12            *x51
        +coeff[ 67]    *x23*x31    *x51
        +coeff[ 68]*x11*x21    *x41*x51
        +coeff[ 69]*x11    *x31*x41*x51
        +coeff[ 70]    *x22    *x42*x51
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 71]*x12*x22
        +coeff[ 72]    *x24*x32
        +coeff[ 73]*x11*x22    *x42
        +coeff[ 74]*x11*x21    *x43
        +coeff[ 75]*x12*x23
        +coeff[ 76]*x11*x23*x31*x41
        +coeff[ 77]            *x43*x53
        +coeff[ 78]    *x21*x32
        +coeff[ 79]    *x21*x31    *x51
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 80]*x11    *x32
        +coeff[ 81]    *x22*x32
        +coeff[ 82]    *x23        *x51
        +coeff[ 83]                *x53
        +coeff[ 84]    *x22*x31    *x51
        +coeff[ 85]*x11        *x41*x51
        +coeff[ 86]*x12*x21
        +coeff[ 87]    *x22        *x52
        +coeff[ 88]    *x24    *x41
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 89]    *x21*x31    *x52
        ;

    return v_p_l5p77_q3ex                            ;
}
float l_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 78;
    float avdat= -0.4338560E-02;
    float xmin[10]={
        -0.14986E-01,-0.46868E-01,-0.14992E-01,-0.32054E-01,-0.47828E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14968E-01, 0.45081E-01, 0.14987E-01, 0.25342E-01, 0.49612E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 79]={
         0.71849283E-02,-0.24639851E+00, 0.14241282E-01,-0.30855881E-01,
         0.41478883E-01,-0.13598329E-01,-0.88087566E-01, 0.98038828E-02,
        -0.20688597E-01, 0.10350027E-01,-0.45829423E-03, 0.53239274E-02,
        -0.32625038E-01,-0.90693235E-02,-0.84969122E-02,-0.79127112E-02,
        -0.22041723E-01,-0.26941504E-01, 0.37724920E-01,-0.16097754E-01,
         0.25898864E-04,-0.17698557E-02,-0.10937732E-01,-0.28302921E-02,
        -0.89718169E-02, 0.24589932E-01,-0.39913547E-02, 0.32921292E-01,
         0.21684222E-01,-0.42237467E-02, 0.13144779E-02, 0.43640877E-02,
        -0.58049872E-02, 0.61414000E-02, 0.51169022E-05,-0.83699539E-04,
         0.15277488E-02,-0.14662374E-02, 0.45826701E-02, 0.74211718E-02,
         0.44287140E-02, 0.62547154E-02,-0.55015826E-03,-0.10489679E-02,
         0.11932333E-02, 0.14857844E-02, 0.10467103E-02, 0.92386169E-03,
         0.88224764E-03,-0.42958171E-02, 0.60369377E-02,-0.11066964E-01,
        -0.94026364E-02,-0.32987741E-02,-0.95361471E-03, 0.16939449E-02,
         0.65062004E-02,-0.17758511E-03,-0.79541438E-03, 0.36694013E-03,
         0.14394934E-02, 0.88449335E-03,-0.34976361E-03, 0.26778673E-03,
         0.31074983E-03,-0.32169279E-03, 0.38621759E-02,-0.41084951E-02,
        -0.76264783E-03, 0.18944256E-02,-0.46442370E-02, 0.26539445E-02,
         0.22178500E-03, 0.11164590E-02, 0.15324955E-02, 0.20307316E-02,
        -0.13923170E-02, 0.37347181E-02,
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
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_l_l5p77_q3ex                            =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]*x11
        +coeff[  5]    *x22
        +coeff[  6]    *x21    *x41
        +coeff[  7]*x11*x21
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x23*x31
        +coeff[ 10]    *x21*x33
        +coeff[ 11]        *x31
        +coeff[ 12]    *x21*x31
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]                *x52
        +coeff[ 15]*x11    *x31
        +coeff[ 16]    *x23
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]*x11*x22
        +coeff[ 20]            *x44
        +coeff[ 21]    *x23*x31*x41
        +coeff[ 22]            *x42
        +coeff[ 23]*x11            *x51
        +coeff[ 24]    *x22*x31
        +coeff[ 25]    *x21*x31*x41
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 26]*x12*x21
        +coeff[ 27]    *x23    *x41
        +coeff[ 28]*x11*x22    *x41
        +coeff[ 29]        *x31*x41
        +coeff[ 30]            *x41*x51
        +coeff[ 31]    *x21*x32
        +coeff[ 32]*x11*x21    *x41
        +coeff[ 33]*x11        *x42
        +coeff[ 34]*x11    *x33*x41
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 35]*x11*x24*x31
        +coeff[ 36]                *x53
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]*x11    *x31*x41
        +coeff[ 39]*x11*x22*x31
        +coeff[ 40]*x12*x21    *x41
        +coeff[ 41]    *x23    *x44
        +coeff[ 42]        *x32
        +coeff[ 43]    *x22        *x51
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 44]    *x21*x31    *x51
        +coeff[ 45]    *x21    *x41*x51
        +coeff[ 46]            *x42*x51
        +coeff[ 47]*x11    *x32
        +coeff[ 48]*x11        *x41*x51
        +coeff[ 49]    *x24
        +coeff[ 50]    *x22    *x42
        +coeff[ 51]    *x21*x31*x42
        +coeff[ 52]    *x21    *x43
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 53]*x11*x23
        +coeff[ 54]*x12*x22
        +coeff[ 55]*x12*x21*x31
        +coeff[ 56]    *x24    *x41
        +coeff[ 57]        *x33*x41*x51
        +coeff[ 58]    *x24*x31*x41
        +coeff[ 59]    *x23*x32*x41
        +coeff[ 60]            *x43
        +coeff[ 61]        *x31*x41*x51
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 62]*x11*x21        *x51
        +coeff[ 63]*x11    *x31    *x51
        +coeff[ 64]*x11            *x52
        +coeff[ 65]*x13
        +coeff[ 66]    *x22*x31*x41
        +coeff[ 67]    *x21*x32*x41
        +coeff[ 68]*x11        *x42*x51
        +coeff[ 69]    *x24*x31
        +coeff[ 70]    *x23    *x42
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 71]    *x22*x31*x42
        +coeff[ 72]        *x33*x42
        +coeff[ 73]        *x31*x44
        +coeff[ 74]*x11*x24
        +coeff[ 75]*x11*x23    *x41
        +coeff[ 76]    *x23    *x41*x52
        +coeff[ 77]    *x22*x32*x43
        ;

    return v_l_l5p77_q3ex                            ;
}
float x_l5p77_sen                             (float *x,int m){
    int ncoeff=  8;
    float avdat=  0.1121842E+00;
    float xmin[10]={
        -0.14999E-01,-0.49847E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.50021E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[  9]={
        -0.59950207E-02, 0.40866030E-06,-0.22915381E-04, 0.85465103E-06,
         0.15126678E-01, 0.32083437E-01, 0.84466577E-04, 0.32757715E-04,
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
    float x21 = x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x51 = x5;

//                 function

    float v_x_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]    *x21
        +coeff[  4]        *x31
        +coeff[  5]            *x41
        +coeff[  6]            *x42
        +coeff[  7]        *x31*x41
        ;

    return v_x_l5p77_sen                             ;
}
float t_l5p77_sen                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1068508E+00;
    float xmin[10]={
        -0.14999E-01,-0.49847E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.50021E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
        -0.46903566E-02,-0.17249627E-04, 0.56506065E-03, 0.30643784E-01,
        -0.22659960E-03,-0.22669460E-03,-0.48594078E-03, 0.10692203E-03,
        -0.23483225E-03,-0.31906617E-03,-0.16288423E-05, 0.13750342E-03,
        -0.53752919E-04, 0.10447811E-03,-0.13101609E-03,-0.48742448E-04,
        -0.15435139E-03,-0.52097024E-04, 0.52003918E-04,-0.21249787E-04,
         0.38720522E-04, 0.45492310E-04,-0.12451738E-03, 0.99168261E-04,
         0.16143646E-03, 0.10039459E-03, 0.18232731E-03, 0.21010923E-04,
         0.10373779E-04,-0.33074968E-04,-0.27404467E-04, 0.34119981E-04,
        -0.42899159E-04, 0.93504663E-04, 0.19868510E-04, 0.23080813E-04,
         0.53977851E-05, 0.13604460E-04,-0.89625373E-05,-0.36887770E-04,
         0.25279989E-04, 0.19165431E-04,-0.46427125E-04,-0.19238279E-04,
         0.24044562E-04,-0.26413069E-04, 0.15448068E-04, 0.39204111E-04,
        -0.15926120E-04, 0.30366366E-04,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_t_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]        *x31*x41
        +coeff[  6]    *x22    *x41
        +coeff[  7]    *x22
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[  8]    *x22*x31
        +coeff[  9]*x11*x21    *x41
        +coeff[ 10]*x13*x21
        +coeff[ 11]*x11*x21
        +coeff[ 12]        *x32
        +coeff[ 13]    *x21    *x41
        +coeff[ 14]            *x42
        +coeff[ 15]            *x41*x51
        +coeff[ 16]*x11*x21*x31
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 17]    *x21
        +coeff[ 18]    *x21*x31
        +coeff[ 19]        *x31    *x51
        +coeff[ 20]*x11*x22
        +coeff[ 21]    *x23
        +coeff[ 22]*x11*x23
        +coeff[ 23]*x11*x21*x31*x41
        +coeff[ 24]    *x22*x31*x41
        +coeff[ 25]*x11*x21    *x42
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 26]    *x22    *x42
        +coeff[ 27]*x11        *x41
        +coeff[ 28]                *x52
        +coeff[ 29]*x12        *x41
        +coeff[ 30]    *x21    *x42
        +coeff[ 31]    *x22*x32
        +coeff[ 32]    *x23    *x41
        +coeff[ 33]*x11*x23    *x41
        +coeff[ 34]    *x23*x31*x41
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 35]*x12
        +coeff[ 36]*x11    *x31
        +coeff[ 37]*x12*x21
        +coeff[ 38]*x12    *x31
        +coeff[ 39]    *x21*x31*x41
        +coeff[ 40]        *x31*x42
        +coeff[ 41]            *x43
        +coeff[ 42]*x12*x22
        +coeff[ 43]    *x23*x31
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 44]*x11*x21*x32
        +coeff[ 45]*x11*x22    *x41
        +coeff[ 46]    *x22    *x41*x51
        +coeff[ 47]*x11*x23*x31
        +coeff[ 48]*x12*x21*x32
        +coeff[ 49]*x12*x22    *x41
        ;

    return v_t_l5p77_sen                             ;
}
float y_l5p77_sen                             (float *x,int m){
    int ncoeff= 11;
    float avdat= -0.1196481E-03;
    float xmin[10]={
        -0.14999E-01,-0.49847E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.50021E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 12]={
         0.20016705E-04,-0.54163944E-01,-0.14932141E-01,-0.21588321E-03,
        -0.10162621E-03, 0.51870611E-05,-0.16113601E-04,-0.76549832E-05,
        -0.75931175E-05, 0.23753817E-05,-0.10414150E-04,
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
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x31 = x3;
    float x41 = x4;
    float x51 = x5;

//                 function

    float v_y_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]*x11
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x21*x31
        +coeff[  5]            *x41
        +coeff[  6]*x11        *x41
        +coeff[  7]*x11    *x31
    ;
    v_y_l5p77_sen                             =v_y_l5p77_sen
        +coeff[  8]    *x21        *x51
        +coeff[  9]        *x31
        +coeff[ 10]    *x23
        ;

    return v_y_l5p77_sen                             ;
}
float p_l5p77_sen                             (float *x,int m){
    int ncoeff= 14;
    float avdat= -0.2192133E-03;
    float xmin[10]={
        -0.14999E-01,-0.49847E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.50021E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 15]={
        -0.26132944E-04,-0.48657309E-01, 0.48017025E-03,-0.96269342E-03,
        -0.37021007E-03,-0.19778428E-03,-0.17181183E-03, 0.56776087E-04,
         0.64195061E-04,-0.64422558E-04,-0.20590927E-03,-0.84063904E-04,
         0.19113971E-03, 0.26940284E-03,
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
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x51 = x5;

//                 function

    float v_p_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]*x11
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x21*x31
        +coeff[  5]*x11        *x41
        +coeff[  6]*x11*x22
        +coeff[  7]            *x41
    ;
    v_p_l5p77_sen                             =v_p_l5p77_sen
        +coeff[  8]    *x22
        +coeff[  9]    *x21        *x51
        +coeff[ 10]    *x23
        +coeff[ 11]*x11    *x31
        +coeff[ 12]    *x21    *x42
        +coeff[ 13]    *x23*x31*x41
        ;

    return v_p_l5p77_sen                             ;
}
float l_l5p77_sen                             (float *x,int m){
    int ncoeff= 10;
    float avdat= -0.8123883E-03;
    float xmin[10]={
        -0.14999E-01,-0.49847E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.50021E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 11]={
         0.11794835E-02,-0.55921623E-05,-0.15217921E-02,-0.31326029E-02,
        -0.13476285E-02,-0.50958124E-05,-0.47669525E-03, 0.23187715E-05,
         0.26669011E-05,-0.37908310E-05,
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
    float x21 = x2;
    float x22 = x21*x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x51 = x5;

//                 function

    float v_l_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]    *x22
        +coeff[  5]        *x31*x41
        +coeff[  6]            *x42
        +coeff[  7]                *x51
    ;
    v_l_l5p77_sen                             =v_l_l5p77_sen
        +coeff[  8]*x11*x21
        +coeff[  9]    *x22    *x41
        ;

    return v_l_l5p77_sen                             ;
}
float x_l5p77_sex                             (float *x,int m){
    int ncoeff= 58;
    float avdat=  0.2705576E+00;
    float xmin[10]={
        -0.14999E-01,-0.49683E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.49843E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 59]={
        -0.11199610E-01,-0.20569874E-03,-0.30679307E-02,-0.83651696E-03,
         0.65685511E-01, 0.51152189E-02,-0.26679403E-03,-0.13245286E-02,
        -0.14512676E-02,-0.18165293E-03,-0.17622302E-03, 0.79576921E-05,
         0.14071431E-04,-0.15136121E-02,-0.60955153E-04,-0.60231599E-04,
         0.11973791E-03,-0.17241333E-04, 0.17602968E-04,-0.15756763E-04,
         0.19691128E-03,-0.26003884E-05,-0.89653244E-04,-0.18620062E-03,
         0.17082481E-01,-0.40119826E-02, 0.24340858E-02,-0.18409061E-02,
         0.58797991E-03,-0.74518024E-03, 0.29354604E-03, 0.15938535E-03,
        -0.20759820E-03, 0.23565094E-03, 0.18127001E-02, 0.13365202E-03,
         0.11449235E-03,-0.67538203E-04,-0.22950124E-03, 0.25710129E-03,
        -0.25234083E-03, 0.78689127E-03, 0.11823756E-02, 0.18499339E-03,
        -0.22713250E-04,-0.47267946E-04, 0.41241343E-04,-0.11179689E-03,
         0.18445232E-03,-0.16436521E-03, 0.15192703E-03,-0.28597272E-03,
         0.60449319E-03, 0.60975646E-04, 0.18755965E-03, 0.13937723E-03,
        -0.17031522E-03, 0.14464698E-03,
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
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;
    float x54 = x53*x5;

//                 function

    float v_x_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]    *x21
        +coeff[  4]            *x41
        +coeff[  5]    *x22
        +coeff[  6]        *x32
        +coeff[  7]        *x31*x41
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[  8]            *x42
        +coeff[  9]*x12    *x31
        +coeff[ 10]*x12        *x41
        +coeff[ 11]        *x31    *x52
        +coeff[ 12]            *x41*x52
        +coeff[ 13]    *x22*x31
        +coeff[ 14]        *x33
        +coeff[ 15]*x13*x21
        +coeff[ 16]*x14    *x31
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 17]*x12    *x31    *x52
        +coeff[ 18]        *x31    *x54
        +coeff[ 19]*x12*x22*x31
        +coeff[ 20]*x12    *x33
        +coeff[ 21]        *x33    *x52
        +coeff[ 22]    *x22*x33
        +coeff[ 23]*x14    *x33
        +coeff[ 24]        *x31
        +coeff[ 25]    *x22    *x41
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 26]*x11*x21
        +coeff[ 27]*x11*x21    *x41
        +coeff[ 28]    *x21    *x41
        +coeff[ 29]*x11*x21*x31
        +coeff[ 30]*x12
        +coeff[ 31]                *x52
        +coeff[ 32]            *x41*x51
        +coeff[ 33]    *x21*x31
        +coeff[ 34]    *x22    *x42
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 35]    *x22*x33*x41
        +coeff[ 36]*x11        *x41
        +coeff[ 37]        *x31    *x51
        +coeff[ 38]    *x22        *x51
        +coeff[ 39]    *x23
        +coeff[ 40]    *x21    *x42
        +coeff[ 41]*x11*x21    *x42
        +coeff[ 42]    *x22*x31*x41
        +coeff[ 43]        *x31*x43
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 44]*x13*x22
        +coeff[ 45]*x13*x21*x31*x41
        +coeff[ 46]*x11    *x31
        +coeff[ 47]*x11*x21        *x51
        +coeff[ 48]*x11*x22
        +coeff[ 49]    *x21*x31*x41
        +coeff[ 50]            *x43
        +coeff[ 51]*x11*x23
        +coeff[ 52]*x11*x21*x31*x41
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 53]    *x23        *x51
        +coeff[ 54]    *x22*x32
        +coeff[ 55]        *x32*x42
        +coeff[ 56]*x14*x22
        +coeff[ 57]        *x33*x42
        ;

    return v_x_l5p77_sex                             ;
}
float t_l5p77_sex                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.2259126E+00;
    float xmin[10]={
        -0.14999E-01,-0.49683E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.49843E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
        -0.71544354E-02, 0.25543002E-02, 0.37448425E-01,-0.62721171E-02,
         0.99104010E-02, 0.55724304E-04, 0.48390368E-03,-0.60237884E-02,
        -0.13513529E-04,-0.14841382E-02, 0.41275397E-02,-0.23652981E-02,
        -0.30198272E-02,-0.26780327E-02,-0.33426957E-03, 0.41659392E-03,
         0.85010962E-03,-0.28997884E-03, 0.33518113E-03,-0.22720557E-02,
         0.37032925E-02,-0.16055748E-03, 0.20358572E-03, 0.35200789E-03,
        -0.42488310E-03, 0.31751391E-03,-0.93742879E-03,-0.47473723E-03,
         0.22572228E-02, 0.14571224E-02,-0.15701979E-04, 0.52322725E-05,
         0.41843414E-04,-0.23798327E-03,-0.44718612E-03, 0.48813265E-03,
        -0.22417474E-03,-0.30532086E-03, 0.33847959E-03, 0.99256262E-03,
        -0.10888204E-03, 0.15562111E-04,-0.10531847E-03, 0.51187177E-04,
         0.17290714E-03, 0.59980997E-04,-0.84350737E-04,-0.25235850E-03,
         0.17058579E-03, 0.50625362E-03,
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
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]*x12*x21
        +coeff[  6]    *x23
        +coeff[  7]    *x22    *x41
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[  8]*x13*x21
        +coeff[  9]    *x21
        +coeff[ 10]*x11*x21
        +coeff[ 11]        *x31*x41
        +coeff[ 12]            *x42
        +coeff[ 13]*x11*x21    *x41
        +coeff[ 14]*x11
        +coeff[ 15]*x12
        +coeff[ 16]    *x21    *x41
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 17]            *x41*x51
        +coeff[ 18]                *x52
        +coeff[ 19]    *x22*x31
        +coeff[ 20]    *x22    *x42
        +coeff[ 21]*x13*x21*x31
        +coeff[ 22]    *x22*x33*x41
        +coeff[ 23]    *x21*x31
        +coeff[ 24]        *x32
        +coeff[ 25]*x11*x22
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 26]*x11*x21*x31
        +coeff[ 27]    *x22        *x51
        +coeff[ 28]    *x22*x31*x41
        +coeff[ 29]*x11*x21    *x42
        +coeff[ 30]*x11        *x43
        +coeff[ 31]        *x31*x43
        +coeff[ 32]*x13*x21*x31*x41
        +coeff[ 33]*x12        *x41
        +coeff[ 34]    *x21    *x42
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 35]            *x43
        +coeff[ 36]*x11*x21        *x51
        +coeff[ 37]*x11*x23
        +coeff[ 38]    *x22*x32
        +coeff[ 39]*x11*x21*x31*x41
        +coeff[ 40]    *x21*x31*x42
        +coeff[ 41]        *x32*x42
        +coeff[ 42]        *x31    *x53
        +coeff[ 43]*x11    *x31
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 44]*x11        *x41
        +coeff[ 45]    *x21        *x51
        +coeff[ 46]*x12    *x31
        +coeff[ 47]    *x21*x31*x41
        +coeff[ 48]        *x32*x41
        +coeff[ 49]        *x31*x42
        ;

    return v_t_l5p77_sex                             ;
}
float y_l5p77_sex                             (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.3898993E-03;
    float xmin[10]={
        -0.14999E-01,-0.49683E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.49843E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.31333332E-03,-0.95927365E-01,-0.13344292E-01,-0.67471084E-02,
        -0.25549044E-02,-0.14667620E-02, 0.48085253E-03,-0.58411493E-03,
         0.25805535E-02, 0.20803334E-02,-0.22021823E-02,-0.15503650E-02,
         0.19275278E-03, 0.54999371E-03, 0.21529088E-02, 0.25742489E-03,
        -0.26413810E-03,-0.48118475E-03, 0.58896979E-03, 0.40860483E-03,
         0.47369048E-03, 0.27950740E-03,-0.37166709E-03, 0.14090784E-02,
        -0.99535020E-04,-0.19080080E-03,-0.14814736E-03,-0.73010378E-04,
        -0.17837624E-03,-0.21350892E-03,-0.70702942E-03,-0.85323345E-03,
         0.79870952E-03, 0.54284959E-03, 0.30458777E-03, 0.25737034E-04,
         0.28523333E-04,-0.82592735E-04, 0.94811905E-04, 0.96128206E-04,
         0.60599607E-04,-0.29039904E-03, 0.12533113E-03,-0.32102951E-04,
        -0.14944194E-04,-0.16162332E-03,-0.18528460E-03,-0.66152599E-04,
        -0.26017255E-04,-0.59771410E-03,-0.49145479E-03,-0.40332490E-03,
        -0.34270596E-03, 0.10783207E-04,-0.78333233E-05,-0.14922450E-04,
        -0.59564682E-05, 0.60847484E-04,-0.19733687E-04,-0.16395396E-04,
         0.18204211E-04, 0.48009802E-04,-0.92465729E-04,-0.24553334E-04,
        -0.21956030E-04,-0.44121116E-04,-0.83816325E-04, 0.90764413E-04,
        -0.20823092E-04, 0.64521293E-04, 0.14818413E-04, 0.99544799E-04,
         0.58225862E-04, 0.77097116E-04, 0.21148420E-04,-0.10422384E-04,
         0.10970716E-03,-0.82377510E-05, 0.95568153E-04, 0.18048286E-04,
         0.35653095E-04,-0.50636820E-04, 0.26424099E-04,-0.14521114E-04,
        -0.18354242E-04, 0.18195195E-04, 0.71814226E-04,-0.89043446E-04,
         0.48556907E-04,-0.77757213E-04,-0.68447727E-04,-0.57616497E-04,
        -0.62678686E-04,-0.45326302E-04, 0.14882456E-03, 0.80348880E-04,
         0.11977318E-03,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;

//                 function

    float v_y_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]*x11
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x21*x31
        +coeff[  5]*x11        *x41
        +coeff[  6]            *x41
        +coeff[  7]*x11    *x31
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[  8]    *x21    *x42
        +coeff[  9]    *x21*x31*x41
        +coeff[ 10]    *x23
        +coeff[ 11]*x11*x22
        +coeff[ 12]        *x31
        +coeff[ 13]    *x22
        +coeff[ 14]    *x23    *x41
        +coeff[ 15]*x11*x21
        +coeff[ 16]    *x21        *x51
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]*x11        *x42
        +coeff[ 19]    *x21*x32
        +coeff[ 20]*x11    *x31*x41
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]*x12*x21
        +coeff[ 23]*x11*x22    *x41
        +coeff[ 24]    *x21*x31*x43
        +coeff[ 25]            *x42
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 26]        *x31*x41
        +coeff[ 27]*x11            *x51
        +coeff[ 28]    *x22*x31
        +coeff[ 29]*x11*x21    *x41
        +coeff[ 30]    *x21    *x43
        +coeff[ 31]    *x21*x31*x42
        +coeff[ 32]    *x23*x31
        +coeff[ 33]*x11*x22*x31
        +coeff[ 34]*x12*x21    *x41
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 35]                *x51
        +coeff[ 36]*x12
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]*x11    *x32
        +coeff[ 39]    *x21*x31    *x51
        +coeff[ 40]*x11        *x41*x51
        +coeff[ 41]    *x21*x32*x41
        +coeff[ 42]*x12*x21*x31
        +coeff[ 43]        *x32
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 44]            *x41*x51
        +coeff[ 45]*x11        *x43
        +coeff[ 46]*x11    *x31*x42
        +coeff[ 47]    *x21    *x42*x51
        +coeff[ 48]*x13
        +coeff[ 49]    *x23    *x42
        +coeff[ 50]    *x23*x31*x41
        +coeff[ 51]*x11*x22    *x42
        +coeff[ 52]*x11*x22*x31*x41
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 53]    *x22    *x44
        +coeff[ 54]    *x22*x31*x43
        +coeff[ 55]*x11    *x34*x41
        +coeff[ 56]        *x31    *x51
        +coeff[ 57]        *x31*x42
        +coeff[ 58]    *x22        *x51
        +coeff[ 59]*x12        *x41
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]*x11*x21    *x42
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 62]    *x24
        +coeff[ 63]    *x21*x33
        +coeff[ 64]            *x45
        +coeff[ 65]    *x21*x31*x41*x51
        +coeff[ 66]*x11*x23
        +coeff[ 67]    *x23        *x51
        +coeff[ 68]*x12*x22
        +coeff[ 69]*x11*x22        *x51
        +coeff[ 70]*x11*x21*x31*x42
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 71]*x11*x24
        +coeff[ 72]*x12*x23
        +coeff[ 73]            *x43
        +coeff[ 74]        *x32*x41
        +coeff[ 75]        *x31*x43
        +coeff[ 76]    *x22    *x42
        +coeff[ 77]*x11*x21        *x51
        +coeff[ 78]    *x22*x31*x41
        +coeff[ 79]    *x22*x32
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 80]*x11*x21*x31*x41
        +coeff[ 81]*x11    *x32*x41
        +coeff[ 82]    *x21    *x44
        +coeff[ 83]*x11        *x42*x51
        +coeff[ 84]    *x22    *x43
        +coeff[ 85]*x13        *x41
        +coeff[ 86]    *x24    *x41
        +coeff[ 87]    *x23*x32
        +coeff[ 88]*x11*x23    *x41
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 89]    *x23    *x41*x51
        +coeff[ 90]*x12*x21    *x42
        +coeff[ 91]*x11*x22*x32
        +coeff[ 92]*x12*x21*x31*x41
        +coeff[ 93]*x11*x22    *x41*x51
        +coeff[ 94]    *x23*x31*x42
        +coeff[ 95]*x11*x22*x31*x42
        +coeff[ 96]    *x21*x31*x45
        ;

    return v_y_l5p77_sex                             ;
}
float p_l5p77_sex                             (float *x,int m){
    int ncoeff= 40;
    float avdat= -0.5540310E-03;
    float xmin[10]={
        -0.14999E-01,-0.49683E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.49843E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 41]={
        -0.30882028E-03,-0.40416207E-01, 0.23580936E-02,-0.10608178E-01,
        -0.22275848E-02, 0.13287743E-02, 0.73454704E-03,-0.38240443E-02,
        -0.34454288E-02, 0.34888810E-02,-0.22434127E-02, 0.73567685E-03,
        -0.82587689E-03, 0.26309830E-02, 0.39949236E-02, 0.23135561E-02,
         0.27063364E-03,-0.34426892E-03, 0.31844154E-03,-0.75651688E-03,
         0.57793548E-03,-0.49160630E-03,-0.18030190E-03,-0.24915001E-03,
         0.47076668E-03, 0.35507811E-03,-0.29173755E-03, 0.39834643E-03,
        -0.12752295E-02,-0.12766715E-02, 0.18166462E-04, 0.80849958E-03,
         0.44213713E-03, 0.26718542E-04,-0.13237182E-03,-0.79171448E-04,
         0.12101854E-03,-0.10261720E-03,-0.40741882E-03, 0.17507910E-03,
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

    float v_p_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]*x11
        +coeff[  3]    *x21    *x41
        +coeff[  4]*x11        *x41
        +coeff[  5]    *x23*x31
        +coeff[  6]            *x41
        +coeff[  7]    *x21*x31
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[  8]    *x23
        +coeff[  9]    *x21    *x42
        +coeff[ 10]*x11*x22
        +coeff[ 11]    *x22
        +coeff[ 12]*x11    *x31
        +coeff[ 13]    *x21*x31*x41
        +coeff[ 14]    *x23    *x41
        +coeff[ 15]*x11*x22    *x41
        +coeff[ 16]        *x31
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 17]    *x21        *x51
        +coeff[ 18]*x11*x21
        +coeff[ 19]    *x22    *x41
        +coeff[ 20]*x11        *x42
        +coeff[ 21]*x12*x21
        +coeff[ 22]            *x42
        +coeff[ 23]    *x22*x31
        +coeff[ 24]    *x21*x32
        +coeff[ 25]    *x21    *x41*x51
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 26]*x11*x21    *x41
        +coeff[ 27]*x11    *x31*x41
        +coeff[ 28]    *x21*x31*x42
        +coeff[ 29]    *x21    *x43
        +coeff[ 30]        *x33*x41
        +coeff[ 31]*x11*x22*x31
        +coeff[ 32]*x12*x21    *x41
        +coeff[ 33]                *x51
        +coeff[ 34]        *x31*x41
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 35]*x11            *x51
        +coeff[ 36]    *x21*x31    *x51
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]    *x21*x32*x41
        +coeff[ 39]*x12*x21*x31
        ;

    return v_p_l5p77_sex                             ;
}
float l_l5p77_sex                             (float *x,int m){
    int ncoeff= 46;
    float avdat= -0.1379200E-02;
    float xmin[10]={
        -0.14999E-01,-0.49683E-01,-0.14992E-01,-0.32620E-01,-0.49954E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.14973E-01, 0.49843E-01, 0.14999E-01, 0.25877E-01, 0.49961E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 47]={
         0.23682129E-02, 0.10633060E-03,-0.18417379E-02,-0.84824935E-02,
         0.54767611E-03,-0.32232900E-02, 0.17544199E-03,-0.75099524E-03,
        -0.36929853E-03, 0.14053189E-03, 0.35184632E-04,-0.63947336E-04,
        -0.36515714E-04,-0.53881380E-04, 0.28421081E-03, 0.19951603E-03,
        -0.14435075E-04,-0.11302230E-04,-0.33642558E-04, 0.46751506E-04,
         0.53582582E-04, 0.16381313E-04,-0.16309905E-04, 0.18402630E-03,
         0.11214067E-03, 0.22964985E-04, 0.26846303E-04,-0.14572700E-03,
        -0.20822475E-03,-0.80565420E-04,-0.74970383E-04,-0.76661927E-05,
        -0.75883854E-05, 0.20409439E-04, 0.42941163E-04,-0.22240800E-04,
        -0.14937894E-04,-0.12747203E-04, 0.11556363E-04,-0.33272438E-04,
         0.51644649E-04,-0.16050777E-04, 0.24630661E-04,-0.13744141E-04,
        -0.42385207E-04, 0.95695868E-04,
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
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_l_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]        *x31*x41
        +coeff[  7]            *x42
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[  8]*x11*x21
        +coeff[  9]            *x41*x51
        +coeff[ 10]*x11
        +coeff[ 11]    *x21    *x41
        +coeff[ 12]                *x52
        +coeff[ 13]*x12
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]*x11*x21    *x41
        +coeff[ 16]    *x24*x31
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 17]*x11*x23*x31
        +coeff[ 18]    *x21*x31
        +coeff[ 19]        *x32
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]        *x31    *x51
        +coeff[ 22]*x11        *x41
        +coeff[ 23]    *x22*x31
        +coeff[ 24]*x11*x21*x31
        +coeff[ 25]*x11*x21        *x51
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 26]*x12        *x41
        +coeff[ 27]    *x22*x31*x41
        +coeff[ 28]    *x22    *x42
        +coeff[ 29]*x11*x21*x31*x41
        +coeff[ 30]*x11*x21    *x42
        +coeff[ 31]    *x21        *x51
        +coeff[ 32]*x11    *x31
        +coeff[ 33]    *x21*x31*x41
        +coeff[ 34]    *x21    *x42
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 35]        *x31*x42
        +coeff[ 36]            *x43
        +coeff[ 37]*x11*x22
        +coeff[ 38]*x12    *x31
        +coeff[ 39]    *x22*x32
        +coeff[ 40]*x11*x23
        +coeff[ 41]*x11*x21*x32
        +coeff[ 42]*x12*x22
        +coeff[ 43]        *x34*x41
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 44]    *x23    *x42
        +coeff[ 45]    *x24    *x42
        ;

    return v_l_l5p77_sex                             ;
}
}
