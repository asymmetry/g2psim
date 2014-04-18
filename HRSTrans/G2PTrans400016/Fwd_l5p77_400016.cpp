#include "Fwd_l5p77_400016.h"

namespace S400016
{
float x_l5p77_den                             (float *x,int m){
    int ncoeff= 60;
    float avdat= -0.5095617E+01;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
        -0.39391336E-03, 0.14594869E-01,-0.98443374E-01,-0.29902651E-02,
        -0.55198395E-02,-0.37073061E-01,-0.33051023E-04, 0.64491917E-03,
         0.10607323E-04, 0.47681932E-02,-0.10381549E-02,-0.12330023E-02,
        -0.30083882E-02,-0.81148390E-02,-0.13310478E-01,-0.96069947E-02,
         0.14922896E-01,-0.11121973E-03,-0.58841892E-03,-0.56498905E-03,
        -0.12737103E-02,-0.60985454E-04,-0.16504895E-02,-0.68726586E-02,
         0.98595377E-02, 0.14864930E-01,-0.79895828E-04, 0.10819768E-02,
         0.26094015E-02, 0.16933709E-02, 0.89315064E-02, 0.17726877E-03,
         0.76293807E-04, 0.47479090E-03, 0.16850334E-02, 0.96860871E-03,
         0.16919863E-02, 0.29824150E-02,-0.40137745E-02,-0.50062803E-03,
         0.18458509E-03, 0.35842805E-03, 0.40627914E-03, 0.30998443E-03,
         0.34187993E-03, 0.88606519E-03, 0.33939481E-03,-0.15368874E-03,
        -0.15021549E-02,-0.38812228E-02, 0.15777450E-04, 0.50777719E-04,
         0.80041507E-04, 0.75129465E-04, 0.15210491E-03, 0.34161718E-03,
        -0.46483910E-03,-0.44922411E-03,-0.76409831E-03,-0.47123319E-03,
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

    float v_x_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]            *x41
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x21    *x41
        +coeff[  6]*x13        *x41
        +coeff[  7]*x12*x21*x31
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[  8]    *x21*x31    *x52
        +coeff[  9]    *x23*x31
        +coeff[ 10]        *x31
        +coeff[ 11]*x11            *x51
        +coeff[ 12]*x11    *x31
        +coeff[ 13]*x11        *x41
        +coeff[ 14]    *x21*x31
        +coeff[ 15]    *x23
        +coeff[ 16]    *x21    *x42
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 17]*x13*x22
        +coeff[ 18]    *x23*x31*x41
        +coeff[ 19]                *x51
        +coeff[ 20]*x11*x21
        +coeff[ 21]*x13
        +coeff[ 22]*x12*x21
        +coeff[ 23]*x11*x22
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]    *x23    *x41
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 26]*x13*x22    *x41
        +coeff[ 27]    *x22
        +coeff[ 28]*x11        *x42
        +coeff[ 29]    *x21*x32
        +coeff[ 30]*x11*x22    *x41
        +coeff[ 31]*x13    *x31*x41
        +coeff[ 32]*x11*x22*x33
        +coeff[ 33]    *x21        *x52
        +coeff[ 34]*x11    *x31*x41
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 35]    *x22    *x41
        +coeff[ 36]*x12*x21    *x41
        +coeff[ 37]*x11*x22*x31
        +coeff[ 38]    *x21    *x43
        +coeff[ 39]    *x21*x33*x42
        +coeff[ 40]            *x41*x51
        +coeff[ 41]        *x31*x41
        +coeff[ 42]            *x42
        +coeff[ 43]*x11    *x32
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 44]    *x22        *x51
        +coeff[ 45]    *x21    *x41*x51
        +coeff[ 46]    *x22*x31
        +coeff[ 47]    *x21*x33
        +coeff[ 48]    *x21*x32*x41
        +coeff[ 49]    *x21*x31*x42
        +coeff[ 50]*x12*x21*x31    *x51
        +coeff[ 51]        *x31    *x51
        +coeff[ 52]        *x32
    ;
    v_x_l5p77_den                             =v_x_l5p77_den
        +coeff[ 53]*x11            *x52
        +coeff[ 54]*x11        *x41*x51
        +coeff[ 55]    *x21*x31    *x51
        +coeff[ 56]*x11        *x43
        +coeff[ 57]    *x21    *x42*x51
        +coeff[ 58]    *x23    *x41*x51
        +coeff[ 59]*x13    *x31*x42
        ;

    return v_x_l5p77_den                             ;
}
float t_l5p77_den                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.3590667E+01;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
        -0.47305722E-01,-0.11644359E+00, 0.58043706E+00, 0.22119245E-01,
         0.11082673E-01,-0.37915781E-01, 0.15746072E+00, 0.26802993E+00,
         0.37279040E-01, 0.99749528E-01,-0.28387114E-03, 0.10364122E+00,
        -0.42604599E-02,-0.22951603E-01, 0.97228121E-03,-0.25243426E-02,
         0.70469556E-02, 0.96835606E-01, 0.52948669E-01, 0.23109105E-02,
         0.35150174E-01, 0.34608461E-01,-0.96405603E-01,-0.48335623E-02,
        -0.41790190E-02, 0.46786959E-02, 0.19076077E-01, 0.28815538E-01,
         0.81837382E-02,-0.92981569E-02,-0.56759991E-01,-0.73623046E-01,
        -0.10648450E-03, 0.15661165E-01,-0.87236864E-02,-0.15708145E-01,
         0.10547618E-01,-0.12027637E-01,-0.59918571E-01,-0.37030693E-01,
        -0.12344376E-02,-0.20673354E-02, 0.61096381E-02,-0.93644634E-02,
         0.40493268E-02,-0.55102701E-04,-0.18091626E-01,-0.16001612E-01,
         0.90503572E-02,-0.69546537E-02,
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

    float v_t_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]*x11*x21
        +coeff[  6]    *x22
        +coeff[  7]    *x21    *x41
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[  8]    *x21        *x51
        +coeff[  9]    *x23
        +coeff[ 10]        *x33
        +coeff[ 11]    *x22    *x41
        +coeff[ 12]*x12*x21*x31
        +coeff[ 13]    *x23*x31
        +coeff[ 14]*x13        *x41
        +coeff[ 15]*x13*x22*x31
        +coeff[ 16]        *x31
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 17]    *x21*x31
        +coeff[ 18]*x11        *x41
        +coeff[ 19]        *x31*x41
        +coeff[ 20]*x11*x22
        +coeff[ 21]    *x22*x31
        +coeff[ 22]    *x21    *x42
        +coeff[ 23]*x12*x21*x31*x41
        +coeff[ 24]    *x23*x31*x41
        +coeff[ 25]*x12
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 26]*x11    *x31
        +coeff[ 27]            *x42
        +coeff[ 28]*x11            *x51
        +coeff[ 29]            *x41*x51
        +coeff[ 30]    *x21*x31*x41
        +coeff[ 31]    *x23    *x41
        +coeff[ 32]*x13*x22    *x41
        +coeff[ 33]*x12*x21
        +coeff[ 34]    *x21*x32
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 35]*x11        *x42
        +coeff[ 36]    *x22        *x51
        +coeff[ 37]*x12*x21    *x41
        +coeff[ 38]*x11*x22    *x41
        +coeff[ 39]    *x22    *x42
        +coeff[ 40]*x13    *x31*x41
        +coeff[ 41]*x12        *x41
        +coeff[ 42]*x11*x21    *x41
        +coeff[ 43]*x11    *x31*x41
    ;
    v_t_l5p77_den                             =v_t_l5p77_den
        +coeff[ 44]            *x42*x51
        +coeff[ 45]*x13    *x31
        +coeff[ 46]*x11*x22*x31
        +coeff[ 47]    *x22*x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x12*x21        *x52
        ;

    return v_t_l5p77_den                             ;
}
float y_l5p77_den                             (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1791847E-01;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.12602125E-01, 0.13989040E+00, 0.43638372E-02, 0.12958188E-01,
        -0.21742608E-01, 0.34156125E-01, 0.16846715E-01, 0.14548969E-01,
        -0.22372723E-01, 0.10209350E-02,-0.12438775E-01, 0.27758086E-02,
        -0.89389598E-02, 0.21306043E-02,-0.79715578E-02,-0.95815230E-02,
        -0.16443621E-02, 0.16384446E-02,-0.11419749E-02,-0.37181349E-02,
         0.22315474E-02, 0.54137176E-03,-0.94645051E-03, 0.26007029E-03,
        -0.13194820E-02, 0.11441636E-01, 0.87407697E-02, 0.47717402E-02,
        -0.56610932E-02, 0.36551401E-02,-0.45036939E-02,-0.25790869E-03,
         0.71281429E-04,-0.76401123E-03,-0.11662419E-02, 0.81825792E-03,
        -0.14357899E-02,-0.28215332E-04, 0.74342859E-03,-0.43700359E-03,
         0.15415063E-02, 0.21732456E-03,-0.91998430E-03, 0.70118561E-03,
         0.55506651E-03, 0.13823096E-02, 0.91891452E-04, 0.27129016E-03,
         0.76814718E-03, 0.31708382E-03,-0.12876265E-03, 0.13137657E-03,
        -0.20712179E-03, 0.11141366E-03,-0.42675322E-03,-0.37539294E-02,
        -0.55407960E-03, 0.11142380E-03,-0.30913444E-02,-0.35628112E-03,
        -0.11735051E-02, 0.15770725E-02, 0.15253260E-02,-0.18787060E-03,
         0.12122162E-03, 0.10057618E-03,-0.26002855E-03, 0.62793668E-04,
         0.42800562E-03, 0.21378913E-03,-0.37250359E-03, 0.55516831E-03,
         0.40897331E-03,-0.21147181E-03,-0.49410749E-03,-0.17396694E-03,
         0.68791298E-03, 0.74927567E-03,-0.32839936E-03,-0.20525046E-02,
        -0.10485058E-02, 0.69484947E-03, 0.35399702E-03, 0.14084046E-03,
         0.33773755E-04,-0.10267392E-03, 0.22354217E-03,-0.74281612E-04,
        -0.39288145E-04,-0.13750546E-03, 0.84353735E-04, 0.45811868E-03,
        -0.18732823E-03,-0.26176113E-03,-0.63498120E-03,-0.35175854E-04,
         0.35869787E-03,
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

    float v_y_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]            *x41*x51
        +coeff[  7]*x11*x21
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[  8]    *x22    *x41
        +coeff[  9]*x11
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21    *x41
        +coeff[ 12]        *x31*x41
        +coeff[ 13]        *x31    *x51
        +coeff[ 14]    *x22*x31
        +coeff[ 15]*x11*x21    *x41
        +coeff[ 16]        *x32
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 17]*x12
        +coeff[ 18]                *x52
        +coeff[ 19]*x11*x21*x31
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]    *x21*x31
        +coeff[ 22]*x11        *x41
        +coeff[ 23]            *x43
        +coeff[ 24]            *x41*x52
        +coeff[ 25]    *x22    *x42
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 26]    *x22*x31*x41
        +coeff[ 27]*x11*x21    *x42
        +coeff[ 28]    *x24
        +coeff[ 29]*x11*x21*x31*x41
        +coeff[ 30]*x11*x23
        +coeff[ 31]*x11    *x31
        +coeff[ 32]            *x42*x51
        +coeff[ 33]    *x21    *x41*x51
        +coeff[ 34]*x12        *x41
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 35]*x11*x21        *x51
        +coeff[ 36]*x12*x22
        +coeff[ 37]    *x21        *x51
        +coeff[ 38]    *x23
        +coeff[ 39]*x12    *x31
        +coeff[ 40]    *x22*x32
        +coeff[ 41]                *x53
        +coeff[ 42]    *x22    *x41*x51
        +coeff[ 43]*x11*x21*x32
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 44]    *x21    *x42
        +coeff[ 45]        *x31*x42
        +coeff[ 46]*x11            *x51
        +coeff[ 47]    *x21*x31*x41
        +coeff[ 48]        *x32*x41
        +coeff[ 49]*x11*x22
        +coeff[ 50]    *x21*x31    *x51
        +coeff[ 51]*x11        *x41*x51
        +coeff[ 52]        *x31    *x52
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 53]*x12            *x51
        +coeff[ 54]    *x22*x31    *x51
        +coeff[ 55]    *x22    *x43
        +coeff[ 56]*x11*x21    *x41*x51
        +coeff[ 57]*x11        *x44
        +coeff[ 58]    *x22*x31*x42
        +coeff[ 59]    *x22        *x52
        +coeff[ 60]    *x22*x32*x41
        +coeff[ 61]    *x22    *x42*x51
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 62]*x11*x23    *x41
        +coeff[ 63]*x13*x21        *x52
        +coeff[ 64]        *x33
        +coeff[ 65]*x11    *x31*x41
        +coeff[ 66]        *x31*x41*x51
        +coeff[ 67]    *x21        *x52
        +coeff[ 68]            *x43*x51
        +coeff[ 69]    *x23    *x41
        +coeff[ 70]            *x45
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 71]*x12        *x42
        +coeff[ 72]*x12    *x31*x41
        +coeff[ 73]*x11*x21*x31    *x51
        +coeff[ 74]*x11*x21*x31*x42
        +coeff[ 75]*x13*x21
        +coeff[ 76]    *x22*x31*x41*x51
        +coeff[ 77]*x11*x23*x31
        +coeff[ 78]    *x21    *x44*x51
        +coeff[ 79]    *x24    *x42
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 80]    *x24*x31*x41
        +coeff[ 81]*x12*x24    *x41
        +coeff[ 82]*x12*x24*x31
        +coeff[ 83]*x11        *x42
        +coeff[ 84]*x12*x21
        +coeff[ 85]*x11        *x43
        +coeff[ 86]        *x31*x42*x51
        +coeff[ 87]    *x23*x31
        +coeff[ 88]*x11            *x52
    ;
    v_y_l5p77_den                             =v_y_l5p77_den
        +coeff[ 89]    *x23        *x51
        +coeff[ 90]*x12    *x32
        +coeff[ 91]    *x24*x31
        +coeff[ 92]    *x22*x33
        +coeff[ 93]*x11*x21*x32*x41
        +coeff[ 94]    *x22    *x44
        +coeff[ 95]        *x34    *x51
        +coeff[ 96]    *x22*x32*x42
        ;

    return v_y_l5p77_den                             ;
}
float p_l5p77_den                             (float *x,int m){
    int ncoeff= 71;
    float avdat= -0.4185514E-02;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 72]={
         0.71497872E-02,-0.10865224E-02,-0.18621515E-01,-0.76492727E-01,
         0.46035284E-02,-0.83984109E-02,-0.38409315E-02,-0.24007121E-01,
         0.20295493E-02, 0.31944602E-02, 0.25548993E-02, 0.16937442E-01,
        -0.33483137E-02, 0.44926456E-02,-0.87231705E-02,-0.22460683E-04,
        -0.28610305E-03,-0.86711897E-02, 0.48961812E-02, 0.57416786E-02,
         0.46305876E-03,-0.21738780E-02, 0.76159358E-03,-0.36146194E-02,
         0.12506649E-02,-0.13366885E-02,-0.10683221E-01, 0.30920975E-03,
        -0.45289472E-02,-0.28718790E-03, 0.11802517E-02,-0.30694963E-03,
        -0.10180955E-02,-0.49620294E-02,-0.17384165E-02, 0.99514716E-03,
        -0.16873764E-04,-0.20088433E-03, 0.31644152E-03, 0.16569567E-02,
        -0.16204406E-03,-0.11483510E-02,-0.31498200E-02, 0.34733280E-03,
        -0.36707736E-03,-0.10854924E-02,-0.11638803E-02,-0.26081642E-03,
        -0.84228005E-03, 0.32775973E-02,-0.48009519E-03,-0.92913740E-03,
        -0.18467831E-02,-0.23035030E-02,-0.67342789E-03, 0.10568277E-03,
         0.82467523E-04, 0.18288986E-02,-0.13343882E-03,-0.24182136E-03,
        -0.47088927E-03, 0.37781309E-03, 0.17342493E-02,-0.30975725E-03,
         0.29352296E-03,-0.44152475E-03,-0.11773903E-02, 0.21130226E-02,
         0.17560000E-03,-0.71998278E-03,-0.16504697E-02,
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
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;
    float x54 = x53*x5;

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
        +coeff[ 10]        *x31    *x51
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]*x11        *x41
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]                *x54
        +coeff[ 16]*x11
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]    *x22        *x51
        +coeff[ 19]*x11*x21    *x41
        +coeff[ 20]    *x25
        +coeff[ 21]                *x52
        +coeff[ 22]*x11    *x31
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]*x11*x21*x31
        +coeff[ 25]            *x41*x52
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 26]    *x22    *x42
        +coeff[ 27]            *x45
        +coeff[ 28]    *x23
        +coeff[ 29]*x11            *x51
        +coeff[ 30]    *x21    *x41*x51
        +coeff[ 31]*x12
        +coeff[ 32]*x11*x22
        +coeff[ 33]    *x22*x31*x41
        +coeff[ 34]*x11        *x42
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 35]*x11*x21        *x51
        +coeff[ 36]        *x33*x41
        +coeff[ 37]*x11*x22*x31*x41
        +coeff[ 38]        *x32
        +coeff[ 39]        *x31*x41
        +coeff[ 40]        *x32*x41
        +coeff[ 41]        *x31*x42
        +coeff[ 42]            *x43
        +coeff[ 43]    *x21*x31    *x51
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 44]        *x31*x41*x51
        +coeff[ 45]    *x24
        +coeff[ 46]    *x23    *x41
        +coeff[ 47]        *x31    *x52
        +coeff[ 48]*x11    *x31*x41
        +coeff[ 49]    *x21    *x43
        +coeff[ 50]*x11        *x41*x51
        +coeff[ 51]    *x22    *x41*x51
        +coeff[ 52]*x11*x22    *x41
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 53]*x11*x21    *x42
        +coeff[ 54]*x11*x21    *x41*x51
        +coeff[ 55]*x11*x23*x31*x41
        +coeff[ 56]*x12*x21*x32*x41
        +coeff[ 57]    *x22    *x42*x53
        +coeff[ 58]        *x32    *x51
        +coeff[ 59]    *x21        *x52
        +coeff[ 60]    *x22*x32
        +coeff[ 61]    *x23        *x51
    ;
    v_p_l5p77_den                             =v_p_l5p77_den
        +coeff[ 62]    *x21*x31*x42
        +coeff[ 63]    *x22        *x52
        +coeff[ 64]    *x24*x31
        +coeff[ 65]    *x23*x32
        +coeff[ 66]*x11*x21*x31*x41
        +coeff[ 67]    *x23    *x42
        +coeff[ 68]*x12            *x51
        +coeff[ 69]*x12*x21    *x41
        +coeff[ 70]*x11*x23    *x41
        ;

    return v_p_l5p77_den                             ;
}
float l_l5p77_den                             (float *x,int m){
    int ncoeff= 65;
    float avdat=  0.1516336E-02;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 66]={
         0.39049501E-02, 0.19059914E+00, 0.11652781E-01,-0.28444843E-01,
        -0.82956357E-02, 0.71511075E-01, 0.10677863E-01, 0.16038816E-01,
        -0.98285321E-02, 0.39120195E-02, 0.25825746E-01,-0.92135658E-02,
         0.45197844E-02, 0.58856755E-02, 0.18460697E-01,-0.27388640E-01,
         0.13266541E-01, 0.11953975E-02, 0.10342068E-02, 0.24273649E-02,
        -0.10159287E-01,-0.19110739E-01,-0.29251581E-01,-0.17537119E-01,
        -0.26467766E-02, 0.16077575E-02,-0.32375283E-02,-0.48599737E-02,
         0.31946138E-02,-0.57673156E-02, 0.56751136E-03,-0.19925354E-04,
        -0.93934895E-03,-0.33215787E-02, 0.74828179E-02,-0.33831804E-02,
         0.18973226E-02,-0.21834600E-03, 0.16013077E-03,-0.19359533E-03,
        -0.29772814E-03,-0.86850673E-03,-0.17658640E-02, 0.70618623E-03,
        -0.58668503E-03,-0.73488191E-03,-0.20892795E-02, 0.28258751E-02,
         0.72583342E-02,-0.10150791E-02,-0.11514783E-02,-0.25894309E-02,
         0.49442961E-03, 0.13549415E-02,-0.27037449E-02,-0.97981712E-03,
        -0.26131229E-03, 0.21150585E-03, 0.11867697E-02, 0.20837728E-02,
         0.12818356E-02, 0.25375069E-02, 0.16346149E-02, 0.84082881E-03,
         0.11766455E-02,
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

//                 function

    float v_l_l5p77_den                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]            *x41
        +coeff[  3]*x11
        +coeff[  4]    *x22
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x21        *x51
        +coeff[  7]*x11        *x41
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[  8]    *x23*x31
        +coeff[  9]        *x31
        +coeff[ 10]    *x21*x31
        +coeff[ 11]            *x42
        +coeff[ 12]*x11*x21
        +coeff[ 13]*x11    *x31
        +coeff[ 14]    *x23
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x11*x22
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 17]    *x23*x31*x41
        +coeff[ 18]                *x51
        +coeff[ 19]*x11            *x51
        +coeff[ 20]    *x22    *x41
        +coeff[ 21]    *x21*x31*x41
        +coeff[ 22]    *x23    *x41
        +coeff[ 23]*x11*x22    *x41
        +coeff[ 24]        *x31*x41
        +coeff[ 25]            *x41*x51
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 26]    *x21*x32
        +coeff[ 27]*x11        *x42
        +coeff[ 28]*x12*x21
        +coeff[ 29]*x11*x22*x31
        +coeff[ 30]    *x24*x31
        +coeff[ 31]*x11    *x33*x41
        +coeff[ 32]    *x21        *x52
        +coeff[ 33]*x11    *x31*x41
        +coeff[ 34]    *x21    *x43
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 35]*x12*x21    *x41
        +coeff[ 36]    *x23*x31*x42
        +coeff[ 37]        *x32
        +coeff[ 38]        *x31    *x51
        +coeff[ 39]                *x52
        +coeff[ 40]*x12
        +coeff[ 41]    *x22        *x51
        +coeff[ 42]    *x21    *x41*x51
        +coeff[ 43]            *x42*x51
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 44]*x11    *x32
        +coeff[ 45]*x11*x21    *x41
        +coeff[ 46]    *x24
        +coeff[ 47]    *x22    *x42
        +coeff[ 48]    *x21*x31*x42
        +coeff[ 49]*x11*x23
        +coeff[ 50]*x12*x21*x31
        +coeff[ 51]    *x21    *x44
        +coeff[ 52]    *x23*x31    *x51
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 53]    *x23*x32*x41
        +coeff[ 54]    *x22*x31
        +coeff[ 55]    *x21*x31    *x51
        +coeff[ 56]*x11        *x41*x51
        +coeff[ 57]*x13
        +coeff[ 58]    *x22*x31*x41
        +coeff[ 59]    *x21*x32*x41
        +coeff[ 60]    *x21    *x42*x51
        +coeff[ 61]    *x24    *x41
    ;
    v_l_l5p77_den                             =v_l_l5p77_den
        +coeff[ 62]    *x23    *x41*x51
        +coeff[ 63]    *x23*x33
        +coeff[ 64]    *x23*x31*x41*x51
        ;

    return v_l_l5p77_den                             ;
}
float x_l5p77_dex                             (float *x,int m){
    int ncoeff= 60;
    float avdat= -0.1208238E-01;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
         0.14241206E-03,-0.52119441E-01, 0.13378075E+00, 0.26367402E+00,
         0.90340171E-02, 0.35096601E-01, 0.11536633E+00,-0.45938077E-03,
        -0.26952445E-02, 0.15541201E-03,-0.15303486E-01, 0.89938603E-02,
         0.25081202E-01, 0.12086685E-01, 0.42180032E-01, 0.31167908E-01,
        -0.45327324E-01, 0.35464656E-03,-0.24394237E-03, 0.32280586E-02,
        -0.35466556E-02,-0.65021445E-02, 0.39859372E-03, 0.49145087E-02,
         0.20371743E-01,-0.29620426E-01,-0.44937134E-01, 0.53749094E-03,
        -0.20446607E-02, 0.55211675E-02,-0.43904167E-02,-0.27490830E-01,
         0.91036293E-03, 0.11919651E-02,-0.72782510E-02, 0.12844797E-02,
         0.31675610E-02,-0.68389405E-02,-0.88722259E-02, 0.57099381E-03,
         0.12342742E-03, 0.14978982E-03,-0.30838489E-03,-0.83172839E-03,
        -0.86357101E-03,-0.52829981E-02, 0.19361933E-02, 0.17620833E-02,
         0.39077196E-02,-0.36506741E-02, 0.49953717E-02, 0.55011180E-02,
         0.60093934E-02, 0.66992780E-02, 0.23110566E-03, 0.51878579E-03,
         0.41569342E-03, 0.35173784E-03, 0.43052284E-03, 0.41108864E-03,
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

    float v_x_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]    *x21
        +coeff[  4]            *x41
        +coeff[  5]    *x21        *x51
        +coeff[  6]    *x21    *x41
        +coeff[  7]*x13        *x41
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[  8]*x12*x21*x31
        +coeff[  9]    *x21*x31    *x52
        +coeff[ 10]    *x23*x31
        +coeff[ 11]*x11    *x31
        +coeff[ 12]*x11        *x41
        +coeff[ 13]    *x22
        +coeff[ 14]    *x21*x31
        +coeff[ 15]    *x23
        +coeff[ 16]    *x21    *x42
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 17]*x13*x22
        +coeff[ 18]    *x23*x31*x41
        +coeff[ 19]        *x31
        +coeff[ 20]                *x52
        +coeff[ 21]            *x42
        +coeff[ 22]*x13
        +coeff[ 23]*x12*x21
        +coeff[ 24]*x11*x22
        +coeff[ 25]    *x21*x31*x41
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 26]    *x23    *x41
        +coeff[ 27]*x13*x22    *x41
        +coeff[ 28]        *x31*x41
        +coeff[ 29]    *x21    *x41*x51
        +coeff[ 30]    *x21*x32
        +coeff[ 31]*x11*x22    *x41
        +coeff[ 32]*x11            *x51
        +coeff[ 33]            *x41*x51
        +coeff[ 34]*x11        *x42
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 35]    *x21*x31    *x51
        +coeff[ 36]    *x22    *x41
        +coeff[ 37]*x12*x21    *x41
        +coeff[ 38]*x11*x22*x31
        +coeff[ 39]*x11    *x31*x42
        +coeff[ 40]*x13    *x31*x41
        +coeff[ 41]    *x22*x33
        +coeff[ 42]        *x32
        +coeff[ 43]    *x21        *x52
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 44]*x11    *x32
        +coeff[ 45]*x11    *x31*x41
        +coeff[ 46]    *x22*x31
        +coeff[ 47]    *x23        *x51
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]    *x23    *x42
        +coeff[ 50]    *x21*x31*x43
        +coeff[ 51]*x12*x21*x31*x42
        +coeff[ 52]*x12*x21    *x43
    ;
    v_x_l5p77_dex                             =v_x_l5p77_dex
        +coeff[ 53]    *x23*x31*x42
        +coeff[ 54]*x12
        +coeff[ 55]*x11*x21        *x51
        +coeff[ 56]*x11        *x41*x51
        +coeff[ 57]*x11*x21*x31
        +coeff[ 58]    *x22        *x51
        +coeff[ 59]        *x31*x41*x51
        ;

    return v_x_l5p77_dex                             ;
}
float t_l5p77_dex                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.5946257E+00;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.11097419E-02, 0.12160958E-01,-0.82930095E-01, 0.30563768E-01,
        -0.32364786E-01,-0.65177893E-02,-0.41248633E-04, 0.62629313E-03,
         0.39303997E-02,-0.10870630E-03,-0.27245556E-02,-0.24895023E-02,
        -0.11472783E-01,-0.68405881E-02,-0.13120668E-02,-0.86155273E-02,
         0.13207668E-01,-0.24932050E-03,-0.78266370E-03,-0.12403782E-02,
         0.68206713E-03,-0.73785573E-03,-0.55977842E-02, 0.73427721E-02,
         0.13132011E-01, 0.25352623E-03,-0.14103214E-02, 0.11523908E-02,
        -0.13302754E-02, 0.74498435E-02,-0.73186704E-03,-0.79468533E-04,
         0.42598558E-03,-0.12706447E-02, 0.19757794E-02, 0.51922578E-03,
        -0.68706099E-03, 0.46895124E-03, 0.57091174E-03, 0.25742729E-02,
         0.14578033E-02, 0.26253963E-03, 0.21959569E-03,-0.14918161E-04,
         0.22174476E-03,-0.20643621E-03,-0.36601067E-03, 0.11754013E-02,
         0.72663993E-03,-0.65153360E-03,
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
        +coeff[  3]                *x51
        +coeff[  4]    *x21    *x41
        +coeff[  5]    *x21        *x51
        +coeff[  6]        *x32*x41
        +coeff[  7]*x12*x21*x31
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[  8]    *x23*x31
        +coeff[  9]*x13        *x41
        +coeff[ 10]            *x41
        +coeff[ 11]*x11    *x31
        +coeff[ 12]    *x21*x31
        +coeff[ 13]*x11        *x41
        +coeff[ 14]                *x52
        +coeff[ 15]    *x23
        +coeff[ 16]    *x21    *x42
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 17]    *x23*x31*x41
        +coeff[ 18]        *x31
        +coeff[ 19]    *x22
        +coeff[ 20]        *x31*x41
        +coeff[ 21]*x11            *x51
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]    *x23    *x41
        +coeff[ 25]*x13*x22    *x41
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 26]*x12*x21
        +coeff[ 27]    *x21*x32
        +coeff[ 28]            *x42*x51
        +coeff[ 29]*x11*x22    *x41
        +coeff[ 30]            *x42
        +coeff[ 31]        *x31    *x51
        +coeff[ 32]            *x41*x51
        +coeff[ 33]    *x22    *x41
        +coeff[ 34]*x11        *x42
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 35]    *x22        *x51
        +coeff[ 36]    *x21    *x41*x51
        +coeff[ 37]    *x21        *x52
        +coeff[ 38]            *x41*x52
        +coeff[ 39]*x11*x22*x31
        +coeff[ 40]*x12*x21    *x41
        +coeff[ 41]*x11    *x32*x41
        +coeff[ 42]*x11    *x31*x42
        +coeff[ 43]*x11        *x43
    ;
    v_t_l5p77_dex                             =v_t_l5p77_dex
        +coeff[ 44]*x13    *x31*x41
        +coeff[ 45]*x11*x21
        +coeff[ 46]*x11*x21    *x41
        +coeff[ 47]*x11    *x31*x41
        +coeff[ 48]    *x22    *x42
        +coeff[ 49]    *x21*x31*x42
        ;

    return v_t_l5p77_dex                             ;
}
float y_l5p77_dex                             (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.2133712E-01;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.71358574E-02, 0.80412649E-01, 0.40040840E-02,-0.14572145E-01,
        -0.26727779E-01,-0.15582691E-01,-0.48446365E-01, 0.39306346E-01,
        -0.59689488E-02, 0.67246458E-02, 0.55606984E-01, 0.16914485E-01,
         0.76811993E-02, 0.77975090E-02,-0.56045121E-02,-0.46015870E-01,
        -0.97526964E-02,-0.11472446E-01, 0.11391870E-01, 0.95222128E-03,
        -0.64950937E-03,-0.80240145E-03,-0.10872807E-01,-0.18504147E-02,
        -0.44776672E-02, 0.20370274E-02,-0.67752395E-02,-0.70553301E-02,
        -0.36321115E-02,-0.42760181E-02, 0.31464726E-02,-0.11163497E-01,
         0.13659063E-02,-0.15730988E-01,-0.71900961E-03,-0.44113849E-02,
        -0.64346887E-03,-0.28438317E-02,-0.26280617E-02,-0.15650576E-02,
        -0.35569270E-03,-0.63038650E-02, 0.63005573E-03,-0.58419106E-03,
        -0.22088869E-02, 0.79174172E-02,-0.69655286E-03, 0.45771324E-02,
         0.76768105E-03,-0.38954762E-02,-0.46320367E-02,-0.20617889E-02,
        -0.86748056E-04,-0.14994604E-02,-0.42676300E-03, 0.39237812E-02,
         0.14326244E-02, 0.52811863E-03, 0.37122057E-02, 0.68271044E-03,
        -0.14898753E-02,-0.11587029E-02, 0.67991926E-02,-0.15452991E-02,
        -0.76127693E-03,-0.10031991E-02, 0.44821785E-02, 0.47106044E-02,
         0.17800258E-02,-0.85727661E-04,-0.35544143E-02, 0.59465109E-03,
        -0.29471819E-03,-0.27391504E-03, 0.51891580E-02, 0.28645532E-03,
         0.27824743E-03, 0.10672365E-02, 0.10475388E-02,-0.17123179E-02,
         0.54042245E-03, 0.82415348E-03,-0.73206425E-03, 0.51230157E-03,
        -0.12680602E-02, 0.27141976E-02,-0.46837874E-03, 0.18056488E-02,
         0.91006710E-04,-0.19393671E-03, 0.65764201E-04, 0.10962040E-02,
        -0.51744503E-03, 0.11276496E-02, 0.77525867E-04, 0.28655017E-03,
         0.53497741E-03,
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

    float v_y_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]            *x42
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x22
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[  8]    *x21*x31
        +coeff[  9]*x11        *x41
        +coeff[ 10]            *x41*x51
        +coeff[ 11]*x11*x21
        +coeff[ 12]    *x21        *x51
        +coeff[ 13]        *x31    *x51
        +coeff[ 14]                *x52
        +coeff[ 15]    *x22    *x41
        +coeff[ 16]    *x23
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 17]    *x22*x31
        +coeff[ 18]    *x22        *x51
        +coeff[ 19]        *x31*x43
        +coeff[ 20]    *x21    *x44
        +coeff[ 21]        *x31*x45
        +coeff[ 22]        *x31*x41
        +coeff[ 23]        *x32
        +coeff[ 24]            *x43
        +coeff[ 25]*x12
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 26]    *x21*x31*x41
        +coeff[ 27]*x11*x21    *x41
        +coeff[ 28]*x11*x21*x31
        +coeff[ 29]            *x41*x52
        +coeff[ 30]*x11*x21        *x51
        +coeff[ 31]    *x24
        +coeff[ 32]*x11
        +coeff[ 33]    *x21    *x42
        +coeff[ 34]*x11            *x51
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 35]*x11        *x42
        +coeff[ 36]            *x42*x51
        +coeff[ 37]    *x21    *x41*x51
        +coeff[ 38]*x11*x22
        +coeff[ 39]*x12        *x41
        +coeff[ 40]    *x22    *x42
        +coeff[ 41]*x11*x23
        +coeff[ 42]*x11    *x31
        +coeff[ 43]        *x31*x42
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 44]*x11    *x31*x41
        +coeff[ 45]    *x21    *x43
        +coeff[ 46]        *x31    *x52
        +coeff[ 47]    *x22*x31*x41
        +coeff[ 48]                *x53
        +coeff[ 49]*x11*x22    *x41
        +coeff[ 50]    *x22    *x41*x51
        +coeff[ 51]*x11*x21    *x41*x51
        +coeff[ 52]    *x23*x31*x42
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 53]        *x31*x41*x51
        +coeff[ 54]*x12    *x31
        +coeff[ 55]*x11*x21    *x42
        +coeff[ 56]    *x22*x32
        +coeff[ 57]*x12            *x51
        +coeff[ 58]*x11*x21*x31*x41
        +coeff[ 59]*x11*x21*x32
        +coeff[ 60]    *x22*x31    *x51
        +coeff[ 61]*x12*x21    *x41
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 62]    *x23    *x42
        +coeff[ 63]*x12*x22
        +coeff[ 64]*x11*x21*x31    *x51
        +coeff[ 65]    *x22        *x52
        +coeff[ 66]*x11*x22    *x42
        +coeff[ 67]    *x22    *x42*x51
        +coeff[ 68]    *x22*x31*x41*x51
        +coeff[ 69]    *x21    *x44*x51
        +coeff[ 70]    *x24    *x43
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 71]        *x32*x41
        +coeff[ 72]    *x21*x32
        +coeff[ 73]        *x32    *x51
        +coeff[ 74]    *x21*x31*x42
        +coeff[ 75]*x12*x21
        +coeff[ 76]    *x21        *x52
        +coeff[ 77]            *x43*x51
        +coeff[ 78]    *x21*x32*x41
        +coeff[ 79]    *x21    *x42*x51
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 80]        *x31*x42*x51
        +coeff[ 81]    *x23*x31
        +coeff[ 82]    *x21*x31*x41*x51
        +coeff[ 83]    *x21    *x41*x52
        +coeff[ 84]    *x24    *x41
        +coeff[ 85]    *x23*x31*x41
        +coeff[ 86]*x11*x21        *x52
        +coeff[ 87]*x11*x22*x31*x41
        +coeff[ 88]        *x33
    ;
    v_y_l5p77_dex                             =v_y_l5p77_dex
        +coeff[ 89]*x11    *x32
        +coeff[ 90]*x11    *x31    *x51
        +coeff[ 91]*x11        *x43
        +coeff[ 92]    *x23    *x41
        +coeff[ 93]*x11    *x31*x42
        +coeff[ 94]*x11            *x52
        +coeff[ 95]*x11    *x32*x41
        +coeff[ 96]        *x31*x44
        ;

    return v_y_l5p77_dex                             ;
}
float p_l5p77_dex                             (float *x,int m){
    int ncoeff= 90;
    float avdat=  0.2452515E-02;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.34324898E-03, 0.28541012E-03,-0.59078657E-02,-0.24187032E-02,
        -0.31464109E-02, 0.43447972E-02,-0.11827397E-02,-0.10927876E-01,
        -0.12504106E-02,-0.19306745E-02, 0.17285836E-02, 0.15238721E-02,
         0.12272272E-01, 0.20361969E-02,-0.23876601E-02,-0.15699255E-02,
        -0.14336584E-02, 0.15756374E-02,-0.81887292E-02,-0.34444332E-02,
         0.30228572E-02,-0.18289343E-02, 0.77796559E-03, 0.14995629E-03,
        -0.12437834E-02,-0.15867222E-03,-0.11213848E-02, 0.23575328E-03,
        -0.61483512E-03,-0.20277505E-02,-0.51860249E-03,-0.16437310E-02,
        -0.15220781E-02,-0.16917687E-04,-0.19973816E-03,-0.31585037E-03,
        -0.36230256E-03, 0.27714900E-03,-0.39086616E-03, 0.54241414E-03,
        -0.85733162E-03,-0.95265463E-03,-0.35946145E-04,-0.27083443E-03,
        -0.10499532E-02,-0.47420544E-03,-0.67495195E-04, 0.17735281E-03,
        -0.72353920E-04, 0.27408917E-03,-0.24565490E-04,-0.30835057E-03,
        -0.25649759E-03, 0.13929978E-03, 0.16211439E-02,-0.40683005E-03,
         0.14917992E-03,-0.51595661E-03,-0.14898699E-03, 0.25939321E-03,
         0.16896493E-02, 0.13834417E-03,-0.20164273E-03, 0.41509437E-03,
         0.11127443E-03,-0.20290181E-03, 0.34978584E-03, 0.17744131E-03,
         0.57411748E-04,-0.30673013E-03,-0.38689122E-04, 0.14216792E-03,
         0.90314797E-03,-0.16098580E-03, 0.20478347E-03,-0.48037455E-04,
         0.48248318E-03, 0.19314179E-03,-0.15116805E-03, 0.89211558E-03,
        -0.17934412E-03,-0.15201708E-03,-0.82555627E-04,-0.55181759E-03,
         0.24050359E-03, 0.54554612E-03, 0.14793289E-03,-0.32085128E-03,
        -0.52882289E-03, 0.41001619E-03,
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
    float x25 = x24*x2;
    float x26 = x25*x2;
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
        +coeff[ 14]    *x23
        +coeff[ 15]                *x52
        +coeff[ 16]    *x22*x31
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 17]*x11        *x41
        +coeff[ 18]    *x22    *x41
        +coeff[ 19]    *x21    *x42
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]*x11*x21        *x51
        +coeff[ 23]*x11
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]*x11            *x51
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 26]            *x43
        +coeff[ 27]*x12
        +coeff[ 28]*x11*x22
        +coeff[ 29]    *x24
        +coeff[ 30]            *x41*x52
        +coeff[ 31]    *x22    *x42
        +coeff[ 32]    *x22    *x41*x51
        +coeff[ 33]*x13*x21*x31
        +coeff[ 34]        *x32
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 35]        *x31*x42
        +coeff[ 36]        *x31*x41*x51
        +coeff[ 37]    *x21        *x52
        +coeff[ 38]*x11*x21*x31
        +coeff[ 39]    *x23    *x41
        +coeff[ 40]*x11        *x42
        +coeff[ 41]*x11*x23
        +coeff[ 42]*x11*x22*x31
        +coeff[ 43]*x12        *x41
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 44]*x11*x22    *x41
        +coeff[ 45]*x11*x21    *x41*x51
        +coeff[ 46]*x13    *x31*x41
        +coeff[ 47]*x11    *x31
        +coeff[ 48]        *x32    *x51
        +coeff[ 49]    *x23*x31
        +coeff[ 50]*x11*x21    *x41
        +coeff[ 51]*x11    *x31*x41
        +coeff[ 52]    *x23        *x51
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 53]                *x53
        +coeff[ 54]    *x21    *x43
        +coeff[ 55]    *x22*x31    *x51
        +coeff[ 56]*x11        *x41*x51
        +coeff[ 57]    *x21    *x42*x51
        +coeff[ 58]    *x21    *x41*x52
        +coeff[ 59]*x11*x21*x31*x41
        +coeff[ 60]    *x23    *x42
        +coeff[ 61]*x12            *x51
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 62]*x11*x21*x31    *x51
        +coeff[ 63]    *x23    *x41*x51
        +coeff[ 64]    *x21*x33    *x51
        +coeff[ 65]*x12*x21    *x41
        +coeff[ 66]    *x23*x32*x41
        +coeff[ 67]    *x23*x31*x42
        +coeff[ 68]    *x21*x33*x42
        +coeff[ 69]            *x42*x51
        +coeff[ 70]        *x31    *x52
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 71]    *x22*x32
        +coeff[ 72]    *x21*x31*x42
        +coeff[ 73]    *x21*x31*x41*x51
        +coeff[ 74]            *x43*x51
        +coeff[ 75]*x12    *x31
        +coeff[ 76]    *x23*x31*x41
        +coeff[ 77]*x11*x21    *x42
        +coeff[ 78]    *x24        *x51
        +coeff[ 79]    *x22    *x42*x51
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 80]*x12*x22
        +coeff[ 81]*x11*x24
        +coeff[ 82]*x11*x21        *x52
        +coeff[ 83]*x11*x23    *x41
        +coeff[ 84]    *x22    *x41*x52
        +coeff[ 85]*x11*x22    *x42
        +coeff[ 86]    *x24    *x42
        +coeff[ 87]    *x26    *x41
        +coeff[ 88]    *x24    *x43
    ;
    v_p_l5p77_dex                             =v_p_l5p77_dex
        +coeff[ 89]    *x24*x31*x41*x51
        ;

    return v_p_l5p77_dex                             ;
}
float l_l5p77_dex                             (float *x,int m){
    int ncoeff= 67;
    float avdat= -0.1868164E-01;
    float xmin[10]={
        -0.14896E-01,-0.45454E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 68]={
         0.61223824E-02,-0.37794760E+00,-0.65720002E-02,-0.10009057E+00,
         0.67633815E-01,-0.26280245E-01,-0.15742211E+00,-0.35015296E-01,
        -0.33916529E-01, 0.18747184E-01,-0.56147948E-01,-0.12430671E-01,
        -0.43263357E-01, 0.61845157E-01,-0.27982956E-01,-0.26350019E-02,
        -0.24675599E-02,-0.31268864E-02,-0.15351003E-01, 0.41321203E-01,
        -0.69689173E-02, 0.60933065E-01, 0.37192378E-01, 0.18186297E-02,
         0.37177682E-02,-0.53761159E-02, 0.70077418E-02, 0.97537665E-02,
         0.12185124E-01,-0.21382559E-03, 0.21483355E-03,-0.26590661E-02,
         0.95835677E-03,-0.19376655E-02, 0.11515807E-02, 0.71106218E-02,
         0.73826932E-02, 0.28632677E-03,-0.81877160E-03,-0.97611680E-03,
         0.14116832E-02,-0.18593052E-02,-0.13083748E-01, 0.24987822E-02,
         0.62099652E-03, 0.40450973E-05, 0.34414945E-03,-0.75762481E-02,
        -0.51654694E-02,-0.17313514E-03,-0.21941699E-03, 0.33539446E-03,
         0.13550968E-02,-0.64169063E-03,-0.54213259E-03,-0.34153853E-02,
         0.32367476E-02,-0.12176144E-01,-0.12934529E-02, 0.12554483E-02,
        -0.25529719E-02, 0.46500412E-03, 0.58236970E-02,-0.56725874E-03,
         0.29857152E-02,-0.49231332E-02,-0.45345575E-02,
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

//                 function

    float v_l_l5p77_dex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]*x11
        +coeff[  5]    *x22
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x21        *x51
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x23*x31
        +coeff[ 10]    *x21*x31
        +coeff[ 11]*x11    *x31
        +coeff[ 12]    *x23
        +coeff[ 13]    *x21    *x42
        +coeff[ 14]*x11*x22
        +coeff[ 15]    *x23*x31*x41
        +coeff[ 16]        *x31
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 17]*x11            *x51
        +coeff[ 18]    *x22    *x41
        +coeff[ 19]    *x21*x31*x41
        +coeff[ 20]*x12*x21
        +coeff[ 21]    *x23    *x41
        +coeff[ 22]*x11*x22    *x41
        +coeff[ 23]                *x52
        +coeff[ 24]*x11*x21
        +coeff[ 25]    *x22*x31
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 26]    *x21*x32
        +coeff[ 27]*x11        *x42
        +coeff[ 28]*x11*x22*x31
        +coeff[ 29]*x11    *x33*x41
        +coeff[ 30]        *x31*x41
        +coeff[ 31]            *x42
        +coeff[ 32]            *x41*x51
        +coeff[ 33]    *x21    *x41*x51
        +coeff[ 34]*x11    *x32
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]*x12*x21    *x41
        +coeff[ 37]        *x31    *x51
        +coeff[ 38]*x12
        +coeff[ 39]    *x22        *x51
        +coeff[ 40]    *x21        *x52
        +coeff[ 41]*x11*x21    *x41
        +coeff[ 42]    *x21    *x43
        +coeff[ 43]*x12*x21*x31
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 44]    *x21*x31*x43
        +coeff[ 45]            *x42*x53
        +coeff[ 46]    *x23*x31*x42
        +coeff[ 47]    *x21*x33*x42
        +coeff[ 48]    *x21*x32*x43
        +coeff[ 49]    *x21*x31*x44
        +coeff[ 50]    *x23*x34*x41
        +coeff[ 51]        *x32
        +coeff[ 52]            *x42*x51
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 53]*x11*x21*x31
        +coeff[ 54]*x13
        +coeff[ 55]    *x24
        +coeff[ 56]    *x22    *x42
        +coeff[ 57]    *x21*x31*x42
        +coeff[ 58]    *x23        *x51
        +coeff[ 59]    *x22    *x41*x51
        +coeff[ 60]*x11*x23
        +coeff[ 61]    *x22*x32*x41
    ;
    v_l_l5p77_dex                             =v_l_l5p77_dex
        +coeff[ 62]    *x21    *x44
        +coeff[ 63]*x11*x23        *x51
        +coeff[ 64]    *x24*x31*x41
        +coeff[ 65]    *x21*x34*x41
        +coeff[ 66]    *x23    *x43*x51
        ;

    return v_l_l5p77_dex                             ;
}
float x_l5p77_fp                              (float *x,int m){
    int ncoeff= 49;
    float avdat=  0.4850814E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 50]={
         0.73665543E-02,-0.37995819E-01, 0.65293777E+00,-0.52604344E-01,
         0.43976791E-01, 0.15137231E-01, 0.36930632E-01,-0.30657919E-02,
         0.77053253E-02, 0.14181478E-01,-0.14557523E-01,-0.29401777E-01,
         0.39109876E-02, 0.56102672E-02, 0.17066566E-01, 0.10962527E-01,
        -0.27189681E-02, 0.27488384E-02, 0.41870973E-02,-0.44432557E-02,
        -0.79075755E-04,-0.32041839E-02, 0.55091535E-02,-0.10267943E-01,
         0.39130161E-02,-0.20187652E-01,-0.13401572E-02, 0.16599483E-03,
         0.11631739E-02, 0.14370932E-02, 0.15681073E-02, 0.20424565E-02,
        -0.14076697E-02,-0.10330845E-01,-0.49683470E-02,-0.43462724E-02,
         0.63118450E-02, 0.64997573E-03,-0.14600570E-02, 0.39820308E-02,
         0.35143178E-02,-0.11653090E-02,-0.25365625E-02, 0.37741067E-02,
        -0.25085355E-02, 0.11700870E-02,-0.48278663E-02, 0.34487743E-02,
         0.16507468E-02,
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
        +coeff[  7]    *x21
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x21*x31
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]            *x41
        +coeff[ 13]                *x53
        +coeff[ 14]    *x21    *x41*x51
        +coeff[ 15]    *x23
        +coeff[ 16]*x11            *x51
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 17]*x11    *x31
        +coeff[ 18]            *x41*x51
        +coeff[ 19]        *x31*x41
        +coeff[ 20]*x12    *x31
        +coeff[ 21]    *x21        *x52
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]*x11*x23
        +coeff[ 25]    *x23    *x41
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 26]    *x23*x31    *x51
        +coeff[ 27]*x13*x22    *x41
        +coeff[ 28]        *x31
        +coeff[ 29]        *x31    *x51
        +coeff[ 30]*x12*x21
        +coeff[ 31]    *x22        *x51
        +coeff[ 32]    *x21*x32
        +coeff[ 33]*x11*x22    *x41
        +coeff[ 34]    *x21    *x42*x51
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 35]    *x23*x31
        +coeff[ 36]    *x22    *x42
        +coeff[ 37]*x11            *x52
        +coeff[ 38]*x11    *x31*x41
        +coeff[ 39]    *x21*x31    *x51
        +coeff[ 40]*x11*x22        *x51
        +coeff[ 41]    *x22        *x52
        +coeff[ 42]*x11*x22*x31
        +coeff[ 43]    *x23        *x51
    ;
    v_x_l5p77_fp                              =v_x_l5p77_fp
        +coeff[ 44]    *x21*x31*x41*x51
        +coeff[ 45]*x13        *x41*x51
        +coeff[ 46]    *x23    *x41*x51
        +coeff[ 47]    *x22    *x42*x51
        +coeff[ 48]*x13    *x31    *x53
        ;

    return v_x_l5p77_fp                              ;
}
float t_l5p77_fp                              (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1246783E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.13610454E-02,-0.34647479E-02,-0.18788554E-01,-0.54564876E-04,
         0.58289759E-04, 0.11297005E+00, 0.53696432E-02,-0.10615933E-01,
         0.11454147E-02,-0.11069913E-02,-0.22524882E-02,-0.22141344E-02,
         0.21624689E-02,-0.78605616E-03, 0.84544026E-03,-0.69176994E-03,
        -0.31655843E-03, 0.63878635E-03, 0.11204408E-02, 0.25421859E-04,
        -0.17859080E-03,-0.27667859E-03,-0.36402428E-03,-0.73233357E-03,
         0.40273461E-03, 0.40990915E-03,-0.63067366E-03, 0.13648443E-02,
        -0.78024612E-04, 0.24767272E-03,-0.28033348E-03, 0.22862288E-03,
         0.19827017E-03, 0.16813961E-03, 0.41615704E-03,-0.68232138E-03,
         0.71254856E-03, 0.67612348E-03,-0.87276765E-03,-0.64878201E-04,
         0.99597441E-04,-0.60150411E-04, 0.10238977E-03, 0.89657100E-04,
         0.11578119E-03,-0.11452480E-03, 0.16251726E-03, 0.14832940E-03,
        -0.24610086E-03,-0.11560744E-03,
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
        +coeff[  9]    *x21    *x41
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]    *x21    *x41*x51
        +coeff[ 13]    *x21        *x52
        +coeff[ 14]*x11*x23
        +coeff[ 15]        *x31*x41
        +coeff[ 16]*x11            *x51
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 17]            *x41*x51
        +coeff[ 18]                *x53
        +coeff[ 19]        *x31    *x53
        +coeff[ 20]    *x21*x31
        +coeff[ 21]*x11        *x41
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x22    *x41
        +coeff[ 24]*x11        *x42
        +coeff[ 25]    *x21*x31    *x51
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 26]    *x23    *x41
        +coeff[ 27]    *x22    *x42
        +coeff[ 28]*x11    *x31
        +coeff[ 29]        *x31    *x51
        +coeff[ 30]    *x22*x31
        +coeff[ 31]    *x22        *x51
        +coeff[ 32]            *x42*x51
        +coeff[ 33]*x11            *x52
        +coeff[ 34]    *x23        *x51
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 35]    *x21    *x42*x51
        +coeff[ 36]    *x22    *x43
        +coeff[ 37]*x13*x22        *x51
        +coeff[ 38]    *x23*x31*x41*x51
        +coeff[ 39]*x12
        +coeff[ 40]*x11*x21
        +coeff[ 41]*x12*x21
        +coeff[ 42]        *x31*x42
        +coeff[ 43]*x11    *x31    *x51
    ;
    v_t_l5p77_fp                              =v_t_l5p77_fp
        +coeff[ 44]*x11        *x41*x51
        +coeff[ 45]            *x41*x52
        +coeff[ 46]*x12*x21    *x41
        +coeff[ 47]*x11*x21*x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x13            *x51
        ;

    return v_t_l5p77_fp                              ;
}
float y_l5p77_fp                              (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.1005435E-01;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
         0.65170373E-02,-0.59205592E-01, 0.13075430E-01, 0.11103833E-01,
        -0.17059630E-01,-0.15923406E-02,-0.32548006E-02,-0.71285306E-02,
        -0.17285218E-02,-0.46245498E-02,-0.33543792E-02, 0.10840210E-01,
        -0.68390802E-02,-0.68930436E-04,-0.64344858E-04, 0.69019529E-02,
        -0.71534398E-03, 0.52608508E-02, 0.24673267E-03, 0.13028190E-02,
        -0.25707176E-02, 0.64810598E-02, 0.45800009E-02, 0.78306231E-03,
         0.26954676E-02,-0.88454189E-03, 0.63687604E-03,-0.19460334E-02,
         0.45873541E-02, 0.16111672E-02, 0.27144498E-02, 0.78555220E-03,
        -0.43049683E-02, 0.10812604E-02, 0.29410238E-02,-0.28960968E-02,
        -0.16620336E-02,-0.56533527E-03, 0.20937492E-03,-0.46297577E-02,
         0.78232371E-03,-0.18276208E-02, 0.53455081E-03, 0.83134189E-03,
         0.28245812E-02,-0.11332378E-02, 0.54305681E-03,-0.21958493E-02,
        -0.14797725E-02, 0.26076841E-02,-0.16080920E-02,-0.62764520E-02,
         0.22293869E-02, 0.12296953E-02, 0.63777401E-03, 0.17202366E-03,
         0.43299180E-03,-0.19844405E-02, 0.44707116E-03,-0.11669641E-02,
        -0.30878978E-03,-0.15550818E-02,-0.43529461E-03, 0.61408174E-03,
        -0.11350361E-02, 0.60451211E-03, 0.12153829E-02,-0.22525168E-02,
        -0.67765207E-03,-0.84672478E-03, 0.23823092E-02, 0.26475871E-03,
         0.10984284E-02,-0.94178133E-04,-0.29518903E-03,-0.34515437E-03,
        -0.14498901E-03,-0.12038877E-03, 0.20628994E-03, 0.99936570E-03,
        -0.41741613E-03,-0.58441347E-03, 0.17096348E-03, 0.53942105E-03,
        -0.35443399E-03,-0.33422135E-03, 0.13789161E-02, 0.29712886E-03,
        -0.52222196E-03,-0.26421202E-03, 0.39295095E-03, 0.40099127E-03,
        -0.27005718E-03, 0.43283586E-03, 0.75291598E-03,-0.13170148E-03,
        -0.59410941E-03,
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

    float v_y_l5p77_fp                              =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]                *x51
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x22
        +coeff[  5]*x11        *x41
        +coeff[  6]            *x41*x51
        +coeff[  7]*x11*x21
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[  8]    *x21        *x51
        +coeff[  9]        *x31    *x51
        +coeff[ 10]                *x52
        +coeff[ 11]    *x22    *x41
        +coeff[ 12]    *x21    *x41*x51
        +coeff[ 13]            *x44
        +coeff[ 14]        *x31*x43
        +coeff[ 15]            *x41*x52
        +coeff[ 16]*x11
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 17]*x11*x21    *x41
        +coeff[ 18]    *x22*x31    *x52
        +coeff[ 19]    *x22    *x41*x53
        +coeff[ 20]        *x31
        +coeff[ 21]            *x42
        +coeff[ 22]        *x31*x41
        +coeff[ 23]        *x32
        +coeff[ 24]    *x21    *x42
        +coeff[ 25]*x12
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 26]*x11            *x51
        +coeff[ 27]            *x42*x51
        +coeff[ 28]    *x22*x31
        +coeff[ 29]*x11*x21*x31
        +coeff[ 30]    *x22        *x51
        +coeff[ 31]*x11        *x41*x51
        +coeff[ 32]    *x22    *x42
        +coeff[ 33]*x11*x21        *x51
        +coeff[ 34]        *x31    *x52
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 35]    *x22*x31*x41
        +coeff[ 36]    *x23        *x51
        +coeff[ 37]    *x23    *x42
        +coeff[ 38]    *x21*x31    *x52
        +coeff[ 39]            *x41*x53
        +coeff[ 40]    *x23        *x52
        +coeff[ 41]    *x21
        +coeff[ 42]            *x43
        +coeff[ 43]*x11        *x42
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 44]    *x23
        +coeff[ 45]        *x31*x41*x51
        +coeff[ 46]*x11*x22
        +coeff[ 47]*x11*x21    *x42
        +coeff[ 48]    *x21    *x42*x51
        +coeff[ 49]    *x24
        +coeff[ 50]*x11*x21*x31*x41
        +coeff[ 51]    *x22    *x41*x51
        +coeff[ 52]*x11*x23
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 53]    *x21    *x41*x52
        +coeff[ 54]    *x21*x31
        +coeff[ 55]*x11    *x31
        +coeff[ 56]    *x21*x31*x41
        +coeff[ 57]    *x21    *x43
        +coeff[ 58]*x12        *x41
        +coeff[ 59]    *x21*x31*x42
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]            *x43*x51
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 62]*x11            *x52
        +coeff[ 63]                *x53
        +coeff[ 64]    *x22*x31    *x51
        +coeff[ 65]*x12*x22
        +coeff[ 66]    *x21    *x43*x51
        +coeff[ 67]    *x22    *x42*x51
        +coeff[ 68]    *x21        *x53
        +coeff[ 69]        *x31    *x53
        +coeff[ 70]    *x22    *x41*x52
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 71]*x11    *x31*x41*x52
        +coeff[ 72]            *x44*x52
        +coeff[ 73]        *x31*x42
        +coeff[ 74]        *x32*x41
        +coeff[ 75]    *x21*x31    *x51
        +coeff[ 76]        *x32    *x51
        +coeff[ 77]*x12*x21
        +coeff[ 78]*x12    *x31
        +coeff[ 79]    *x21        *x52
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 80]    *x23*x31
        +coeff[ 81]    *x22*x32
        +coeff[ 82]*x12            *x51
        +coeff[ 83]*x11*x22    *x41
        +coeff[ 84]*x11        *x42*x51
        +coeff[ 85]*x11*x21*x32
        +coeff[ 86]    *x22    *x43
        +coeff[ 87]*x12*x21    *x41
        +coeff[ 88]*x11*x21    *x41*x51
    ;
    v_y_l5p77_fp                              =v_y_l5p77_fp
        +coeff[ 89]*x11*x21*x31    *x51
        +coeff[ 90]    *x22        *x52
        +coeff[ 91]*x11    *x31*x43
        +coeff[ 92]*x12        *x41*x51
        +coeff[ 93]*x11        *x41*x52
        +coeff[ 94]    *x21*x31*x42*x51
        +coeff[ 95]    *x21*x34
        +coeff[ 96]*x11*x23    *x41
        ;

    return v_y_l5p77_fp                              ;
}
float p_l5p77_fp                              (float *x,int m){
    int ncoeff= 88;
    float avdat= -0.8726283E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 89]={
         0.31087459E-02,-0.17550972E-02, 0.60103172E-02,-0.32550465E-01,
         0.11132461E-01,-0.15728759E-01, 0.21492946E-02, 0.20090608E-01,
         0.65019373E-02,-0.32881526E-02,-0.36899175E-02,-0.20557264E-01,
        -0.68714907E-02, 0.13136443E-02,-0.27759343E-02, 0.17848238E-01,
         0.55295043E-02,-0.36096491E-02,-0.18700310E-02, 0.16883461E-03,
        -0.16874974E-03,-0.44804325E-03,-0.66132176E-04,-0.14528082E-03,
         0.78978599E-03, 0.44170893E-02, 0.48174555E-02,-0.86248352E-03,
         0.43175304E-02, 0.30869315E-02,-0.49025076E-03,-0.63911907E-03,
         0.43229964E-02, 0.19765997E-02, 0.36523244E-03,-0.39645040E-03,
         0.10976845E-02, 0.16878002E-02, 0.94880757E-03, 0.29695451E-02,
         0.14757375E-02,-0.88586583E-03, 0.25625774E-02, 0.12724109E-02,
        -0.12439745E-02, 0.23487753E-03, 0.12303308E-02, 0.66318561E-03,
        -0.25460871E-02, 0.21825028E-05, 0.62858849E-03, 0.13259274E-02,
         0.68047480E-03, 0.51072915E-03,-0.11788579E-02,-0.34821380E-04,
        -0.20769864E-03,-0.29133211E-03, 0.14745751E-02,-0.61941094E-03,
        -0.13578311E-03, 0.17639925E-03,-0.14418671E-02,-0.14228341E-02,
        -0.14806696E-03, 0.49855199E-03,-0.25576004E-03,-0.16686295E-02,
        -0.21363680E-03,-0.34475076E-03, 0.68151834E-03, 0.85092586E-04,
        -0.10972444E-02,-0.11117007E-03,-0.86817985E-04, 0.56695804E-04,
         0.46252844E-04,-0.11806986E-02, 0.13195942E-03,-0.11505177E-03,
         0.35964517E-03,-0.46569301E-03,-0.28071320E-03,-0.12628649E-02,
        -0.22171719E-03, 0.55413943E-03, 0.17941769E-03,-0.43722929E-03,
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
    float x25 = x24*x2;
    float x26 = x25*x2;
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
        +coeff[  8]            *x42
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
        +coeff[ 18]    *x22*x31*x41
        +coeff[ 19]        *x33*x41
        +coeff[ 20]    *x25
        +coeff[ 21]    *x24*x31
        +coeff[ 22]            *x43*x52
        +coeff[ 23]        *x35*x41
        +coeff[ 24]        *x32
        +coeff[ 25]        *x31*x41
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 26]    *x22*x31
        +coeff[ 27]*x12
        +coeff[ 28]    *x24
        +coeff[ 29]*x11*x21    *x41
        +coeff[ 30]*x11*x23*x31
        +coeff[ 31]*x11
        +coeff[ 32]    *x23
        +coeff[ 33]    *x21*x31*x41
        +coeff[ 34]*x11            *x51
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 35]            *x42*x51
        +coeff[ 36]*x11*x22
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]        *x31    *x52
        +coeff[ 39]            *x41*x52
        +coeff[ 40]*x11        *x42
        +coeff[ 41]*x11*x21        *x51
        +coeff[ 42]*x11*x23
        +coeff[ 43]    *x22    *x43
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 44]            *x41*x53
        +coeff[ 45]            *x45
        +coeff[ 46]*x11*x22    *x43
        +coeff[ 47]*x11    *x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x11*x22*x31
        +coeff[ 50]*x12        *x41
        +coeff[ 51]*x11*x22    *x41
        +coeff[ 52]*x12*x22
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 53]*x12*x21    *x41
        +coeff[ 54]    *x23*x31*x42
        +coeff[ 55]*x13    *x31
        +coeff[ 56]*x11    *x31
        +coeff[ 57]        *x32*x41
        +coeff[ 58]            *x43
        +coeff[ 59]    *x22*x32
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]*x12    *x31
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 62]*x11*x21*x31*x41
        +coeff[ 63]*x11*x21    *x42
        +coeff[ 64]*x12            *x51
        +coeff[ 65]*x11*x21    *x41*x51
        +coeff[ 66]        *x31    *x53
        +coeff[ 67]    *x22    *x42*x51
        +coeff[ 68]    *x25        *x51
        +coeff[ 69]*x11*x22        *x52
        +coeff[ 70]    *x26    *x41
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 71]        *x31*x44*x51
        +coeff[ 72]    *x25*x32*x41
        +coeff[ 73]    *x21*x31    *x51
        +coeff[ 74]    *x21    *x41*x51
        +coeff[ 75]        *x32    *x51
        +coeff[ 76]*x11    *x32
        +coeff[ 77]    *x21*x31*x42
        +coeff[ 78]*x11        *x41*x51
        +coeff[ 79]*x12*x21
    ;
    v_p_l5p77_fp                              =v_p_l5p77_fp
        +coeff[ 80]    *x22        *x52
        +coeff[ 81]            *x43*x51
        +coeff[ 82]*x11*x21*x32
        +coeff[ 83]    *x23    *x42
        +coeff[ 84]*x11        *x43
        +coeff[ 85]    *x21    *x44
        +coeff[ 86]    *x23        *x52
        +coeff[ 87]    *x25*x31
        ;

    return v_p_l5p77_fp                              ;
}
float l_l5p77_fp                              (float *x,int m){
    int ncoeff= 81;
    float avdat= -0.1987057E-01;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 82]={
         0.83441073E-02,-0.24360231E+00,-0.34637661E-02,-0.33813886E-01,
         0.41575912E-01,-0.26284715E-01,-0.99122055E-01,-0.26512446E-01,
        -0.21665953E-01, 0.11416291E-01,-0.35530519E-01, 0.39644921E-02,
        -0.79143019E-02,-0.27368957E-01,-0.20160317E-01, 0.39728206E-01,
        -0.17377678E-01, 0.40605819E-03,-0.30036096E-03,-0.37966294E-02,
        -0.79438174E-02,-0.34970627E-02, 0.25771398E-02,-0.16564309E-02,
         0.25794279E-01,-0.23089738E-02, 0.49719093E-02,-0.44403509E-02,
         0.38621280E-01, 0.23833012E-01,-0.92207291E-03,-0.56213476E-02,
         0.42304178E-02, 0.59777969E-02, 0.77205044E-02,-0.77581900E-03,
        -0.24938094E-03,-0.82184188E-03, 0.92981802E-03,-0.27437564E-02,
         0.47346670E-02,-0.41362992E-02, 0.45841848E-02,-0.21131453E-03,
         0.13493160E-02, 0.14714708E-02, 0.71005872E-03, 0.93072787E-03,
         0.67209657E-02,-0.10662257E-01,-0.94233248E-02,-0.74733078E-03,
        -0.27614215E-02, 0.15268204E-02, 0.71182516E-02, 0.50195251E-02,
        -0.13378052E-03, 0.73608979E-04, 0.12467285E-02,-0.45532800E-03,
         0.42610307E-03, 0.24744534E-03, 0.33570733E-03, 0.10494553E-03,
        -0.32920748E-03, 0.30218664E-03, 0.31781604E-02,-0.33803007E-02,
        -0.48886548E-03, 0.10078787E-02,-0.10469309E-02,-0.80185715E-03,
         0.51515625E-03,-0.81025198E-03,-0.52791811E-03,-0.73398004E-03,
         0.18900851E-02, 0.66593068E-03, 0.25862025E-02, 0.96391345E-03,
         0.41168286E-02,
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
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x23*x31
        +coeff[ 10]    *x21*x31
        +coeff[ 11]*x11*x21
        +coeff[ 12]*x11    *x31
        +coeff[ 13]    *x23
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x11*x22
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 17]            *x44
        +coeff[ 18]    *x22*x33
        +coeff[ 19]    *x23*x31*x41
        +coeff[ 20]            *x42
        +coeff[ 21]    *x21        *x51
        +coeff[ 22]            *x41*x51
        +coeff[ 23]*x11            *x51
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]    *x21        *x52
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 26]                *x53
        +coeff[ 27]*x12*x21
        +coeff[ 28]    *x23    *x41
        +coeff[ 29]*x11*x22    *x41
        +coeff[ 30]        *x31
        +coeff[ 31]    *x22*x31
        +coeff[ 32]    *x21*x32
        +coeff[ 33]*x11        *x42
        +coeff[ 34]*x11*x22*x31
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 35]*x11*x22*x31*x41
        +coeff[ 36]*x11    *x33*x41
        +coeff[ 37]*x12
        +coeff[ 38]        *x31*x41*x51
        +coeff[ 39]*x11*x21    *x41
        +coeff[ 40]*x11    *x31*x41
        +coeff[ 41]    *x24
        +coeff[ 42]*x12*x21    *x41
        +coeff[ 43]        *x31*x41
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 44]    *x21*x31    *x51
        +coeff[ 45]    *x21    *x41*x51
        +coeff[ 46]*x11    *x32
        +coeff[ 47]*x11        *x41*x51
        +coeff[ 48]    *x22    *x42
        +coeff[ 49]    *x21*x31*x42
        +coeff[ 50]    *x21    *x43
        +coeff[ 51]    *x21*x31*x41*x51
        +coeff[ 52]*x11*x23
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 53]*x12*x21*x31
        +coeff[ 54]    *x24    *x41
        +coeff[ 55]    *x21    *x44
        +coeff[ 56]    *x24*x31*x41
        +coeff[ 57]    *x23*x32*x41
        +coeff[ 58]            *x43
        +coeff[ 59]    *x22        *x51
        +coeff[ 60]            *x41*x52
        +coeff[ 61]*x11    *x31    *x51
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 62]*x11            *x52
        +coeff[ 63]*x12    *x31
        +coeff[ 64]*x13
        +coeff[ 65]    *x22*x32
        +coeff[ 66]    *x22*x31*x41
        +coeff[ 67]    *x21*x32*x41
        +coeff[ 68]    *x23        *x51
        +coeff[ 69]    *x21    *x42*x51
        +coeff[ 70]    *x21    *x41*x52
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 71]            *x42*x52
        +coeff[ 72]    *x21        *x53
        +coeff[ 73]                *x54
        +coeff[ 74]*x11*x21    *x41*x51
        +coeff[ 75]*x11        *x42*x51
        +coeff[ 76]    *x24*x31
        +coeff[ 77]        *x33*x42
        +coeff[ 78]*x11*x23    *x41
        +coeff[ 79]*x13*x21    *x42
    ;
    v_l_l5p77_fp                              =v_l_l5p77_fp
        +coeff[ 80]    *x23*x31*x43
        ;

    return v_l_l5p77_fp                              ;
}
float x_l5p77_q1en                            (float *x,int m){
    int ncoeff= 42;
    float avdat=  0.7227922E-03;
    float xmin[10]={
        -0.14934E-01,-0.45642E-01,-0.14989E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26979E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 43]={
         0.80256379E-03, 0.12972687E-01, 0.83134888E-04, 0.98311067E-01,
         0.82191620E-02, 0.29101090E-02, 0.67206105E-03, 0.20416840E-02,
         0.24221349E-03, 0.73011231E-03, 0.25607753E-02,-0.31309172E-02,
         0.41625994E-04,-0.15509959E-04,-0.66711502E-04, 0.31749462E-03,
         0.19516352E-04, 0.17708282E-02,-0.26377642E-02, 0.17661961E-03,
         0.24789502E-03, 0.15327513E-03, 0.20961775E-03, 0.41888733E-03,
        -0.63985621E-03,-0.29388021E-03,-0.40021687E-03,-0.22324207E-02,
        -0.17940442E-02,-0.79171848E-03, 0.96587371E-03, 0.39200537E-03,
         0.59295464E-04, 0.94379728E-04,-0.10770930E-03,-0.43872980E-03,
        -0.97679258E-04,-0.33770688E-03,-0.54267747E-03,-0.98956865E-04,
         0.77067304E-03, 0.64745627E-03,
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

    float v_x_l5p77_q1en                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]    *x21
        +coeff[  4]    *x21    *x41
        +coeff[  5]    *x21*x31
        +coeff[  6]            *x41
        +coeff[  7]*x11        *x41
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[  8]        *x31
        +coeff[  9]*x11    *x31
        +coeff[ 10]    *x23
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]*x13*x22
        +coeff[ 13]*x12*x21*x31*x41
        +coeff[ 14]    *x23*x31*x41
        +coeff[ 15]    *x21        *x51
        +coeff[ 16]*x13
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[ 17]*x11*x22
        +coeff[ 18]    *x23    *x41
        +coeff[ 19]*x11*x22*x31*x41
        +coeff[ 20]*x11*x22    *x42
        +coeff[ 21]*x11*x22    *x43
        +coeff[ 22]    *x22
        +coeff[ 23]*x12*x21
        +coeff[ 24]*x11        *x42
        +coeff[ 25]    *x21    *x41*x51
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[ 26]    *x21*x32
        +coeff[ 27]    *x21*x31*x41
        +coeff[ 28]*x11*x22    *x41
        +coeff[ 29]    *x23*x31
        +coeff[ 30]    *x21    *x43
        +coeff[ 31]    *x23*x31*x42
        +coeff[ 32]*x11            *x51
        +coeff[ 33]*x11*x21
        +coeff[ 34]            *x42
    ;
    v_x_l5p77_q1en                            =v_x_l5p77_q1en
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]    *x21*x31    *x51
        +coeff[ 37]*x12*x21    *x41
        +coeff[ 38]*x11*x22*x31
        +coeff[ 39]        *x33*x41
        +coeff[ 40]    *x21*x31*x42
        +coeff[ 41]    *x23*x32*x41
        ;

    return v_x_l5p77_q1en                            ;
}
float t_l5p77_q1en                            (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.1070047E-02;
    float xmin[10]={
        -0.14934E-01,-0.45642E-01,-0.14989E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26979E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.16119963E-03,-0.23088937E-02, 0.36846112E-01, 0.10805798E-01,
        -0.16123441E-03,-0.12341746E-02, 0.19103953E-04, 0.85140805E-03,
         0.83960372E-03, 0.37708865E-02, 0.23460118E-02, 0.29090596E-02,
        -0.36508117E-02, 0.19374516E-03, 0.28919536E-03, 0.38208021E-03,
         0.20546769E-02,-0.38103047E-02, 0.17083794E-05, 0.99784011E-04,
         0.37815346E-03, 0.47848304E-03,-0.26091600E-02,-0.66329731E-03,
        -0.39901151E-03,-0.22676212E-02, 0.14893463E-02, 0.15263744E-03,
        -0.18571963E-03, 0.76793847E-04,-0.44372273E-03,-0.51352323E-03,
        -0.75114053E-03,-0.47010730E-03, 0.81317529E-06, 0.13518051E-02,
        -0.11334268E-03, 0.21206593E-03,-0.12930934E-03, 0.40489892E-03,
        -0.19676929E-04,-0.27072678E-04, 0.59159789E-04,-0.71109644E-04,
         0.73979179E-04,-0.70175309E-04, 0.16560641E-03, 0.30135788E-03,
         0.10742922E-03, 0.14458403E-03,
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
        +coeff[ 11]    *x23
        +coeff[ 12]    *x21    *x42
        +coeff[ 13]    *x23*x31*x41
        +coeff[ 14]        *x31
        +coeff[ 15]    *x21        *x51
        +coeff[ 16]*x11*x22
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 17]    *x23    *x41
        +coeff[ 18]*x13*x22    *x41
        +coeff[ 19]                *x51
        +coeff[ 20]    *x22
        +coeff[ 21]*x12*x21
        +coeff[ 22]    *x21*x31*x41
        +coeff[ 23]*x11        *x42
        +coeff[ 24]    *x21    *x41*x51
        +coeff[ 25]*x11*x22    *x41
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 26]    *x21    *x43
        +coeff[ 27]*x11*x21
        +coeff[ 28]            *x42
        +coeff[ 29]*x11            *x51
        +coeff[ 30]    *x21*x32
        +coeff[ 31]*x11    *x31*x41
        +coeff[ 32]*x11*x22*x31
        +coeff[ 33]*x12*x21    *x41
        +coeff[ 34]        *x33*x41
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 35]    *x21*x31*x42
        +coeff[ 36]        *x31*x41
        +coeff[ 37]    *x22    *x41
        +coeff[ 38]    *x21*x31    *x51
        +coeff[ 39]    *x21*x32*x41
        +coeff[ 40]        *x32
        +coeff[ 41]            *x41*x51
        +coeff[ 42]    *x22*x31
        +coeff[ 43]*x11    *x32
    ;
    v_t_l5p77_q1en                            =v_t_l5p77_q1en
        +coeff[ 44]*x11*x21    *x41
        +coeff[ 45]*x11        *x41*x51
        +coeff[ 46]*x11        *x43
        +coeff[ 47]    *x23    *x42
        +coeff[ 48]*x11    *x31*x43
        +coeff[ 49]*x13    *x31*x42
        ;

    return v_t_l5p77_q1en                            ;
}
float y_l5p77_q1en                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.8726323E-02;
    float xmin[10]={
        -0.14934E-01,-0.45642E-01,-0.14989E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26979E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.99162757E-02, 0.73314637E-01, 0.16859217E-01,-0.41868044E-02,
         0.70912768E-02, 0.32332088E-02, 0.88447391E-03,-0.29507740E-02,
        -0.53010066E-02,-0.21017937E-02,-0.25065481E-02, 0.20428335E-03,
         0.39329237E-03,-0.19801604E-02,-0.96367742E-03, 0.31444731E-02,
        -0.22597609E-03,-0.38319640E-03, 0.21446146E-03, 0.22473524E-02,
         0.14120407E-02, 0.10353411E-02,-0.63800755E-04,-0.17850954E-03,
         0.55803260E-03, 0.55929046E-03,-0.29515885E-03,-0.30685210E-03,
        -0.13741088E-03,-0.97911293E-03,-0.88839419E-03,-0.39427116E-04,
        -0.35826160E-04,-0.72870425E-04, 0.20523799E-03,-0.12258330E-03,
         0.37406743E-03,-0.32798050E-03,-0.14370842E-04,-0.11208567E-03,
        -0.82795683E-04, 0.91568174E-04, 0.95379983E-04, 0.72246934E-04,
         0.68374815E-04, 0.15781494E-03, 0.11872471E-03, 0.16902655E-03,
        -0.57056092E-03, 0.85300482E-04,-0.73567574E-03, 0.71623555E-03,
        -0.25399387E-03, 0.65238279E-03, 0.19456253E-03, 0.28879142E-04,
        -0.14108260E-04, 0.76129538E-04,-0.25915087E-03,-0.30485954E-03,
        -0.47251116E-04, 0.30484985E-03,-0.96356314E-04, 0.29172795E-03,
         0.94147297E-04,-0.82646757E-05,-0.13379679E-04, 0.13189948E-04,
         0.13139383E-04,-0.26218846E-04,-0.61732331E-04,-0.11905923E-04,
        -0.35729998E-04,-0.16365297E-04,-0.74161171E-05,-0.53622422E-04,
        -0.23592067E-04, 0.64603402E-04, 0.18756953E-04, 0.30137173E-05,
        -0.27986227E-04,-0.46783278E-03,-0.38089367E-03,-0.38100299E-03,
        -0.30592456E-03, 0.11600891E-03,-0.90200498E-04, 0.55008426E-04,
        -0.78134326E-04,-0.39047827E-05, 0.77869527E-05, 0.49372138E-05,
         0.29756677E-04, 0.23735562E-04, 0.35460111E-04, 0.19291509E-04,
         0.12223386E-04,
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
        +coeff[  7]            *x42
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[  8]    *x22    *x41
        +coeff[  9]        *x31*x41
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]*x11
        +coeff[ 12]*x12
        +coeff[ 13]    *x22*x31
        +coeff[ 14]*x11*x21*x31
        +coeff[ 15]    *x22    *x42
        +coeff[ 16]    *x21    *x41
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 17]        *x32
        +coeff[ 18]                *x52
        +coeff[ 19]    *x22*x31*x41
        +coeff[ 20]*x11*x21    *x42
        +coeff[ 21]*x11*x21*x31*x41
        +coeff[ 22]    *x21*x31
        +coeff[ 23]            *x41*x51
        +coeff[ 24]            *x43
        +coeff[ 25]        *x31*x42
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 26]    *x22        *x51
        +coeff[ 27]*x12        *x41
        +coeff[ 28]*x11*x21        *x51
        +coeff[ 29]    *x24
        +coeff[ 30]*x11*x23
        +coeff[ 31]*x11        *x41
        +coeff[ 32]    *x21        *x51
        +coeff[ 33]        *x31    *x51
        +coeff[ 34]        *x32*x41
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 35]*x12    *x31
        +coeff[ 36]    *x22*x32
        +coeff[ 37]*x12*x22
        +coeff[ 38]*x11    *x31
        +coeff[ 39]    *x21    *x42
        +coeff[ 40]    *x21*x31*x41
        +coeff[ 41]            *x42*x51
        +coeff[ 42]    *x23
        +coeff[ 43]        *x31*x41*x51
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 44]*x11*x22
        +coeff[ 45]    *x22    *x41*x51
        +coeff[ 46]*x12        *x42
        +coeff[ 47]*x11*x21*x32
        +coeff[ 48]    *x22    *x43
        +coeff[ 49]*x12    *x31*x41
        +coeff[ 50]    *x22*x31*x42
        +coeff[ 51]    *x24    *x41
        +coeff[ 52]    *x22*x32*x41
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 53]*x11*x23    *x41
        +coeff[ 54]*x12*x22    *x41
        +coeff[ 55]        *x33
        +coeff[ 56]*x12            *x51
        +coeff[ 57]*x11*x21    *x41*x51
        +coeff[ 58]*x11*x21    *x43
        +coeff[ 59]*x11*x21*x31*x42
        +coeff[ 60]*x13*x21
        +coeff[ 61]    *x24*x31
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 62]*x11*x21*x32*x41
        +coeff[ 63]*x11*x23*x31
        +coeff[ 64]*x12*x22*x31
        +coeff[ 65]*x11            *x51
        +coeff[ 66]    *x21*x32
        +coeff[ 67]        *x32    *x51
        +coeff[ 68]*x12*x21
        +coeff[ 69]*x11        *x43
        +coeff[ 70]    *x23    *x41
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 71]*x11    *x31*x42
        +coeff[ 72]    *x23*x31
        +coeff[ 73]            *x45
        +coeff[ 74]                *x53
        +coeff[ 75]*x11*x22    *x41
        +coeff[ 76]*x11*x22*x31
        +coeff[ 77]    *x22*x31    *x51
        +coeff[ 78]*x12    *x32
        +coeff[ 79]*x11*x21*x31    *x51
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 80]    *x22*x33
        +coeff[ 81]    *x24    *x42
        +coeff[ 82]    *x24*x31*x41
        +coeff[ 83]*x11*x23    *x42
        +coeff[ 84]*x11*x23*x31*x41
        +coeff[ 85]    *x22*x31*x44
        +coeff[ 86]    *x22    *x44*x51
        +coeff[ 87]*x11*x23*x31    *x51
        +coeff[ 88]    *x22*x31*x43*x51
    ;
    v_y_l5p77_q1en                            =v_y_l5p77_q1en
        +coeff[ 89]*x11    *x31*x41
        +coeff[ 90]    *x21    *x41*x51
        +coeff[ 91]    *x21*x31    *x51
        +coeff[ 92]    *x21    *x43
        +coeff[ 93]        *x31*x43
        +coeff[ 94]    *x21*x31*x42
        +coeff[ 95]        *x32*x42
        +coeff[ 96]    *x21*x32*x41
        ;

    return v_y_l5p77_q1en                            ;
}
float p_l5p77_q1en                            (float *x,int m){
    int ncoeff= 47;
    float avdat=  0.8072414E-02;
    float xmin[10]={
        -0.14934E-01,-0.45642E-01,-0.14989E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26979E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 48]={
        -0.55424706E-02, 0.12215672E-02, 0.22678166E-02, 0.36800280E-01,
        -0.62580779E-02, 0.96053528E-02,-0.36301564E-02, 0.40409467E-02,
        -0.57084318E-02,-0.24752934E-02,-0.26294240E-02, 0.27406219E-03,
        -0.45561083E-03, 0.31980101E-03,-0.20589984E-02, 0.49046276E-03,
        -0.96219749E-03, 0.36647108E-02, 0.32081132E-03,-0.40744347E-03,
        -0.40927451E-03, 0.14955814E-02, 0.27798835E-03,-0.13366486E-03,
        -0.17450341E-03,-0.86257693E-04, 0.47690806E-03,-0.14798625E-02,
         0.21910954E-02,-0.16821596E-03, 0.72945062E-04,-0.12013969E-02,
        -0.26321458E-03, 0.88324118E-03,-0.37596826E-03,-0.27189524E-04,
         0.17821870E-03,-0.36655969E-04,-0.56083205E-04, 0.79556630E-04,
        -0.87137305E-04, 0.45861481E-03, 0.77505218E-04, 0.32226401E-03,
        -0.88652632E-04, 0.16227867E-03, 0.12870852E-03,
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
        +coeff[  9]        *x31*x41
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]*x11
        +coeff[ 12]    *x21    *x41
        +coeff[ 13]                *x52
        +coeff[ 14]    *x22*x31
        +coeff[ 15]*x12
        +coeff[ 16]*x11*x21*x31
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 17]    *x22    *x42
        +coeff[ 18]    *x24*x31*x41
        +coeff[ 19]        *x32
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]*x11*x21    *x42
        +coeff[ 22]*x11*x23*x31*x41
        +coeff[ 23]    *x21*x31
        +coeff[ 24]            *x41*x51
        +coeff[ 25]*x11        *x41
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 26]            *x43
        +coeff[ 27]    *x24
        +coeff[ 28]    *x22*x31*x41
        +coeff[ 29]*x11*x21        *x51
        +coeff[ 30]        *x31*x43
        +coeff[ 31]*x11*x23
        +coeff[ 32]*x12        *x41
        +coeff[ 33]*x11*x21*x31*x41
        +coeff[ 34]*x12*x22
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 35]        *x31*x44
        +coeff[ 36]    *x22*x32*x42
        +coeff[ 37]    *x21        *x51
        +coeff[ 38]        *x31    *x51
        +coeff[ 39]    *x23
        +coeff[ 40]    *x21    *x42
        +coeff[ 41]        *x31*x42
        +coeff[ 42]            *x42*x51
        +coeff[ 43]    *x22*x32
    ;
    v_p_l5p77_q1en                            =v_p_l5p77_q1en
        +coeff[ 44]*x12    *x31
        +coeff[ 45]*x11*x21*x32
        +coeff[ 46]        *x34*x41
        ;

    return v_p_l5p77_q1en                            ;
}
float l_l5p77_q1en                            (float *x,int m){
    int ncoeff= 58;
    float avdat= -0.7516675E-03;
    float xmin[10]={
        -0.14934E-01,-0.45642E-01,-0.14989E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26979E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 59]={
         0.11631451E-03, 0.59104888E-02,-0.19520973E-02,-0.97936008E-05,
        -0.12933025E-03,-0.12886937E-02, 0.10167957E-05,-0.18387733E-03,
        -0.22396102E-06, 0.67192625E-04, 0.28087737E-06,-0.36788197E-08,
        -0.74087512E-06,-0.25828251E-05, 0.25215080E-04, 0.88400634E-06,
         0.14945741E-04, 0.34620944E-05, 0.81371303E-06,-0.15280737E-06,
        -0.33173530E-05,-0.20038508E-04,-0.46579495E-04, 0.18284371E-02,
        -0.10794187E-03, 0.17835488E-03,-0.74954593E-03,-0.26847186E-03,
        -0.57133842E-04, 0.14305365E-03, 0.99303194E-04, 0.34375035E-03,
         0.15061707E-04,-0.50201514E-04,-0.16170656E-03,-0.11433541E-03,
         0.10729399E-03, 0.14311474E-04, 0.52511559E-05,-0.10573845E-04,
        -0.57212901E-05,-0.69438279E-05, 0.16103261E-03,-0.11470186E-03,
         0.11484281E-03, 0.38194880E-05, 0.12506042E-04,-0.80558802E-05,
         0.13074266E-04, 0.80367327E-05,-0.19396957E-04, 0.18410292E-03,
         0.43060136E-04,-0.27983751E-04,-0.86535132E-04,-0.99056097E-05,
         0.30804746E-04, 0.27614586E-04,
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

    float v_l_l5p77_q1en                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x22
        +coeff[  3]        *x32
        +coeff[  4]        *x31*x41
        +coeff[  5]            *x42
        +coeff[  6]*x11            *x51
        +coeff[  7]    *x22*x31
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[  8]        *x33
        +coeff[  9]        *x31*x42
        +coeff[ 10]        *x32    *x51
        +coeff[ 11]        *x31    *x52
        +coeff[ 12]*x12    *x31
        +coeff[ 13]            *x43*x51
        +coeff[ 14]    *x24*x31
        +coeff[ 15]    *x22*x33
        +coeff[ 16]        *x34*x41
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 17]        *x33*x42
        +coeff[ 18]        *x33    *x52
        +coeff[ 19]*x12    *x33
        +coeff[ 20]    *x24*x33
        +coeff[ 21]    *x22*x34*x41
        +coeff[ 22]    *x21
        +coeff[ 23]        *x31
        +coeff[ 24]                *x51
        +coeff[ 25]*x11*x21
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 26]    *x22    *x41
        +coeff[ 27]*x11*x21    *x41
        +coeff[ 28]    *x21    *x41
        +coeff[ 29]            *x41*x51
        +coeff[ 30]            *x43
        +coeff[ 31]    *x22    *x42
        +coeff[ 32]    *x24*x31*x41
        +coeff[ 33]*x11*x21*x31
        +coeff[ 34]    *x24
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 35]*x11*x23
        +coeff[ 36]*x11*x21    *x42
        +coeff[ 37]*x11*x23*x31*x41
        +coeff[ 38]*x11
        +coeff[ 39]    *x21*x31
        +coeff[ 40]                *x52
        +coeff[ 41]*x11        *x41
        +coeff[ 42]    *x24    *x41
        +coeff[ 43]    *x22    *x43
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 44]*x11*x23    *x41
        +coeff[ 45]        *x31    *x51
        +coeff[ 46]*x12
        +coeff[ 47]    *x23
        +coeff[ 48]    *x22        *x51
        +coeff[ 49]*x11*x21        *x51
        +coeff[ 50]*x12        *x41
        +coeff[ 51]    *x22*x31*x41
        +coeff[ 52]*x11*x21*x31*x41
    ;
    v_l_l5p77_q1en                            =v_l_l5p77_q1en
        +coeff[ 53]*x12*x22
        +coeff[ 54]    *x22*x31*x42
        +coeff[ 55]            *x43*x52
        +coeff[ 56]*x12*x22    *x41
        +coeff[ 57]    *x24*x32
        ;

    return v_l_l5p77_q1en                            ;
}
float x_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 48;
    float avdat= -0.5084220E-03;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 49]={
         0.89667126E-04, 0.45305328E-02, 0.11675917E+00, 0.30337854E-02,
         0.21381993E-01,-0.35925451E-03, 0.16905329E-02, 0.48775827E-02,
         0.75667324E-02, 0.61490422E-03, 0.17638762E-02, 0.57726931E-02,
        -0.75924457E-02, 0.88173045E-04,-0.91345028E-05, 0.70145697E-03,
         0.64191976E-04, 0.99775498E-03, 0.40827938E-02,-0.53038653E-02,
        -0.79275118E-02, 0.22904988E-03, 0.39285942E-03,-0.95495733E-03,
        -0.48723575E-02,-0.25206199E-02, 0.13897605E-04,-0.11716212E-04,
        -0.21687623E-03, 0.28406788E-03,-0.37303404E-03,-0.10717461E-02,
        -0.15651096E-02,-0.43657093E-03, 0.25333351E-03,-0.11008385E-02,
        -0.15428085E-02, 0.20897323E-04, 0.25425681E-02, 0.11775506E-02,
        -0.44111039E-04,-0.25684136E-03,-0.12019317E-03,-0.18139792E-03,
        -0.18113227E-03, 0.19372735E-02, 0.21094440E-03, 0.14106415E-02,
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

    float v_x_l5p77_q1ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]    *x21        *x51
        +coeff[  4]    *x21    *x41
        +coeff[  5]*x12*x21*x31
        +coeff[  6]            *x41
        +coeff[  7]*x11        *x41
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[  8]    *x21*x31
        +coeff[  9]        *x31
        +coeff[ 10]*x11    *x31
        +coeff[ 11]    *x23
        +coeff[ 12]    *x21    *x42
        +coeff[ 13]*x13*x22
        +coeff[ 14]    *x23*x31*x41
        +coeff[ 15]    *x22
        +coeff[ 16]*x13
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 17]*x12*x21
        +coeff[ 18]*x11*x22
        +coeff[ 19]    *x21*x31*x41
        +coeff[ 20]    *x23    *x41
        +coeff[ 21]                *x51
        +coeff[ 22]*x11            *x51
        +coeff[ 23]    *x21*x32
        +coeff[ 24]*x11*x22    *x41
        +coeff[ 25]    *x23*x31
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 26]*x13        *x42
        +coeff[ 27]*x11    *x33*x41
        +coeff[ 28]*x11*x22*x33
        +coeff[ 29]*x11*x21
        +coeff[ 30]            *x42
        +coeff[ 31]*x11    *x31*x41
        +coeff[ 32]*x11        *x42
        +coeff[ 33]    *x21    *x41*x51
        +coeff[ 34]    *x22    *x41
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 35]*x12*x21    *x41
        +coeff[ 36]*x11*x22*x31
        +coeff[ 37]        *x33*x41
        +coeff[ 38]    *x21    *x43
        +coeff[ 39]    *x23*x31*x42
        +coeff[ 40]        *x32
        +coeff[ 41]        *x31*x41
        +coeff[ 42]    *x21        *x52
        +coeff[ 43]*x11    *x32
    ;
    v_x_l5p77_q1ex                            =v_x_l5p77_q1ex
        +coeff[ 44]    *x21*x31    *x51
        +coeff[ 45]    *x21*x31*x42
        +coeff[ 46]*x12*x21*x32*x41
        +coeff[ 47]    *x23*x32*x41
        ;

    return v_x_l5p77_q1ex                            ;
}
float t_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.1177512E-02;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
        -0.53511136E-04,-0.62986016E-02,-0.63281581E-02, 0.42235677E-03,
         0.11151687E-02, 0.53156484E-02, 0.22688343E-02, 0.11753939E-02,
        -0.11765899E-03,-0.79402601E-03,-0.17245058E-05, 0.30667892E-04,
        -0.24854651E-05,-0.34417881E-03, 0.78964858E-05, 0.14389788E-03,
         0.21527421E-03, 0.39157338E-03, 0.18812015E-02, 0.85975387E-03,
        -0.12912990E-02,-0.18889785E-02,-0.22808381E-02, 0.64783831E-04,
         0.23833156E-03, 0.20027281E-03,-0.18708757E-03,-0.13521927E-02,
         0.79352743E-04,-0.10815091E-03,-0.36847338E-03, 0.18112535E-03,
        -0.97784992E-04,-0.42883464E-03,-0.22955333E-03, 0.45363419E-06,
        -0.42382282E-04,-0.65726104E-04,-0.20622244E-03, 0.33800079E-04,
         0.34279100E-03,-0.51083352E-04, 0.26292357E-03, 0.14898460E-03,
         0.55751088E-03,-0.14202880E-04, 0.29413854E-04, 0.35698111E-04,
         0.26666383E-04, 0.78248406E-04,
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
        +coeff[  3]            *x41
        +coeff[  4]*x11        *x41
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x21        *x51
        +coeff[  7]    *x23
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[  8]*x12*x21*x31
        +coeff[  9]    *x23*x31
        +coeff[ 10]*x11    *x33
        +coeff[ 11]    *x21*x33
        +coeff[ 12]*x13            *x51
        +coeff[ 13]    *x23    *x42
        +coeff[ 14]    *x23*x33
        +coeff[ 15]        *x31
        +coeff[ 16]    *x22
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 17]*x11    *x31
        +coeff[ 18]    *x21*x31
        +coeff[ 19]*x11*x22
        +coeff[ 20]    *x21*x31*x41
        +coeff[ 21]    *x21    *x42
        +coeff[ 22]    *x23    *x41
        +coeff[ 23]                *x51
        +coeff[ 24]*x11            *x51
        +coeff[ 25]*x12*x21
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 26]    *x21*x32
        +coeff[ 27]*x11*x22    *x41
        +coeff[ 28]*x11*x21
        +coeff[ 29]            *x42
        +coeff[ 30]*x11        *x42
        +coeff[ 31]    *x21    *x41*x51
        +coeff[ 32]    *x21        *x52
        +coeff[ 33]*x11*x22*x31
        +coeff[ 34]*x12*x21    *x41
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 35]        *x33*x41
        +coeff[ 36]*x13    *x31*x41
        +coeff[ 37]        *x31*x41
        +coeff[ 38]*x11    *x31*x41
        +coeff[ 39]    *x21*x31    *x51
        +coeff[ 40]    *x21    *x43
        +coeff[ 41]*x13    *x32
        +coeff[ 42]    *x21*x31*x43
        +coeff[ 43]*x12*x21*x31*x42
    ;
    v_t_l5p77_q1ex                            =v_t_l5p77_q1ex
        +coeff[ 44]    *x23*x31*x42
        +coeff[ 45]        *x32
        +coeff[ 46]*x13
        +coeff[ 47]    *x22    *x41
        +coeff[ 48]*x11        *x41*x51
        +coeff[ 49]    *x21*x32*x41
        ;

    return v_t_l5p77_q1ex                            ;
}
float y_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1440121E-01;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.10022627E-01, 0.17543498E+00, 0.36949359E-02, 0.28414050E-01,
        -0.18457061E-01,-0.11674689E-01, 0.29147293E-01, 0.12685195E-01,
        -0.83160531E-02,-0.18992381E-01,-0.88637024E-02, 0.81490370E-03,
        -0.27667978E-02, 0.14635699E-02,-0.70016016E-02,-0.33664184E-02,
        -0.11962499E-02,-0.14927927E-02, 0.11233740E-02, 0.11944790E-01,
         0.81291338E-02, 0.16488932E-03,-0.33776314E-03,-0.66234486E-03,
        -0.11001750E-02,-0.62647864E-03, 0.50842445E-02,-0.44403607E-02,
        -0.62767597E-03, 0.35549954E-02,-0.13145832E-02,-0.37718420E-02,
        -0.26429558E-03,-0.18139194E-03,-0.15879545E-02,-0.43007999E-03,
         0.14362308E-02,-0.89139980E-03,-0.12656800E-02,-0.90235015E-04,
        -0.50582684E-03, 0.24301633E-02,-0.37097832E-03, 0.37932582E-03,
         0.26619452E-03, 0.13027426E-03, 0.62646071E-03,-0.42208759E-03,
        -0.17845302E-03, 0.24729497E-02,-0.46727429E-04, 0.74679480E-03,
         0.33468747E-03, 0.56315708E-04,-0.73792908E-04, 0.65363856E-03,
         0.54526405E-03, 0.38739890E-03, 0.34724450E-03,-0.17101335E-05,
         0.21015282E-02,-0.20375724E-03, 0.80538582E-03, 0.20521458E-02,
         0.82582375E-03, 0.75251324E-03, 0.33807382E-03, 0.40336946E-03,
         0.12179973E-03, 0.69268528E-04, 0.22821686E-03, 0.23163017E-03,
        -0.16744452E-03,-0.85607433E-04,-0.54696036E-04, 0.25320772E-03,
        -0.21749267E-02,-0.23434730E-02, 0.56539284E-04, 0.13808609E-03,
        -0.10049208E-02,-0.96441712E-03,-0.10919577E-02,-0.16863314E-03,
        -0.63931358E-04,-0.37961168E-03, 0.15251478E-03, 0.53731535E-04,
         0.52306795E-03, 0.31851698E-04, 0.11594896E-02, 0.81847282E-03,
         0.14128796E-04, 0.27157266E-04, 0.70289127E-04, 0.21191926E-03,
         0.63246829E-04,
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
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]            *x42
        +coeff[  6]    *x22
        +coeff[  7]*x11*x21
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[  8]        *x31*x41
        +coeff[  9]    *x22    *x41
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]*x11
        +coeff[ 12]            *x41*x51
        +coeff[ 13]*x12
        +coeff[ 14]    *x22*x31
        +coeff[ 15]*x11*x21*x31
        +coeff[ 16]    *x21    *x41
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 17]        *x32
        +coeff[ 18]                *x52
        +coeff[ 19]    *x22    *x42
        +coeff[ 20]    *x22*x31*x41
        +coeff[ 21]    *x24        *x51
        +coeff[ 22]    *x21*x31
        +coeff[ 23]        *x31    *x51
        +coeff[ 24]*x12        *x41
        +coeff[ 25]*x11*x21        *x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 26]*x11*x21    *x42
        +coeff[ 27]    *x24
        +coeff[ 28]            *x45
        +coeff[ 29]*x11*x21*x31*x41
        +coeff[ 30]        *x31*x44
        +coeff[ 31]*x11*x23
        +coeff[ 32]*x11        *x41
        +coeff[ 33]    *x21        *x51
        +coeff[ 34]    *x22        *x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 35]*x12    *x31
        +coeff[ 36]    *x22*x32
        +coeff[ 37]        *x32*x43
        +coeff[ 38]*x12*x22
        +coeff[ 39]*x11    *x31
        +coeff[ 40]    *x21    *x42
        +coeff[ 41]        *x31*x42
        +coeff[ 42]    *x21*x31*x41
        +coeff[ 43]    *x23
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 44]*x11*x22
        +coeff[ 45]            *x41*x52
        +coeff[ 46]*x11*x21*x32
        +coeff[ 47]        *x33*x42
        +coeff[ 48]        *x32*x45
        +coeff[ 49]            *x43
        +coeff[ 50]*x11            *x51
        +coeff[ 51]        *x32*x41
        +coeff[ 52]        *x31*x41*x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 53]*x12*x21
        +coeff[ 54]*x12            *x51
        +coeff[ 55]    *x22    *x41*x51
        +coeff[ 56]*x12        *x42
        +coeff[ 57]*x12    *x31*x41
        +coeff[ 58]*x11*x21    *x41*x51
        +coeff[ 59]            *x44*x51
        +coeff[ 60]    *x24    *x41
        +coeff[ 61]*x13*x21
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 62]    *x24*x31
        +coeff[ 63]*x11*x23    *x41
        +coeff[ 64]*x11*x23*x31
        +coeff[ 65]*x12*x22    *x41
        +coeff[ 66]*x12*x22*x31
        +coeff[ 67]            *x42*x51
        +coeff[ 68]        *x33
        +coeff[ 69]        *x32    *x51
        +coeff[ 70]    *x21    *x43
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 71]    *x21*x31*x42
        +coeff[ 72]    *x23    *x41
        +coeff[ 73]    *x23*x31
        +coeff[ 74]                *x53
        +coeff[ 75]    *x22*x31    *x51
        +coeff[ 76]    *x22    *x43
        +coeff[ 77]    *x22*x31*x42
        +coeff[ 78]*x12    *x32
        +coeff[ 79]*x11*x21*x31    *x51
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 80]*x11*x21    *x43
        +coeff[ 81]    *x22*x32*x41
        +coeff[ 82]*x11*x21*x31*x42
        +coeff[ 83]    *x22*x33
        +coeff[ 84]    *x21*x34
        +coeff[ 85]*x11*x21*x32*x41
        +coeff[ 86]    *x23*x31*x41*x51
        +coeff[ 87]    *x21    *x41*x51
        +coeff[ 88]            *x44
    ;
    v_y_l5p77_q1ex                            =v_y_l5p77_q1ex
        +coeff[ 89]    *x21*x31    *x51
        +coeff[ 90]        *x31*x43
        +coeff[ 91]        *x32*x42
        +coeff[ 92]*x11    *x31    *x51
        +coeff[ 93]        *x31    *x52
        +coeff[ 94]    *x21*x32*x41
        +coeff[ 95]        *x33*x41
        +coeff[ 96]    *x21    *x42*x51
        ;

    return v_y_l5p77_q1ex                            ;
}
float p_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 44;
    float avdat=  0.7687971E-02;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 45]={
        -0.45596543E-02, 0.94885835E-02, 0.74260212E-01,-0.95395073E-02,
         0.14814159E-01,-0.22607229E-02, 0.63278377E-02, 0.13263938E-03,
         0.18681972E-02,-0.56298459E-02,-0.91800485E-02,-0.42138994E-02,
         0.38433545E-02,-0.30400666E-04, 0.43376273E-03,-0.69415389E-03,
        -0.38904864E-02, 0.67411043E-03,-0.33499941E-02, 0.77165460E-03,
        -0.15638145E-02, 0.59423954E-02,-0.65426138E-03,-0.46462484E-03,
        -0.89728174E-03, 0.24092312E-02,-0.19505055E-03,-0.12153704E-03,
         0.10448048E-02,-0.23285530E-02,-0.37872684E-03,-0.19622173E-02,
        -0.46528454E-03, 0.16128626E-02,-0.62409311E-03, 0.25166548E-03,
        -0.98858574E-04, 0.10538647E-04, 0.77666028E-03, 0.13031877E-03,
         0.65589463E-03,-0.16553215E-03, 0.29849555E-03, 0.36125074E-03,
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
        +coeff[ 15]        *x32
        +coeff[ 16]        *x31*x41
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[ 17]                *x52
        +coeff[ 18]    *x22*x31
        +coeff[ 19]*x12
        +coeff[ 20]*x11*x21*x31
        +coeff[ 21]    *x22    *x42
        +coeff[ 22]    *x21    *x41
        +coeff[ 23]        *x31    *x51
        +coeff[ 24]    *x22        *x51
        +coeff[ 25]*x11*x21    *x42
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[ 26]    *x21*x31
        +coeff[ 27]*x11        *x41
        +coeff[ 28]            *x43
        +coeff[ 29]    *x24
        +coeff[ 30]*x11*x21        *x51
        +coeff[ 31]*x11*x23
        +coeff[ 32]*x12        *x41
        +coeff[ 33]*x11*x21*x31*x41
        +coeff[ 34]*x12*x22
    ;
    v_p_l5p77_q1ex                            =v_p_l5p77_q1ex
        +coeff[ 35]        *x33*x42
        +coeff[ 36]    *x21        *x51
        +coeff[ 37]        *x33
        +coeff[ 38]        *x31*x42
        +coeff[ 39]            *x41*x52
        +coeff[ 40]    *x22*x32
        +coeff[ 41]*x12    *x31
        +coeff[ 42]*x11*x21*x32
        +coeff[ 43]        *x34*x41
        ;

    return v_p_l5p77_q1ex                            ;
}
float l_l5p77_q1ex                            (float *x,int m){
    int ncoeff= 86;
    float avdat= -0.1754399E-02;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 87]={
         0.12203738E-02,-0.39713708E-04, 0.55235615E-02,-0.23757643E-02,
        -0.26268277E-04,-0.86406234E-03,-0.44509228E-02, 0.10149255E-02,
         0.18696385E-03, 0.73239262E-05,-0.41254426E-03, 0.58956143E-05,
        -0.22622878E-02, 0.10394066E-03, 0.47833219E-03, 0.60844439E-03,
        -0.10537148E-04,-0.80097729E-03,-0.48696497E-05, 0.68225905E-04,
         0.90995763E-06,-0.34551380E-04,-0.11950768E-03, 0.15227505E-03,
         0.69208591E-05,-0.50831864E-05,-0.17172883E-05, 0.83711944E-06,
        -0.25667108E-04, 0.17952451E-02,-0.22556966E-03, 0.10447425E-03,
        -0.69040805E-04, 0.12817987E-02, 0.18764863E-03, 0.98566697E-05,
         0.99726931E-04,-0.30071868E-03, 0.62121253E-03,-0.22828943E-03,
         0.52079745E-03,-0.64572800E-04, 0.26197285E-04,-0.34163244E-04,
        -0.13585410E-04, 0.80430000E-05, 0.10145702E-03,-0.85862928E-04,
        -0.12138014E-03,-0.77865290E-04, 0.23865041E-03, 0.61080465E-03,
         0.44878601E-03,-0.61930809E-05,-0.57135891E-04,-0.50445313E-04,
         0.63500775E-04,-0.23834456E-04,-0.64506086E-04,-0.70060120E-03,
        -0.13965434E-03, 0.22582972E-04, 0.11181350E-04,-0.12286393E-03,
        -0.25038089E-04,-0.14365691E-04, 0.40817453E-04,-0.44574618E-03,
         0.11353014E-03,-0.26161553E-03, 0.11651635E-03,-0.20547492E-03,
        -0.13582040E-03,-0.53497852E-05, 0.86797036E-05, 0.13749663E-04,
        -0.27763934E-04, 0.21900822E-04,-0.92837136E-05,-0.13833393E-04,
         0.16164324E-04,-0.15736476E-03,-0.17068961E-03,-0.11709079E-03,
        -0.13648436E-03, 0.36145451E-04,
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
        +coeff[  1]    *x21
        +coeff[  2]            *x41
        +coeff[  3]    *x22
        +coeff[  4]    *x21*x31
        +coeff[  5]        *x31*x41
        +coeff[  6]            *x42
        +coeff[  7]            *x41*x51
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[  8]*x11*x21
        +coeff[  9]*x11            *x51
        +coeff[ 10]    *x22*x31
        +coeff[ 11]        *x33
        +coeff[ 12]    *x22    *x41
        +coeff[ 13]        *x32*x41
        +coeff[ 14]        *x31*x42
        +coeff[ 15]            *x43
        +coeff[ 16]        *x31    *x52
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 17]*x11*x21    *x41
        +coeff[ 18]*x12    *x31
        +coeff[ 19]    *x22*x32
        +coeff[ 20]        *x34
        +coeff[ 21]        *x32*x42
        +coeff[ 22]        *x31*x43
        +coeff[ 23]    *x24*x31
        +coeff[ 24]    *x22*x33
        +coeff[ 25]        *x33*x42
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 26]        *x33    *x52
        +coeff[ 27]*x12    *x33
        +coeff[ 28]    *x24*x33
        +coeff[ 29]        *x31
        +coeff[ 30]    *x21    *x41
        +coeff[ 31]        *x31    *x51
        +coeff[ 32]                *x52
        +coeff[ 33]    *x22    *x42
        +coeff[ 34]    *x22        *x51
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 35]                *x53
        +coeff[ 36]*x11*x21        *x51
        +coeff[ 37]    *x24
        +coeff[ 38]    *x22*x31*x41
        +coeff[ 39]*x11*x23
        +coeff[ 40]*x11*x21    *x42
        +coeff[ 41]                *x51
        +coeff[ 42]    *x21        *x51
        +coeff[ 43]*x11        *x41
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 44]*x12
        +coeff[ 45]        *x31*x41*x51
        +coeff[ 46]            *x42*x51
        +coeff[ 47]            *x41*x52
        +coeff[ 48]*x11*x21*x31
        +coeff[ 49]*x12        *x41
        +coeff[ 50]*x11*x21*x31*x41
        +coeff[ 51]    *x24    *x41
        +coeff[ 52]*x11*x23    *x41
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 53]*x11
        +coeff[ 54]        *x32
        +coeff[ 55]    *x23
        +coeff[ 56]    *x21    *x42
        +coeff[ 57]*x11*x22
        +coeff[ 58]*x12*x22
        +coeff[ 59]    *x22    *x43
        +coeff[ 60]    *x24*x31*x42
        +coeff[ 61]    *x21*x31*x41
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 62]*x12            *x51
        +coeff[ 63]            *x44
        +coeff[ 64]    *x22*x31    *x51
        +coeff[ 65]*x11*x21*x31    *x51
        +coeff[ 66]*x12        *x42
        +coeff[ 67]    *x22*x31*x42
        +coeff[ 68]*x11*x23*x31
        +coeff[ 69]*x11*x21    *x43
        +coeff[ 70]*x12*x22    *x41
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 71]    *x24*x32*x41
        +coeff[ 72]*x11*x23*x31*x42
        +coeff[ 73]        *x32    *x51
        +coeff[ 74]*x11        *x42
        +coeff[ 75]    *x21    *x43
        +coeff[ 76]    *x22        *x52
        +coeff[ 77]*x11*x21*x32
        +coeff[ 78]*x11*x21    *x41*x51
        +coeff[ 79]*x11*x21        *x52
    ;
    v_l_l5p77_q1ex                            =v_l_l5p77_q1ex
        +coeff[ 80]*x12    *x31*x41
        +coeff[ 81]*x11*x21*x31*x42
        +coeff[ 82]    *x24    *x42
        +coeff[ 83]*x11*x23    *x42
        +coeff[ 84]*x11*x23*x32*x41
        +coeff[ 85]*x12*x24*x31
        ;

    return v_l_l5p77_q1ex                            ;
}
float x_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.1309680E+01;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
        -0.10935051E-03,-0.22814563E-01, 0.19085447E+00, 0.52652992E-02,
         0.10765613E-01, 0.67431606E-01,-0.10105114E-03,-0.10692475E-02,
        -0.34385554E-04,-0.85247038E-02, 0.17906999E-02, 0.52959900E-02,
         0.14888884E-01, 0.23657860E-01, 0.16612660E-01,-0.26579140E-01,
         0.12633014E-03, 0.12413692E-02, 0.85938093E-03, 0.20795774E-02,
         0.23877434E-02, 0.17235684E-03, 0.12004319E-01,-0.18199248E-01,
        -0.27267504E-01,-0.19085393E-03, 0.28800236E-02,-0.46997583E-02,
        -0.30518032E-02,-0.16063238E-01, 0.17236034E-04, 0.95099368E-03,
        -0.11634597E-02,-0.31487765E-02,-0.33112578E-02,-0.52505694E-02,
        -0.12257276E-04,-0.74123865E-03,-0.78951987E-03,-0.55035076E-03,
         0.91353565E-03, 0.10866260E-02, 0.61999131E-02, 0.30688141E-02,
         0.38859004E-02, 0.54409094E-04, 0.31126896E-02,-0.16469405E-03,
        -0.15270387E-03, 0.34025215E-03,-0.60636312E-03,-0.14442399E-02,
         0.33914330E-03, 0.59964461E-03, 0.45459853E-02,-0.30826361E-03,
         0.12258007E-02, 0.73065888E-03, 0.41706889E-03, 0.28069543E-02,
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

    float v_x_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]            *x41
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x21    *x41
        +coeff[  6]*x13        *x41
        +coeff[  7]*x12*x21*x31
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[  8]    *x21*x31    *x52
        +coeff[  9]    *x23*x31
        +coeff[ 10]        *x31
        +coeff[ 11]*x11    *x31
        +coeff[ 12]*x11        *x41
        +coeff[ 13]    *x21*x31
        +coeff[ 14]    *x23
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x13*x22
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 17]    *x23*x31*x41
        +coeff[ 18]                *x51
        +coeff[ 19]*x11            *x51
        +coeff[ 20]    *x22
        +coeff[ 21]*x13
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]    *x23    *x41
        +coeff[ 25]*x13*x22    *x41
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 26]*x12*x21
        +coeff[ 27]*x11        *x42
        +coeff[ 28]    *x21*x32
        +coeff[ 29]*x11*x22    *x41
        +coeff[ 30]*x13    *x31*x41
        +coeff[ 31]*x11*x21
        +coeff[ 32]            *x42
        +coeff[ 33]*x11    *x31*x41
        +coeff[ 34]*x12*x21    *x41
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 35]*x11*x22*x31
        +coeff[ 36]        *x33*x41
        +coeff[ 37]        *x31*x41
        +coeff[ 38]    *x21        *x52
        +coeff[ 39]*x11    *x32
        +coeff[ 40]    *x22    *x41
        +coeff[ 41]    *x21    *x42*x51
        +coeff[ 42]    *x21    *x43
        +coeff[ 43]    *x21*x32*x42
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 44]    *x21*x31*x43
        +coeff[ 45]*x12*x21*x31*x42
        +coeff[ 46]    *x23*x31*x42
        +coeff[ 47]            *x41*x51
        +coeff[ 48]        *x32
        +coeff[ 49]*x11*x21    *x41
        +coeff[ 50]    *x21*x31    *x51
        +coeff[ 51]    *x21    *x41*x51
        +coeff[ 52]    *x22*x31
    ;
    v_x_l5p77_q2ex                            =v_x_l5p77_q2ex
        +coeff[ 53]    *x21*x31*x41*x51
        +coeff[ 54]    *x21*x31*x42
        +coeff[ 55]*x13        *x41*x51
        +coeff[ 56]    *x23    *x41*x51
        +coeff[ 57]    *x21*x33*x41
        +coeff[ 58]*x12*x21*x32*x41
        +coeff[ 59]    *x23*x32*x41
        ;

    return v_x_l5p77_q2ex                            ;
}
float t_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.5806373E+00;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.30308435E-03,-0.11053805E-01, 0.56641079E-01, 0.18744939E-02,
         0.24615837E-01, 0.31139341E-02,-0.44039727E-03,-0.31184002E-02,
         0.12744533E-04,-0.24973002E-03, 0.58777485E-03, 0.83024683E-03,
         0.86055454E-02, 0.53689033E-02, 0.86642557E-03, 0.59557818E-02,
        -0.10150120E-01, 0.14171154E-03, 0.36696805E-04, 0.28595523E-03,
         0.35101364E-03, 0.18746324E-02, 0.43312157E-04, 0.42457273E-02,
        -0.62205247E-02,-0.10255328E-01,-0.35652331E-04,-0.63654472E-03,
         0.10459892E-02,-0.98884932E-03,-0.16688217E-02,-0.12184107E-02,
        -0.59372336E-02,-0.65993539E-04,-0.53644623E-03,-0.10777795E-02,
        -0.33547406E-03, 0.62086285E-04,-0.17983117E-02, 0.21946947E-02,
        -0.36838310E-03,-0.32066109E-03,-0.84738953E-04,-0.18880509E-03,
        -0.20482924E-03,-0.15461126E-03,-0.28735501E-03,-0.81155093E-04,
        -0.19330223E-03, 0.18285535E-02,
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
        +coeff[  3]            *x41
        +coeff[  4]    *x21    *x41
        +coeff[  5]    *x21        *x51
        +coeff[  6]*x12*x21*x31
        +coeff[  7]    *x23*x31
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[  8]*x13        *x41
        +coeff[  9]*x13*x22*x31
        +coeff[ 10]        *x31
        +coeff[ 11]*x11*x21
        +coeff[ 12]    *x21*x31
        +coeff[ 13]*x11        *x41
        +coeff[ 14]*x11            *x51
        +coeff[ 15]    *x23
        +coeff[ 16]    *x21    *x42
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 17]*x13*x22
        +coeff[ 18]*x12*x21*x31*x41
        +coeff[ 19]    *x23*x31*x41
        +coeff[ 20]                *x51
        +coeff[ 21]*x11    *x31
        +coeff[ 22]*x13
        +coeff[ 23]*x11*x22
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]    *x23    *x41
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 26]*x13*x22    *x41
        +coeff[ 27]    *x22
        +coeff[ 28]*x12*x21
        +coeff[ 29]    *x21*x32
        +coeff[ 30]*x11        *x42
        +coeff[ 31]*x12*x21    *x41
        +coeff[ 32]*x11*x22    *x41
        +coeff[ 33]*x13    *x31*x41
        +coeff[ 34]    *x22    *x41
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]    *x21        *x52
        +coeff[ 37]*x13    *x31
        +coeff[ 38]*x11*x22*x31
        +coeff[ 39]    *x21    *x43
        +coeff[ 40]    *x23*x31*x42
        +coeff[ 41]            *x42
        +coeff[ 42]            *x41*x51
        +coeff[ 43]*x11    *x32
    ;
    v_t_l5p77_q2ex                            =v_t_l5p77_q2ex
        +coeff[ 44]    *x22        *x51
        +coeff[ 45]*x11        *x41*x51
        +coeff[ 46]    *x21    *x41*x51
        +coeff[ 47]*x12    *x31*x41
        +coeff[ 48]        *x33*x41
        +coeff[ 49]    *x21*x31*x42
        ;

    return v_t_l5p77_q2ex                            ;
}
float y_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.2170862E-01;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.12681827E-01, 0.21400788E+00, 0.54143686E-02, 0.27470741E-01,
        -0.27022172E-01,-0.16939441E-01, 0.42321868E-01, 0.18260771E-01,
        -0.28071335E-01,-0.11897201E-01, 0.57013682E-02,-0.10120691E-01,
        -0.12761756E-01, 0.11882407E-02,-0.21295608E-02, 0.20845993E-02,
        -0.47302828E-02,-0.17575145E-02, 0.17103273E-01, 0.12181210E-01,
         0.68651028E-02,-0.66679302E-02, 0.48876265E-02,-0.56040818E-02,
        -0.50545903E-03, 0.38146184E-03, 0.41813526E-03,-0.15407923E-02,
        -0.18285610E-02,-0.36689281E-03,-0.12431979E-03, 0.25366005E-02,
         0.63529838E-03,-0.52149949E-03,-0.59731322E-03, 0.20223458E-02,
        -0.50818594E-03, 0.89876616E-03,-0.27776250E-03, 0.37720893E-04,
        -0.85610326E-03,-0.55301812E-03, 0.50908537E-03, 0.34981765E-03,
        -0.33993351E-02, 0.26202807E-02, 0.43656206E-03,-0.93724339E-04,
         0.12264134E-03, 0.26877705E-03,-0.17577127E-03,-0.33513055E-03,
        -0.72017318E-03, 0.73963549E-03, 0.53561089E-03,-0.38076155E-02,
         0.20999184E-02,-0.27545777E-03, 0.11440574E-02, 0.10736848E-02,
         0.14595393E-03, 0.83972333E-03, 0.41005839E-03, 0.97544561E-03,
         0.18813750E-02, 0.91921201E-03,-0.26174387E-03, 0.13916792E-02,
         0.10351740E-02, 0.69369831E-04,-0.71714654E-04,-0.93295421E-04,
         0.54662043E-04,-0.18370504E-03, 0.85287087E-04,-0.12379807E-02,
        -0.10496870E-02,-0.19017019E-03, 0.59789489E-03,-0.83181326E-03,
        -0.45650697E-03, 0.62866166E-03,-0.94659376E-03,-0.19935612E-02,
        -0.14318840E-02,-0.69480372E-03,-0.11607055E-04,-0.92084971E-04,
        -0.74384101E-04, 0.39226248E-04, 0.53656305E-03, 0.20376263E-04,
         0.40971936E-04, 0.10766292E-03, 0.83157276E-04,-0.20505458E-04,
         0.74223644E-04,
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

    float v_y_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]            *x42
        +coeff[  6]    *x22
        +coeff[  7]*x11*x21
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[  8]    *x22    *x41
        +coeff[  9]        *x31*x41
        +coeff[ 10]            *x41*x51
        +coeff[ 11]    *x22*x31
        +coeff[ 12]*x11*x21    *x41
        +coeff[ 13]*x11
        +coeff[ 14]        *x32
        +coeff[ 15]*x12
        +coeff[ 16]*x11*x21*x31
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 17]    *x21    *x41
        +coeff[ 18]    *x22    *x42
        +coeff[ 19]    *x22*x31*x41
        +coeff[ 20]*x11*x21    *x42
        +coeff[ 21]    *x24
        +coeff[ 22]*x11*x21*x31*x41
        +coeff[ 23]*x11*x23
        +coeff[ 24]    *x21*x31
        +coeff[ 25]        *x31    *x51
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 26]                *x52
        +coeff[ 27]*x12        *x41
        +coeff[ 28]*x12*x22
        +coeff[ 29]*x11        *x41
        +coeff[ 30]*x11    *x31
        +coeff[ 31]        *x31*x42
        +coeff[ 32]            *x42*x51
        +coeff[ 33]            *x41*x52
        +coeff[ 34]*x12    *x31
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 35]    *x22*x32
        +coeff[ 36]            *x45
        +coeff[ 37]*x11*x21*x32
        +coeff[ 38]        *x32*x43
        +coeff[ 39]        *x32*x45
        +coeff[ 40]    *x21    *x42
        +coeff[ 41]    *x21*x31*x41
        +coeff[ 42]    *x23
        +coeff[ 43]*x11*x22
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 44]    *x22    *x43
        +coeff[ 45]*x11*x23    *x41
        +coeff[ 46]    *x24    *x43
        +coeff[ 47]    *x21*x32
        +coeff[ 48]        *x33
        +coeff[ 49]        *x31*x41*x51
        +coeff[ 50]*x11*x21        *x51
        +coeff[ 51]    *x23    *x41
        +coeff[ 52]        *x31*x44
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 53]*x12        *x42
        +coeff[ 54]*x12    *x31*x41
        +coeff[ 55]    *x22*x31*x42
        +coeff[ 56]    *x24    *x41
        +coeff[ 57]*x13*x21
        +coeff[ 58]    *x24*x31
        +coeff[ 59]*x11*x23*x31
        +coeff[ 60]    *x24        *x51
        +coeff[ 61]*x12*x22    *x41
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 62]*x12*x22*x31
        +coeff[ 63]    *x22*x31*x44
        +coeff[ 64]            *x43
        +coeff[ 65]        *x32*x41
        +coeff[ 66]    *x22        *x51
        +coeff[ 67]        *x31*x43
        +coeff[ 68]        *x32*x42
        +coeff[ 69]*x12*x21
        +coeff[ 70]        *x31    *x52
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 71]    *x23*x31
        +coeff[ 72]                *x53
        +coeff[ 73]*x11*x22    *x41
        +coeff[ 74]*x12    *x32
        +coeff[ 75]    *x22*x32*x41
        +coeff[ 76]*x11*x21*x31*x42
        +coeff[ 77]    *x22*x33
        +coeff[ 78]    *x21    *x45
        +coeff[ 79]        *x31*x45
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 80]*x11*x21*x32*x41
        +coeff[ 81]    *x21*x31*x44
        +coeff[ 82]        *x32*x44
        +coeff[ 83]    *x24    *x42
        +coeff[ 84]    *x24*x31*x41
        +coeff[ 85]*x11*x21    *x45
        +coeff[ 86]    *x21        *x51
        +coeff[ 87]*x11        *x42
        +coeff[ 88]*x11    *x31*x41
    ;
    v_y_l5p77_q2ex                            =v_y_l5p77_q2ex
        +coeff[ 89]    *x21    *x41*x51
        +coeff[ 90]            *x44
        +coeff[ 91]    *x21*x31    *x51
        +coeff[ 92]        *x32    *x51
        +coeff[ 93]    *x21*x32*x41
        +coeff[ 94]        *x33*x41
        +coeff[ 95]*x12            *x51
        +coeff[ 96]    *x22    *x41*x51
        ;

    return v_y_l5p77_q2ex                            ;
}
float p_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 54;
    float avdat= -0.1879460E-02;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 55]={
         0.17982154E-02,-0.64724875E-02,-0.31554405E-01, 0.25422145E-02,
        -0.39336733E-02, 0.15522774E-02, 0.75518561E-03, 0.50959541E-02,
        -0.17655463E-02,-0.46968870E-03, 0.94052387E-03,-0.73491904E-03,
         0.11198396E-02,-0.11019913E-03, 0.20901990E-03, 0.10939458E-02,
         0.22524886E-02, 0.12591619E-02, 0.43046440E-03,-0.12796177E-03,
         0.18630516E-03, 0.95516630E-03,-0.80601935E-03,-0.22124514E-03,
         0.44100778E-03,-0.33785685E-03,-0.22319593E-02,-0.11258398E-03,
        -0.10690645E-03,-0.99023775E-04,-0.44548005E-03,-0.13592542E-03,
        -0.47698623E-03, 0.46666447E-03, 0.10395477E-03,-0.72663761E-03,
        -0.99648778E-05, 0.48397629E-04, 0.20116808E-03, 0.43430631E-04,
        -0.11135263E-03,-0.97876786E-04, 0.42477922E-03,-0.60622013E-04,
        -0.14125537E-03,-0.11800629E-02, 0.72312498E-04,-0.13149317E-03,
        -0.46331962E-03,-0.41586376E-03, 0.20084729E-03, 0.54921034E-04,
        -0.22713297E-03, 0.17992995E-03,
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
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]            *x42
        +coeff[  6]        *x31    *x51
        +coeff[  7]            *x41*x51
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[  8]*x11*x21
        +coeff[  9]    *x21
        +coeff[ 10]    *x21    *x41
        +coeff[ 11]                *x52
        +coeff[ 12]    *x22        *x51
        +coeff[ 13]*x11
        +coeff[ 14]    *x21*x31
        +coeff[ 15]        *x31*x41
        +coeff[ 16]    *x22    *x41
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 17]*x11*x21    *x41
        +coeff[ 18]*x11*x21        *x51
        +coeff[ 19]    *x24*x31
        +coeff[ 20]        *x32
        +coeff[ 21]    *x22*x31
        +coeff[ 22]            *x43
        +coeff[ 23]*x12
        +coeff[ 24]*x11*x21*x31
        +coeff[ 25]            *x41*x52
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 26]    *x22    *x42
        +coeff[ 27]    *x24*x31*x41
        +coeff[ 28]*x11        *x41
        +coeff[ 29]        *x32*x41
        +coeff[ 30]        *x31*x42
        +coeff[ 31]        *x31*x41*x51
        +coeff[ 32]    *x22    *x41*x51
        +coeff[ 33]*x11*x23
        +coeff[ 34]*x12        *x41
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 35]*x11*x21    *x42
        +coeff[ 36]*x11*x23*x31*x41
        +coeff[ 37]    *x21        *x51
        +coeff[ 38]    *x21    *x42
        +coeff[ 39]*x11            *x51
        +coeff[ 40]    *x21    *x41*x51
        +coeff[ 41]            *x42*x51
        +coeff[ 42]    *x24
        +coeff[ 43]        *x31    *x52
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 44]    *x22*x32
        +coeff[ 45]    *x22*x31*x41
        +coeff[ 46]                *x53
        +coeff[ 47]    *x22*x31    *x51
        +coeff[ 48]    *x24    *x41
        +coeff[ 49]*x11*x21*x31*x41
        +coeff[ 50]    *x23*x31*x41
        +coeff[ 51]*x12            *x51
        +coeff[ 52]*x11*x21    *x41*x51
    ;
    v_p_l5p77_q2ex                            =v_p_l5p77_q2ex
        +coeff[ 53]*x12*x22
        ;

    return v_p_l5p77_q2ex                            ;
}
float l_l5p77_q2ex                            (float *x,int m){
    int ncoeff= 90;
    float avdat= -0.2944912E-02;
    float xmin[10]={
        -0.14896E-01,-0.45642E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.26779E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.24858036E-02, 0.17680849E-02, 0.52641598E-02,-0.44883317E-02,
        -0.60625527E-04,-0.12929135E-03,-0.42068830E-03,-0.18421689E-02,
        -0.79897493E-02, 0.22039158E-03, 0.19014189E-02, 0.99985686E-03,
        -0.12405177E-03,-0.12057694E-02,-0.53490442E-02, 0.41603966E-03,
        -0.13300639E-02, 0.24920097E-02, 0.17458177E-03,-0.66410212E-04,
        -0.13059151E-03,-0.12604406E-03, 0.11422273E-02,-0.83169213E-03,
         0.13047387E-02,-0.52328245E-03, 0.70043618E-03,-0.71146519E-05,
         0.10777087E-04, 0.90608047E-03,-0.19609531E-03,-0.17838590E-03,
         0.16868010E-03, 0.15107588E-02, 0.84779633E-03,-0.27720289E-04,
         0.34044297E-04, 0.22675455E-04, 0.22403598E-03, 0.59264814E-04,
         0.98397497E-04, 0.10823738E-03, 0.14142971E-03, 0.32251107E-03,
         0.39265311E-03,-0.11917639E-02,-0.21882859E-04, 0.12104004E-04,
        -0.16981927E-04, 0.11211408E-04,-0.29911607E-04,-0.39769729E-04,
         0.26610629E-04, 0.34605026E-04,-0.10073531E-03,-0.22843790E-03,
        -0.53959189E-04,-0.95648415E-04,-0.87362336E-03, 0.19275055E-03,
        -0.35186231E-03,-0.34851382E-04, 0.36176461E-04, 0.19443436E-04,
        -0.21741515E-04,-0.25077674E-04, 0.77429002E-04,-0.25803328E-03,
        -0.70226815E-04,-0.29305937E-04,-0.31624702E-04,-0.35506670E-04,
        -0.30680494E-04,-0.19575348E-03,-0.15978319E-04,-0.38183178E-04,
        -0.22295000E-03, 0.11419816E-03,-0.20132848E-04,-0.61504106E-04,
        -0.61624356E-04,-0.56539488E-05,-0.68738309E-05,-0.29398820E-04,
         0.50684041E-04,-0.82335835E-04, 0.14635933E-04,-0.22698458E-04,
        -0.21490368E-04,-0.43171342E-04,
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
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_l_l5p77_q2ex                            =avdat
        +coeff[  0]
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]    *x22
        +coeff[  4]    *x21*x31
        +coeff[  5]        *x32
        +coeff[  6]    *x21    *x41
        +coeff[  7]        *x31*x41
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[  8]            *x42
        +coeff[  9]        *x31    *x51
        +coeff[ 10]            *x41*x51
        +coeff[ 11]*x11*x21
        +coeff[ 12]    *x23
        +coeff[ 13]    *x22*x31
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]            *x42*x51
        +coeff[ 16]*x11*x21    *x41
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 17]    *x22    *x42
        +coeff[ 18]    *x21
        +coeff[ 19]*x11
        +coeff[ 20]                *x52
        +coeff[ 21]*x12
        +coeff[ 22]            *x43
        +coeff[ 23]    *x24
        +coeff[ 24]    *x22*x31*x41
        +coeff[ 25]*x11*x23
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 26]*x11*x21    *x42
        +coeff[ 27]        *x33*x42
        +coeff[ 28]        *x33
        +coeff[ 29]        *x31*x42
        +coeff[ 30]            *x41*x52
        +coeff[ 31]*x11*x21*x31
        +coeff[ 32]*x11*x21        *x51
        +coeff[ 33]    *x24    *x41
        +coeff[ 34]*x11*x23    *x41
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 35]                *x51
        +coeff[ 36]    *x21        *x51
        +coeff[ 37]*x11            *x51
        +coeff[ 38]        *x32*x41
        +coeff[ 39]    *x21    *x42
        +coeff[ 40]    *x22        *x51
        +coeff[ 41]        *x31*x41*x51
        +coeff[ 42]    *x22*x32
        +coeff[ 43]*x11*x21*x31*x41
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 44]    *x24*x31
        +coeff[ 45]    *x22    *x43
        +coeff[ 46]    *x24*x31*x42
        +coeff[ 47]*x11    *x31
        +coeff[ 48]*x11        *x41
        +coeff[ 49]    *x21*x31*x41
        +coeff[ 50]        *x31    *x52
        +coeff[ 51]*x11*x22
        +coeff[ 52]*x12    *x31
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 53]*x12            *x51
        +coeff[ 54]    *x23    *x41
        +coeff[ 55]            *x44
        +coeff[ 56]    *x22*x31    *x51
        +coeff[ 57]*x12*x22
        +coeff[ 58]    *x22*x31*x42
        +coeff[ 59]*x11*x23*x31
        +coeff[ 60]*x11*x21    *x43
        +coeff[ 61]*x12        *x43
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 62]        *x33*x43
        +coeff[ 63]                *x53
        +coeff[ 64]    *x23*x31
        +coeff[ 65]        *x33*x41
        +coeff[ 66]    *x21    *x43
        +coeff[ 67]        *x31*x43
        +coeff[ 68]            *x43*x51
        +coeff[ 69]            *x42*x52
        +coeff[ 70]*x11*x22    *x41
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 71]*x11*x21*x31    *x51
        +coeff[ 72]*x11*x21    *x41*x51
        +coeff[ 73]    *x22*x32*x41
        +coeff[ 74]        *x32*x43
        +coeff[ 75]*x11*x21*x32*x41
        +coeff[ 76]*x11*x21*x31*x42
        +coeff[ 77]*x12*x22    *x41
        +coeff[ 78]        *x33*x42*x51
        +coeff[ 79]    *x24        *x52
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 80]*x11*x23        *x52
        +coeff[ 81]    *x21*x31    *x51
        +coeff[ 82]    *x21        *x52
        +coeff[ 83]*x12        *x41
        +coeff[ 84]    *x21*x31*x42
        +coeff[ 85]        *x32*x42
        +coeff[ 86]    *x23        *x51
        +coeff[ 87]    *x22    *x41*x51
        +coeff[ 88]        *x32*x41*x51
    ;
    v_l_l5p77_q2ex                            =v_l_l5p77_q2ex
        +coeff[ 89]        *x31*x42*x51
        ;

    return v_l_l5p77_q2ex                            ;
}
float x_l5p77_q3en                            (float *x,int m){
    int ncoeff= 60;
    float avdat= -0.7118862E-02;
    float xmin[10]={
        -0.14896E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
         0.23707168E-02,-0.35852060E-01, 0.13909730E+00, 0.16397388E+00,
         0.68509295E-02, 0.27263606E-01, 0.15980406E-01, 0.74595012E-01,
        -0.13559041E-03,-0.15223293E-02, 0.21631886E-03,-0.96889725E-02,
         0.24494126E-02,-0.59482716E-02, 0.59182295E-02, 0.16130695E-01,
         0.27626589E-01,-0.58269948E-02, 0.20191226E-01,-0.31397849E-01,
         0.12925268E-01,-0.19294187E-01,-0.27755074E-01,-0.17123516E-02,
        -0.18471009E-02, 0.32617161E-02, 0.53758058E-02,-0.28569303E-02,
         0.47417139E-02,-0.17481508E-01, 0.50038165E-04, 0.14770703E-02,
        -0.43792762E-02, 0.21260686E-02, 0.13414639E-02, 0.22475270E-02,
        -0.32001447E-02,-0.57611475E-02, 0.16103239E-04, 0.37715174E-03,
        -0.27822590E-03, 0.56867930E-03,-0.45200737E-03,-0.33080070E-02,
        -0.61913603E-03, 0.17994574E-02, 0.35785104E-02,-0.53303109E-04,
         0.29522546E-02, 0.16100568E-02, 0.45007570E-02, 0.23741086E-03,
         0.28023799E-03, 0.31024317E-03,-0.29276649E-03, 0.19962188E-03,
        -0.30923402E-03,-0.53778029E-03, 0.60502591E-03, 0.68856147E-03,
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
        +coeff[  4]            *x41
        +coeff[  5]    *x21        *x51
        +coeff[  6]    *x22
        +coeff[  7]    *x21    *x41
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[  8]*x13        *x41
        +coeff[  9]*x12*x21*x31
        +coeff[ 10]    *x21*x31    *x52
        +coeff[ 11]    *x23*x31
        +coeff[ 12]        *x31
        +coeff[ 13]                *x52
        +coeff[ 14]*x11    *x31
        +coeff[ 15]*x11        *x41
        +coeff[ 16]    *x21*x31
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 17]            *x42
        +coeff[ 18]    *x23
        +coeff[ 19]    *x21    *x42
        +coeff[ 20]*x11*x22
        +coeff[ 21]    *x21*x31*x41
        +coeff[ 22]    *x23    *x41
        +coeff[ 23]*x11*x21
        +coeff[ 24]        *x31*x41
        +coeff[ 25]*x12*x21
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 26]    *x21    *x41*x51
        +coeff[ 27]    *x21*x32
        +coeff[ 28]    *x22    *x41
        +coeff[ 29]*x11*x22    *x41
        +coeff[ 30]    *x22*x33
        +coeff[ 31]            *x41*x51
        +coeff[ 32]*x11        *x42
        +coeff[ 33]    *x22        *x51
        +coeff[ 34]    *x21*x31    *x51
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 35]    *x22*x31
        +coeff[ 36]*x12*x21    *x41
        +coeff[ 37]*x11*x22*x31
        +coeff[ 38]*x13    *x31*x41
        +coeff[ 39]*x12
        +coeff[ 40]        *x32
        +coeff[ 41]*x11        *x41*x51
        +coeff[ 42]*x11*x21    *x41
        +coeff[ 43]*x11    *x31*x41
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 44]            *x42*x51
        +coeff[ 45]    *x23        *x51
        +coeff[ 46]    *x21    *x43
        +coeff[ 47]*x13    *x32
        +coeff[ 48]    *x21*x31*x43
        +coeff[ 49]*x12*x21*x31*x42
        +coeff[ 50]    *x23*x31*x42
        +coeff[ 51]        *x31    *x51
        +coeff[ 52]*x13
    ;
    v_x_l5p77_q3en                            =v_x_l5p77_q3en
        +coeff[ 53]                *x53
        +coeff[ 54]*x12        *x41
        +coeff[ 55]*x11    *x31    *x51
        +coeff[ 56]    *x21        *x52
        +coeff[ 57]*x11    *x32
        +coeff[ 58]*x11*x22        *x51
        +coeff[ 59]*x11*x23
        ;

    return v_x_l5p77_q3en                            ;
}
float t_l5p77_q3en                            (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.7813910E-03;
    float xmin[10]={
        -0.14896E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.77841186E-03, 0.90991827E-02,-0.61016120E-01, 0.22359289E-01,
        -0.43262946E-02,-0.23707325E-01,-0.41066133E-02, 0.11803888E-03,
         0.41144917E-03, 0.28785742E-02,-0.28602530E-04,-0.21274930E-02,
        -0.18524814E-02,-0.86755510E-02,-0.50456077E-02,-0.13999171E-02,
        -0.62628821E-02, 0.66738132E-04, 0.85072042E-02, 0.34893113E-04,
         0.92965965E-04,-0.78713015E-03, 0.68920880E-03, 0.54522301E-03,
        -0.58979762E-03,-0.73168521E-04,-0.41985903E-02,-0.17487147E-02,
         0.54968609E-02, 0.84517328E-02,-0.10523635E-02,-0.85751526E-03,
         0.83148747E-03, 0.52433126E-02,-0.12207310E-04,-0.17480619E-03,
         0.10306796E-02, 0.15310673E-02, 0.32257609E-03, 0.17467557E-02,
         0.11634573E-02,-0.89409086E-03, 0.15679369E-03, 0.82229380E-04,
         0.19766617E-03, 0.16725600E-03,-0.15231362E-03, 0.29500359E-03,
         0.17976816E-03, 0.66101126E-03,
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
        +coeff[  3]                *x51
        +coeff[  4]    *x22
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x21        *x51
        +coeff[  7]*x12        *x41
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[  8]*x12*x21*x31
        +coeff[  9]    *x23*x31
        +coeff[ 10]*x13        *x41
        +coeff[ 11]            *x41
        +coeff[ 12]*x11    *x31
        +coeff[ 13]    *x21*x31
        +coeff[ 14]*x11        *x41
        +coeff[ 15]                *x52
        +coeff[ 16]    *x23
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 17]*x12    *x31
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]*x13*x22
        +coeff[ 20]    *x23*x31*x41
        +coeff[ 21]        *x31
        +coeff[ 22]*x11*x21
        +coeff[ 23]            *x42
        +coeff[ 24]*x11            *x51
        +coeff[ 25]*x13
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 26]*x11*x22
        +coeff[ 27]    *x22    *x41
        +coeff[ 28]    *x21*x31*x41
        +coeff[ 29]    *x23    *x41
        +coeff[ 30]*x12*x21
        +coeff[ 31]    *x22*x31
        +coeff[ 32]    *x21*x32
        +coeff[ 33]*x11*x22    *x41
        +coeff[ 34]*x13        *x42
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 35]*x12
        +coeff[ 36]*x11    *x31*x41
        +coeff[ 37]*x11        *x42
        +coeff[ 38]            *x42*x51
        +coeff[ 39]*x11*x22*x31
        +coeff[ 40]*x12*x21    *x41
        +coeff[ 41]    *x21    *x43
        +coeff[ 42]        *x31*x41
        +coeff[ 43]        *x31    *x51
    ;
    v_t_l5p77_q3en                            =v_t_l5p77_q3en
        +coeff[ 44]*x11    *x32
        +coeff[ 45]*x11*x21    *x41
        +coeff[ 46]    *x22        *x51
        +coeff[ 47]    *x21        *x52
        +coeff[ 48]                *x53
        +coeff[ 49]    *x22    *x42
        ;

    return v_t_l5p77_q3en                            ;
}
float y_l5p77_q3en                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.2444004E-01;
    float xmin[10]={
        -0.14896E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.70955786E-02, 0.80954097E-01, 0.49266065E-02,-0.20432927E-01,
        -0.30641388E-01,-0.18433414E-01,-0.60126923E-01, 0.43401759E-01,
        -0.69807838E-02,-0.21288607E-02, 0.83093913E-02, 0.69217518E-01,
         0.19082643E-01, 0.99974982E-02, 0.98720109E-02,-0.70238514E-02,
        -0.53033657E-01,-0.12703007E-01,-0.12695294E-01,-0.47971536E-02,
         0.13902059E-01, 0.16569743E-03,-0.57339314E-02, 0.95233218E-05,
        -0.64184633E-03,-0.12556320E-01,-0.55752662E-02, 0.21922756E-02,
        -0.73410808E-02,-0.84267948E-02,-0.41123289E-02, 0.38315845E-02,
        -0.12770748E-01, 0.18148316E-02,-0.18777758E-01,-0.92990475E-03,
        -0.52905921E-02,-0.80660399E-03,-0.36582751E-02,-0.20449220E-02,
        -0.62209490E-03,-0.51862490E-02,-0.71870661E-02, 0.65874384E-03,
        -0.82902657E-03,-0.24634483E-02, 0.93195997E-02, 0.57345266E-02,
        -0.87921217E-03, 0.98820869E-03,-0.52797925E-02,-0.26419831E-02,
         0.85400203E-02,-0.17775672E-02,-0.59640250E-03, 0.61149604E-03,
         0.52498011E-02, 0.46983776E-02,-0.23082425E-02, 0.13410576E-02,
         0.16862889E-02, 0.64041559E-03, 0.42047640E-02,-0.77707466E-03,
         0.77914633E-03,-0.18724216E-02,-0.14141535E-02,-0.16783186E-02,
        -0.95347763E-03,-0.13492015E-02, 0.14978420E-02, 0.51528434E-02,
         0.56644427E-02, 0.21271773E-02, 0.19958052E-02,-0.32561193E-02,
         0.20884150E-02, 0.71614998E-03, 0.22756249E-03,-0.34281827E-03,
         0.32609250E-03, 0.29088819E-03, 0.11511892E-02, 0.89044578E-03,
         0.74528286E-03, 0.84168056E-03,-0.51045808E-03, 0.52602391E-03,
        -0.14430516E-02, 0.30486158E-02,-0.59345039E-03, 0.21479612E-02,
         0.87697018E-03,-0.14743858E-03, 0.10499596E-03,-0.12163113E-03,
         0.12811265E-03,
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

    float v_y_l5p77_q3en                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]            *x42
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x22
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[  8]    *x21*x31
        +coeff[  9]        *x32
        +coeff[ 10]*x11        *x41
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]        *x31    *x51
        +coeff[ 15]                *x52
        +coeff[ 16]    *x22    *x41
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 17]    *x23
        +coeff[ 18]    *x22*x31
        +coeff[ 19]    *x21    *x41*x51
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]        *x31*x43
        +coeff[ 22]            *x41*x52
        +coeff[ 23]        *x33*x41
        +coeff[ 24]    *x21    *x44
        +coeff[ 25]        *x31*x41
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 26]            *x43
        +coeff[ 27]*x12
        +coeff[ 28]    *x21*x31*x41
        +coeff[ 29]*x11*x21    *x41
        +coeff[ 30]*x11*x21*x31
        +coeff[ 31]*x11*x21        *x51
        +coeff[ 32]    *x24
        +coeff[ 33]*x11
        +coeff[ 34]    *x21    *x42
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 35]*x11            *x51
        +coeff[ 36]*x11        *x42
        +coeff[ 37]            *x42*x51
        +coeff[ 38]*x11*x22
        +coeff[ 39]*x12        *x41
        +coeff[ 40]    *x22    *x42
        +coeff[ 41]*x11*x22    *x41
        +coeff[ 42]*x11*x23
        +coeff[ 43]*x11    *x31
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 44]        *x31*x42
        +coeff[ 45]*x11    *x31*x41
        +coeff[ 46]    *x21    *x43
        +coeff[ 47]    *x21*x31*x42
        +coeff[ 48]        *x31    *x52
        +coeff[ 49]                *x53
        +coeff[ 50]    *x22    *x41*x51
        +coeff[ 51]*x11*x21    *x41*x51
        +coeff[ 52]    *x23    *x42
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 53]        *x31*x41*x51
        +coeff[ 54]*x12    *x31
        +coeff[ 55]    *x21        *x52
        +coeff[ 56]    *x22*x31*x41
        +coeff[ 57]*x11*x21    *x42
        +coeff[ 58]    *x21    *x42*x51
        +coeff[ 59]    *x23*x31
        +coeff[ 60]    *x22*x32
        +coeff[ 61]*x12            *x51
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 62]*x11*x21*x31*x41
        +coeff[ 63]    *x21*x31*x41*x51
        +coeff[ 64]*x11*x21*x32
        +coeff[ 65]    *x22*x31    *x51
        +coeff[ 66]*x12*x21    *x41
        +coeff[ 67]*x12*x22
        +coeff[ 68]*x11*x21*x31    *x51
        +coeff[ 69]    *x22        *x52
        +coeff[ 70]*x11*x21    *x43
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 71]*x11*x22    *x42
        +coeff[ 72]    *x22    *x42*x51
        +coeff[ 73]    *x22*x31*x41*x51
        +coeff[ 74]    *x23*x32*x41
        +coeff[ 75]    *x24    *x43
        +coeff[ 76]            *x45*x53
        +coeff[ 77]        *x32*x41
        +coeff[ 78]    *x21*x31    *x51
        +coeff[ 79]        *x32    *x51
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 80]*x11        *x41*x51
        +coeff[ 81]*x12*x21
        +coeff[ 82]*x11        *x43
        +coeff[ 83]*x11    *x31*x42
        +coeff[ 84]        *x31*x44
        +coeff[ 85]*x12        *x42
        +coeff[ 86]    *x23        *x51
        +coeff[ 87]*x12    *x31*x41
        +coeff[ 88]    *x24    *x41
    ;
    v_y_l5p77_q3en                            =v_y_l5p77_q3en
        +coeff[ 89]    *x23*x31*x41
        +coeff[ 90]*x11*x21        *x52
        +coeff[ 91]*x11*x22*x31*x41
        +coeff[ 92]        *x31*x42*x53
        +coeff[ 93]    *x21*x32
        +coeff[ 94]        *x33
        +coeff[ 95]*x11    *x32
        +coeff[ 96]*x11    *x31    *x51
        ;

    return v_y_l5p77_q3en                            ;
}
float p_l5p77_q3en                            (float *x,int m){
    int ncoeff= 90;
    float avdat=  0.3161908E-02;
    float xmin[10]={
        -0.14896E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
        -0.23281977E-04, 0.71653660E-03,-0.58645322E-02, 0.17414524E-02,
        -0.39901412E-02, 0.54968284E-02,-0.15241711E-02,-0.11612339E-01,
        -0.28576501E-03,-0.16656176E-02,-0.25035893E-02, 0.17156854E-02,
         0.16817837E-02, 0.13864052E-01, 0.25558537E-02,-0.17417731E-02,
        -0.18266230E-02, 0.16513511E-02,-0.94946846E-02,-0.39151199E-02,
         0.33674522E-02,-0.12896501E-02, 0.85107691E-03, 0.23526196E-03,
        -0.24644400E-02,-0.15749652E-02,-0.11704691E-02, 0.30422065E-03,
        -0.22289110E-02,-0.54817559E-03,-0.15991679E-02,-0.14876323E-02,
         0.38552946E-04, 0.24249156E-03,-0.10099507E-03,-0.40676942E-03,
        -0.69472345E-03, 0.23010485E-03,-0.90623094E-03,-0.51458640E-03,
        -0.11140850E-02,-0.33222823E-03,-0.10254388E-02, 0.81332073E-05,
        -0.51274325E-03,-0.13788475E-04,-0.44871822E-05, 0.20066729E-03,
         0.83136933E-04,-0.84544154E-04,-0.51913247E-03,-0.23433423E-03,
        -0.36633376E-03, 0.16209349E-03, 0.15476682E-02, 0.16248751E-02,
         0.15212147E-03,-0.55815658E-03, 0.38450054E-03,-0.22182862E-03,
        -0.26243244E-04, 0.40500172E-03, 0.28509196E-03,-0.88648027E-04,
        -0.16520577E-03,-0.26201949E-03,-0.28292146E-04, 0.18129693E-03,
         0.90121449E-03,-0.23980583E-03,-0.67380606E-04,-0.74424746E-03,
         0.42843280E-03, 0.41347698E-03, 0.35502741E-03,-0.22554035E-03,
         0.81184250E-03,-0.23725080E-03, 0.62334999E-04, 0.45413204E-03,
         0.69655705E-03, 0.39484943E-03, 0.18370226E-03,-0.11232050E-03,
        -0.10739580E-03, 0.38826474E-03,-0.41296027E-03, 0.97090298E-04,
        -0.33065495E-04, 0.15996407E-03,
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
    float x25 = x24*x2;
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
        +coeff[  8]        *x32
        +coeff[  9]        *x31*x41
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21        *x51
        +coeff[ 12]        *x31    *x51
        +coeff[ 13]            *x41*x51
        +coeff[ 14]*x11*x21
        +coeff[ 15]                *x52
        +coeff[ 16]    *x22*x31
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 17]*x11        *x41
        +coeff[ 18]    *x22    *x41
        +coeff[ 19]    *x21    *x42
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]*x11*x21        *x51
        +coeff[ 23]    *x25
        +coeff[ 24]    *x23
        +coeff[ 25]    *x21*x31*x41
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 26]            *x43
        +coeff[ 27]*x12
        +coeff[ 28]    *x24
        +coeff[ 29]            *x41*x52
        +coeff[ 30]    *x22    *x42
        +coeff[ 31]    *x22    *x41*x51
        +coeff[ 32]*x11*x23*x31
        +coeff[ 33]*x11    *x31
        +coeff[ 34]*x11            *x51
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 35]        *x31*x41*x51
        +coeff[ 36]*x11*x22
        +coeff[ 37]    *x21        *x52
        +coeff[ 38]*x11        *x42
        +coeff[ 39]    *x22*x31    *x51
        +coeff[ 40]*x11*x23
        +coeff[ 41]*x12        *x41
        +coeff[ 42]*x11*x22    *x41
        +coeff[ 43]*x13
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 44]*x11*x23    *x41
        +coeff[ 45]*x11    *x33*x41
        +coeff[ 46]*x13*x22
        +coeff[ 47]*x11
        +coeff[ 48]        *x32*x41
        +coeff[ 49]        *x32    *x51
        +coeff[ 50]*x11*x21*x31
        +coeff[ 51]*x11*x21    *x41
        +coeff[ 52]*x11    *x31*x41
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 53]                *x53
        +coeff[ 54]    *x21    *x43
        +coeff[ 55]    *x23    *x42
        +coeff[ 56]*x12            *x51
        +coeff[ 57]*x11*x21    *x41*x51
        +coeff[ 58]    *x23    *x41*x51
        +coeff[ 59]*x12*x21    *x41
        +coeff[ 60]    *x22*x31    *x52
        +coeff[ 61]    *x23*x31*x42
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 62]    *x23    *x43
        +coeff[ 63]    *x21*x32
        +coeff[ 64]        *x31*x42
        +coeff[ 65]            *x42*x51
        +coeff[ 66]*x11    *x32
        +coeff[ 67]    *x22*x32
        +coeff[ 68]    *x21*x31*x42
        +coeff[ 69]    *x21    *x42*x51
        +coeff[ 70]*x12    *x31
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 71]    *x24    *x41
        +coeff[ 72]*x11*x21*x31*x41
        +coeff[ 73]    *x23*x31*x41
        +coeff[ 74]*x11*x21    *x42
        +coeff[ 75]*x11*x21*x31    *x51
        +coeff[ 76]    *x22    *x42*x51
        +coeff[ 77]*x12*x22
        +coeff[ 78]    *x21*x33    *x51
        +coeff[ 79]    *x22    *x41*x52
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 80]*x11*x22    *x42
        +coeff[ 81]    *x23*x32*x41
        +coeff[ 82]*x12*x21    *x42
        +coeff[ 83]*x11*x22*x33
        +coeff[ 84]    *x21*x33    *x52
        +coeff[ 85]    *x24*x31*x41*x51
        +coeff[ 86]    *x22    *x45
        +coeff[ 87]    *x23*x31
        +coeff[ 88]        *x31    *x52
    ;
    v_p_l5p77_q3en                            =v_p_l5p77_q3en
        +coeff[ 89]    *x22*x31*x41
        ;

    return v_p_l5p77_q3en                            ;
}
float l_l5p77_q3en                            (float *x,int m){
    int ncoeff= 67;
    float avdat= -0.1414646E-01;
    float xmin[10]={
        -0.14896E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 68]={
         0.33203999E-02,-0.24372169E+00,-0.30230922E-02,-0.33368412E-01,
         0.41490786E-01,-0.22715457E-01,-0.99016473E-01,-0.16004261E-01,
        -0.21564076E-01, 0.11283266E-01,-0.11743783E-02,-0.35416599E-01,
        -0.78710457E-02,-0.27055467E-01, 0.37919965E-01,-0.17460516E-01,
        -0.18147826E-02, 0.37046173E-02,-0.28937408E-02,-0.15412476E-01,
         0.25519298E-01, 0.37149809E-01, 0.26807238E-04, 0.23186600E-01,
        -0.57672746E-02, 0.18271315E-02,-0.44227177E-02, 0.43020016E-02,
         0.62687281E-02,-0.43684407E-02, 0.75508789E-02,-0.13175332E-03,
         0.28364614E-03,-0.74464292E-03,-0.10165828E-02, 0.45286808E-02,
        -0.34403799E-02,-0.97514903E-02, 0.46776030E-02,-0.14761944E-02,
         0.18311052E-03, 0.11398190E-02, 0.95510256E-03, 0.76716422E-03,
        -0.12123750E-02,-0.37238435E-02, 0.46278480E-02,-0.98032681E-02,
        -0.21472615E-02, 0.15597002E-02, 0.14644226E-02,-0.79654570E-03,
        -0.20275594E-03, 0.45870978E-03, 0.96147886E-03, 0.40482014E-03,
        -0.32735339E-03, 0.47520106E-03, 0.21784454E-03,-0.32450087E-03,
         0.24363673E-02, 0.33763796E-02, 0.14141529E-02, 0.48019164E-02,
        -0.13033790E-02,-0.33633035E-03, 0.12821726E-02,
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
    float x54 = x53*x5;

//                 function

    float v_l_l5p77_q3en                            =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]*x11
        +coeff[  5]    *x22
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x21        *x51
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x23*x31
        +coeff[ 10]        *x31
        +coeff[ 11]    *x21*x31
        +coeff[ 12]*x11    *x31
        +coeff[ 13]    *x23
        +coeff[ 14]    *x21    *x42
        +coeff[ 15]*x11*x22
        +coeff[ 16]    *x23*x31*x41
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 17]*x11*x21
        +coeff[ 18]*x11            *x51
        +coeff[ 19]    *x22    *x41
        +coeff[ 20]    *x21*x31*x41
        +coeff[ 21]    *x23    *x41
        +coeff[ 22]            *x44
        +coeff[ 23]*x11*x22    *x41
        +coeff[ 24]            *x42
        +coeff[ 25]            *x41*x51
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 26]    *x22*x31
        +coeff[ 27]    *x21*x32
        +coeff[ 28]*x11        *x42
        +coeff[ 29]*x12*x21
        +coeff[ 30]*x11*x22*x31
        +coeff[ 31]*x11    *x33*x41
        +coeff[ 32]        *x31    *x51
        +coeff[ 33]*x12
        +coeff[ 34]    *x22        *x51
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 35]*x11    *x31*x41
        +coeff[ 36]    *x24
        +coeff[ 37]    *x21    *x43
        +coeff[ 38]*x12*x21    *x41
        +coeff[ 39]    *x21*x33*x42
        +coeff[ 40]        *x32
        +coeff[ 41]            *x42*x51
        +coeff[ 42]    *x21        *x52
        +coeff[ 43]*x11    *x32
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 44]*x11*x21    *x41
        +coeff[ 45]    *x21*x32*x41
        +coeff[ 46]    *x22    *x42
        +coeff[ 47]    *x21*x31*x42
        +coeff[ 48]*x11*x23
        +coeff[ 49]*x12*x21*x31
        +coeff[ 50]    *x23*x31    *x51
        +coeff[ 51]        *x31*x41
        +coeff[ 52]                *x52
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 53]            *x43
        +coeff[ 54]    *x21    *x41*x51
        +coeff[ 55]        *x31*x41*x51
        +coeff[ 56]*x11*x21*x31
        +coeff[ 57]*x11        *x41*x51
        +coeff[ 58]*x11            *x52
        +coeff[ 59]*x13
        +coeff[ 60]    *x22*x31*x41
        +coeff[ 61]    *x24    *x41
    ;
    v_l_l5p77_q3en                            =v_l_l5p77_q3en
        +coeff[ 62]    *x21*x31*x43
        +coeff[ 63]    *x21    *x44
        +coeff[ 64]    *x21    *x43*x51
        +coeff[ 65]            *x41*x54
        +coeff[ 66]*x11*x24*x31    *x51
        ;

    return v_l_l5p77_q3en                            ;
}
float x_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 51;
    float avdat=  0.1263800E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 52]={
         0.35508536E-02, 0.32765919E+00, 0.50983094E-01,-0.22299033E-01,
         0.28597495E-01, 0.40520322E-01, 0.11455074E-03,-0.27846839E-01,
         0.85493391E-02, 0.11646447E-01, 0.14847383E-01, 0.38004464E-02,
        -0.80836583E-02, 0.11139254E-01,-0.23190700E-01, 0.12878602E-02,
         0.31043862E-02,-0.24740342E-02, 0.26352585E-02, 0.64402279E-02,
         0.93004331E-02,-0.10273652E-01,-0.11362261E-03,-0.18065583E-01,
        -0.15067108E-02, 0.23916140E-02, 0.16752903E-02, 0.21800429E-02,
         0.19806230E-02,-0.15768469E-02, 0.16617653E-02,-0.10088260E-01,
        -0.47438033E-02, 0.69961016E-03,-0.18236373E-02,-0.16109756E-02,
         0.66769717E-03, 0.22769561E-02,-0.26557141E-02,-0.15323921E-03,
         0.25103865E-02, 0.24422293E-02,-0.76546293E-03, 0.77266450E-03,
        -0.65129192E-03,-0.42841097E-03,-0.15288236E-02, 0.18570543E-02,
        -0.16505859E-02, 0.55503106E-03,-0.19112455E-02,
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
        +coeff[  1]                *x51
        +coeff[  2]    *x21
        +coeff[  3]                *x52
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x21    *x41
        +coeff[  6]*x13
        +coeff[  7]*x11
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x22
        +coeff[ 10]    *x21*x31
        +coeff[ 11]            *x41
        +coeff[ 12]            *x42
        +coeff[ 13]    *x23
        +coeff[ 14]    *x21    *x42
        +coeff[ 15]        *x31
        +coeff[ 16]*x11    *x31
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 17]        *x31*x41
        +coeff[ 18]                *x53
        +coeff[ 19]*x11*x22
        +coeff[ 20]    *x21    *x41*x51
        +coeff[ 21]    *x21*x31*x41
        +coeff[ 22]            *x43*x51
        +coeff[ 23]    *x23    *x41
        +coeff[ 24]*x11            *x51
        +coeff[ 25]            *x41*x51
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 26]*x12*x21
        +coeff[ 27]    *x22        *x51
        +coeff[ 28]    *x21*x31    *x51
        +coeff[ 29]    *x21*x32
        +coeff[ 30]    *x22    *x41
        +coeff[ 31]*x11*x22    *x41
        +coeff[ 32]    *x23*x31
        +coeff[ 33]        *x31    *x51
        +coeff[ 34]*x11    *x31*x41
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 35]*x11        *x42
        +coeff[ 36]    *x22*x31
        +coeff[ 37]*x11*x23
        +coeff[ 38]*x11*x22*x31
        +coeff[ 39]*x11*x21    *x42
        +coeff[ 40]    *x23        *x51
        +coeff[ 41]    *x22    *x42
        +coeff[ 42]*x11*x21
        +coeff[ 43]*x11        *x41*x51
    ;
    v_x_l5p77_q3ex                            =v_x_l5p77_q3ex
        +coeff[ 44]    *x21        *x52
        +coeff[ 45]*x11    *x32
        +coeff[ 46]*x12*x21    *x41
        +coeff[ 47]*x11*x22        *x51
        +coeff[ 48]    *x21    *x42*x51
        +coeff[ 49]*x13    *x31    *x51
        +coeff[ 50]    *x23    *x41*x51
        ;

    return v_x_l5p77_q3ex                            ;
}
float t_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1246837E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.13608640E-02,-0.34646273E-02,-0.18788569E-01,-0.54402746E-04,
         0.58262685E-04, 0.11296885E+00, 0.53697024E-02,-0.10615812E-01,
         0.11454481E-02,-0.11060281E-02,-0.22518050E-02,-0.22147051E-02,
         0.21618523E-02,-0.78484538E-03, 0.84518746E-03,-0.69185277E-03,
        -0.31664295E-03, 0.63866854E-03, 0.11240422E-02, 0.25394205E-04,
        -0.17871444E-03,-0.27659346E-03,-0.36401671E-03,-0.73168788E-03,
         0.40247562E-03, 0.40985883E-03,-0.63140126E-03, 0.13638838E-02,
        -0.77996650E-04, 0.24726734E-03,-0.28028592E-03, 0.22916816E-03,
         0.20002149E-03, 0.16718329E-03, 0.41653728E-03,-0.68317848E-03,
         0.71123429E-03, 0.67644997E-03,-0.87192270E-03,-0.64880056E-04,
         0.99643170E-04,-0.60118640E-04, 0.10201860E-03, 0.89475609E-04,
         0.11619672E-03,-0.11535591E-03, 0.16234341E-03, 0.14805578E-03,
        -0.24682831E-03,-0.11562378E-03,
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
        +coeff[  9]    *x21    *x41
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]    *x21    *x41*x51
        +coeff[ 13]    *x21        *x52
        +coeff[ 14]*x11*x23
        +coeff[ 15]        *x31*x41
        +coeff[ 16]*x11            *x51
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 17]            *x41*x51
        +coeff[ 18]                *x53
        +coeff[ 19]        *x31    *x53
        +coeff[ 20]    *x21*x31
        +coeff[ 21]*x11        *x41
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x22    *x41
        +coeff[ 24]*x11        *x42
        +coeff[ 25]    *x21*x31    *x51
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 26]    *x23    *x41
        +coeff[ 27]    *x22    *x42
        +coeff[ 28]*x11    *x31
        +coeff[ 29]        *x31    *x51
        +coeff[ 30]    *x22*x31
        +coeff[ 31]    *x22        *x51
        +coeff[ 32]            *x42*x51
        +coeff[ 33]*x11            *x52
        +coeff[ 34]    *x23        *x51
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 35]    *x21    *x42*x51
        +coeff[ 36]    *x22    *x43
        +coeff[ 37]*x13*x22        *x51
        +coeff[ 38]    *x23*x31*x41*x51
        +coeff[ 39]*x12
        +coeff[ 40]*x11*x21
        +coeff[ 41]*x12*x21
        +coeff[ 42]        *x31*x42
        +coeff[ 43]*x11    *x31    *x51
    ;
    v_t_l5p77_q3ex                            =v_t_l5p77_q3ex
        +coeff[ 44]*x11        *x41*x51
        +coeff[ 45]            *x41*x52
        +coeff[ 46]*x12*x21    *x41
        +coeff[ 47]*x11*x21*x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x13            *x51
        ;

    return v_t_l5p77_q3ex                            ;
}
float y_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.1505348E-01;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
        -0.24318937E-02, 0.34582000E-01, 0.33583138E-02,-0.19851178E-01,
        -0.19005712E-01,-0.12210285E-01,-0.47466330E-01, 0.28051175E-01,
        -0.58326642E-02,-0.13808965E-02, 0.65096538E-02, 0.56042615E-01,
         0.12655817E-01, 0.73956405E-02, 0.61942372E-02,-0.15540482E-01,
        -0.71735797E-02,-0.40856589E-01,-0.94777802E-02,-0.85003218E-02,
        -0.66410648E-02, 0.13115576E-01, 0.31751842E-03, 0.36251980E-02,
        -0.82066930E-04, 0.31949385E-03,-0.82397284E-02,-0.47646752E-02,
         0.15418458E-02,-0.62629147E-02,-0.35826748E-02,-0.25792280E-02,
        -0.10014201E-01,-0.71032112E-02,-0.52006636E-03,-0.37112413E-02,
        -0.19466734E-02, 0.72124242E-02,-0.13703362E-02,-0.11702012E-02,
        -0.42522610E-02, 0.56042798E-04,-0.38546189E-02,-0.53644595E-02,
         0.27378183E-02, 0.91182947E-05,-0.98289871E-04,-0.28733615E-03,
        -0.58238226E-03, 0.10079734E-02, 0.77705825E-03,-0.75980835E-03,
        -0.14617818E-02,-0.26005434E-02, 0.48965233E-03, 0.46021636E-02,
         0.12017688E-02,-0.18152593E-02, 0.59248036E-03, 0.76869014E-03,
        -0.98715688E-03,-0.81481895E-03,-0.21581522E-02, 0.73852963E-02,
        -0.11906981E-02,-0.14444712E-02,-0.10133984E-03,-0.33437790E-03,
        -0.34901797E-03, 0.11080576E-02, 0.27540310E-02, 0.80926158E-03,
         0.10597454E-02, 0.23503779E-02, 0.43162212E-03,-0.18409149E-02,
        -0.83402242E-03, 0.23621221E-02,-0.99295401E-03, 0.41963425E-02,
        -0.59293414E-03, 0.21207742E-02, 0.83006558E-03, 0.67757751E-03,
         0.17985696E-02, 0.14077540E-02, 0.86803705E-03,-0.27261714E-02,
         0.86455047E-03, 0.27652533E-03,-0.29674204E-03, 0.57143945E-03,
         0.15928052E-03, 0.16865229E-03, 0.43669544E-03, 0.10000198E-02,
         0.14284424E-02,
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

    float v_y_l5p77_q3ex                            =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]                *x51
        +coeff[  5]            *x42
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x22
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[  8]    *x21*x31
        +coeff[  9]        *x32
        +coeff[ 10]*x11        *x41
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]        *x31    *x51
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]                *x52
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]    *x23
        +coeff[ 19]    *x22*x31
        +coeff[ 20]    *x21    *x41*x51
        +coeff[ 21]    *x22        *x51
        +coeff[ 22]        *x31*x43
        +coeff[ 23]*x11*x21        *x51
        +coeff[ 24]        *x33*x41
        +coeff[ 25]        *x33*x43
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 26]        *x31*x41
        +coeff[ 27]            *x43
        +coeff[ 28]*x12
        +coeff[ 29]    *x21*x31*x41
        +coeff[ 30]*x11*x21    *x41
        +coeff[ 31]*x11*x21*x31
        +coeff[ 32]    *x24
        +coeff[ 33]    *x22    *x41*x51
        +coeff[ 34]*x11            *x51
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 35]*x11        *x42
        +coeff[ 36]        *x31*x41*x51
        +coeff[ 37]    *x21    *x43
        +coeff[ 38]*x12        *x41
        +coeff[ 39]            *x41*x52
        +coeff[ 40]    *x22    *x42
        +coeff[ 41]*x13
        +coeff[ 42]*x11*x22    *x41
        +coeff[ 43]*x11*x23
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 44]*x11*x22    *x42
        +coeff[ 45]*x13    *x31
        +coeff[ 46]*x11    *x33*x41
        +coeff[ 47]*x11*x24
        +coeff[ 48]    *x24*x31    *x51
        +coeff[ 49]*x11
        +coeff[ 50]*x11    *x31
        +coeff[ 51]        *x31*x42
        +coeff[ 52]*x11    *x31*x41
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 53]*x11*x22
        +coeff[ 54]*x11        *x41*x51
        +coeff[ 55]    *x21*x31*x42
        +coeff[ 56]    *x21        *x52
        +coeff[ 57]    *x21    *x42*x51
        +coeff[ 58]*x12            *x51
        +coeff[ 59]                *x53
        +coeff[ 60]    *x23        *x51
        +coeff[ 61]*x12*x21    *x41
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 62]*x11*x21    *x41*x51
        +coeff[ 63]    *x23    *x42
        +coeff[ 64]*x12*x22
        +coeff[ 65]            *x42*x51
        +coeff[ 66]*x11    *x32
        +coeff[ 67]        *x32    *x51
        +coeff[ 68]*x12    *x31
        +coeff[ 69]    *x23    *x41
        +coeff[ 70]    *x22*x31*x41
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 71]    *x23*x31
        +coeff[ 72]    *x22*x32
        +coeff[ 73]*x11*x21*x31*x41
        +coeff[ 74]*x11*x21*x32
        +coeff[ 75]    *x22*x31    *x51
        +coeff[ 76]*x11*x21*x31    *x51
        +coeff[ 77]    *x23*x31*x41
        +coeff[ 78]            *x41*x53
        +coeff[ 79]    *x22    *x42*x51
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 80]*x11*x23    *x41
        +coeff[ 81]    *x22*x31*x41*x51
        +coeff[ 82]        *x31*x42*x52
        +coeff[ 83]*x11*x21    *x44
        +coeff[ 84]    *x24    *x42
        +coeff[ 85]*x11*x23    *x42
        +coeff[ 86]*x11*x21    *x45
        +coeff[ 87]    *x24    *x43
        +coeff[ 88]        *x32*x43*x52
    ;
    v_y_l5p77_q3ex                            =v_y_l5p77_q3ex
        +coeff[ 89]        *x32*x41
        +coeff[ 90]    *x21*x32
        +coeff[ 91]            *x44
        +coeff[ 92]*x12*x21
        +coeff[ 93]        *x31    *x52
        +coeff[ 94]*x11        *x43
        +coeff[ 95]    *x21*x32*x41
        +coeff[ 96]*x11*x21    *x42
        ;

    return v_y_l5p77_q3ex                            ;
}
float p_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 88;
    float avdat= -0.8725910E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 89]={
         0.31084111E-02,-0.17547710E-02, 0.60107457E-02,-0.32550678E-01,
         0.11130784E-01,-0.15729778E-01, 0.21491549E-02, 0.20090874E-01,
         0.65015019E-02,-0.32873510E-02,-0.36900502E-02,-0.20560391E-01,
        -0.68719201E-02, 0.13152147E-02,-0.27765601E-02, 0.17847329E-01,
         0.55311257E-02,-0.36079611E-02,-0.18697224E-02, 0.16796653E-03,
        -0.16786235E-03,-0.44927406E-03,-0.66494067E-04,-0.14498005E-03,
         0.78957481E-03, 0.44171745E-02, 0.48179985E-02,-0.86263171E-03,
         0.43185633E-02, 0.30867313E-02,-0.48966508E-03,-0.63896802E-03,
         0.43213353E-02, 0.19762651E-02, 0.36503316E-03,-0.39626792E-03,
         0.10974137E-02, 0.16869592E-02, 0.94435859E-03, 0.29716534E-02,
         0.14756710E-02,-0.88570087E-03, 0.25636246E-02, 0.12732509E-02,
        -0.12297677E-02, 0.23358315E-03, 0.12255072E-02, 0.66288887E-03,
        -0.25476434E-02, 0.13386421E-05, 0.62885397E-03, 0.13278543E-02,
         0.68063155E-03, 0.51030930E-03,-0.11810203E-02,-0.34573648E-04,
        -0.20765839E-03,-0.29129747E-03, 0.14757407E-02,-0.61898580E-03,
        -0.13508914E-03, 0.17667572E-03,-0.14412284E-02,-0.14215291E-02,
        -0.14764450E-03, 0.50019997E-03,-0.25398377E-03,-0.16710362E-02,
        -0.21322744E-03,-0.34421100E-03, 0.68040064E-03, 0.84123567E-04,
        -0.11004800E-02,-0.11076540E-03,-0.84862993E-04, 0.56856657E-04,
         0.46260750E-04,-0.11808337E-02, 0.13058676E-03,-0.11506011E-03,
         0.36264682E-03,-0.46644878E-03,-0.28073456E-03,-0.12637689E-02,
        -0.22124774E-03, 0.55284455E-03, 0.18259471E-03,-0.43779021E-03,
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
    float x25 = x24*x2;
    float x26 = x25*x2;
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
        +coeff[  8]            *x42
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
        +coeff[ 18]    *x22*x31*x41
        +coeff[ 19]        *x33*x41
        +coeff[ 20]    *x25
        +coeff[ 21]    *x24*x31
        +coeff[ 22]            *x43*x52
        +coeff[ 23]        *x35*x41
        +coeff[ 24]        *x32
        +coeff[ 25]        *x31*x41
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 26]    *x22*x31
        +coeff[ 27]*x12
        +coeff[ 28]    *x24
        +coeff[ 29]*x11*x21    *x41
        +coeff[ 30]*x11*x23*x31
        +coeff[ 31]*x11
        +coeff[ 32]    *x23
        +coeff[ 33]    *x21*x31*x41
        +coeff[ 34]*x11            *x51
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 35]            *x42*x51
        +coeff[ 36]*x11*x22
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]        *x31    *x52
        +coeff[ 39]            *x41*x52
        +coeff[ 40]*x11        *x42
        +coeff[ 41]*x11*x21        *x51
        +coeff[ 42]*x11*x23
        +coeff[ 43]    *x22    *x43
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 44]            *x41*x53
        +coeff[ 45]            *x45
        +coeff[ 46]*x11*x22    *x43
        +coeff[ 47]*x11    *x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x11*x22*x31
        +coeff[ 50]*x12        *x41
        +coeff[ 51]*x11*x22    *x41
        +coeff[ 52]*x12*x22
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 53]*x12*x21    *x41
        +coeff[ 54]    *x23*x31*x42
        +coeff[ 55]*x13    *x31
        +coeff[ 56]*x11    *x31
        +coeff[ 57]        *x32*x41
        +coeff[ 58]            *x43
        +coeff[ 59]    *x22*x32
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]*x12    *x31
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 62]*x11*x21*x31*x41
        +coeff[ 63]*x11*x21    *x42
        +coeff[ 64]*x12            *x51
        +coeff[ 65]*x11*x21    *x41*x51
        +coeff[ 66]        *x31    *x53
        +coeff[ 67]    *x22    *x42*x51
        +coeff[ 68]    *x25        *x51
        +coeff[ 69]*x11*x22        *x52
        +coeff[ 70]    *x26    *x41
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 71]        *x31*x44*x51
        +coeff[ 72]    *x25*x32*x41
        +coeff[ 73]    *x21*x31    *x51
        +coeff[ 74]    *x21    *x41*x51
        +coeff[ 75]        *x32    *x51
        +coeff[ 76]*x11    *x32
        +coeff[ 77]    *x21*x31*x42
        +coeff[ 78]*x11        *x41*x51
        +coeff[ 79]*x12*x21
    ;
    v_p_l5p77_q3ex                            =v_p_l5p77_q3ex
        +coeff[ 80]    *x22        *x52
        +coeff[ 81]            *x43*x51
        +coeff[ 82]*x11*x21*x32
        +coeff[ 83]    *x23    *x42
        +coeff[ 84]*x11        *x43
        +coeff[ 85]    *x21    *x44
        +coeff[ 86]    *x23        *x52
        +coeff[ 87]    *x25*x31
        ;

    return v_p_l5p77_q3ex                            ;
}
float l_l5p77_q3ex                            (float *x,int m){
    int ncoeff= 67;
    float avdat= -0.1433178E-01;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 68]={
         0.26995754E-02,-0.24379745E+00,-0.31171765E-02,-0.32926898E-01,
         0.41462358E-01,-0.25345350E-01,-0.99193439E-01,-0.21686895E-01,
         0.11527403E-01,-0.35253096E-01,-0.98069226E-02,-0.83938567E-02,
         0.43195281E-02,-0.79289488E-02,-0.27324988E-01, 0.26424617E-04,
        -0.17502174E-01, 0.38537458E-01,-0.17075201E-01,-0.19945044E-02,
        -0.61267796E-02,-0.28843130E-02,-0.49089766E-02, 0.25535718E-01,
        -0.43953098E-02, 0.37359435E-01, 0.23372313E-01, 0.22435079E-02,
         0.42520417E-02, 0.66923085E-02, 0.76416740E-02,-0.14799997E-02,
        -0.23675174E-03,-0.11001963E-02,-0.79833175E-03, 0.14482556E-02,
        -0.14734596E-02, 0.48777158E-02,-0.39548851E-02,-0.93346499E-02,
         0.45333435E-02,-0.26959695E-02, 0.33647160E-03,-0.91004750E-03,
         0.11847747E-02, 0.16399116E-02, 0.68457995E-03, 0.56879583E-03,
         0.48192721E-02,-0.93724634E-02,-0.88525692E-03,-0.26643937E-02,
         0.14531041E-02, 0.34102059E-02,-0.14473068E-02, 0.13834273E-03,
         0.64315839E-03, 0.40845005E-03,-0.30813989E-03,-0.52073027E-03,
         0.21160184E-02,-0.29609345E-02,-0.62908645E-03,-0.61075186E-03,
         0.40179370E-02,-0.73171675E-03,-0.22333269E-02,
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
        +coeff[  7]*x11        *x41
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[  8]    *x23*x31
        +coeff[  9]    *x21*x31
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]                *x52
        +coeff[ 12]*x11*x21
        +coeff[ 13]*x11    *x31
        +coeff[ 14]    *x23
        +coeff[ 15]        *x33
        +coeff[ 16]    *x22    *x41
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]*x11*x22
        +coeff[ 19]    *x23*x31*x41
        +coeff[ 20]            *x42
        +coeff[ 21]*x11            *x51
        +coeff[ 22]    *x22*x31
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]*x12*x21
        +coeff[ 25]    *x23    *x41
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 26]*x11*x22    *x41
        +coeff[ 27]            *x41*x51
        +coeff[ 28]    *x21*x32
        +coeff[ 29]*x11        *x42
        +coeff[ 30]*x11*x22*x31
        +coeff[ 31]*x11*x22*x31*x41
        +coeff[ 32]*x11    *x33*x41
        +coeff[ 33]        *x31
        +coeff[ 34]*x12
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 35]                *x53
        +coeff[ 36]*x11*x21    *x41
        +coeff[ 37]*x11    *x31*x41
        +coeff[ 38]    *x24
        +coeff[ 39]    *x21    *x43
        +coeff[ 40]*x12*x21    *x41
        +coeff[ 41]    *x23*x31*x42
        +coeff[ 42]        *x31    *x51
        +coeff[ 43]    *x22        *x51
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 44]    *x21*x31    *x51
        +coeff[ 45]    *x21    *x41*x51
        +coeff[ 46]*x11    *x32
        +coeff[ 47]*x11        *x41*x51
        +coeff[ 48]    *x22    *x42
        +coeff[ 49]    *x21*x31*x42
        +coeff[ 50]    *x21*x31*x41*x51
        +coeff[ 51]*x11*x23
        +coeff[ 52]*x12*x21*x31
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 53]    *x21    *x44
        +coeff[ 54]    *x23*x32*x41
        +coeff[ 55]        *x32
        +coeff[ 56]            *x43
        +coeff[ 57]*x11            *x52
        +coeff[ 58]*x13
        +coeff[ 59]    *x21*x33
        +coeff[ 60]    *x22*x31*x41
        +coeff[ 61]    *x21*x32*x41
    ;
    v_l_l5p77_q3ex                            =v_l_l5p77_q3ex
        +coeff[ 62]        *x33*x41
        +coeff[ 63]*x11*x21    *x41*x51
        +coeff[ 64]    *x24    *x41
        +coeff[ 65]*x11*x23*x31
        +coeff[ 66]*x11*x22    *x42
        ;

    return v_l_l5p77_q3ex                            ;
}
float x_l5p77_sen                             (float *x,int m){
    int ncoeff=  6;
    float avdat=  0.1116779E+00;
    float xmin[10]={
        -0.14934E-01,-0.45984E-01,-0.14997E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.50023E-01, 0.14990E-01, 0.27012E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[  7]={
        -0.57888664E-02,-0.45730890E-05, 0.15074139E-01, 0.33173349E-01,
         0.96990203E-04, 0.42454321E-04,
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
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
//  set up monomials   functions
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x51 = x5;

//                 function

    float v_x_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]                *x51
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]            *x42
        +coeff[  5]        *x31*x41
        ;

    return v_x_l5p77_sen                             ;
}
float t_l5p77_sen                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1042097E+00;
    float xmin[10]={
        -0.14934E-01,-0.45984E-01,-0.14997E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.50023E-01, 0.14990E-01, 0.27012E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
        -0.40551065E-02, 0.28165974E-03, 0.31231135E-01,-0.12303058E-03,
         0.10585090E-03,-0.16085958E-03,-0.28927505E-03, 0.32048221E-04,
         0.77892320E-04, 0.68470428E-04,-0.13620673E-03,-0.17624490E-03,
         0.55618793E-05,-0.36521920E-04,-0.26010799E-04,-0.79315170E-04,
        -0.18811852E-04,-0.59623624E-04,-0.48355473E-05,-0.11996951E-04,
        -0.21214006E-04,-0.52019368E-04,-0.21562444E-05, 0.12347523E-03,
         0.63484775E-04, 0.10664433E-03, 0.37288148E-04,-0.60134362E-05,
         0.15545123E-04,-0.65146346E-05, 0.65366398E-05,-0.23832488E-04,
         0.23715418E-04, 0.56828881E-04, 0.19482113E-04,-0.48377275E-04,
        -0.36607128E-05,-0.83422301E-05,-0.17638994E-05,-0.95028827E-05,
         0.82670385E-05,-0.47732847E-05, 0.96378535E-05, 0.67120827E-05,
         0.79406564E-05, 0.67786550E-05, 0.12955617E-04, 0.22861755E-04,
        -0.73966266E-05, 0.11602300E-04,
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
        +coeff[  1]        *x31
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]    *x21
        +coeff[  5]        *x31*x41
        +coeff[  6]    *x22    *x41
        +coeff[  7]*x11
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[  8]*x11*x21
        +coeff[  9]    *x22
        +coeff[ 10]    *x22*x31
        +coeff[ 11]*x11*x21    *x41
        +coeff[ 12]*x13*x21*x31
        +coeff[ 13]        *x32
        +coeff[ 14]            *x41*x51
        +coeff[ 15]*x11*x21*x31
        +coeff[ 16]    *x21    *x41
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 17]            *x42
        +coeff[ 18]    *x21        *x51
        +coeff[ 19]        *x31    *x51
        +coeff[ 20]*x12        *x41
        +coeff[ 21]*x11*x23
        +coeff[ 22]    *x23*x31
        +coeff[ 23]    *x22*x31*x41
        +coeff[ 24]*x11*x21    *x42
        +coeff[ 25]    *x22    *x42
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 26]        *x31*x43
        +coeff[ 27]*x13*x21*x31*x41
        +coeff[ 28]*x12
        +coeff[ 29]*x11        *x41
        +coeff[ 30]                *x52
        +coeff[ 31]*x12*x22
        +coeff[ 32]    *x22*x32
        +coeff[ 33]*x11*x21*x31*x41
        +coeff[ 34]        *x32*x42
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 35]    *x22*x31*x43
        +coeff[ 36]*x11    *x31
        +coeff[ 37]    *x21*x31
        +coeff[ 38]*x11            *x51
        +coeff[ 39]*x12    *x31
        +coeff[ 40]            *x43
        +coeff[ 41]*x13*x21
        +coeff[ 42]*x11*x21*x32
        +coeff[ 43]*x11*x21    *x41*x51
    ;
    v_t_l5p77_sen                             =v_t_l5p77_sen
        +coeff[ 44]    *x22    *x41*x51
        +coeff[ 45]        *x31*x42*x51
        +coeff[ 46]*x12*x22    *x41
        +coeff[ 47]*x11*x23    *x41
        +coeff[ 48]    *x22*x32*x41
        +coeff[ 49]*x12    *x31*x42
        ;

    return v_t_l5p77_sen                             ;
}
float y_l5p77_sen                             (float *x,int m){
    int ncoeff= 11;
    float avdat= -0.2395591E-02;
    float xmin[10]={
        -0.14934E-01,-0.45984E-01,-0.14997E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.50023E-01, 0.14990E-01, 0.27012E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 12]={
         0.17931301E-03,-0.52200783E-01,-0.14957558E-01,-0.17749335E-03,
        -0.80072743E-04,-0.97305829E-05,-0.45018787E-05,-0.49424752E-05,
        -0.19665558E-05,-0.40813661E-05,-0.36329393E-05,
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
//  set up monomials   functions
    float x11 = x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x31 = x3;
    float x41 = x4;

//                 function

    float v_y_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]*x11
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x21*x31
        +coeff[  5]            *x41
        +coeff[  6]        *x31
        +coeff[  7]*x11        *x41
    ;
    v_y_l5p77_sen                             =v_y_l5p77_sen
        +coeff[  8]*x11    *x31
        +coeff[  9]    *x23
        +coeff[ 10]*x11*x22
        ;

    return v_y_l5p77_sen                             ;
}
float p_l5p77_sen                             (float *x,int m){
    int ncoeff= 10;
    float avdat= -0.1965474E-02;
    float xmin[10]={
        -0.14934E-01,-0.45984E-01,-0.14997E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.50023E-01, 0.14990E-01, 0.27012E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 11]={
         0.26183916E-03,-0.47432907E-01,-0.62035001E-03, 0.23146243E-03,
        -0.19282349E-03,-0.77415367E-04,-0.11682450E-03,-0.34024211E-04,
        -0.12898463E-03,-0.10973937E-03,
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
//  set up monomials   functions
    float x11 = x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x31 = x3;
    float x41 = x4;

//                 function

    float v_p_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]    *x21    *x41
        +coeff[  3]*x11
        +coeff[  4]    *x21*x31
        +coeff[  5]            *x41
        +coeff[  6]*x11        *x41
        +coeff[  7]        *x31
    ;
    v_p_l5p77_sen                             =v_p_l5p77_sen
        +coeff[  8]    *x23
        +coeff[  9]*x11*x22
        ;

    return v_p_l5p77_sen                             ;
}
float l_l5p77_sen                             (float *x,int m){
    int ncoeff=  8;
    float avdat= -0.7653352E-03;
    float xmin[10]={
        -0.14934E-01,-0.45984E-01,-0.14997E-01,-0.33638E-01,-0.49937E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.50023E-01, 0.14990E-01, 0.27012E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[  9]={
         0.11223612E-02,-0.10518804E-03,-0.15161977E-02,-0.32294181E-02,
        -0.12520733E-02,-0.48099678E-05,-0.51060133E-03,-0.34084824E-05,
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
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
//  set up monomials   functions
    float x21 = x2;
    float x22 = x21*x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_l_l5p77_sen                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]    *x22
        +coeff[  5]        *x31*x41
        +coeff[  6]            *x42
        +coeff[  7]    *x22    *x41
        ;

    return v_l_l5p77_sen                             ;
}
float x_l5p77_sex                             (float *x,int m){
    int ncoeff= 49;
    float avdat=  0.4850814E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 50]={
         0.73665543E-02,-0.37995819E-01, 0.65293777E+00,-0.52604344E-01,
         0.43976791E-01, 0.15137231E-01, 0.36930632E-01,-0.30657919E-02,
         0.77053253E-02, 0.14181478E-01,-0.14557523E-01,-0.29401777E-01,
         0.39109876E-02, 0.56102672E-02, 0.17066566E-01, 0.10962527E-01,
        -0.27189681E-02, 0.27488384E-02, 0.41870973E-02,-0.44432557E-02,
        -0.79075755E-04,-0.32041839E-02, 0.55091535E-02,-0.10267943E-01,
         0.39130161E-02,-0.20187652E-01,-0.13401572E-02, 0.16599483E-03,
         0.11631739E-02, 0.14370932E-02, 0.15681073E-02, 0.20424565E-02,
        -0.14076697E-02,-0.10330845E-01,-0.49683470E-02,-0.43462724E-02,
         0.63118450E-02, 0.64997573E-03,-0.14600570E-02, 0.39820308E-02,
         0.35143178E-02,-0.11653090E-02,-0.25365625E-02, 0.37741067E-02,
        -0.25085355E-02, 0.11700870E-02,-0.48278663E-02, 0.34487743E-02,
         0.16507468E-02,
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
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_x_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]                *x51
        +coeff[  3]                *x52
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21    *x41
        +coeff[  7]    *x21
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x21*x31
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]            *x41
        +coeff[ 13]                *x53
        +coeff[ 14]    *x21    *x41*x51
        +coeff[ 15]    *x23
        +coeff[ 16]*x11            *x51
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 17]*x11    *x31
        +coeff[ 18]            *x41*x51
        +coeff[ 19]        *x31*x41
        +coeff[ 20]*x12    *x31
        +coeff[ 21]    *x21        *x52
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x21*x31*x41
        +coeff[ 24]*x11*x23
        +coeff[ 25]    *x23    *x41
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 26]    *x23*x31    *x51
        +coeff[ 27]*x13*x22    *x41
        +coeff[ 28]        *x31
        +coeff[ 29]        *x31    *x51
        +coeff[ 30]*x12*x21
        +coeff[ 31]    *x22        *x51
        +coeff[ 32]    *x21*x32
        +coeff[ 33]*x11*x22    *x41
        +coeff[ 34]    *x21    *x42*x51
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 35]    *x23*x31
        +coeff[ 36]    *x22    *x42
        +coeff[ 37]*x11            *x52
        +coeff[ 38]*x11    *x31*x41
        +coeff[ 39]    *x21*x31    *x51
        +coeff[ 40]*x11*x22        *x51
        +coeff[ 41]    *x22        *x52
        +coeff[ 42]*x11*x22*x31
        +coeff[ 43]    *x23        *x51
    ;
    v_x_l5p77_sex                             =v_x_l5p77_sex
        +coeff[ 44]    *x21*x31*x41*x51
        +coeff[ 45]*x13        *x41*x51
        +coeff[ 46]    *x23    *x41*x51
        +coeff[ 47]    *x22    *x42*x51
        +coeff[ 48]*x13    *x31    *x53
        ;

    return v_x_l5p77_sex                             ;
}
float t_l5p77_sex                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1246783E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.13610454E-02,-0.34647479E-02,-0.18788554E-01,-0.54564876E-04,
         0.58289759E-04, 0.11297005E+00, 0.53696432E-02,-0.10615933E-01,
         0.11454147E-02,-0.11069913E-02,-0.22524882E-02,-0.22141344E-02,
         0.21624689E-02,-0.78605616E-03, 0.84544026E-03,-0.69176994E-03,
        -0.31655843E-03, 0.63878635E-03, 0.11204408E-02, 0.25421859E-04,
        -0.17859080E-03,-0.27667859E-03,-0.36402428E-03,-0.73233357E-03,
         0.40273461E-03, 0.40990915E-03,-0.63067366E-03, 0.13648443E-02,
        -0.78024612E-04, 0.24767272E-03,-0.28033348E-03, 0.22862288E-03,
         0.19827017E-03, 0.16813961E-03, 0.41615704E-03,-0.68232138E-03,
         0.71254856E-03, 0.67612348E-03,-0.87276765E-03,-0.64878201E-04,
         0.99597441E-04,-0.60150411E-04, 0.10238977E-03, 0.89657100E-04,
         0.11578119E-03,-0.11452480E-03, 0.16251726E-03, 0.14832940E-03,
        -0.24610086E-03,-0.11560744E-03,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]*x11
        +coeff[  2]    *x21
        +coeff[  3]        *x31
        +coeff[  4]            *x41
        +coeff[  5]                *x51
        +coeff[  6]    *x21        *x51
        +coeff[  7]                *x52
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[  8]    *x22
        +coeff[  9]    *x21    *x41
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]    *x21    *x41*x51
        +coeff[ 13]    *x21        *x52
        +coeff[ 14]*x11*x23
        +coeff[ 15]        *x31*x41
        +coeff[ 16]*x11            *x51
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 17]            *x41*x51
        +coeff[ 18]                *x53
        +coeff[ 19]        *x31    *x53
        +coeff[ 20]    *x21*x31
        +coeff[ 21]*x11        *x41
        +coeff[ 22]*x11*x22
        +coeff[ 23]    *x22    *x41
        +coeff[ 24]*x11        *x42
        +coeff[ 25]    *x21*x31    *x51
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 26]    *x23    *x41
        +coeff[ 27]    *x22    *x42
        +coeff[ 28]*x11    *x31
        +coeff[ 29]        *x31    *x51
        +coeff[ 30]    *x22*x31
        +coeff[ 31]    *x22        *x51
        +coeff[ 32]            *x42*x51
        +coeff[ 33]*x11            *x52
        +coeff[ 34]    *x23        *x51
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 35]    *x21    *x42*x51
        +coeff[ 36]    *x22    *x43
        +coeff[ 37]*x13*x22        *x51
        +coeff[ 38]    *x23*x31*x41*x51
        +coeff[ 39]*x12
        +coeff[ 40]*x11*x21
        +coeff[ 41]*x12*x21
        +coeff[ 42]        *x31*x42
        +coeff[ 43]*x11    *x31    *x51
    ;
    v_t_l5p77_sex                             =v_t_l5p77_sex
        +coeff[ 44]*x11        *x41*x51
        +coeff[ 45]            *x41*x52
        +coeff[ 46]*x12*x21    *x41
        +coeff[ 47]*x11*x21*x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x13            *x51
        ;

    return v_t_l5p77_sex                             ;
}
float y_l5p77_sex                             (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.1005435E-01;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
         0.65170373E-02,-0.59205592E-01, 0.13075430E-01, 0.11103833E-01,
        -0.17059630E-01,-0.15923406E-02,-0.32548006E-02,-0.71285306E-02,
        -0.17285218E-02,-0.46245498E-02,-0.33543792E-02, 0.10840210E-01,
        -0.68390802E-02,-0.68930436E-04,-0.64344858E-04, 0.69019529E-02,
        -0.71534398E-03, 0.52608508E-02, 0.24673267E-03, 0.13028190E-02,
        -0.25707176E-02, 0.64810598E-02, 0.45800009E-02, 0.78306231E-03,
         0.26954676E-02,-0.88454189E-03, 0.63687604E-03,-0.19460334E-02,
         0.45873541E-02, 0.16111672E-02, 0.27144498E-02, 0.78555220E-03,
        -0.43049683E-02, 0.10812604E-02, 0.29410238E-02,-0.28960968E-02,
        -0.16620336E-02,-0.56533527E-03, 0.20937492E-03,-0.46297577E-02,
         0.78232371E-03,-0.18276208E-02, 0.53455081E-03, 0.83134189E-03,
         0.28245812E-02,-0.11332378E-02, 0.54305681E-03,-0.21958493E-02,
        -0.14797725E-02, 0.26076841E-02,-0.16080920E-02,-0.62764520E-02,
         0.22293869E-02, 0.12296953E-02, 0.63777401E-03, 0.17202366E-03,
         0.43299180E-03,-0.19844405E-02, 0.44707116E-03,-0.11669641E-02,
        -0.30878978E-03,-0.15550818E-02,-0.43529461E-03, 0.61408174E-03,
        -0.11350361E-02, 0.60451211E-03, 0.12153829E-02,-0.22525168E-02,
        -0.67765207E-03,-0.84672478E-03, 0.23823092E-02, 0.26475871E-03,
         0.10984284E-02,-0.94178133E-04,-0.29518903E-03,-0.34515437E-03,
        -0.14498901E-03,-0.12038877E-03, 0.20628994E-03, 0.99936570E-03,
        -0.41741613E-03,-0.58441347E-03, 0.17096348E-03, 0.53942105E-03,
        -0.35443399E-03,-0.33422135E-03, 0.13789161E-02, 0.29712886E-03,
        -0.52222196E-03,-0.26421202E-03, 0.39295095E-03, 0.40099127E-03,
        -0.27005718E-03, 0.43283586E-03, 0.75291598E-03,-0.13170148E-03,
        -0.59410941E-03,
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

    float v_y_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]            *x41
        +coeff[  2]                *x51
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x22
        +coeff[  5]*x11        *x41
        +coeff[  6]            *x41*x51
        +coeff[  7]*x11*x21
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[  8]    *x21        *x51
        +coeff[  9]        *x31    *x51
        +coeff[ 10]                *x52
        +coeff[ 11]    *x22    *x41
        +coeff[ 12]    *x21    *x41*x51
        +coeff[ 13]            *x44
        +coeff[ 14]        *x31*x43
        +coeff[ 15]            *x41*x52
        +coeff[ 16]*x11
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 17]*x11*x21    *x41
        +coeff[ 18]    *x22*x31    *x52
        +coeff[ 19]    *x22    *x41*x53
        +coeff[ 20]        *x31
        +coeff[ 21]            *x42
        +coeff[ 22]        *x31*x41
        +coeff[ 23]        *x32
        +coeff[ 24]    *x21    *x42
        +coeff[ 25]*x12
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 26]*x11            *x51
        +coeff[ 27]            *x42*x51
        +coeff[ 28]    *x22*x31
        +coeff[ 29]*x11*x21*x31
        +coeff[ 30]    *x22        *x51
        +coeff[ 31]*x11        *x41*x51
        +coeff[ 32]    *x22    *x42
        +coeff[ 33]*x11*x21        *x51
        +coeff[ 34]        *x31    *x52
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 35]    *x22*x31*x41
        +coeff[ 36]    *x23        *x51
        +coeff[ 37]    *x23    *x42
        +coeff[ 38]    *x21*x31    *x52
        +coeff[ 39]            *x41*x53
        +coeff[ 40]    *x23        *x52
        +coeff[ 41]    *x21
        +coeff[ 42]            *x43
        +coeff[ 43]*x11        *x42
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 44]    *x23
        +coeff[ 45]        *x31*x41*x51
        +coeff[ 46]*x11*x22
        +coeff[ 47]*x11*x21    *x42
        +coeff[ 48]    *x21    *x42*x51
        +coeff[ 49]    *x24
        +coeff[ 50]*x11*x21*x31*x41
        +coeff[ 51]    *x22    *x41*x51
        +coeff[ 52]*x11*x23
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 53]    *x21    *x41*x52
        +coeff[ 54]    *x21*x31
        +coeff[ 55]*x11    *x31
        +coeff[ 56]    *x21*x31*x41
        +coeff[ 57]    *x21    *x43
        +coeff[ 58]*x12        *x41
        +coeff[ 59]    *x21*x31*x42
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]            *x43*x51
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 62]*x11            *x52
        +coeff[ 63]                *x53
        +coeff[ 64]    *x22*x31    *x51
        +coeff[ 65]*x12*x22
        +coeff[ 66]    *x21    *x43*x51
        +coeff[ 67]    *x22    *x42*x51
        +coeff[ 68]    *x21        *x53
        +coeff[ 69]        *x31    *x53
        +coeff[ 70]    *x22    *x41*x52
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 71]*x11    *x31*x41*x52
        +coeff[ 72]            *x44*x52
        +coeff[ 73]        *x31*x42
        +coeff[ 74]        *x32*x41
        +coeff[ 75]    *x21*x31    *x51
        +coeff[ 76]        *x32    *x51
        +coeff[ 77]*x12*x21
        +coeff[ 78]*x12    *x31
        +coeff[ 79]    *x21        *x52
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 80]    *x23*x31
        +coeff[ 81]    *x22*x32
        +coeff[ 82]*x12            *x51
        +coeff[ 83]*x11*x22    *x41
        +coeff[ 84]*x11        *x42*x51
        +coeff[ 85]*x11*x21*x32
        +coeff[ 86]    *x22    *x43
        +coeff[ 87]*x12*x21    *x41
        +coeff[ 88]*x11*x21    *x41*x51
    ;
    v_y_l5p77_sex                             =v_y_l5p77_sex
        +coeff[ 89]*x11*x21*x31    *x51
        +coeff[ 90]    *x22        *x52
        +coeff[ 91]*x11    *x31*x43
        +coeff[ 92]*x12        *x41*x51
        +coeff[ 93]*x11        *x41*x52
        +coeff[ 94]    *x21*x31*x42*x51
        +coeff[ 95]    *x21*x34
        +coeff[ 96]*x11*x23    *x41
        ;

    return v_y_l5p77_sex                             ;
}
float p_l5p77_sex                             (float *x,int m){
    int ncoeff= 88;
    float avdat= -0.8726283E-02;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 89]={
         0.31087459E-02,-0.17550972E-02, 0.60103172E-02,-0.32550465E-01,
         0.11132461E-01,-0.15728759E-01, 0.21492946E-02, 0.20090608E-01,
         0.65019373E-02,-0.32881526E-02,-0.36899175E-02,-0.20557264E-01,
        -0.68714907E-02, 0.13136443E-02,-0.27759343E-02, 0.17848238E-01,
         0.55295043E-02,-0.36096491E-02,-0.18700310E-02, 0.16883461E-03,
        -0.16874974E-03,-0.44804325E-03,-0.66132176E-04,-0.14528082E-03,
         0.78978599E-03, 0.44170893E-02, 0.48174555E-02,-0.86248352E-03,
         0.43175304E-02, 0.30869315E-02,-0.49025076E-03,-0.63911907E-03,
         0.43229964E-02, 0.19765997E-02, 0.36523244E-03,-0.39645040E-03,
         0.10976845E-02, 0.16878002E-02, 0.94880757E-03, 0.29695451E-02,
         0.14757375E-02,-0.88586583E-03, 0.25625774E-02, 0.12724109E-02,
        -0.12439745E-02, 0.23487753E-03, 0.12303308E-02, 0.66318561E-03,
        -0.25460871E-02, 0.21825028E-05, 0.62858849E-03, 0.13259274E-02,
         0.68047480E-03, 0.51072915E-03,-0.11788579E-02,-0.34821380E-04,
        -0.20769864E-03,-0.29133211E-03, 0.14745751E-02,-0.61941094E-03,
        -0.13578311E-03, 0.17639925E-03,-0.14418671E-02,-0.14228341E-02,
        -0.14806696E-03, 0.49855199E-03,-0.25576004E-03,-0.16686295E-02,
        -0.21363680E-03,-0.34475076E-03, 0.68151834E-03, 0.85092586E-04,
        -0.10972444E-02,-0.11117007E-03,-0.86817985E-04, 0.56695804E-04,
         0.46252844E-04,-0.11806986E-02, 0.13195942E-03,-0.11505177E-03,
         0.35964517E-03,-0.46569301E-03,-0.28071320E-03,-0.12628649E-02,
        -0.22171719E-03, 0.55413943E-03, 0.17941769E-03,-0.43722929E-03,
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
    float x25 = x24*x2;
    float x26 = x25*x2;
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
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]        *x31
        +coeff[  3]            *x41
        +coeff[  4]                *x51
        +coeff[  5]    *x22
        +coeff[  6]    *x21*x31
        +coeff[  7]    *x21    *x41
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[  8]            *x42
        +coeff[  9]    *x21        *x51
        +coeff[ 10]        *x31    *x51
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21
        +coeff[ 13]                *x52
        +coeff[ 14]*x11        *x41
        +coeff[ 15]    *x22    *x41
        +coeff[ 16]    *x21    *x42
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 17]    *x22        *x51
        +coeff[ 18]    *x22*x31*x41
        +coeff[ 19]        *x33*x41
        +coeff[ 20]    *x25
        +coeff[ 21]    *x24*x31
        +coeff[ 22]            *x43*x52
        +coeff[ 23]        *x35*x41
        +coeff[ 24]        *x32
        +coeff[ 25]        *x31*x41
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 26]    *x22*x31
        +coeff[ 27]*x12
        +coeff[ 28]    *x24
        +coeff[ 29]*x11*x21    *x41
        +coeff[ 30]*x11*x23*x31
        +coeff[ 31]*x11
        +coeff[ 32]    *x23
        +coeff[ 33]    *x21*x31*x41
        +coeff[ 34]*x11            *x51
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 35]            *x42*x51
        +coeff[ 36]*x11*x22
        +coeff[ 37]*x11*x21*x31
        +coeff[ 38]        *x31    *x52
        +coeff[ 39]            *x41*x52
        +coeff[ 40]*x11        *x42
        +coeff[ 41]*x11*x21        *x51
        +coeff[ 42]*x11*x23
        +coeff[ 43]    *x22    *x43
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 44]            *x41*x53
        +coeff[ 45]            *x45
        +coeff[ 46]*x11*x22    *x43
        +coeff[ 47]*x11    *x31*x41
        +coeff[ 48]    *x21    *x43
        +coeff[ 49]*x11*x22*x31
        +coeff[ 50]*x12        *x41
        +coeff[ 51]*x11*x22    *x41
        +coeff[ 52]*x12*x22
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 53]*x12*x21    *x41
        +coeff[ 54]    *x23*x31*x42
        +coeff[ 55]*x13    *x31
        +coeff[ 56]*x11    *x31
        +coeff[ 57]        *x32*x41
        +coeff[ 58]            *x43
        +coeff[ 59]    *x22*x32
        +coeff[ 60]*x11    *x31    *x51
        +coeff[ 61]*x12    *x31
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 62]*x11*x21*x31*x41
        +coeff[ 63]*x11*x21    *x42
        +coeff[ 64]*x12            *x51
        +coeff[ 65]*x11*x21    *x41*x51
        +coeff[ 66]        *x31    *x53
        +coeff[ 67]    *x22    *x42*x51
        +coeff[ 68]    *x25        *x51
        +coeff[ 69]*x11*x22        *x52
        +coeff[ 70]    *x26    *x41
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 71]        *x31*x44*x51
        +coeff[ 72]    *x25*x32*x41
        +coeff[ 73]    *x21*x31    *x51
        +coeff[ 74]    *x21    *x41*x51
        +coeff[ 75]        *x32    *x51
        +coeff[ 76]*x11    *x32
        +coeff[ 77]    *x21*x31*x42
        +coeff[ 78]*x11        *x41*x51
        +coeff[ 79]*x12*x21
    ;
    v_p_l5p77_sex                             =v_p_l5p77_sex
        +coeff[ 80]    *x22        *x52
        +coeff[ 81]            *x43*x51
        +coeff[ 82]*x11*x21*x32
        +coeff[ 83]    *x23    *x42
        +coeff[ 84]*x11        *x43
        +coeff[ 85]    *x21    *x44
        +coeff[ 86]    *x23        *x52
        +coeff[ 87]    *x25*x31
        ;

    return v_p_l5p77_sex                             ;
}
float l_l5p77_sex                             (float *x,int m){
    int ncoeff= 81;
    float avdat= -0.1987057E-01;
    float xmin[10]={
        -0.14998E-01,-0.44508E-01,-0.14989E-01,-0.33113E-01,-0.47801E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.48696E-01, 0.14990E-01, 0.25296E-01, 0.49919E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 82]={
         0.83441073E-02,-0.24360231E+00,-0.34637661E-02,-0.33813886E-01,
         0.41575912E-01,-0.26284715E-01,-0.99122055E-01,-0.26512446E-01,
        -0.21665953E-01, 0.11416291E-01,-0.35530519E-01, 0.39644921E-02,
        -0.79143019E-02,-0.27368957E-01,-0.20160317E-01, 0.39728206E-01,
        -0.17377678E-01, 0.40605819E-03,-0.30036096E-03,-0.37966294E-02,
        -0.79438174E-02,-0.34970627E-02, 0.25771398E-02,-0.16564309E-02,
         0.25794279E-01,-0.23089738E-02, 0.49719093E-02,-0.44403509E-02,
         0.38621280E-01, 0.23833012E-01,-0.92207291E-03,-0.56213476E-02,
         0.42304178E-02, 0.59777969E-02, 0.77205044E-02,-0.77581900E-03,
        -0.24938094E-03,-0.82184188E-03, 0.92981802E-03,-0.27437564E-02,
         0.47346670E-02,-0.41362992E-02, 0.45841848E-02,-0.21131453E-03,
         0.13493160E-02, 0.14714708E-02, 0.71005872E-03, 0.93072787E-03,
         0.67209657E-02,-0.10662257E-01,-0.94233248E-02,-0.74733078E-03,
        -0.27614215E-02, 0.15268204E-02, 0.71182516E-02, 0.50195251E-02,
        -0.13378052E-03, 0.73608979E-04, 0.12467285E-02,-0.45532800E-03,
         0.42610307E-03, 0.24744534E-03, 0.33570733E-03, 0.10494553E-03,
        -0.32920748E-03, 0.30218664E-03, 0.31781604E-02,-0.33803007E-02,
        -0.48886548E-03, 0.10078787E-02,-0.10469309E-02,-0.80185715E-03,
         0.51515625E-03,-0.81025198E-03,-0.52791811E-03,-0.73398004E-03,
         0.18900851E-02, 0.66593068E-03, 0.25862025E-02, 0.96391345E-03,
         0.41168286E-02,
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
    float x54 = x53*x5;

//                 function

    float v_l_l5p77_sex                             =avdat
        +coeff[  0]
        +coeff[  1]    *x21
        +coeff[  2]            *x41
        +coeff[  3]                *x51
        +coeff[  4]*x11
        +coeff[  5]    *x22
        +coeff[  6]    *x21    *x41
        +coeff[  7]                *x52
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[  8]*x11        *x41
        +coeff[  9]    *x23*x31
        +coeff[ 10]    *x21*x31
        +coeff[ 11]*x11*x21
        +coeff[ 12]*x11    *x31
        +coeff[ 13]    *x23
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x11*x22
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 17]            *x44
        +coeff[ 18]    *x22*x33
        +coeff[ 19]    *x23*x31*x41
        +coeff[ 20]            *x42
        +coeff[ 21]    *x21        *x51
        +coeff[ 22]            *x41*x51
        +coeff[ 23]*x11            *x51
        +coeff[ 24]    *x21*x31*x41
        +coeff[ 25]    *x21        *x52
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 26]                *x53
        +coeff[ 27]*x12*x21
        +coeff[ 28]    *x23    *x41
        +coeff[ 29]*x11*x22    *x41
        +coeff[ 30]        *x31
        +coeff[ 31]    *x22*x31
        +coeff[ 32]    *x21*x32
        +coeff[ 33]*x11        *x42
        +coeff[ 34]*x11*x22*x31
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 35]*x11*x22*x31*x41
        +coeff[ 36]*x11    *x33*x41
        +coeff[ 37]*x12
        +coeff[ 38]        *x31*x41*x51
        +coeff[ 39]*x11*x21    *x41
        +coeff[ 40]*x11    *x31*x41
        +coeff[ 41]    *x24
        +coeff[ 42]*x12*x21    *x41
        +coeff[ 43]        *x31*x41
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 44]    *x21*x31    *x51
        +coeff[ 45]    *x21    *x41*x51
        +coeff[ 46]*x11    *x32
        +coeff[ 47]*x11        *x41*x51
        +coeff[ 48]    *x22    *x42
        +coeff[ 49]    *x21*x31*x42
        +coeff[ 50]    *x21    *x43
        +coeff[ 51]    *x21*x31*x41*x51
        +coeff[ 52]*x11*x23
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 53]*x12*x21*x31
        +coeff[ 54]    *x24    *x41
        +coeff[ 55]    *x21    *x44
        +coeff[ 56]    *x24*x31*x41
        +coeff[ 57]    *x23*x32*x41
        +coeff[ 58]            *x43
        +coeff[ 59]    *x22        *x51
        +coeff[ 60]            *x41*x52
        +coeff[ 61]*x11    *x31    *x51
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 62]*x11            *x52
        +coeff[ 63]*x12    *x31
        +coeff[ 64]*x13
        +coeff[ 65]    *x22*x32
        +coeff[ 66]    *x22*x31*x41
        +coeff[ 67]    *x21*x32*x41
        +coeff[ 68]    *x23        *x51
        +coeff[ 69]    *x21    *x42*x51
        +coeff[ 70]    *x21    *x41*x52
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 71]            *x42*x52
        +coeff[ 72]    *x21        *x53
        +coeff[ 73]                *x54
        +coeff[ 74]*x11*x21    *x41*x51
        +coeff[ 75]*x11        *x42*x51
        +coeff[ 76]    *x24*x31
        +coeff[ 77]        *x33*x42
        +coeff[ 78]*x11*x23    *x41
        +coeff[ 79]*x13*x21    *x42
    ;
    v_l_l5p77_sex                             =v_l_l5p77_sex
        +coeff[ 80]    *x23*x31*x43
        ;

    return v_l_l5p77_sex                             ;
}
}
