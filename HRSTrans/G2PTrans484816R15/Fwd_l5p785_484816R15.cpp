#include "Fwd_l5p785_484816R15.h"

namespace S484816R15
{
float x_l5p785_fp                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.2267227E-01;
    float xmin[10]={
        -0.14988E-01,-0.47247E-01,-0.14993E-01,-0.32667E-01,-0.46394E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.46621E-01, 0.14972E-01, 0.25411E-01, 0.49909E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.34130765E-02,-0.37364822E-01, 0.64359063E+00,-0.51526140E-01,
         0.43762371E-01, 0.11521450E-01, 0.35487361E-01, 0.13356934E-01,
        -0.13595128E-01,-0.29254250E-01, 0.48478203E-04,-0.28343536E-02,
         0.74754567E-02,-0.33870300E-02, 0.55735903E-02, 0.13365731E-01,
         0.89096548E-02,-0.10647596E-01,-0.91178226E-03,-0.27111326E-02,
         0.26007928E-02, 0.26802076E-02,-0.32504895E-02, 0.50107338E-02,
         0.48910370E-02,-0.20084010E-01, 0.13849534E-02, 0.28694002E-02,
        -0.17019579E-02,-0.10775101E-02,-0.99359183E-02,-0.26208111E-02,
        -0.41235262E-02,-0.48270723E-03, 0.70351874E-02, 0.41066427E-02,
         0.61823253E-03, 0.11000105E-02,-0.11046140E-02, 0.65015466E-03,
         0.11425747E-02, 0.34384371E-02,-0.20811602E-02,-0.53026583E-02,
         0.32771074E-02, 0.14176104E-02, 0.80179080E-03, 0.20275726E-02,
         0.89105703E-02,-0.54695243E-02,
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

    float v_x_l5p785_fp                             =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]                *x51
        +coeff[  3]                *x52
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21    *x41    
        +coeff[  7]    *x21*x31        
    ;
    v_x_l5p785_fp                             =v_x_l5p785_fp                             
        +coeff[  8]            *x42    
        +coeff[  9]    *x21    *x42    
        +coeff[ 10]    *x21            
        +coeff[ 11]            *x41    
        +coeff[ 12]*x11        *x41    
        +coeff[ 13]        *x31*x41    
        +coeff[ 14]                *x53
        +coeff[ 15]    *x21    *x41*x51
        +coeff[ 16]    *x23            
    ;
    v_x_l5p785_fp                             =v_x_l5p785_fp                             
        +coeff[ 17]    *x21*x31*x41    
        +coeff[ 18]        *x31        
        +coeff[ 19]*x11            *x51
        +coeff[ 20]*x11    *x31        
        +coeff[ 21]            *x41*x51
        +coeff[ 22]    *x21        *x52
        +coeff[ 23]*x11*x22            
        +coeff[ 24]    *x22    *x41    
        +coeff[ 25]    *x23    *x41    
    ;
    v_x_l5p785_fp                             =v_x_l5p785_fp                             
        +coeff[ 26]*x12*x21            
        +coeff[ 27]    *x21*x31    *x51
        +coeff[ 28]    *x21*x32        
        +coeff[ 29]*x12*x22            
        +coeff[ 30]*x11*x22    *x41    
        +coeff[ 31]    *x21    *x42*x51
        +coeff[ 32]    *x23*x31        
        +coeff[ 33]    *x22*x31*x41    
        +coeff[ 34]    *x22    *x42    
    ;
    v_x_l5p785_fp                             =v_x_l5p785_fp                             
        +coeff[ 35]    *x22    *x42*x51
        +coeff[ 36]        *x31    *x51
        +coeff[ 37]*x11*x21    *x41    
        +coeff[ 38]*x11    *x31*x41    
        +coeff[ 39]    *x22*x31        
        +coeff[ 40]            *x43    
        +coeff[ 41]*x11*x22        *x51
        +coeff[ 42]*x11*x22*x31        
        +coeff[ 43]*x11*x21    *x42    
    ;
    v_x_l5p785_fp                             =v_x_l5p785_fp                             
        +coeff[ 44]    *x23        *x51
        +coeff[ 45]*x13        *x41*x51
        +coeff[ 46]*x13*x23            
        +coeff[ 47]*x11*x23        *x52
        +coeff[ 48]*x11*x23    *x42    
        +coeff[ 49]*x12*x22*x32*x41    
        ;

    return v_x_l5p785_fp                             ;
}
float t_l5p785_fp                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.3411545E-02;
    float xmin[10]={
        -0.14988E-01,-0.47247E-01,-0.14993E-01,-0.32667E-01,-0.46394E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.46621E-01, 0.14972E-01, 0.25411E-01, 0.49909E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.53207792E-03,-0.33823878E-02,-0.18817483E-01, 0.89439374E-04,
         0.11137834E+00, 0.52485834E-02,-0.10345989E-01, 0.12956716E-02,
        -0.24398796E-02, 0.16271793E-03,-0.96347288E-03,-0.22775205E-02,
         0.22427717E-02,-0.81448851E-03,-0.79633592E-03,-0.41826366E-03,
         0.36883421E-03, 0.10929212E-02,-0.19773291E-03,-0.25646854E-03,
        -0.34457262E-03, 0.34254187E-03, 0.34819733E-03, 0.33169266E-03,
         0.90879790E-03, 0.15985129E-02,-0.47422483E-04,-0.85862848E-03,
         0.21232643E-04,-0.76181983E-04, 0.16089047E-03,-0.19768106E-03,
        -0.32895879E-03,-0.33200317E-03, 0.34670750E-03, 0.14138711E-03,
         0.53107005E-03, 0.53275953E-03,-0.80059626E-03, 0.24585001E-03,
        -0.88362541E-03,-0.65614768E-04,-0.91787850E-04, 0.22223241E-03,
        -0.84502215E-04, 0.10387886E-03,-0.11682043E-03,-0.42757165E-03,
         0.17359023E-03, 0.28083549E-03,
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
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_l5p785_fp                             =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x21        *x51
        +coeff[  6]                *x52
        +coeff[  7]    *x22            
    ;
    v_t_l5p785_fp                             =v_t_l5p785_fp                             
        +coeff[  8]            *x42    
        +coeff[  9]*x11*x21            
        +coeff[ 10]    *x21    *x41    
        +coeff[ 11]    *x21    *x42    
        +coeff[ 12]    *x21    *x41*x51
        +coeff[ 13]    *x21        *x52
        +coeff[ 14]        *x31*x41    
        +coeff[ 15]*x11            *x51
        +coeff[ 16]            *x41*x51
    ;
    v_t_l5p785_fp                             =v_t_l5p785_fp                             
        +coeff[ 17]                *x53
        +coeff[ 18]    *x21*x31        
        +coeff[ 19]*x11        *x41    
        +coeff[ 20]*x11*x22            
        +coeff[ 21]*x11        *x42    
        +coeff[ 22]    *x21*x31    *x51
        +coeff[ 23]            *x42*x51
        +coeff[ 24]*x11*x23            
        +coeff[ 25]    *x22    *x42    
    ;
    v_t_l5p785_fp                             =v_t_l5p785_fp                             
        +coeff[ 26]        *x33    *x51
        +coeff[ 27]    *x23    *x43    
        +coeff[ 28]        *x31        
        +coeff[ 29]*x11    *x31        
        +coeff[ 30]        *x31    *x51
        +coeff[ 31]    *x23            
        +coeff[ 32]    *x22*x31        
        +coeff[ 33]    *x22    *x41    
        +coeff[ 34]            *x43    
    ;
    v_t_l5p785_fp                             =v_t_l5p785_fp                             
        +coeff[ 35]*x11            *x52
        +coeff[ 36]*x11*x22        *x51
        +coeff[ 37]    *x23        *x51
        +coeff[ 38]    *x23    *x41*x51
        +coeff[ 39]        *x31*x41*x53
        +coeff[ 40]    *x23    *x42*x51
        +coeff[ 41]*x12                
        +coeff[ 42]*x12*x21            
        +coeff[ 43]        *x31*x42    
    ;
    v_t_l5p785_fp                             =v_t_l5p785_fp                             
        +coeff[ 44]    *x22        *x51
        +coeff[ 45]*x11        *x41*x51
        +coeff[ 46]            *x41*x52
        +coeff[ 47]    *x23    *x41    
        +coeff[ 48]*x11*x21*x31*x41    
        +coeff[ 49]    *x22*x31*x41    
        ;

    return v_t_l5p785_fp                             ;
}
float y_l5p785_fp                             (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.8176396E-02;
    float xmin[10]={
        -0.14988E-01,-0.47247E-01,-0.14993E-01,-0.32667E-01,-0.46394E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.46621E-01, 0.14972E-01, 0.25411E-01, 0.49909E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 98]={
         0.57107909E-02,-0.58257539E-01, 0.30012408E-02, 0.49144565E-03,
         0.11917459E-01, 0.55577625E-02, 0.80952933E-02,-0.16498631E-01,
        -0.21205309E-02,-0.37700420E-02,-0.66721155E-02,-0.25925501E-02,
        -0.44445661E-02,-0.27213611E-02, 0.10485073E-01,-0.55081556E-02,
         0.64244606E-02,-0.27003074E-02, 0.42423932E-02, 0.44810432E-02,
         0.44514365E-02, 0.26847434E-02, 0.14028304E-03, 0.98650251E-03,
         0.80175087E-03,-0.77749783E-03,-0.17637697E-02, 0.14608323E-02,
         0.28656246E-02, 0.11956474E-02, 0.10327016E-02, 0.14611003E-02,
         0.88499990E-04,-0.20868274E-03, 0.40101112E-03, 0.84283249E-03,
         0.10038463E-02, 0.20899465E-02,-0.11548385E-02,-0.39862711E-02,
        -0.29907087E-02, 0.75218518E-03,-0.14463580E-02,-0.42990330E-02,
        -0.96345681E-03, 0.25210995E-03,-0.12659912E-02,-0.97297409E-04,
         0.26841543E-02,-0.26423816E-03, 0.40903370E-03,-0.27899629E-04,
        -0.17925751E-03, 0.45921674E-03,-0.26185455E-03, 0.21602948E-03,
        -0.13180878E-02, 0.67381695E-03,-0.19350646E-02, 0.25269394E-02,
        -0.36888584E-03, 0.66152139E-03,-0.16321390E-02,-0.52520158E-02,
         0.22457668E-02,-0.10656066E-02, 0.61929127E-03,-0.52373734E-03,
        -0.76621416E-03, 0.24578103E-03,-0.19356264E-02, 0.78690326E-03,
        -0.61249601E-04,-0.12579975E-02,-0.86272508E-03,-0.16462979E-03,
        -0.29996521E-03,-0.49910275E-03,-0.57246396E-03, 0.36837859E-03,
        -0.79935370E-03, 0.85172546E-03,-0.34987673E-03, 0.43201610E-03,
        -0.36269205E-03, 0.86242624E-03,-0.60925062E-03,-0.69959363E-03,
        -0.29991911E-03, 0.19412909E-02,-0.58366085E-03, 0.16546226E-03,
        -0.24080810E-02,-0.13701945E-02, 0.20883977E-02, 0.18659842E-02,
         0.82562445E-03,
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
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_l5p785_fp                             =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]*x11                
        +coeff[  4]                *x51
        +coeff[  5]            *x42    
        +coeff[  6]    *x21    *x41    
        +coeff[  7]    *x22            
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[  8]*x11        *x41    
        +coeff[  9]            *x41*x51
        +coeff[ 10]*x11*x21            
        +coeff[ 11]    *x21        *x51
        +coeff[ 12]        *x31    *x51
        +coeff[ 13]                *x52
        +coeff[ 14]    *x22    *x41    
        +coeff[ 15]    *x21    *x41*x51
        +coeff[ 16]            *x41*x52
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 17]        *x31        
        +coeff[ 18]        *x31*x41    
        +coeff[ 19]    *x22*x31        
        +coeff[ 20]*x11*x21    *x41    
        +coeff[ 21]        *x31    *x52
        +coeff[ 22]    *x23    *x42    
        +coeff[ 23]    *x22    *x41*x53
        +coeff[ 24]        *x32        
        +coeff[ 25]*x12                
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 26]            *x42*x51
        +coeff[ 27]*x11*x21*x31        
        +coeff[ 28]    *x22        *x51
        +coeff[ 29]*x11        *x41*x51
        +coeff[ 30]*x11*x21        *x51
        +coeff[ 31]    *x21        *x52
        +coeff[ 32]            *x45    
        +coeff[ 33]    *x21*x31        
        +coeff[ 34]*x11            *x51
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 35]    *x21*x31*x41    
        +coeff[ 36]*x11        *x42    
        +coeff[ 37]    *x23            
        +coeff[ 38]        *x31*x41*x51
        +coeff[ 39]    *x22    *x42    
        +coeff[ 40]    *x22*x31*x41    
        +coeff[ 41]    *x21    *x44    
        +coeff[ 42]    *x23        *x51
        +coeff[ 43]            *x41*x53
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 44]    *x21    *x44*x51
        +coeff[ 45]*x11    *x31*x41*x52
        +coeff[ 46]    *x23    *x44    
        +coeff[ 47]*x11    *x31        
        +coeff[ 48]    *x21    *x42    
        +coeff[ 49]        *x32*x41    
        +coeff[ 50]*x11    *x31*x41    
        +coeff[ 51]            *x44    
        +coeff[ 52]        *x32    *x51
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 53]*x12        *x41    
        +coeff[ 54]*x12*x21            
        +coeff[ 55]*x12    *x31        
        +coeff[ 56]            *x43*x51
        +coeff[ 57]    *x23    *x41    
        +coeff[ 58]*x11*x21    *x42    
        +coeff[ 59]    *x24            
        +coeff[ 60]*x11            *x52
        +coeff[ 61]*x11*x22    *x41    
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 62]*x11*x21*x31*x41    
        +coeff[ 63]    *x22    *x41*x51
        +coeff[ 64]*x11*x23            
        +coeff[ 65]    *x22*x31    *x51
        +coeff[ 66]*x12*x22            
        +coeff[ 67]    *x21        *x53
        +coeff[ 68]        *x31    *x53
        +coeff[ 69]            *x42*x53
        +coeff[ 70]    *x21    *x41*x53
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 71]            *x43    
        +coeff[ 72]    *x21*x32        
        +coeff[ 73]    *x21    *x43    
        +coeff[ 74]    *x21*x31*x42    
        +coeff[ 75]*x11    *x31    *x51
        +coeff[ 76]    *x21    *x42*x51
        +coeff[ 77]    *x23*x31        
        +coeff[ 78]    *x22*x32        
        +coeff[ 79]                *x53
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 80]*x11        *x42*x51
        +coeff[ 81]            *x42*x52
        +coeff[ 82]*x11*x21*x32        
        +coeff[ 83]*x12*x21    *x41    
        +coeff[ 84]*x11    *x31*x41*x51
        +coeff[ 85]    *x21    *x41*x52
        +coeff[ 86]*x11*x22    *x42    
        +coeff[ 87]*x11*x21    *x42*x51
        +coeff[ 88]    *x23*x31    *x51
    ;
    v_y_l5p785_fp                             =v_y_l5p785_fp                             
        +coeff[ 89]    *x22    *x41*x52
        +coeff[ 90]    *x22    *x44*x51
        +coeff[ 91]*x12            *x53
        +coeff[ 92]    *x24    *x42*x51
        +coeff[ 93]*x11*x22    *x43*x51
        +coeff[ 94]    *x23    *x41*x53
        +coeff[ 95]    *x24    *x45    
        +coeff[ 96]*x11*x22*x32    *x52
        ;

    return v_y_l5p785_fp                             ;
}
float p_l5p785_fp                             (float *x,int m){
    int ncoeff= 90;
    float avdat= -0.6752405E-02;
    float xmin[10]={
        -0.14988E-01,-0.47247E-01,-0.14993E-01,-0.32667E-01,-0.46394E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E-01, 0.46621E-01, 0.14972E-01, 0.25411E-01, 0.49909E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 91]={
         0.31118165E-02, 0.29739486E-02, 0.60560764E-02,-0.31288404E-01,
         0.10556458E-01,-0.15422394E-01, 0.14951651E-02, 0.18569820E-01,
         0.37267830E-02, 0.49918015E-02,-0.25324312E-02,-0.37073093E-02,
        -0.20262284E-01,-0.63094245E-02, 0.14590473E-02,-0.34529406E-02,
         0.16884951E-01, 0.51195920E-02,-0.34447764E-02, 0.94627387E-04,
        -0.85537955E-04, 0.43848911E-03, 0.71786856E-03, 0.28909331E-02,
         0.44393968E-02, 0.19855807E-02, 0.56142011E-03, 0.22967786E-02,
        -0.71447267E-03, 0.22196702E-02, 0.28257503E-02,-0.45455783E-03,
        -0.37532230E-03, 0.39082337E-02, 0.12881315E-02, 0.83943136E-03,
        -0.16043334E-02, 0.18440231E-02,-0.85456611E-03, 0.22449722E-02,
         0.43416523E-04,-0.11522343E-03, 0.68487995E-03,-0.10077141E-02,
         0.82713593E-03,-0.25368433E-02, 0.55771822E-03, 0.13407511E-02,
        -0.12687207E-02,-0.12243782E-02, 0.70547924E-03,-0.10542022E-02,
        -0.11371035E-03, 0.79960172E-03, 0.50784665E-03,-0.25554406E-03,
        -0.36790193E-03,-0.64175040E-03,-0.10185582E-03,-0.19857714E-03,
        -0.18369303E-03,-0.11914542E-03, 0.18115919E-03,-0.29845303E-03,
        -0.12140220E-03, 0.29187277E-03, 0.50716772E-03, 0.43913510E-03,
        -0.11428407E-02, 0.11664178E-02, 0.79545673E-04,-0.54117129E-03,
        -0.26331126E-03,-0.14448229E-02, 0.32302813E-03, 0.33143978E-03,
         0.49591664E-03,-0.61188279E-04,-0.50757214E-03, 0.70025888E-03,
         0.16477953E-03, 0.16129929E-03, 0.17419006E-03,-0.15987656E-02,
        -0.13824989E-03,-0.18401413E-03, 0.23767413E-03,-0.14778969E-03,
         0.46284746E-04,-0.90235693E-03,
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

    float v_p_l5p785_fp                             =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21*x31        
        +coeff[  7]    *x21    *x41    
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
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
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 17]    *x21    *x42    
        +coeff[ 18]    *x22        *x51
        +coeff[ 19]*x11    *x32        
        +coeff[ 20]    *x24*x31        
        +coeff[ 21]*x11                
        +coeff[ 22]        *x32        
        +coeff[ 23]    *x23            
        +coeff[ 24]    *x22*x31        
        +coeff[ 25]    *x21*x31*x41    
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 26]*x11            *x51
        +coeff[ 27]            *x43    
        +coeff[ 28]*x12                
        +coeff[ 29]*x11*x21    *x41    
        +coeff[ 30]            *x41*x52
        +coeff[ 31]*x11    *x31        
        +coeff[ 32]            *x42*x51
        +coeff[ 33]    *x24            
        +coeff[ 34]*x11*x21*x31        
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 35]        *x31    *x52
        +coeff[ 36]    *x22*x31*x41    
        +coeff[ 37]*x11        *x42    
        +coeff[ 38]*x11*x21        *x51
        +coeff[ 39]*x11*x23            
        +coeff[ 40]*x11    *x33*x41    
        +coeff[ 41]        *x32*x41    
        +coeff[ 42]*x11*x22            
        +coeff[ 43]    *x23    *x41    
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 44]*x11    *x31*x41    
        +coeff[ 45]    *x21    *x43    
        +coeff[ 46]*x12        *x41    
        +coeff[ 47]*x11*x22    *x41    
        +coeff[ 48]*x11*x21*x31*x41    
        +coeff[ 49]*x11*x21    *x42    
        +coeff[ 50]*x11*x21    *x41*x51
        +coeff[ 51]            *x41*x53
        +coeff[ 52]    *x23*x31*x42    
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 53]    *x24    *x41*x51
        +coeff[ 54]        *x31*x42    
        +coeff[ 55]    *x21*x31    *x51
        +coeff[ 56]    *x21    *x41*x51
        +coeff[ 57]    *x22*x32        
        +coeff[ 58]                *x53
        +coeff[ 59]*x11    *x31    *x51
        +coeff[ 60]*x12*x21            
        +coeff[ 61]*x11            *x52
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 62]*x12    *x31        
        +coeff[ 63]*x11*x21*x32        
        +coeff[ 64]*x12            *x51
        +coeff[ 65]    *x21    *x44    
        +coeff[ 66]*x12*x22            
        +coeff[ 67]*x12*x21    *x41    
        +coeff[ 68]*x11*x22    *x42    
        +coeff[ 69]    *x22    *x42*x52
        +coeff[ 70]    *x21        *x52
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 71]    *x23*x31        
        +coeff[ 72]    *x23        *x51
        +coeff[ 73]    *x21*x31*x42    
        +coeff[ 74]    *x22*x31    *x51
        +coeff[ 75]    *x21*x31*x41*x51
        +coeff[ 76]    *x21    *x42*x51
        +coeff[ 77]        *x31*x42*x51
        +coeff[ 78]            *x43*x51
        +coeff[ 79]    *x24    *x41    
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 80]    *x21*x31    *x52
        +coeff[ 81]    *x21    *x41*x52
        +coeff[ 82]    *x23*x32        
        +coeff[ 83]    *x23    *x42    
        +coeff[ 84]*x11*x22        *x51
        +coeff[ 85]*x11        *x43    
        +coeff[ 86]*x11*x21*x31    *x51
        +coeff[ 87]        *x31    *x53
        +coeff[ 88]    *x22*x31*x41*x51
    ;
    v_p_l5p785_fp                             =v_p_l5p785_fp                             
        +coeff[ 89]    *x22    *x42*x51
        ;

    return v_p_l5p785_fp                             ;
}
// float l_r6_h90s_fp                            (float *x,int m){
//     int ncoeff= 83;
//     float avdat= -0.7604773E-02;
//     float xmin[10]={
//         -0.14988E-01,-0.47247E-01,-0.14993E-01,-0.32667E-01,-0.46394E-01,
//          0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
//     float xmax[10]={
//          0.15000E-01, 0.46621E-01, 0.14972E-01, 0.25411E-01, 0.49909E-01,
//          0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
//     float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//     float coeff[ 84]={
//          0.68144388E-02,-0.25147295E+00, 0.13013301E-01,-0.32102790E-01,
//          0.41630790E-01,-0.15869312E-01,-0.91396213E-01,-0.26322620E-01,
//          0.93306694E-02,-0.20866483E-01, 0.10946541E-01,-0.41654904E-03,
//          0.49361167E-02,-0.33570770E-01,-0.11965007E-01,-0.78005111E-02,
//         -0.23962090E-01,-0.28495587E-01, 0.37867703E-01,-0.16569924E-01,
//         -0.34104118E-02,-0.80615273E-02, 0.25176777E-01,-0.23253227E-02,
//          0.47656521E-02,-0.59454031E-02,-0.42534028E-02, 0.36639344E-01,
//          0.22658460E-01,-0.34726425E-02,-0.33983865E-02, 0.18648598E-02,
//         -0.18020272E-02, 0.44247140E-02, 0.58843070E-02, 0.76899603E-02,
//         -0.94240306E-04,-0.11735513E-02, 0.45820489E-02, 0.45919982E-02,
//         -0.55336143E-03, 0.13948147E-02, 0.16036105E-02, 0.98936877E-03,
//          0.31686260E-03, 0.87351276E-03, 0.99146215E-03,-0.49428311E-02,
//          0.67440164E-02,-0.11250781E-01,-0.94760219E-02,-0.93664450E-03,
//         -0.37883404E-02,-0.71263983E-03, 0.16104080E-02, 0.81911758E-02,
//          0.48357090E-02,-0.86280610E-03,-0.25664672E-03, 0.58916041E-04,
//         -0.19834729E-03,-0.18694664E-03, 0.19223026E-02,-0.84178703E-03,
//         -0.29576255E-03, 0.34883665E-03,-0.26544795E-03, 0.39870809E-02,
//         -0.40147430E-02,-0.11385175E-02, 0.76333148E-03, 0.63950277E-03,
//         -0.78098057E-03, 0.15753348E-02, 0.42739735E-03, 0.60231483E-03,
//          0.63895789E-03, 0.11816076E-02, 0.12050889E-02, 0.28367771E-02,
//         -0.82084088E-03,-0.23198575E-02, 0.48742821E-02,
//               0.      };
//     int ientry=0;

//     if (ientry==0){
//         ientry=1;
//         for(int i=0;i<m;i++){
//             if(xmin[i]==xmax[i]) continue;
//             scale[i]=2./(xmax[i]-xmin[i]);
//         }
//     }
// //  normalize variables between -1 and +1
//     float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
//     float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
//     float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
//     float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
//     float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
// //  set up monomials   functions
//     float x11 = x1;
//     float x12 = x11*x1;
//     float x13 = x12*x1;
//     float x21 = x2;
//     float x22 = x21*x2;
//     float x23 = x22*x2;
//     float x24 = x23*x2;
//     float x31 = x3;
//     float x32 = x31*x3;
//     float x33 = x32*x3;
//     float x34 = x33*x3;
//     float x41 = x4;
//     float x42 = x41*x4;
//     float x43 = x42*x4;
//     float x44 = x43*x4;
//     float x51 = x5;
//     float x52 = x51*x5;
//     float x53 = x52*x5;
//     float x54 = x53*x5;

// //                 function

//     float v_l_r6_h90s_fp                            =avdat
//         +coeff[  0]                    
//         +coeff[  1]    *x21            
//         +coeff[  2]            *x41    
//         +coeff[  3]                *x51
//         +coeff[  4]*x11                
//         +coeff[  5]    *x22            
//         +coeff[  6]    *x21    *x41    
//         +coeff[  7]                *x52
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[  8]*x11*x21            
//         +coeff[  9]*x11        *x41    
//         +coeff[ 10]    *x23*x31        
//         +coeff[ 11]    *x21*x33        
//         +coeff[ 12]        *x31        
//         +coeff[ 13]    *x21*x31        
//         +coeff[ 14]            *x42    
//         +coeff[ 15]*x11    *x31        
//         +coeff[ 16]    *x23            
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 17]    *x22    *x41    
//         +coeff[ 18]    *x21    *x42    
//         +coeff[ 19]*x11*x22            
//         +coeff[ 20]    *x23*x31*x41    
//         +coeff[ 21]    *x22*x31        
//         +coeff[ 22]    *x21*x31*x41    
//         +coeff[ 23]    *x21        *x52
//         +coeff[ 24]                *x53
//         +coeff[ 25]*x11*x21    *x41    
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 26]*x12*x21            
//         +coeff[ 27]    *x23    *x41    
//         +coeff[ 28]*x11*x22    *x41    
//         +coeff[ 29]        *x31*x41    
//         +coeff[ 30]    *x21        *x51
//         +coeff[ 31]            *x41*x51
//         +coeff[ 32]*x11            *x51
//         +coeff[ 33]    *x21*x32        
//         +coeff[ 34]*x11        *x42    
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 35]*x11*x22*x31        
//         +coeff[ 36]*x11    *x33*x41    
//         +coeff[ 37]*x11*x21*x31        
//         +coeff[ 38]*x11    *x31*x41    
//         +coeff[ 39]*x12*x21    *x41    
//         +coeff[ 40]        *x32        
//         +coeff[ 41]    *x21*x31    *x51
//         +coeff[ 42]    *x21    *x41*x51
//         +coeff[ 43]        *x31*x41*x51
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 44]            *x42*x51
//         +coeff[ 45]*x11    *x32        
//         +coeff[ 46]*x11        *x41*x51
//         +coeff[ 47]    *x24            
//         +coeff[ 48]    *x22    *x42    
//         +coeff[ 49]    *x21*x31*x42    
//         +coeff[ 50]    *x21    *x43    
//         +coeff[ 51]    *x21*x31*x41*x51
//         +coeff[ 52]*x11*x23            
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 53]*x12*x22            
//         +coeff[ 54]*x12*x21*x31        
//         +coeff[ 55]    *x24    *x41    
//         +coeff[ 56]    *x21    *x44    
//         +coeff[ 57]    *x24*x31*x41    
//         +coeff[ 58]    *x23*x32*x41    
//         +coeff[ 59]    *x21*x34*x41    
//         +coeff[ 60]        *x31    *x51
//         +coeff[ 61]*x12                
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 62]            *x43    
//         +coeff[ 63]    *x22        *x51
//         +coeff[ 64]*x11*x21        *x51
//         +coeff[ 65]*x11            *x52
//         +coeff[ 66]*x13                
//         +coeff[ 67]    *x22*x31*x41    
//         +coeff[ 68]    *x21*x32*x41    
//         +coeff[ 69]    *x21    *x41*x52
//         +coeff[ 70]    *x21        *x53
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 71]*x11*x21    *x42    
//         +coeff[ 72]*x11        *x42*x51
//         +coeff[ 73]    *x24*x31        
//         +coeff[ 74]        *x34*x41    
//         +coeff[ 75]    *x22*x31*x42    
//         +coeff[ 76]        *x33*x42    
//         +coeff[ 77]        *x31*x44    
//         +coeff[ 78]*x11*x24            
//         +coeff[ 79]*x11*x23    *x41    
//     ;
//     v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
//         +coeff[ 80]            *x42*x54
//         +coeff[ 81]*x11*x22*x31*x42    
//         +coeff[ 82]    *x23*x31*x43    
//         ;

//     return v_l_r6_h90s_fp                            ;
// }
}
