#include "Bwd_l5p77_484816.h"

namespace S484816
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
    int ncoeff= 28;
    float avdat=  0.1806338E-02;
    float xmin[10]={
        -0.70676E+00,-0.24649E-01,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.20444E-01, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 29]={
        -0.55276840E-02, 0.47911994E-01, 0.43851542E-02, 0.25593217E-02,
         0.38472114E-02,-0.50137630E-02, 0.44005872E-02,-0.12836595E-02,
        -0.29512662E-02, 0.10639990E-03,-0.82995590E-04, 0.13929608E-02,
        -0.44814622E-03,-0.50798483E-03, 0.28227379E-02,-0.38339228E-02,
         0.60364150E-03, 0.40065660E-03, 0.53261570E-03, 0.34851281E-03,
         0.19447175E-02, 0.47912743E-03, 0.51077513E-03,-0.23204344E-02,
        -0.10371626E-02, 0.18809825E-02,-0.90336660E-03, 0.50889765E-03,
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
    float x51 = x5;

//                 function

    float v_delta_l5p77                             =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]*x12                
        +coeff[  3]                *x51
        +coeff[  4]*x11*x21            
        +coeff[  5]    *x21*x31        
        +coeff[  6]    *x21    *x41    
        +coeff[  7]    *x22            
    ;
    v_delta_l5p77                             =v_delta_l5p77                             
        +coeff[  8]    *x21*x32        
        +coeff[  9]        *x34        
        +coeff[ 10]*x13        *x41    
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]        *x31        
        +coeff[ 13]*x11    *x31        
        +coeff[ 14]        *x32        
        +coeff[ 15]        *x31*x41    
        +coeff[ 16]    *x24            
    ;
    v_delta_l5p77                             =v_delta_l5p77                             
        +coeff[ 17]    *x21        *x51
        +coeff[ 18]    *x21            
        +coeff[ 19]            *x41    
        +coeff[ 20]            *x42    
        +coeff[ 21]*x13                
        +coeff[ 22]*x11*x22            
        +coeff[ 23]    *x22*x31        
        +coeff[ 24]*x11*x21    *x41    
        +coeff[ 25]    *x22    *x41    
    ;
    v_delta_l5p77                             =v_delta_l5p77                             
        +coeff[ 26]            *x41*x51
        +coeff[ 27]        *x32    *x51
        ;

    return v_delta_l5p77                             ;
}
float theta_l5p77                             (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.1079361E-02;
    float xmin[10]={
        -0.70676E+00,-0.24649E-01,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.20444E-01, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.73364414E-02,-0.71991766E-02,-0.52481033E-01,-0.69722626E-03,
        -0.24773958E-02, 0.61953538E-02,-0.24543686E-02,-0.33725344E-01,
         0.28755948E-01, 0.26158299E-02, 0.46845856E-02,-0.96575241E-03,
        -0.97825285E-02, 0.68466491E-02, 0.20929264E-01,-0.19430449E-02,
         0.10018193E-01, 0.94507600E-03,-0.40358347E-02,-0.35768040E-01,
         0.60541850E-01,-0.28068395E-01,-0.66178949E-02,-0.41906796E-02,
         0.26411945E-02, 0.19989067E-02,-0.18468504E-02,-0.74122488E-02,
         0.26923807E-02, 0.97303577E-02,-0.99370570E-03,-0.82478533E-02,
         0.24442570E-02,-0.26629157E-02, 0.20338618E-02,-0.78537880E-03,
        -0.89369379E-02,-0.26880408E-01,-0.17632040E-02,-0.50871670E-02,
        -0.73664859E-02, 0.44570349E-01, 0.13053592E-01,-0.17787386E-01,
        -0.23871791E-03,-0.64172298E-02,-0.74703866E-02,-0.29596246E-02,
         0.57698943E-03, 0.18889846E-02,
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
        +coeff[  3]            *x41    
        +coeff[  4]*x12                
        +coeff[  5]*x11*x21            
        +coeff[  6]    *x22            
        +coeff[  7]    *x21*x31        
    ;
    v_theta_l5p77                             =v_theta_l5p77                             
        +coeff[  8]    *x21    *x41    
        +coeff[  9]*x11*x22            
        +coeff[ 10]    *x23            
        +coeff[ 11]*x11*x21*x31        
        +coeff[ 12]    *x22*x31        
        +coeff[ 13]                *x51
        +coeff[ 14]    *x23*x31        
        +coeff[ 15]*x11            *x51
        +coeff[ 16]        *x31    *x51
    ;
    v_theta_l5p77                             =v_theta_l5p77                             
        +coeff[ 17]        *x31        
        +coeff[ 18]*x11    *x31        
        +coeff[ 19]    *x21*x32        
        +coeff[ 20]    *x21*x31*x41    
        +coeff[ 21]    *x21    *x42    
        +coeff[ 22]            *x41*x51
        +coeff[ 23]    *x22        *x51
        +coeff[ 24]        *x32    *x51
        +coeff[ 25]        *x32        
    ;
    v_theta_l5p77                             =v_theta_l5p77                             
        +coeff[ 26]*x11    *x32        
        +coeff[ 27]*x11*x21    *x41    
        +coeff[ 28]    *x22    *x41    
        +coeff[ 29]*x11*x22*x31        
        +coeff[ 30]*x11*x21        *x51
        +coeff[ 31]    *x22*x31    *x51
        +coeff[ 32]*x11        *x41    
        +coeff[ 33]        *x31*x41    
        +coeff[ 34]            *x42    
    ;
    v_theta_l5p77                             =v_theta_l5p77                             
        +coeff[ 35]*x12*x21            
        +coeff[ 36]*x11*x21*x32        
        +coeff[ 37]    *x22*x32        
        +coeff[ 38]    *x21*x33        
        +coeff[ 39]*x11*x22    *x41    
        +coeff[ 40]    *x23    *x41    
        +coeff[ 41]    *x22*x31*x41    
        +coeff[ 42]*x11*x21    *x42    
        +coeff[ 43]    *x22    *x42    
    ;
    v_theta_l5p77                             =v_theta_l5p77                             
        +coeff[ 44]    *x21        *x51
        +coeff[ 45]*x12*x21*x32        
        +coeff[ 46]*x11*x22*x32        
        +coeff[ 47]    *x21*x33*x41    
        +coeff[ 48]*x12            *x51
        +coeff[ 49]    *x21*x31    *x51
        ;

    return v_theta_l5p77                             ;
}
float phi_l5p77                               (float *x,int m){
    int ncoeff= 70;
    float avdat= -0.2389198E-02;
    float xmin[10]={
        -0.70676E+00,-0.24649E-01,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.20444E-01, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 71]={
         0.43710956E-03,-0.28982565E-01,-0.40501598E-02, 0.53884727E-02,
         0.30220046E-04, 0.97498903E-02,-0.12211645E-01, 0.15377986E-01,
         0.10519944E-01,-0.15544492E-01,-0.89674173E-02,-0.19325825E-02,
        -0.30355619E-01, 0.63038631E-02, 0.40771966E-02,-0.36201043E-04,
         0.60399296E-02, 0.16653709E-01,-0.19458594E-02,-0.68724044E-02,
         0.33892693E-02, 0.28309342E-01, 0.27564175E-02,-0.41526485E-01,
         0.80039341E-03,-0.12727531E-03,-0.14933942E-01,-0.10071074E-02,
        -0.45904610E-02,-0.47122300E-03, 0.31552629E-02, 0.26867124E-02,
        -0.57251425E-02, 0.55083267E-01,-0.13264680E-02,-0.77932016E-02,
        -0.15112128E-01, 0.45957183E-02, 0.94947629E-02,-0.42963950E-02,
         0.15736868E-02, 0.21331001E-01, 0.23997135E-02,-0.30180311E-01,
         0.83291689E-02, 0.55938396E-02,-0.83683822E-02, 0.61003454E-02,
         0.13964762E-02, 0.22550174E-02,-0.18014201E-02,-0.70316783E-04,
        -0.45486097E-02, 0.58491663E-02,-0.43193396E-03,-0.18359589E-02,
        -0.20903321E-02, 0.83398953E-03,-0.15164639E-01, 0.69670123E-02,
        -0.12749496E-02,-0.37496660E-02,-0.88142483E-02,-0.43048211E-02,
         0.86755231E-02, 0.32005962E-02, 0.20387203E-02,-0.20084667E-02,
         0.87503844E-03,-0.11048340E-02,
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
        +coeff[ 10]    *x22            
        +coeff[ 11]*x12    *x31        
        +coeff[ 12]    *x22*x31        
        +coeff[ 13]*x12        *x41    
        +coeff[ 14]    *x24            
        +coeff[ 15]                *x51
        +coeff[ 16]    *x21        *x51
    ;
    v_phi_l5p77                               =v_phi_l5p77                               
        +coeff[ 17]    *x21*x31    *x51
        +coeff[ 18]*x11*x21            
        +coeff[ 19]*x11*x21*x31        
        +coeff[ 20]*x11*x21    *x41    
        +coeff[ 21]    *x22    *x41    
        +coeff[ 22]*x11*x22            
        +coeff[ 23]    *x22*x32        
        +coeff[ 24]*x11            *x51
        +coeff[ 25]        *x31*x41*x51
    ;
    v_phi_l5p77                               =v_phi_l5p77                               
        +coeff[ 26]    *x21    *x41*x51
        +coeff[ 27]*x11*x21        *x51
        +coeff[ 28]    *x23        *x51
        +coeff[ 29]*x11    *x33    *x51
        +coeff[ 30]*x11    *x32        
        +coeff[ 31]    *x21*x32        
        +coeff[ 32]*x11*x21*x32        
        +coeff[ 33]    *x22*x31*x41    
        +coeff[ 34]*x13    *x31        
    ;
    v_phi_l5p77                               =v_phi_l5p77                               
        +coeff[ 35]    *x23*x31        
        +coeff[ 36]*x11*x22    *x41    
        +coeff[ 37]    *x23    *x41    
        +coeff[ 38]    *x24*x31        
        +coeff[ 39]*x11*x24*x31        
        +coeff[ 40]*x11    *x31    *x51
        +coeff[ 41]    *x21*x32    *x51
        +coeff[ 42]*x11    *x31*x41*x51
        +coeff[ 43]    *x21*x31*x41*x51
    ;
    v_phi_l5p77                               =v_phi_l5p77                               
        +coeff[ 44]    *x21    *x42*x51
        +coeff[ 45]*x11*x21    *x41*x51
        +coeff[ 46]    *x23*x31    *x51
        +coeff[ 47]    *x22*x33    *x51
        +coeff[ 48]    *x23*x35*x42    
        +coeff[ 49]            *x42    
        +coeff[ 50]    *x21    *x41    
        +coeff[ 51]*x12                
        +coeff[ 52]        *x33        
    ;
    v_phi_l5p77                               =v_phi_l5p77                               
        +coeff[ 53]        *x32*x41    
        +coeff[ 54]    *x23            
        +coeff[ 55]*x11    *x33        
        +coeff[ 56]    *x21*x33        
        +coeff[ 57]*x11*x21    *x42    
        +coeff[ 58]    *x22    *x42    
        +coeff[ 59]*x11*x22*x31        
        +coeff[ 60]*x12*x22            
        +coeff[ 61]    *x22*x33        
    ;
    v_phi_l5p77                               =v_phi_l5p77                               
        +coeff[ 62]    *x23*x32        
        +coeff[ 63]*x11*x22*x31*x41    
        +coeff[ 64]    *x23*x31*x41    
        +coeff[ 65]*x11*x23*x31        
        +coeff[ 66]*x14        *x41    
        +coeff[ 67]*x13*x21    *x41    
        +coeff[ 68]        *x31    *x51
        +coeff[ 69]            *x41*x51
        ;

    return v_phi_l5p77                               ;
}
float y00_l5p77                               (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.8433543E-03;
    float xmin[10]={
        -0.70676E+00,-0.24649E-01,-0.64821E-01,-0.48367E-01,-0.14986E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61008E+00, 0.20444E-01, 0.49615E-01, 0.41878E-01, 0.14968E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
         0.56533720E-02,-0.47242589E-01, 0.84736817E-01,-0.44322810E-02,
         0.34817269E-02,-0.16404880E-01, 0.20632038E-01,-0.56160716E-02,
        -0.48149325E-01,-0.50194129E-01, 0.13796888E-01, 0.18094197E-01,
         0.38736613E-03, 0.34246158E-02, 0.10469529E-01, 0.23417281E-01,
        -0.53168170E-01, 0.50643288E-01,-0.22584628E-01, 0.26624396E-02,
         0.99628814E-03,-0.89699682E-02, 0.12552742E-02, 0.14637413E-01,
         0.51937472E-01,-0.12898319E-01,-0.17032564E-01,-0.32416537E-01,
        -0.52508009E-02, 0.18329092E-04,-0.16808581E-01, 0.16552359E-01,
         0.11341276E-01, 0.54904398E-01, 0.13823218E-01, 0.38047207E-02,
        -0.12366975E-02, 0.14133148E-02,-0.73240204E-02,-0.22297854E-01,
         0.11712445E-02, 0.19627826E-01, 0.19356111E-02,-0.34878405E-02,
        -0.22597721E-01, 0.22142373E-01, 0.48253727E-02,-0.76237557E-04,
        -0.52351383E-02, 0.76057301E-02, 0.39067497E-02,-0.58884077E-01,
         0.65571233E-02, 0.98300623E-02,-0.57457155E-03,-0.84202533E-03,
        -0.28901470E-02,-0.20488648E-01, 0.16783204E-01,-0.52032107E-02,
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
        +coeff[ 19]    *x21*x32        
        +coeff[ 20]*x11    *x31*x41    
        +coeff[ 21]    *x21*x31*x41    
        +coeff[ 22]*x12    *x31        
        +coeff[ 23]*x11*x21*x31        
        +coeff[ 24]    *x22*x31        
        +coeff[ 25]*x12        *x41    
    ;
    v_y00_l5p77                               =v_y00_l5p77                               
        +coeff[ 26]*x11*x21    *x41    
        +coeff[ 27]    *x22    *x41    
        +coeff[ 28]*x11*x22            
        +coeff[ 29]                *x51
        +coeff[ 30]*x11    *x31*x42    
        +coeff[ 31]*x11        *x43    
        +coeff[ 32]*x11*x21*x32        
        +coeff[ 33]    *x22*x32        
        +coeff[ 34]    *x22    *x42    
    ;
    v_y00_l5p77                               =v_y00_l5p77                               
        +coeff[ 35]*x13    *x31        
        +coeff[ 36]        *x31    *x51
        +coeff[ 37]*x12*x22            
        +coeff[ 38]    *x21        *x51
        +coeff[ 39]    *x21*x31    *x51
        +coeff[ 40]*x11        *x41*x51
        +coeff[ 41]    *x21    *x41*x51
        +coeff[ 42]*x11*x21        *x51
        +coeff[ 43]*x11    *x32    *x51
    ;
    v_y00_l5p77                               =v_y00_l5p77                               
        +coeff[ 44]    *x21*x32    *x51
        +coeff[ 45]    *x21*x31*x41*x51
        +coeff[ 46]*x12*x23*x31        
        +coeff[ 47]*x13            *x51
        +coeff[ 48]*x11    *x32        
        +coeff[ 49]*x11    *x33        
        +coeff[ 50]    *x21*x33        
        +coeff[ 51]    *x22*x31*x41    
        +coeff[ 52]    *x23*x31        
    ;
    v_y00_l5p77                               =v_y00_l5p77                               
        +coeff[ 53]*x11*x22    *x41    
        +coeff[ 54]*x11*x23            
        +coeff[ 55]*x11            *x51
        +coeff[ 56]*x11    *x31    *x51
        +coeff[ 57]*x13*x22*x31        
        +coeff[ 58]*x13*x22    *x41    
        +coeff[ 59]*x11*x21    *x41*x51
        ;

    return v_y00_l5p77                               ;
}
}
