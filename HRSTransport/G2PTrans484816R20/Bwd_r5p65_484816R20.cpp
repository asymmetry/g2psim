#include "Bwd_r5p65_484816R20.h"

namespace S484816R20
{
    float txfit_r5p65                             (float *x,int m){
        int ncoeff=  2;
        float avdat=  0.3558419E-02;
        float xmin[10]={
            -0.67894E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
            0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
            float xmax[10]={
                0.60087E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
                0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
                float scale[10]={0};
                float coeff[  3]={
                    -0.10712370E-01, 0.11366472E+00,
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
     
                    float v_txfit_r5p65                             =avdat
                                +coeff[  0]    
                                +coeff[  1]*x11
                                ;
     
                    return v_txfit_r5p65                             ;
    }
     
     
    float delta_r5p65                             (float *x,int m){
        int ncoeff= 21;
        float avdat=  0.1900347E-02;
        float xmin[10]={
            -0.70610E+00,-0.24697E-01,-0.50228E-01,-0.47683E-01,-0.20000E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.60639E+00, 0.20490E-01, 0.48409E-01, 0.33001E-01, 0.19967E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0};
        float coeff[ 22]={
            -0.55465251E-02, 0.47821123E-01, 0.44066645E-02, 0.42855269E-02,
            -0.13069338E-02, 0.31658183E-02,-0.32333382E-02, 0.59591464E-04,
             0.50774684E-04, 0.37877620E-03,-0.22398306E-03,-0.62973198E-03,
             0.16217638E-02, 0.80424116E-03,-0.12591368E-02, 0.84210047E-03,
            -0.17070004E-02, 0.16723020E-02, 0.36748691E-03,-0.85193518E-03,
            -0.10496072E-03,
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
        float x33 = x32*x3;
        float x34 = x33*x3;
        float x41 = x4;
        float x42 = x41*x4;
        float x51 = x5;
     
    //                 function
     
        float v_delta_r5p65                             =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]*x12                
            +coeff[  3]*x11*x21            
            +coeff[  4]    *x22            
            +coeff[  5]    *x21*x31        
            +coeff[  6]    *x21*x32        
            +coeff[  7]        *x34        
        ;
        v_delta_r5p65                             =v_delta_r5p65                             
            +coeff[  8]        *x33*x41    
            +coeff[  9]                *x51
            +coeff[ 10]        *x31        
            +coeff[ 11]            *x41    
            +coeff[ 12]        *x32        
            +coeff[ 13]*x11*x22            
            +coeff[ 14]*x11*x21*x31        
            +coeff[ 15]*x11    *x31        
            +coeff[ 16]        *x31*x41    
        ;
        v_delta_r5p65                             =v_delta_r5p65                             
            +coeff[ 17]            *x42    
            +coeff[ 18]*x13                
            +coeff[ 19]*x11    *x32        
            +coeff[ 20]        *x31    *x51
            ;
     
        return v_delta_r5p65                             ;
    }
    float theta_r5p65                             (float *x,int m){
        int ncoeff= 50;
        float avdat= -0.6241306E-03;
        float xmin[10]={
            -0.70610E+00,-0.24697E-01,-0.50228E-01,-0.47683E-01,-0.20000E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.60639E+00, 0.20490E-01, 0.48409E-01, 0.33001E-01, 0.19967E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0};
        float coeff[ 51]={
             0.72575402E-02,-0.77559552E-02,-0.51556468E-01, 0.63029169E-02,
            -0.39994782E-02, 0.26576439E-01,-0.16576106E-01,-0.48418245E-02,
             0.68669597E-03, 0.87551301E-03,-0.12988755E-01,-0.39049045E-02,
            -0.10971467E-03,-0.21243852E-02, 0.29677469E-02,-0.33116448E-02,
            -0.12798026E-01,-0.15659145E-02, 0.51188380E-02, 0.75786696E-02,
            -0.14400104E-02, 0.23212556E-02,-0.98654376E-02,-0.10797083E-02,
             0.47586067E-02,-0.72340906E-03, 0.47300849E-02, 0.11943238E-01,
            -0.94826886E-03,-0.82797660E-02,-0.83205393E-02, 0.81339572E-02,
             0.14726567E-01,-0.30333095E-02, 0.66413931E-02,-0.34368769E-03,
             0.15164766E-01,-0.55514988E-02, 0.46647305E-03, 0.28422652E-03,
            -0.57497050E-03,-0.56520698E-03,-0.41559819E-03, 0.89346860E-02,
            -0.43455847E-02, 0.47203433E-02, 0.52603181E-02,-0.11324953E-01,
            -0.61002146E-02,-0.22634456E-02,
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
        float x51 = x5;
     
    //                 function
     
        float v_theta_r5p65                             =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]*x11*x21            
            +coeff[  4]    *x22            
            +coeff[  5]    *x21*x31        
            +coeff[  6]    *x21    *x41    
            +coeff[  7]        *x31*x41    
        ;
        v_theta_r5p65                             =v_theta_r5p65                             
            +coeff[  8]*x11    *x32        
            +coeff[  9]                *x51
            +coeff[ 10]*x12    *x31*x41    
            +coeff[ 11]        *x31        
            +coeff[ 12]            *x41    
            +coeff[ 13]*x12                
            +coeff[ 14]        *x32        
            +coeff[ 15]*x11*x21*x31        
            +coeff[ 16]    *x21*x32        
        ;
        v_theta_r5p65                             =v_theta_r5p65                             
            +coeff[ 17]*x11*x21    *x41    
            +coeff[ 18]*x11    *x31        
            +coeff[ 19]            *x42    
            +coeff[ 20]*x12*x21            
            +coeff[ 21]*x11*x22            
            +coeff[ 22]    *x22    *x41    
            +coeff[ 23]    *x21*x31*x41    
            +coeff[ 24]    *x21*x33        
            +coeff[ 25]        *x31    *x51
        ;
        v_theta_r5p65                             =v_theta_r5p65                             
            +coeff[ 26]    *x23            
            +coeff[ 27]    *x22*x31        
            +coeff[ 28]*x12        *x41    
            +coeff[ 29]*x11    *x31*x41    
            +coeff[ 30]    *x23*x31        
            +coeff[ 31]*x12    *x32        
            +coeff[ 32]*x11*x21*x32        
            +coeff[ 33]*x11*x21*x31*x41    
            +coeff[ 34]*x12        *x42    
        ;
        v_theta_r5p65                             =v_theta_r5p65                             
            +coeff[ 35]*x11            *x51
            +coeff[ 36]*x12*x21*x31*x41    
            +coeff[ 37]*x12*x21    *x42    
            +coeff[ 38]        *x32    *x51
            +coeff[ 39]*x13                
            +coeff[ 40]        *x33        
            +coeff[ 41]    *x21    *x42    
            +coeff[ 42]*x12*x22            
            +coeff[ 43]*x12*x21*x31        
        ;
        v_theta_r5p65                             =v_theta_r5p65                             
            +coeff[ 44]*x12*x21    *x41    
            +coeff[ 45]*x11*x22    *x41    
            +coeff[ 46]    *x22*x31*x41    
            +coeff[ 47]*x12*x21*x32        
            +coeff[ 48]*x11*x22*x32        
            +coeff[ 49]    *x21*x33*x41    
            ;
     
        return v_theta_r5p65                             ;
    }
    float phi_r5p65                               (float *x,int m){
        int ncoeff= 48;
        float avdat= -0.3636865E-03;
        float xmin[10]={
            -0.70610E+00,-0.24697E-01,-0.50228E-01,-0.47683E-01,-0.20000E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.60639E+00, 0.20490E-01, 0.48409E-01, 0.33001E-01, 0.19967E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0};
        float coeff[ 49]={
             0.17570065E-02,-0.23052601E-01,-0.56504733E-02,-0.34306443E-02,
            -0.37547320E-02, 0.44382755E-02,-0.27597049E-01, 0.23992192E-01,
             0.12921245E-01, 0.14384366E-01,-0.14404798E-01,-0.62213428E-02,
             0.98849032E-02,-0.28347704E-02, 0.34903183E-02, 0.44693379E-02,
            -0.83015644E-03, 0.29754038E-02,-0.68466594E-02, 0.11031725E-01,
             0.78285960E-02, 0.10616847E-01,-0.10481354E-01,-0.26658332E-01,
             0.11528999E-01, 0.47951713E-02,-0.10212207E-01, 0.15092775E-02,
            -0.17071633E-02,-0.76694572E-02,-0.69685635E-03,-0.30132022E-02,
             0.30866736E-02, 0.26304638E-02,-0.79994131E-03,-0.14833496E-03,
            -0.33574447E-02, 0.79165082E-02,-0.43209707E-02, 0.20707939E-02,
            -0.25839768E-02,-0.14420595E-02,-0.84233627E-03, 0.73129719E-03,
             0.12339714E-03, 0.28012057E-02, 0.65540220E-03,-0.31889719E-03,
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
        float x51 = x5;
     
    //                 function
     
        float v_phi_r5p65                               =avdat
            +coeff[  0]                    
            +coeff[  1]        *x31        
            +coeff[  2]            *x41    
            +coeff[  3]*x11                
            +coeff[  4]    *x21            
            +coeff[  5]        *x32        
            +coeff[  6]        *x31*x41    
            +coeff[  7]            *x42    
        ;
        v_phi_r5p65                               =v_phi_r5p65                               
            +coeff[  8]*x11    *x31        
            +coeff[  9]    *x21*x31        
            +coeff[ 10]*x11        *x41    
            +coeff[ 11]    *x21    *x41    
            +coeff[ 12]    *x22            
            +coeff[ 13]*x11*x22            
            +coeff[ 14]    *x23            
            +coeff[ 15]*x11    *x31*x42    
            +coeff[ 16]    *x21        *x51
        ;
        v_phi_r5p65                               =v_phi_r5p65                               
            +coeff[ 17]*x11*x21            
            +coeff[ 18]        *x33        
            +coeff[ 19]        *x32*x41    
            +coeff[ 20]*x11    *x32        
            +coeff[ 21]    *x21*x32        
            +coeff[ 22]*x11    *x31*x41    
            +coeff[ 23]    *x21*x31*x41    
            +coeff[ 24]    *x21    *x42    
            +coeff[ 25]*x12    *x31        
        ;
        v_phi_r5p65                               =v_phi_r5p65                               
            +coeff[ 26]    *x22*x31        
            +coeff[ 27]*x12*x21            
            +coeff[ 28]*x11    *x33        
            +coeff[ 29]*x11    *x32*x41    
            +coeff[ 30]*x14                
            +coeff[ 31]*x11*x23            
            +coeff[ 32]*x13    *x32        
            +coeff[ 33]*x14        *x41    
            +coeff[ 34]*x12                
        ;
        v_phi_r5p65                               =v_phi_r5p65                               
            +coeff[ 35]*x11*x21*x31        
            +coeff[ 36]    *x21*x33        
            +coeff[ 37]    *x22*x32        
            +coeff[ 38]*x12    *x31*x41    
            +coeff[ 39]*x11*x21*x31*x41    
            +coeff[ 40]    *x23*x31        
            +coeff[ 41]*x13        *x41    
            +coeff[ 42]*x13*x21            
            +coeff[ 43]*x12*x22            
        ;
        v_phi_r5p65                               =v_phi_r5p65                               
            +coeff[ 44]                *x51
            +coeff[ 45]*x11*x23*x32        
            +coeff[ 46]    *x21*x31    *x51
            +coeff[ 47]    *x22        *x51
            ;
     
        return v_phi_r5p65                               ;
    }
    float y00_r5p65                               (float *x,int m){
        int ncoeff= 60;
        float avdat=  0.3318688E-04;
        float xmin[10]={
            -0.70610E+00,-0.24697E-01,-0.50228E-01,-0.47683E-01,-0.20000E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.60639E+00, 0.20490E-01, 0.48409E-01, 0.33001E-01, 0.19967E-02,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0};
        float coeff[ 61]={
            -0.72100614E-02,-0.26094975E-01, 0.39174877E-01, 0.82583766E-03,
            -0.10656346E-01,-0.44135146E-01,-0.42671315E-01, 0.33389318E-02,
             0.47339130E-01, 0.29388091E-01,-0.34760959E-01, 0.55737019E-01,
             0.40146140E-02,-0.76450303E-03,-0.16396437E-02,-0.35327359E-03,
            -0.24496207E-01,-0.11102409E-01, 0.25352066E-01, 0.25144916E-02,
             0.32573449E-02, 0.47549335E-02, 0.13709274E-01,-0.12966079E-01,
            -0.15894126E-01,-0.15847275E-01,-0.14954400E-01,-0.28869556E-01,
             0.13274276E-02, 0.51305041E-01,-0.47020137E-01, 0.12306320E-02,
            -0.96451881E-03, 0.44393013E-02, 0.36638693E-02,-0.46246327E-02,
             0.27123336E-02, 0.12868509E-01, 0.24467986E-01, 0.47500581E-02,
             0.31496848E-02,-0.57433564E-02, 0.61532976E-02,-0.44876882E-02,
            -0.36370193E-02, 0.33094953E-02,-0.20784924E-01, 0.19795083E-01,
             0.82760174E-02,-0.83090309E-02,-0.18497730E-01,-0.21003305E-02,
             0.72054295E-02, 0.86920466E-02,-0.13078620E-02,-0.16447081E-02,
            -0.14784225E-02, 0.52982899E-02, 0.42168493E-02, 0.30432774E-02,
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
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
     
    //                 function
     
        float v_y00_r5p65                               =avdat
            +coeff[  0]                    
            +coeff[  1]        *x31        
            +coeff[  2]            *x41    
            +coeff[  3]*x11*x21            
            +coeff[  4]    *x22            
            +coeff[  5]*x11    *x31        
            +coeff[  6]    *x21*x31        
            +coeff[  7]        *x32        
        ;
        v_y00_r5p65                               =v_y00_r5p65                               
            +coeff[  8]*x11        *x41    
            +coeff[  9]    *x21    *x41    
            +coeff[ 10]        *x31*x41    
            +coeff[ 11]            *x42    
            +coeff[ 12]*x13                
            +coeff[ 13]*x12*x21            
            +coeff[ 14]*x11*x22            
            +coeff[ 15]    *x23            
            +coeff[ 16]*x12    *x31        
        ;
        v_y00_r5p65                               =v_y00_r5p65                               
            +coeff[ 17]*x11*x21*x31        
            +coeff[ 18]    *x22*x31        
            +coeff[ 19]*x11    *x32        
            +coeff[ 20]    *x21*x32        
            +coeff[ 21]        *x33        
            +coeff[ 22]*x12        *x41    
            +coeff[ 23]*x11*x21    *x41    
            +coeff[ 24]    *x22    *x41    
            +coeff[ 25]    *x21*x31*x41    
        ;
        v_y00_r5p65                               =v_y00_r5p65                               
            +coeff[ 26]        *x32*x41    
            +coeff[ 27]*x11        *x42    
            +coeff[ 28]    *x21    *x42    
            +coeff[ 29]        *x31*x42    
            +coeff[ 30]            *x43    
            +coeff[ 31]*x14                
            +coeff[ 32]*x13*x21            
            +coeff[ 33]*x11*x23            
            +coeff[ 34]    *x24            
        ;
        v_y00_r5p65                               =v_y00_r5p65                               
            +coeff[ 35]*x13    *x31        
            +coeff[ 36]*x12*x21*x31        
            +coeff[ 37]*x11*x22*x31        
            +coeff[ 38]    *x23*x31        
            +coeff[ 39]*x12    *x32        
            +coeff[ 40]*x11*x21*x32        
            +coeff[ 41]    *x22*x32        
            +coeff[ 42]        *x34        
            +coeff[ 43]*x13        *x41    
        ;
        v_y00_r5p65                               =v_y00_r5p65                               
            +coeff[ 44]*x12*x21    *x41    
            +coeff[ 45]*x11*x22    *x41    
            +coeff[ 46]    *x23    *x41    
            +coeff[ 47]*x12    *x31*x41    
            +coeff[ 48]*x11*x21*x31*x41    
            +coeff[ 49]        *x33*x41    
            +coeff[ 50]*x12        *x42    
            +coeff[ 51]    *x22    *x42    
            +coeff[ 52]    *x21*x31*x42    
        ;
        v_y00_r5p65                               =v_y00_r5p65                               
            +coeff[ 53]*x11        *x43    
            +coeff[ 54]*x13*x22            
            +coeff[ 55]*x12*x23            
            +coeff[ 56]*x11*x24            
            +coeff[ 57]*x13*x21*x31        
            +coeff[ 58]*x12*x22*x31        
            +coeff[ 59]    *x24*x31        
            ;
     
        return v_y00_r5p65                               ;
    }
}
