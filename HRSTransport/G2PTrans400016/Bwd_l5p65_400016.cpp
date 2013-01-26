#include "Bwd_l5p65_400016.h"

namespace S400016
{
    float txfit_l5p65                     (float *x,int m){
        int ncoeff=  2;
        float avdat=  0.3074509E-02;
        float xmin[10]={
            -0.68262E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.62245E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[  3]={
            -0.77049020E-02, 0.11616763E+00,
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
     
        float v_txfit_l5p65                     =avdat
            +coeff[  0]    
            +coeff[  1]*x11
            ;
     
        return v_txfit_l5p65                     ;
    }
    float delta_l5p65                     (float *x,int m){
        int ncoeff= 27;
        float avdat=  0.2636606E-02;
        float xmin[10]={
            -0.69130E+00,-0.25676E-01,-0.54269E-01,-0.41176E-01,-0.14985E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.62245E+00, 0.19559E-01, 0.56361E-01, 0.49778E-01, 0.14999E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 28]={
            -0.51159216E-02, 0.48094507E-01, 0.44774259E-02, 0.36355774E-02,
             0.26125412E-02,-0.17138434E-02,-0.48939181E-02, 0.36582286E-02,
             0.50390488E-03, 0.13136938E-02, 0.45609544E-03,-0.55146968E-03,
             0.29546588E-02,-0.42893556E-02,-0.29436182E-02, 0.56111661E-03,
             0.26434993E-02, 0.41622904E-03, 0.47739499E-03,-0.16741237E-02,
            -0.95380720E-03, 0.12749824E-02,-0.57398021E-03, 0.78731927E-03,
            -0.75228838E-03, 0.42600371E-03, 0.18956560E-03,
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
     
        float v_delta_l5p65                     =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]*x12                
            +coeff[  3]*x11*x21            
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]    *x21*x31        
            +coeff[  7]    *x21    *x41    
        ;
        v_delta_l5p65                     =v_delta_l5p65                     
            +coeff[  8]    *x21        *x51
            +coeff[  9]        *x31    *x51
            +coeff[ 10]        *x31        
            +coeff[ 11]*x11    *x31        
            +coeff[ 12]        *x32        
            +coeff[ 13]        *x31*x41    
            +coeff[ 14]    *x21*x32        
            +coeff[ 15]    *x21            
            +coeff[ 16]            *x42    
        ;
        v_delta_l5p65                     =v_delta_l5p65                     
            +coeff[ 17]*x13                
            +coeff[ 18]*x11*x22            
            +coeff[ 19]    *x22*x31        
            +coeff[ 20]*x11*x21    *x41    
            +coeff[ 21]    *x22    *x41    
            +coeff[ 22]*x11    *x31*x41    
            +coeff[ 23]    *x24            
            +coeff[ 24]            *x41*x51
            +coeff[ 25]        *x32    *x51
        ;
        v_delta_l5p65                     =v_delta_l5p65                     
            +coeff[ 26]*x13            *x51
            ;
     
        return v_delta_l5p65                     ;
    }
    float theta_l5p65                     (float *x,int m){
        int ncoeff= 50;
        float avdat=  0.7915922E-04;
        float xmin[10]={
            -0.69130E+00,-0.25676E-01,-0.54269E-01,-0.41176E-01,-0.14985E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.62245E+00, 0.19559E-01, 0.56361E-01, 0.49778E-01, 0.14999E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 51]={
             0.84537389E-02,-0.76094363E-02,-0.52707892E-01,-0.24669352E-02,
             0.37432837E-02,-0.50889612E-02,-0.30187258E-01, 0.26163066E-01,
             0.57187206E-02,-0.33168094E-02, 0.68692784E-02, 0.95498972E-02,
             0.70728483E-02,-0.63399398E-02,-0.35289838E-02, 0.22485794E-02,
            -0.20079708E-01,-0.26718237E-01,-0.63467957E-02, 0.14572031E-01,
             0.38108170E-01,-0.15969982E-01,-0.10955429E-02,-0.67428346E-02,
             0.21694165E-02, 0.71928478E-02, 0.43697190E-02,-0.97615598E-02,
            -0.76353748E-03, 0.12820273E-01,-0.27778402E-02,-0.48252959E-02,
             0.17583586E-02,-0.71269586E-02, 0.70343274E-02,-0.67543001E-02,
             0.34273462E-02,-0.88134623E-03, 0.48838404E-03,-0.80286304E-03,
             0.53031333E-02,-0.19764565E-01,-0.74797794E-02,-0.69250627E-02,
             0.22648960E-01, 0.69363732E-02, 0.12550668E-02,-0.58358638E-02,
            -0.96223160E-03, 0.19361848E-02,
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
     
        float v_theta_l5p65                     =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]*x12                
            +coeff[  4]*x11*x21            
            +coeff[  5]    *x22            
            +coeff[  6]    *x21*x31        
            +coeff[  7]    *x21    *x41    
        ;
        v_theta_l5p65                     =v_theta_l5p65                     
            +coeff[  8]    *x23            
            +coeff[  9]*x11*x21*x31        
            +coeff[ 10]                *x51
            +coeff[ 11]        *x31    *x51
            +coeff[ 12]        *x31        
            +coeff[ 13]            *x41    
            +coeff[ 14]*x11    *x31        
            +coeff[ 15]*x11*x22            
            +coeff[ 16]    *x22*x31        
        ;
        v_theta_l5p65                     =v_theta_l5p65                     
            +coeff[ 17]    *x21*x32        
            +coeff[ 18]*x11*x21    *x41    
            +coeff[ 19]    *x22    *x41    
            +coeff[ 20]    *x21*x31*x41    
            +coeff[ 21]    *x21    *x42    
            +coeff[ 22]*x11            *x51
            +coeff[ 23]            *x41*x51
            +coeff[ 24]        *x32    *x51
            +coeff[ 25]        *x32        
        ;
        v_theta_l5p65                     =v_theta_l5p65                     
            +coeff[ 26]*x11        *x41    
            +coeff[ 27]        *x31*x41    
            +coeff[ 28]*x11    *x32        
            +coeff[ 29]    *x23*x31        
            +coeff[ 30]    *x21*x33        
            +coeff[ 31]    *x22        *x51
            +coeff[ 32]*x11        *x41*x51
            +coeff[ 33]    *x22*x31    *x51
            +coeff[ 34]    *x23*x31    *x51
        ;
        v_theta_l5p65                     =v_theta_l5p65                     
            +coeff[ 35]    *x23    *x41*x51
            +coeff[ 36]            *x42    
            +coeff[ 37]*x12*x21            
            +coeff[ 38]        *x33        
            +coeff[ 39]*x12        *x41    
            +coeff[ 40]*x11*x22*x31        
            +coeff[ 41]    *x22*x32        
            +coeff[ 42]*x11*x22    *x41    
            +coeff[ 43]*x11*x21*x31*x41    
        ;
        v_theta_l5p65                     =v_theta_l5p65                     
            +coeff[ 44]    *x22*x31*x41    
            +coeff[ 45]*x11*x21    *x42    
            +coeff[ 46]    *x21        *x51
            +coeff[ 47]*x11*x22*x32        
            +coeff[ 48]*x11*x21        *x51
            +coeff[ 49]    *x21*x31    *x51
            ;
     
        return v_theta_l5p65                     ;
    }
    float phi_l5p65                       (float *x,int m){
        int ncoeff= 70;
        float avdat=  0.1725391E-02;
        float xmin[10]={
            -0.69130E+00,-0.25676E-01,-0.54269E-01,-0.41176E-01,-0.14985E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.62245E+00, 0.19559E-01, 0.56361E-01, 0.49778E-01, 0.14999E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 71]={
            -0.28719886E-02,-0.28119942E-01,-0.32492650E-02, 0.53481581E-02,
             0.41905609E-02, 0.64109219E-02,-0.73789312E-02, 0.22297373E-03,
             0.15141395E-01, 0.13943932E-01,-0.15352657E-01,-0.86096162E-02,
            -0.16620832E-02,-0.24622740E-01, 0.64157797E-02, 0.13201669E-01,
            -0.25650826E-02,-0.32223773E-02, 0.49820938E-03,-0.37835583E-02,
             0.41463152E-02,-0.13196831E-03, 0.20157952E-01, 0.19502009E-02,
            -0.17664387E-02,-0.93547686E-03,-0.12131199E-01,-0.89024380E-02,
             0.53401017E-02, 0.34555460E-02, 0.13498546E-02,-0.98183649E-02,
            -0.13629668E-02, 0.92929238E-02,-0.58581461E-02, 0.61634998E-02,
            -0.10088404E-01,-0.12470415E-02, 0.18607294E-03, 0.61641731E-02,
             0.21695171E-02, 0.72429460E-02,-0.21329373E-02,-0.17132909E-02,
             0.23517902E-03,-0.11361873E-01, 0.34336790E-01,-0.39383590E-01,
             0.15556064E-01,-0.15868390E-01, 0.15079169E-01, 0.11013871E-02,
             0.13039842E-03,-0.36964577E-03,-0.63947206E-02, 0.88999253E-02,
            -0.16331910E-02, 0.57744808E-02,-0.15364131E-02,-0.30942534E-02,
            -0.78414902E-02,-0.55066925E-02,-0.78885898E-03,-0.10100672E-02,
             0.39632008E-02,-0.79565163E-03, 0.11015135E-02, 0.13578867E-02,
            -0.31078155E-02, 0.18448508E-02,
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
        float x34 = x33*x3;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x51 = x5;
     
    //                 function
     
        float v_phi_l5p65                       =avdat
            +coeff[  0]                    
            +coeff[  1]        *x31        
            +coeff[  2]            *x41    
            +coeff[  3]*x11                
            +coeff[  4]    *x21            
            +coeff[  5]        *x32        
            +coeff[  6]        *x31*x41    
            +coeff[  7]            *x42    
        ;
        v_phi_l5p65                       =v_phi_l5p65                       
            +coeff[  8]*x11    *x31        
            +coeff[  9]    *x21*x31        
            +coeff[ 10]*x11        *x41    
            +coeff[ 11]    *x22            
            +coeff[ 12]*x12    *x31        
            +coeff[ 13]    *x22*x31        
            +coeff[ 14]    *x21        *x51
            +coeff[ 15]    *x21*x31    *x51
            +coeff[ 16]*x12*x23            
        ;
        v_phi_l5p65                       =v_phi_l5p65                       
            +coeff[ 17]    *x21    *x41    
            +coeff[ 18]*x12                
            +coeff[ 19]*x11*x21            
            +coeff[ 20]*x11    *x31*x41    
            +coeff[ 21]*x11        *x42    
            +coeff[ 22]    *x22    *x41    
            +coeff[ 23]*x11*x22            
            +coeff[ 24]    *x23            
            +coeff[ 25]                *x51
        ;
        v_phi_l5p65                       =v_phi_l5p65                       
            +coeff[ 26]    *x22*x32        
            +coeff[ 27]    *x23*x31        
            +coeff[ 28]    *x23    *x41    
            +coeff[ 29]*x11*x23            
            +coeff[ 30]*x11            *x51
            +coeff[ 31]    *x21    *x41*x51
            +coeff[ 32]*x11*x21        *x51
            +coeff[ 33]    *x21*x32        
            +coeff[ 34]*x11*x21*x31        
        ;
        v_phi_l5p65                       =v_phi_l5p65                       
            +coeff[ 35]*x12        *x41    
            +coeff[ 36]*x11*x21    *x42    
            +coeff[ 37]        *x31    *x51
            +coeff[ 38]            *x41*x51
            +coeff[ 39]*x11*x22    *x42    
            +coeff[ 40]*x11    *x31    *x51
            +coeff[ 41]    *x21*x32    *x51
            +coeff[ 42]*x11*x22        *x51
            +coeff[ 43]    *x23        *x51
        ;
        v_phi_l5p65                       =v_phi_l5p65                       
            +coeff[ 44]        *x34    *x51
            +coeff[ 45]        *x33        
            +coeff[ 46]        *x32*x41    
            +coeff[ 47]        *x31*x42    
            +coeff[ 48]            *x43    
            +coeff[ 49]    *x21*x31*x41    
            +coeff[ 50]    *x21    *x42    
            +coeff[ 51]*x11*x21    *x41    
            +coeff[ 52]*x12*x21            
        ;
        v_phi_l5p65                       =v_phi_l5p65                       
            +coeff[ 53]        *x33*x41    
            +coeff[ 54]*x11    *x33        
            +coeff[ 55]*x11    *x32*x41    
            +coeff[ 56]*x11        *x43    
            +coeff[ 57]    *x21    *x43    
            +coeff[ 58]*x12*x21*x31        
            +coeff[ 59]*x11*x22    *x41    
            +coeff[ 60]    *x22*x31*x42    
            +coeff[ 61]*x11*x21    *x43    
        ;
        v_phi_l5p65                       =v_phi_l5p65                       
            +coeff[ 62]        *x32    *x51
            +coeff[ 63]            *x42*x51
            +coeff[ 64]*x11*x23*x31        
            +coeff[ 65]*x11        *x41*x51
            +coeff[ 66]    *x22        *x51
            +coeff[ 67]*x11    *x32    *x51
            +coeff[ 68]    *x21*x31*x41*x51
            +coeff[ 69]    *x22*x31    *x51
            ;
     
        return v_phi_l5p65                       ;
    }
    float y00_l5p65                       (float *x,int m){
        int ncoeff= 60;
        float avdat=  0.7328165E-03;
        float xmin[10]={
            -0.69130E+00,-0.25676E-01,-0.54269E-01,-0.41176E-01,-0.14985E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.62245E+00, 0.19559E-01, 0.56361E-01, 0.49778E-01, 0.14999E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 61]={
             0.53526638E-02,-0.47400024E-01, 0.84597647E-01,-0.10292803E-01,
            -0.25778611E-02,-0.85732685E-02,-0.44427224E-01,-0.51539827E-01,
             0.91904439E-02, 0.20621480E-01, 0.26908312E-02, 0.14711594E-01,
             0.28122162E-01,-0.72037257E-01, 0.56664880E-01,-0.14147057E-01,
             0.12054422E-01,-0.19986168E-01,-0.53446036E-03, 0.11060343E-01,
             0.38177755E-01,-0.17866375E-01, 0.11385725E-02, 0.40903282E-02,
             0.65648678E-03,-0.58420231E-02, 0.15422775E-01, 0.10356030E-02,
             0.34673312E-02, 0.14013611E-01, 0.11733919E-02,-0.11494552E-02,
            -0.10700230E-01,-0.65248855E-02, 0.34397745E-02,-0.56657111E-04,
            -0.80033978E-02,-0.13102611E-02, 0.14842374E-02,-0.25863241E-03,
            -0.25358321E-02,-0.12770157E-02,-0.30787089E-02,-0.41352821E-03,
            -0.15801092E-02, 0.29393956E-02,-0.61334101E-02, 0.14607376E-01,
            -0.59705833E-02,-0.14852655E-01,-0.17786623E-02, 0.14703891E-02,
            -0.15195629E-01, 0.86391373E-02, 0.37339111E-02,-0.40840339E-02,
             0.82329446E-02, 0.64510126E-02, 0.44786427E-02,-0.18883625E-02,
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
     
        float v_y00_l5p65                       =avdat
            +coeff[  0]                    
            +coeff[  1]        *x31        
            +coeff[  2]            *x41    
            +coeff[  3]*x11                
            +coeff[  4]        *x32        
            +coeff[  5]        *x31*x41    
            +coeff[  6]*x11    *x31        
            +coeff[  7]    *x21*x31        
        ;
        v_y00_l5p65                       =v_y00_l5p65                       
            +coeff[  8]*x11        *x41    
            +coeff[  9]    *x21    *x41    
            +coeff[ 10]*x11*x21            
            +coeff[ 11]    *x22            
            +coeff[ 12]        *x33        
            +coeff[ 13]        *x32*x41    
            +coeff[ 14]        *x31*x42    
            +coeff[ 15]            *x43    
            +coeff[ 16]*x11    *x32        
        ;
        v_y00_l5p65                       =v_y00_l5p65                       
            +coeff[ 17]*x11    *x31*x41    
            +coeff[ 18]*x12    *x31        
            +coeff[ 19]*x11*x21*x31        
            +coeff[ 20]    *x22*x31        
            +coeff[ 21]*x11*x21    *x41    
            +coeff[ 22]*x12*x21            
            +coeff[ 23]    *x23            
            +coeff[ 24]        *x32*x42    
            +coeff[ 25]        *x31*x43    
        ;
        v_y00_l5p65                       =v_y00_l5p65                       
            +coeff[ 26]*x11    *x33        
            +coeff[ 27]*x11    *x31*x42    
            +coeff[ 28]*x12        *x42    
            +coeff[ 29]*x11*x22    *x41    
            +coeff[ 30]        *x31    *x51
            +coeff[ 31]*x11            *x51
            +coeff[ 32]    *x21        *x51
            +coeff[ 33]    *x23    *x42    
            +coeff[ 34]        *x32*x44    
        ;
        v_y00_l5p65                       =v_y00_l5p65                       
            +coeff[ 35]        *x32    *x51
            +coeff[ 36]    *x21*x31    *x51
            +coeff[ 37]*x13*x22            
            +coeff[ 38]*x12*x23            
            +coeff[ 39]*x12            *x51
            +coeff[ 40]    *x22        *x51
            +coeff[ 41]    *x21*x31*x41*x51
            +coeff[ 42]*x13*x22*x31        
            +coeff[ 43]*x12*x23*x31        
        ;
        v_y00_l5p65                       =v_y00_l5p65                       
            +coeff[ 44]*x11*x21*x31    *x51
            +coeff[ 45]    *x23        *x51
            +coeff[ 46]    *x21            
            +coeff[ 47]            *x42    
            +coeff[ 48]*x12        *x41    
            +coeff[ 49]    *x22    *x41    
            +coeff[ 50]*x11*x22            
            +coeff[ 51]                *x51
            +coeff[ 52]*x11    *x32*x41    
        ;
        v_y00_l5p65                       =v_y00_l5p65                       
            +coeff[ 53]    *x22*x32        
            +coeff[ 54]*x11*x21*x31*x41    
            +coeff[ 55]*x11*x22*x31        
            +coeff[ 56]    *x23*x31        
            +coeff[ 57]*x12*x21    *x41    
            +coeff[ 58]        *x34*x41    
            +coeff[ 59]            *x41*x51
            ;
     
        return v_y00_l5p65                       ;
    }
}

