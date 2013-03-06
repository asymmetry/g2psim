#include "Bwd_l5p77_484816R00.h"

namespace S484816R00
{
    float txfit_l5p77          (float *x,int m){
        int ncoeff=  2;
        float avdat= -0.5051652E-02;
        float xmin[10]={
            -0.69016E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.59782E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[  3]={
            -0.37142362E-02, 0.11289136E+00,
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
     
        float v_txfit_l5p77          =avdat
            +coeff[  0]    
            +coeff[  1]*x11
            ;
     
        return v_txfit_l5p77          ;
    }
    float delta_l5p77          (float *x,int m){
        int ncoeff= 18;
        float avdat=  0.5034296E-03;
        float xmin[10]={
            -0.70542E+00,-0.22107E-01,-0.46979E-01,-0.30069E-01, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.61246E+00, 0.18344E-01, 0.49705E-01, 0.44151E-01, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 19]={
            -0.45872242E-02, 0.48122853E-01, 0.41402006E-02, 0.39591291E-02,
            -0.29271997E-02,-0.21437474E-02, 0.56736398E-03,-0.15408965E-02,
             0.15305240E-02, 0.54070377E-03, 0.47061394E-03, 0.51013345E-03,
             0.21917562E-03,-0.44351269E-03,-0.89930749E-03, 0.63209207E-03,
            -0.49693318E-03, 0.58011909E-03,
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
        float x21 = x2;
        float x22 = x21*x2;
        float x23 = x22*x2;
        float x24 = x23*x2;
        float x31 = x3;
        float x32 = x31*x3;
        float x41 = x4;
     
    //                 function
     
        float v_delta_l5p77          =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]*x12                
            +coeff[  3]*x11*x21            
            +coeff[  4]    *x21*x31        
            +coeff[  5]    *x21*x32        
            +coeff[  6]    *x24            
            +coeff[  7]    *x22            
        ;
        v_delta_l5p77          =v_delta_l5p77          
            +coeff[  8]        *x32        
            +coeff[  9]    *x21    *x41    
            +coeff[ 10]*x13                
            +coeff[ 11]*x11*x22            
            +coeff[ 12]        *x31        
            +coeff[ 13]            *x41    
            +coeff[ 14]        *x31*x41    
            +coeff[ 15]*x11*x21*x31        
            +coeff[ 16]    *x22*x31        
        ;
        v_delta_l5p77          =v_delta_l5p77          
            +coeff[ 17]    *x23*x31        
            ;
     
        return v_delta_l5p77          ;
    }
    float theta_l5p77          (float *x,int m){
        int ncoeff= 40;
        float avdat= -0.3971255E-03;
        float xmin[10]={
            -0.70542E+00,-0.22107E-01,-0.46979E-01,-0.30069E-01, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.61246E+00, 0.18344E-01, 0.49705E-01, 0.44151E-01, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 41]={
             0.49985838E-02,-0.39173942E-02,-0.49279194E-01,-0.15320665E-02,
            -0.22362464E-02, 0.67868209E-02,-0.26056264E-02, 0.34901053E-02,
            -0.12814080E-01, 0.17116034E-01,-0.14988331E-01, 0.19374939E-03,
             0.11031644E-02, 0.20544867E-02,-0.48480986E-03, 0.12846899E-02,
            -0.62206476E-02, 0.11973709E-01,-0.18150895E-02, 0.21469037E-03,
             0.88219572E-03,-0.36735804E-04,-0.19103611E-02,-0.10156128E-01,
            -0.32532492E-03,-0.56575858E-02,-0.17861429E-02,-0.11881483E-01,
             0.17284817E-02, 0.10124472E-01, 0.63683727E-03, 0.56275859E-03,
            -0.19294531E-02, 0.79208268E-02,-0.11841056E-01, 0.52631791E-02,
            -0.17624695E-02, 0.22386918E-02, 0.25065732E-02,-0.61466527E-03,
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
        float x21 = x2;
        float x22 = x21*x2;
        float x23 = x22*x2;
        float x31 = x3;
        float x32 = x31*x3;
        float x33 = x32*x3;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
     
    //                 function
     
        float v_theta_l5p77          =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]            *x41    
            +coeff[  4]*x12                
            +coeff[  5]*x11*x21            
            +coeff[  6]    *x22            
            +coeff[  7]*x11    *x31        
        ;
        v_theta_l5p77          =v_theta_l5p77          
            +coeff[  8]    *x21*x31        
            +coeff[  9]*x11*x21*x31        
            +coeff[ 10]    *x21*x32        
            +coeff[ 11]*x11*x22            
            +coeff[ 12]        *x31        
            +coeff[ 13]        *x32        
            +coeff[ 14]*x12*x21            
            +coeff[ 15]            *x43    
            +coeff[ 16]*x12*x21*x31        
        ;
        v_theta_l5p77          =v_theta_l5p77          
            +coeff[ 17]*x11*x21*x32        
            +coeff[ 18]*x11    *x33        
            +coeff[ 19]*x13                
            +coeff[ 20]    *x21    *x42    
            +coeff[ 21]*x12*x22            
            +coeff[ 22]    *x22*x32        
            +coeff[ 23]    *x21*x33        
            +coeff[ 24]        *x31*x43    
            +coeff[ 25]*x11        *x41    
        ;
        v_theta_l5p77          =v_theta_l5p77          
            +coeff[ 26]        *x31*x41    
            +coeff[ 27]*x11*x21    *x41    
            +coeff[ 28]    *x22    *x41    
            +coeff[ 29]    *x21*x31*x41    
            +coeff[ 30]*x13*x21            
            +coeff[ 31]*x11*x23            
            +coeff[ 32]*x11*x22*x31        
            +coeff[ 33]*x12*x21    *x41    
            +coeff[ 34]*x11*x21*x31*x41    
        ;
        v_theta_l5p77          =v_theta_l5p77          
            +coeff[ 35]    *x21*x32*x41    
            +coeff[ 36]    *x21*x33*x41    
            +coeff[ 37]*x12*x22*x32        
            +coeff[ 38]*x12    *x31        
            +coeff[ 39]*x12        *x41    
            ;
     
        return v_theta_l5p77          ;
    }
    float phi_l5p77            (float *x,int m){
        int ncoeff= 70;
        float avdat=  0.1560849E-02;
        float xmin[10]={
            -0.70542E+00,-0.22107E-01,-0.46979E-01,-0.30069E-01, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.61246E+00, 0.18344E-01, 0.49705E-01, 0.44151E-01, 0.00000E+00,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 71]={
            -0.10325093E-02,-0.28951598E-01, 0.72670332E-02, 0.50966284E-03,
             0.18644988E-02,-0.61076325E-02, 0.24393076E-01,-0.15637545E-01,
             0.11524010E-01,-0.47155377E-02,-0.22305641E-01, 0.10126584E-01,
            -0.69011748E-02,-0.37376321E-03, 0.11513538E-01,-0.10189350E-02,
            -0.79365373E-02,-0.51137088E-02,-0.11956488E-02, 0.10336140E-02,
            -0.54444429E-02,-0.10315485E-01,-0.97639477E-02, 0.16939973E-02,
            -0.20107662E-02, 0.20463143E-02, 0.58668097E-02, 0.24301040E-02,
             0.25572518E-02,-0.18079837E-02, 0.15153234E-01,-0.76862625E-02,
            -0.48405625E-03,-0.56420406E-03, 0.37372197E-03,-0.15778200E-02,
             0.53532491E-02, 0.16870649E-02, 0.23375158E-03, 0.64078794E-03,
             0.22368843E-02,-0.50313068E-02, 0.58924081E-02,-0.78466712E-02,
             0.33821433E-02,-0.28096014E-03,-0.12984859E-02,-0.88210485E-03,
            -0.25649914E-02, 0.22205788E-03, 0.24366602E-02, 0.82398532E-02,
             0.62744836E-02, 0.39985153E-03, 0.64522732E-03,-0.25688516E-03,
            -0.32964512E-02,-0.47652270E-02, 0.39307624E-02,-0.42505935E-02,
             0.32469090E-02,-0.65668472E-02,-0.97512389E-02, 0.24716503E-02,
             0.32808071E-02, 0.19086858E-03, 0.23924487E-03, 0.17258001E-02,
            -0.11351557E-02, 0.85677608E-03,
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
     
    //                 function
     
        float v_phi_l5p77            =avdat
            +coeff[  0]                    
            +coeff[  1]        *x31        
            +coeff[  2]            *x41    
            +coeff[  3]*x11                
            +coeff[  4]    *x21            
            +coeff[  5]        *x32        
            +coeff[  6]        *x31*x41    
            +coeff[  7]            *x42    
        ;
        v_phi_l5p77            =v_phi_l5p77            
            +coeff[  8]*x11    *x31        
            +coeff[  9]    *x21*x31        
            +coeff[ 10]*x11        *x41    
            +coeff[ 11]    *x21    *x41    
            +coeff[ 12]    *x22            
            +coeff[ 13]*x11*x22            
            +coeff[ 14]    *x21*x31*x41    
            +coeff[ 15]    *x21    *x42    
            +coeff[ 16]*x11*x21*x31*x41    
        ;
        v_phi_l5p77            =v_phi_l5p77            
            +coeff[ 17]*x11        *x42    
            +coeff[ 18]*x12*x21            
            +coeff[ 19]    *x23            
            +coeff[ 20]        *x33        
            +coeff[ 21]*x11    *x32        
            +coeff[ 22]    *x21*x32        
            +coeff[ 23]*x13    *x31        
            +coeff[ 24]    *x23*x32        
            +coeff[ 25]*x12                
        ;
        v_phi_l5p77            =v_phi_l5p77            
            +coeff[ 26]*x11*x21    *x42    
            +coeff[ 27]*x11*x21*x33        
            +coeff[ 28]*x13*x22*x32        
            +coeff[ 29]*x11*x21            
            +coeff[ 30]*x11    *x31*x41    
            +coeff[ 31]    *x22*x31        
            +coeff[ 32]*x12        *x41    
            +coeff[ 33]*x12    *x32        
            +coeff[ 34]*x12    *x31*x41    
        ;
        v_phi_l5p77            =v_phi_l5p77            
            +coeff[ 35]*x12        *x42    
            +coeff[ 36]    *x23*x31        
            +coeff[ 37]*x13        *x41    
            +coeff[ 38]        *x34*x41    
            +coeff[ 39]*x12*x22            
            +coeff[ 40]*x11*x21*x32*x41    
            +coeff[ 41]*x12    *x31*x42    
            +coeff[ 42]*x12        *x43    
            +coeff[ 43]*x11*x22*x31*x41    
        ;
        v_phi_l5p77            =v_phi_l5p77            
            +coeff[ 44]    *x23*x31*x41    
            +coeff[ 45]*x13*x21*x31        
            +coeff[ 46]*x13*x21    *x41    
            +coeff[ 47]*x12    *x32*x42    
            +coeff[ 48]*x13*x21*x31*x41    
            +coeff[ 49]*x13*x21*x33        
            +coeff[ 50]*x12*x22*x33        
            +coeff[ 51]        *x32*x41    
            +coeff[ 52]*x12    *x31        
        ;
        v_phi_l5p77            =v_phi_l5p77            
            +coeff[ 53]*x11*x21    *x41    
            +coeff[ 54]*x13                
            +coeff[ 55]        *x34        
            +coeff[ 56]*x11    *x33        
            +coeff[ 57]    *x21*x32*x41    
            +coeff[ 58]    *x21    *x43    
            +coeff[ 59]    *x22*x32        
            +coeff[ 60]    *x22*x31*x41    
            +coeff[ 61]*x11*x22    *x41    
        ;
        v_phi_l5p77            =v_phi_l5p77            
            +coeff[ 62]*x11*x21*x31*x42    
            +coeff[ 63]*x13    *x31*x41    
            +coeff[ 64]*x12*x21*x31*x41    
            +coeff[ 65]*x12*x21    *x42    
            +coeff[ 66]        *x33*x43    
            +coeff[ 67]*x12*x22*x31        
            +coeff[ 68]*x13    *x33        
            +coeff[ 69]*x13*x22    *x41    
            ;
     
        return v_phi_l5p77            ;
    }
}
