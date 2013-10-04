#include "Bwd_l5p785_484816R15.h"

namespace S484816R15
{
float txfit_l5p785                            (float *x,int m){
    int ncoeff=  2;
    float avdat=  0.2701266E-02;
    float xmin[10]={
        -0.67363E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.60496E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[  3]={
        -0.92148632E-02, 0.11361258E+00,
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

    float v_txfit_l5p785                            =avdat
        +coeff[  0]    
        +coeff[  1]*x11
        ;

    return v_txfit_l5p785                            ;
}
float delta_l5p785                            (float *x,int m){
    int ncoeff= 27;
    float avdat=  0.2606124E-02;
    float xmin[10]={
        -0.67363E+00,-0.26141E-01,-0.56814E-01,-0.49659E-01,-0.14988E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61352E+00, 0.20733E-01, 0.52640E-01, 0.49812E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 28]={
        -0.50260937E-02, 0.47071032E-01, 0.42860452E-02, 0.38573628E-02,
         0.26577062E-02,-0.55121994E-02, 0.48650387E-02, 0.48045421E-03,
        -0.13955381E-02, 0.21870015E-02, 0.13800788E-03,-0.25796844E-02,
        -0.53548762E-02, 0.14840937E-02,-0.51839481E-03,-0.21450797E-02,
         0.25947245E-02, 0.89572795E-03,-0.14464570E-02, 0.33941524E-03,
        -0.10655499E-02, 0.12196632E-02, 0.42751915E-03,-0.92004589E-03,
         0.59020044E-02,-0.31909649E-02, 0.45118894E-03,
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

    float v_delta_l5p785                            =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]*x12                
        +coeff[  3]*x11*x21            
        +coeff[  4]                *x51
        +coeff[  5]    *x21*x31        
        +coeff[  6]    *x21    *x41    
        +coeff[  7]    *x21            
    ;
    v_delta_l5p785                            =v_delta_l5p785                            
        +coeff[  8]    *x22            
        +coeff[  9]        *x32        
        +coeff[ 10]*x11        *x41    
        +coeff[ 11]        *x31*x41    
        +coeff[ 12]    *x21*x32        
        +coeff[ 13]        *x31    *x51
        +coeff[ 14]*x11    *x31        
        +coeff[ 15]    *x22*x31        
        +coeff[ 16]    *x22    *x41    
    ;
    v_delta_l5p785                            =v_delta_l5p785                            
        +coeff[ 17]    *x24            
        +coeff[ 18]*x11*x22    *x41    
        +coeff[ 19]    *x21        *x51
        +coeff[ 20]            *x41*x51
        +coeff[ 21]            *x42    
        +coeff[ 22]*x13                
        +coeff[ 23]*x11*x21    *x41    
        +coeff[ 24]    *x21*x31*x41    
        +coeff[ 25]    *x21    *x42    
    ;
    v_delta_l5p785                            =v_delta_l5p785                            
        +coeff[ 26]        *x32    *x51
        ;

    return v_delta_l5p785                            ;
}
float theta_l5p785                            (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.1636986E-02;
    float xmin[10]={
        -0.67363E+00,-0.26141E-01,-0.56814E-01,-0.49659E-01,-0.14988E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61352E+00, 0.20733E-01, 0.52640E-01, 0.49812E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 51]={
         0.78286892E-02,-0.73323213E-02,-0.55560276E-01, 0.71071897E-03,
        -0.49856288E-03,-0.22636335E-02, 0.49943356E-02,-0.36840879E-02,
        -0.36243882E-01, 0.36547631E-01, 0.33306149E-02, 0.71391817E-02,
        -0.20479560E-02,-0.13160788E-01,-0.11733328E-01, 0.46150889E-02,
         0.71698809E-02,-0.17508256E-02, 0.11283785E-01,-0.99892002E-02,
         0.10634024E-01,-0.47970810E-02,-0.36095891E-01, 0.54026525E-01,
        -0.15205430E-01,-0.53038220E-02, 0.25145870E-02, 0.19304002E-01,
        -0.69567855E-02,-0.91247084E-02,-0.75899321E-03,-0.91843950E-02,
        -0.19832468E-02, 0.19475123E-02, 0.36826848E-04, 0.40133805E-02,
         0.16018380E-02,-0.20945184E-02, 0.58692932E-03,-0.17176535E-01,
         0.96876826E-02,-0.78398231E-02, 0.85458858E-02, 0.11290393E-01,
        -0.17913632E-02,-0.14914316E-02,-0.21786711E-02, 0.51662320E-03,
         0.14985079E-02, 0.11612388E-02,
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

    float v_theta_l5p785                            =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]            *x41    
        +coeff[  5]*x12                
        +coeff[  6]*x11*x21            
        +coeff[  7]    *x22            
    ;
    v_theta_l5p785                            =v_theta_l5p785                            
        +coeff[  8]    *x21*x31        
        +coeff[  9]    *x21    *x41    
        +coeff[ 10]*x11*x22            
        +coeff[ 11]    *x23            
        +coeff[ 12]*x11*x21*x31        
        +coeff[ 13]    *x22*x31        
        +coeff[ 14]*x11*x21    *x41    
        +coeff[ 15]    *x22    *x41    
        +coeff[ 16]                *x51
    ;
    v_theta_l5p785                            =v_theta_l5p785                            
        +coeff[ 17]*x11            *x51
        +coeff[ 18]        *x31    *x51
        +coeff[ 19]            *x41*x51
        +coeff[ 20]    *x23*x32        
        +coeff[ 21]*x11    *x31        
        +coeff[ 22]    *x21*x32        
        +coeff[ 23]    *x21*x31*x41    
        +coeff[ 24]*x11*x21*x31*x41    
        +coeff[ 25]    *x22        *x51
    ;
    v_theta_l5p785                            =v_theta_l5p785                            
        +coeff[ 26]        *x32    *x51
        +coeff[ 27]    *x23*x31        
        +coeff[ 28]    *x22*x32        
        +coeff[ 29]    *x21*x33        
        +coeff[ 30]*x11*x21        *x51
        +coeff[ 31]    *x22*x31    *x51
        +coeff[ 32]    *x21*x32    *x51
        +coeff[ 33]            *x43*x51
        +coeff[ 34]        *x32        
    ;
    v_theta_l5p785                            =v_theta_l5p785                            
        +coeff[ 35]*x11        *x41    
        +coeff[ 36]        *x31*x41    
        +coeff[ 37]*x11    *x32        
        +coeff[ 38]        *x33        
        +coeff[ 39]    *x21    *x42    
        +coeff[ 40]*x11*x22*x31        
        +coeff[ 41]*x11*x22    *x41    
        +coeff[ 42]    *x21*x32*x41    
        +coeff[ 43]*x11*x21    *x42    
    ;
    v_theta_l5p785                            =v_theta_l5p785                            
        +coeff[ 44]*x11    *x31*x42    
        +coeff[ 45]*x12*x23            
        +coeff[ 46]*x12*x21*x32        
        +coeff[ 47]*x12            *x51
        +coeff[ 48]    *x21*x31    *x51
        +coeff[ 49]*x11        *x41*x51
        ;

    return v_theta_l5p785                            ;
}
float phi_l5p785                              (float *x,int m){
    int ncoeff= 70;
    float avdat= -0.2417408E-02;
    float xmin[10]={
        -0.67363E+00,-0.26141E-01,-0.56814E-01,-0.49659E-01,-0.14988E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61352E+00, 0.20733E-01, 0.52640E-01, 0.49812E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 71]={
        -0.25688335E-02,-0.26419979E-01,-0.58889086E-02, 0.55891583E-02,
         0.68461453E-03, 0.78945234E-02,-0.97377095E-02, 0.91557875E-02,
        -0.15715191E-01,-0.11717387E-01, 0.20853935E-02, 0.34476558E-02,
         0.36085006E-02,-0.49676007E-03, 0.95741041E-02,-0.44235639E-01,
         0.61411474E-01,-0.27410580E-01, 0.54491585E-03, 0.50428114E-02,
         0.93651644E-04, 0.70927492E-02,-0.46711056E-02,-0.79375779E-03,
         0.14638882E-01, 0.11253374E-02,-0.27083503E-02,-0.12627065E-02,
        -0.32448459E-01, 0.29071146E-02, 0.31936463E-01, 0.32994314E-02,
        -0.37127405E-01,-0.20689513E-02, 0.11281308E-02, 0.16790515E-01,
        -0.18115127E-01,-0.18838823E-02, 0.20999843E-02, 0.33715230E-02,
        -0.74903946E-02, 0.62378403E-02,-0.67234840E-02, 0.46112541E-01,
        -0.53186915E-02, 0.88099577E-03,-0.95586403E-03, 0.25678361E-02,
        -0.22919159E-02, 0.14288709E-01,-0.12942438E-01, 0.11896965E-01,
        -0.11169466E-02,-0.13433628E-01,-0.28921566E-02, 0.21621329E-02,
        -0.43250751E-02,-0.17658126E-02,-0.69858000E-03, 0.13251979E-02,
        -0.70413374E-02, 0.15481170E-02,-0.51768995E-02,-0.68354760E-02,
        -0.73136333E-02, 0.11555344E-02, 0.71442272E-02,-0.12916987E-02,
        -0.55783913E-02, 0.39242753E-02,
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

//                 function

    float v_phi_l5p785                              =avdat
        +coeff[  0]                    
        +coeff[  1]        *x31        
        +coeff[  2]            *x41    
        +coeff[  3]*x11                
        +coeff[  4]    *x21            
        +coeff[  5]        *x32        
        +coeff[  6]        *x31*x41    
        +coeff[  7]    *x21*x31        
    ;
    v_phi_l5p785                              =v_phi_l5p785                              
        +coeff[  8]*x11        *x41    
        +coeff[  9]    *x22            
        +coeff[ 10]*x11    *x32        
        +coeff[ 11]    *x21*x32        
        +coeff[ 12]*x11        *x42    
        +coeff[ 13]*x12*x21            
        +coeff[ 14]*x11    *x33        
        +coeff[ 15]*x11    *x32*x41    
        +coeff[ 16]*x11    *x31*x42    
    ;
    v_phi_l5p785                              =v_phi_l5p785                              
        +coeff[ 17]*x11        *x43    
        +coeff[ 18]*x12    *x31*x41    
        +coeff[ 19]    *x24            
        +coeff[ 20]                *x51
        +coeff[ 21]    *x21        *x51
        +coeff[ 22]    *x23        *x51
        +coeff[ 23]            *x42    
        +coeff[ 24]*x11    *x31        
        +coeff[ 25]    *x21    *x41    
    ;
    v_phi_l5p785                              =v_phi_l5p785                              
        +coeff[ 26]*x11*x21            
        +coeff[ 27]*x12    *x31        
        +coeff[ 28]    *x22*x31        
        +coeff[ 29]*x11*x21    *x41    
        +coeff[ 30]    *x22    *x41    
        +coeff[ 31]*x11*x22            
        +coeff[ 32]    *x22*x32        
        +coeff[ 33]*x12        *x42    
        +coeff[ 34]        *x31    *x51
    ;
    v_phi_l5p785                              =v_phi_l5p785                              
        +coeff[ 35]    *x21*x31    *x51
        +coeff[ 36]    *x21    *x41*x51
        +coeff[ 37]*x11*x21        *x51
        +coeff[ 38]*x11    *x32    *x51
        +coeff[ 39]        *x32*x41    
        +coeff[ 40]*x11*x21*x31        
        +coeff[ 41]*x12        *x41    
        +coeff[ 42]*x11*x21*x32        
        +coeff[ 43]    *x22*x31*x41    
    ;
    v_phi_l5p785                              =v_phi_l5p785                              
        +coeff[ 44]    *x23*x31        
        +coeff[ 45]*x11            *x51
        +coeff[ 46]        *x31*x41*x51
        +coeff[ 47]*x11    *x31    *x51
        +coeff[ 48]        *x31*x42*x51
        +coeff[ 49]    *x21*x32    *x51
        +coeff[ 50]    *x21*x31*x41*x51
        +coeff[ 51]    *x21*x32*x41*x51
        +coeff[ 52]    *x24        *x51
    ;
    v_phi_l5p785                              =v_phi_l5p785                              
        +coeff[ 53]    *x23*x33    *x51
        +coeff[ 54]        *x33        
        +coeff[ 55]    *x21*x31*x41    
        +coeff[ 56]    *x22    *x42    
        +coeff[ 57]*x13    *x31        
        +coeff[ 58]*x12*x21*x31        
        +coeff[ 59]*x13        *x41    
        +coeff[ 60]*x11*x22    *x41    
        +coeff[ 61]*x11*x23            
    ;
    v_phi_l5p785                              =v_phi_l5p785                              
        +coeff[ 62]    *x22*x33        
        +coeff[ 63]    *x23*x32        
        +coeff[ 64]*x11*x22*x31*x41    
        +coeff[ 65]    *x24*x31        
        +coeff[ 66]    *x24    *x41    
        +coeff[ 67]            *x41*x51
        +coeff[ 68]    *x22*x33*x41    
        +coeff[ 69]*x11*x21*x32*x42    
        ;

    return v_phi_l5p785                              ;
}
float y00_l5p785                              (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.7281741E-03;
    float xmin[10]={
        -0.67363E+00,-0.26141E-01,-0.56814E-01,-0.49659E-01,-0.14988E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.61352E+00, 0.20733E-01, 0.52640E-01, 0.49812E-01, 0.15000E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float coeff[ 61]={
         0.71167271E-02,-0.48264183E-01, 0.95211320E-01,-0.76308232E-02,
         0.58573484E-03,-0.75104027E-02, 0.10845303E-02, 0.87214205E-02,
        -0.46774078E-01,-0.49886972E-01, 0.15576288E-01, 0.20699032E-01,
        -0.48935931E-03, 0.31883956E-02, 0.13152758E-01,-0.71322615E-02,
         0.10025625E-02, 0.15171115E-01, 0.59663199E-01,-0.13141427E-01,
        -0.19020895E-01,-0.42797089E-01,-0.57081538E-02,-0.20740191E-03,
         0.10866559E-01, 0.40659169E-02,-0.40999260E-02,-0.18317895E-02,
        -0.20006674E-02,-0.91940826E-02,-0.22734294E-01, 0.20664511E-01,
        -0.12508994E-02, 0.12793585E-02,-0.36180010E-02, 0.10008450E-01,
         0.17467994E-02,-0.16845954E-02, 0.32980398E-02, 0.17066229E-01,
        -0.31708520E-01, 0.16380522E-01,-0.46866871E-02, 0.23733963E-03,
         0.36621136E-02, 0.47132712E-01,-0.37364606E-01, 0.95106708E-02,
         0.14900285E-01, 0.12236992E-02, 0.30204942E-02, 0.60281926E-02,
         0.10628526E-02,-0.27393554E-02, 0.74515573E-03,-0.25754811E-01,
         0.38055029E-01,-0.14624398E-01,-0.95096435E-02,-0.75640394E-02,
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
    float x51 = x5;

//                 function

    float v_y00_l5p785                              =avdat
        +coeff[  0]                    
        +coeff[  1]        *x31        
        +coeff[  2]            *x41    
        +coeff[  3]*x11                
        +coeff[  4]    *x21            
        +coeff[  5]        *x32        
        +coeff[  6]        *x31*x41    
        +coeff[  7]            *x42    
    ;
    v_y00_l5p785                              =v_y00_l5p785                              
        +coeff[  8]*x11    *x31        
        +coeff[  9]    *x21*x31        
        +coeff[ 10]*x11        *x41    
        +coeff[ 11]    *x21    *x41    
        +coeff[ 12]*x12                
        +coeff[ 13]*x11*x21            
        +coeff[ 14]    *x22            
        +coeff[ 15]    *x21*x31*x41    
        +coeff[ 16]*x12    *x31        
    ;
    v_y00_l5p785                              =v_y00_l5p785                              
        +coeff[ 17]*x11*x21*x31        
        +coeff[ 18]    *x22*x31        
        +coeff[ 19]*x12        *x41    
        +coeff[ 20]*x11*x21    *x41    
        +coeff[ 21]    *x22    *x41    
        +coeff[ 22]*x11*x22            
        +coeff[ 23]                *x51
        +coeff[ 24]*x11*x21*x32        
        +coeff[ 25]*x13    *x31        
    ;
    v_y00_l5p785                              =v_y00_l5p785                              
        +coeff[ 26]    *x23    *x41    
        +coeff[ 27]        *x34*x41    
        +coeff[ 28]        *x31    *x51
        +coeff[ 29]    *x21        *x51
        +coeff[ 30]    *x21*x31    *x51
        +coeff[ 31]    *x21    *x41*x51
        +coeff[ 32]    *x22*x34        
        +coeff[ 33]*x11*x21        *x51
        +coeff[ 34]*x11    *x32    *x51
    ;
    v_y00_l5p785                              =v_y00_l5p785                              
        +coeff[ 35]*x12*x23*x31        
        +coeff[ 36]*x11*x21*x31    *x51
        +coeff[ 37]*x13            *x51
        +coeff[ 38]    *x23        *x51
        +coeff[ 39]        *x33        
        +coeff[ 40]        *x32*x41    
        +coeff[ 41]        *x31*x42    
        +coeff[ 42]*x11    *x31*x41    
        +coeff[ 43]        *x34        
    ;
    v_y00_l5p785                              =v_y00_l5p785                              
        +coeff[ 44]*x11    *x33        
        +coeff[ 45]    *x22*x32        
        +coeff[ 46]    *x22*x31*x41    
        +coeff[ 47]    *x23*x31        
        +coeff[ 48]*x11*x22    *x41    
        +coeff[ 49]            *x41*x51
        +coeff[ 50]*x12*x22            
        +coeff[ 51]    *x22*x33        
        +coeff[ 52]        *x31*x41*x51
    ;
    v_y00_l5p785                              =v_y00_l5p785                              
        +coeff[ 53]*x11    *x31    *x51
        +coeff[ 54]        *x33    *x51
        +coeff[ 55]    *x21*x32    *x51
        +coeff[ 56]    *x21*x31*x41*x51
        +coeff[ 57]    *x21    *x42*x51
        +coeff[ 58]*x13*x22*x31        
        +coeff[ 59]*x11*x21    *x41*x51
        ;

    return v_y00_l5p785                              ;
}
}
