#include "Fwd_l5p65_400016.h"

namespace S400016
{    
    float x_l5p65_ep10_q1en               (float *x,int m){
        int ncoeff= 45;
        float avdat= -0.2266574E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46818E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 46]={
             0.44502190E-03, 0.97995184E-01, 0.65494905E-03, 0.20548312E-02,
             0.82005793E-02, 0.41945550E-04, 0.15429505E-04, 0.17023881E-02,
            -0.69776055E-04,-0.21771541E-04, 0.49322556E-04, 0.13076128E-01,
             0.26645914E-02, 0.35461644E-03, 0.26132620E-02,-0.32753560E-02,
            -0.43642565E-04,-0.26827745E-03, 0.14571303E-03, 0.31371901E-03,
            -0.25114040E-02,-0.11521666E-03, 0.66147193E-04, 0.41237104E-03,
            -0.63457672E-03,-0.29394106E-03, 0.36157604E-03,-0.20555244E-02,
            -0.51400275E-04,-0.14814588E-02, 0.10364782E-02, 0.58169189E-04,
            -0.62809681E-03, 0.51764800E-04, 0.60550869E-04,-0.67912792E-04,
             0.13071590E-03,-0.45739266E-03,-0.93200986E-04,-0.34236570E-03,
            -0.50224410E-03, 0.13333359E-04, 0.84964861E-03,-0.45774083E-03,
             0.45726742E-03,
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
     
        float v_x_l5p65_ep10_q1en               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]*x11    *x31        
            +coeff[  3]*x11        *x41    
            +coeff[  4]    *x21    *x41    
            +coeff[  5]*x13                
            +coeff[  6]*x11            *x52
            +coeff[  7]*x11*x22            
        ;
        v_x_l5p65_ep10_q1en               =v_x_l5p65_ep10_q1en               
            +coeff[  8]*x11    *x32        
            +coeff[  9]*x13            *x52
            +coeff[ 10]*x13*x22            
            +coeff[ 11]*x11                
            +coeff[ 12]    *x21*x31        
            +coeff[ 13]            *x41    
            +coeff[ 14]    *x23            
            +coeff[ 15]    *x21    *x42    
            +coeff[ 16]*x12*x21*x31*x41    
        ;
        v_x_l5p65_ep10_q1en               =v_x_l5p65_ep10_q1en               
            +coeff[ 17]    *x23*x31*x41    
            +coeff[ 18]        *x31        
            +coeff[ 19]    *x21        *x51
            +coeff[ 20]    *x23    *x41    
            +coeff[ 21]*x13*x22    *x41    
            +coeff[ 22]                *x51
            +coeff[ 23]*x12*x21            
            +coeff[ 24]*x11        *x42    
            +coeff[ 25]    *x21    *x41*x51
        ;
        v_x_l5p65_ep10_q1en               =v_x_l5p65_ep10_q1en               
            +coeff[ 26]    *x22    *x41    
            +coeff[ 27]    *x21*x31*x41    
            +coeff[ 28]*x13        *x41    
            +coeff[ 29]*x11*x22    *x41    
            +coeff[ 30]    *x21    *x43    
            +coeff[ 31]*x13    *x31*x41    
            +coeff[ 32]    *x23*x32        
            +coeff[ 33]    *x21*x33*x42    
            +coeff[ 34]*x11            *x51
        ;
        v_x_l5p65_ep10_q1en               =v_x_l5p65_ep10_q1en               
            +coeff[ 35]    *x22            
            +coeff[ 36]*x11*x21    *x41    
            +coeff[ 37]*x11    *x31*x41    
            +coeff[ 38]    *x21*x31    *x51
            +coeff[ 39]*x12*x21    *x41    
            +coeff[ 40]    *x23*x31        
            +coeff[ 41]    *x21*x33        
            +coeff[ 42]    *x21*x31*x42    
            +coeff[ 43]*x13*x22*x31        
        ;
        v_x_l5p65_ep10_q1en               =v_x_l5p65_ep10_q1en               
            +coeff[ 44]    *x23*x32*x41    
            ;
     
        return v_x_l5p65_ep10_q1en               ;
    }
    float t_l5p65_ep10_q1en               (float *x,int m){
        int ncoeff= 50;
        float avdat= -0.2100820E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46818E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 51]={
             0.22642150E-05,-0.22277180E-02, 0.37041429E-01, 0.10842867E-01,
            -0.15020823E-03,-0.11351769E-02,-0.21855782E-04, 0.47502090E-03,
             0.35774857E-02,-0.39197654E-02, 0.16354054E-03, 0.78601227E-03,
             0.23702160E-02, 0.36861422E-03, 0.19674087E-02, 0.27738123E-02,
            -0.37810348E-02,-0.95690655E-06, 0.15250422E-03, 0.46709427E-03,
             0.59084652E-03,-0.26346773E-02,-0.79819042E-03,-0.22566512E-02,
             0.16758174E-02, 0.81367114E-04,-0.43265490E-03,-0.50910149E-03,
            -0.40246872E-03,-0.68837445E-03,-0.47684315E-03, 0.14816623E-02,
            -0.12362784E-04, 0.57599482E-04,-0.10422175E-03, 0.74342315E-04,
             0.16954828E-03, 0.25415173E-03,-0.11585683E-03,-0.78151541E-04,
             0.41806433E-03, 0.17037921E-04, 0.31007949E-03,-0.65513443E-04,
             0.65130764E-04,-0.72313662E-04, 0.22085884E-03, 0.29883371E-03,
            -0.12189488E-04, 0.32380114E-04,
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
     
        float v_t_l5p65_ep10_q1en               =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]    *x21    *x41    
            +coeff[  4]*x12*x21*x31        
            +coeff[  5]    *x23*x31        
            +coeff[  6]*x13        *x41    
            +coeff[  7]            *x41    
        ;
        v_t_l5p65_ep10_q1en               =v_t_l5p65_ep10_q1en               
            +coeff[  8]    *x21*x31        
            +coeff[  9]    *x21    *x42    
            +coeff[ 10]    *x23*x31*x41    
            +coeff[ 11]*x11    *x31        
            +coeff[ 12]*x11        *x41    
            +coeff[ 13]    *x21        *x51
            +coeff[ 14]*x11*x22            
            +coeff[ 15]    *x23            
            +coeff[ 16]    *x23    *x41    
        ;
        v_t_l5p65_ep10_q1en               =v_t_l5p65_ep10_q1en               
            +coeff[ 17]*x13*x22    *x41    
            +coeff[ 18]        *x31        
            +coeff[ 19]*x12*x21            
            +coeff[ 20]    *x22    *x41    
            +coeff[ 21]    *x21*x31*x41    
            +coeff[ 22]*x11        *x42    
            +coeff[ 23]*x11*x22    *x41    
            +coeff[ 24]    *x21    *x43    
            +coeff[ 25]                *x51
        ;
        v_t_l5p65_ep10_q1en               =v_t_l5p65_ep10_q1en               
            +coeff[ 26]    *x21*x32        
            +coeff[ 27]*x11    *x31*x41    
            +coeff[ 28]    *x21    *x41*x51
            +coeff[ 29]*x11*x22*x31        
            +coeff[ 30]*x12*x21    *x41    
            +coeff[ 31]    *x21*x31*x42    
            +coeff[ 32]*x13*x21    *x41    
            +coeff[ 33]    *x22            
            +coeff[ 34]            *x42    
        ;
        v_t_l5p65_ep10_q1en               =v_t_l5p65_ep10_q1en               
            +coeff[ 35]*x11            *x51
            +coeff[ 36]    *x22*x31        
            +coeff[ 37]*x11*x21    *x41    
            +coeff[ 38]    *x21*x31    *x51
            +coeff[ 39]*x11        *x41*x51
            +coeff[ 40]    *x21*x32*x41    
            +coeff[ 41]        *x33*x41    
            +coeff[ 42]*x11        *x43    
            +coeff[ 43]        *x31*x41    
        ;
        v_t_l5p65_ep10_q1en               =v_t_l5p65_ep10_q1en               
            +coeff[ 44]*x11*x21*x31        
            +coeff[ 45]*x11    *x32        
            +coeff[ 46]*x11    *x31*x42    
            +coeff[ 47]    *x23    *x42    
            +coeff[ 48]            *x41*x51
            +coeff[ 49]*x13                
            ;
     
        return v_t_l5p65_ep10_q1en               ;
    }
    float y_l5p65_ep10_q1en               (float *x,int m){
        int ncoeff= 97;
        float avdat=  0.1196650E-01;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46818E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 98]={
            -0.13586352E-01, 0.77074923E-01, 0.16790839E-01,-0.40373355E-02,
             0.67090713E-02, 0.30752134E-02,-0.31596117E-02,-0.52938135E-02,
             0.37711742E-03,-0.21336873E-02,-0.24793546E-02, 0.36983201E-03,
            -0.18802813E-02,-0.90744451E-03, 0.33165668E-02, 0.90001347E-04,
            -0.36122114E-03, 0.20641554E-03, 0.22169419E-02, 0.14667261E-02,
             0.10144163E-02,-0.18322584E-03, 0.59825252E-03,-0.33666086E-03,
             0.59476972E-03,-0.26807826E-03,-0.31836724E-03,-0.12315180E-03,
            -0.90916565E-03,-0.82709355E-03,-0.70024398E-04,-0.23088184E-03,
             0.20718470E-03, 0.22754305E-03, 0.16176779E-03, 0.91354348E-04,
            -0.11898506E-03, 0.35185443E-03, 0.17267431E-03,-0.28595957E-03,
             0.14271284E-03, 0.59125439E-04, 0.12360098E-03, 0.71800656E-04,
             0.17904559E-03, 0.16501734E-03,-0.63196855E-03, 0.11025585E-03,
            -0.67014457E-03, 0.65268175E-03, 0.59728208E-03, 0.23279415E-03,
            -0.12390156E-04,-0.63082851E-04,-0.42147949E-04, 0.99046185E-04,
             0.37003651E-04, 0.12882942E-04,-0.15244279E-03,-0.13500484E-04,
            -0.11324848E-03, 0.85697320E-04,-0.31426255E-03, 0.22362949E-04,
            -0.23548759E-03,-0.31912304E-03,-0.46970927E-04,-0.12675759E-03,
             0.39099024E-04, 0.16215407E-04, 0.21921445E-04,-0.49881877E-04,
             0.28646005E-04,-0.47207992E-04,-0.31519503E-04, 0.64905362E-04,
            -0.30170118E-04, 0.17876477E-04, 0.31383810E-04,-0.13223830E-03,
             0.24724210E-03,-0.51296178E-04,-0.47879377E-04, 0.23313625E-03,
             0.90972999E-05,-0.41501591E-03,-0.24762037E-04, 0.92575930E-04,
            -0.26490085E-03,-0.32341795E-03,-0.20771293E-03,-0.15664616E-03,
            -0.87231914E-04,-0.28639552E-04,-0.33892240E-05, 0.98578939E-05,
             0.32703250E-04,
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
     
        float v_y_l5p65_ep10_q1en               =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]        *x31        
            +coeff[  3]                *x51
            +coeff[  4]    *x22            
            +coeff[  5]*x11*x21            
            +coeff[  6]            *x42    
            +coeff[  7]    *x22    *x41    
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[  8]    *x21            
            +coeff[  9]        *x31*x41    
            +coeff[ 10]*x11*x21    *x41    
            +coeff[ 11]*x12                
            +coeff[ 12]    *x22*x31        
            +coeff[ 13]*x11*x21*x31        
            +coeff[ 14]    *x22    *x42    
            +coeff[ 15]*x11                
            +coeff[ 16]        *x32        
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 17]                *x52
            +coeff[ 18]    *x22*x31*x41    
            +coeff[ 19]*x11*x21    *x42    
            +coeff[ 20]*x11*x21*x31*x41    
            +coeff[ 21]            *x41*x51
            +coeff[ 22]            *x43    
            +coeff[ 23]    *x21    *x42    
            +coeff[ 24]        *x31*x42    
            +coeff[ 25]    *x22        *x51
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 26]*x12        *x41    
            +coeff[ 27]*x11*x21        *x51
            +coeff[ 28]    *x24            
            +coeff[ 29]*x11*x23            
            +coeff[ 30]        *x31    *x51
            +coeff[ 31]    *x21*x31*x41    
            +coeff[ 32]        *x32*x41    
            +coeff[ 33]    *x23            
            +coeff[ 34]*x11*x22            
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 35]    *x21    *x43    
            +coeff[ 36]*x12    *x31        
            +coeff[ 37]    *x22*x32        
            +coeff[ 38]*x11*x21*x32        
            +coeff[ 39]*x12*x22            
            +coeff[ 40]    *x21    *x41    
            +coeff[ 41]    *x21*x31        
            +coeff[ 42]            *x42*x51
            +coeff[ 43]        *x31*x41*x51
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 44]    *x22    *x41*x51
            +coeff[ 45]*x12        *x42    
            +coeff[ 46]    *x22    *x43    
            +coeff[ 47]*x12    *x31*x41    
            +coeff[ 48]    *x22*x31*x42    
            +coeff[ 49]    *x24    *x41    
            +coeff[ 50]*x11*x23    *x41    
            +coeff[ 51]*x12*x22    *x41    
            +coeff[ 52]    *x21        *x51
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 53]*x11        *x42    
            +coeff[ 54]    *x21*x32        
            +coeff[ 55]    *x21*x31*x42    
            +coeff[ 56]*x12*x21            
            +coeff[ 57]*x11        *x43    
            +coeff[ 58]    *x23    *x41    
            +coeff[ 59]*x12            *x51
            +coeff[ 60]*x11*x22    *x41    
            +coeff[ 61]*x11*x21    *x41*x51
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 62]*x11*x21    *x43    
            +coeff[ 63]*x11    *x31*x43    
            +coeff[ 64]    *x22*x32*x41    
            +coeff[ 65]*x11*x21*x31*x42    
            +coeff[ 66]*x13*x21            
            +coeff[ 67]*x11*x21*x32*x41    
            +coeff[ 68]*x11        *x41    
            +coeff[ 69]*x11    *x31        
            +coeff[ 70]        *x33        
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 71]*x11    *x31*x41    
            +coeff[ 72]    *x21*x32*x41    
            +coeff[ 73]    *x23*x31        
            +coeff[ 74]*x11*x22*x31        
            +coeff[ 75]    *x22*x31    *x51
            +coeff[ 76]*x12*x21    *x41    
            +coeff[ 77]*x12    *x32        
            +coeff[ 78]*x11*x21*x31    *x51
            +coeff[ 79]    *x22    *x42*x51
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 80]    *x24*x31        
            +coeff[ 81]    *x22*x31*x41*x51
            +coeff[ 82]*x11*x21    *x42*x51
            +coeff[ 83]*x11*x23*x31        
            +coeff[ 84]        *x34    *x51
            +coeff[ 85]    *x24    *x42    
            +coeff[ 86]            *x42*x53
            +coeff[ 87]*x12*x22*x31        
            +coeff[ 88]    *x24*x31*x41    
        ;
        v_y_l5p65_ep10_q1en               =v_y_l5p65_ep10_q1en               
            +coeff[ 89]*x11*x23    *x42    
            +coeff[ 90]*x11*x23*x31*x41    
            +coeff[ 91]*x12*x22    *x42    
            +coeff[ 92]*x12*x22*x31*x41    
            +coeff[ 93]*x11*x22*x34        
            +coeff[ 94]*x11            *x51
            +coeff[ 95]            *x44    
            +coeff[ 96]        *x31*x43    
            ;
     
        return v_y_l5p65_ep10_q1en               ;
    }
    float p_l5p65_ep10_q1en               (float *x,int m){
        int ncoeff= 46;
        float avdat=  0.5410907E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46818E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 47]={
            -0.71830088E-02, 0.53594355E-03, 0.21862586E-02, 0.38543701E-01,
            -0.60068988E-02, 0.91120880E-02,-0.38519406E-02, 0.38406814E-02,
            -0.57410570E-02,-0.24523146E-02,-0.26473901E-02, 0.30354297E-03,
            -0.19468390E-02, 0.47140117E-03,-0.91799052E-03, 0.38878908E-02,
             0.10141003E-03, 0.12314295E-03,-0.39365102E-03,-0.40711756E-03,
            -0.37718014E-03,-0.25198818E-03, 0.15849640E-02, 0.43460681E-04,
            -0.47061487E-04,-0.16586373E-03,-0.21342079E-03, 0.47108982E-03,
            -0.14000183E-02, 0.22135493E-02,-0.16031123E-03, 0.20585039E-03,
            -0.11258014E-02, 0.95219049E-03,-0.35217608E-03, 0.15891223E-04,
            -0.43254891E-04, 0.30109068E-03,-0.10713492E-04, 0.30204694E-03,
             0.17762381E-03, 0.31256941E-03,-0.70792230E-04, 0.20115965E-03,
            -0.72219525E-04, 0.14883802E-03,
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
        float x51 = x5;
        float x52 = x51*x5;
     
    //                 function
     
        float v_p_l5p65_ep10_q1en               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]            *x42    
            +coeff[  7]*x11*x21            
        ;
        v_p_l5p65_ep10_q1en               =v_p_l5p65_ep10_q1en               
            +coeff[  8]    *x22    *x41    
            +coeff[  9]        *x31*x41    
            +coeff[ 10]*x11*x21    *x41    
            +coeff[ 11]                *x52
            +coeff[ 12]    *x22*x31        
            +coeff[ 13]*x12                
            +coeff[ 14]*x11*x21*x31        
            +coeff[ 15]    *x22    *x42    
            +coeff[ 16]    *x24*x31*x41    
        ;
        v_p_l5p65_ep10_q1en               =v_p_l5p65_ep10_q1en               
            +coeff[ 17]*x11                
            +coeff[ 18]        *x32        
            +coeff[ 19]    *x21    *x42    
            +coeff[ 20]    *x22        *x51
            +coeff[ 21]*x12        *x41    
            +coeff[ 22]*x11*x21    *x42    
            +coeff[ 23]*x11*x23*x31*x41    
            +coeff[ 24]    *x21    *x41    
            +coeff[ 25]            *x41*x51
        ;
        v_p_l5p65_ep10_q1en               =v_p_l5p65_ep10_q1en               
            +coeff[ 26]    *x21*x31*x41    
            +coeff[ 27]            *x43    
            +coeff[ 28]    *x24            
            +coeff[ 29]    *x22*x31*x41    
            +coeff[ 30]*x11*x21        *x51
            +coeff[ 31]        *x31*x43    
            +coeff[ 32]*x11*x23            
            +coeff[ 33]*x11*x21*x31*x41    
            +coeff[ 34]*x12*x22            
        ;
        v_p_l5p65_ep10_q1en               =v_p_l5p65_ep10_q1en               
            +coeff[ 35]        *x33*x42    
            +coeff[ 36]        *x31    *x51
            +coeff[ 37]    *x23            
            +coeff[ 38]        *x33        
            +coeff[ 39]        *x31*x42    
            +coeff[ 40]*x11*x22            
            +coeff[ 41]    *x22*x32        
            +coeff[ 42]*x11        *x42    
            +coeff[ 43]        *x32*x42    
        ;
        v_p_l5p65_ep10_q1en               =v_p_l5p65_ep10_q1en               
            +coeff[ 44]*x12    *x31        
            +coeff[ 45]*x11*x21*x32        
            ;
     
        return v_p_l5p65_ep10_q1en               ;
    }
    float l_l5p65_ep10_q1en               (float *x,int m){
        int ncoeff= 44;
        float avdat= -0.1808039E-04;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46818E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 45]={
            -0.17040178E-03, 0.89816953E-04, 0.18607065E-02, 0.63168122E-02,
            -0.11683669E-03,-0.19323246E-02,-0.12496021E-03,-0.14221091E-02,
            -0.75452082E-03, 0.14285410E-03, 0.18001682E-03,-0.16473743E-03,
            -0.26056534E-03, 0.10179533E-03, 0.35806306E-03, 0.24609199E-04,
            -0.50883191E-05, 0.58263860E-04,-0.52431911E-04,-0.14414925E-03,
            -0.11184438E-03, 0.12316115E-03, 0.15311378E-04, 0.14356751E-04,
            -0.20793379E-04,-0.84986013E-05,-0.31013551E-04,-0.12094562E-03,
             0.14875299E-05, 0.36067704E-05,-0.12313013E-04, 0.91075572E-05,
            -0.84330832E-05, 0.54508901E-05,-0.25272970E-04, 0.16227391E-03,
             0.48369344E-04, 0.14827811E-03,-0.12601745E-04,-0.57076268E-04,
             0.10520942E-03,-0.34237317E-04, 0.32582422E-04,-0.17923441E-04,
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
     
        float v_l_l5p65_ep10_q1en               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]        *x31*x41    
            +coeff[  7]            *x42    
        ;
        v_l_l5p65_ep10_q1en               =v_l_l5p65_ep10_q1en               
            +coeff[  8]    *x22    *x41    
            +coeff[  9]            *x41*x51
            +coeff[ 10]*x11*x21            
            +coeff[ 11]    *x22*x31        
            +coeff[ 12]*x11*x21    *x41    
            +coeff[ 13]            *x43    
            +coeff[ 14]    *x22    *x42    
            +coeff[ 15]    *x24*x31*x41    
            +coeff[ 16]    *x21    *x41    
        ;
        v_l_l5p65_ep10_q1en               =v_l_l5p65_ep10_q1en               
            +coeff[ 17]        *x31*x42    
            +coeff[ 18]*x11*x21*x31        
            +coeff[ 19]    *x24            
            +coeff[ 20]*x11*x23            
            +coeff[ 21]*x11*x21    *x42    
            +coeff[ 22]*x11*x23*x31*x41    
            +coeff[ 23]*x12                
            +coeff[ 24]*x12        *x41    
            +coeff[ 25]        *x32    *x52
        ;
        v_l_l5p65_ep10_q1en               =v_l_l5p65_ep10_q1en               
            +coeff[ 26]*x12*x22            
            +coeff[ 27]    *x22    *x43    
            +coeff[ 28]*x11                
            +coeff[ 29]        *x31    *x51
            +coeff[ 30]    *x21    *x42    
            +coeff[ 31]    *x22        *x51
            +coeff[ 32]            *x41*x52
            +coeff[ 33]*x11*x21        *x51
            +coeff[ 34]    *x23    *x41    
        ;
        v_l_l5p65_ep10_q1en               =v_l_l5p65_ep10_q1en               
            +coeff[ 35]    *x22*x31*x41    
            +coeff[ 36]*x11*x21*x31*x41    
            +coeff[ 37]    *x24    *x41    
            +coeff[ 38]    *x23*x31*x41    
            +coeff[ 39]    *x22*x31*x42    
            +coeff[ 40]*x11*x23    *x41    
            +coeff[ 41]*x11*x21    *x43    
            +coeff[ 42]*x12*x22    *x41    
            +coeff[ 43]*x11*x24    *x41    
            ;
     
        return v_l_l5p65_ep10_q1en               ;
    }
    float x_l5p65_ep13_q1ex               (float *x,int m){
        int ncoeff= 45;
        float avdat= -0.3973688E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46641E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 46]={
            -0.54735900E-03, 0.48707849E-02, 0.10872995E-03, 0.11734181E+00,
             0.19871345E-01,-0.30795528E-03, 0.45566307E-02, 0.29733391E-02,
             0.70107910E-02, 0.81624696E-03, 0.16413331E-02, 0.52566710E-02,
            -0.70142462E-02, 0.42306754E-04, 0.74225631E-04, 0.28024853E-03,
             0.38722792E-03, 0.62819687E-04, 0.91385690E-03, 0.37700757E-02,
            -0.49076919E-02,-0.73324507E-02,-0.22585586E-04,-0.14148365E-02,
            -0.86606387E-03, 0.11928389E-02,-0.44364305E-02,-0.23114202E-02,
            -0.38750801E-04, 0.49485138E-03,-0.96245256E-03,-0.41590593E-03,
             0.39434165E-03,-0.11066620E-02,-0.14909542E-02, 0.23043484E-02,
             0.85599319E-03, 0.96475807E-04,-0.12741220E-03, 0.17009824E-03,
            -0.17061408E-03,-0.14760786E-03, 0.17969561E-02, 0.28804704E-03,
             0.11184139E-02,
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
     
        float v_x_l5p65_ep13_q1ex               =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]                *x51
            +coeff[  3]    *x21            
            +coeff[  4]    *x21    *x41    
            +coeff[  5]*x12*x21*x31        
            +coeff[  6]*x11        *x41    
            +coeff[  7]    *x21        *x51
        ;
        v_x_l5p65_ep13_q1ex               =v_x_l5p65_ep13_q1ex               
            +coeff[  8]    *x21*x31        
            +coeff[  9]            *x41    
            +coeff[ 10]*x11    *x31        
            +coeff[ 11]    *x23            
            +coeff[ 12]    *x21    *x42    
            +coeff[ 13]*x13*x22            
            +coeff[ 14]    *x23*x31*x41    
            +coeff[ 15]        *x31        
            +coeff[ 16]*x11            *x51
        ;
        v_x_l5p65_ep13_q1ex               =v_x_l5p65_ep13_q1ex               
            +coeff[ 17]*x13                
            +coeff[ 18]*x12*x21            
            +coeff[ 19]*x11*x22            
            +coeff[ 20]    *x21*x31*x41    
            +coeff[ 21]    *x23    *x41    
            +coeff[ 22]*x13*x22    *x41    
            +coeff[ 23]*x11        *x42    
            +coeff[ 24]    *x21*x32        
            +coeff[ 25]    *x22    *x41    
        ;
        v_x_l5p65_ep13_q1ex               =v_x_l5p65_ep13_q1ex               
            +coeff[ 26]*x11*x22    *x41    
            +coeff[ 27]    *x23*x31        
            +coeff[ 28]*x13    *x31*x41    
            +coeff[ 29]*x11*x21    *x41    
            +coeff[ 30]*x11    *x31*x41    
            +coeff[ 31]    *x21    *x41*x51
            +coeff[ 32]    *x22*x31        
            +coeff[ 33]*x12*x21    *x41    
            +coeff[ 34]*x11*x22*x31        
        ;
        v_x_l5p65_ep13_q1ex               =v_x_l5p65_ep13_q1ex               
            +coeff[ 35]    *x21    *x43    
            +coeff[ 36]    *x23*x31*x42    
            +coeff[ 37]    *x22            
            +coeff[ 38]    *x21        *x52
            +coeff[ 39]*x11*x21*x31        
            +coeff[ 40]*x11    *x32        
            +coeff[ 41]    *x21*x31    *x51
            +coeff[ 42]    *x21*x31*x42    
            +coeff[ 43]*x12*x21*x32*x41    
        ;
        v_x_l5p65_ep13_q1ex               =v_x_l5p65_ep13_q1ex               
            +coeff[ 44]    *x23*x32*x41    
            ;
     
        return v_x_l5p65_ep13_q1ex               ;
    }
    float t_l5p65_ep13_q1ex               (float *x,int m){
        int ncoeff= 50;
        float avdat= -0.7830739E-03;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46641E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 51]={
            -0.16004466E-03,-0.62297224E-02,-0.59132529E-02, 0.20540204E-03,
             0.95883530E-04, 0.16920891E-02, 0.49466793E-02, 0.22482872E-02,
             0.10177942E-02,-0.19865311E-02,-0.26805001E-04,-0.24633562E-05,
             0.36468505E-03, 0.10421140E-02, 0.78829663E-03,-0.11733904E-02,
            -0.21505489E-02, 0.41483196E-04, 0.63542138E-04, 0.23806319E-03,
             0.19288764E-03, 0.30498439E-03,-0.33695722E-03,-0.62498235E-03,
            -0.12578481E-02, 0.23542038E-04,-0.19870242E-04, 0.10355390E-03,
             0.12701849E-03,-0.21168792E-03, 0.15386246E-03,-0.95980824E-04,
            -0.37769898E-03,-0.24617134E-03, 0.59279997E-03,-0.14685524E-03,
             0.26724871E-04, 0.41652744E-04,-0.32678778E-04,-0.20206340E-03,
             0.30644085E-04,-0.75230237E-04, 0.59042149E-03, 0.21985008E-04,
            -0.20807802E-04, 0.14783157E-04, 0.17906161E-04,-0.25191970E-04,
             0.15412530E-03, 0.45995643E-04,
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
     
        float v_t_l5p65_ep13_q1ex               =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]            *x41    
            +coeff[  4]    *x22            
            +coeff[  5]    *x21*x31        
            +coeff[  6]    *x21    *x41    
            +coeff[  7]    *x21        *x51
        ;
        v_t_l5p65_ep13_q1ex               =v_t_l5p65_ep13_q1ex               
            +coeff[  8]    *x23            
            +coeff[  9]    *x21    *x42    
            +coeff[ 10]*x13        *x41    
            +coeff[ 11]*x13            *x51
            +coeff[ 12]*x11    *x31        
            +coeff[ 13]*x11        *x41    
            +coeff[ 14]*x11*x22            
            +coeff[ 15]    *x21*x31*x41    
            +coeff[ 16]    *x23    *x41    
        ;
        v_t_l5p65_ep13_q1ex               =v_t_l5p65_ep13_q1ex               
            +coeff[ 17]*x13*x22    *x41    
            +coeff[ 18]        *x31        
            +coeff[ 19]*x11            *x51
            +coeff[ 20]*x12*x21            
            +coeff[ 21]    *x22    *x41    
            +coeff[ 22]*x11        *x42    
            +coeff[ 23]    *x23*x31        
            +coeff[ 24]*x11*x22    *x41    
            +coeff[ 25]    *x23*x32        
        ;
        v_t_l5p65_ep13_q1ex               =v_t_l5p65_ep13_q1ex               
            +coeff[ 26]                *x51
            +coeff[ 27]    *x22*x31        
            +coeff[ 28]*x11*x21    *x41    
            +coeff[ 29]*x11    *x31*x41    
            +coeff[ 30]    *x21    *x41*x51
            +coeff[ 31]    *x21        *x52
            +coeff[ 32]*x11*x22*x31        
            +coeff[ 33]*x12*x21    *x41    
            +coeff[ 34]    *x21    *x43    
        ;
        v_t_l5p65_ep13_q1ex               =v_t_l5p65_ep13_q1ex               
            +coeff[ 35]    *x23*x31*x42    
            +coeff[ 36]*x11*x21            
            +coeff[ 37]*x11*x21*x31        
            +coeff[ 38]*x11    *x32        
            +coeff[ 39]    *x21*x32        
            +coeff[ 40]    *x21*x31    *x51
            +coeff[ 41]*x12*x21*x31        
            +coeff[ 42]    *x21*x31*x42    
            +coeff[ 43]    *x23*x32*x41    
        ;
        v_t_l5p65_ep13_q1ex               =v_t_l5p65_ep13_q1ex               
            +coeff[ 44]            *x42    
            +coeff[ 45]*x13                
            +coeff[ 46]*x11        *x41*x51
            +coeff[ 47]*x12    *x31*x41    
            +coeff[ 48]    *x21*x32*x41    
            +coeff[ 49]*x11        *x43    
            ;
     
        return v_t_l5p65_ep13_q1ex               ;
    }
    float y_l5p65_ep13_q1ex               (float *x,int m){
        int ncoeff= 97;
        float avdat=  0.1137931E-01;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46641E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 98]={
            -0.10541648E-01, 0.17392212E+00, 0.15156624E-02, 0.27881708E-01,
            -0.17809393E-01,-0.10781345E-01, 0.26964515E-01, 0.11786405E-01,
            -0.76162913E-02,-0.27014297E-02,-0.17415486E-01,-0.81376228E-02,
             0.13724159E-02,-0.59872414E-02,-0.30255464E-02, 0.34993561E-03,
            -0.13645780E-02, 0.10742816E-02, 0.11059014E-01, 0.74198050E-02,
            -0.64782356E-03,-0.12827786E-02, 0.46878326E-03,-0.13384436E-02,
            -0.10239559E-02,-0.60973968E-03, 0.47600879E-02,-0.41224221E-02,
             0.33204614E-02,-0.34774721E-02, 0.20464193E-02, 0.22778560E-02,
            -0.89489471E-03, 0.95083989E-03, 0.62103663E-03,-0.35557602E-03,
             0.12769379E-02,-0.11772776E-02, 0.46635876E-04, 0.62204275E-03,
            -0.25010199E-03, 0.36071212E-03, 0.62093948E-03, 0.13252023E-03,
            -0.66429748E-04,-0.15423895E-03, 0.66664034E-04, 0.24819665E-03,
             0.10912504E-03, 0.13458160E-03,-0.66588451E-04, 0.63786790E-03,
             0.48559374E-03, 0.28874938E-03,-0.21002085E-02, 0.35141277E-03,
             0.33536655E-03,-0.22690697E-02,-0.10744487E-02, 0.20985641E-02,
            -0.10451307E-02,-0.18320688E-03, 0.19910240E-02, 0.21232891E-04,
             0.66557957E-03,-0.13889424E-02,-0.12116837E-02,-0.15523152E-03,
             0.41452073E-04, 0.37682895E-03, 0.40665720E-03, 0.17616239E-03,
            -0.46906414E-03, 0.14951284E-03,-0.18589986E-03,-0.59148973E-04,
            -0.35970649E-03, 0.28715953E-04,-0.55351277E-03, 0.66812019E-04,
             0.13912914E-03, 0.34267239E-04,-0.42158112E-03, 0.12781281E-04,
            -0.27840291E-03, 0.43946563E-03, 0.15482556E-04,-0.23768084E-03,
            -0.28363219E-03,-0.16614480E-03, 0.10740517E-03,-0.47592319E-04,
             0.65240479E-03, 0.83132392E-04, 0.74673211E-04, 0.23993049E-04,
            -0.31413776E-04,
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
     
        float v_y_l5p65_ep13_q1ex               =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]        *x31        
            +coeff[  4]                *x51
            +coeff[  5]            *x42    
            +coeff[  6]    *x22            
            +coeff[  7]*x11*x21            
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[  8]        *x31*x41    
            +coeff[  9]            *x41*x51
            +coeff[ 10]    *x22    *x41    
            +coeff[ 11]*x11*x21    *x41    
            +coeff[ 12]*x12                
            +coeff[ 13]    *x22*x31        
            +coeff[ 14]*x11*x21*x31        
            +coeff[ 15]*x11                
            +coeff[ 16]        *x32        
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 17]                *x52
            +coeff[ 18]    *x22    *x42    
            +coeff[ 19]    *x22*x31*x41    
            +coeff[ 20]        *x31    *x51
            +coeff[ 21]    *x21    *x42    
            +coeff[ 22]            *x42*x51
            +coeff[ 23]    *x22        *x51
            +coeff[ 24]*x12        *x41    
            +coeff[ 25]*x11*x21        *x51
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 26]*x11*x21    *x42    
            +coeff[ 27]    *x24            
            +coeff[ 28]*x11*x21*x31*x41    
            +coeff[ 29]*x11*x23            
            +coeff[ 30]            *x43    
            +coeff[ 31]        *x31*x42    
            +coeff[ 32]    *x21*x31*x41    
            +coeff[ 33]    *x23            
            +coeff[ 34]*x11*x22            
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 35]*x12    *x31        
            +coeff[ 36]    *x22*x32        
            +coeff[ 37]*x12*x22            
            +coeff[ 38]        *x34*x41    
            +coeff[ 39]        *x32*x41    
            +coeff[ 40]*x11        *x42    
            +coeff[ 41]        *x31*x41*x51
            +coeff[ 42]*x11*x21*x32        
            +coeff[ 43]    *x21*x31        
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 44]    *x21        *x51
            +coeff[ 45]    *x21*x32        
            +coeff[ 46]        *x33        
            +coeff[ 47]        *x31*x43    
            +coeff[ 48]            *x41*x52
            +coeff[ 49]*x12*x21            
            +coeff[ 50]*x12            *x51
            +coeff[ 51]    *x22    *x41*x51
            +coeff[ 52]*x12        *x42    
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 53]    *x22*x31    *x51
            +coeff[ 54]    *x22    *x43    
            +coeff[ 55]*x12    *x31*x41    
            +coeff[ 56]*x11*x21    *x41*x51
            +coeff[ 57]    *x22*x31*x42    
            +coeff[ 58]*x11*x21    *x43    
            +coeff[ 59]    *x24    *x41    
            +coeff[ 60]*x11*x21*x31*x42    
            +coeff[ 61]*x13*x21            
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 62]*x11*x23    *x41    
            +coeff[ 63]*x11*x22*x31*x41    
            +coeff[ 64]*x12*x22    *x41    
            +coeff[ 65]    *x24*x32*x41    
            +coeff[ 66]*x11*x23*x32*x41    
            +coeff[ 67]*x11    *x31*x41    
            +coeff[ 68]        *x32    *x51
            +coeff[ 69]    *x21    *x43    
            +coeff[ 70]    *x21*x31*x42    
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 71]        *x32*x42    
            +coeff[ 72]    *x23    *x41    
            +coeff[ 73]    *x21*x32*x41    
            +coeff[ 74]    *x23*x31        
            +coeff[ 75]                *x53
            +coeff[ 76]*x11*x22    *x41    
            +coeff[ 77]*x11    *x32*x41    
            +coeff[ 78]        *x31*x44    
            +coeff[ 79]*x12    *x32        
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 80]*x11*x21*x31    *x51
            +coeff[ 81]*x13        *x41    
            +coeff[ 82]    *x22    *x42*x51
            +coeff[ 83]*x13    *x31        
            +coeff[ 84]    *x22*x31*x41*x51
            +coeff[ 85]*x11*x23*x31        
            +coeff[ 86]*x13            *x51
            +coeff[ 87]*x12*x22*x31        
            +coeff[ 88]        *x32*x45    
        ;
        v_y_l5p65_ep13_q1ex               =v_y_l5p65_ep13_q1ex               
            +coeff[ 89]*x11*x24*x31        
            +coeff[ 90]    *x24        *x52
            +coeff[ 91]*x13    *x32*x41    
            +coeff[ 92]*x12*x24*x31        
            +coeff[ 93]    *x21    *x41    
            +coeff[ 94]*x11        *x41    
            +coeff[ 95]*x11    *x31        
            +coeff[ 96]*x11            *x51
            ;
     
        return v_y_l5p65_ep13_q1ex               ;
    }
    float p_l5p65_ep13_q1ex               (float *x,int m){
        int ncoeff= 41;
        float avdat=  0.4687944E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46641E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 42]={
            -0.47518993E-02, 0.77725365E-03, 0.92262486E-02, 0.73699042E-01,
            -0.91919070E-02, 0.13730549E-01,-0.22153950E-02, 0.59519084E-02,
            -0.52755359E-02,-0.83670253E-02,-0.38448866E-02, 0.34299183E-02,
            -0.27927357E-04,-0.36145386E-02, 0.64548693E-03,-0.29998780E-02,
             0.72642323E-03,-0.14165735E-02, 0.54027331E-02, 0.16209661E-03,
            -0.61976677E-03,-0.43609191E-03,-0.84193528E-03, 0.44941396E-03,
            -0.20935251E-02,-0.35705703E-03,-0.17243336E-02,-0.41819963E-03,
             0.21005012E-02, 0.10594374E-02, 0.52238158E-02, 0.52439712E-03,
            -0.43124743E-02,-0.48062956E-03, 0.64902700E-03, 0.30058849E-03,
             0.44259045E-03,-0.14048163E-03,-0.57015778E-03,-0.53966703E-03,
             0.27980862E-03,
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
        float x34 = x33*x3;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x44 = x43*x4;
        float x45 = x44*x4;
        float x51 = x5;
        float x52 = x51*x5;
     
    //                 function
     
        float v_p_l5p65_ep13_q1ex               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]            *x41*x51
            +coeff[  7]*x11*x21            
        ;
        v_p_l5p65_ep13_q1ex               =v_p_l5p65_ep13_q1ex               
            +coeff[  8]            *x42    
            +coeff[  9]    *x22    *x41    
            +coeff[ 10]*x11*x21    *x41    
            +coeff[ 11]    *x22*x31*x41    
            +coeff[ 12]        *x33*x41    
            +coeff[ 13]        *x31*x41    
            +coeff[ 14]                *x52
            +coeff[ 15]    *x22*x31        
            +coeff[ 16]*x12                
        ;
        v_p_l5p65_ep13_q1ex               =v_p_l5p65_ep13_q1ex               
            +coeff[ 17]*x11*x21*x31        
            +coeff[ 18]    *x22    *x42    
            +coeff[ 19]*x11                
            +coeff[ 20]        *x32        
            +coeff[ 21]        *x31    *x51
            +coeff[ 22]    *x22        *x51
            +coeff[ 23]    *x23            
            +coeff[ 24]    *x24            
            +coeff[ 25]*x11*x21        *x51
        ;
        v_p_l5p65_ep13_q1ex               =v_p_l5p65_ep13_q1ex               
            +coeff[ 26]*x11*x23            
            +coeff[ 27]*x12        *x41    
            +coeff[ 28]*x11*x21    *x42    
            +coeff[ 29]            *x45    
            +coeff[ 30]*x11*x23*x31*x41    
            +coeff[ 31]        *x33*x44    
            +coeff[ 32]*x11*x25*x31*x41    
            +coeff[ 33]    *x21    *x42    
            +coeff[ 34]        *x31*x42    
        ;
        v_p_l5p65_ep13_q1ex               =v_p_l5p65_ep13_q1ex               
            +coeff[ 35]*x11*x22            
            +coeff[ 36]    *x22*x32        
            +coeff[ 37]*x12    *x31        
            +coeff[ 38]    *x23*x31*x41    
            +coeff[ 39]*x12*x22            
            +coeff[ 40]        *x34*x41    
            ;
     
        return v_p_l5p65_ep13_q1ex               ;
    }
    float l_l5p65_ep13_q1ex               (float *x,int m){
        int ncoeff= 86;
        float avdat= -0.1160080E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.46641E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 87]={
             0.12021636E-02, 0.12989946E-03, 0.18602116E-02, 0.58440911E-02,
            -0.10818230E-03,-0.23095510E-02,-0.86793018E-03,-0.44107279E-02,
             0.97336415E-04, 0.96587883E-03, 0.20935181E-03,-0.21211738E-02,
            -0.76404930E-03,-0.38760999E-03,-0.64683598E-04,-0.60728111E-04,
             0.54506958E-03, 0.16202667E-03, 0.99315759E-04,-0.11574693E-03,
             0.86489155E-04,-0.27746754E-03, 0.12338255E-02,-0.21096805E-03,
             0.50153740E-03, 0.58493923E-04,-0.85806772E-04,-0.55711182E-04,
             0.37016362E-05, 0.43489056E-03,-0.81041951E-04,-0.73922849E-04,
             0.66616235E-03, 0.22301970E-03,-0.59842190E-04, 0.53143507E-03,
             0.36429628E-04, 0.40214945E-03,-0.11123943E-04, 0.11646686E-04,
            -0.80281374E-04,-0.43062792E-04,-0.62272907E-03, 0.10517170E-03,
            -0.74440024E-04, 0.81652266E-04,-0.11087807E-04, 0.64432737E-04,
             0.61405110E-04, 0.40285104E-04, 0.13264794E-03,-0.31386787E-05,
            -0.13584321E-03,-0.46365760E-03, 0.90417692E-04,-0.24781708E-03,
             0.16307573E-04,-0.21672776E-03, 0.12023464E-04,-0.53346939E-04,
            -0.33586018E-05,-0.37970242E-05, 0.37577486E-05,-0.10470465E-04,
            -0.55446476E-05,-0.10104870E-04, 0.48242159E-05, 0.34859426E-04,
            -0.68085610E-04,-0.31872602E-04, 0.74060599E-05,-0.14819970E-04,
            -0.21658581E-04, 0.19302079E-04,-0.11714329E-04,-0.10865716E-04,
             0.12878949E-04,-0.51525316E-04, 0.28136743E-04,-0.16216464E-03,
            -0.54608667E-04,-0.76581397E-04,-0.83218503E-04,-0.13468949E-03,
            -0.25085101E-04,-0.85736334E-04,
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
     
        float v_l_l5p65_ep13_q1ex               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]        *x31*x41    
            +coeff[  7]            *x42    
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[  8]        *x31    *x51
            +coeff[  9]            *x41*x51
            +coeff[ 10]*x11*x21            
            +coeff[ 11]    *x22    *x41    
            +coeff[ 12]*x11*x21    *x41    
            +coeff[ 13]    *x22*x31        
            +coeff[ 14]    *x21    *x41    
            +coeff[ 15]                *x52
            +coeff[ 16]            *x43    
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 17]    *x22        *x51
            +coeff[ 18]            *x42*x51
            +coeff[ 19]*x11*x21*x31        
            +coeff[ 20]*x11*x21        *x51
            +coeff[ 21]    *x24            
            +coeff[ 22]    *x22    *x42    
            +coeff[ 23]*x11*x23            
            +coeff[ 24]*x11*x21    *x42    
            +coeff[ 25]        *x33*x42    
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 26]    *x24*x31*x41    
            +coeff[ 27]        *x32        
            +coeff[ 28]        *x33        
            +coeff[ 29]        *x31*x42    
            +coeff[ 30]            *x41*x52
            +coeff[ 31]*x12        *x41    
            +coeff[ 32]    *x22*x31*x41    
            +coeff[ 33]*x11*x21*x31*x41    
            +coeff[ 34]*x12*x22            
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 35]    *x24    *x41    
            +coeff[ 36]        *x34*x41    
            +coeff[ 37]*x11*x23    *x41    
            +coeff[ 38]*x11                
            +coeff[ 39]    *x21        *x51
            +coeff[ 40]    *x23    *x41    
            +coeff[ 41]*x11*x22    *x41    
            +coeff[ 42]    *x22    *x43    
            +coeff[ 43]*x12*x22    *x41    
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 44]    *x22*x33*x42    
            +coeff[ 45]        *x32*x41    
            +coeff[ 46]        *x31    *x52
            +coeff[ 47]    *x22*x32        
            +coeff[ 48]    *x21    *x43    
            +coeff[ 49]*x12        *x42    
            +coeff[ 50]    *x24*x31        
            +coeff[ 51]    *x22*x33        
            +coeff[ 52]    *x22*x32*x41    
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 53]    *x22*x31*x42    
            +coeff[ 54]*x11*x23*x31        
            +coeff[ 55]*x11*x21    *x43    
            +coeff[ 56]*x12            *x53
            +coeff[ 57]    *x24    *x42    
            +coeff[ 58]    *x24*x31    *x51
            +coeff[ 59]*x11*x23*x31*x42    
            +coeff[ 60]    *x21*x31        
            +coeff[ 61]*x11        *x41    
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 62]*x11            *x51
            +coeff[ 63]*x12                
            +coeff[ 64]        *x32    *x51
            +coeff[ 65]    *x23*x31        
            +coeff[ 66]    *x21*x32*x41    
            +coeff[ 67]    *x21*x31*x42    
            +coeff[ 68]            *x44    
            +coeff[ 69]    *x22*x31    *x51
            +coeff[ 70]        *x32*x41*x51
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 71]            *x43*x51
            +coeff[ 72]    *x22        *x52
            +coeff[ 73]*x11*x21*x32        
            +coeff[ 74]*x11*x21*x31    *x51
            +coeff[ 75]*x11*x21        *x52
            +coeff[ 76]*x12    *x31*x41    
            +coeff[ 77]        *x31*x44    
            +coeff[ 78]    *x22*x31*x41*x51
            +coeff[ 79]*x11*x21*x31*x42    
        ;
        v_l_l5p65_ep13_q1ex               =v_l_l5p65_ep13_q1ex               
            +coeff[ 80]        *x34*x42    
            +coeff[ 81]    *x22*x31*x43    
            +coeff[ 82]        *x33*x43    
            +coeff[ 83]*x11*x23    *x42    
            +coeff[ 84]    *x23*x33*x41    
            +coeff[ 85]*x11*x23*x32*x41    
            ;
     
        return v_l_l5p65_ep13_q1ex               ;
    }
    float x_l5p65_ep24_dex                (float *x,int m){
        int ncoeff= 60;
        float avdat= -0.1975655E-01;
        float xmin[10]={
            -0.14981E-01,-0.47334E-01,-0.14977E-01,-0.32322E-01,-0.49995E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.47289E-01, 0.14991E-01, 0.30724E-01, 0.49882E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 61]={
            -0.46704631E-03,-0.50400790E-01, 0.13257769E+00, 0.27184528E+00,
             0.35015885E-01, 0.11692772E+00,-0.38376451E-03,-0.18139366E-02,
             0.58537080E-05,-0.13454748E-01, 0.58842571E-02, 0.85563054E-02,
             0.25461592E-01, 0.96878838E-02, 0.38765475E-01, 0.29684685E-01,
            -0.55267200E-01, 0.74029755E-03, 0.12464526E-02, 0.17200217E-02,
            -0.34662818E-02,-0.55778073E-02, 0.19035220E-01,-0.31291954E-01,
            -0.47159564E-01, 0.44911663E-03, 0.47432710E-02,-0.75617358E-02,
             0.57917056E-02,-0.48913998E-02, 0.81758136E-02,-0.28432578E-01,
            -0.34500493E-03,-0.50450879E-03, 0.70387329E-03,-0.16403265E-02,
             0.16008230E-02,-0.11518778E-02,-0.48116422E-02, 0.11263527E-02,
             0.32918158E-02,-0.54751746E-02,-0.81403181E-02, 0.13552934E-01,
            -0.12150554E-03,-0.47570943E-05,-0.10023019E-02, 0.10557071E-02,
            -0.85832441E-03, 0.22204749E-02, 0.23989205E-02, 0.13076355E-01,
             0.49683836E-03, 0.36923279E-03, 0.49668684E-03,-0.81039290E-03,
             0.10898081E-02, 0.11081074E-02,-0.10184695E-02, 0.39587626E-02,
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
     
        float v_x_l5p65_ep24_dex                =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]                *x51
            +coeff[  3]    *x21            
            +coeff[  4]    *x21        *x51
            +coeff[  5]    *x21    *x41    
            +coeff[  6]*x13        *x41    
            +coeff[  7]*x12*x21*x31        
        ;
        v_x_l5p65_ep24_dex                =v_x_l5p65_ep24_dex                
            +coeff[  8]    *x21*x31    *x52
            +coeff[  9]    *x23*x31        
            +coeff[ 10]            *x41    
            +coeff[ 11]*x11    *x31        
            +coeff[ 12]*x11        *x41    
            +coeff[ 13]    *x22            
            +coeff[ 14]    *x21*x31        
            +coeff[ 15]    *x23            
            +coeff[ 16]    *x21    *x42    
        ;
        v_x_l5p65_ep24_dex                =v_x_l5p65_ep24_dex                
            +coeff[ 17]*x13*x22            
            +coeff[ 18]    *x23*x31*x41    
            +coeff[ 19]        *x31        
            +coeff[ 20]                *x52
            +coeff[ 21]            *x42    
            +coeff[ 22]*x11*x22            
            +coeff[ 23]    *x21*x31*x41    
            +coeff[ 24]    *x23    *x41    
            +coeff[ 25]*x13*x22    *x41    
        ;
        v_x_l5p65_ep24_dex                =v_x_l5p65_ep24_dex                
            +coeff[ 26]*x12*x21            
            +coeff[ 27]*x11        *x42    
            +coeff[ 28]    *x21    *x41*x51
            +coeff[ 29]    *x21*x32        
            +coeff[ 30]    *x22    *x41    
            +coeff[ 31]*x11*x22    *x41    
            +coeff[ 32]*x13    *x31*x41    
            +coeff[ 33]*x11*x22*x33        
            +coeff[ 34]*x11            *x51
        ;
        v_x_l5p65_ep24_dex                =v_x_l5p65_ep24_dex                
            +coeff[ 35]*x11*x21            
            +coeff[ 36]            *x41*x51
            +coeff[ 37]        *x31*x41    
            +coeff[ 38]*x11    *x31*x41    
            +coeff[ 39]    *x21*x31    *x51
            +coeff[ 40]    *x22*x31        
            +coeff[ 41]*x12*x21    *x41    
            +coeff[ 42]*x11*x22*x31        
            +coeff[ 43]    *x21    *x43    
        ;
        v_x_l5p65_ep24_dex                =v_x_l5p65_ep24_dex                
            +coeff[ 44]*x13*x21    *x41    
            +coeff[ 45]    *x23*x31*x42    
            +coeff[ 46]    *x21        *x52
            +coeff[ 47]*x11*x21*x31        
            +coeff[ 48]*x11    *x32        
            +coeff[ 49]*x11*x21    *x41    
            +coeff[ 50]    *x23        *x51
            +coeff[ 51]    *x21*x31*x42    
            +coeff[ 52]*x12*x21*x32*x41    
        ;
        v_x_l5p65_ep24_dex                =v_x_l5p65_ep24_dex                
            +coeff[ 53]    *x23*x32*x41    
            +coeff[ 54]*x11        *x41*x51
            +coeff[ 55]            *x42*x51
            +coeff[ 56]*x11*x22        *x51
            +coeff[ 57]*x11*x23            
            +coeff[ 58]    *x22    *x41*x51
            +coeff[ 59]    *x21*x32*x41    
            ;
     
        return v_x_l5p65_ep24_dex                ;
    }
    float t_l5p65_ep24_dex                (float *x,int m){
        int ncoeff= 80;
        float avdat=  0.5976589E+00;
        float xmin[10]={
            -0.14981E-01,-0.47334E-01,-0.14977E-01,-0.32322E-01,-0.49995E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.47289E-01, 0.14991E-01, 0.30724E-01, 0.49882E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 81]={
             0.69969735E-03, 0.11757700E-01,-0.85240833E-01,-0.12401777E-02,
             0.30670952E-01,-0.32318313E-01,-0.64806310E-02, 0.47453155E-03,
             0.36431730E-02, 0.17626150E-03,-0.23793683E-02,-0.10679368E-01,
            -0.70100329E-02,-0.84430724E-02, 0.14914857E-01, 0.81976556E-04,
            -0.61060977E-03,-0.47716484E-03,-0.75301470E-03,-0.13187689E-02,
            -0.54009426E-02, 0.82397442E-02, 0.12996060E-01,-0.18769646E-04,
             0.29496729E-03, 0.11565730E-02,-0.13251409E-02, 0.13074091E-02,
            -0.29624698E-02, 0.21168441E-02,-0.17267818E-02, 0.78560552E-02,
            -0.55183441E-03, 0.48686139E-04, 0.11566786E-03, 0.22484188E-03,
            -0.14535222E-02,-0.10190824E-02, 0.14427203E-02, 0.48644131E-03,
            -0.32043798E-03, 0.48284561E-03, 0.65399444E-03, 0.23049866E-02,
             0.15435790E-02,-0.48119546E-03, 0.22656220E-03,-0.37185799E-02,
            -0.37912829E-02,-0.11533130E-02,-0.14651089E-03,-0.15690511E-03,
            -0.51492840E-04,-0.99072524E-04,-0.13722894E-03, 0.66182890E-03,
             0.13833669E-03, 0.14107778E-03, 0.14681152E-03, 0.40578988E-03,
            -0.10741001E-02, 0.13690686E-02,-0.41776188E-03, 0.38385784E-03,
             0.77368648E-04, 0.16430632E-03,-0.51009678E-03, 0.79900544E-03,
            -0.10283161E-03,-0.65109378E-03,-0.20463945E-03,-0.54543594E-04,
            -0.11143721E-03,-0.10548545E-04,-0.83522762E-04, 0.30172314E-03,
             0.39435650E-04,-0.10544739E-03, 0.43007924E-03,-0.22659668E-03,
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
     
        float v_t_l5p65_ep24_dex                =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x21    *x41    
            +coeff[  6]    *x21        *x51
            +coeff[  7]*x12*x21*x31        
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[  8]    *x23*x31        
            +coeff[  9]*x13        *x41    
            +coeff[ 10]*x11    *x31        
            +coeff[ 11]    *x21*x31        
            +coeff[ 12]*x11        *x41    
            +coeff[ 13]    *x23            
            +coeff[ 14]    *x21    *x42    
            +coeff[ 15]*x13*x22            
            +coeff[ 16]    *x23*x31*x41    
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[ 17]        *x31        
            +coeff[ 18]*x11            *x51
            +coeff[ 19]                *x52
            +coeff[ 20]*x11*x22            
            +coeff[ 21]    *x21*x31*x41    
            +coeff[ 22]    *x23    *x41    
            +coeff[ 23]*x13*x22    *x41    
            +coeff[ 24]        *x31*x41    
            +coeff[ 25]            *x41*x51
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[ 26]*x12*x21            
            +coeff[ 27]    *x21*x32        
            +coeff[ 28]    *x22    *x41    
            +coeff[ 29]*x11        *x42    
            +coeff[ 30]            *x42*x51
            +coeff[ 31]*x11*x22    *x41    
            +coeff[ 32]            *x42*x52
            +coeff[ 33]*x13    *x31*x41    
            +coeff[ 34]*x11*x22*x33        
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[ 35]*x11*x21            
            +coeff[ 36]            *x42    
            +coeff[ 37]*x11*x21    *x41    
            +coeff[ 38]*x11    *x31*x41    
            +coeff[ 39]    *x22        *x51
            +coeff[ 40]    *x21    *x41*x51
            +coeff[ 41]    *x21        *x52
            +coeff[ 42]            *x41*x52
            +coeff[ 43]*x11*x22*x31        
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[ 44]*x12*x21    *x41    
            +coeff[ 45]    *x22*x31        
            +coeff[ 46]*x11    *x32        
            +coeff[ 47]    *x21*x31*x42    
            +coeff[ 48]    *x21    *x43    
            +coeff[ 49]    *x22    *x41*x51
            +coeff[ 50]*x12*x21*x32*x41    
            +coeff[ 51]    *x23*x32*x41    
            +coeff[ 52]*x12                
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[ 53]        *x31    *x51
            +coeff[ 54]*x11*x21*x31        
            +coeff[ 55]            *x43    
            +coeff[ 56]*x11*x21        *x51
            +coeff[ 57]    *x21*x31    *x51
            +coeff[ 58]        *x31*x41*x51
            +coeff[ 59]    *x22*x31*x41    
            +coeff[ 60]    *x21*x32*x41    
            +coeff[ 61]    *x22    *x42    
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[ 62]*x11*x21    *x41*x51
            +coeff[ 63]    *x21    *x42*x51
            +coeff[ 64]        *x33*x42    
            +coeff[ 65]*x13        *x41*x51
            +coeff[ 66]    *x23    *x41*x51
            +coeff[ 67]    *x22    *x42*x51
            +coeff[ 68]*x12*x22*x32        
            +coeff[ 69]    *x22*x31*x41*x52
            +coeff[ 70]    *x22            
        ;
        v_t_l5p65_ep24_dex                =v_t_l5p65_ep24_dex                
            +coeff[ 71]        *x32        
            +coeff[ 72]*x13                
            +coeff[ 73]        *x33        
            +coeff[ 74]*x12        *x41    
            +coeff[ 75]        *x31*x42    
            +coeff[ 76]*x11            *x52
            +coeff[ 77]*x11*x23            
            +coeff[ 78]*x11*x21    *x42    
            +coeff[ 79]*x11        *x43    
            ;
     
        return v_t_l5p65_ep24_dex                ;
    }
    float y_l5p65_ep24_dex                (float *x,int m){
        int ncoeff= 97;
        float avdat=  0.2622084E-02;
        float xmin[10]={
            -0.14981E-01,-0.47334E-01,-0.14977E-01,-0.32322E-01,-0.49995E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.47289E-01, 0.14991E-01, 0.30724E-01, 0.49882E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 98]={
            -0.92777582E-02, 0.86630195E-01, 0.40679942E-02,-0.14972846E-01,
            -0.28402276E-01,-0.15690865E-01,-0.48651297E-01,-0.10527425E-01,
             0.39503250E-01,-0.47673937E-02, 0.75270878E-02, 0.59516110E-01,
             0.15895586E-01, 0.70414520E-02, 0.76062884E-02,-0.18756671E-01,
            -0.52645663E-02,-0.48336565E-01,-0.80109444E-02,-0.10704917E-01,
             0.11054131E-01,-0.16883319E-02,-0.60906769E-02, 0.18724524E-02,
            -0.76180608E-02,-0.49216142E-02,-0.63446909E-02,-0.32369238E-02,
            -0.45572105E-02, 0.29988869E-02, 0.86854358E-03,-0.73680026E-03,
            -0.27066185E-02, 0.98827500E-02,-0.10413237E-01,-0.51960228E-02,
            -0.58457726E-02,-0.57828153E-03,-0.73282904E-03,-0.23698837E-02,
            -0.17843838E-02,-0.17316904E-02,-0.11629983E-02, 0.58932696E-02,
            -0.66672900E-03, 0.41348203E-02, 0.76512928E-03,-0.41525057E-02,
             0.33319946E-02, 0.82498202E-02, 0.13381924E-02, 0.32248229E-03,
            -0.35226572E-03,-0.15351673E-02, 0.35708811E-03,-0.48803125E-03,
             0.30935647E-02,-0.24174859E-02, 0.11579137E-02, 0.42290086E-03,
            -0.11038897E-02,-0.15211446E-02,-0.12004341E-02,-0.20946409E-02,
            -0.14553916E-02,-0.81893173E-03,-0.78583282E-03, 0.49161315E-02,
             0.57897223E-02, 0.61678194E-03,-0.64139436E-02, 0.18908427E-02,
             0.45284363E-04, 0.37615685E-03,-0.22990133E-03, 0.16973655E-03,
            -0.29120871E-03, 0.11805024E-03, 0.25494839E-03, 0.13485403E-02,
             0.96398720E-03, 0.68898406E-03, 0.51248923E-03, 0.73646993E-03,
             0.48831961E-03, 0.54196035E-03, 0.30921015E-02, 0.21186369E-02,
             0.21249854E-02, 0.97528267E-04, 0.96642826E-03,-0.49950567E-03,
             0.10036482E-02, 0.89916524E-04, 0.30211220E-03,-0.27836888E-03,
            -0.37615321E-03,
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
     
        float v_y_l5p65_ep24_dex                =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]        *x31        
            +coeff[  4]                *x51
            +coeff[  5]            *x42    
            +coeff[  6]    *x21    *x41    
            +coeff[  7]        *x31*x41    
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
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
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]    *x23            
            +coeff[ 19]    *x22*x31        
            +coeff[ 20]    *x22        *x51
            +coeff[ 21]        *x32        
            +coeff[ 22]            *x43    
            +coeff[ 23]*x12                
            +coeff[ 24]    *x21*x31*x41    
            +coeff[ 25]*x11        *x42    
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 26]*x11*x21    *x41    
            +coeff[ 27]*x11*x21*x31        
            +coeff[ 28]            *x41*x52
            +coeff[ 29]*x11*x21        *x51
            +coeff[ 30]*x11    *x31        
            +coeff[ 31]*x11            *x51
            +coeff[ 32]    *x21    *x41*x51
            +coeff[ 33]    *x21    *x43    
            +coeff[ 34]    *x24            
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 35]    *x22    *x41*x51
            +coeff[ 36]*x11*x23            
            +coeff[ 37]        *x31*x42    
            +coeff[ 38]            *x42*x51
            +coeff[ 39]*x11    *x31*x41    
            +coeff[ 40]*x11*x22            
            +coeff[ 41]*x12        *x41    
            +coeff[ 42]    *x22    *x42    
            +coeff[ 43]    *x21*x31*x42    
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 44]        *x31    *x52
            +coeff[ 45]    *x22*x31*x41    
            +coeff[ 46]                *x53
            +coeff[ 47]*x11*x22    *x41    
            +coeff[ 48]*x11*x21*x31*x41    
            +coeff[ 49]    *x23    *x42    
            +coeff[ 50]*x11*x21    *x44    
            +coeff[ 51]*x11                
            +coeff[ 52]    *x21*x32        
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 53]        *x31*x41*x51
            +coeff[ 54]*x12*x21            
            +coeff[ 55]*x12    *x31        
            +coeff[ 56]*x11*x21    *x42    
            +coeff[ 57]    *x21    *x42*x51
            +coeff[ 58]    *x22*x32        
            +coeff[ 59]*x12            *x51
            +coeff[ 60]    *x21*x31*x41*x51
            +coeff[ 61]    *x22*x31    *x51
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 62]*x12*x21    *x41    
            +coeff[ 63]*x11*x21    *x41*x51
            +coeff[ 64]*x12*x22            
            +coeff[ 65]*x11*x21*x31    *x51
            +coeff[ 66]    *x22        *x52
            +coeff[ 67]*x11*x22    *x42    
            +coeff[ 68]    *x22    *x42*x51
            +coeff[ 69]*x11*x21*x34        
            +coeff[ 70]    *x24    *x43    
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 71]    *x22*x32*x43    
            +coeff[ 72]    *x22*x33*x41*x51
            +coeff[ 73]        *x32*x41    
            +coeff[ 74]*x11    *x32        
            +coeff[ 75]    *x21*x31    *x51
            +coeff[ 76]        *x32    *x51
            +coeff[ 77]*x11    *x31    *x51
            +coeff[ 78]    *x21        *x52
            +coeff[ 79]            *x43*x51
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 80]    *x21*x32*x41    
            +coeff[ 81]        *x31*x42*x51
            +coeff[ 82]    *x23*x31        
            +coeff[ 83]*x12        *x42    
            +coeff[ 84]*x12    *x31*x41    
            +coeff[ 85]    *x21    *x41*x52
            +coeff[ 86]    *x23*x31*x41    
            +coeff[ 87]*x11*x22*x31*x41    
            +coeff[ 88]    *x22*x31*x41*x51
        ;
        v_y_l5p65_ep24_dex                =v_y_l5p65_ep24_dex                
            +coeff[ 89]*x11        *x41*x51
            +coeff[ 90]*x11        *x43    
            +coeff[ 91]    *x23    *x41    
            +coeff[ 92]*x11    *x31*x42    
            +coeff[ 93]*x11            *x52
            +coeff[ 94]*x11    *x32*x41    
            +coeff[ 95]            *x42*x52
            +coeff[ 96]*x11*x22*x31        
            ;
     
        return v_y_l5p65_ep24_dex                ;
    }
    float p_l5p65_ep24_dex                (float *x,int m){
        int ncoeff= 90;
        float avdat= -0.8897337E-03;
        float xmin[10]={
            -0.14981E-01,-0.47334E-01,-0.14977E-01,-0.32322E-01,-0.49995E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.47289E-01, 0.14991E-01, 0.30724E-01, 0.49882E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 91]={
             0.11231305E-03, 0.68509136E-03,-0.59476192E-02,-0.25091174E-02,
            -0.38171138E-02,-0.15575382E-05, 0.48435903E-02,-0.98805199E-03,
            -0.11172793E-01,-0.12040918E-02,-0.18402211E-02, 0.17040236E-02,
             0.14690236E-02, 0.13291125E-01, 0.18438930E-02,-0.20657941E-02,
            -0.14803694E-02, 0.19070592E-03, 0.16836372E-02,-0.85477103E-02,
            -0.39125192E-02, 0.29896488E-02,-0.17924338E-02,-0.17070892E-02,
             0.63493659E-04,-0.13212262E-02,-0.14069373E-02,-0.18551646E-03,
            -0.14113373E-02, 0.23576501E-03,-0.51887514E-03,-0.92558918E-03,
            -0.18842502E-02, 0.73171820E-03, 0.20104331E-03,-0.16166644E-03,
             0.43614597E-04,-0.34591509E-03,-0.38238979E-03,-0.51159674E-03,
            -0.20906022E-02, 0.28871052E-03,-0.33255131E-03, 0.70819689E-03,
            -0.35463701E-03, 0.20284296E-02,-0.80589799E-03,-0.10156983E-02,
            -0.35732685E-03, 0.13093621E-03, 0.17842466E-03, 0.11463617E-02,
            -0.38691139E-03, 0.14960459E-03,-0.61581651E-03,-0.26508197E-03,
             0.21604404E-02, 0.11759177E-03,-0.41181652E-03,-0.29778879E-03,
             0.71824608E-04, 0.25251810E-03,-0.32848558E-04, 0.15564926E-03,
            -0.21585594E-03, 0.45124507E-04,-0.52808813E-04,-0.18031189E-03,
             0.16449299E-03, 0.61579503E-03,-0.19519034E-03, 0.34908371E-03,
            -0.14137071E-03,-0.14430455E-03,-0.14924753E-03,-0.65441767E-04,
            -0.17844740E-03,-0.49400097E-03, 0.21749252E-03,-0.51897921E-03,
             0.60110306E-03, 0.36215762E-03,-0.10621082E-03,-0.39507510E-03,
            -0.13113902E-04, 0.29068082E-03, 0.20350425E-02, 0.51077816E-03,
             0.45885047E-03,-0.11198429E-02,
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
        float x45 = x44*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_p_l5p65_ep24_dex                =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]*x11                
            +coeff[  6]    *x22            
            +coeff[  7]    *x21*x31        
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[  8]    *x21    *x41    
            +coeff[  9]        *x31*x41    
            +coeff[ 10]            *x42    
            +coeff[ 11]    *x21        *x51
            +coeff[ 12]        *x31    *x51
            +coeff[ 13]            *x41*x51
            +coeff[ 14]*x11*x21            
            +coeff[ 15]    *x23            
            +coeff[ 16]                *x52
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 17]*x11    *x31        
            +coeff[ 18]*x11        *x41    
            +coeff[ 19]    *x22    *x41    
            +coeff[ 20]    *x21    *x42    
            +coeff[ 21]    *x22        *x51
            +coeff[ 22]    *x21    *x41*x51
            +coeff[ 23]    *x22    *x41*x51
            +coeff[ 24]    *x24*x31        
            +coeff[ 25]    *x22*x31        
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 26]    *x21*x31*x41    
            +coeff[ 27]*x11            *x51
            +coeff[ 28]            *x43    
            +coeff[ 29]*x12                
            +coeff[ 30]            *x41*x52
            +coeff[ 31]*x11        *x42    
            +coeff[ 32]    *x22    *x42    
            +coeff[ 33]*x11*x21        *x51
            +coeff[ 34]    *x26            
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 35]        *x32        
            +coeff[ 36]        *x32*x41    
            +coeff[ 37]        *x31*x42    
            +coeff[ 38]        *x31*x41*x51
            +coeff[ 39]*x11*x22            
            +coeff[ 40]    *x24            
            +coeff[ 41]    *x21        *x52
            +coeff[ 42]*x11*x21*x31        
            +coeff[ 43]    *x23    *x41    
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 44]*x11    *x31*x41    
            +coeff[ 45]    *x21    *x43    
            +coeff[ 46]*x11*x23            
            +coeff[ 47]*x11*x22    *x41    
            +coeff[ 48]    *x25        *x51
            +coeff[ 49]    *x21*x31    *x51
            +coeff[ 50]                *x53
            +coeff[ 51]    *x21*x31*x42    
            +coeff[ 52]    *x22*x31    *x51
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 53]*x11        *x41*x51
            +coeff[ 54]    *x21    *x42*x51
            +coeff[ 55]*x12        *x41    
            +coeff[ 56]    *x23    *x42    
            +coeff[ 57]*x12            *x51
            +coeff[ 58]*x11*x21    *x41*x51
            +coeff[ 59]        *x32*x42*x51
            +coeff[ 60]*x11*x23*x31    *x51
            +coeff[ 61]    *x23*x31        
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 62]*x11    *x32        
            +coeff[ 63]    *x21*x32*x41    
            +coeff[ 64]    *x21*x31*x41*x51
            +coeff[ 65]*x12*x21            
            +coeff[ 66]*x12    *x31        
            +coeff[ 67]    *x21    *x41*x52
            +coeff[ 68]*x11*x21*x31*x41    
            +coeff[ 69]    *x23*x31*x41    
            +coeff[ 70]*x11*x21*x31    *x51
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 71]    *x23    *x41*x51
            +coeff[ 72]            *x41*x53
            +coeff[ 73]    *x22    *x42*x51
            +coeff[ 74]*x12*x22            
            +coeff[ 75]*x11*x21        *x52
            +coeff[ 76]*x12*x21    *x41    
            +coeff[ 77]*x11*x23    *x41    
            +coeff[ 78]    *x22    *x41*x52
            +coeff[ 79]            *x44*x51
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 80]*x11*x22    *x42    
            +coeff[ 81]*x11*x21    *x43    
            +coeff[ 82]        *x32    *x53
            +coeff[ 83]    *x26    *x41    
            +coeff[ 84]        *x31*x45    
            +coeff[ 85]            *x45*x51
            +coeff[ 86]    *x22    *x44*x51
            +coeff[ 87]    *x22*x31*x41*x53
            +coeff[ 88]        *x32*x42*x53
        ;
        v_p_l5p65_ep24_dex                =v_p_l5p65_ep24_dex                
            +coeff[ 89]    *x24    *x45    
            ;
     
        return v_p_l5p65_ep24_dex                ;
    }
    float l_l5p65_ep24_dex                (float *x,int m){
        int ncoeff= 61;
        float avdat= -0.5696807E-02;
        float xmin[10]={
            -0.14981E-01,-0.47334E-01,-0.14977E-01,-0.32322E-01,-0.49995E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.47289E-01, 0.14991E-01, 0.30724E-01, 0.49882E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 62]={
             0.59275101E-02,-0.38796178E+00,-0.20619483E-03,-0.98868333E-01,
             0.65449528E-01,-0.22811055E-01,-0.15987945E+00,-0.34879036E-01,
            -0.34478791E-01, 0.18162131E-01,-0.53052314E-01,-0.11656574E-01,
            -0.41221160E-01, 0.68984099E-01,-0.26425026E-01,-0.21381620E-02,
            -0.29404552E-03, 0.56738281E-02,-0.31320408E-02,-0.26032941E-01,
             0.41549191E-01, 0.62743455E-01, 0.38219154E-01,-0.56030126E-02,
             0.17765535E-02,-0.73334640E-02, 0.66950074E-02, 0.10460925E-01,
            -0.65696491E-02, 0.11477116E-01,-0.12763208E-03, 0.12311048E-02,
            -0.30149966E-02, 0.15125558E-02,-0.44812071E-02, 0.70301355E-02,
            -0.20615220E-01, 0.74524386E-02, 0.37818265E-03,-0.15775076E-02,
            -0.14484866E-02,-0.60641591E-03,-0.89609367E-03, 0.17074113E-02,
             0.11304225E-02, 0.47889748E-02,-0.18991793E-01, 0.25440766E-02,
             0.35093699E-02,-0.10068533E-02, 0.28527816E-03,-0.45574998E-03,
            -0.18006003E-02,-0.40565457E-03,-0.37496090E-02, 0.28358272E-02,
            -0.54745222E-02,-0.12239910E-02,-0.28997331E-02, 0.51390459E-02,
             0.87077478E-02,
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
     
        float v_l_l5p65_ep24_dex                =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]            *x41    
            +coeff[  3]                *x51
            +coeff[  4]*x11                
            +coeff[  5]    *x22            
            +coeff[  6]    *x21    *x41    
            +coeff[  7]    *x21        *x51
        ;
        v_l_l5p65_ep24_dex                =v_l_l5p65_ep24_dex                
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
        v_l_l5p65_ep24_dex                =v_l_l5p65_ep24_dex                
            +coeff[ 17]*x11*x21            
            +coeff[ 18]*x11            *x51
            +coeff[ 19]    *x22    *x41    
            +coeff[ 20]    *x21*x31*x41    
            +coeff[ 21]    *x23    *x41    
            +coeff[ 22]*x11*x22    *x41    
            +coeff[ 23]            *x42    
            +coeff[ 24]                *x52
            +coeff[ 25]    *x22*x31        
        ;
        v_l_l5p65_ep24_dex                =v_l_l5p65_ep24_dex                
            +coeff[ 26]    *x21*x32        
            +coeff[ 27]*x11        *x42    
            +coeff[ 28]*x12*x21            
            +coeff[ 29]*x11*x22*x31        
            +coeff[ 30]*x11    *x33*x41    
            +coeff[ 31]            *x41*x51
            +coeff[ 32]    *x21    *x41*x51
            +coeff[ 33]    *x21        *x52
            +coeff[ 34]*x11*x21    *x41    
        ;
        v_l_l5p65_ep24_dex                =v_l_l5p65_ep24_dex                
            +coeff[ 35]*x11    *x31*x41    
            +coeff[ 36]    *x21    *x43    
            +coeff[ 37]*x12*x21    *x41    
            +coeff[ 38]*x11*x21*x33        
            +coeff[ 39]    *x23*x31*x42    
            +coeff[ 40]        *x31*x41    
            +coeff[ 41]*x12                
            +coeff[ 42]    *x22        *x51
            +coeff[ 43]            *x42*x51
        ;
        v_l_l5p65_ep24_dex                =v_l_l5p65_ep24_dex                
            +coeff[ 44]*x11    *x32        
            +coeff[ 45]    *x22    *x42    
            +coeff[ 46]    *x21*x31*x42    
            +coeff[ 47]*x12*x21*x31        
            +coeff[ 48]    *x21*x31*x43    
            +coeff[ 49]    *x23*x32*x41    
            +coeff[ 50]        *x31    *x51
            +coeff[ 51]            *x41*x52
            +coeff[ 52]*x11*x21*x31        
        ;
        v_l_l5p65_ep24_dex                =v_l_l5p65_ep24_dex                
            +coeff[ 53]*x13                
            +coeff[ 54]    *x24            
            +coeff[ 55]    *x22*x31*x41    
            +coeff[ 56]    *x21*x32*x41    
            +coeff[ 57]    *x23        *x51
            +coeff[ 58]*x11*x23            
            +coeff[ 59]    *x24    *x41    
            +coeff[ 60]    *x21    *x44    
            ;
     
        return v_l_l5p65_ep24_dex                ;
    }
    float x_l5p65_ep26_q3en               (float *x,int m){
        int ncoeff= 58;
        float avdat= -0.5767067E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 59]={
            -0.80824988E-02,-0.34104515E-01, 0.13752362E+00, 0.16842097E+00,
             0.27105821E-01, 0.13120853E-01, 0.69171317E-01, 0.95135150E-02,
            -0.15109820E-03,-0.93937060E-03,-0.77534332E-02,-0.59652450E-02,
             0.54480694E-02, 0.15343712E-01, 0.24070609E-01, 0.18052470E-01,
            -0.32125220E-01, 0.24638313E-02,-0.30393286E-02,-0.41416618E-02,
             0.11351391E-01, 0.37076098E-02,-0.18219987E-01,-0.26974753E-01,
             0.29340603E-02, 0.50748261E-02,-0.30066385E-02,-0.16866289E-01,
             0.77406591E-03, 0.13137716E-02,-0.67733455E-03,-0.42860624E-02,
             0.18016354E-02, 0.12714298E-02,-0.36280772E-02,-0.51749516E-02,
            -0.35861670E-03, 0.44157749E-03, 0.84803608E-03,-0.52686717E-03,
             0.18129948E-02,-0.29289285E-02,-0.82118780E-03, 0.66443416E-02,
            -0.90404856E-03,-0.82331459E-03, 0.19379627E-03, 0.75684155E-04,
             0.29752494E-03,-0.40286893E-03, 0.25061698E-03, 0.86957595E-03,
             0.16783761E-02,-0.70126652E-03, 0.22402096E-02, 0.73119984E-02,
             0.57043979E-03, 0.97908871E-03,
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
     
        float v_x_l5p65_ep26_q3en               =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]                *x51
            +coeff[  3]    *x21            
            +coeff[  4]    *x21        *x51
            +coeff[  5]    *x22            
            +coeff[  6]    *x21    *x41    
            +coeff[  7]    *x22    *x41    
        ;
        v_x_l5p65_ep26_q3en               =v_x_l5p65_ep26_q3en               
            +coeff[  8]*x13        *x41    
            +coeff[  9]*x12*x21*x31        
            +coeff[ 10]    *x23*x31        
            +coeff[ 11]                *x52
            +coeff[ 12]*x11    *x31        
            +coeff[ 13]*x11        *x41    
            +coeff[ 14]    *x21*x31        
            +coeff[ 15]    *x23            
            +coeff[ 16]    *x21    *x42    
        ;
        v_x_l5p65_ep26_q3en               =v_x_l5p65_ep26_q3en               
            +coeff[ 17]            *x41    
            +coeff[ 18]*x11*x21            
            +coeff[ 19]            *x42    
            +coeff[ 20]*x11*x22            
            +coeff[ 21]    *x22*x31        
            +coeff[ 22]    *x21*x31*x41    
            +coeff[ 23]    *x23    *x41    
            +coeff[ 24]*x12*x21            
            +coeff[ 25]    *x21    *x41*x51
        ;
        v_x_l5p65_ep26_q3en               =v_x_l5p65_ep26_q3en               
            +coeff[ 26]    *x21*x32        
            +coeff[ 27]*x11*x22    *x41    
            +coeff[ 28]        *x31        
            +coeff[ 29]            *x41*x51
            +coeff[ 30]        *x31*x41    
            +coeff[ 31]*x11        *x42    
            +coeff[ 32]    *x22        *x51
            +coeff[ 33]    *x21*x31    *x51
            +coeff[ 34]*x12*x21    *x41    
        ;
        v_x_l5p65_ep26_q3en               =v_x_l5p65_ep26_q3en               
            +coeff[ 35]*x11*x22*x31        
            +coeff[ 36]*x13    *x31*x41    
            +coeff[ 37]*x11        *x41*x51
            +coeff[ 38]*x11*x21*x31        
            +coeff[ 39]*x11    *x32        
            +coeff[ 40]*x11*x21    *x41    
            +coeff[ 41]*x11    *x31*x41    
            +coeff[ 42]            *x42*x51
            +coeff[ 43]    *x21    *x43    
        ;
        v_x_l5p65_ep26_q3en               =v_x_l5p65_ep26_q3en               
            +coeff[ 44]*x12*x21*x31*x42    
            +coeff[ 45]    *x23*x31*x42    
            +coeff[ 46]*x12                
            +coeff[ 47]*x11            *x51
            +coeff[ 48]                *x53
            +coeff[ 49]    *x21        *x52
            +coeff[ 50]            *x41*x52
            +coeff[ 51]*x11*x23            
            +coeff[ 52]    *x23        *x51
        ;
        v_x_l5p65_ep26_q3en               =v_x_l5p65_ep26_q3en               
            +coeff[ 53]    *x22*x31*x41    
            +coeff[ 54]    *x21*x32*x41    
            +coeff[ 55]    *x21*x31*x42    
            +coeff[ 56]*x11*x22*x31    *x51
            +coeff[ 57]*x13*x22        *x51
            ;
     
        return v_x_l5p65_ep26_q3en               ;
    }
    float t_l5p65_ep26_q3en               (float *x,int m){
        int ncoeff= 80;
        float avdat=  0.3015241E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 81]={
             0.19268958E-02, 0.85407514E-02,-0.62420078E-01, 0.22689067E-01,
            -0.35066006E-02,-0.75292587E-02,-0.21447213E-01,-0.38962387E-02,
             0.85195024E-04,-0.75047236E-03,-0.16941904E-02,-0.47006966E-02,
            -0.14265346E-02,-0.55252011E-02, 0.84681371E-02,-0.46740723E-03,
            -0.24247401E-03, 0.11214594E-02,-0.58396871E-03,-0.36030393E-02,
            -0.33107395E-02, 0.57116337E-02, 0.78564864E-02,-0.19977169E-04,
            -0.92607894E-03,-0.12592026E-02, 0.92941424E-03, 0.18281513E-02,
             0.24890967E-02, 0.50504818E-02, 0.14271392E-05,-0.23159230E-03,
            -0.50708366E-03, 0.44913933E-03, 0.35865980E-03, 0.16535006E-02,
             0.11127815E-02,-0.25971823E-02,-0.10298346E-03,-0.18181239E-04,
             0.42156094E-04,-0.19524198E-03,-0.23030835E-03, 0.21888383E-03,
             0.12512943E-02,-0.14718472E-03,-0.18470912E-03, 0.32610400E-03,
             0.92848466E-03,-0.25315434E-02,-0.45693276E-03,-0.18125944E-03,
            -0.14516676E-03, 0.23410194E-03, 0.13475603E-03, 0.36681831E-03,
            -0.79775177E-03,-0.57765207E-03, 0.41595654E-03,-0.84216363E-03,
             0.11067509E-03, 0.91544702E-04,-0.67851848E-04, 0.16564094E-03,
             0.14373004E-03, 0.45484656E-04, 0.66189896E-04, 0.27510221E-04,
             0.73097282E-04,-0.23285423E-03,-0.13670874E-03,-0.63741801E-03,
             0.51528197E-04,-0.11225511E-03, 0.12338971E-03, 0.11828647E-03,
            -0.96341304E-04, 0.93660725E-04,-0.51431771E-03, 0.20647590E-03,
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
     
        float v_t_l5p65_ep26_q3en               =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]                *x51
            +coeff[  4]    *x22            
            +coeff[  5]    *x21*x31        
            +coeff[  6]    *x21    *x41    
            +coeff[  7]    *x21        *x51
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[  8]*x13        *x41    
            +coeff[  9]            *x41    
            +coeff[ 10]*x11    *x31        
            +coeff[ 11]*x11        *x41    
            +coeff[ 12]                *x52
            +coeff[ 13]    *x23            
            +coeff[ 14]    *x21    *x42    
            +coeff[ 15]    *x23*x31*x41    
            +coeff[ 16]        *x31        
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[ 17]*x11*x21            
            +coeff[ 18]*x11            *x51
            +coeff[ 19]*x11*x22            
            +coeff[ 20]    *x22    *x41    
            +coeff[ 21]    *x21*x31*x41    
            +coeff[ 22]    *x23    *x41    
            +coeff[ 23]*x13*x22    *x41    
            +coeff[ 24]*x12*x21            
            +coeff[ 25]    *x22*x31        
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[ 26]    *x21*x32        
            +coeff[ 27]*x11        *x42    
            +coeff[ 28]    *x23*x31        
            +coeff[ 29]*x11*x22    *x41    
            +coeff[ 30]*x13    *x31*x41    
            +coeff[ 31]        *x31*x41    
            +coeff[ 32]*x11*x21    *x41    
            +coeff[ 33]            *x42*x51
            +coeff[ 34]    *x21        *x52
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[ 35]*x11*x22*x31        
            +coeff[ 36]*x12*x21    *x41    
            +coeff[ 37]    *x21    *x43    
            +coeff[ 38]*x12                
            +coeff[ 39]            *x42    
            +coeff[ 40]        *x31    *x51
            +coeff[ 41]            *x41*x51
            +coeff[ 42]*x11*x21*x31        
            +coeff[ 43]*x11    *x32        
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[ 44]*x11    *x31*x41    
            +coeff[ 45]*x11*x21        *x51
            +coeff[ 46]            *x41*x52
            +coeff[ 47]*x12*x21*x31        
            +coeff[ 48]    *x22    *x42    
            +coeff[ 49]    *x21*x31*x42    
            +coeff[ 50]    *x21    *x42*x51
            +coeff[ 51]*x12*x21*x32*x41    
            +coeff[ 52]    *x22        *x51
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[ 53]    *x21    *x41*x51
            +coeff[ 54]                *x53
            +coeff[ 55]    *x22*x31*x41    
            +coeff[ 56]    *x21*x32*x41    
            +coeff[ 57]*x11        *x43    
            +coeff[ 58]    *x22    *x41*x51
            +coeff[ 59]*x11*x22    *x42    
            +coeff[ 60]    *x21*x31    *x53
            +coeff[ 61]*x13    *x31*x42    
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[ 62]*x13                
            +coeff[ 63]        *x31*x42    
            +coeff[ 64]            *x43    
            +coeff[ 65]*x11    *x31    *x51
            +coeff[ 66]*x11        *x41*x51
            +coeff[ 67]*x13    *x31        
            +coeff[ 68]    *x22*x32        
            +coeff[ 69]*x11    *x32*x41    
            +coeff[ 70]*x11*x21    *x42    
        ;
        v_t_l5p65_ep26_q3en               =v_t_l5p65_ep26_q3en               
            +coeff[ 71]*x11    *x31*x42    
            +coeff[ 72]*x11*x22        *x51
            +coeff[ 73]    *x23        *x51
            +coeff[ 74]    *x22*x31    *x51
            +coeff[ 75]*x11*x21    *x41*x51
            +coeff[ 76]    *x21*x31*x41*x51
            +coeff[ 77]            *x42*x52
            +coeff[ 78]*x11*x22*x31*x41    
            +coeff[ 79]    *x22*x32*x41    
            ;
     
        return v_t_l5p65_ep26_q3en               ;
    }
    float y_l5p65_ep26_q3en               (float *x,int m){
        int ncoeff= 97;
        float avdat=  0.1141901E-03;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 98]={
            -0.49567721E-02, 0.83266646E-01, 0.16425070E-02,-0.21110952E-01,
             0.37643535E-03,-0.30789057E-01,-0.56865130E-01,-0.11591160E-01,
             0.42084839E-01,-0.51598325E-02, 0.84650908E-02, 0.70466258E-01,
             0.17764362E-01, 0.83500640E-02, 0.95674004E-02,-0.17322991E-01,
            -0.68025859E-02,-0.50356284E-01,-0.68558576E-02,-0.90834871E-02,
            -0.11643404E-01, 0.39574763E-03, 0.13281249E-01,-0.16716108E-01,
            -0.19077807E-02,-0.69103455E-02, 0.20635317E-02,-0.62454394E-02,
            -0.45523443E-02,-0.36628712E-02,-0.58245156E-02, 0.35007806E-02,
            -0.11269955E-02,-0.52183545E-02,-0.10965826E-02,-0.18177869E-02,
            -0.11581951E-01,-0.41613802E-02,-0.49647335E-02,-0.68710693E-02,
            -0.45467788E-03, 0.87899889E-03,-0.12094872E-02,-0.22412420E-02,
            -0.20832550E-02, 0.11094757E-01, 0.61413813E-02,-0.81520522E-03,
             0.10329125E-02,-0.23201199E-02,-0.26047771E-03,-0.15432430E-03,
             0.45471545E-03,-0.54858200E-03, 0.62142976E-03, 0.46214475E-02,
             0.39116149E-02,-0.25725639E-02, 0.10894862E-02, 0.16191860E-02,
             0.62799180E-03, 0.37211038E-02,-0.27235844E-02, 0.86443470E-03,
            -0.15313213E-02,-0.11353639E-02, 0.76314751E-02,-0.16535533E-02,
            -0.83042122E-03, 0.46186437E-03, 0.35824657E-02, 0.37596147E-02,
             0.52593631E-03, 0.10824084E-02,-0.14725466E-02, 0.15534010E-02,
            -0.33888154E-03,-0.20089111E-03,-0.21510715E-03, 0.76989020E-03,
             0.33798299E-03,-0.22740418E-02, 0.10865346E-02, 0.27895984E-02,
             0.11448434E-02, 0.86560939E-03,-0.74157672E-03,-0.25153060E-02,
            -0.19985363E-02, 0.61145762E-03, 0.57263568E-03, 0.54061325E-03,
            -0.10352780E-02,-0.27514305E-02, 0.19406661E-02,-0.40960283E-03,
            -0.46901969E-03,
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
     
        float v_y_l5p65_ep26_q3en               =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]        *x31        
            +coeff[  4]*x11                
            +coeff[  5]                *x51
            +coeff[  6]    *x21    *x41    
            +coeff[  7]        *x31*x41    
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
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
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]    *x21*x31*x41    
            +coeff[ 19]    *x23            
            +coeff[ 20]    *x22*x31        
            +coeff[ 21]            *x44    
            +coeff[ 22]    *x22        *x51
            +coeff[ 23]            *x42    
            +coeff[ 24]        *x32        
            +coeff[ 25]            *x43    
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 26]*x12                
            +coeff[ 27]*x11*x21    *x41    
            +coeff[ 28]    *x21    *x41*x51
            +coeff[ 29]*x11*x21*x31        
            +coeff[ 30]            *x41*x52
            +coeff[ 31]*x11*x21        *x51
            +coeff[ 32]*x11            *x51
            +coeff[ 33]*x11        *x42    
            +coeff[ 34]        *x31*x41*x51
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 35]*x12        *x41    
            +coeff[ 36]    *x24            
            +coeff[ 37]*x11*x22    *x41    
            +coeff[ 38]    *x22    *x41*x51
            +coeff[ 39]*x11*x23            
            +coeff[ 40]    *x22*x31*x43    
            +coeff[ 41]*x11    *x31        
            +coeff[ 42]        *x31*x42    
            +coeff[ 43]*x11    *x31*x41    
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 44]*x11*x22            
            +coeff[ 45]    *x21    *x43    
            +coeff[ 46]    *x21*x31*x42    
            +coeff[ 47]        *x31    *x52
            +coeff[ 48]                *x53
            +coeff[ 49]*x11*x21    *x41*x51
            +coeff[ 50]*x12        *x42*x51
            +coeff[ 51]*x12*x22    *x42    
            +coeff[ 52]*x12*x21            
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 53]*x12    *x31        
            +coeff[ 54]    *x21        *x52
            +coeff[ 55]    *x22*x31*x41    
            +coeff[ 56]*x11*x21    *x42    
            +coeff[ 57]    *x21    *x42*x51
            +coeff[ 58]    *x23*x31        
            +coeff[ 59]    *x22*x32        
            +coeff[ 60]*x12            *x51
            +coeff[ 61]*x11*x21*x31*x41    
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 62]    *x21    *x44    
            +coeff[ 63]*x12        *x42    
            +coeff[ 64]    *x22*x31    *x51
            +coeff[ 65]*x12*x21    *x41    
            +coeff[ 66]    *x23    *x42    
            +coeff[ 67]*x12*x22            
            +coeff[ 68]*x11*x21*x31    *x51
            +coeff[ 69]        *x34*x41    
            +coeff[ 70]*x11*x22    *x42    
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 71]    *x22    *x42*x51
            +coeff[ 72]    *x23*x31    *x51
            +coeff[ 73]            *x45*x51
            +coeff[ 74]    *x23    *x43    
            +coeff[ 75]*x11*x23*x32        
            +coeff[ 76]    *x21*x32        
            +coeff[ 77]*x11    *x32        
            +coeff[ 78]        *x32    *x51
            +coeff[ 79]        *x31*x43    
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 80]*x11        *x41*x51
            +coeff[ 81]    *x22    *x42    
            +coeff[ 82]*x11        *x43    
            +coeff[ 83]    *x23    *x41    
            +coeff[ 84]    *x21*x32*x41    
            +coeff[ 85]*x11    *x31*x42    
            +coeff[ 86]    *x21*x31*x41*x51
            +coeff[ 87]    *x22    *x43    
            +coeff[ 88]    *x21*x31*x43    
        ;
        v_y_l5p65_ep26_q3en               =v_y_l5p65_ep26_q3en               
            +coeff[ 89]*x12    *x31*x41    
            +coeff[ 90]    *x21    *x41*x52
            +coeff[ 91]*x11*x22        *x51
            +coeff[ 92]    *x22        *x52
            +coeff[ 93]    *x24    *x41    
            +coeff[ 94]    *x23*x31*x41    
            +coeff[ 95]*x12        *x41*x51
            +coeff[ 96]*x11*x21        *x52
            ;
     
        return v_y_l5p65_ep26_q3en               ;
    }
    float p_l5p65_ep26_q3en               (float *x,int m){
        int ncoeff= 90;
        float avdat= -0.1075884E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 91]={
             0.13994661E-03, 0.20466933E-03,-0.59229331E-02, 0.20426896E-02,
            -0.42572268E-02, 0.44815155E-04, 0.56037996E-02,-0.13209431E-02,
            -0.11056893E-01,-0.14675951E-02,-0.21374302E-02, 0.14072302E-02,
             0.16174641E-02, 0.14067519E-01, 0.23182342E-02,-0.16572677E-02,
            -0.16721940E-02, 0.16264531E-02,-0.96236244E-02,-0.35269721E-02,
             0.31278469E-02,-0.12330235E-02, 0.78937213E-03,-0.15326238E-02,
             0.84561971E-05,-0.19832645E-04,-0.16581308E-02,-0.15705326E-03,
            -0.14477075E-02, 0.29379813E-03,-0.21220413E-02,-0.44470571E-03,
            -0.21780776E-02,-0.26925580E-04,-0.25728601E-03, 0.25484696E-03,
            -0.15642849E-02,-0.28701025E-03,-0.36522013E-03,-0.43842717E-03,
            -0.92896010E-03, 0.17359649E-02,-0.10281764E-02,-0.29592996E-03,
            -0.88282430E-03, 0.39833784E-03,-0.30565214E-04,-0.73486197E-04,
             0.10290484E-03, 0.18706574E-03, 0.42704778E-04,-0.52570464E-03,
            -0.42989099E-03, 0.10807751E-02,-0.51738520E-03, 0.15579186E-02,
            -0.47746632E-03, 0.55400183E-03,-0.17903489E-03, 0.46159895E-03,
             0.78660867E-03, 0.21723714E-04,-0.10617686E-03,-0.73553274E-04,
            -0.24237241E-03, 0.45749257E-03, 0.22705493E-03, 0.17152219E-03,
             0.22681252E-03,-0.31123919E-03, 0.50225819E-04,-0.55343677E-04,
             0.36391645E-03, 0.63364342E-03, 0.22394997E-03, 0.17117892E-03,
             0.10308740E-03,-0.19791175E-03, 0.73519000E-03,-0.20712242E-03,
            -0.61742426E-03, 0.45938065E-03,-0.55005378E-03,-0.19339631E-03,
             0.22974331E-03, 0.15847851E-03, 0.63132953E-04, 0.11164095E-03,
            -0.32564156E-04,-0.28375278E-04,
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
     
        float v_p_l5p65_ep26_q3en               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]*x11                
            +coeff[  6]    *x22            
            +coeff[  7]    *x21*x31        
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[  8]    *x21    *x41    
            +coeff[  9]        *x31*x41    
            +coeff[ 10]            *x42    
            +coeff[ 11]    *x21        *x51
            +coeff[ 12]        *x31    *x51
            +coeff[ 13]            *x41*x51
            +coeff[ 14]*x11*x21            
            +coeff[ 15]    *x23            
            +coeff[ 16]                *x52
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 17]*x11        *x41    
            +coeff[ 18]    *x22    *x41    
            +coeff[ 19]    *x21    *x42    
            +coeff[ 20]    *x22        *x51
            +coeff[ 21]    *x21    *x41*x51
            +coeff[ 22]*x11*x21        *x51
            +coeff[ 23]    *x22    *x41*x51
            +coeff[ 24]    *x24*x31        
            +coeff[ 25]    *x21*x33*x41    
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 26]    *x22*x31        
            +coeff[ 27]*x11            *x51
            +coeff[ 28]            *x43    
            +coeff[ 29]*x12                
            +coeff[ 30]    *x24            
            +coeff[ 31]*x11*x21*x31        
            +coeff[ 32]    *x22    *x42    
            +coeff[ 33]            *x41*x54
            +coeff[ 34]        *x32        
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 35]*x11    *x31        
            +coeff[ 36]    *x21*x31*x41    
            +coeff[ 37]        *x31*x42    
            +coeff[ 38]        *x31*x41*x51
            +coeff[ 39]*x11*x22            
            +coeff[ 40]*x11        *x42    
            +coeff[ 41]    *x21    *x43    
            +coeff[ 42]*x11*x23            
            +coeff[ 43]*x12        *x41    
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 44]*x11*x22    *x41    
            +coeff[ 45]*x11*x22*x31*x41    
            +coeff[ 46]*x11    *x33*x41    
            +coeff[ 47]*x12            *x53
            +coeff[ 48]    *x21*x31    *x51
            +coeff[ 49]    *x21        *x52
            +coeff[ 50]*x11*x21    *x41    
            +coeff[ 51]            *x41*x52
            +coeff[ 52]*x11    *x31*x41    
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 53]    *x21*x31*x42    
            +coeff[ 54]    *x22*x31    *x51
            +coeff[ 55]    *x23    *x42    
            +coeff[ 56]*x11*x21    *x41*x51
            +coeff[ 57]    *x23    *x41*x51
            +coeff[ 58]*x12*x21    *x41    
            +coeff[ 59]    *x22    *x41*x52
            +coeff[ 60]*x11*x22    *x42    
            +coeff[ 61]*x11*x21*x33    *x51
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 62]    *x21*x32        
            +coeff[ 63]        *x32    *x51
            +coeff[ 64]            *x42*x51
            +coeff[ 65]    *x23    *x41    
            +coeff[ 66]    *x22*x32        
            +coeff[ 67]                *x53
            +coeff[ 68]    *x21*x32*x41    
            +coeff[ 69]    *x21    *x42*x51
            +coeff[ 70]*x12*x21            
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 71]*x12    *x31        
            +coeff[ 72]*x11*x21*x31*x41    
            +coeff[ 73]    *x23*x31*x41    
            +coeff[ 74]*x11*x21    *x42    
            +coeff[ 75]*x12            *x51
            +coeff[ 76]*x11*x22        *x51
            +coeff[ 77]*x11*x21*x31    *x51
            +coeff[ 78]    *x22    *x42*x51
            +coeff[ 79]*x12*x22            
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 80]*x11*x23    *x41    
            +coeff[ 81]    *x24    *x42    
            +coeff[ 82]    *x26    *x41    
            +coeff[ 83]    *x21*x33*x41*x51
            +coeff[ 84]*x11*x23*x32        
            +coeff[ 85]            *x45*x51
            +coeff[ 86]        *x32*x41    
            +coeff[ 87]    *x23*x31        
            +coeff[ 88]        *x31    *x52
        ;
        v_p_l5p65_ep26_q3en               =v_p_l5p65_ep26_q3en               
            +coeff[ 89]*x11    *x32        
            ;
     
        return v_p_l5p65_ep26_q3en               ;
    }
    float l_l5p65_ep26_q3en               (float *x,int m){
        int ncoeff= 68;
        float avdat= -0.8692704E-02;
        float xmin[10]={
            -0.14979E-01,-0.47010E-01,-0.14993E-01,-0.29713E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14992E-01, 0.30136E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 69]={
             0.12477567E-01,-0.24901684E+00,-0.32336015E-01, 0.39415825E-01,
            -0.19154182E-01,-0.31808693E-01,-0.91185354E-01,-0.15452548E-01,
            -0.20253446E-01, 0.57267752E-02,-0.72334446E-02,-0.23117704E-01,
             0.35720646E-01,-0.16038232E-01,-0.11501692E-02,-0.85428450E-02,
            -0.28487353E-02,-0.21698195E-01, 0.22962503E-01, 0.13172952E-02,
            -0.38475720E-02, 0.34228958E-01, 0.21785924E-01, 0.30074050E-02,
             0.16999217E-02,-0.60459818E-02, 0.39596376E-02,-0.39799078E-02,
             0.60709678E-02, 0.10218356E-01, 0.69474559E-02, 0.24267808E-03,
             0.82591502E-03,-0.25032135E-02, 0.82944211E-03, 0.40899678E-02,
             0.45012161E-02, 0.20124756E-02, 0.24669548E-03,-0.52099384E-03,
            -0.95823372E-03,-0.45052732E-03,-0.10766366E-02, 0.73520618E-03,
            -0.33566519E-02, 0.51384275E-02,-0.10937200E-01,-0.10316891E-01,
            -0.25769733E-02, 0.15392984E-02, 0.55746292E-02,-0.76818615E-04,
            -0.51762111E-03,-0.56795012E-04,-0.27774426E-03, 0.10814241E-02,
             0.70283725E-03, 0.54015632E-03, 0.47194876E-03, 0.39235508E-03,
            -0.30459822E-03,-0.16758383E-03, 0.29780851E-02,-0.32989057E-02,
             0.56046381E-03, 0.40899618E-02, 0.25835470E-02, 0.13654457E-02,
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
     
        float v_l_l5p65_ep26_q3en               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]                *x51
            +coeff[  3]*x11                
            +coeff[  4]    *x22            
            +coeff[  5]    *x21*x31        
            +coeff[  6]    *x21    *x41    
            +coeff[  7]    *x21        *x51
        ;
        v_l_l5p65_ep26_q3en               =v_l_l5p65_ep26_q3en               
            +coeff[  8]*x11        *x41    
            +coeff[  9]*x11*x21            
            +coeff[ 10]*x11    *x31        
            +coeff[ 11]    *x23            
            +coeff[ 12]    *x21    *x42    
            +coeff[ 13]*x11*x22            
            +coeff[ 14]    *x23*x31*x41    
            +coeff[ 15]            *x42    
            +coeff[ 16]*x11            *x51
        ;
        v_l_l5p65_ep26_q3en               =v_l_l5p65_ep26_q3en               
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]    *x21*x31*x41    
            +coeff[ 19]            *x43    
            +coeff[ 20]*x12*x21            
            +coeff[ 21]    *x23    *x41    
            +coeff[ 22]*x11*x22    *x41    
            +coeff[ 23]            *x41    
            +coeff[ 24]            *x41*x51
            +coeff[ 25]    *x22*x31        
        ;
        v_l_l5p65_ep26_q3en               =v_l_l5p65_ep26_q3en               
            +coeff[ 26]    *x21*x32        
            +coeff[ 27]*x11*x21    *x41    
            +coeff[ 28]*x11        *x42    
            +coeff[ 29]    *x23*x31        
            +coeff[ 30]*x11*x22*x31        
            +coeff[ 31]*x11    *x33*x41    
            +coeff[ 32]        *x31        
            +coeff[ 33]        *x31*x41    
            +coeff[ 34]    *x21        *x52
        ;
        v_l_l5p65_ep26_q3en               =v_l_l5p65_ep26_q3en               
            +coeff[ 35]*x11    *x31*x41    
            +coeff[ 36]*x12*x21    *x41    
            +coeff[ 37]            *x44*x51
            +coeff[ 38]        *x31    *x51
            +coeff[ 39]*x12                
            +coeff[ 40]    *x22        *x51
            +coeff[ 41]            *x41*x52
            +coeff[ 42]*x11*x21*x31        
            +coeff[ 43]*x11    *x32        
        ;
        v_l_l5p65_ep26_q3en               =v_l_l5p65_ep26_q3en               
            +coeff[ 44]    *x24            
            +coeff[ 45]    *x22    *x42    
            +coeff[ 46]    *x21*x31*x42    
            +coeff[ 47]    *x21    *x43    
            +coeff[ 48]*x11*x23            
            +coeff[ 49]*x12*x21*x31        
            +coeff[ 50]    *x21    *x44    
            +coeff[ 51]    *x24*x31*x41    
            +coeff[ 52]    *x23*x32*x41    
        ;
        v_l_l5p65_ep26_q3en               =v_l_l5p65_ep26_q3en               
            +coeff[ 53]    *x23*x33    *x51
            +coeff[ 54]                *x52
            +coeff[ 55]        *x31*x42    
            +coeff[ 56]    *x21*x31    *x51
            +coeff[ 57]    *x21    *x41*x51
            +coeff[ 58]        *x31*x41*x51
            +coeff[ 59]*x11        *x41*x51
            +coeff[ 60]*x13                
            +coeff[ 61]        *x34        
        ;
        v_l_l5p65_ep26_q3en               =v_l_l5p65_ep26_q3en               
            +coeff[ 62]    *x22*x31*x41    
            +coeff[ 63]    *x21*x32*x41    
            +coeff[ 64]    *x22    *x41*x51
            +coeff[ 65]    *x24    *x41    
            +coeff[ 66]    *x21*x31*x43    
            +coeff[ 67]*x11*x24            
            ;
     
        return v_l_l5p65_ep26_q3en               ;
    }
    float x_l5p65_ep29_q3ex               (float *x,int m){
        int ncoeff= 54;
        float avdat=  0.1438850E-01;
        float xmin[10]={
            -0.14979E-01,-0.46767E-01,-0.14993E-01,-0.29713E-01,-0.46582E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14981E-01, 0.30136E-01, 0.49667E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 55]={
            -0.89000650E-02,-0.27006982E-01, 0.32085121E+00, 0.53453211E-01,
            -0.21539886E-01, 0.28067745E-01, 0.38292062E-01, 0.98708561E-02,
             0.12929000E-01,-0.23646757E-01, 0.78949183E-02,-0.70239492E-02,
             0.10360048E-01, 0.15826626E-02, 0.29141074E-02,-0.19674606E-02,
             0.54770429E-02, 0.82065230E-02,-0.10153783E-01,-0.17655987E-01,
             0.50667080E-03,-0.13451432E-02, 0.17828146E-02, 0.24572255E-02,
             0.15112623E-02, 0.18055289E-02, 0.40886761E-02,-0.96193627E-02,
            -0.14849224E-02, 0.44672971E-03,-0.15351380E-02,-0.14000449E-02,
             0.17763120E-02, 0.15174528E-02,-0.16208162E-02, 0.25492853E-02,
             0.21565037E-02,-0.19477272E-02,-0.39918642E-02, 0.26079384E-02,
            -0.18216629E-02, 0.44520217E-03, 0.84514701E-03, 0.44708213E-03,
             0.77868352E-03,-0.72133297E-03, 0.56928518E-03, 0.10560332E-02,
            -0.14492052E-03, 0.17322799E-02,-0.27612306E-02, 0.30923099E-02,
             0.24175243E-02, 0.18735357E-02,
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
        float x53 = x52*x5;
     
    //                 function
     
        float v_x_l5p65_ep29_q3ex               =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]                *x51
            +coeff[  3]    *x21            
            +coeff[  4]                *x52
            +coeff[  5]    *x21        *x51
            +coeff[  6]    *x21    *x41    
            +coeff[  7]    *x22            
        ;
        v_x_l5p65_ep29_q3ex               =v_x_l5p65_ep29_q3ex               
            +coeff[  8]    *x21*x31        
            +coeff[  9]    *x21    *x42    
            +coeff[ 10]*x11        *x41    
            +coeff[ 11]            *x42    
            +coeff[ 12]    *x23            
            +coeff[ 13]            *x41    
            +coeff[ 14]*x11    *x31        
            +coeff[ 15]        *x31*x41    
            +coeff[ 16]*x11*x22            
        ;
        v_x_l5p65_ep29_q3ex               =v_x_l5p65_ep29_q3ex               
            +coeff[ 17]    *x21    *x41*x51
            +coeff[ 18]    *x21*x31*x41    
            +coeff[ 19]    *x23    *x41    
            +coeff[ 20]        *x31        
            +coeff[ 21]*x11            *x51
            +coeff[ 22]            *x41*x51
            +coeff[ 23]                *x53
            +coeff[ 24]*x12*x21            
            +coeff[ 25]    *x22        *x51
        ;
        v_x_l5p65_ep29_q3ex               =v_x_l5p65_ep29_q3ex               
            +coeff[ 26]    *x22    *x41    
            +coeff[ 27]*x11*x22    *x41    
            +coeff[ 28]*x11*x21            
            +coeff[ 29]        *x31    *x51
            +coeff[ 30]*x11    *x31*x41    
            +coeff[ 31]*x11        *x42    
            +coeff[ 32]    *x21*x31    *x51
            +coeff[ 33]    *x22*x31        
            +coeff[ 34]    *x21*x32        
        ;
        v_x_l5p65_ep29_q3ex               =v_x_l5p65_ep29_q3ex               
            +coeff[ 35]*x11*x23            
            +coeff[ 36]    *x23        *x51
            +coeff[ 37]    *x21    *x42*x51
            +coeff[ 38]    *x23*x31        
            +coeff[ 39]    *x22    *x42    
            +coeff[ 40]*x14*x21    *x41    
            +coeff[ 41]*x13*x22*x31        
            +coeff[ 42]*x13*x23    *x41    
            +coeff[ 43]*x11    *x31    *x51
        ;
        v_x_l5p65_ep29_q3ex               =v_x_l5p65_ep29_q3ex               
            +coeff[ 44]*x11        *x41*x51
            +coeff[ 45]    *x21        *x52
            +coeff[ 46]*x11*x21*x31        
            +coeff[ 47]*x11*x21    *x41    
            +coeff[ 48]*x13    *x31        
            +coeff[ 49]*x11*x22        *x51
            +coeff[ 50]*x11*x22*x31        
            +coeff[ 51]    *x21*x31*x42    
            +coeff[ 52]    *x21    *x43    
        ;
        v_x_l5p65_ep29_q3ex               =v_x_l5p65_ep29_q3ex               
            +coeff[ 53]    *x23*x32*x41    
            ;
     
        return v_x_l5p65_ep29_q3ex               ;
    }
    float t_l5p65_ep29_q3ex               (float *x,int m){
        int ncoeff= 80;
        float avdat=  0.6323135E-02;
        float xmin[10]={
            -0.14979E-01,-0.46767E-01,-0.14993E-01,-0.29713E-01,-0.46582E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14981E-01, 0.30136E-01, 0.49667E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 81]={
            -0.15522081E-02,-0.34886049E-02,-0.18747732E-01, 0.11083505E+00,
             0.52433326E-02,-0.10172439E-01, 0.10611741E-02,-0.21129772E-02,
             0.17529317E-03,-0.42694987E-03,-0.74328569E-03,-0.25144108E-02,
             0.24468289E-02,-0.71558770E-03, 0.10676727E-02,-0.32891289E-03,
             0.43957264E-03,-0.59407495E-03,-0.70649519E-04,-0.26440871E-03,
            -0.32769586E-03, 0.14170106E-03,-0.37456208E-03,-0.19619823E-03,
             0.14296289E-03, 0.45999489E-03, 0.22626256E-03, 0.46982875E-03,
             0.22330541E-03, 0.78235212E-03, 0.12323883E-02, 0.32002534E-03,
            -0.95644667E-04, 0.19692905E-03,-0.27000206E-03, 0.13696119E-03,
            -0.82199019E-03,-0.10017589E-03, 0.81869555E-04, 0.11644472E-03,
             0.14621479E-03, 0.11497988E-03, 0.13708865E-03, 0.12090330E-03,
             0.22972001E-03,-0.75356103E-03,-0.44723402E-04, 0.45243566E-03,
            -0.29783320E-03, 0.10214149E-03,-0.23215584E-03, 0.13184900E-03,
             0.19362029E-04, 0.27176808E-04, 0.51999697E-03,-0.79946214E-03,
             0.42064066E-03,-0.38213917E-03, 0.32636477E-03, 0.67798106E-03,
             0.16053877E-04, 0.16754286E-03, 0.58474638E-04,-0.49592440E-04,
            -0.51595918E-04, 0.13255597E-03,-0.50240371E-03, 0.53498443E-04,
             0.96041018E-04, 0.10074986E-03,-0.90300586E-04, 0.81214792E-04,
            -0.12250173E-03,-0.14482051E-03, 0.10484871E-03,-0.21698452E-03,
            -0.15916047E-03,-0.17210760E-03,-0.21174617E-03,-0.13564165E-03,
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
     
        float v_t_l5p65_ep29_q3ex               =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]                *x51
            +coeff[  4]    *x21        *x51
            +coeff[  5]                *x52
            +coeff[  6]    *x22            
            +coeff[  7]            *x42    
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[  8]*x11*x21            
            +coeff[  9]    *x21    *x41    
            +coeff[ 10]        *x31*x41    
            +coeff[ 11]    *x21    *x42    
            +coeff[ 12]    *x21    *x41*x51
            +coeff[ 13]    *x21        *x52
            +coeff[ 14]                *x53
            +coeff[ 15]*x11            *x51
            +coeff[ 16]            *x41*x51
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]*x12                
            +coeff[ 19]    *x21*x31        
            +coeff[ 20]*x11        *x41    
            +coeff[ 21]        *x31    *x51
            +coeff[ 22]*x11*x22            
            +coeff[ 23]    *x22*x31        
            +coeff[ 24]*x11    *x31*x41    
            +coeff[ 25]*x11        *x42    
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[ 26]    *x22        *x51
            +coeff[ 27]    *x21*x31    *x51
            +coeff[ 28]            *x42*x51
            +coeff[ 29]*x11*x23            
            +coeff[ 30]    *x22    *x42    
            +coeff[ 31]    *x23    *x43    
            +coeff[ 32]*x11    *x31        
            +coeff[ 33]    *x23            
            +coeff[ 34]*x11*x21    *x42    
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[ 35]    *x23        *x51
            +coeff[ 36]    *x21    *x42*x51
            +coeff[ 37]    *x21        *x53
            +coeff[ 38]*x13*x22        *x51
            +coeff[ 39]*x11*x21    *x41    
            +coeff[ 40]            *x43    
            +coeff[ 41]*x11    *x31    *x51
            +coeff[ 42]*x11        *x41*x51
            +coeff[ 43]*x11            *x52
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[ 44]    *x23*x31        
            +coeff[ 45]    *x23    *x41    
            +coeff[ 46]*x13            *x51
            +coeff[ 47]*x11*x22        *x51
            +coeff[ 48]    *x21*x31*x41*x51
            +coeff[ 49]        *x31*x41*x52
            +coeff[ 50]            *x42*x52
            +coeff[ 51]*x11*x23*x31        
            +coeff[ 52]*x12    *x31*x42    
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[ 53]        *x33*x42    
            +coeff[ 54]    *x22    *x43    
            +coeff[ 55]    *x23    *x41*x51
            +coeff[ 56]    *x22    *x42*x51
            +coeff[ 57]    *x22    *x41*x52
            +coeff[ 58]*x13*x21    *x41*x51
            +coeff[ 59]    *x23        *x53
            +coeff[ 60]*x12*x21            
            +coeff[ 61]        *x31*x42    
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[ 62]*x12            *x51
            +coeff[ 63]        *x32    *x51
            +coeff[ 64]        *x31    *x52
            +coeff[ 65]*x11*x22*x31        
            +coeff[ 66]    *x21    *x43    
            +coeff[ 67]*x12*x21        *x51
            +coeff[ 68]    *x22*x31    *x51
            +coeff[ 69]*x11*x21        *x52
            +coeff[ 70]    *x21    *x41*x52
        ;
        v_t_l5p65_ep29_q3ex               =v_t_l5p65_ep29_q3ex               
            +coeff[ 71]            *x41*x53
            +coeff[ 72]*x12*x23            
            +coeff[ 73]*x12*x22    *x41    
            +coeff[ 74]    *x22*x32*x41    
            +coeff[ 75]*x11*x22    *x42    
            +coeff[ 76]*x11    *x31*x43    
            +coeff[ 77]*x12*x22        *x51
            +coeff[ 78]    *x23*x31    *x51
            +coeff[ 79]*x11*x21    *x42*x51
            ;
     
        return v_t_l5p65_ep29_q3ex               ;
    }
    float y_l5p65_ep29_q3ex               (float *x,int m){
        int ncoeff= 97;
        float avdat= -0.2989810E-02;
        float xmin[10]={
            -0.14979E-01,-0.46767E-01,-0.14993E-01,-0.29713E-01,-0.46582E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14981E-01, 0.30136E-01, 0.49667E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 98]={
            -0.14336201E-02, 0.36300026E-01, 0.13025569E-02,-0.20198669E-01,
            -0.20182509E-01,-0.44524509E-01,-0.75232899E-02, 0.27673883E-01,
            -0.47023809E-02, 0.64263963E-02, 0.55838015E-01, 0.11730332E-01,
             0.63276510E-02, 0.58060647E-02,-0.14023820E-01, 0.13666422E-02,
            -0.63244659E-02,-0.37705906E-01,-0.49919626E-02,-0.67093717E-02,
            -0.76174554E-02,-0.10910170E-03, 0.12318593E-01, 0.33456839E-02,
             0.45037019E-03, 0.21887472E-03,-0.10840774E-01,-0.12813548E-02,
            -0.51562926E-02,-0.39513749E-02,-0.58898358E-02,-0.23353384E-02,
            -0.88553475E-02,-0.63023400E-02, 0.80890337E-03,-0.70719951E-03,
            -0.18524365E-02,-0.17991684E-02,-0.16383857E-02,-0.13713912E-02,
            -0.47870213E-02,-0.33901627E-02, 0.43688258E-03,-0.49788910E-02,
            -0.15116854E-02, 0.74576656E-02,-0.11027348E-02, 0.40247659E-02,
             0.10078853E-02, 0.54130529E-03, 0.68884151E-03,-0.18537010E-02,
            -0.18030216E-02, 0.50407378E-02,-0.10886515E-02,-0.10720539E-02,
             0.49805228E-03, 0.31389398E-03,-0.35366148E-03, 0.21736156E-02,
             0.20787225E-02,-0.25054861E-02, 0.68307173E-03, 0.11385452E-02,
             0.22952564E-02,-0.10522903E-02,-0.85012993E-03,-0.10199698E-02,
            -0.75366441E-03,-0.35464216E-02,-0.93211367E-03, 0.24799621E-02,
             0.24494224E-02,-0.61626243E-03,-0.24000583E-02, 0.18603008E-03,
            -0.19366809E-03,-0.30112659E-03,-0.13779050E-03, 0.98678493E-03,
            -0.13440595E-02,-0.49033825E-03, 0.35914624E-03,-0.29934518E-03,
            -0.14492590E-03, 0.43995774E-03, 0.23733011E-03, 0.96285390E-03,
             0.53298479E-03, 0.62686042E-03, 0.19402724E-02, 0.74641936E-03,
             0.65085030E-03,-0.16922235E-02,-0.76753792E-03, 0.41499548E-03,
            -0.42526962E-03,
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
     
        float v_y_l5p65_ep29_q3ex               =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]        *x31        
            +coeff[  4]                *x51
            +coeff[  5]    *x21    *x41    
            +coeff[  6]        *x31*x41    
            +coeff[  7]    *x22            
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[  8]    *x21*x31        
            +coeff[  9]*x11        *x41    
            +coeff[ 10]            *x41*x51
            +coeff[ 11]*x11*x21            
            +coeff[ 12]    *x21        *x51
            +coeff[ 13]        *x31    *x51
            +coeff[ 14]    *x21    *x42    
            +coeff[ 15]*x12                
            +coeff[ 16]                *x52
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]    *x21*x31*x41    
            +coeff[ 19]    *x23            
            +coeff[ 20]    *x22*x31        
            +coeff[ 21]            *x44    
            +coeff[ 22]    *x22        *x51
            +coeff[ 23]*x11*x21        *x51
            +coeff[ 24]    *x21    *x43*x51
            +coeff[ 25]*x11                
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 26]            *x42    
            +coeff[ 27]        *x32        
            +coeff[ 28]            *x43    
            +coeff[ 29]*x11        *x42    
            +coeff[ 30]    *x21    *x41*x51
            +coeff[ 31]*x11*x21*x31        
            +coeff[ 32]    *x24            
            +coeff[ 33]    *x22    *x41*x51
            +coeff[ 34]*x11    *x31        
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 35]*x11            *x51
            +coeff[ 36]*x11*x21    *x41    
            +coeff[ 37]*x11    *x31*x41    
            +coeff[ 38]*x11*x22            
            +coeff[ 39]*x12        *x41    
            +coeff[ 40]    *x22    *x42    
            +coeff[ 41]*x11*x22    *x41    
            +coeff[ 42]        *x31*x44    
            +coeff[ 43]*x11*x23            
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 44]        *x31*x41*x51
            +coeff[ 45]    *x21    *x43    
            +coeff[ 46]            *x41*x52
            +coeff[ 47]    *x21*x31*x42    
            +coeff[ 48]    *x21        *x52
            +coeff[ 49]*x12            *x51
            +coeff[ 50]                *x53
            +coeff[ 51]    *x22*x31    *x51
            +coeff[ 52]*x11*x21    *x41*x51
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 53]    *x23    *x42    
            +coeff[ 54]        *x31*x42    
            +coeff[ 55]            *x42*x51
            +coeff[ 56]*x11        *x41*x51
            +coeff[ 57]*x12*x21            
            +coeff[ 58]*x12    *x31        
            +coeff[ 59]    *x23    *x41    
            +coeff[ 60]    *x22*x31*x41    
            +coeff[ 61]    *x21    *x42*x51
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 62]    *x23*x31        
            +coeff[ 63]    *x22*x32        
            +coeff[ 64]*x11*x21*x31*x41    
            +coeff[ 65]    *x23        *x51
            +coeff[ 66]*x12*x21    *x41    
            +coeff[ 67]*x12*x22            
            +coeff[ 68]*x11*x21*x31    *x51
            +coeff[ 69]    *x24    *x41    
            +coeff[ 70]            *x41*x53
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 71]*x11*x22    *x42    
            +coeff[ 72]    *x22    *x42*x51
            +coeff[ 73]*x11*x21        *x52
            +coeff[ 74]*x11*x23    *x41    
            +coeff[ 75]        *x31*x42*x52
            +coeff[ 76]    *x23*x31    *x51
            +coeff[ 77]        *x34    *x51
            +coeff[ 78]*x11*x21    *x44    
            +coeff[ 79]*x11*x23*x32        
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 80]    *x21*x31*x41*x53
            +coeff[ 81]*x12*x22        *x52
            +coeff[ 82]        *x32*x41    
            +coeff[ 83]    *x21*x32        
            +coeff[ 84]*x11    *x32        
            +coeff[ 85]    *x21*x31    *x51
            +coeff[ 86]        *x31    *x52
            +coeff[ 87]*x11        *x43    
            +coeff[ 88]            *x43*x51
        ;
        v_y_l5p65_ep29_q3ex               =v_y_l5p65_ep29_q3ex               
            +coeff[ 89]    *x21*x32*x41    
            +coeff[ 90]*x11*x21    *x42    
            +coeff[ 91]*x11    *x31*x42    
            +coeff[ 92]*x12        *x42    
            +coeff[ 93]    *x22    *x43    
            +coeff[ 94]    *x21*x31*x43    
            +coeff[ 95]*x12    *x31*x41    
            +coeff[ 96]    *x22        *x52
            ;
     
        return v_y_l5p65_ep29_q3ex               ;
    }
    float p_l5p65_ep29_q3ex               (float *x,int m){
        int ncoeff= 90;
        float avdat= -0.2037797E-03;
        float xmin[10]={
            -0.14979E-01,-0.46767E-01,-0.14993E-01,-0.29713E-01,-0.46582E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14981E-01, 0.30136E-01, 0.49667E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 91]={
             0.21447118E-02, 0.62491517E-02,-0.33592813E-01, 0.10495326E-01,
            -0.11112604E-03,-0.14906316E-01, 0.14462954E-02, 0.18739611E-01,
             0.41195368E-02,-0.26083533E-02,-0.35977615E-02,-0.20761052E-01,
            -0.63436837E-02, 0.30377605E-02,-0.28784492E-02, 0.17296594E-01,
             0.60862256E-02,-0.35203518E-02,-0.17324793E-03,-0.13633736E-03,
            -0.25086457E-03,-0.23648512E-03, 0.74407413E-04, 0.60278415E-02,
             0.15306278E-02, 0.42952369E-02, 0.22749184E-02, 0.22283669E-02,
            -0.80306362E-03, 0.38391065E-02, 0.24939915E-02, 0.27507499E-02,
            -0.76139148E-03, 0.70715969E-03, 0.32906909E-03, 0.12761109E-02,
             0.81825122E-03, 0.17397581E-02,-0.79708797E-03, 0.24530026E-02,
             0.56776736E-03,-0.33879926E-03, 0.40160626E-03, 0.15471141E-03,
            -0.41389582E-03, 0.65055996E-03, 0.70920936E-03,-0.19069517E-02,
            -0.30036187E-02, 0.10366493E-03, 0.15149629E-02, 0.67177898E-03,
            -0.99742704E-03,-0.11897220E-02,-0.66831498E-03,-0.60551963E-03,
            -0.71877264E-04, 0.66048250E-03,-0.17193357E-03,-0.43204793E-03,
             0.17095705E-03,-0.13730824E-02,-0.15277719E-02,-0.22716050E-02,
            -0.15233012E-03, 0.66683040E-03, 0.36890156E-04, 0.35531732E-03,
            -0.13205977E-02, 0.72088116E-03,-0.59260614E-03,-0.18827226E-03,
             0.15824460E-03, 0.70492853E-04,-0.25742917E-03,-0.10570374E-02,
            -0.10388219E-03, 0.23024516E-03,-0.34379429E-03,-0.12476397E-03,
             0.20749598E-03, 0.25409099E-03, 0.13979757E-03,-0.49911829E-03,
            -0.16180841E-03,-0.23677867E-03, 0.73952001E-03, 0.21549068E-03,
            -0.18881859E-03,-0.24367377E-03,
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
        float x44 = x43*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_p_l5p65_ep29_q3ex               =avdat
            +coeff[  0]                    
            +coeff[  1]        *x31        
            +coeff[  2]            *x41    
            +coeff[  3]                *x51
            +coeff[  4]*x11                
            +coeff[  5]    *x22            
            +coeff[  6]    *x21*x31        
            +coeff[  7]    *x21    *x41    
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[  8]        *x31*x41    
            +coeff[  9]    *x21        *x51
            +coeff[ 10]        *x31    *x51
            +coeff[ 11]            *x41*x51
            +coeff[ 12]*x11*x21            
            +coeff[ 13]    *x23            
            +coeff[ 14]*x11        *x41    
            +coeff[ 15]    *x22    *x41    
            +coeff[ 16]    *x21    *x42    
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 17]    *x22        *x51
            +coeff[ 18]    *x22    *x42    
            +coeff[ 19]        *x32*x42    
            +coeff[ 20]            *x44    
            +coeff[ 21]    *x24*x31        
            +coeff[ 22]            *x42*x52
            +coeff[ 23]            *x42    
            +coeff[ 24]                *x52
            +coeff[ 25]    *x22*x31        
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 26]    *x21*x31*x41    
            +coeff[ 27]            *x43    
            +coeff[ 28]*x12                
            +coeff[ 29]    *x24            
            +coeff[ 30]*x11*x21    *x41    
            +coeff[ 31]            *x41*x52
            +coeff[ 32]    *x21            
            +coeff[ 33]        *x32        
            +coeff[ 34]*x11            *x51
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 35]*x11*x21*x31        
            +coeff[ 36]        *x31    *x52
            +coeff[ 37]*x11        *x42    
            +coeff[ 38]*x11*x21        *x51
            +coeff[ 39]*x11*x23            
            +coeff[ 40]*x12        *x41    
            +coeff[ 41]*x11    *x31        
            +coeff[ 42]        *x31*x42    
            +coeff[ 43]        *x31*x41*x51
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 44]            *x42*x51
            +coeff[ 45]*x11*x22            
            +coeff[ 46]*x11    *x31*x41    
            +coeff[ 47]    *x22*x31*x41    
            +coeff[ 48]    *x21    *x43    
            +coeff[ 49]*x11        *x41*x51
            +coeff[ 50]*x11*x22    *x41    
            +coeff[ 51]*x11*x21    *x41*x51
            +coeff[ 52]            *x41*x53
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 53]    *x23*x31*x42    
            +coeff[ 54]    *x23    *x41    
            +coeff[ 55]    *x22*x32        
            +coeff[ 56]                *x53
            +coeff[ 57]    *x22    *x41*x51
            +coeff[ 58]*x12*x21            
            +coeff[ 59]            *x43*x51
            +coeff[ 60]*x12    *x31        
            +coeff[ 61]*x11*x21*x31*x41    
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 62]*x11*x21    *x42    
            +coeff[ 63]    *x23    *x42    
            +coeff[ 64]*x12            *x51
            +coeff[ 65]*x12*x22            
            +coeff[ 66]        *x34*x41    
            +coeff[ 67]*x12*x21    *x41    
            +coeff[ 68]*x11*x22    *x42    
            +coeff[ 69]    *x26    *x41    
            +coeff[ 70]*x11*x23*x32        
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 71]        *x32*x41    
            +coeff[ 72]    *x21        *x52
            +coeff[ 73]*x11    *x32        
            +coeff[ 74]    *x23        *x51
            +coeff[ 75]    *x21*x31*x42    
            +coeff[ 76]*x11    *x31    *x51
            +coeff[ 77]    *x22*x31    *x51
            +coeff[ 78]        *x31*x43    
            +coeff[ 79]*x11            *x52
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 80]*x11*x22*x31        
            +coeff[ 81]    *x21    *x41*x52
            +coeff[ 82]    *x23*x32        
            +coeff[ 83]    *x23*x31*x41    
            +coeff[ 84]    *x21        *x53
            +coeff[ 85]*x11        *x43    
            +coeff[ 86]    *x22    *x43    
            +coeff[ 87]*x11*x21*x31    *x51
            +coeff[ 88]    *x23*x31    *x51
        ;
        v_p_l5p65_ep29_q3ex               =v_p_l5p65_ep29_q3ex               
            +coeff[ 89]        *x31    *x53
            ;
     
        return v_p_l5p65_ep29_q3ex               ;
    }
    float l_l5p65_ep29_q3ex               (float *x,int m){
        int ncoeff= 72;
        float avdat= -0.1098074E-01;
        float xmin[10]={
            -0.14979E-01,-0.46767E-01,-0.14993E-01,-0.29713E-01,-0.46582E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.45617E-01, 0.14981E-01, 0.30136E-01, 0.49667E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 73]={
             0.13070384E-01,-0.24838650E+00,-0.32277022E-01, 0.39042164E-01,
            -0.21587769E-01,-0.91116183E-01,-0.20213049E-01, 0.99459905E-02,
            -0.31685527E-01,-0.87507861E-02,-0.89013176E-02,-0.81140390E-02,
             0.61292681E-02,-0.72188322E-02,-0.23192311E-01,-0.23630623E-01,
             0.36054358E-01,-0.15307121E-01,-0.14941675E-02, 0.27951028E-02,
            -0.26865562E-02,-0.65497248E-02, 0.22811227E-01, 0.34118481E-01,
             0.21799425E-01, 0.21986288E-02, 0.39630504E-02,-0.47971481E-02,
             0.59732869E-02,-0.38509001E-02, 0.69224765E-02, 0.14793170E-03,
             0.84401539E-03,-0.19461086E-02, 0.13841289E-02, 0.41257041E-02,
            -0.41856985E-02, 0.44342000E-02, 0.18846625E-03,-0.56227687E-03,
             0.85914269E-03,-0.10351174E-02, 0.73052006E-03, 0.48916857E-03,
             0.23561285E-02,-0.10790171E-01,-0.10255579E-01,-0.27271870E-02,
             0.14731887E-02, 0.53107301E-02,-0.19259978E-03,-0.18083270E-02,
             0.41356060E-03, 0.38713482E-02,-0.64003840E-03,-0.18426683E-03,
             0.14359327E-02,-0.73480001E-03, 0.11449030E-02, 0.19993782E-02,
             0.47098013E-03, 0.31188226E-03,-0.31420847E-02, 0.54926784E-02,
             0.32351571E-02, 0.23650758E-02, 0.92130329E-03, 0.14779010E-02,
             0.19294524E-02,-0.90478628E-03, 0.43988326E-02, 0.34612236E-02,
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
     
        float v_l_l5p65_ep29_q3ex               =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]                *x51
            +coeff[  3]*x11                
            +coeff[  4]    *x22            
            +coeff[  5]    *x21    *x41    
            +coeff[  6]*x11        *x41    
            +coeff[  7]    *x23*x31        
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[  8]    *x21*x31        
            +coeff[  9]            *x42    
            +coeff[ 10]    *x21        *x51
            +coeff[ 11]                *x52
            +coeff[ 12]*x11*x21            
            +coeff[ 13]*x11    *x31        
            +coeff[ 14]    *x23            
            +coeff[ 15]    *x22    *x41    
            +coeff[ 16]    *x21    *x42    
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[ 17]*x11*x22            
            +coeff[ 18]    *x23*x31*x41    
            +coeff[ 19]            *x41    
            +coeff[ 20]*x11            *x51
            +coeff[ 21]    *x22*x31        
            +coeff[ 22]    *x21*x31*x41    
            +coeff[ 23]    *x23    *x41    
            +coeff[ 24]*x11*x22    *x41    
            +coeff[ 25]            *x41*x51
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[ 26]    *x21*x32        
            +coeff[ 27]*x11*x21    *x41    
            +coeff[ 28]*x11        *x42    
            +coeff[ 29]*x12*x21            
            +coeff[ 30]*x11*x22*x31        
            +coeff[ 31]*x11    *x33*x41    
            +coeff[ 32]        *x31        
            +coeff[ 33]        *x31*x41    
            +coeff[ 34]                *x53
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[ 35]*x11    *x31*x41    
            +coeff[ 36]    *x24            
            +coeff[ 37]*x12*x21    *x41    
            +coeff[ 38]        *x31    *x51
            +coeff[ 39]*x12                
            +coeff[ 40]            *x42*x51
            +coeff[ 41]*x11*x21*x31        
            +coeff[ 42]*x11    *x32        
            +coeff[ 43]*x11        *x41*x51
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[ 44]    *x22    *x42    
            +coeff[ 45]    *x21*x31*x42    
            +coeff[ 46]    *x21    *x43    
            +coeff[ 47]*x11*x23            
            +coeff[ 48]*x12*x21*x31        
            +coeff[ 49]    *x21    *x44    
            +coeff[ 50]    *x23*x31    *x51
            +coeff[ 51]    *x23    *x41*x51
            +coeff[ 52]    *x22*x31*x41*x51
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[ 53]    *x24*x31*x41    
            +coeff[ 54]    *x23*x32*x41    
            +coeff[ 55]        *x32        
            +coeff[ 56]            *x43    
            +coeff[ 57]    *x22        *x51
            +coeff[ 58]    *x21*x31    *x51
            +coeff[ 59]    *x21    *x41*x51
            +coeff[ 60]        *x31*x41*x51
            +coeff[ 61]*x11            *x52
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[ 62]    *x21*x32*x41    
            +coeff[ 63]    *x24    *x41    
            +coeff[ 64]    *x22*x31*x42    
            +coeff[ 65]    *x21*x31*x43    
            +coeff[ 66]        *x31*x44    
            +coeff[ 67]*x11*x24            
            +coeff[ 68]*x11*x23    *x41    
            +coeff[ 69]*x13*x22            
            +coeff[ 70]    *x24    *x42    
        ;
        v_l_l5p65_ep29_q3ex               =v_l_l5p65_ep29_q3ex               
            +coeff[ 71]    *x22*x32*x43    
            ;
     
        return v_l_l5p65_ep29_q3ex               ;
    }
    float x_l5p65_ep5                     (float *x,int m){
        int ncoeff=  7;
        float avdat=  0.1131176E+00;
        float xmin[10]={
            -0.14999E-01,-0.50000E-01,-0.14977E-01,-0.32453E-01,-0.49984E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.50023E-01, 0.14991E-01, 0.30973E-01, 0.49958E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[  8]={
            -0.66975919E-02, 0.99990712E-06, 0.37790876E-05, 0.15064772E-01,
             0.34694370E-01, 0.10470957E-03, 0.43554624E-04,
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
        float x31 = x3;
        float x41 = x4;
        float x42 = x41*x4;
     
    //                 function
     
        float v_x_l5p65_ep5                     =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]        *x31        
            +coeff[  4]            *x41    
            +coeff[  5]            *x42    
            +coeff[  6]        *x31*x41    
            ;
     
        return v_x_l5p65_ep5                     ;
    }
    float t_l5p65_ep5                     (float *x,int m){
        int ncoeff= 50;
        float avdat=  0.1057637E+00;
        float xmin[10]={
            -0.14999E-01,-0.50000E-01,-0.14977E-01,-0.32453E-01,-0.49984E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.50023E-01, 0.14991E-01, 0.30973E-01, 0.49958E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 51]={
            -0.52106064E-02, 0.25526482E-04, 0.92808608E-04, 0.26578573E-03,
             0.32633398E-01,-0.11834662E-03,-0.15561591E-03,-0.31227592E-03,
             0.65363594E-04,-0.13727882E-03,-0.19829896E-03,-0.46568252E-05,
             0.78755409E-04,-0.34646775E-04,-0.24054316E-04,-0.83435698E-04,
            -0.65716253E-04,-0.95291307E-05, 0.10227678E-03, 0.75973658E-04,
             0.11991811E-03, 0.23311479E-05, 0.55175265E-05,-0.32692913E-06,
            -0.50293206E-05, 0.15177608E-04,-0.66670807E-04,-0.48406237E-05,
             0.55841347E-05, 0.10423521E-04, 0.11384315E-04,-0.29500334E-04,
             0.32890297E-04, 0.21415708E-04,-0.27511831E-04, 0.69806061E-04,
             0.44853281E-04, 0.66981443E-05, 0.66323873E-05,-0.64574683E-05,
             0.13370210E-04,-0.38028063E-05, 0.14059205E-04, 0.17556500E-04,
            -0.55983728E-05, 0.21455013E-04, 0.20994084E-04,-0.60180437E-05,
            -0.98074752E-05,-0.26186008E-04,
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
     
        float v_t_l5p65_ep5                     =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]        *x31        
            +coeff[  4]            *x41    
            +coeff[  5]                *x51
            +coeff[  6]        *x31*x41    
            +coeff[  7]    *x22    *x41    
        ;
        v_t_l5p65_ep5                     =v_t_l5p65_ep5                     
            +coeff[  8]    *x22            
            +coeff[  9]    *x22*x31        
            +coeff[ 10]*x11*x21    *x41    
            +coeff[ 11]*x13*x21            
            +coeff[ 12]*x11*x21            
            +coeff[ 13]        *x32        
            +coeff[ 14]            *x41*x51
            +coeff[ 15]*x11*x21*x31        
            +coeff[ 16]*x11*x23            
        ;
        v_t_l5p65_ep5                     =v_t_l5p65_ep5                     
            +coeff[ 17]        *x31    *x51
            +coeff[ 18]    *x22*x31*x41    
            +coeff[ 19]*x11*x21    *x42    
            +coeff[ 20]    *x22    *x42    
            +coeff[ 21]        *x32*x42    
            +coeff[ 22]*x12        *x43    
            +coeff[ 23]*x13*x21*x31*x41    
            +coeff[ 24]    *x22*x32*x42    
            +coeff[ 25]*x12                
        ;
        v_t_l5p65_ep5                     =v_t_l5p65_ep5                     
            +coeff[ 26]            *x42    
            +coeff[ 27]    *x21        *x51
            +coeff[ 28]                *x52
            +coeff[ 29]*x11*x22            
            +coeff[ 30]    *x23            
            +coeff[ 31]*x12        *x41    
            +coeff[ 32]        *x31*x42    
            +coeff[ 33]            *x43    
            +coeff[ 34]*x12*x22            
        ;
        v_t_l5p65_ep5                     =v_t_l5p65_ep5                     
            +coeff[ 35]*x11*x21*x31*x41    
            +coeff[ 36]*x11*x23    *x41    
            +coeff[ 37]*x12    *x32*x41    
            +coeff[ 38]*x12*x21            
            +coeff[ 39]*x12    *x31        
            +coeff[ 40]        *x32*x41    
            +coeff[ 41]*x11*x21        *x51
            +coeff[ 42]*x11*x21*x32        
            +coeff[ 43]    *x22*x32        
        ;
        v_t_l5p65_ep5                     =v_t_l5p65_ep5                     
            +coeff[ 44]*x11*x22    *x41    
            +coeff[ 45]*x11*x23*x31        
            +coeff[ 46]*x12*x22    *x41    
            +coeff[ 47]    *x23*x31*x41    
            +coeff[ 48]*x12*x21    *x42    
            +coeff[ 49]*x11*x21*x31*x42    
            ;
     
        return v_t_l5p65_ep5                     ;
    }
    float y_l5p65_ep5                     (float *x,int m){
        int ncoeff=  6;
        float avdat=  0.1132580E-02;
        float xmin[10]={
            -0.14999E-01,-0.50000E-01,-0.14977E-01,-0.32453E-01,-0.49984E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.50023E-01, 0.14991E-01, 0.30973E-01, 0.49958E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[  7]={
            -0.11139134E-02,-0.20130881E-05,-0.54391790E-01,-0.14972868E-01,
            -0.18853023E-03,-0.81170117E-04,
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
        float x31 = x3;
        float x41 = x4;
     
    //                 function
     
        float v_y_l5p65_ep5                     =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]*x11                
            +coeff[  4]    *x21    *x41    
            +coeff[  5]    *x21*x31        
            ;
     
        return v_y_l5p65_ep5                     ;
    }
    float p_l5p65_ep5                     (float *x,int m){
        int ncoeff= 10;
        float avdat=  0.1297381E-02;
        float xmin[10]={
            -0.14999E-01,-0.50000E-01,-0.14977E-01,-0.32453E-01,-0.49984E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.50023E-01, 0.14991E-01, 0.30973E-01, 0.49958E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 11]={
            -0.10283185E-02,-0.49456440E-01,-0.64673298E-03, 0.23015760E-03,
            -0.18912545E-03,-0.12906549E-03,-0.43643908E-04,-0.14367150E-03,
            -0.49337075E-04,-0.11970787E-03,
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
     
        float v_p_l5p65_ep5                     =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]    *x21    *x41    
            +coeff[  3]*x11                
            +coeff[  4]    *x21*x31        
            +coeff[  5]*x11        *x41    
            +coeff[  6]            *x41    
            +coeff[  7]    *x23            
        ;
        v_p_l5p65_ep5                     =v_p_l5p65_ep5                     
            +coeff[  8]*x11    *x31        
            +coeff[  9]*x11*x22            
            ;
     
        return v_p_l5p65_ep5                     ;
    }
    float l_l5p65_ep5                     (float *x,int m){
        int ncoeff=  8;
        float avdat= -0.1172051E-02;
        float xmin[10]={
            -0.14999E-01,-0.50000E-01,-0.14977E-01,-0.32453E-01,-0.49984E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14960E-01, 0.50023E-01, 0.14991E-01, 0.30973E-01, 0.49958E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[  9]={
             0.12507946E-02,-0.49669978E-06,-0.14839326E-02,-0.33939756E-02,
            -0.13590114E-02,-0.49256541E-05,-0.55832946E-03,-0.36828210E-05,
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
     
        float v_l_l5p65_ep5                     =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]    *x22            
            +coeff[  5]        *x31*x41    
            +coeff[  6]            *x42    
            +coeff[  7]    *x22    *x41    
            ;
     
        return v_l_l5p65_ep5                     ;
    }
    float x_l5p65_ep7                     (float *x,int m){
        int ncoeff= 40;
        float avdat=  0.2707396E+00;
        float xmin[10]={
            -0.14979E-01,-0.49703E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.49762E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 41]={
            -0.13375221E-01,-0.29272670E-02, 0.35733142E-03, 0.16751107E-01,
             0.71168683E-01, 0.51578628E-02, 0.23748705E-02,-0.41578296E-02,
            -0.15970954E-02,-0.18701813E-02,-0.21523774E-04, 0.27313476E-03,
            -0.13087519E-02,-0.14953233E-02, 0.96151663E-04, 0.14967607E-03,
            -0.19542921E-03,-0.23384450E-03,-0.66843227E-03, 0.19815390E-02,
            -0.69553658E-04,-0.18766460E-03,-0.23634567E-03, 0.83144271E-03,
             0.13324155E-02, 0.19915486E-03,-0.27812334E-05,-0.42349602E-04,
            -0.55620989E-04,-0.11160655E-03, 0.83200961E-04,-0.15212891E-03,
             0.17492268E-03,-0.18796601E-03, 0.53715560E-03, 0.15655186E-03,
             0.12499611E-03, 0.12082088E-03,-0.16075092E-03, 0.11979001E-03,
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
     
        float v_x_l5p65_ep7                     =avdat
            +coeff[  0]                    
            +coeff[  1]                *x51
            +coeff[  2]    *x21            
            +coeff[  3]        *x31        
            +coeff[  4]            *x41    
            +coeff[  5]    *x22            
            +coeff[  6]*x11*x21            
            +coeff[  7]    *x22    *x41    
        ;
        v_x_l5p65_ep7                     =v_x_l5p65_ep7                     
            +coeff[  8]            *x42    
            +coeff[  9]*x11*x21    *x41    
            +coeff[ 10]        *x33*x41    
            +coeff[ 11]*x12                
            +coeff[ 12]        *x31*x41    
            +coeff[ 13]    *x22*x31        
            +coeff[ 14]*x11                
            +coeff[ 15]                *x52
            +coeff[ 16]            *x41*x51
        ;
        v_x_l5p65_ep7                     =v_x_l5p65_ep7                     
            +coeff[ 17]        *x32        
            +coeff[ 18]*x11*x21*x31        
            +coeff[ 19]    *x22    *x42    
            +coeff[ 20]    *x22*x33*x41    
            +coeff[ 21]*x12        *x41    
            +coeff[ 22]    *x22        *x51
            +coeff[ 23]*x11*x21    *x42    
            +coeff[ 24]    *x22*x31*x41    
            +coeff[ 25]        *x31*x43    
        ;
        v_x_l5p65_ep7                     =v_x_l5p65_ep7                     
            +coeff[ 26]*x13*x21*x31*x41    
            +coeff[ 27]        *x31    *x51
            +coeff[ 28]*x12    *x31        
            +coeff[ 29]*x11*x21        *x51
            +coeff[ 30]    *x23            
            +coeff[ 31]    *x21    *x42    
            +coeff[ 32]            *x43    
            +coeff[ 33]*x11*x23            
            +coeff[ 34]*x11*x21*x31*x41    
        ;
        v_x_l5p65_ep7                     =v_x_l5p65_ep7                     
            +coeff[ 35]    *x22*x32        
            +coeff[ 36]    *x23    *x41    
            +coeff[ 37]        *x32*x42    
            +coeff[ 38]    *x23*x31*x41    
            +coeff[ 39]        *x33*x42    
            ;
     
        return v_x_l5p65_ep7                     ;
    }
    float t_l5p65_ep7                     (float *x,int m){
        int ncoeff= 50;
        float avdat=  0.2248077E+00;
        float xmin[10]={
            -0.14979E-01,-0.49703E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.49762E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 51]={
            -0.83040716E-02, 0.50986884E-03, 0.23339810E-02, 0.40491834E-01,
            -0.61374730E-02, 0.10135975E-01,-0.64502461E-02,-0.60507791E-04,
             0.45202469E-03, 0.41900803E-02,-0.35653634E-02,-0.28153777E-02,
             0.52220421E-03, 0.14374698E-04,-0.24879905E-02, 0.31388548E-03,
            -0.22651365E-02, 0.42013009E-02, 0.29262883E-04,-0.22904034E-03,
             0.13130707E-03,-0.40872110E-03,-0.31297732E-03,-0.46466049E-03,
             0.56115381E-03,-0.50536689E-03, 0.25443784E-02, 0.16157739E-02,
             0.13804850E-03, 0.35392601E-03,-0.10157872E-02,-0.31952429E-03,
            -0.25618187E-03,-0.20719408E-03, 0.92397392E-03,-0.46749923E-04,
             0.20246464E-03,-0.11349138E-03, 0.17543166E-03,-0.95631476E-04,
            -0.35890628E-03, 0.34539495E-03, 0.45147334E-03, 0.24125902E-04,
             0.51535837E-04, 0.70318587E-04,-0.15179499E-03, 0.92375019E-04,
             0.14734104E-03,-0.38122398E-04,
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
     
        float v_t_l5p65_ep7                     =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]    *x22    *x41    
            +coeff[  7]*x13*x21            
        ;
        v_t_l5p65_ep7                     =v_t_l5p65_ep7                     
            +coeff[  8]*x12                
            +coeff[  9]*x11*x21            
            +coeff[ 10]            *x42    
            +coeff[ 11]*x11*x21    *x41    
            +coeff[ 12]        *x31*x42    
            +coeff[ 13]        *x33*x41    
            +coeff[ 14]        *x31*x41    
            +coeff[ 15]                *x52
            +coeff[ 16]    *x22*x31        
        ;
        v_t_l5p65_ep7                     =v_t_l5p65_ep7                     
            +coeff[ 17]    *x22    *x42    
            +coeff[ 18]*x13*x21*x31        
            +coeff[ 19]    *x22*x33*x41    
            +coeff[ 20]*x11                
            +coeff[ 21]        *x32        
            +coeff[ 22]            *x41*x51
            +coeff[ 23]    *x21    *x42    
            +coeff[ 24]            *x43    
            +coeff[ 25]    *x22        *x51
        ;
        v_t_l5p65_ep7                     =v_t_l5p65_ep7                     
            +coeff[ 26]    *x22*x31*x41    
            +coeff[ 27]*x11*x21    *x42    
            +coeff[ 28]*x13*x21*x31*x41    
            +coeff[ 29]    *x23            
            +coeff[ 30]*x11*x21*x31        
            +coeff[ 31]*x12        *x41    
            +coeff[ 32]    *x21*x31*x41    
            +coeff[ 33]*x11*x21        *x51
            +coeff[ 34]*x11*x21*x31*x41    
        ;
        v_t_l5p65_ep7                     =v_t_l5p65_ep7                     
            +coeff[ 35]        *x31    *x51
            +coeff[ 36]*x11*x22            
            +coeff[ 37]*x12    *x31        
            +coeff[ 38]        *x32*x41    
            +coeff[ 39]*x11        *x42    
            +coeff[ 40]*x11*x23            
            +coeff[ 41]    *x22*x32        
            +coeff[ 42]*x12*x22*x31*x41    
            +coeff[ 43]    *x21*x31        
        ;
        v_t_l5p65_ep7                     =v_t_l5p65_ep7                     
            +coeff[ 44]*x12*x21            
            +coeff[ 45]            *x42*x51
            +coeff[ 46]*x12*x22            
            +coeff[ 47]*x11*x21*x32        
            +coeff[ 48]*x12        *x42    
            +coeff[ 49]    *x23        *x51
            ;
     
        return v_t_l5p65_ep7                     ;
    }
    float y_l5p65_ep7                     (float *x,int m){
        int ncoeff= 92;
        float avdat=  0.2196197E-02;
        float xmin[10]={
            -0.14979E-01,-0.49703E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.49762E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 93]={
            -0.10215237E-02,-0.25261901E-03,-0.96336514E-01,-0.13535958E-01,
            -0.70869029E-02,-0.24047629E-02,-0.15162870E-02, 0.28534664E-02,
             0.20392726E-02,-0.21884316E-02,-0.15185145E-02,-0.93926152E-04,
            -0.53770054E-03, 0.22353036E-02, 0.14788499E-02,-0.59045666E-04,
            -0.25092441E-03,-0.35330819E-03, 0.63058268E-03, 0.35846914E-03,
             0.45171843E-03, 0.28729034E-03,-0.34737811E-03,-0.21770742E-04,
             0.81849757E-04,-0.66445224E-04,-0.88059535E-03,-0.91525808E-03,
            -0.28882720E-03, 0.75278111E-03, 0.50756935E-03, 0.31467315E-03,
            -0.12201146E-03,-0.14195246E-03, 0.80772450E-04, 0.88361820E-04,
             0.56535719E-04,-0.43268727E-04, 0.11103502E-03,-0.19269366E-05,
             0.38003924E-04,-0.47666283E-04, 0.12155854E-03,-0.16626924E-03,
             0.77750490E-04,-0.16821449E-03, 0.13968022E-04,-0.54023060E-03,
            -0.34982961E-03,-0.35965428E-03,-0.25588615E-03, 0.50581971E-05,
             0.42166597E-04,-0.52831415E-05,-0.10124463E-04, 0.18077380E-04,
            -0.31246561E-05, 0.48515667E-04,-0.65678134E-04,-0.23215567E-04,
            -0.58642763E-04,-0.27356833E-04, 0.28469081E-04,-0.55022236E-04,
            -0.38182407E-04,-0.39146807E-05, 0.86967812E-04, 0.62578452E-04,
            -0.14304611E-03,-0.37079317E-04,-0.29503432E-04, 0.80737536E-05,
             0.41724281E-04, 0.12500555E-04, 0.11533815E-04,-0.62027875E-05,
             0.66484870E-04,-0.15730851E-04, 0.37761121E-04, 0.16559929E-04,
            -0.57448218E-04, 0.14033322E-04,-0.57298450E-04,-0.87838140E-04,
             0.93053088E-04,-0.36397119E-04,-0.35604186E-04, 0.87762448E-04,
             0.42135380E-04, 0.15518157E-04, 0.89917950E-04,-0.65956519E-04,
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
     
    //                 function
     
        float v_y_l5p65_ep7                     =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]*x11                
            +coeff[  4]    *x21    *x41    
            +coeff[  5]    *x21*x31        
            +coeff[  6]*x11        *x41    
            +coeff[  7]    *x21    *x42    
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[  8]    *x21*x31*x41    
            +coeff[  9]    *x23            
            +coeff[ 10]*x11*x22            
            +coeff[ 11]        *x31        
            +coeff[ 12]*x11    *x31        
            +coeff[ 13]    *x23    *x41    
            +coeff[ 14]*x11*x22    *x41    
            +coeff[ 15]                *x51
            +coeff[ 16]    *x21        *x51
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]*x11        *x42    
            +coeff[ 19]    *x21*x32        
            +coeff[ 20]*x11    *x31*x41    
            +coeff[ 21]    *x21    *x41*x51
            +coeff[ 22]*x12*x21            
            +coeff[ 23]    *x21*x31*x43    
            +coeff[ 24]    *x22            
            +coeff[ 25]*x11            *x51
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 26]    *x21    *x43    
            +coeff[ 27]    *x21*x31*x42    
            +coeff[ 28]    *x21*x32*x41    
            +coeff[ 29]    *x23*x31        
            +coeff[ 30]*x11*x22*x31        
            +coeff[ 31]*x12*x21    *x41    
            +coeff[ 32]    *x22*x31        
            +coeff[ 33]*x11*x21    *x41    
            +coeff[ 34]*x11    *x32        
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 35]    *x21*x31    *x51
            +coeff[ 36]*x11        *x41*x51
            +coeff[ 37]*x11*x23            
            +coeff[ 38]*x12*x21*x31        
            +coeff[ 39]*x11    *x31*x43    
            +coeff[ 40]*x11*x21            
            +coeff[ 41]*x11*x21*x31        
            +coeff[ 42]    *x22    *x42    
            +coeff[ 43]*x11        *x43    
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 44]    *x22*x31*x41    
            +coeff[ 45]*x11    *x31*x42    
            +coeff[ 46]    *x22    *x41*x51
            +coeff[ 47]    *x23    *x42    
            +coeff[ 48]    *x23*x31*x41    
            +coeff[ 49]*x11*x22    *x42    
            +coeff[ 50]*x11*x22*x31*x41    
            +coeff[ 51]        *x31    *x51
            +coeff[ 52]            *x43    
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 53]            *x44    
            +coeff[ 54]*x12        *x41    
            +coeff[ 55]*x11    *x31    *x51
            +coeff[ 56]        *x33*x41    
            +coeff[ 57]*x11*x21    *x42    
            +coeff[ 58]    *x21    *x42*x51
            +coeff[ 59]*x13                
            +coeff[ 60]    *x24            
            +coeff[ 61]    *x21*x33        
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 62]*x11*x21*x31*x41    
            +coeff[ 63]*x11    *x32*x41    
            +coeff[ 64]    *x21*x31*x41*x51
            +coeff[ 65]        *x31*x44    
            +coeff[ 66]    *x23        *x51
            +coeff[ 67]*x11*x22        *x51
            +coeff[ 68]*x12*x23*x31*x41    
            +coeff[ 69]            *x42    
            +coeff[ 70]        *x31*x41    
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 71]            *x41*x51
            +coeff[ 72]        *x31*x42    
            +coeff[ 73]        *x32*x41    
            +coeff[ 74]    *x21        *x52
            +coeff[ 75]        *x34        
            +coeff[ 76]    *x21    *x44    
            +coeff[ 77]    *x21    *x41*x52
            +coeff[ 78]    *x21*x32*x42    
            +coeff[ 79]*x13        *x41    
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 80]    *x23*x32        
            +coeff[ 81]*x12*x21        *x51
            +coeff[ 82]    *x23    *x41*x51
            +coeff[ 83]*x12*x21    *x42    
            +coeff[ 84]*x11*x24            
            +coeff[ 85]*x11*x22*x32        
            +coeff[ 86]*x11*x22    *x41*x51
            +coeff[ 87]    *x23*x31*x42    
            +coeff[ 88]*x12*x23            
        ;
        v_y_l5p65_ep7                     =v_y_l5p65_ep7                     
            +coeff[ 89]    *x24*x32        
            +coeff[ 90]    *x21*x31*x45    
            +coeff[ 91]*x11*x24    *x41    
            ;
     
        return v_y_l5p65_ep7                     ;
    }
    float p_l5p65_ep7                     (float *x,int m){
        int ncoeff= 41;
        float avdat=  0.2167722E-02;
        float xmin[10]={
            -0.14979E-01,-0.49703E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.49762E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 42]={
            -0.24963633E-03,-0.40513057E-01, 0.22347961E-02,-0.11573684E-01,
            -0.23059223E-02, 0.12291002E-02,-0.47076115E-03,-0.36740927E-02,
            -0.35240056E-02, 0.40743779E-02,-0.22875539E-02,-0.60563576E-06,
            -0.79485303E-03, 0.43659797E-02, 0.24731243E-02,-0.15746537E-03,
            -0.34028362E-03,-0.60437014E-03, 0.27673622E-02, 0.71557431E-03,
            -0.16697308E-02,-0.48438681E-03,-0.80785809E-04, 0.45367351E-03,
             0.39917041E-03, 0.39729004E-03,-0.14977552E-02, 0.70696278E-03,
             0.44553718E-03, 0.61755680E-04,-0.15075933E-03,-0.72318442E-04,
             0.84050116E-04,-0.14993679E-03,-0.77242490E-04, 0.11242334E-03,
            -0.25504935E-03,-0.39547950E-03, 0.72593211E-04, 0.69275426E-04,
            -0.20812583E-03,
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
     
        float v_p_l5p65_ep7                     =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]*x11                
            +coeff[  3]    *x21    *x41    
            +coeff[  4]*x11        *x41    
            +coeff[  5]    *x23*x31        
            +coeff[  6]            *x41    
            +coeff[  7]    *x21*x31        
        ;
        v_p_l5p65_ep7                     =v_p_l5p65_ep7                     
            +coeff[  8]    *x23            
            +coeff[  9]    *x21    *x42    
            +coeff[ 10]*x11*x22            
            +coeff[ 11]    *x23*x31*x41    
            +coeff[ 12]*x11    *x31        
            +coeff[ 13]    *x23    *x41    
            +coeff[ 14]*x11*x22    *x41    
            +coeff[ 15]        *x31        
            +coeff[ 16]    *x21        *x51
        ;
        v_p_l5p65_ep7                     =v_p_l5p65_ep7                     
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]    *x21*x31*x41    
            +coeff[ 19]*x11        *x42    
            +coeff[ 20]    *x21    *x43    
            +coeff[ 21]*x12*x21            
            +coeff[ 22]                *x51
            +coeff[ 23]    *x21*x32        
            +coeff[ 24]    *x21    *x41*x51
            +coeff[ 25]*x11    *x31*x41    
        ;
        v_p_l5p65_ep7                     =v_p_l5p65_ep7                     
            +coeff[ 26]    *x21*x31*x42    
            +coeff[ 27]*x11*x22*x31        
            +coeff[ 28]*x12*x21    *x41    
            +coeff[ 29]*x11*x23    *x41    
            +coeff[ 30]    *x23*x32*x41    
            +coeff[ 31]    *x22            
            +coeff[ 32]            *x42    
            +coeff[ 33]    *x22*x31        
            +coeff[ 34]*x11            *x51
        ;
        v_p_l5p65_ep7                     =v_p_l5p65_ep7                     
            +coeff[ 35]    *x21*x31    *x51
            +coeff[ 36]*x11*x21    *x41    
            +coeff[ 37]    *x21*x32*x41    
            +coeff[ 38]*x11        *x41*x51
            +coeff[ 39]        *x33*x41    
            +coeff[ 40]*x11        *x43    
            ;
     
        return v_p_l5p65_ep7                     ;
    }
    float l_l5p65_ep7                     (float *x,int m){
        int ncoeff= 42;
        float avdat= -0.2380195E-02;
        float xmin[10]={
            -0.14979E-01,-0.49703E-01,-0.14993E-01,-0.32870E-01,-0.49980E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14977E-01, 0.49762E-01, 0.14992E-01, 0.31015E-01, 0.49921E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 43]={
             0.26497513E-02,-0.44797080E-05,-0.17629267E-02,-0.91087734E-02,
             0.51990111E-03,-0.32529803E-02, 0.16754253E-03,-0.91545627E-03,
            -0.37472585E-03, 0.14254193E-03,-0.51799645E-04, 0.29570743E-03,
             0.19676035E-03,-0.15010693E-04, 0.34188331E-04,-0.33049120E-04,
             0.16594240E-03, 0.57139161E-04, 0.10019003E-03,-0.79215075E-04,
            -0.27283279E-04, 0.12929451E-04, 0.59130610E-06, 0.26225958E-04,
             0.26330010E-04,-0.11694447E-03,-0.86257736E-04,-0.37865284E-04,
            -0.64532891E-04, 0.27911809E-04,-0.71494674E-05, 0.68767276E-05,
             0.15586462E-04,-0.27862923E-04,-0.83288705E-05,-0.15204421E-04,
             0.10403538E-04,-0.97014105E-04, 0.39807124E-04,-0.80388985E-04,
             0.14253159E-04, 0.32375112E-04,
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
     
    //                 function
     
        float v_l_l5p65_ep7                     =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]        *x31*x41    
            +coeff[  7]            *x42    
        ;
        v_l_l5p65_ep7                     =v_l_l5p65_ep7                     
            +coeff[  8]*x11*x21            
            +coeff[  9]            *x41*x51
            +coeff[ 10]*x12                
            +coeff[ 11]    *x22    *x41    
            +coeff[ 12]*x11*x21    *x41    
            +coeff[ 13]*x11                
            +coeff[ 14]        *x32        
            +coeff[ 15]                *x52
            +coeff[ 16]    *x22*x31        
        ;
        v_l_l5p65_ep7                     =v_l_l5p65_ep7                     
            +coeff[ 17]    *x22        *x51
            +coeff[ 18]*x11*x21*x31        
            +coeff[ 19]    *x22*x31*x42    
            +coeff[ 20]    *x21    *x41    
            +coeff[ 21]        *x31    *x51
            +coeff[ 22]            *x43    
            +coeff[ 23]*x11*x21        *x51
            +coeff[ 24]*x12        *x41    
            +coeff[ 25]    *x22    *x42    
        ;
        v_l_l5p65_ep7                     =v_l_l5p65_ep7                     
            +coeff[ 26]*x11*x21    *x42    
            +coeff[ 27]    *x22*x32*x41    
            +coeff[ 28]    *x22    *x43    
            +coeff[ 29]    *x21    *x44    
            +coeff[ 30]*x11*x21*x31*x42    
            +coeff[ 31]    *x22*x33*x41    
            +coeff[ 32]*x11*x21*x33*x41    
            +coeff[ 33]    *x23            
            +coeff[ 34]            *x41*x52
        ;
        v_l_l5p65_ep7                     =v_l_l5p65_ep7                     
            +coeff[ 35]*x11*x22            
            +coeff[ 36]*x12    *x31        
            +coeff[ 37]    *x22*x31*x41    
            +coeff[ 38]*x11*x23            
            +coeff[ 39]*x11*x21*x31*x41    
            +coeff[ 40]    *x21*x33*x41    
            +coeff[ 41]*x12*x24            
            ;
     
        return v_l_l5p65_ep7                     ;
    }
    float x_l5p65_fp                      (float *x,int m){
        int ncoeff= 46;
        float avdat=  0.1947775E-01;
        float xmin[10]={
            -0.14985E-01,-0.46325E-01,-0.14986E-01,-0.29661E-01,-0.47283E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14999E-01, 0.47334E-01, 0.14991E-01, 0.30861E-01, 0.49993E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 47]={
            -0.27709529E-02,-0.36931135E-01, 0.64789128E+00,-0.52326437E-01,
             0.44880353E-01, 0.13993715E-01, 0.36335804E-01, 0.12421987E-01,
            -0.14167706E-01, 0.69180001E-02, 0.58017350E-02, 0.11374928E-01,
            -0.99262595E-02,-0.30775055E-01,-0.76341158E-03, 0.22080487E-02,
            -0.19657579E-02,-0.87501446E-03, 0.23985868E-02, 0.36033113E-02,
            -0.43195356E-02,-0.27793224E-02, 0.47715893E-02, 0.16336216E-01,
             0.63351974E-04, 0.51856111E-02,-0.19006683E-01, 0.10455914E-02,
             0.17627977E-02, 0.10334588E-03, 0.29285159E-02, 0.27075813E-02,
            -0.14674159E-02,-0.84977923E-02,-0.38623214E-02, 0.66705104E-02,
             0.95292146E-03,-0.10859545E-02, 0.12588683E-02, 0.25814422E-02,
            -0.31323261E-02,-0.45394008E-02, 0.37351504E-02, 0.24903342E-02,
            -0.24440812E-02, 0.29540169E-02,
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
     
        float v_x_l5p65_fp                      =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]                *x51
            +coeff[  3]                *x52
            +coeff[  4]    *x21        *x51
            +coeff[  5]    *x22            
            +coeff[  6]    *x21    *x41    
            +coeff[  7]    *x21*x31        
        ;
        v_x_l5p65_fp                      =v_x_l5p65_fp                      
            +coeff[  8]            *x42    
            +coeff[  9]*x11        *x41    
            +coeff[ 10]                *x53
            +coeff[ 11]    *x23            
            +coeff[ 12]    *x21*x31*x41    
            +coeff[ 13]    *x21    *x42    
            +coeff[ 14]*x12*x21    *x41*x51
            +coeff[ 15]            *x41    
            +coeff[ 16]*x11            *x51
        ;
        v_x_l5p65_fp                      =v_x_l5p65_fp                      
            +coeff[ 17]*x11*x21            
            +coeff[ 18]*x11    *x31        
            +coeff[ 19]            *x41*x51
            +coeff[ 20]        *x31*x41    
            +coeff[ 21]    *x21        *x52
            +coeff[ 22]*x11*x22            
            +coeff[ 23]    *x21    *x41*x51
            +coeff[ 24]        *x33        
            +coeff[ 25]*x11*x23            
        ;
        v_x_l5p65_fp                      =v_x_l5p65_fp                      
            +coeff[ 26]    *x23    *x41    
            +coeff[ 27]        *x31    *x51
            +coeff[ 28]*x12*x21            
            +coeff[ 29]*x11        *x42    
            +coeff[ 30]    *x22        *x51
            +coeff[ 31]    *x21*x31    *x51
            +coeff[ 32]    *x21*x32        
            +coeff[ 33]*x11*x22    *x41    
            +coeff[ 34]    *x21    *x42*x51
        ;
        v_x_l5p65_fp                      =v_x_l5p65_fp                      
            +coeff[ 35]    *x22    *x42    
            +coeff[ 36]        *x31        
            +coeff[ 37]*x11    *x31*x41    
            +coeff[ 38]            *x42*x51
            +coeff[ 39]    *x23        *x51
            +coeff[ 40]    *x23*x31        
            +coeff[ 41]    *x23    *x41*x51
            +coeff[ 42]    *x22    *x43    
            +coeff[ 43]*x13*x22        *x51
        ;
        v_x_l5p65_fp                      =v_x_l5p65_fp                      
            +coeff[ 44]*x13*x22*x31        
            +coeff[ 45]*x13*x23    *x41    
            ;
     
        return v_x_l5p65_fp                      ;
    }
    float t_l5p65_fp                      (float *x,int m){
        int ncoeff= 50;
        float avdat=  0.3970148E-02;
        float xmin[10]={
            -0.14985E-01,-0.46325E-01,-0.14986E-01,-0.29661E-01,-0.47283E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14999E-01, 0.47334E-01, 0.14991E-01, 0.30861E-01, 0.49993E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 51]={
            -0.48515787E-04,-0.34947386E-02,-0.18967204E-01,-0.37465707E-04,
             0.11223308E+00, 0.53727585E-02,-0.10495048E-01, 0.10492973E-02,
            -0.22506521E-02, 0.95202260E-04,-0.64148766E-03,-0.25228390E-02,
             0.21423001E-02,-0.75959263E-03,-0.84019493E-05,-0.71415392E-03,
            -0.35064289E-03, 0.51272241E-03, 0.10875443E-02, 0.14177888E-02,
            -0.24923580E-03,-0.31485662E-03, 0.16391608E-03,-0.36667599E-03,
            -0.85410202E-03, 0.44974277E-03, 0.39484055E-03, 0.31834125E-03,
             0.51962623E-04, 0.94150036E-03,-0.80100574E-04,-0.65145927E-04,
            -0.23280593E-03, 0.17712478E-03, 0.23363048E-03, 0.12392229E-03,
             0.13017576E-03,-0.68455504E-03,-0.37651046E-03, 0.44636181E-03,
             0.51613996E-03,-0.69147092E-03, 0.98105124E-03,-0.59521606E-03,
             0.92861068E-04, 0.76253215E-04, 0.53489552E-04, 0.13542773E-03,
             0.95588744E-04, 0.11758942E-03,
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
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_t_l5p65_fp                      =avdat
            +coeff[  0]                    
            +coeff[  1]*x11                
            +coeff[  2]    *x21            
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x21        *x51
            +coeff[  6]                *x52
            +coeff[  7]    *x22            
        ;
        v_t_l5p65_fp                      =v_t_l5p65_fp                      
            +coeff[  8]            *x42    
            +coeff[  9]*x11*x21            
            +coeff[ 10]    *x21    *x41    
            +coeff[ 11]    *x21    *x42    
            +coeff[ 12]    *x21    *x41*x51
            +coeff[ 13]    *x21        *x52
            +coeff[ 14]        *x31        
            +coeff[ 15]        *x31*x41    
            +coeff[ 16]*x11            *x51
        ;
        v_t_l5p65_fp                      =v_t_l5p65_fp                      
            +coeff[ 17]            *x41*x51
            +coeff[ 18]                *x53
            +coeff[ 19]    *x22    *x42    
            +coeff[ 20]    *x21*x31        
            +coeff[ 21]*x11        *x41    
            +coeff[ 22]        *x31    *x51
            +coeff[ 23]*x11*x22            
            +coeff[ 24]    *x22    *x41    
            +coeff[ 25]*x11        *x42    
        ;
        v_t_l5p65_fp                      =v_t_l5p65_fp                      
            +coeff[ 26]    *x21*x31    *x51
            +coeff[ 27]            *x42*x51
            +coeff[ 28]*x12*x22            
            +coeff[ 29]*x11*x23            
            +coeff[ 30]*x12                
            +coeff[ 31]*x11    *x31        
            +coeff[ 32]    *x22*x31        
            +coeff[ 33]*x11*x21    *x41    
            +coeff[ 34]    *x22        *x51
        ;
        v_t_l5p65_fp                      =v_t_l5p65_fp                      
            +coeff[ 35]*x11            *x52
            +coeff[ 36]    *x23*x31        
            +coeff[ 37]    *x23    *x41    
            +coeff[ 38]    *x21    *x43    
            +coeff[ 39]*x11*x22        *x51
            +coeff[ 40]    *x23        *x51
            +coeff[ 41]    *x21    *x42*x51
            +coeff[ 42]    *x22    *x43    
            +coeff[ 43]    *x21*x31*x41*x53
        ;
        v_t_l5p65_fp                      =v_t_l5p65_fp                      
            +coeff[ 44]    *x23            
            +coeff[ 45]*x11*x21*x31        
            +coeff[ 46]*x11    *x31*x41    
            +coeff[ 47]        *x31*x42    
            +coeff[ 48]*x11    *x31    *x51
            +coeff[ 49]*x11        *x41*x51
            ;
     
        return v_t_l5p65_fp                      ;
    }
    float y_l5p65_fp                      (float *x,int m){
        int ncoeff= 97;
        float avdat= -0.5008442E-02;
        float xmin[10]={
            -0.14985E-01,-0.46325E-01,-0.14986E-01,-0.29661E-01,-0.47283E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14999E-01, 0.47334E-01, 0.14991E-01, 0.30861E-01, 0.49993E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 98]={
             0.53559770E-02,-0.60466979E-01,-0.15070972E-02,-0.26900083E-03,
             0.99107279E-02, 0.10631410E-01,-0.15593206E-01,-0.16860241E-02,
            -0.38495716E-02,-0.67891818E-02,-0.15531536E-02,-0.46003638E-02,
            -0.21233500E-02, 0.10782406E-01,-0.67840638E-02,-0.36186184E-03,
            -0.20896179E-03, 0.67643044E-02, 0.43651098E-02,-0.12687203E-04,
             0.49674469E-02, 0.26094422E-02,-0.17682809E-04, 0.10397878E-02,
            -0.21423427E-02, 0.65353871E-02, 0.43316032E-02, 0.75211580E-03,
             0.30660632E-02,-0.81703917E-03, 0.26946270E-03,-0.19098996E-02,
             0.23059987E-02, 0.13683953E-02,-0.48528300E-02, 0.94547309E-03,
             0.16364052E-02, 0.28394128E-02,-0.98092816E-04, 0.10750551E-03,
             0.10195917E-02, 0.95776521E-03, 0.42180327E-03, 0.94787637E-03,
            -0.29768436E-02, 0.24778510E-02,-0.37695499E-03, 0.28395743E-03,
             0.21018342E-02,-0.15381481E-02,-0.43388433E-02,-0.33192639E-03,
            -0.18804696E-04, 0.14604877E-03, 0.69679314E-03, 0.31408475E-03,
            -0.11797792E-02,-0.25979034E-02,-0.10005489E-02,-0.16534941E-03,
            -0.26106508E-03,-0.23205215E-02,-0.95677719E-03,-0.43132188E-03,
             0.48710470E-03, 0.59461861E-03,-0.16559195E-02,-0.55933124E-02,
             0.40915303E-03, 0.86143549E-03, 0.70776790E-03,-0.87251735E-03,
             0.26513196E-02,-0.28148065E-05,-0.28972217E-03, 0.37314842E-03,
            -0.15072485E-03, 0.17236728E-03,-0.12575730E-02,-0.13081463E-03,
            -0.53717225E-03, 0.14403887E-03, 0.11128023E-02,-0.33211455E-03,
            -0.95911388E-03, 0.26326301E-03, 0.13714128E-02,-0.49956288E-03,
            -0.12639879E-02, 0.16977638E-03,-0.11864791E-03,-0.45523015E-03,
            -0.10141250E-02, 0.66890521E-03,-0.11067818E-02, 0.18777716E-02,
            -0.72117167E-03,
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
     
        float v_y_l5p65_fp                      =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]*x11                
            +coeff[  4]                *x51
            +coeff[  5]    *x21    *x41    
            +coeff[  6]    *x22            
            +coeff[  7]*x11        *x41    
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[  8]            *x41*x51
            +coeff[  9]*x11*x21            
            +coeff[ 10]    *x21        *x51
            +coeff[ 11]        *x31    *x51
            +coeff[ 12]                *x52
            +coeff[ 13]    *x22    *x41    
            +coeff[ 14]    *x21    *x41*x51
            +coeff[ 15]            *x44    
            +coeff[ 16]        *x31*x43    
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 17]            *x41*x52
            +coeff[ 18]    *x22*x31        
            +coeff[ 19]        *x33        
            +coeff[ 20]*x11*x21    *x41    
            +coeff[ 21]    *x22        *x51
            +coeff[ 22]        *x33    *x52
            +coeff[ 23]    *x22    *x41*x53
            +coeff[ 24]        *x31        
            +coeff[ 25]            *x42    
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 26]        *x31*x41    
            +coeff[ 27]        *x32        
            +coeff[ 28]    *x21    *x42    
            +coeff[ 29]*x12                
            +coeff[ 30]*x11            *x51
            +coeff[ 31]            *x42*x51
            +coeff[ 32]    *x23            
            +coeff[ 33]*x11*x21*x31        
            +coeff[ 34]    *x22    *x42    
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 35]*x11*x21        *x51
            +coeff[ 36]    *x21        *x52
            +coeff[ 37]        *x31    *x52
            +coeff[ 38]        *x33*x41*x51
            +coeff[ 39]*x11    *x31        
            +coeff[ 40]            *x43    
            +coeff[ 41]*x11        *x42    
            +coeff[ 42]*x12        *x41    
            +coeff[ 43]*x11        *x41*x51
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 44]    *x22*x31*x41    
            +coeff[ 45]    *x24            
            +coeff[ 46]*x11            *x52
            +coeff[ 47]    *x21*x31*x41*x51
            +coeff[ 48]*x11*x23            
            +coeff[ 49]    *x23        *x51
            +coeff[ 50]            *x41*x53
            +coeff[ 51]    *x24    *x41*x51
            +coeff[ 52]    *x22*x31    *x53
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 53]    *x21*x31        
            +coeff[ 54]    *x21*x31*x41    
            +coeff[ 55]*x11    *x31*x41    
            +coeff[ 56]        *x31*x41*x51
            +coeff[ 57]    *x21    *x43    
            +coeff[ 58]    *x21*x31*x42    
            +coeff[ 59]*x12*x21            
            +coeff[ 60]*x11    *x31    *x51
            +coeff[ 61]*x11*x21    *x42    
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 62]    *x21    *x42*x51
            +coeff[ 63]    *x23*x31        
            +coeff[ 64]                *x53
            +coeff[ 65]*x11*x22    *x41    
            +coeff[ 66]*x11*x21*x31*x41    
            +coeff[ 67]    *x22    *x41*x51
            +coeff[ 68]*x12*x21    *x41    
            +coeff[ 69]    *x21    *x41*x52
            +coeff[ 70]*x12*x22            
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 71]        *x31    *x53
            +coeff[ 72]    *x22    *x41*x52
            +coeff[ 73]        *x31*x42    
            +coeff[ 74]        *x32*x41    
            +coeff[ 75]*x11*x22            
            +coeff[ 76]        *x32    *x51
            +coeff[ 77]*x12    *x31        
            +coeff[ 78]            *x43*x51
            +coeff[ 79]    *x23    *x41    
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 80]    *x22*x32        
            +coeff[ 81]*x12            *x51
            +coeff[ 82]            *x42*x52
            +coeff[ 83]*x11*x21*x32        
            +coeff[ 84]    *x22*x31    *x51
            +coeff[ 85]    *x22        *x52
            +coeff[ 86]    *x21    *x43*x51
            +coeff[ 87]*x11*x22    *x42    
            +coeff[ 88]    *x22    *x42*x51
        ;
        v_y_l5p65_fp                      =v_y_l5p65_fp                      
            +coeff[ 89]*x13*x21            
            +coeff[ 90]    *x21*x34        
            +coeff[ 91]    *x21        *x53
            +coeff[ 92]*x11*x23    *x41    
            +coeff[ 93]    *x23    *x41*x51
            +coeff[ 94]    *x21    *x42*x52
            +coeff[ 95]    *x23    *x43    
            +coeff[ 96]*x11*x22    *x41*x51
            ;
     
        return v_y_l5p65_fp                      ;
    }
    float p_l5p65_fp                      (float *x,int m){
        int ncoeff= 87;
        float avdat= -0.1520557E-02;
        float xmin[10]={
            -0.14985E-01,-0.46325E-01,-0.14986E-01,-0.29661E-01,-0.47283E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14999E-01, 0.47334E-01, 0.14991E-01, 0.30861E-01, 0.49993E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 88]={
             0.29462771E-02,-0.11356351E-02, 0.63721994E-02,-0.33295032E-01,
             0.10168304E-01,-0.14792960E-01, 0.17743715E-02, 0.20147569E-01,
             0.41142623E-02, 0.65080901E-02,-0.30774139E-02,-0.21355079E-01,
            -0.63022673E-02, 0.34092199E-02, 0.15788957E-02,-0.28835537E-02,
             0.18032772E-01, 0.55835792E-02,-0.35129692E-02, 0.34800125E-03,
            -0.66599889E-04,-0.10675369E-03, 0.41781881E-04, 0.72359829E-03,
            -0.36186269E-02, 0.42951899E-02, 0.22520442E-02,-0.78234728E-03,
             0.25475996E-02, 0.29325848E-02,-0.84849721E-03, 0.22553196E-02,
             0.30915102E-03, 0.69307443E-03, 0.39947904E-02, 0.12600687E-02,
             0.87676052E-03, 0.23325644E-02, 0.20592909E-03,-0.24808562E-03,
            -0.24779711E-03,-0.21727887E-03, 0.12066135E-03,-0.41297777E-03,
             0.67068095E-03,-0.17512335E-02, 0.13349081E-02,-0.30427927E-02,
             0.60017529E-03, 0.15974359E-02,-0.13650217E-02,-0.15375924E-02,
             0.70994749E-03,-0.11603239E-02,-0.80912729E-03, 0.23356857E-03,
             0.16527242E-03,-0.56153379E-03, 0.62041648E-03,-0.11994361E-03,
            -0.15759388E-03,-0.18512390E-02,-0.12955125E-03,-0.29003879E-03,
             0.64006395E-03, 0.60149294E-03, 0.39502268E-03, 0.12474407E-03,
             0.55455257E-04,-0.52648474E-03, 0.73086259E-04,-0.39059023E-03,
            -0.14894110E-02,-0.80997946E-04, 0.12789910E-03,-0.37858065E-03,
             0.38584636E-03, 0.39222973E-03, 0.21658614E-03, 0.14992419E-03,
             0.60494745E-03,-0.20676470E-03, 0.13140215E-03,-0.74698770E-03,
             0.19977642E-03,-0.66799950E-03,-0.36775635E-03,
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
        float x34 = x33*x3;
        float x35 = x34*x3;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x44 = x43*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_p_l5p65_fp                      =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]        *x31        
            +coeff[  3]            *x41    
            +coeff[  4]                *x51
            +coeff[  5]    *x22            
            +coeff[  6]    *x21*x31        
            +coeff[  7]    *x21    *x41    
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[  8]        *x31*x41    
            +coeff[  9]            *x42    
            +coeff[ 10]    *x21        *x51
            +coeff[ 11]            *x41*x51
            +coeff[ 12]*x11*x21            
            +coeff[ 13]    *x23            
            +coeff[ 14]                *x52
            +coeff[ 15]*x11        *x41    
            +coeff[ 16]    *x22    *x41    
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 17]    *x21    *x42    
            +coeff[ 18]    *x22        *x51
            +coeff[ 19]    *x22*x31    *x51
            +coeff[ 20]        *x33    *x51
            +coeff[ 21]    *x24*x31        
            +coeff[ 22]        *x35    *x51
            +coeff[ 23]        *x32        
            +coeff[ 24]        *x31    *x51
            +coeff[ 25]    *x22*x31        
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 26]            *x43    
            +coeff[ 27]*x12                
            +coeff[ 28]*x11*x21    *x41    
            +coeff[ 29]            *x41*x52
            +coeff[ 30]*x11*x21        *x51
            +coeff[ 31]    *x21*x31*x41    
            +coeff[ 32]*x11            *x51
            +coeff[ 33]*x11*x22            
            +coeff[ 34]    *x24            
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 35]*x11*x21*x31        
            +coeff[ 36]        *x31    *x52
            +coeff[ 37]*x11*x23            
            +coeff[ 38]*x11        *x44    
            +coeff[ 39]*x11                
            +coeff[ 40]*x11    *x31        
            +coeff[ 41]        *x32*x41    
            +coeff[ 42]        *x31*x41*x51
            +coeff[ 43]            *x42*x51
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 44]*x11    *x31*x41    
            +coeff[ 45]    *x22*x31*x41    
            +coeff[ 46]*x11        *x42    
            +coeff[ 47]    *x21    *x43    
            +coeff[ 48]*x12        *x41    
            +coeff[ 49]*x11*x22    *x41    
            +coeff[ 50]*x11*x21*x31*x41    
            +coeff[ 51]*x11*x21    *x42    
            +coeff[ 52]*x11*x21    *x41*x51
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 53]            *x41*x53
            +coeff[ 54]    *x23*x31*x42    
            +coeff[ 55]        *x31*x42    
            +coeff[ 56]    *x21        *x52
            +coeff[ 57]    *x22*x32        
            +coeff[ 58]    *x22    *x41*x51
            +coeff[ 59]*x12*x21            
            +coeff[ 60]*x11            *x52
            +coeff[ 61]    *x23    *x42    
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 62]*x12            *x51
            +coeff[ 63]        *x31    *x53
            +coeff[ 64]    *x21    *x44    
            +coeff[ 65]*x12*x22            
            +coeff[ 66]*x12*x21    *x41    
            +coeff[ 67]    *x21    *x41*x51
            +coeff[ 68]        *x32    *x51
            +coeff[ 69]    *x23    *x41    
            +coeff[ 70]*x11    *x32        
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 71]    *x22    *x42    
            +coeff[ 72]    *x21*x31*x42    
            +coeff[ 73]*x11    *x31    *x51
            +coeff[ 74]*x11        *x41*x51
            +coeff[ 75]            *x44    
            +coeff[ 76]    *x21*x31*x41*x51
            +coeff[ 77]    *x21    *x42*x51
            +coeff[ 78]    *x22        *x52
            +coeff[ 79]*x12    *x31        
        ;
        v_p_l5p65_fp                      =v_p_l5p65_fp                      
            +coeff[ 80]    *x24    *x41    
            +coeff[ 81]*x11*x21*x32        
            +coeff[ 82]    *x23*x32        
            +coeff[ 83]    *x23*x31*x41    
            +coeff[ 84]*x11*x21*x31    *x51
            +coeff[ 85]    *x22    *x42*x51
            +coeff[ 86]    *x25*x31        
            ;
     
        return v_p_l5p65_fp                      ;
    }
    float l_l5p65_fp                      (float *x,int m){
        int ncoeff= 76;
        float avdat= -0.1693699E-01;
        float xmin[10]={
            -0.14985E-01,-0.46325E-01,-0.14986E-01,-0.29661E-01,-0.47283E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.14999E-01, 0.47334E-01, 0.14991E-01, 0.30861E-01, 0.49993E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 77]={
             0.13271747E-01,-0.25367579E+00,-0.34108263E-01, 0.39095115E-01,
            -0.24554597E-01,-0.93674853E-01,-0.26349677E-01,-0.20245083E-01,
             0.10566272E-01, 0.64536566E-02, 0.58721058E-03,-0.32167610E-01,
            -0.94843162E-02, 0.53921952E-02,-0.71488172E-02,-0.24338955E-01,
            -0.22585683E-01, 0.22924401E-01, 0.38619366E-01,-0.15101566E-01,
             0.31100863E-02,-0.55013220E-02, 0.48788763E-02,-0.39644591E-02,
             0.37617225E-01, 0.22709930E-01,-0.28535852E-02,-0.14912982E-02,
             0.39777192E-02,-0.36925646E-02, 0.57903756E-02, 0.69775353E-02,
             0.14193852E-03,-0.57186722E-03,-0.24311752E-02, 0.41181324E-02,
             0.58547622E-02, 0.42322930E-02, 0.11146805E-02, 0.12802882E-02,
             0.48139811E-03,-0.70919702E-03, 0.69989177E-03,-0.44829953E-02,
            -0.11468450E-01,-0.10274542E-01,-0.33313967E-02, 0.12139782E-02,
             0.43453658E-02, 0.18492268E-02, 0.30075316E-03,-0.10189854E-02,
            -0.79377217E-03, 0.14667237E-02,-0.46030077E-03, 0.38972305E-03,
             0.95108384E-03, 0.28300061E-03,-0.30070750E-03, 0.25450671E-02,
            -0.33893958E-02,-0.71624084E-03,-0.59093256E-03, 0.11663972E-02,
            -0.12748616E-02, 0.59790979E-03,-0.70125051E-03,-0.13809204E-02,
            -0.62569370E-03,-0.13483760E-02, 0.94966916E-03, 0.30611323E-02,
             0.56620984E-03,-0.52207697E-03,-0.78260200E-03,-0.10114067E-02,
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
     
        float v_l_l5p65_fp                      =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]                *x51
            +coeff[  3]*x11                
            +coeff[  4]    *x22            
            +coeff[  5]    *x21    *x41    
            +coeff[  6]                *x52
            +coeff[  7]*x11        *x41    
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[  8]    *x23*x31        
            +coeff[  9]    *x21    *x44    
            +coeff[ 10]            *x41    
            +coeff[ 11]    *x21*x31        
            +coeff[ 12]            *x42    
            +coeff[ 13]*x11*x21            
            +coeff[ 14]*x11    *x31        
            +coeff[ 15]    *x23            
            +coeff[ 16]    *x22    *x41    
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[ 17]    *x21*x31*x41    
            +coeff[ 18]    *x21    *x42    
            +coeff[ 19]*x11*x22            
            +coeff[ 20]            *x41*x51
            +coeff[ 21]    *x22*x31        
            +coeff[ 22]                *x53
            +coeff[ 23]*x12*x21            
            +coeff[ 24]    *x23    *x41    
            +coeff[ 25]*x11*x22    *x41    
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[ 26]    *x21        *x51
            +coeff[ 27]*x11            *x51
            +coeff[ 28]    *x21*x32        
            +coeff[ 29]*x11*x21    *x41    
            +coeff[ 30]*x11        *x42    
            +coeff[ 31]*x11*x22*x31        
            +coeff[ 32]*x11    *x33*x41    
            +coeff[ 33]*x12                
            +coeff[ 34]    *x21        *x52
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[ 35]*x11    *x31*x41    
            +coeff[ 36]    *x22    *x42    
            +coeff[ 37]*x12*x21    *x41    
            +coeff[ 38]    *x21*x31    *x51
            +coeff[ 39]        *x31*x41*x51
            +coeff[ 40]            *x41*x52
            +coeff[ 41]*x11*x21*x31        
            +coeff[ 42]*x11    *x32        
            +coeff[ 43]    *x24            
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[ 44]    *x21*x31*x42    
            +coeff[ 45]    *x21    *x43    
            +coeff[ 46]*x11*x23            
            +coeff[ 47]*x12*x21*x31        
            +coeff[ 48]    *x24    *x41    
            +coeff[ 49]    *x21    *x41*x53
            +coeff[ 50]*x11        *x43*x51
            +coeff[ 51]    *x23*x32*x41    
            +coeff[ 52]        *x31*x41    
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[ 53]            *x43    
            +coeff[ 54]    *x22        *x51
            +coeff[ 55]*x11    *x31    *x51
            +coeff[ 56]*x11        *x41*x51
            +coeff[ 57]*x11            *x52
            +coeff[ 58]*x13                
            +coeff[ 59]    *x22*x31*x41    
            +coeff[ 60]    *x21*x32*x41    
            +coeff[ 61]    *x23        *x51
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[ 62]    *x21*x31*x41*x51
            +coeff[ 63]    *x21    *x42*x51
            +coeff[ 64]    *x21    *x41*x52
            +coeff[ 65]    *x21        *x53
            +coeff[ 66]                *x54
            +coeff[ 67]*x11        *x42*x51
            +coeff[ 68]*x12*x22            
            +coeff[ 69]    *x23*x31*x41    
            +coeff[ 70]        *x33*x42    
        ;
        v_l_l5p65_fp                      =v_l_l5p65_fp                      
            +coeff[ 71]    *x21*x31*x43    
            +coeff[ 72]*x12    *x32    *x51
            +coeff[ 73]        *x34    *x52
            +coeff[ 74]        *x33*x41*x52
            +coeff[ 75]*x11    *x33*x41*x51
            ;
     
        return v_l_l5p65_fp                      ;
    }
}
