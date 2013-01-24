#include "Fwd_l5p77_484816R00.h"

namespace S484816R00
{
    float x_l5p77_fp           (float *x,int m){
        int ncoeff= 24;
        float avdat=  0.1389984E-01;
        float xmin[10]={
             0.00000E+00,-0.42303E-01, 0.00000E+00,-0.21263E-01,-0.46243E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.43132E-01, 0.00000E+00, 0.27251E-01, 0.49919E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 25]={
             0.18673921E-01, 0.64265275E+00,-0.51232368E-01, 0.42292271E-01,
             0.10371306E-01, 0.21321056E-01, 0.55447067E-02, 0.49898922E-02,
            -0.18915461E-01,-0.54419767E-02,-0.94244666E-02, 0.54466925E-02,
             0.11432515E-01,-0.13507863E-01, 0.25003385E-02,-0.29455372E-02,
             0.42883833E-02, 0.17648686E-02, 0.46923584E-02, 0.88578672E-03,
             0.59361779E-03, 0.12125837E-02,-0.27403352E-02,-0.40278230E-02,
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
        float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
        float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
    //  set up monomials   functions
        float x21 = x2;
        float x22 = x21*x2;
        float x23 = x22*x2;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_x_l5p77_fp           =avdat
            +coeff[  0]                    
            +coeff[  1]                *x51
            +coeff[  2]                *x52
            +coeff[  3]    *x21        *x51
            +coeff[  4]    *x22            
            +coeff[  5]    *x21    *x41    
            +coeff[  6]    *x21            
            +coeff[  7]    *x23            
        ;
        v_x_l5p77_fp           =v_x_l5p77_fp           
            +coeff[  8]    *x21    *x42    
            +coeff[  9]            *x41    
            +coeff[ 10]            *x42    
            +coeff[ 11]                *x53
            +coeff[ 12]    *x21    *x41*x51
            +coeff[ 13]    *x23    *x41    
            +coeff[ 14]            *x41*x51
            +coeff[ 15]    *x21        *x52
            +coeff[ 16]    *x22    *x41    
        ;
        v_x_l5p77_fp           =v_x_l5p77_fp           
            +coeff[ 17]    *x23        *x51
            +coeff[ 18]    *x22    *x42    
            +coeff[ 19]    *x22        *x51
            +coeff[ 20]            *x42*x51
            +coeff[ 21]            *x43    
            +coeff[ 22]    *x21    *x42*x51
            +coeff[ 23]    *x23    *x41*x51
            ;
     
        return v_x_l5p77_fp           ;
    }
    float t_l5p77_fp           (float *x,int m){
        int ncoeff= 46;
        float avdat=  0.1639012E-02;
        float xmin[10]={
             0.00000E+00,-0.42303E-01, 0.00000E+00,-0.21263E-01,-0.46243E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.43132E-01, 0.00000E+00, 0.27251E-01, 0.49919E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 47]={
             0.32005305E-02,-0.17290648E-01,-0.51969930E-03, 0.11120266E+00,
             0.52461312E-02,-0.10204817E-01, 0.10785257E-02,-0.15935412E-02,
            -0.95726678E-03, 0.16611993E-02,-0.71489939E-03, 0.39865795E-03,
            -0.15939062E-02, 0.10953646E-02, 0.10139035E-02, 0.15052613E-03,
             0.18623174E-03,-0.77362514E-04,-0.38485092E-03, 0.28701042E-03,
            -0.78350172E-03,-0.19991436E-03,-0.10782171E-03,-0.64244773E-03,
            -0.33310783E-03,-0.36752378E-03,-0.24002579E-03,-0.39741182E-03,
             0.33275736E-03,-0.17199728E-03,-0.21429495E-03,-0.86694206E-04,
            -0.41560786E-04,-0.14364668E-03, 0.65528904E-03, 0.37719338E-03,
             0.33025545E-03, 0.32328401E-03, 0.71650073E-04,-0.10830008E-03,
             0.81513412E-04, 0.98616176E-04, 0.60858314E-04, 0.71570154E-04,
             0.95820025E-04,-0.65658372E-04,
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
        float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
        float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
    //  set up monomials   functions
        float x21 = x2;
        float x22 = x21*x2;
        float x23 = x22*x2;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_t_l5p77_fp           =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]            *x41    
            +coeff[  3]                *x51
            +coeff[  4]    *x21        *x51
            +coeff[  5]                *x52
            +coeff[  6]    *x22            
            +coeff[  7]            *x42    
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[  8]    *x21    *x41    
            +coeff[  9]    *x21    *x41*x51
            +coeff[ 10]    *x21        *x52
            +coeff[ 11]            *x41*x51
            +coeff[ 12]    *x21    *x42    
            +coeff[ 13]                *x53
            +coeff[ 14]    *x22    *x42    
            +coeff[ 15]            *x43    
            +coeff[ 16]            *x42*x51
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 17]    *x23        *x51
            +coeff[ 18]    *x21    *x42*x51
            +coeff[ 19]    *x22    *x43    
            +coeff[ 20]    *x23    *x41*x51
            +coeff[ 21]    *x23            
            +coeff[ 22]            *x41*x52
            +coeff[ 23]    *x23    *x41    
            +coeff[ 24]    *x21    *x43    
            +coeff[ 25]    *x22        *x52
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 26]    *x23        *x52
            +coeff[ 27]    *x22    *x41*x52
            +coeff[ 28]    *x21    *x42*x52
            +coeff[ 29]            *x42*x53
            +coeff[ 30]    *x22        *x51
            +coeff[ 31]    *x22    *x41*x51
            +coeff[ 32]    *x21    *x41*x52
            +coeff[ 33]            *x42*x52
            +coeff[ 34]    *x22    *x42*x51
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 35]    *x22    *x43*x51
            +coeff[ 36]    *x22    *x42*x52
            +coeff[ 37]    *x23        *x53
            +coeff[ 38]    *x22    *x41    
            +coeff[ 39]    *x21        *x53
            +coeff[ 40]            *x41*x53
            +coeff[ 41]            *x43*x52
            +coeff[ 42]    *x21    *x43*x51
            +coeff[ 43]    *x22        *x53
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 44]    *x23    *x43    
            +coeff[ 45]            *x43*x53
            ;
     
        return v_t_l5p77_fp           ;
    }
    float y_l5p77_fp           (float *x,int m){
        int ncoeff= 85;
        float avdat= -0.7036306E-02;
        float xmin[10]={
             0.00000E+00,-0.42303E-01, 0.00000E+00,-0.21263E-01,-0.46243E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.43132E-01, 0.00000E+00, 0.27251E-01, 0.49919E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 86]={
            -0.25976351E-02,-0.46636157E-01, 0.32649206E-02, 0.10384003E-01,
             0.41599432E-02, 0.70719370E-02,-0.12898093E-01,-0.34869441E-02,
            -0.30341588E-02,-0.15614951E-02, 0.52378988E-02,-0.54033301E-02,
             0.52671782E-02, 0.29310980E-02,-0.31993433E-02,-0.34611395E-02,
            -0.25331436E-03, 0.16931664E-02, 0.19412530E-02,-0.37409014E-02,
             0.19221196E-04,-0.13711492E-02, 0.16081283E-02, 0.15581662E-02,
            -0.12152784E-02,-0.18746992E-04, 0.31840824E-04,-0.38285396E-03,
            -0.85020787E-03, 0.64956897E-03, 0.64870476E-03,-0.96950046E-03,
             0.18451497E-02, 0.41585714E-04,-0.10802536E-02,-0.24093266E-04,
            -0.16158704E-02,-0.53068309E-03, 0.17874944E-02,-0.10626921E-02,
             0.27875745E-03,-0.10422415E-02,-0.13930172E-02,-0.54544862E-03,
             0.30986881E-02, 0.23434457E-03, 0.11129740E-02, 0.11219335E-02,
            -0.95618120E-03, 0.25749995E-03, 0.16186264E-02, 0.16464748E-02,
            -0.80029812E-03,-0.55160741E-04, 0.73016237E-03, 0.31643172E-04,
             0.20880897E-03,-0.34209265E-03, 0.10210322E-02,-0.12468221E-02,
             0.19542927E-03,-0.16645189E-03, 0.99442329E-03, 0.96906541E-03,
            -0.82235405E-04, 0.10215097E-03, 0.43541388E-03,-0.54178579E-03,
             0.93493232E-04, 0.53810614E-04, 0.73327485E-03, 0.96230692E-03,
            -0.15166582E-02, 0.23201606E-04, 0.73647866E-03,-0.45456507E-03,
            -0.76747511E-03,-0.11291754E-02, 0.67651720E-03, 0.39304176E-03,
            -0.68134803E-03, 0.37640581E-03, 0.21960004E-03, 0.25630964E-03,
             0.66998487E-04,
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
        float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
        float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
    //  set up monomials   functions
        float x21 = x2;
        float x22 = x21*x2;
        float x23 = x22*x2;
        float x24 = x23*x2;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x44 = x43*x4;
        float x45 = x44*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_y_l5p77_fp           =avdat
            +coeff[  0]                    
            +coeff[  1]            *x41    
            +coeff[  2]    *x21            
            +coeff[  3]                *x51
            +coeff[  4]            *x42    
            +coeff[  5]    *x21    *x41    
            +coeff[  6]    *x22            
            +coeff[  7]            *x41*x51
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[  8]    *x21        *x51
            +coeff[  9]                *x52
            +coeff[ 10]    *x22    *x41    
            +coeff[ 11]    *x21    *x41*x51
            +coeff[ 12]            *x41*x52
            +coeff[ 13]    *x22        *x51
            +coeff[ 14]    *x22    *x42    
            +coeff[ 15]            *x41*x53
            +coeff[ 16]    *x23        *x52
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 17]    *x21    *x42    
            +coeff[ 18]    *x23            
            +coeff[ 19]    *x22    *x41*x51
            +coeff[ 20]            *x44*x51
            +coeff[ 21]            *x42*x51
            +coeff[ 22]    *x21        *x52
            +coeff[ 23]    *x24            
            +coeff[ 24]    *x23        *x51
            +coeff[ 25]    *x21    *x45    
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 26]            *x43    
            +coeff[ 27]            *x43*x51
            +coeff[ 28]    *x21    *x42*x51
            +coeff[ 29]            *x42*x52
            +coeff[ 30]    *x21    *x41*x52
            +coeff[ 31]    *x23    *x42    
            +coeff[ 32]    *x22    *x41*x52
            +coeff[ 33]    *x24        *x52
            +coeff[ 34]    *x21    *x43    
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 35]                *x53
            +coeff[ 36]    *x22    *x42*x51
            +coeff[ 37]    *x21        *x53
            +coeff[ 38]    *x23    *x41*x51
            +coeff[ 39]    *x24        *x51
            +coeff[ 40]    *x21    *x44*x51
            +coeff[ 41]    *x21    *x41*x53
            +coeff[ 42]    *x22    *x45    
            +coeff[ 43]    *x22        *x53
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 44]    *x22    *x43    
            +coeff[ 45]    *x22        *x52
            +coeff[ 46]    *x21    *x43*x51
            +coeff[ 47]    *x24    *x41    
            +coeff[ 48]    *x21    *x42*x52
            +coeff[ 49]    *x23    *x43    
            +coeff[ 50]    *x24    *x42    
            +coeff[ 51]    *x23    *x42*x51
            +coeff[ 52]    *x24    *x41*x51
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 53]    *x23    *x41*x52
            +coeff[ 54]    *x22    *x41*x53
            +coeff[ 55]    *x23    *x41    
            +coeff[ 56]            *x45    
            +coeff[ 57]            *x45*x51
            +coeff[ 58]    *x22    *x42*x52
            +coeff[ 59]            *x43*x53
            +coeff[ 60]    *x21    *x44*x52
            +coeff[ 61]    *x23        *x53
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 62]    *x23    *x42*x52
            +coeff[ 63]            *x45*x53
            +coeff[ 64]            *x44    
            +coeff[ 65]            *x43*x52
            +coeff[ 66]    *x22    *x44    
            +coeff[ 67]    *x22    *x43*x51
            +coeff[ 68]            *x44*x52
            +coeff[ 69]    *x21    *x43*x52
            +coeff[ 70]    *x21    *x42*x53
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 71]    *x21    *x43*x53
            +coeff[ 72]    *x23    *x44*x51
            +coeff[ 73]    *x21    *x44    
            +coeff[ 74]    *x23    *x44    
            +coeff[ 75]    *x21    *x45*x51
            +coeff[ 76]    *x24    *x43    
            +coeff[ 77]    *x23    *x43*x51
            +coeff[ 78]    *x23    *x45    
            +coeff[ 79]    *x22    *x43*x52
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 80]    *x24    *x44    
            +coeff[ 81]    *x22    *x45*x51
            +coeff[ 82]    *x24    *x41*x52
            +coeff[ 83]    *x23    *x41*x53
            +coeff[ 84]            *x42*x53
            ;
     
        return v_y_l5p77_fp           ;
    }
    float p_l5p77_fp           (float *x,int m){
        int ncoeff= 36;
        float avdat= -0.5212206E-02;
        float xmin[10]={
             0.00000E+00,-0.42303E-01, 0.00000E+00,-0.21263E-01,-0.46243E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.43132E-01, 0.00000E+00, 0.27251E-01, 0.49919E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 37]={
            -0.10062685E-02, 0.47768536E-02,-0.24497107E-01, 0.74334247E-02,
            -0.11544489E-01, 0.15763357E-01,-0.24463811E-02,-0.16722651E-01,
             0.25051001E-02, 0.18284488E-02, 0.11857779E-01, 0.29587157E-02,
            -0.27623686E-02, 0.11700285E-02, 0.23239984E-02, 0.58169833E-04,
            -0.11140255E-03, 0.42260094E-02, 0.29616689E-02,-0.10452799E-03,
            -0.25269401E-03,-0.16230006E-02, 0.13514395E-04,-0.10174111E-02,
            -0.20803536E-03,-0.11754299E-02,-0.25135899E-03,-0.22430689E-03,
             0.34520283E-03, 0.62453549E-03,-0.16864050E-02,-0.10595344E-02,
             0.70917340E-04, 0.20325230E-03,-0.50142757E-03, 0.17308090E-03,
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
        float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
        float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
    //  set up monomials   functions
        float x21 = x2;
        float x22 = x21*x2;
        float x23 = x22*x2;
        float x24 = x23*x2;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x44 = x43*x4;
        float x45 = x44*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
     
    //                 function
     
        float v_p_l5p77_fp           =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]            *x41    
            +coeff[  3]                *x51
            +coeff[  4]    *x22            
            +coeff[  5]    *x21    *x41    
            +coeff[  6]    *x21        *x51
            +coeff[  7]            *x41*x51
        ;
        v_p_l5p77_fp           =v_p_l5p77_fp           
            +coeff[  8]    *x23            
            +coeff[  9]                *x52
            +coeff[ 10]    *x22    *x41    
            +coeff[ 11]    *x21    *x42    
            +coeff[ 12]    *x22        *x51
            +coeff[ 13]            *x43    
            +coeff[ 14]            *x41*x52
            +coeff[ 15]    *x22    *x42    
            +coeff[ 16]            *x44    
        ;
        v_p_l5p77_fp           =v_p_l5p77_fp           
            +coeff[ 17]            *x42    
            +coeff[ 18]    *x24            
            +coeff[ 19]    *x21    *x41*x51
            +coeff[ 20]            *x42*x51
            +coeff[ 21]    *x21    *x43    
            +coeff[ 22]    *x22    *x41*x51
            +coeff[ 23]            *x41*x53
            +coeff[ 24]    *x22        *x53
            +coeff[ 25]    *x23    *x41    
        ;
        v_p_l5p77_fp           =v_p_l5p77_fp           
            +coeff[ 26]    *x23        *x51
            +coeff[ 27]                *x53
            +coeff[ 28]    *x21    *x42*x51
            +coeff[ 29]    *x24    *x41    
            +coeff[ 30]    *x23    *x42    
            +coeff[ 31]    *x22    *x42*x51
            +coeff[ 32]            *x45*x51
            +coeff[ 33]    *x22        *x52
            +coeff[ 34]            *x43*x51
        ;
        v_p_l5p77_fp           =v_p_l5p77_fp           
            +coeff[ 35]            *x42*x52
            ;
     
        return v_p_l5p77_fp           ;
    }
    float l_l5p77_fp           (float *x,int m){
        int ncoeff= 31;
        float avdat= -0.9646633E-02;
        float xmin[10]={
             0.00000E+00,-0.42303E-01, 0.00000E+00,-0.21263E-01,-0.46243E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.43132E-01, 0.00000E+00, 0.27251E-01, 0.49919E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 32]={
             0.71470574E-02,-0.24359792E+00, 0.70879129E-02,-0.32579765E-01,
            -0.17616162E-01,-0.60604174E-01,-0.25576519E-01,-0.14712851E-01,
             0.24701739E-01,-0.17694095E-01, 0.23387646E-01, 0.28066363E-03,
            -0.25372456E-02, 0.47503631E-02,-0.75330110E-02, 0.17586849E-02,
            -0.26919860E-02,-0.22582305E-02,-0.49126702E-02, 0.13446944E-02,
             0.43809912E-02, 0.51880009E-02,-0.77715347E-03, 0.35529703E-03,
             0.10507976E-02,-0.44094832E-03,-0.79061236E-03,-0.65777346E-03,
             0.49621559E-03,-0.65502850E-03,-0.20575426E-02,
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
        float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
        float x5 =1.+(x[  4]-xmax[  4])*scale[  4];
    //  set up monomials   functions
        float x21 = x2;
        float x22 = x21*x2;
        float x23 = x22*x2;
        float x24 = x23*x2;
        float x41 = x4;
        float x42 = x41*x4;
        float x43 = x42*x4;
        float x44 = x43*x4;
        float x51 = x5;
        float x52 = x51*x5;
        float x53 = x52*x5;
        float x54 = x53*x5;
     
    //                 function
     
        float v_l_l5p77_fp           =avdat
            +coeff[  0]                    
            +coeff[  1]    *x21            
            +coeff[  2]            *x41    
            +coeff[  3]                *x51
            +coeff[  4]    *x22            
            +coeff[  5]    *x21    *x41    
            +coeff[  6]                *x52
            +coeff[  7]    *x23            
        ;
        v_l_l5p77_fp           =v_l_l5p77_fp           
            +coeff[  8]    *x21    *x42    
            +coeff[  9]    *x22    *x41    
            +coeff[ 10]    *x23    *x41    
            +coeff[ 11]            *x44    
            +coeff[ 12]    *x21        *x51
            +coeff[ 13]                *x53
            +coeff[ 14]            *x42    
            +coeff[ 15]            *x41*x51
            +coeff[ 16]    *x24            
        ;
        v_l_l5p77_fp           =v_l_l5p77_fp           
            +coeff[ 17]    *x21        *x52
            +coeff[ 18]    *x21    *x43    
            +coeff[ 19]    *x21    *x41*x51
            +coeff[ 20]    *x22    *x42    
            +coeff[ 21]    *x24    *x41    
            +coeff[ 22]    *x24        *x51
            +coeff[ 23]            *x43*x52
            +coeff[ 24]            *x43    
            +coeff[ 25]    *x23        *x51
        ;
        v_l_l5p77_fp           =v_l_l5p77_fp           
            +coeff[ 26]    *x21    *x41*x52
            +coeff[ 27]            *x42*x52
            +coeff[ 28]    *x21        *x53
            +coeff[ 29]                *x54
            +coeff[ 30]    *x23    *x42    
            ;
     
        return v_l_l5p77_fp           ;
    }
}
