#include "Fwd_l5p77_484816R00.h"

namespace S484816R00
{
    float x_l5p77_fp           (float *x,int m){
        int ncoeff= 23;
        float avdat= -0.4057213E-03;
        float xmin[10]={
             0.00000E+00,-0.42213E-01, 0.00000E+00,-0.20331E-01,-0.46868E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.42434E-01, 0.00000E+00, 0.28594E-01, 0.49812E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 24]={
             0.27927356E-01, 0.64675194E+00,-0.51907498E-01, 0.42320278E-01,
             0.99067399E-02, 0.20589093E-01, 0.69121439E-02,-0.92437631E-02,
            -0.18613795E-01, 0.13248004E-02,-0.55274409E-02, 0.55988682E-02,
             0.10481548E-01, 0.49499404E-02,-0.12606413E-01, 0.22531396E-02,
            -0.29221540E-02, 0.37154970E-02, 0.43875864E-02, 0.10083147E-02,
             0.15206177E-02,-0.27517381E-02,-0.29351772E-02,
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
            +coeff[  7]            *x42    
        ;
        v_x_l5p77_fp           =v_x_l5p77_fp           
            +coeff[  8]    *x21    *x42    
            +coeff[  9]            *x43    
            +coeff[ 10]            *x41    
            +coeff[ 11]                *x53
            +coeff[ 12]    *x21    *x41*x51
            +coeff[ 13]    *x23            
            +coeff[ 14]    *x23    *x41    
            +coeff[ 15]            *x41*x51
            +coeff[ 16]    *x21        *x52
        ;
        v_x_l5p77_fp           =v_x_l5p77_fp           
            +coeff[ 17]    *x22    *x41    
            +coeff[ 18]    *x22    *x42    
            +coeff[ 19]    *x22        *x51
            +coeff[ 20]    *x23        *x51
            +coeff[ 21]    *x21    *x42*x51
            +coeff[ 22]    *x23    *x41*x51
            ;
     
        return v_x_l5p77_fp           ;
    }
    float t_l5p77_fp           (float *x,int m){
        int ncoeff= 46;
        float avdat= -0.9912001E-03;
        float xmin[10]={
             0.00000E+00,-0.42213E-01, 0.00000E+00,-0.20331E-01,-0.46868E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.42434E-01, 0.00000E+00, 0.28594E-01, 0.49812E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 47]={
             0.51136850E-02,-0.17186437E-01,-0.54432201E-03, 0.11193797E+00,
             0.53213462E-02,-0.10349713E-01, 0.95745834E-03,-0.15845710E-02,
            -0.77139970E-03,-0.16375001E-02, 0.15442029E-02,-0.69137476E-03,
             0.35441341E-03, 0.10791699E-02,-0.44518885E-04, 0.10309876E-02,
             0.15370976E-03, 0.21574319E-03,-0.66919718E-03,-0.16340498E-03,
            -0.48533361E-03,-0.68984012E-03,-0.20949714E-03,-0.42114366E-03,
            -0.27739193E-03, 0.38362516E-03,-0.42167338E-03, 0.29393076E-03,
            -0.13564230E-03,-0.91231683E-04,-0.20788222E-04,-0.18678306E-03,
             0.12640761E-03, 0.48325339E-03,-0.45366782E-04,-0.25569092E-03,
             0.70281225E-04, 0.36675998E-03, 0.13821252E-03, 0.32458329E-03,
             0.23734341E-03,-0.12609136E-03, 0.20932262E-03,-0.17536852E-03,
             0.21441134E-03, 0.16928454E-03,
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
            +coeff[  9]    *x21    *x42    
            +coeff[ 10]    *x21    *x41*x51
            +coeff[ 11]    *x21        *x52
            +coeff[ 12]            *x41*x51
            +coeff[ 13]                *x53
            +coeff[ 14]    *x22    *x41    
            +coeff[ 15]    *x22    *x42    
            +coeff[ 16]            *x43    
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 17]            *x42*x51
            +coeff[ 18]    *x23    *x41    
            +coeff[ 19]    *x23        *x51
            +coeff[ 20]    *x21    *x42*x51
            +coeff[ 21]    *x23    *x41*x51
            +coeff[ 22]    *x23        *x52
            +coeff[ 23]    *x21    *x43    
            +coeff[ 24]    *x22        *x52
            +coeff[ 25]    *x22    *x43    
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 26]    *x22    *x41*x52
            +coeff[ 27]    *x21    *x42*x52
            +coeff[ 28]            *x42*x53
            +coeff[ 29]    *x23            
            +coeff[ 30]    *x22    *x41*x51
            +coeff[ 31]            *x42*x52
            +coeff[ 32]            *x41*x53
            +coeff[ 33]    *x23        *x53
            +coeff[ 34]    *x21    *x41*x52
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 35]    *x21        *x53
            +coeff[ 36]    *x23    *x42    
            +coeff[ 37]    *x22    *x42*x51
            +coeff[ 38]    *x21    *x43*x51
            +coeff[ 39]    *x23    *x43    
            +coeff[ 40]    *x22    *x42*x52
            +coeff[ 41]    *x22    *x41*x53
            +coeff[ 42]    *x21    *x42*x53
            +coeff[ 43]    *x22        *x51
        ;
        v_t_l5p77_fp           =v_t_l5p77_fp           
            +coeff[ 44]    *x22        *x53
            +coeff[ 45]    *x22    *x43*x51
            ;
     
        return v_t_l5p77_fp           ;
    }
    float y_l5p77_fp           (float *x,int m){
        int ncoeff= 84;
        float avdat= -0.5412241E-02;
        float xmin[10]={
             0.00000E+00,-0.42213E-01, 0.00000E+00,-0.20331E-01,-0.46868E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.42434E-01, 0.00000E+00, 0.28594E-01, 0.49812E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 85]={
            -0.16194639E-02,-0.46509203E-01, 0.28970672E-02, 0.85074054E-02,
             0.39926860E-02, 0.69209388E-02,-0.11781494E-01,-0.37042622E-02,
            -0.27663426E-02, 0.46644160E-02,-0.51575922E-02, 0.53774025E-02,
             0.28088323E-02,-0.34969666E-02,-0.35629747E-02, 0.11319256E-03,
             0.16109240E-02,-0.83432742E-03, 0.16139270E-02, 0.13520198E-02,
            -0.31770347E-02, 0.33817854E-04, 0.49833205E-04,-0.14165747E-02,
             0.14134701E-02,-0.99670934E-03,-0.70667936E-03,-0.54037466E-03,
            -0.50980103E-03, 0.91024925E-03, 0.48612483E-03,-0.77000511E-03,
             0.26016985E-02,-0.14426673E-02,-0.71289937E-03, 0.58864523E-03,
             0.30863748E-02,-0.15571704E-02,-0.32013498E-03, 0.14771263E-02,
            -0.97798475E-03,-0.58833865E-03,-0.13825801E-02, 0.81966928E-03,
            -0.40615257E-03,-0.32971299E-03, 0.21102573E-02, 0.10493076E-02,
            -0.98815060E-03,-0.64805849E-03,-0.95724652E-03, 0.52473496E-03,
             0.88926079E-03,-0.58698701E-03, 0.15946242E-03, 0.13998246E-02,
            -0.31664287E-03, 0.90902223E-03, 0.67739049E-03,-0.63015491E-03,
            -0.28968728E-03, 0.13929838E-02,-0.12727812E-04,-0.99349309E-04,
            -0.12006356E-03,-0.69493492E-03, 0.49658887E-04, 0.33153282E-03,
             0.47796682E-03, 0.37374029E-04, 0.11979853E-02,-0.42491218E-04,
            -0.29967868E-03,-0.33044678E-03,-0.53633016E-03, 0.78529376E-03,
             0.15273596E-03, 0.96785807E-05, 0.36861323E-03,-0.95608753E-04,
             0.56879758E-03, 0.19971219E-03, 0.10350346E-03,-0.33051541E-03,
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
            +coeff[  9]    *x22    *x41    
            +coeff[ 10]    *x21    *x41*x51
            +coeff[ 11]            *x41*x52
            +coeff[ 12]    *x22        *x51
            +coeff[ 13]    *x22    *x42    
            +coeff[ 14]            *x41*x53
            +coeff[ 15]    *x23        *x52
            +coeff[ 16]    *x21    *x42    
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 17]                *x52
            +coeff[ 18]    *x23            
            +coeff[ 19]    *x24            
            +coeff[ 20]    *x22    *x41*x51
            +coeff[ 21]            *x44*x51
            +coeff[ 22]            *x43    
            +coeff[ 23]            *x42*x51
            +coeff[ 24]    *x21        *x52
            +coeff[ 25]    *x23        *x51
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 26]    *x21    *x43    
            +coeff[ 27]            *x43*x51
            +coeff[ 28]    *x21    *x42*x51
            +coeff[ 29]            *x42*x52
            +coeff[ 30]    *x21    *x41*x52
            +coeff[ 31]    *x23    *x42    
            +coeff[ 32]    *x22    *x41*x52
            +coeff[ 33]    *x22    *x43*x51
            +coeff[ 34]    *x22        *x53
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 35]    *x24        *x52
            +coeff[ 36]    *x22    *x43    
            +coeff[ 37]    *x22    *x42*x51
            +coeff[ 38]    *x21        *x53
            +coeff[ 39]    *x23    *x41*x51
            +coeff[ 40]    *x24        *x51
            +coeff[ 41]    *x21    *x41*x53
            +coeff[ 42]    *x24    *x41*x51
            +coeff[ 43]    *x21    *x43*x51
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 44]    *x21    *x42*x52
            +coeff[ 45]            *x45*x51
            +coeff[ 46]    *x24    *x42    
            +coeff[ 47]    *x23    *x42*x51
            +coeff[ 48]    *x22    *x45    
            +coeff[ 49]    *x24    *x43    
            +coeff[ 50]            *x43*x53
            +coeff[ 51]            *x45*x52
            +coeff[ 52]    *x22    *x41*x53
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 53]    *x24    *x44    
            +coeff[ 54]    *x24    *x45    
            +coeff[ 55]    *x24    *x41    
            +coeff[ 56]            *x43*x52
            +coeff[ 57]    *x22    *x44    
            +coeff[ 58]    *x22    *x42*x52
            +coeff[ 59]    *x24    *x42*x51
            +coeff[ 60]    *x23        *x53
            +coeff[ 61]    *x22    *x45*x51
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 62]            *x44    
            +coeff[ 63]                *x53
            +coeff[ 64]    *x22        *x52
            +coeff[ 65]    *x23    *x43    
            +coeff[ 66]    *x21    *x44*x51
            +coeff[ 67]            *x42*x53
            +coeff[ 68]    *x21    *x43*x52
            +coeff[ 69]    *x21    *x44*x52
            +coeff[ 70]    *x23    *x45    
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 71]            *x44*x53
            +coeff[ 72]    *x24    *x41*x52
            +coeff[ 73]    *x22    *x42*x53
            +coeff[ 74]    *x21    *x45*x52
            +coeff[ 75]            *x45*x53
            +coeff[ 76]            *x45    
            +coeff[ 77]    *x21    *x44    
            +coeff[ 78]    *x23    *x44    
            +coeff[ 79]    *x21    *x45*x51
        ;
        v_y_l5p77_fp           =v_y_l5p77_fp           
            +coeff[ 80]    *x22    *x44*x51
            +coeff[ 81]    *x24        *x53
            +coeff[ 82]    *x23    *x41    
            +coeff[ 83]    *x21    *x45    
            ;
     
        return v_y_l5p77_fp           ;
    }
    float p_l5p77_fp           (float *x,int m){
        int ncoeff= 39;
        float avdat= -0.1528836E-02;
        float xmin[10]={
             0.00000E+00,-0.42213E-01, 0.00000E+00,-0.20331E-01,-0.46868E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.42434E-01, 0.00000E+00, 0.28594E-01, 0.49812E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 40]={
            -0.25555800E-03, 0.43955781E-02,-0.24177471E-01, 0.68141883E-02,
            -0.10936971E-01, 0.15582929E-01,-0.23055440E-02,-0.16874116E-01,
             0.19291786E-02, 0.11110643E-01,-0.25677199E-02, 0.11473084E-02,
             0.44657223E-03,-0.14726779E-02,-0.26563031E-04, 0.46797340E-04,
             0.39654179E-02, 0.22484625E-02, 0.28143423E-02, 0.28600094E-02,
             0.22708126E-02,-0.21281917E-03,-0.11388632E-02,-0.77964884E-04,
            -0.15228873E-02,-0.99388976E-03,-0.22724760E-03,-0.17626054E-03,
             0.28760952E-03, 0.80250227E-03,-0.12316373E-02, 0.19664598E-03,
             0.85690896E-04,-0.18514864E-03, 0.30937817E-03,-0.66065014E-03,
             0.24509092E-03, 0.53002150E-03, 0.35733075E-03,
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
            +coeff[  8]                *x52
            +coeff[  9]    *x22    *x41    
            +coeff[ 10]    *x22        *x51
            +coeff[ 11]            *x43    
            +coeff[ 12]    *x22    *x42    
            +coeff[ 13]    *x21    *x43    
            +coeff[ 14]            *x44    
            +coeff[ 15]            *x43*x52
            +coeff[ 16]            *x42    
        ;
        v_p_l5p77_fp           =v_p_l5p77_fp           
            +coeff[ 17]    *x23            
            +coeff[ 18]    *x21    *x42    
            +coeff[ 19]    *x24            
            +coeff[ 20]            *x41*x52
            +coeff[ 21]            *x42*x51
            +coeff[ 22]    *x23    *x41    
            +coeff[ 23]    *x22    *x41*x51
            +coeff[ 24]    *x23    *x42    
            +coeff[ 25]            *x41*x53
        ;
        v_p_l5p77_fp           =v_p_l5p77_fp           
            +coeff[ 26]    *x22        *x53
            +coeff[ 27]    *x23        *x51
            +coeff[ 28]    *x21    *x42*x51
            +coeff[ 29]    *x24    *x41    
            +coeff[ 30]    *x22    *x42*x51
            +coeff[ 31]            *x45*x51
            +coeff[ 32]    *x21        *x52
            +coeff[ 33]                *x53
            +coeff[ 34]    *x22        *x52
        ;
        v_p_l5p77_fp           =v_p_l5p77_fp           
            +coeff[ 35]            *x43*x51
            +coeff[ 36]            *x42*x52
            +coeff[ 37]    *x22    *x43    
            +coeff[ 38]    *x22    *x41*x52
            ;
     
        return v_p_l5p77_fp           ;
    }
    float l_l5p77_fp           (float *x,int m){
        int ncoeff= 31;
        float avdat= -0.1067273E-01;
        float xmin[10]={
             0.00000E+00,-0.42213E-01, 0.00000E+00,-0.20331E-01,-0.46868E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float xmax[10]={
             0.00000E+00, 0.42434E-01, 0.00000E+00, 0.28594E-01, 0.49812E-01,
             0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
        float scale[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        float coeff[ 32]={
             0.10541497E-01,-0.24550283E+00, 0.75290087E-02,-0.32619953E-01,
            -0.17541528E-01,-0.57822336E-01,-0.25877362E-01,-0.13452678E-01,
             0.22531474E-01,-0.17184380E-01, 0.22171482E-01, 0.10324860E-03,
            -0.22631194E-02, 0.51613008E-02,-0.75701410E-02, 0.21147006E-02,
            -0.24106612E-02,-0.25993087E-02,-0.40057334E-02, 0.13928221E-02,
             0.12321208E-02, 0.10608631E-02, 0.46387408E-02,-0.11861685E-02,
             0.35123087E-02,-0.98140747E-03,-0.71536703E-03,-0.91541611E-03,
            -0.12129236E-02, 0.13680777E-02,-0.89149358E-03,
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
            +coeff[ 16]    *x21        *x52
        ;
        v_l_l5p77_fp           =v_l_l5p77_fp           
            +coeff[ 17]    *x24            
            +coeff[ 18]    *x21    *x43    
            +coeff[ 19]    *x22    *x44    
            +coeff[ 20]            *x43    
            +coeff[ 21]    *x21    *x41*x51
            +coeff[ 22]    *x24    *x41    
            +coeff[ 23]    *x22        *x53
            +coeff[ 24]    *x22    *x42    
            +coeff[ 25]    *x21    *x41*x52
        ;
        v_l_l5p77_fp           =v_l_l5p77_fp           
            +coeff[ 26]            *x42*x52
            +coeff[ 27]                *x54
            +coeff[ 28]    *x23    *x42    
            +coeff[ 29]    *x21    *x44    
            +coeff[ 30]    *x22    *x43*x51
            ;
     
        return v_l_l5p77_fp           ;
    }
}
