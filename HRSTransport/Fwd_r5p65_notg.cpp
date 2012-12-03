float x_r5p65_ep11_q1en                       (float *x,int m){
    int ncoeff= 32;
    float avdat= -0.4991575E-03;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 33]={
        -0.27636946E-04, 0.89440299E-02, 0.10536432E+00,-0.65390519E-02,
        -0.10723757E-02,-0.18883918E-02, 0.24730754E-02,-0.23124486E-02,
        -0.53920265E-03, 0.11172722E-02,-0.12308083E-02, 0.23885707E-02,
        -0.28216289E-03, 0.25426058E-03, 0.57361816E-03, 0.49491430E-03,
         0.10191972E-02,-0.14418154E-03, 0.16407679E-03,-0.14925346E-03,
        -0.30582375E-03, 0.22184663E-03, 0.86174626E-03,-0.77350478E-03,
        -0.89629943E-03,-0.69729860E-04, 0.15737257E-03, 0.12735001E-03,
        -0.14342443E-03, 0.11200690E-03,-0.15110751E-03, 0.24724077E-03,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;

//                 function

    float v_x_r5p65_ep11_q1en                       =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]    *x21    *x41    
        +coeff[  4]*x11        *x41    
        +coeff[  5]    *x21*x31        
        +coeff[  6]    *x23            
        +coeff[  7]    *x21    *x42    
    ;
    v_x_r5p65_ep11_q1en                       =v_x_r5p65_ep11_q1en                       
        +coeff[  8]            *x41    
        +coeff[  9]*x11*x22            
        +coeff[ 10]    *x21*x31*x41    
        +coeff[ 11]    *x23    *x41    
        +coeff[ 12]*x11    *x31        
        +coeff[ 13]    *x21        *x51
        +coeff[ 14]    *x22            
        +coeff[ 15]    *x22    *x41    
        +coeff[ 16]*x11*x22    *x41    
    ;
    v_x_r5p65_ep11_q1en                       =v_x_r5p65_ep11_q1en                       
        +coeff[ 17]        *x31        
        +coeff[ 18]*x11*x21            
        +coeff[ 19]            *x42    
        +coeff[ 20]*x11        *x42    
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]    *x23*x31        
        +coeff[ 23]    *x21    *x43    
        +coeff[ 24]    *x23*x31*x42    
        +coeff[ 25]        *x31*x41    
    ;
    v_x_r5p65_ep11_q1en                       =v_x_r5p65_ep11_q1en                       
        +coeff[ 26]*x12*x21            
        +coeff[ 27]*x11*x21    *x41    
        +coeff[ 28]*x11    *x31*x41    
        +coeff[ 29]    *x22*x31        
        +coeff[ 30]    *x21*x32        
        +coeff[ 31]*x11*x22*x31        
        ;

    return v_x_r5p65_ep11_q1en                       ;
}
float t_r5p65_ep11_q1en                       (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.4347269E-03;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 51]={
        -0.17659836E-03,-0.12205893E-02, 0.40857319E-01,-0.86565930E-02,
         0.37774324E-04,-0.61795500E-03,-0.21885920E-02,-0.11939370E-02,
         0.25676957E-02,-0.26388890E-02, 0.37009914E-02, 0.57507568E-03,
         0.12088970E-02,-0.13248386E-02,-0.16260751E-03, 0.16983130E-03,
        -0.33230861E-03, 0.31740047E-03, 0.72471512E-03,-0.32838309E-03,
         0.14080643E-02,-0.11277000E-02,-0.19689686E-03,-0.18054199E-03,
         0.18789855E-03, 0.33495590E-03, 0.86491497E-03,-0.79870847E-03,
        -0.10315656E-03, 0.18205191E-03, 0.17957864E-03,-0.17292256E-03,
         0.80427752E-04, 0.40193781E-03,-0.18523735E-03, 0.20617912E-04,
         0.39536579E-04, 0.47978414E-04,-0.61240018E-04, 0.19796740E-03,
         0.13446820E-03,-0.14741081E-03, 0.39113811E-03, 0.58460686E-04,
        -0.20218776E-03, 0.13188583E-04, 0.17810591E-04,-0.24696246E-04,
         0.69783608E-04, 0.41593765E-04,
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

//                 function

    float v_t_r5p65_ep11_q1en                       =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]    *x21    *x41    
        +coeff[  4]*x12*x21*x31        
        +coeff[  5]            *x41    
        +coeff[  6]    *x21*x31        
        +coeff[  7]*x11        *x41    
    ;
    v_t_r5p65_ep11_q1en                       =v_t_r5p65_ep11_q1en                       
        +coeff[  8]    *x23            
        +coeff[  9]    *x21    *x42    
        +coeff[ 10]    *x23    *x41    
        +coeff[ 11]    *x22            
        +coeff[ 12]*x11*x22            
        +coeff[ 13]    *x21*x31*x41    
        +coeff[ 14]        *x31        
        +coeff[ 15]*x11*x21            
        +coeff[ 16]*x11    *x31        
    ;
    v_t_r5p65_ep11_q1en                       =v_t_r5p65_ep11_q1en                       
        +coeff[ 17]    *x21        *x51
        +coeff[ 18]    *x22    *x41    
        +coeff[ 19]*x11        *x42    
        +coeff[ 20]*x11*x22    *x41    
        +coeff[ 21]    *x21    *x43    
        +coeff[ 22]            *x42    
        +coeff[ 23]    *x21*x32        
        +coeff[ 24]*x11*x21    *x41    
        +coeff[ 25]    *x21    *x41*x51
    ;
    v_t_r5p65_ep11_q1en                       =v_t_r5p65_ep11_q1en                       
        +coeff[ 26]    *x23*x31        
        +coeff[ 27]    *x21*x31*x42    
        +coeff[ 28]        *x31*x41    
        +coeff[ 29]*x12*x21            
        +coeff[ 30]    *x22*x31        
        +coeff[ 31]*x11    *x31*x41    
        +coeff[ 32]    *x21*x31    *x51
        +coeff[ 33]*x11*x22*x31        
        +coeff[ 34]    *x21*x32*x41    
    ;
    v_t_r5p65_ep11_q1en                       =v_t_r5p65_ep11_q1en                       
        +coeff[ 35]                *x51
        +coeff[ 36]*x11            *x51
        +coeff[ 37]*x11*x21*x31        
        +coeff[ 38]            *x43    
        +coeff[ 39]*x12*x21    *x41    
        +coeff[ 40]    *x22    *x42    
        +coeff[ 41]*x11        *x43    
        +coeff[ 42]    *x23    *x42    
        +coeff[ 43]*x13        *x41*x51
    ;
    v_t_r5p65_ep11_q1en                       =v_t_r5p65_ep11_q1en                       
        +coeff[ 44]*x11*x22*x31*x42    
        +coeff[ 45]*x12                
        +coeff[ 46]            *x41*x51
        +coeff[ 47]        *x31*x42    
        +coeff[ 48]    *x22*x31*x41    
        +coeff[ 49]*x11*x21    *x42    
        ;

    return v_t_r5p65_ep11_q1en                       ;
}
float y_r5p65_ep11_q1en                       (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.1470180E-01;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
         0.69188224E-02, 0.65996148E-01,-0.10138815E-02, 0.10943611E-01,
         0.41649677E-02,-0.68541830E-02, 0.22044503E-02,-0.19680823E-02,
        -0.42993492E-02, 0.11774013E-02,-0.13082414E-02,-0.15473836E-03,
        -0.67890243E-03,-0.11736796E-02,-0.27979908E-02,-0.15022063E-02,
        -0.18772703E-03,-0.15806234E-03,-0.21115389E-03,-0.37005809E-03,
        -0.82857144E-03, 0.10087470E-02, 0.57253090E-03, 0.14563602E-03,
        -0.10894102E-03,-0.13669262E-03, 0.45904523E-03,-0.35067531E-03,
         0.31206224E-03,-0.17278288E-03, 0.30194485E-03, 0.13147609E-03,
         0.27817284E-03,-0.44692747E-03,-0.31648528E-04, 0.42035677E-04,
         0.73744028E-04,-0.10402793E-03, 0.77346871E-04,-0.16069881E-03,
        -0.82957488E-03, 0.14309301E-03, 0.39228701E-03,-0.43567525E-04,
        -0.52516607E-04,-0.88269553E-04,-0.53402877E-04, 0.17673918E-03,
         0.85216401E-04, 0.14422236E-03,-0.66333957E-03,-0.28403589E-03,
         0.22276714E-03,-0.18763538E-04,-0.27613551E-04, 0.21963919E-04,
        -0.17818364E-04,-0.33992124E-04,-0.46261077E-04, 0.57878049E-04,
         0.48621150E-04, 0.11438600E-03, 0.31976333E-05,-0.84411644E-04,
        -0.18026687E-03, 0.10674591E-03, 0.73170282E-04, 0.75033800E-04,
         0.71827241E-03, 0.40858219E-03, 0.40971473E-03,-0.42458778E-04,
         0.25439833E-03, 0.16984213E-04, 0.65705890E-05,-0.41578573E-04,
        -0.15222395E-03,-0.29257741E-04, 0.46059944E-04,-0.38481878E-04,
         0.26138312E-04, 0.15263093E-04, 0.43851003E-04, 0.60766590E-04,
         0.97158692E-04, 0.10792167E-03, 0.32522436E-04, 0.42955951E-04,
        -0.98573430E-04, 0.35533867E-04, 0.65447621E-05, 0.18484816E-04,
         0.52463450E-03, 0.42995511E-03, 0.11574897E-04, 0.38724396E-03,
         0.22784510E-03,
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
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x45 = x44*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_r5p65_ep11_q1en                       =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]            *x42    
        +coeff[  7]*x11*x21            
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[  8]    *x22    *x41    
        +coeff[  9]        *x31*x41    
        +coeff[ 10]*x11*x21    *x41    
        +coeff[ 11]*x11                
        +coeff[ 12]    *x21    *x41    
        +coeff[ 13]    *x22*x31        
        +coeff[ 14]    *x22    *x42    
        +coeff[ 15]    *x22*x31*x41    
        +coeff[ 16]    *x21*x31        
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 17]*x12                
        +coeff[ 18]                *x52
        +coeff[ 19]*x11*x21*x31        
        +coeff[ 20]*x11*x21    *x42    
        +coeff[ 21]    *x24            
        +coeff[ 22]*x11*x23            
        +coeff[ 23]        *x32        
        +coeff[ 24]*x11        *x41    
        +coeff[ 25]            *x41*x51
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 26]            *x43    
        +coeff[ 27]    *x21    *x42    
        +coeff[ 28]        *x31*x42    
        +coeff[ 29]    *x21*x31*x41    
        +coeff[ 30]    *x23            
        +coeff[ 31]*x11*x22            
        +coeff[ 32]    *x22        *x51
        +coeff[ 33]*x11*x21*x31*x41    
        +coeff[ 34]*x11    *x31        
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 35]    *x21        *x51
        +coeff[ 36]        *x32*x41    
        +coeff[ 37]*x12        *x41    
        +coeff[ 38]*x11*x21        *x51
        +coeff[ 39]    *x22*x32        
        +coeff[ 40]    *x22    *x43    
        +coeff[ 41]*x12*x22            
        +coeff[ 42]    *x24    *x41    
        +coeff[ 43]        *x31    *x51
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 44]*x11        *x42    
        +coeff[ 45]            *x42*x51
        +coeff[ 46]        *x31*x41*x51
        +coeff[ 47]    *x23    *x41    
        +coeff[ 48]*x11*x22    *x41    
        +coeff[ 49]    *x22    *x41*x51
        +coeff[ 50]    *x22*x31*x42    
        +coeff[ 51]*x11*x21    *x43    
        +coeff[ 52]*x11*x23    *x41    
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 53]    *x21*x32        
        +coeff[ 54]*x11    *x31*x41    
        +coeff[ 55]*x12*x21            
        +coeff[ 56]*x12    *x31        
        +coeff[ 57]*x12        *x42    
        +coeff[ 58]*x11*x21*x32        
        +coeff[ 59]    *x22*x31    *x51
        +coeff[ 60]*x11*x21    *x41*x51
        +coeff[ 61]    *x23    *x42    
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 62]*x11    *x32    *x51
        +coeff[ 63]    *x22*x32*x41    
        +coeff[ 64]*x11*x21*x31*x42    
        +coeff[ 65]    *x24*x31        
        +coeff[ 66]*x11*x23*x31        
        +coeff[ 67]*x12*x22    *x41    
        +coeff[ 68]    *x24    *x42    
        +coeff[ 69]    *x24*x31*x41    
        +coeff[ 70]*x11*x23    *x42    
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 71]    *x22*x31*x42*x51
        +coeff[ 72]*x11*x23*x31*x41    
        +coeff[ 73]    *x21    *x41*x51
        +coeff[ 74]    *x21*x31    *x51
        +coeff[ 75]    *x21    *x43    
        +coeff[ 76]        *x31*x43    
        +coeff[ 77]    *x21*x31*x42    
        +coeff[ 78]    *x23*x31        
        +coeff[ 79]            *x45    
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 80]*x11*x22*x31        
        +coeff[ 81]*x11*x21*x31    *x51
        +coeff[ 82]    *x23*x31*x41    
        +coeff[ 83]*x11*x22    *x42    
        +coeff[ 84]    *x22    *x42*x51
        +coeff[ 85]        *x31*x45    
        +coeff[ 86]*x11*x22*x31*x41    
        +coeff[ 87]    *x22*x31*x41*x51
        +coeff[ 88]    *x22    *x44    
    ;
    v_y_r5p65_ep11_q1en                       =v_y_r5p65_ep11_q1en                       
        +coeff[ 89]*x11*x21    *x42*x51
        +coeff[ 90]*x13            *x51
        +coeff[ 91]    *x22    *x41*x52
        +coeff[ 92]    *x24    *x43    
        +coeff[ 93]    *x24*x31*x42    
        +coeff[ 94]*x12            *x53
        +coeff[ 95]*x11*x23    *x43    
        +coeff[ 96]*x11*x23*x31*x42    
        ;

    return v_y_r5p65_ep11_q1en                       ;
}
float p_r5p65_ep11_q1en                       (float *x,int m){
    int ncoeff= 44;
    float avdat= -0.7248341E-02;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 45]={
         0.42968234E-02,-0.12764777E-02, 0.11934669E-02, 0.32134432E-01,
         0.60629277E-02,-0.93016047E-02, 0.26913173E-02,-0.25046705E-02,
        -0.43729208E-02,-0.66004077E-03, 0.13416461E-02,-0.13794425E-02,
        -0.32694440E-02,-0.19769509E-03,-0.30602998E-03,-0.12477614E-02,
         0.16859433E-02,-0.15011561E-02,-0.17795716E-03,-0.14428659E-03,
         0.44533677E-03,-0.38321837E-03, 0.39732215E-03,-0.19374111E-03,
        -0.37533755E-03, 0.87959610E-03,-0.86209265E-03, 0.17250342E-03,
        -0.90952839E-04,-0.17084784E-03, 0.41418956E-03, 0.18574885E-03,
         0.10278732E-03,-0.41291671E-03, 0.38043712E-04, 0.45983979E-04,
        -0.40851479E-04, 0.66850567E-04, 0.22101261E-03,-0.85336913E-04,
        -0.17506888E-03,-0.10267644E-03,-0.33053794E-03, 0.17466230E-03,
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

    float v_p_r5p65_ep11_q1en                       =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]            *x42    
        +coeff[  7]*x11*x21            
    ;
    v_p_r5p65_ep11_q1en                       =v_p_r5p65_ep11_q1en                       
        +coeff[  8]    *x22    *x41    
        +coeff[  9]    *x21    *x41    
        +coeff[ 10]        *x31*x41    
        +coeff[ 11]*x11*x21    *x41    
        +coeff[ 12]    *x22    *x42    
        +coeff[ 13]*x11                
        +coeff[ 14]                *x52
        +coeff[ 15]    *x22*x31        
        +coeff[ 16]    *x24            
    ;
    v_p_r5p65_ep11_q1en                       =v_p_r5p65_ep11_q1en                       
        +coeff[ 17]    *x22*x31*x41    
        +coeff[ 18]    *x21*x31        
        +coeff[ 19]            *x41*x51
        +coeff[ 20]    *x23            
        +coeff[ 21]    *x21    *x42    
        +coeff[ 22]    *x22        *x51
        +coeff[ 23]*x12                
        +coeff[ 24]*x11*x21*x31        
        +coeff[ 25]*x11*x23            
    ;
    v_p_r5p65_ep11_q1en                       =v_p_r5p65_ep11_q1en                       
        +coeff[ 26]*x11*x21    *x42    
        +coeff[ 27]        *x32        
        +coeff[ 28]*x11        *x41    
        +coeff[ 29]    *x21*x31*x41    
        +coeff[ 30]            *x43    
        +coeff[ 31]*x11*x22            
        +coeff[ 32]*x11*x21        *x51
        +coeff[ 33]*x11*x21*x31*x41    
        +coeff[ 34]        *x33*x42    
    ;
    v_p_r5p65_ep11_q1en                       =v_p_r5p65_ep11_q1en                       
        +coeff[ 35]    *x21        *x51
        +coeff[ 36]        *x31    *x51
        +coeff[ 37]        *x32*x41    
        +coeff[ 38]        *x31*x42    
        +coeff[ 39]            *x42*x51
        +coeff[ 40]    *x22*x32        
        +coeff[ 41]*x12        *x41    
        +coeff[ 42]    *x22    *x43    
        +coeff[ 43]*x12*x22            
        ;

    return v_p_r5p65_ep11_q1en                       ;
}
float l_r6_h90s_ep11_q1en                     (float *x,int m){
    int ncoeff= 43;
    float avdat= -0.1338794E-04;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 44]={
         0.49421360E-03, 0.52247513E-04,-0.12261404E-02,-0.51474012E-02,
        -0.22790423E-02,-0.64769425E-04,-0.10309139E-02, 0.64058736E-03,
        -0.10019424E-03,-0.12523020E-03, 0.10314114E-03, 0.75925098E-04,
         0.15452586E-03, 0.27322609E-03,-0.49298909E-04,-0.15666243E-03,
         0.11403237E-03,-0.38882015E-04, 0.13522753E-03, 0.26216745E-04,
        -0.66594504E-04,-0.76010656E-04, 0.57797446E-04,-0.16647365E-03,
        -0.77794175E-04, 0.96141812E-05,-0.55204764E-05, 0.79315114E-05,
         0.97955162E-05,-0.27813448E-04, 0.13426770E-04, 0.19253968E-04,
        -0.37566610E-04, 0.19693243E-04, 0.10010575E-03, 0.69792557E-04,
         0.52160976E-05, 0.24722356E-05, 0.61071059E-05,-0.12907622E-04,
         0.47993462E-05, 0.55421224E-05,-0.14235044E-04,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_l_r6_h90s_ep11_q1en                     =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]    *x22            
        +coeff[  5]        *x31*x41    
        +coeff[  6]            *x42    
        +coeff[  7]    *x22    *x41    
    ;
    v_l_r6_h90s_ep11_q1en                     =v_l_r6_h90s_ep11_q1en                     
        +coeff[  8]                *x51
        +coeff[  9]            *x41*x51
        +coeff[ 10]*x11*x21            
        +coeff[ 11]    *x21    *x41    
        +coeff[ 12]*x11*x21    *x41    
        +coeff[ 13]    *x22    *x42    
        +coeff[ 14]    *x24*x31        
        +coeff[ 15]    *x24            
        +coeff[ 16]    *x22*x31*x41    
    ;
    v_l_r6_h90s_ep11_q1en                     =v_l_r6_h90s_ep11_q1en                     
        +coeff[ 17]    *x23            
        +coeff[ 18]    *x22*x31        
        +coeff[ 19]    *x21    *x42    
        +coeff[ 20]            *x43    
        +coeff[ 21]*x11*x23            
        +coeff[ 22]*x11*x21    *x42    
        +coeff[ 23]    *x24    *x41    
        +coeff[ 24]*x11*x23    *x41    
        +coeff[ 25]    *x21*x31        
    ;
    v_l_r6_h90s_ep11_q1en                     =v_l_r6_h90s_ep11_q1en                     
        +coeff[ 26]                *x52
        +coeff[ 27]*x11        *x41    
        +coeff[ 28]    *x21*x31*x41    
        +coeff[ 29]        *x31*x42    
        +coeff[ 30]    *x22        *x51
        +coeff[ 31]*x11*x21*x31        
        +coeff[ 32]    *x23    *x41    
        +coeff[ 33]*x11*x21*x31*x41    
        +coeff[ 34]    *x22    *x43    
    ;
    v_l_r6_h90s_ep11_q1en                     =v_l_r6_h90s_ep11_q1en                     
        +coeff[ 35]    *x24*x31*x42    
        +coeff[ 36]*x11                
        +coeff[ 37]    *x21        *x51
        +coeff[ 38]            *x41*x52
        +coeff[ 39]*x11*x22            
        +coeff[ 40]*x11*x21        *x51
        +coeff[ 41]*x12        *x41    
        +coeff[ 42]*x11*x22    *x41    
        ;

    return v_l_r6_h90s_ep11_q1en                     ;
}
float x_r5p65_ep13_q1                         (float *x,int m){
    int ncoeff= 42;
    float avdat= -0.1719703E-02;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.29972E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 43]={
         0.63523307E-03, 0.13040115E+00,-0.14533915E-01, 0.16407259E-04,
         0.61329091E-02,-0.37706804E-02,-0.11302715E-02,-0.22382056E-02,
         0.15949087E-02, 0.45506353E-02,-0.49657580E-02, 0.17455835E-04,
        -0.58211875E-03, 0.10968018E-02,-0.24782077E-02, 0.60271476E-02,
        -0.29200318E-03, 0.31937141E-03, 0.21561421E-02, 0.12337228E-02,
         0.25043874E-02, 0.15113550E-02,-0.49106049E-04, 0.14550782E-03,
        -0.30364541E-03, 0.32839156E-03, 0.29727098E-03,-0.30841891E-03,
         0.50041347E-03,-0.33727131E-03, 0.61737938E-03,-0.18733672E-02,
        -0.10682935E-03, 0.30540352E-04,-0.15062316E-03,-0.59256941E-03,
         0.29194346E-03, 0.31130720E-03,-0.31768979E-03,-0.12400271E-02,
         0.19364129E-03, 0.86911808E-03,
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

    float v_x_r5p65_ep13_q1                         =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]    *x21    *x41    
        +coeff[  3]*x13                
        +coeff[  4]*x11                
        +coeff[  5]    *x21*x31        
        +coeff[  6]            *x41    
        +coeff[  7]*x11        *x41    
    ;
    v_x_r5p65_ep13_q1                         =v_x_r5p65_ep13_q1                         
        +coeff[  8]    *x21        *x51
        +coeff[  9]    *x23            
        +coeff[ 10]    *x21    *x42    
        +coeff[ 11]*x13*x22            
        +coeff[ 12]*x11    *x31        
        +coeff[ 13]    *x22            
        +coeff[ 14]    *x21*x31*x41    
        +coeff[ 15]    *x23    *x41    
        +coeff[ 16]        *x31        
    ;
    v_x_r5p65_ep13_q1                         =v_x_r5p65_ep13_q1                         
        +coeff[ 17]*x11*x21            
        +coeff[ 18]*x11*x22            
        +coeff[ 19]    *x22    *x41    
        +coeff[ 20]*x11*x22    *x41    
        +coeff[ 21]    *x23*x31        
        +coeff[ 22]*x13        *x42    
        +coeff[ 23]*x11            *x51
        +coeff[ 24]            *x42    
        +coeff[ 25]*x12*x21            
    ;
    v_x_r5p65_ep13_q1                         =v_x_r5p65_ep13_q1                         
        +coeff[ 26]*x11*x21    *x41    
        +coeff[ 27]*x11    *x31*x41    
        +coeff[ 28]    *x21    *x41*x51
        +coeff[ 29]    *x21*x32        
        +coeff[ 30]*x11*x22*x31        
        +coeff[ 31]    *x21    *x43    
        +coeff[ 32]    *x23*x31*x42    
        +coeff[ 33]                *x51
        +coeff[ 34]        *x31*x41    
    ;
    v_x_r5p65_ep13_q1                         =v_x_r5p65_ep13_q1                         
        +coeff[ 35]*x11        *x42    
        +coeff[ 36]    *x22*x31        
        +coeff[ 37]*x12*x21    *x41    
        +coeff[ 38]    *x21*x32*x41    
        +coeff[ 39]    *x21*x31*x42    
        +coeff[ 40]    *x21*x31    *x53
        +coeff[ 41]    *x23    *x42    
        ;

    return v_x_r5p65_ep13_q1                         ;
}
float t_r5p65_ep13_q1                         (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.1773553E-03;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.29972E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 51]={
        -0.11559037E-03,-0.39681550E-02,-0.39451439E-02,-0.79321006E-04,
        -0.31063097E-03, 0.79251819E-04, 0.26461732E-03,-0.11146740E-02,
        -0.59756363E-03,-0.43330016E-02, 0.22922247E-02, 0.46570515E-03,
        -0.25308509E-02,-0.27861321E-03, 0.28045166E-02,-0.26941036E-05,
         0.43311364E-04,-0.15397141E-03, 0.16904048E-03, 0.41147461E-03,
         0.43480378E-03,-0.96678885E-03, 0.61980874E-03, 0.99377357E-03,
        -0.22936013E-03,-0.13679903E-03,-0.64619802E-04, 0.71273680E-04,
         0.95882846E-04,-0.86945380E-04, 0.99472672E-04,-0.86094384E-04,
        -0.54347234E-04, 0.21438798E-03, 0.10781989E-03,-0.17350816E-03,
        -0.28929111E-03, 0.17465092E-03, 0.21254798E-03,-0.20470254E-04,
        -0.34180694E-03,-0.23773877E-04, 0.25550309E-04,-0.29036810E-04,
        -0.65134911E-04, 0.56475674E-04, 0.61199418E-04,-0.19145479E-03,
        -0.14356128E-03, 0.47688416E-03,
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

//                 function

    float v_t_r5p65_ep13_q1                         =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]            *x41    
        +coeff[  5]*x11*x21            
        +coeff[  6]    *x22            
        +coeff[  7]    *x21*x31        
    ;
    v_t_r5p65_ep13_q1                         =v_t_r5p65_ep13_q1                         
        +coeff[  8]*x11        *x41    
        +coeff[  9]    *x21    *x41    
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]    *x23            
        +coeff[ 12]    *x21    *x42    
        +coeff[ 13]    *x21    *x41*x51
        +coeff[ 14]    *x23    *x41    
        +coeff[ 15]*x12*x21*x31*x41    
        +coeff[ 16]    *x21*x33*x41    
    ;
    v_t_r5p65_ep13_q1                         =v_t_r5p65_ep13_q1                         
        +coeff[ 17]*x11    *x31        
        +coeff[ 18]*x11            *x51
        +coeff[ 19]*x11*x22            
        +coeff[ 20]    *x22    *x41    
        +coeff[ 21]    *x21*x31*x41    
        +coeff[ 22]    *x23*x31        
        +coeff[ 23]*x11*x22    *x41    
        +coeff[ 24]*x11        *x42    
        +coeff[ 25]    *x21        *x52
    ;
    v_t_r5p65_ep13_q1                         =v_t_r5p65_ep13_q1                         
        +coeff[ 26]            *x42    
        +coeff[ 27]*x12*x21            
        +coeff[ 28]    *x22*x31        
        +coeff[ 29]    *x21*x32        
        +coeff[ 30]*x11*x21    *x41    
        +coeff[ 31]*x11    *x31*x41    
        +coeff[ 32]    *x21*x31    *x51
        +coeff[ 33]*x11*x22*x31        
        +coeff[ 34]*x12*x21    *x41    
    ;
    v_t_r5p65_ep13_q1                         =v_t_r5p65_ep13_q1                         
        +coeff[ 35]    *x21*x31*x42    
        +coeff[ 36]    *x21    *x43    
        +coeff[ 37]    *x23        *x51
        +coeff[ 38]    *x21    *x42*x51
        +coeff[ 39]        *x31*x41*x52
        +coeff[ 40]    *x23    *x42    
        +coeff[ 41]        *x31*x41    
        +coeff[ 42]*x11*x21*x31        
        +coeff[ 43]*x11        *x41*x51
    ;
    v_t_r5p65_ep13_q1                         =v_t_r5p65_ep13_q1                         
        +coeff[ 44]    *x21*x32*x41    
        +coeff[ 45]    *x21*x31*x41*x51
        +coeff[ 46]    *x21    *x41*x52
        +coeff[ 47]    *x23    *x41*x51
        +coeff[ 48]    *x21    *x43*x51
        +coeff[ 49]    *x23    *x43    
        ;

    return v_t_r5p65_ep13_q1                         ;
}
float y_r5p65_ep13_q1                         (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.2085432E-01;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.29972E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
         0.81276540E-02, 0.11713970E+00,-0.27363990E-02, 0.14366286E-01,
         0.12434795E-01, 0.59384950E-02,-0.19298708E-01,-0.53248378E-02,
        -0.10794300E-01,-0.16924401E-02, 0.30514768E-02,-0.29386231E-02,
        -0.33304414E-02,-0.82315067E-02,-0.40603839E-03,-0.11395684E-02,
        -0.68769825E-03,-0.95371203E-03, 0.10398955E-02,-0.40847003E-02,
         0.30428977E-02, 0.16480306E-02,-0.46752571E-03, 0.38770685E-03,
         0.13159076E-02,-0.99905883E-03,-0.40816268E-03,-0.47879206E-03,
         0.85703202E-03,-0.22554281E-02,-0.26850111E-03, 0.13490868E-03,
        -0.20373384E-03, 0.91198238E-03, 0.36879108E-03,-0.26983244E-03,
         0.32493708E-03,-0.11136365E-02,-0.79216814E-04, 0.18940870E-03,
        -0.25742009E-03,-0.11075842E-03,-0.44627851E-03,-0.14700467E-02,
        -0.12433368E-02, 0.38406925E-03, 0.10392684E-02, 0.73203834E-03,
        -0.12920801E-03,-0.53302338E-04,-0.49226321E-04, 0.57652240E-04,
        -0.55151086E-04, 0.41761951E-03, 0.20257405E-03, 0.35197736E-03,
        -0.10687909E-03,-0.12647035E-03, 0.14084105E-03,-0.23806846E-03,
         0.26458612E-03, 0.21397816E-02, 0.10425275E-02, 0.10977150E-02,
         0.21717716E-03, 0.32952783E-03, 0.18981045E-04, 0.47033976E-04,
         0.58394435E-04,-0.53111296E-04, 0.11727665E-03, 0.69622460E-04,
         0.11070457E-03, 0.35873477E-03, 0.48282200E-04,-0.61797378E-04,
        -0.37779240E-03,-0.24072772E-03, 0.30213027E-03,-0.22668055E-03,
         0.20093347E-03,-0.15826761E-03, 0.48352094E-03, 0.41680985E-04,
        -0.16876083E-04,-0.10877405E-03,-0.64060412E-04, 0.32324908E-04,
        -0.19260241E-03, 0.29171370E-04,-0.13909569E-03,-0.33060354E-04,
         0.48040729E-04,-0.38019196E-04,-0.60719580E-04, 0.12932418E-03,
         0.10973142E-03,
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

    float v_y_r5p65_ep13_q1                         =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]                *x51
        +coeff[  5]            *x42    
        +coeff[  6]    *x22            
        +coeff[  7]*x11*x21            
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[  8]    *x22    *x41    
        +coeff[  9]    *x21    *x41    
        +coeff[ 10]        *x31*x41    
        +coeff[ 11]    *x22*x31        
        +coeff[ 12]*x11*x21    *x41    
        +coeff[ 13]    *x22    *x42    
        +coeff[ 14]*x11                
        +coeff[ 15]            *x41*x51
        +coeff[ 16]                *x52
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 17]*x11*x21*x31        
        +coeff[ 18]    *x22        *x51
        +coeff[ 19]    *x22*x31*x41    
        +coeff[ 20]    *x24            
        +coeff[ 21]*x11*x23            
        +coeff[ 22]    *x21*x31        
        +coeff[ 23]        *x32        
        +coeff[ 24]            *x43    
        +coeff[ 25]    *x21    *x42    
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 26]*x12                
        +coeff[ 27]    *x21*x31*x41    
        +coeff[ 28]    *x23            
        +coeff[ 29]*x11*x21    *x42    
        +coeff[ 30]*x11        *x41    
        +coeff[ 31]    *x21        *x51
        +coeff[ 32]        *x31    *x51
        +coeff[ 33]        *x31*x42    
        +coeff[ 34]*x11*x22            
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 35]*x12        *x41    
        +coeff[ 36]*x11*x21        *x51
        +coeff[ 37]*x11*x21*x31*x41    
        +coeff[ 38]*x11    *x31        
        +coeff[ 39]        *x32*x41    
        +coeff[ 40]            *x42*x51
        +coeff[ 41]        *x31*x41*x51
        +coeff[ 42]    *x22*x32        
        +coeff[ 43]    *x22    *x43    
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 44]    *x22*x31*x42    
        +coeff[ 45]*x12*x22            
        +coeff[ 46]    *x24    *x41    
        +coeff[ 47]*x11*x23    *x41    
        +coeff[ 48]*x11        *x42    
        +coeff[ 49]    *x21*x32        
        +coeff[ 50]*x11    *x31*x41    
        +coeff[ 51]            *x41*x52
        +coeff[ 52]*x12    *x31        
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 53]    *x23    *x41    
        +coeff[ 54]*x11*x22    *x41    
        +coeff[ 55]    *x22    *x41*x51
        +coeff[ 56]*x12        *x42    
        +coeff[ 57]*x11*x21*x32        
        +coeff[ 58]*x11*x21    *x41*x51
        +coeff[ 59]    *x22*x32*x41    
        +coeff[ 60]*x11*x23*x31        
        +coeff[ 61]    *x24    *x42    
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 62]    *x24*x31*x41    
        +coeff[ 63]*x11*x23    *x42    
        +coeff[ 64]    *x24*x31*x42    
        +coeff[ 65]    *x22    *x42*x53
        +coeff[ 66]*x11            *x51
        +coeff[ 67]    *x21    *x41*x51
        +coeff[ 68]*x12*x21            
        +coeff[ 69]        *x31*x42*x51
        +coeff[ 70]    *x23*x31        
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 71]*x11*x22*x31        
        +coeff[ 72]    *x22*x31    *x51
        +coeff[ 73]    *x23    *x42    
        +coeff[ 74]*x11*x21*x31    *x51
        +coeff[ 75]    *x22        *x52
        +coeff[ 76]*x11*x21    *x43    
        +coeff[ 77]*x11*x21*x31*x42    
        +coeff[ 78]    *x24*x31        
        +coeff[ 79]    *x24        *x51
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 80]*x12*x22    *x41    
        +coeff[ 81]*x11*x23        *x51
        +coeff[ 82]*x11*x23*x31*x41    
        +coeff[ 83]            *x44    
        +coeff[ 84]        *x32    *x51
        +coeff[ 85]    *x21    *x43    
        +coeff[ 86]    *x21*x31*x42    
        +coeff[ 87]    *x21    *x42*x51
        +coeff[ 88]            *x45    
    ;
    v_y_r5p65_ep13_q1                         =v_y_r5p65_ep13_q1                         
        +coeff[ 89]    *x21*x31*x41*x51
        +coeff[ 90]        *x31*x44    
        +coeff[ 91]    *x23        *x51
        +coeff[ 92]*x12*x21    *x41    
        +coeff[ 93]*x12    *x31*x41    
        +coeff[ 94]        *x31*x43*x51
        +coeff[ 95]    *x23*x31*x41    
        +coeff[ 96]*x11*x22    *x42    
        ;

    return v_y_r5p65_ep13_q1                         ;
}
float p_r5p65_ep13_q1                         (float *x,int m){
    int ncoeff= 44;
    float avdat= -0.1169697E-01;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.29972E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 45]={
         0.51516276E-02,-0.18384609E-02, 0.56299628E-02, 0.64936489E-01,
         0.92817694E-02,-0.14276197E-01,-0.38925826E-02, 0.40458543E-02,
        -0.19697922E-02,-0.66858977E-02,-0.32979224E-04,-0.10693129E-02,
        -0.64314407E-03,-0.18543134E-02, 0.10015281E-02, 0.23273488E-02,
        -0.20635179E-02,-0.54623629E-02,-0.31468019E-03, 0.20781187E-02,
        -0.29794549E-03, 0.94494969E-03,-0.30829263E-03,-0.58059767E-03,
         0.25683157E-04,-0.24858150E-02, 0.13110320E-02, 0.66567992E-03,
        -0.14161897E-02, 0.25151257E-03,-0.15565800E-03,-0.64804725E-03,
         0.49855892E-03, 0.32565213E-03, 0.25152095E-03,-0.66332868E-03,
        -0.30034196E-03, 0.10884847E-03,-0.30154595E-03, 0.12386420E-03,
        -0.25098986E-03,-0.17207098E-03, 0.78306079E-03, 0.28303181E-03,
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
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_p_r5p65_ep13_q1                         =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]*x11*x21            
        +coeff[  7]            *x42    
    ;
    v_p_r5p65_ep13_q1                         =v_p_r5p65_ep13_q1                         
        +coeff[  8]            *x41*x51
        +coeff[  9]    *x22    *x41    
        +coeff[ 10]        *x33*x41    
        +coeff[ 11]    *x21    *x41    
        +coeff[ 12]                *x52
        +coeff[ 13]    *x22*x31        
        +coeff[ 14]    *x22        *x51
        +coeff[ 15]    *x24            
        +coeff[ 16]*x11*x21    *x41    
    ;
    v_p_r5p65_ep13_q1                         =v_p_r5p65_ep13_q1                         
        +coeff[ 17]    *x22    *x42    
        +coeff[ 18]*x11                
        +coeff[ 19]        *x31*x41    
        +coeff[ 20]        *x31    *x51
        +coeff[ 21]            *x43    
        +coeff[ 22]*x12                
        +coeff[ 23]*x11*x21*x31        
        +coeff[ 24]    *x23*x31        
        +coeff[ 25]    *x22*x31*x41    
    ;
    v_p_r5p65_ep13_q1                         =v_p_r5p65_ep13_q1                         
        +coeff[ 26]*x11*x23            
        +coeff[ 27]    *x25            
        +coeff[ 28]*x11*x21    *x42    
        +coeff[ 29]        *x32        
        +coeff[ 30]*x11        *x41    
        +coeff[ 31]    *x21    *x42    
        +coeff[ 32]        *x31*x42    
        +coeff[ 33]*x11*x22            
        +coeff[ 34]*x11*x21        *x51
    ;
    v_p_r5p65_ep13_q1                         =v_p_r5p65_ep13_q1                         
        +coeff[ 35]*x11*x21*x31*x41    
        +coeff[ 36]    *x21*x31        
        +coeff[ 37]    *x21        *x51
        +coeff[ 38]    *x21*x31*x41    
        +coeff[ 39]            *x41*x52
        +coeff[ 40]    *x22*x32        
        +coeff[ 41]*x12        *x41    
        +coeff[ 42]    *x24    *x41    
        +coeff[ 43]*x12*x22            
        ;

    return v_p_r5p65_ep13_q1                         ;
}
float l_r6_h90s_ep13_q1                       (float *x,int m){
    int ncoeff= 68;
    float avdat= -0.6216912E-03;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99630E-02,-0.29972E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 69]={
         0.10880487E-02, 0.54646523E-04,-0.12085693E-02,-0.49089510E-02,
        -0.27493085E-02,-0.24099537E-03,-0.22403996E-02,-0.51192613E-03,
         0.14145464E-02, 0.16142956E-03, 0.12563098E-03, 0.30716334E-03,
        -0.62310719E-04,-0.38917467E-04, 0.21756369E-03,-0.23310415E-03,
        -0.23739868E-03, 0.66165504E-03,-0.25780910E-04,-0.57372195E-04,
        -0.12613989E-03, 0.83699488E-04, 0.31764721E-04, 0.26008993E-03,
        -0.10925283E-03, 0.17359803E-03,-0.42010483E-03,-0.19743561E-03,
         0.25032388E-04, 0.15092628E-04, 0.18259427E-04, 0.81021608E-04,
         0.22867247E-04, 0.36466561E-04,-0.19425957E-04, 0.41568539E-04,
        -0.97355136E-04, 0.52912619E-04, 0.32112247E-03, 0.30551906E-04,
         0.18078142E-05, 0.22379041E-04, 0.18631143E-04,-0.21647653E-04,
        -0.33763132E-04,-0.21250225E-04,-0.94486728E-04, 0.15243563E-03,
        -0.44069468E-04, 0.78368612E-04,-0.10958896E-04,-0.16246242E-04,
         0.49584933E-05, 0.69706748E-05, 0.18260405E-04, 0.28940423E-04,
        -0.28926164E-04,-0.11959187E-04, 0.69456000E-05, 0.17170894E-04,
        -0.33850294E-04, 0.36700537E-05, 0.20941265E-04, 0.29804920E-04,
        -0.32790944E-04,-0.11567964E-03,-0.55815624E-04,-0.75780416E-04,
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
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_l_r6_h90s_ep13_q1                       =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]    *x22            
        +coeff[  5]        *x31*x41    
        +coeff[  6]            *x42    
        +coeff[  7]            *x41*x51
    ;
    v_l_r6_h90s_ep13_q1                       =v_l_r6_h90s_ep13_q1                       
        +coeff[  8]    *x22    *x41    
        +coeff[  9]    *x21    *x41    
        +coeff[ 10]*x11*x21            
        +coeff[ 11]*x11*x21    *x41    
        +coeff[ 12]                *x51
        +coeff[ 13]                *x52
        +coeff[ 14]    *x22*x31        
        +coeff[ 15]            *x43    
        +coeff[ 16]    *x24            
    ;
    v_l_r6_h90s_ep13_q1                       =v_l_r6_h90s_ep13_q1                       
        +coeff[ 17]    *x22    *x42    
        +coeff[ 18]        *x31    *x51
        +coeff[ 19]    *x23            
        +coeff[ 20]        *x31*x42    
        +coeff[ 21]    *x22        *x51
        +coeff[ 22]*x11*x21        *x51
        +coeff[ 23]    *x22*x31*x41    
        +coeff[ 24]*x11*x23            
        +coeff[ 25]*x11*x21    *x42    
    ;
    v_l_r6_h90s_ep13_q1                       =v_l_r6_h90s_ep13_q1                       
        +coeff[ 26]    *x24    *x41    
        +coeff[ 27]*x11*x23    *x41    
        +coeff[ 28]    *x21*x31        
        +coeff[ 29]    *x21        *x51
        +coeff[ 30]*x11        *x41    
        +coeff[ 31]    *x21    *x42    
        +coeff[ 32]            *x42*x51
        +coeff[ 33]            *x41*x52
        +coeff[ 34]*x11*x22            
    ;
    v_l_r6_h90s_ep13_q1                       =v_l_r6_h90s_ep13_q1                       
        +coeff[ 35]*x11*x21*x31        
        +coeff[ 36]    *x23    *x41    
        +coeff[ 37]*x11*x21*x31*x41    
        +coeff[ 38]    *x22    *x43    
        +coeff[ 39]    *x24*x31*x42    
        +coeff[ 40]*x11                
        +coeff[ 41]    *x21*x31*x41    
        +coeff[ 42]*x12        *x41    
        +coeff[ 43]    *x23*x31        
    ;
    v_l_r6_h90s_ep13_q1                       =v_l_r6_h90s_ep13_q1                       
        +coeff[ 44]*x11*x22    *x41    
        +coeff[ 45]*x12*x22            
        +coeff[ 46]    *x24*x31        
        +coeff[ 47]    *x22*x31*x42    
        +coeff[ 48]*x11*x23*x31        
        +coeff[ 49]*x11*x21    *x43    
        +coeff[ 50]        *x32        
        +coeff[ 51]        *x32*x41    
        +coeff[ 52]                *x53
    ;
    v_l_r6_h90s_ep13_q1                       =v_l_r6_h90s_ep13_q1                       
        +coeff[ 53]*x11        *x42    
        +coeff[ 54]    *x22*x32        
        +coeff[ 55]    *x21    *x43    
        +coeff[ 56]            *x44    
        +coeff[ 57]    *x22        *x52
        +coeff[ 58]*x12        *x42    
        +coeff[ 59]    *x22*x32*x41    
        +coeff[ 60]    *x23    *x42    
        +coeff[ 61]*x13            *x51
    ;
    v_l_r6_h90s_ep13_q1                       =v_l_r6_h90s_ep13_q1                       
        +coeff[ 62]    *x22    *x42*x51
        +coeff[ 63]*x11*x21*x31*x42    
        +coeff[ 64]*x12*x22    *x41    
        +coeff[ 65]    *x24    *x42    
        +coeff[ 66]    *x22*x31*x43    
        +coeff[ 67]*x11*x23    *x42    
        ;

    return v_l_r6_h90s_ep13_q1                       ;
}
float x_r5p65_ep14_q1ex                       (float *x,int m){
    int ncoeff= 44;
    float avdat= -0.1908322E-02;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 45]={
         0.46482793E-03, 0.26429130E-04, 0.12722461E+00, 0.30684951E-02,
        -0.17248025E-01, 0.32507207E-04, 0.36612926E-02,-0.25953900E-02,
        -0.44315909E-02, 0.51010568E-02,-0.12720716E-02, 0.12271157E-02,
         0.24419993E-02,-0.60899551E-02, 0.74265962E-02,-0.33184994E-03,
        -0.67880005E-03,-0.29463330E-02, 0.35744745E-03, 0.14522211E-02,
         0.31076083E-02, 0.17925790E-02,-0.12907966E-03, 0.24679414E-03,
        -0.34218637E-03, 0.36926530E-03, 0.36298414E-03,-0.40051673E-03,
        -0.41254173E-03, 0.75137540E-03,-0.23064812E-03,-0.21997979E-02,
         0.57020392E-04, 0.89521018E-04, 0.36423054E-03,-0.11683206E-03,
        -0.72840921E-03, 0.37049115E-03, 0.28090519E-03, 0.64994046E-03,
        -0.16196112E-02, 0.17620370E-03, 0.10977441E-02,-0.71791251E-03,
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

    float v_x_r5p65_ep14_q1ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]                *x51
        +coeff[  2]    *x21            
        +coeff[  3]    *x21        *x51
        +coeff[  4]    *x21    *x41    
        +coeff[  5]*x13                
        +coeff[  6]*x11                
        +coeff[  7]*x11        *x41    
    ;
    v_x_r5p65_ep14_q1ex                       =v_x_r5p65_ep14_q1ex                       
        +coeff[  8]    *x21*x31        
        +coeff[  9]    *x23            
        +coeff[ 10]            *x41    
        +coeff[ 11]    *x22            
        +coeff[ 12]*x11*x22            
        +coeff[ 13]    *x21    *x42    
        +coeff[ 14]    *x23    *x41    
        +coeff[ 15]        *x31        
        +coeff[ 16]*x11    *x31        
    ;
    v_x_r5p65_ep14_q1ex                       =v_x_r5p65_ep14_q1ex                       
        +coeff[ 17]    *x21*x31*x41    
        +coeff[ 18]*x11*x21            
        +coeff[ 19]    *x22    *x41    
        +coeff[ 20]*x11*x22    *x41    
        +coeff[ 21]    *x23*x31        
        +coeff[ 22]*x13        *x42    
        +coeff[ 23]*x11            *x51
        +coeff[ 24]            *x42    
        +coeff[ 25]*x12*x21            
    ;
    v_x_r5p65_ep14_q1ex                       =v_x_r5p65_ep14_q1ex                       
        +coeff[ 26]*x11*x21    *x41    
        +coeff[ 27]*x11    *x31*x41    
        +coeff[ 28]    *x21*x32        
        +coeff[ 29]*x11*x22*x31        
        +coeff[ 30]        *x33*x41    
        +coeff[ 31]    *x21    *x43    
        +coeff[ 32]    *x21    *x41*x53
        +coeff[ 33]    *x22*x33        
        +coeff[ 34]    *x23*x31*x42    
    ;
    v_x_r5p65_ep14_q1ex                       =v_x_r5p65_ep14_q1ex                       
        +coeff[ 35]    *x21        *x52
        +coeff[ 36]*x11        *x42    
        +coeff[ 37]    *x21    *x41*x51
        +coeff[ 38]    *x22*x31        
        +coeff[ 39]*x12*x21    *x41    
        +coeff[ 40]    *x21*x31*x42    
        +coeff[ 41]    *x21*x31    *x53
        +coeff[ 42]    *x23    *x42    
        +coeff[ 43]*x12*x21*x32*x41    
        ;

    return v_x_r5p65_ep14_q1ex                       ;
}
float t_r5p65_ep14_q1ex                       (float *x,int m){
    int ncoeff= 50;
    float avdat= -0.2675710E-03;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 51]={
        -0.42783745E-05,-0.40385923E-02,-0.55826460E-02,-0.29456150E-03,
         0.74556505E-04,-0.10794898E-02,-0.58689981E-03,-0.42904164E-02,
         0.23592981E-02, 0.92557707E-03,-0.15920254E-02, 0.24041752E-02,
         0.24339420E-03,-0.15257824E-03, 0.15528640E-03, 0.47566029E-03,
        -0.68753521E-03, 0.86651288E-03,-0.75102413E-04, 0.37362723E-03,
        -0.19470100E-03,-0.17292023E-03, 0.51112258E-03,-0.78060453E-04,
         0.74528230E-04, 0.89742483E-04,-0.80348887E-04, 0.92787770E-04,
        -0.87044580E-04,-0.11023203E-03, 0.19838366E-03,-0.43444548E-03,
        -0.78605764E-04,-0.33794084E-04, 0.10562970E-03,-0.29064209E-03,
         0.10448531E-03,-0.16522866E-03,-0.60604480E-05, 0.22111288E-04,
        -0.22116968E-04,-0.18616258E-04, 0.22168209E-04, 0.62611907E-04,
         0.36307451E-04,-0.11645705E-03, 0.50605745E-05,-0.46653868E-05,
        -0.11576993E-04, 0.81645085E-05,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_t_r5p65_ep14_q1ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]            *x41    
        +coeff[  4]*x11*x21            
        +coeff[  5]    *x21*x31        
        +coeff[  6]*x11        *x41    
        +coeff[  7]    *x21    *x41    
    ;
    v_t_r5p65_ep14_q1ex                       =v_t_r5p65_ep14_q1ex                       
        +coeff[  8]    *x21        *x51
        +coeff[  9]    *x23            
        +coeff[ 10]    *x21    *x42    
        +coeff[ 11]    *x23    *x41    
        +coeff[ 12]    *x22            
        +coeff[ 13]*x11    *x31        
        +coeff[ 14]*x11            *x51
        +coeff[ 15]*x11*x22            
        +coeff[ 16]    *x21*x31*x41    
    ;
    v_t_r5p65_ep14_q1ex                       =v_t_r5p65_ep14_q1ex                       
        +coeff[ 17]*x11*x22    *x41    
        +coeff[ 18]        *x31        
        +coeff[ 19]    *x22    *x41    
        +coeff[ 20]*x11        *x42    
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]    *x23*x31        
        +coeff[ 23]            *x42    
        +coeff[ 24]*x12*x21            
        +coeff[ 25]    *x22*x31        
    ;
    v_t_r5p65_ep14_q1ex                       =v_t_r5p65_ep14_q1ex                       
        +coeff[ 26]    *x21*x32        
        +coeff[ 27]*x11*x21    *x41    
        +coeff[ 28]*x11    *x31*x41    
        +coeff[ 29]    *x21        *x52
        +coeff[ 30]*x11*x22*x31        
        +coeff[ 31]    *x21    *x43    
        +coeff[ 32]    *x23*x31*x42    
        +coeff[ 33]        *x31*x41    
        +coeff[ 34]*x12*x21    *x41    
    ;
    v_t_r5p65_ep14_q1ex                       =v_t_r5p65_ep14_q1ex                       
        +coeff[ 35]    *x21*x31*x42    
        +coeff[ 36]    *x23        *x51
        +coeff[ 37]    *x23*x32*x41    
        +coeff[ 38]                *x51
        +coeff[ 39]*x11*x21*x31        
        +coeff[ 40]    *x21*x31    *x51
        +coeff[ 41]*x11        *x41*x51
        +coeff[ 42]*x12*x21*x31        
        +coeff[ 43]    *x22    *x42    
    ;
    v_t_r5p65_ep14_q1ex                       =v_t_r5p65_ep14_q1ex                       
        +coeff[ 44]*x11*x22        *x51
        +coeff[ 45]    *x23    *x42    
        +coeff[ 46]*x12                
        +coeff[ 47]        *x32        
        +coeff[ 48]*x11    *x32        
        +coeff[ 49]*x12        *x41    
        ;

    return v_t_r5p65_ep14_q1ex                       ;
}
float y_r5p65_ep14_q1ex                       (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.1311591E-01;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
        -0.29680489E-02, 0.15709081E+00,-0.38722430E-02, 0.17936099E-01,
         0.18164776E-01, 0.85489750E-02,-0.28146796E-01,-0.76750196E-02,
        -0.15185781E-01, 0.43587759E-02,-0.24278050E-02,-0.24002714E-02,
        -0.10771025E-02,-0.40632701E-02,-0.48919097E-02,-0.11730817E-01,
        -0.57758659E-03,-0.66081778E-03,-0.58402563E-03, 0.12229452E-02,
        -0.13833115E-02, 0.15906041E-02,-0.57483800E-02, 0.45474903E-02,
         0.23892778E-02, 0.55266079E-03,-0.37754903E-03,-0.40615047E-03,
        -0.14223468E-02, 0.54190366E-03,-0.33022081E-02,-0.42435949E-03,
        -0.16913034E-02,-0.71696850E-03, 0.17965045E-03,-0.68127154E-03,
        -0.37909765E-03, 0.48155067E-03,-0.27752307E-03, 0.54084306E-03,
        -0.10745190E-03,-0.15934234E-03, 0.66770770E-03,-0.61393320E-03,
         0.59689820E-03, 0.13467909E-02, 0.11041559E-02, 0.20745329E-02,
         0.16182987E-02, 0.30406733E-04,-0.40771475E-03,-0.84699190E-04,
        -0.76046890E-04,-0.16510709E-03, 0.11964135E-03,-0.79003097E-04,
         0.45867182E-04, 0.35800799E-03,-0.14296171E-03,-0.19834303E-03,
         0.19304742E-03, 0.19275128E-03,-0.44093374E-02, 0.36815493E-03,
         0.29519462E-03,-0.40148853E-03, 0.34792619E-03, 0.14110403E-02,
         0.13236257E-02, 0.23199902E-02, 0.35953167E-03, 0.49882619E-04,
        -0.54709235E-04, 0.89224857E-04, 0.10681547E-03, 0.17420223E-03,
         0.79542704E-04,-0.35683415E-02, 0.40769219E-03,-0.10810664E-03,
        -0.12757306E-03, 0.15319693E-03,-0.88539190E-03,-0.57854300E-04,
        -0.31834855E-03, 0.58860405E-04, 0.25691587E-03, 0.22627660E-02,
        -0.19950267E-03, 0.16608707E-02, 0.83085371E-03, 0.37383034E-02,
        -0.58099360E-03, 0.23437063E-02,-0.44937863E-03, 0.14949327E-02,
         0.50881493E-03,
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

    float v_y_r5p65_ep14_q1ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]                *x51
        +coeff[  5]            *x42    
        +coeff[  6]    *x22            
        +coeff[  7]*x11*x21            
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[  8]    *x22    *x41    
        +coeff[  9]        *x31*x41    
        +coeff[ 10]    *x21    *x41    
        +coeff[ 11]            *x41*x51
        +coeff[ 12]                *x52
        +coeff[ 13]    *x22*x31        
        +coeff[ 14]*x11*x21    *x41    
        +coeff[ 15]    *x22    *x42    
        +coeff[ 16]*x11                
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 17]    *x21*x31        
        +coeff[ 18]*x12                
        +coeff[ 19]    *x23            
        +coeff[ 20]*x11*x21*x31        
        +coeff[ 21]    *x22        *x51
        +coeff[ 22]    *x22*x31*x41    
        +coeff[ 23]    *x24            
        +coeff[ 24]*x11*x23            
        +coeff[ 25]        *x32        
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 26]*x11        *x41    
        +coeff[ 27]        *x31    *x51
        +coeff[ 28]    *x21    *x42    
        +coeff[ 29]*x11*x22            
        +coeff[ 30]*x11*x21    *x42    
        +coeff[ 31]            *x45    
        +coeff[ 32]*x11*x21*x31*x41    
        +coeff[ 33]        *x31*x44    
        +coeff[ 34]    *x21        *x51
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 35]    *x21*x31*x41    
        +coeff[ 36]*x12        *x41    
        +coeff[ 37]*x11*x21        *x51
        +coeff[ 38]        *x32*x43    
        +coeff[ 39]*x12*x22            
        +coeff[ 40]*x11    *x31        
        +coeff[ 41]*x11        *x42    
        +coeff[ 42]    *x23    *x41    
        +coeff[ 43]    *x22*x32        
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 44]    *x22    *x41*x51
        +coeff[ 45]    *x24    *x41    
        +coeff[ 46]*x11*x23    *x41    
        +coeff[ 47]            *x43    
        +coeff[ 48]        *x31*x42    
        +coeff[ 49]*x11            *x51
        +coeff[ 50]            *x42*x51
        +coeff[ 51]    *x21*x32        
        +coeff[ 52]*x11    *x31*x41    
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 53]        *x31*x41*x51
        +coeff[ 54]            *x41*x52
        +coeff[ 55]*x12    *x31        
        +coeff[ 56]                *x53
        +coeff[ 57]*x11*x22    *x41    
        +coeff[ 58]*x12        *x42    
        +coeff[ 59]*x11*x21*x32        
        +coeff[ 60]    *x22*x31    *x51
        +coeff[ 61]*x11*x21    *x41*x51
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 62]    *x22*x31*x42    
        +coeff[ 63]    *x22    *x42*x51
        +coeff[ 64]    *x24*x31        
        +coeff[ 65]    *x22    *x44    
        +coeff[ 66]*x11*x23*x31        
        +coeff[ 67]    *x24*x31*x41    
        +coeff[ 68]    *x22    *x45    
        +coeff[ 69]    *x24    *x44    
        +coeff[ 70]        *x32*x41    
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 71]    *x21    *x41*x51
        +coeff[ 72]    *x21*x31*x42    
        +coeff[ 73]*x12*x21            
        +coeff[ 74]    *x21    *x42*x51
        +coeff[ 75]    *x23*x31        
        +coeff[ 76]*x11*x22*x31        
        +coeff[ 77]    *x22    *x43    
        +coeff[ 78]    *x23    *x42    
        +coeff[ 79]    *x22        *x52
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 80]        *x31*x43*x51
        +coeff[ 81]    *x23*x31*x41    
        +coeff[ 82]    *x22*x32*x41    
        +coeff[ 83]*x11*x21        *x52
        +coeff[ 84]    *x24        *x51
        +coeff[ 85]    *x23*x31    *x51
        +coeff[ 86]*x12*x22    *x41    
        +coeff[ 87]    *x24    *x42    
        +coeff[ 88]*x11*x23        *x51
    ;
    v_y_r5p65_ep14_q1ex                       =v_y_r5p65_ep14_q1ex                       
        +coeff[ 89]*x11*x23    *x42    
        +coeff[ 90]*x11*x23*x31*x41    
        +coeff[ 91]    *x22*x31*x44    
        +coeff[ 92]*x11*x21    *x45    
        +coeff[ 93]    *x22*x32*x43    
        +coeff[ 94]*x11*x21*x31*x44    
        +coeff[ 95]    *x24*x31*x42    
        +coeff[ 96]    *x22*x33*x42    
        ;

    return v_y_r5p65_ep14_q1ex                       ;
}
float p_r5p65_ep14_q1ex                       (float *x,int m){
    int ncoeff= 46;
    float avdat= -0.5558616E-02;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 47]={
        -0.84786135E-03,-0.19723768E-02, 0.58207284E-02, 0.66207819E-01,
         0.94096307E-02,-0.14404127E-01,-0.20155488E-02,-0.39678407E-02,
        -0.80929827E-02, 0.42186133E-02,-0.22256933E-02,-0.24427006E-05,
         0.32895184E-03,-0.31166457E-03,-0.10502219E-02, 0.21540702E-02,
        -0.65456110E-03, 0.25267797E-02,-0.54912460E-02,-0.29183389E-03,
        -0.30082819E-03, 0.69331791E-03,-0.22362010E-02, 0.91127807E-03,
        -0.31492667E-03,-0.60126488E-03,-0.26013516E-02, 0.14109304E-02,
         0.28731933E-03,-0.15511984E-03,-0.60817320E-03, 0.76072861E-03,
         0.30881618E-03, 0.23361534E-03,-0.13736141E-02, 0.79950398E-04,
         0.10781730E-03, 0.10640471E-03,-0.27598103E-03, 0.49475057E-03,
        -0.29252397E-03,-0.18393878E-03, 0.94926916E-03,-0.69608114E-03,
         0.31887155E-03, 0.16685751E-03,
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

    float v_p_r5p65_ep14_q1ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]            *x41*x51
        +coeff[  7]*x11*x21            
    ;
    v_p_r5p65_ep14_q1ex                       =v_p_r5p65_ep14_q1ex                       
        +coeff[  8]    *x22    *x41    
        +coeff[  9]            *x42    
        +coeff[ 10]*x11*x21    *x41    
        +coeff[ 11]        *x33*x41    
        +coeff[ 12]    *x24*x31        
        +coeff[ 13]*x11                
        +coeff[ 14]    *x21    *x41    
        +coeff[ 15]        *x31*x41    
        +coeff[ 16]                *x52
    ;
    v_p_r5p65_ep14_q1ex                       =v_p_r5p65_ep14_q1ex                       
        +coeff[ 17]    *x24            
        +coeff[ 18]    *x22    *x42    
        +coeff[ 19]    *x21*x31        
        +coeff[ 20]        *x31    *x51
        +coeff[ 21]    *x23            
        +coeff[ 22]    *x22*x31        
        +coeff[ 23]    *x22        *x51
        +coeff[ 24]*x12                
        +coeff[ 25]*x11*x21*x31        
    ;
    v_p_r5p65_ep14_q1ex                       =v_p_r5p65_ep14_q1ex                       
        +coeff[ 26]    *x22*x31*x41    
        +coeff[ 27]*x11*x23            
        +coeff[ 28]        *x32        
        +coeff[ 29]*x11        *x41    
        +coeff[ 30]    *x21    *x42    
        +coeff[ 31]            *x43    
        +coeff[ 32]*x11*x22            
        +coeff[ 33]*x11*x21        *x51
        +coeff[ 34]*x11*x21    *x42    
    ;
    v_p_r5p65_ep14_q1ex                       =v_p_r5p65_ep14_q1ex                       
        +coeff[ 35]        *x33*x42    
        +coeff[ 36]*x11*x23*x31*x41    
        +coeff[ 37]    *x21        *x51
        +coeff[ 38]    *x21*x31*x41    
        +coeff[ 39]        *x31*x42    
        +coeff[ 40]    *x22*x32        
        +coeff[ 41]*x12        *x41    
        +coeff[ 42]    *x24    *x41    
        +coeff[ 43]*x11*x21*x31*x41    
    ;
    v_p_r5p65_ep14_q1ex                       =v_p_r5p65_ep14_q1ex                       
        +coeff[ 44]*x12*x22            
        +coeff[ 45]        *x34*x41    
        ;

    return v_p_r5p65_ep14_q1ex                       ;
}
float l_r6_h90s_ep14_q1ex                     (float *x,int m){
    int ncoeff= 82;
    float avdat= -0.1248902E-02;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 83]={
         0.16775030E-02,-0.11897199E-02,-0.46404316E-02,-0.29157307E-04,
        -0.28274965E-02,-0.46467796E-03,-0.35454063E-02,-0.90380339E-03,
         0.18826526E-02, 0.44033828E-03, 0.62267631E-04, 0.94204406E-04,
         0.24683599E-03,-0.26089908E-03, 0.10599403E-02, 0.52516167E-04,
         0.22673818E-03,-0.63089654E-04,-0.68696558E-04,-0.60488506E-04,
        -0.42977137E-03, 0.17596067E-03, 0.61681341E-04,-0.11806922E-03,
        -0.49424794E-03, 0.28420798E-05,-0.14332513E-03, 0.24363369E-04,
         0.25518364E-04, 0.29320865E-04, 0.13962261E-03,-0.23682837E-03,
         0.87359644E-04, 0.72682953E-04, 0.47826128E-04,-0.12564629E-03,
         0.42691923E-03, 0.28886425E-03, 0.68001007E-03,-0.25548576E-03,
        -0.23053068E-04, 0.45356785E-04,-0.38697035E-04,-0.21807797E-04,
         0.10849835E-03,-0.29779389E-04,-0.10190898E-03, 0.32522224E-03,
         0.45616041E-04,-0.64452490E-06, 0.42471870E-05, 0.95916603E-05,
         0.14737376E-04,-0.22879527E-04, 0.29700575E-04,-0.68131580E-04,
        -0.23960265E-04,-0.53196363E-04,-0.48119160E-04, 0.14201875E-03,
        -0.30583510E-03, 0.14050333E-05,-0.20512851E-03, 0.72106104E-05,
         0.88529705E-05, 0.50802323E-05, 0.45907286E-05, 0.33328160E-04,
        -0.44753560E-04, 0.14743877E-04, 0.23726063E-04, 0.73948854E-05,
         0.47863537E-04,-0.82387218E-04, 0.62490042E-04,-0.97169377E-05,
         0.85160191E-05, 0.23467430E-04,-0.76370467E-04,-0.42764361E-04,
        -0.34795139E-04,-0.15713302E-03,
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

    float v_l_r6_h90s_ep14_q1ex                     =avdat
        +coeff[  0]                    
        +coeff[  1]        *x31        
        +coeff[  2]            *x41    
        +coeff[  3]                *x51
        +coeff[  4]    *x22            
        +coeff[  5]        *x31*x41    
        +coeff[  6]            *x42    
        +coeff[  7]            *x41*x51
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[  8]    *x22    *x41    
        +coeff[  9]*x11*x21    *x41    
        +coeff[ 10]    *x21    *x43    
        +coeff[ 11]*x11*x21            
        +coeff[ 12]    *x22*x31        
        +coeff[ 13]    *x24            
        +coeff[ 14]    *x22    *x42    
        +coeff[ 15]    *x21            
        +coeff[ 16]    *x21    *x41    
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 17]        *x31    *x51
        +coeff[ 18]                *x52
        +coeff[ 19]    *x23            
        +coeff[ 20]            *x43    
        +coeff[ 21]    *x22        *x51
        +coeff[ 22]*x11*x21        *x51
        +coeff[ 23]*x11*x23            
        +coeff[ 24]    *x24    *x41    
        +coeff[ 25]        *x31*x44    
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 26]    *x24*x31*x41    
        +coeff[ 27]    *x21*x31        
        +coeff[ 28]    *x21        *x51
        +coeff[ 29]*x11        *x41    
        +coeff[ 30]    *x21    *x42    
        +coeff[ 31]        *x31*x42    
        +coeff[ 32]            *x42*x51
        +coeff[ 33]            *x41*x52
        +coeff[ 34]*x11*x21*x31        
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 35]    *x23    *x41    
        +coeff[ 36]    *x22*x31*x41    
        +coeff[ 37]*x11*x21    *x42    
        +coeff[ 38]    *x22    *x43    
        +coeff[ 39]*x11*x23    *x41    
        +coeff[ 40]        *x32        
        +coeff[ 41]    *x21*x31*x41    
        +coeff[ 42]        *x32*x41    
        +coeff[ 43]*x11*x22            
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 44]*x11*x21*x31*x41    
        +coeff[ 45]*x12*x22            
        +coeff[ 46]    *x24*x31        
        +coeff[ 47]    *x22*x31*x42    
        +coeff[ 48]*x12        *x43    
        +coeff[ 49]*x11                
        +coeff[ 50]*x11            *x51
        +coeff[ 51]        *x31*x41*x51
        +coeff[ 52]*x11        *x42    
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 53]    *x23*x31        
        +coeff[ 54]    *x22*x32        
        +coeff[ 55]            *x44    
        +coeff[ 56]    *x22        *x52
        +coeff[ 57]*x11*x22    *x41    
        +coeff[ 58]*x11*x23*x31        
        +coeff[ 59]*x11*x21    *x43    
        +coeff[ 60]    *x24    *x42    
        +coeff[ 61]    *x22*x31*x43    
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 62]*x11*x23    *x42    
        +coeff[ 63]        *x31    *x52
        +coeff[ 64]                *x53
        +coeff[ 65]*x11    *x31*x41    
        +coeff[ 66]*x12            *x51
        +coeff[ 67]    *x21*x31*x42    
        +coeff[ 68]        *x31*x43    
        +coeff[ 69]    *x22*x31    *x51
        +coeff[ 70]            *x43*x51
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 71]*x11*x21*x31    *x51
        +coeff[ 72]    *x22*x32*x41    
        +coeff[ 73]    *x23    *x42    
        +coeff[ 74]*x11*x21*x31*x42    
        +coeff[ 75]*x12*x21        *x52
        +coeff[ 76]    *x22*x31*x42*x51
        +coeff[ 77]        *x31*x44*x51
        +coeff[ 78]*x11*x23*x31*x41    
        +coeff[ 79]*x12*x22    *x42    
    ;
    v_l_r6_h90s_ep14_q1ex                     =v_l_r6_h90s_ep14_q1ex                     
        +coeff[ 80]    *x23*x33*x41    
        +coeff[ 81]    *x24    *x43    
        ;

    return v_l_r6_h90s_ep14_q1ex                     ;
}
float x_r5p65_ep23_den                        (float *x,int m){
    int ncoeff= 46;
    float avdat= -0.6068480E+01;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 47]={
         0.33565375E-03,-0.60099359E-02, 0.46765803E-04, 0.84761113E-01,
         0.42115902E-02,-0.22434514E-01, 0.58970931E-02, 0.11507285E-03,
         0.25383115E-02,-0.16086991E-02,-0.32419260E-02, 0.14274145E-02,
        -0.57531064E-02, 0.28780475E-02,-0.86635100E-02, 0.11428011E-01,
        -0.84745360E-03,-0.38679310E-02,-0.42126433E-03, 0.57306525E-03,
         0.44145895E-03, 0.19251268E-02, 0.43169279E-02,-0.45789336E-03,
         0.45311209E-03,-0.34606439E-03,-0.10359816E-02,-0.46050109E-03,
         0.10219528E-02,-0.23579253E-02,-0.10678021E-03,-0.58243860E-03,
        -0.19805600E-03, 0.50065742E-03,-0.41508250E-03, 0.44480758E-03,
         0.48554179E-03, 0.53808664E-03,-0.15818329E-02, 0.23044937E-03,
        -0.93647587E-03, 0.48106260E-04, 0.10928881E-03, 0.33389009E-03,
         0.33585980E-03,-0.52526034E-03,
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

    float v_x_r5p65_ep23_den                        =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]                *x51
        +coeff[  3]    *x21            
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x21    *x41    
        +coeff[  6]    *x23            
        +coeff[  7]*x12*x21*x31        
    ;
    v_x_r5p65_ep23_den                        =v_x_r5p65_ep23_den                        
        +coeff[  8]    *x23*x31        
        +coeff[  9]            *x41    
        +coeff[ 10]*x11        *x41    
        +coeff[ 11]    *x22            
        +coeff[ 12]    *x21*x31        
        +coeff[ 13]*x11*x22            
        +coeff[ 14]    *x21    *x42    
        +coeff[ 15]    *x23    *x41    
        +coeff[ 16]*x11    *x31        
    ;
    v_x_r5p65_ep23_den                        =v_x_r5p65_ep23_den                        
        +coeff[ 17]    *x21*x31*x41    
        +coeff[ 18]        *x31        
        +coeff[ 19]*x11            *x51
        +coeff[ 20]*x11*x21            
        +coeff[ 21]    *x22    *x41    
        +coeff[ 22]*x11*x22    *x41    
        +coeff[ 23]            *x42    
        +coeff[ 24]*x12*x21            
        +coeff[ 25]    *x21        *x52
    ;
    v_x_r5p65_ep23_den                        =v_x_r5p65_ep23_den                        
        +coeff[ 26]*x11        *x42    
        +coeff[ 27]    *x21*x32        
        +coeff[ 28]*x11*x22*x31        
        +coeff[ 29]    *x21    *x43    
        +coeff[ 30]*x13    *x31*x41    
        +coeff[ 31]    *x23*x31*x42    
        +coeff[ 32]        *x31*x41    
        +coeff[ 33]*x11*x21    *x41    
        +coeff[ 34]*x11    *x31*x41    
    ;
    v_x_r5p65_ep23_den                        =v_x_r5p65_ep23_den                        
        +coeff[ 35]    *x21    *x41*x51
        +coeff[ 36]    *x22*x31        
        +coeff[ 37]*x12*x21    *x41    
        +coeff[ 38]    *x21*x31*x42    
        +coeff[ 39]    *x21*x31    *x53
        +coeff[ 40]    *x23*x32*x41    
        +coeff[ 41]            *x41*x51
        +coeff[ 42]*x11*x21*x31        
        +coeff[ 43]    *x21    *x42*x51
    ;
    v_x_r5p65_ep23_den                        =v_x_r5p65_ep23_den                        
        +coeff[ 44]    *x22    *x42    
        +coeff[ 45]    *x23    *x41*x51
        ;

    return v_x_r5p65_ep23_den                        ;
}
float t_r5p65_ep23_den                        (float *x,int m){
    int ncoeff= 80;
    float avdat=  0.3747658E+01;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 81]={
        -0.32241516E-01, 0.71682400E+00,-0.73084205E-01, 0.14536658E+00,
        -0.22284898E+00, 0.33189919E-01,-0.29987996E-01, 0.90344600E-01,
        -0.57535321E-01, 0.21897595E-01,-0.41428567E-02,-0.15584398E-01,
        -0.21441124E-01,-0.57482995E-01, 0.24107844E-01,-0.90570174E-01,
         0.10145764E+00, 0.29531235E-02, 0.59210381E-03,-0.81543643E-02,
         0.60573947E-02,-0.15400211E-01,-0.40566668E-01, 0.24547547E-01,
         0.43129515E-01, 0.13220662E-02,-0.37117517E-02, 0.59093814E-02,
        -0.10006613E-01, 0.19193963E-02,-0.21599797E-02, 0.11998849E-01,
        -0.47993218E-02, 0.39054500E-02,-0.38690867E-02,-0.49948925E-02,
         0.10621985E-01,-0.22504613E-01,-0.22197828E-01, 0.26571080E-01,
        -0.19473692E-02,-0.15706801E-02, 0.14838125E-02, 0.57317461E-02,
         0.45293649E-02, 0.50994367E-02,-0.16983494E-01, 0.56654084E-02,
        -0.19922375E-02, 0.48324819E-04,-0.27902337E-03, 0.25894344E-02,
        -0.67900232E-03, 0.16173406E-02, 0.78291015E-03, 0.49037836E-02,
        -0.78749256E-02,-0.46672244E-02, 0.36659234E-02, 0.54545980E-02,
        -0.42348756E-02, 0.43537188E-02,-0.93916105E-02, 0.70437086E-02,
         0.10203967E-02,-0.19708932E-02, 0.15168227E-03, 0.44836191E-03,
        -0.70475804E-03,-0.30272503E-03,-0.21371748E-02,-0.12967580E-02,
         0.13050828E-02, 0.12503869E-02, 0.22994513E-02,-0.16905166E-02,
         0.19227387E-03, 0.57758525E-03,-0.23221935E-02, 0.25120860E-02,
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

    float v_t_r5p65_ep23_den                        =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]*x11                
        +coeff[  3]    *x22            
        +coeff[  4]    *x21    *x41    
        +coeff[  5]    *x21        *x51
        +coeff[  6]*x11        *x41    
        +coeff[  7]    *x23            
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[  8]    *x22    *x41    
        +coeff[  9]    *x23*x31        
        +coeff[ 10]        *x31        
        +coeff[ 11]            *x41    
        +coeff[ 12]*x11*x21            
        +coeff[ 13]    *x21*x31        
        +coeff[ 14]*x11*x22            
        +coeff[ 15]    *x21    *x42    
        +coeff[ 16]    *x23    *x41    
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[ 17]    *x23*x31*x41    
        +coeff[ 18]                *x51
        +coeff[ 19]*x11    *x31        
        +coeff[ 20]*x11            *x51
        +coeff[ 21]    *x22*x31        
        +coeff[ 22]    *x21*x31*x41    
        +coeff[ 23]    *x24            
        +coeff[ 24]*x11*x22    *x41    
        +coeff[ 25]    *x24        *x51
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[ 26]            *x42    
        +coeff[ 27]*x12*x21            
        +coeff[ 28]*x11        *x42    
        +coeff[ 29]*x12                
        +coeff[ 30]        *x31*x41    
        +coeff[ 31]    *x22        *x51
        +coeff[ 32]    *x21*x32        
        +coeff[ 33]*x11*x21    *x41    
        +coeff[ 34]    *x21        *x52
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[ 35]*x11    *x31*x41    
        +coeff[ 36]*x11*x22*x31        
        +coeff[ 37]    *x22    *x42    
        +coeff[ 38]    *x21    *x43    
        +coeff[ 39]    *x24    *x41    
        +coeff[ 40]    *x24*x31*x41    
        +coeff[ 41]    *x23*x31*x42    
        +coeff[ 42]*x11*x21*x31        
        +coeff[ 43]    *x21    *x41*x51
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[ 44]*x11*x23            
        +coeff[ 45]*x12*x21    *x41    
        +coeff[ 46]    *x21*x31*x42    
        +coeff[ 47]    *x24*x31        
        +coeff[ 48]    *x23*x31    *x51
        +coeff[ 49]    *x23*x32*x41    
        +coeff[ 50]        *x32        
        +coeff[ 51]    *x21*x31    *x51
        +coeff[ 52]*x11    *x32        
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[ 53]*x12        *x41    
        +coeff[ 54]*x11        *x41*x51
        +coeff[ 55]    *x23        *x51
        +coeff[ 56]    *x22*x31*x41    
        +coeff[ 57]    *x21*x32*x41    
        +coeff[ 58]*x11*x21    *x42    
        +coeff[ 59]    *x21    *x42*x51
        +coeff[ 60]*x11*x24            
        +coeff[ 61]*x11*x23    *x41    
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[ 62]    *x23    *x41*x51
        +coeff[ 63]*x11*x22    *x42    
        +coeff[ 64]*x11*x21        *x53
        +coeff[ 65]    *x24*x32        
        +coeff[ 66]                *x52
        +coeff[ 67]*x12    *x31        
        +coeff[ 68]            *x43    
        +coeff[ 69]*x11            *x52
        +coeff[ 70]    *x22    *x41*x51
    ;
    v_t_r5p65_ep23_den                        =v_t_r5p65_ep23_den                        
        +coeff[ 71]    *x22        *x52
        +coeff[ 72]*x12*x21*x31        
        +coeff[ 73]*x11*x21*x31*x41    
        +coeff[ 74]    *x21*x31*x41*x51
        +coeff[ 75]*x11        *x43    
        +coeff[ 76]        *x32*x41*x51
        +coeff[ 77]            *x41*x53
        +coeff[ 78]*x12*x23            
        +coeff[ 79]*x11*x22*x31*x41    
        ;

    return v_t_r5p65_ep23_den                        ;
}
float y_r5p65_ep23_den                        (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.1276265E-01;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
        -0.13407141E-02, 0.14985606E+00,-0.48283576E-02, 0.10909664E-01,
         0.21989416E-01, 0.10955148E-01,-0.35526115E-01, 0.11606515E-01,
        -0.95804520E-02,-0.20381924E-01, 0.54970505E-02,-0.51930174E-02,
        -0.58714221E-02,-0.69691107E-03,-0.27963063E-02,-0.12527249E-01,
        -0.77838008E-03, 0.67565159E-03, 0.93194394E-03,-0.16122537E-02,
        -0.71527436E-03, 0.61771547E-03, 0.16335084E-02,-0.16496085E-02,
        -0.12662022E-02,-0.66292030E-02,-0.37993726E-02, 0.65492475E-02,
         0.31531996E-02,-0.43573993E-03,-0.83499745E-03, 0.64863078E-03,
        -0.90521050E-03,-0.19542635E-02,-0.11860563E-03,-0.11283578E-03,
         0.99302141E-03,-0.31179207E-03,-0.45030416E-03,-0.24475204E-03,
        -0.72267198E-03, 0.68466092E-03, 0.28333606E-03,-0.20353262E-03,
         0.61096845E-03, 0.18116490E-03,-0.15696987E-03,-0.27171057E-03,
        -0.33433496E-02, 0.44637713E-04, 0.10953705E-02, 0.51794814E-04,
         0.36266526E-05, 0.46464347E-03,-0.91093360E-04,-0.82949373E-04,
        -0.11848044E-03,-0.21983789E-04, 0.27090372E-03,-0.19346102E-03,
        -0.14334705E-03,-0.18880518E-03,-0.37374760E-02, 0.20241074E-03,
        -0.91017544E-03, 0.27252999E-03, 0.38158856E-03, 0.28449327E-02,
        -0.62178056E-04, 0.15812636E-02, 0.17906755E-02, 0.10832681E-03,
         0.43697856E-03,-0.23654637E-04,-0.12354717E-04, 0.13054804E-03,
        -0.98391814E-04,-0.10733956E-03, 0.16607578E-03,-0.26566317E-03,
        -0.88474575E-04, 0.52408088E-03,-0.60259223E-04,-0.66962966E-03,
         0.76698037E-04,-0.21581849E-03,-0.29436708E-03, 0.77783637E-03,
         0.22758914E-02, 0.15818581E-02, 0.16176722E-02, 0.54483843E-03,
         0.29582743E-03, 0.37654009E-03, 0.26665451E-02,-0.30609523E-03,
        -0.11783203E-03,
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

    float v_y_r5p65_ep23_den                        =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]                *x51
        +coeff[  5]            *x42    
        +coeff[  6]    *x22            
        +coeff[  7]            *x41*x51
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[  8]*x11*x21            
        +coeff[  9]    *x22    *x41    
        +coeff[ 10]        *x31*x41    
        +coeff[ 11]    *x22*x31        
        +coeff[ 12]*x11*x21    *x41    
        +coeff[ 13]*x11                
        +coeff[ 14]    *x21    *x41    
        +coeff[ 15]    *x22    *x42    
        +coeff[ 16]    *x21*x31        
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 17]        *x32        
        +coeff[ 18]        *x31    *x51
        +coeff[ 19]    *x21    *x42    
        +coeff[ 20]*x12                
        +coeff[ 21]                *x52
        +coeff[ 22]    *x23            
        +coeff[ 23]*x11*x21*x31        
        +coeff[ 24]    *x22        *x51
        +coeff[ 25]    *x22*x31*x41    
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 26]*x11*x21    *x42    
        +coeff[ 27]    *x24            
        +coeff[ 28]*x11*x23            
        +coeff[ 29]*x11        *x41    
        +coeff[ 30]    *x21*x31*x41    
        +coeff[ 31]*x11*x22            
        +coeff[ 32]            *x41*x52
        +coeff[ 33]*x11*x21*x31*x41    
        +coeff[ 34]*x11    *x31        
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 35]    *x21        *x51
        +coeff[ 36]        *x31*x42    
        +coeff[ 37]            *x42*x51
        +coeff[ 38]*x12        *x41    
        +coeff[ 39]*x11*x21        *x51
        +coeff[ 40]    *x22*x32        
        +coeff[ 41]*x12*x22            
        +coeff[ 42]        *x32*x41    
        +coeff[ 43]*x11        *x42    
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 44]    *x23    *x41    
        +coeff[ 45]    *x23*x31        
        +coeff[ 46]                *x53
        +coeff[ 47]    *x22    *x41*x51
        +coeff[ 48]    *x22    *x43    
        +coeff[ 49]    *x24    *x41    
        +coeff[ 50]*x11*x23    *x41    
        +coeff[ 51]*x11*x23*x32        
        +coeff[ 52]*x12    *x31    *x52
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 53]            *x43    
        +coeff[ 54]    *x21*x32        
        +coeff[ 55]*x11    *x31*x41    
        +coeff[ 56]    *x21    *x41*x51
        +coeff[ 57]        *x31*x41*x51
        +coeff[ 58]*x11*x22    *x41    
        +coeff[ 59]*x12        *x42    
        +coeff[ 60]    *x22*x31    *x51
        +coeff[ 61]*x11*x21    *x41*x51
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 62]    *x22*x31*x42    
        +coeff[ 63]    *x22        *x52
        +coeff[ 64]    *x22    *x42*x51
        +coeff[ 65]    *x24*x31        
        +coeff[ 66]*x11*x23*x31        
        +coeff[ 67]    *x24    *x42    
        +coeff[ 68]*x12*x23            
        +coeff[ 69]    *x24*x31*x41    
        +coeff[ 70]*x11*x23    *x42    
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 71]*x11*x22*x33        
        +coeff[ 72]*x12*x24    *x41    
        +coeff[ 73]*x11            *x51
        +coeff[ 74]    *x21*x31    *x51
        +coeff[ 75]*x12*x21            
        +coeff[ 76]*x12    *x31        
        +coeff[ 77]        *x31    *x52
        +coeff[ 78]            *x43*x51
        +coeff[ 79]*x11*x21*x32        
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 80]*x12    *x31*x41    
        +coeff[ 81]    *x23    *x42    
        +coeff[ 82]*x11*x21*x31    *x51
        +coeff[ 83]    *x22*x32*x41    
        +coeff[ 84]*x13*x21            
        +coeff[ 85]    *x22*x31*x41*x51
        +coeff[ 86]    *x22    *x44    
        +coeff[ 87]*x11*x23*x31*x41    
        +coeff[ 88]    *x22*x31*x44    
    ;
    v_y_r5p65_ep23_den                        =v_y_r5p65_ep23_den                        
        +coeff[ 89]    *x22*x32*x43    
        +coeff[ 90]    *x24*x31*x42    
        +coeff[ 91]    *x22*x33*x42    
        +coeff[ 92]    *x23*x33*x41    
        +coeff[ 93]*x11*x24    *x42    
        +coeff[ 94]    *x24    *x44    
        +coeff[ 95]*x12*x21*x31*x42*x51
        +coeff[ 96]        *x31*x43    
        ;

    return v_y_r5p65_ep23_den                        ;
}
float p_r5p65_ep23_den                        (float *x,int m){
    int ncoeff= 60;
    float avdat=  0.8674432E-02;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 61]={
         0.15554455E-02, 0.35743448E-02,-0.13924614E-01,-0.94233006E-01,
        -0.91408826E-02, 0.13428261E-01,-0.16539840E-01,-0.31519688E-02,
         0.15188347E-01, 0.35147325E-02,-0.24219034E-02,-0.22457782E-02,
         0.16068616E-02, 0.23381256E-02, 0.19664895E-02,-0.41102041E-02,
        -0.17617004E-02, 0.26005262E-02, 0.52351770E-02, 0.30340634E-02,
         0.74088937E-02, 0.11059816E-03, 0.20058453E-02, 0.11964390E-02,
         0.18202128E-02,-0.17774763E-02, 0.19335238E-02, 0.69674285E-03,
        -0.15928471E-02,-0.10333421E-02, 0.29004402E-02,-0.84178027E-03,
        -0.15876084E-03, 0.29899113E-03,-0.60557964E-03, 0.22062361E-03,
        -0.11734110E-02, 0.38343098E-03, 0.67363831E-03, 0.23652106E-02,
        -0.39293958E-03,-0.11705741E-02,-0.96851290E-03,-0.25844099E-02,
         0.13880514E-02, 0.22905113E-03, 0.73705116E-04, 0.17234625E-03,
         0.26101150E-03, 0.28141367E-03,-0.36220692E-03,-0.21496847E-03,
         0.12485348E-02,-0.41504548E-03, 0.40868705E-03,-0.59308234E-03,
         0.52969455E-03,-0.17662843E-02,-0.49895805E-03,-0.97924320E-03,
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

    float v_p_r5p65_ep23_den                        =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21    *x41    
        +coeff[  7]            *x42    
    ;
    v_p_r5p65_ep23_den                        =v_p_r5p65_ep23_den                        
        +coeff[  8]            *x41*x51
        +coeff[  9]*x11*x21            
        +coeff[ 10]    *x21*x31        
        +coeff[ 11]    *x21        *x51
        +coeff[ 12]        *x31    *x51
        +coeff[ 13]                *x52
        +coeff[ 14]*x11        *x41    
        +coeff[ 15]    *x22        *x51
        +coeff[ 16]        *x31*x41    
    ;
    v_p_r5p65_ep23_den                        =v_p_r5p65_ep23_den                        
        +coeff[ 17]    *x22    *x41    
        +coeff[ 18]    *x21    *x42    
        +coeff[ 19]*x11*x21    *x41    
        +coeff[ 20]    *x22    *x42    
        +coeff[ 21]*x11                
        +coeff[ 22]    *x23            
        +coeff[ 23]    *x22*x31        
        +coeff[ 24]    *x21*x31*x41    
        +coeff[ 25]            *x43    
    ;
    v_p_r5p65_ep23_den                        =v_p_r5p65_ep23_den                        
        +coeff[ 26]    *x21    *x41*x51
        +coeff[ 27]*x11*x21*x31        
        +coeff[ 28]    *x23    *x41    
        +coeff[ 29]            *x41*x52
        +coeff[ 30]    *x22*x31*x41    
        +coeff[ 31]*x11*x21        *x51
        +coeff[ 32]        *x32        
        +coeff[ 33]*x11    *x31        
        +coeff[ 34]        *x31*x42    
    ;
    v_p_r5p65_ep23_den                        =v_p_r5p65_ep23_den                        
        +coeff[ 35]*x12                
        +coeff[ 36]    *x24            
        +coeff[ 37]    *x21        *x52
        +coeff[ 38]*x11        *x42    
        +coeff[ 39]    *x21    *x43    
        +coeff[ 40]*x11        *x41*x51
        +coeff[ 41]    *x22    *x41*x51
        +coeff[ 42]*x11*x23            
        +coeff[ 43]    *x24    *x41    
    ;
    v_p_r5p65_ep23_den                        =v_p_r5p65_ep23_den                        
        +coeff[ 44]*x11*x21    *x42    
        +coeff[ 45]    *x21*x32        
        +coeff[ 46]*x11            *x51
        +coeff[ 47]    *x21*x31    *x51
        +coeff[ 48]        *x31*x41*x51
        +coeff[ 49]*x11    *x31*x41    
        +coeff[ 50]    *x23        *x51
        +coeff[ 51]                *x53
        +coeff[ 52]    *x21*x31*x42    
    ;
    v_p_r5p65_ep23_den                        =v_p_r5p65_ep23_den                        
        +coeff[ 53]    *x22*x31    *x51
        +coeff[ 54]    *x22        *x52
        +coeff[ 55]*x11*x22    *x41    
        +coeff[ 56]*x11*x21*x31*x41    
        +coeff[ 57]    *x23    *x42    
        +coeff[ 58]*x11*x21    *x41*x51
        +coeff[ 59]*x11*x23    *x41    
        ;

    return v_p_r5p65_ep23_den                        ;
}
float l_r6_h90s_ep23_den                      (float *x,int m){
    int ncoeff= 90;
    float avdat= -0.3138593E-02;
    float xmin[10]={
        -0.99915E-02,-0.49497E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 91]={
         0.35318267E-02, 0.24584687E-03,-0.11207685E-02,-0.39965520E-02,
         0.59733273E-04,-0.27798354E-04,-0.75788437E-02,-0.11386698E-02,
        -0.67949095E-02,-0.17992558E-02, 0.10281451E-02, 0.59943777E-02,
         0.48913667E-03, 0.83720719E-03,-0.10782413E-02,-0.36151364E-03,
         0.23552144E-02,-0.44209295E-03, 0.51776500E-03,-0.15363198E-03,
        -0.13601416E-03,-0.24575274E-03, 0.10918644E-02,-0.78168634E-03,
         0.21766585E-03, 0.94879814E-03,-0.43780578E-03,-0.21371944E-02,
        -0.11049508E-04, 0.50356775E-04,-0.70029237E-04, 0.21819961E-03,
        -0.45899538E-03,-0.12035358E-03, 0.10112408E-03, 0.86937369E-04,
         0.57510147E-06,-0.70036290E-03,-0.10287237E-03,-0.60178365E-04,
         0.67556932E-04,-0.81637168E-04,-0.54723121E-04, 0.10650127E-03,
         0.38758319E-03, 0.10598761E-02, 0.48587854E-05, 0.75392178E-04,
         0.10280613E-04, 0.21895248E-04, 0.14939616E-04,-0.78239180E-04,
         0.90594251E-04,-0.16512968E-03, 0.55787958E-04, 0.12533882E-03,
         0.72791831E-05,-0.48407299E-04, 0.62793063E-03,-0.14511628E-03,
        -0.50969469E-04,-0.10433007E-03,-0.40079904E-05, 0.60311627E-05,
         0.22787657E-04, 0.10349250E-03,-0.12285663E-03, 0.47774385E-04,
        -0.54480952E-04,-0.22315258E-04, 0.95016447E-04,-0.85853615E-04,
         0.21432282E-03,-0.35707555E-04,-0.13090132E-04, 0.38214464E-04,
         0.33770641E-04, 0.13001743E-03,-0.82265105E-04, 0.11712056E-04,
        -0.97376424E-05,-0.89039631E-05, 0.96516087E-05,-0.11353890E-04,
         0.51873023E-04,-0.23632443E-04,-0.21474991E-04,-0.21129856E-04,
        -0.14832413E-04, 0.11278070E-04,
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

    float v_l_r6_h90s_ep23_den                      =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]*x11                
        +coeff[  6]    *x22            
        +coeff[  7]        *x31*x41    
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[  8]            *x42    
        +coeff[  9]            *x41*x51
        +coeff[ 10]*x11*x21            
        +coeff[ 11]    *x22    *x41    
        +coeff[ 12]            *x42*x51
        +coeff[ 13]*x11*x21    *x41    
        +coeff[ 14]    *x24            
        +coeff[ 15]    *x23    *x41    
        +coeff[ 16]    *x22    *x42    
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 17]    *x24*x31        
        +coeff[ 18]    *x21    *x41    
        +coeff[ 19]        *x31    *x51
        +coeff[ 20]                *x52
        +coeff[ 21]    *x23            
        +coeff[ 22]    *x22*x31        
        +coeff[ 23]            *x43    
        +coeff[ 24]            *x41*x52
        +coeff[ 25]    *x22*x31*x41    
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 26]*x11*x23            
        +coeff[ 27]    *x24    *x41    
        +coeff[ 28]        *x33*x42    
        +coeff[ 29]    *x21        *x51
        +coeff[ 30]*x12                
        +coeff[ 31]    *x21    *x42    
        +coeff[ 32]        *x31*x42    
        +coeff[ 33]    *x22        *x51
        +coeff[ 34]        *x31*x41*x51
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 35]*x11*x21        *x51
        +coeff[ 36]    *x21*x33        
        +coeff[ 37]*x11*x23    *x41    
        +coeff[ 38]*x11*x23    *x42    
        +coeff[ 39]        *x32        
        +coeff[ 40]    *x21*x31*x41    
        +coeff[ 41]        *x32*x41    
        +coeff[ 42]*x11*x22            
        +coeff[ 43]*x11*x21*x31        
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 44]*x11*x21    *x42    
        +coeff[ 45]    *x22    *x43    
        +coeff[ 46]    *x24*x31*x42    
        +coeff[ 47]    *x21*x31        
        +coeff[ 48]*x11            *x51
        +coeff[ 49]        *x31    *x52
        +coeff[ 50]*x12            *x51
        +coeff[ 51]    *x23*x31        
        +coeff[ 52]    *x22*x32        
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 53]            *x44    
        +coeff[ 54]            *x43*x51
        +coeff[ 55]*x11*x21*x31*x41    
        +coeff[ 56]*x11        *x43    
        +coeff[ 57]*x12*x22            
        +coeff[ 58]    *x22*x31*x42    
        +coeff[ 59]*x11*x23*x31        
        +coeff[ 60]    *x22*x31*x43    
        +coeff[ 61]*x11*x24    *x41    
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 62]*x11    *x31        
        +coeff[ 63]    *x21*x31    *x51
        +coeff[ 64]                *x53
        +coeff[ 65]    *x21    *x43    
        +coeff[ 66]        *x31*x43    
        +coeff[ 67]        *x31*x42*x51
        +coeff[ 68]            *x42*x52
        +coeff[ 69]*x11*x21        *x52
        +coeff[ 70]    *x22*x32*x41    
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 71]    *x23    *x42    
        +coeff[ 72]*x11*x21    *x43    
        +coeff[ 73]        *x34*x42    
        +coeff[ 74]    *x21*x31*x44    
        +coeff[ 75]    *x24*x31    *x51
        +coeff[ 76]*x11*x23*x31    *x51
        +coeff[ 77]*x11*x21*x31*x44    
        +coeff[ 78]*x12*x24    *x41    
        +coeff[ 79]*x11        *x41    
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 80]    *x21    *x41*x51
        +coeff[ 81]    *x21        *x52
        +coeff[ 82]*x11        *x42    
        +coeff[ 83]*x12    *x31        
        +coeff[ 84]    *x21*x31*x42    
        +coeff[ 85]    *x22    *x41*x51
        +coeff[ 86]    *x22        *x52
        +coeff[ 87]            *x41*x53
        +coeff[ 88]*x11*x22*x31        
    ;
    v_l_r6_h90s_ep23_den                      =v_l_r6_h90s_ep23_den                      
        +coeff[ 89]*x11*x21*x32        
        ;

    return v_l_r6_h90s_ep23_den                      ;
}
float x_r5p65_ep25_dex                        (float *x,int m){
    int ncoeff= 55;
    float avdat= -0.2566867E-02;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 56]={
        -0.29103137E-02,-0.30694257E-01, 0.13283221E+00, 0.30549100E+00,
         0.37515704E-01, 0.15738467E-01,-0.90945736E-01, 0.26174223E-01,
         0.52419410E-03, 0.10118104E-01,-0.60176849E-02,-0.13228908E-01,
        -0.23619832E-01, 0.11487924E-01,-0.41116841E-01, 0.46192303E-01,
         0.22195096E-02,-0.18375329E-02,-0.34225758E-02,-0.34913600E-02,
        -0.59058601E-02,-0.17723426E-01, 0.17739287E-01,-0.50254771E-02,
         0.61889812E-02, 0.86802518E-03,-0.14009245E-02,-0.13581016E-02,
         0.18444364E-02, 0.24705913E-02,-0.37975232E-02,-0.20421634E-02,
         0.42411992E-02,-0.24180561E-03, 0.69271430E-03,-0.87561738E-03,
        -0.17871023E-02,-0.71865949E-03, 0.12496805E-02, 0.20975978E-02,
         0.21932914E-02,-0.73747598E-02,-0.93253488E-02, 0.29467302E-03,
         0.19546573E-03, 0.49199187E-03,-0.27000092E-03, 0.27490340E-03,
        -0.81037363E-03,-0.81273139E-03, 0.14550254E-02,-0.20590117E-02,
         0.23191327E-02, 0.53379568E-03, 0.30656236E-02,
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

    float v_x_r5p65_ep25_dex                        =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]                *x51
        +coeff[  3]    *x21            
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21    *x41    
        +coeff[  7]    *x23            
    ;
    v_x_r5p65_ep25_dex                        =v_x_r5p65_ep25_dex                        
        +coeff[  8]*x12*x21*x31        
        +coeff[  9]    *x23*x31        
        +coeff[ 10]            *x41    
        +coeff[ 11]*x11        *x41    
        +coeff[ 12]    *x21*x31        
        +coeff[ 13]*x11*x22            
        +coeff[ 14]    *x21    *x42    
        +coeff[ 15]    *x23    *x41    
        +coeff[ 16]    *x23*x31*x41    
    ;
    v_x_r5p65_ep25_dex                        =v_x_r5p65_ep25_dex                        
        +coeff[ 17]        *x31        
        +coeff[ 18]                *x52
        +coeff[ 19]*x11    *x31        
        +coeff[ 20]            *x42    
        +coeff[ 21]    *x21*x31*x41    
        +coeff[ 22]*x11*x22    *x41    
        +coeff[ 23]    *x21    *x41*x51
        +coeff[ 24]    *x22    *x41    
        +coeff[ 25]*x11*x21            
    ;
    v_x_r5p65_ep25_dex                        =v_x_r5p65_ep25_dex                        
        +coeff[ 26]            *x41*x51
        +coeff[ 27]        *x31*x41    
        +coeff[ 28]*x12*x21            
        +coeff[ 29]*x11*x21    *x41    
        +coeff[ 30]*x11        *x42    
        +coeff[ 31]    *x21*x32        
        +coeff[ 32]*x11*x22*x31        
        +coeff[ 33]*x13    *x31*x41    
        +coeff[ 34]*x11            *x51
    ;
    v_x_r5p65_ep25_dex                        =v_x_r5p65_ep25_dex                        
        +coeff[ 35]    *x21        *x52
        +coeff[ 36]*x11    *x31*x41    
        +coeff[ 37]    *x21*x31    *x51
        +coeff[ 38]    *x22*x31        
        +coeff[ 39]*x12*x21    *x41    
        +coeff[ 40]    *x23        *x51
        +coeff[ 41]    *x21*x31*x42    
        +coeff[ 42]    *x21    *x43    
        +coeff[ 43]    *x23*x32*x41    
    ;
    v_x_r5p65_ep25_dex                        =v_x_r5p65_ep25_dex                        
        +coeff[ 44]*x12                
        +coeff[ 45]*x11*x21*x31        
        +coeff[ 46]*x11    *x32        
        +coeff[ 47]    *x22        *x51
        +coeff[ 48]            *x42*x51
        +coeff[ 49]            *x43    
        +coeff[ 50]    *x22    *x41*x51
        +coeff[ 51]    *x21*x32*x41    
        +coeff[ 52]    *x22    *x42    
    ;
    v_x_r5p65_ep25_dex                        =v_x_r5p65_ep25_dex                        
        +coeff[ 53]*x13*x21        *x51
        +coeff[ 54]    *x23    *x42    
        ;

    return v_x_r5p65_ep25_dex                        ;
}
float t_r5p65_ep25_dex                        (float *x,int m){
    int ncoeff= 80;
    float avdat=  0.5942838E+00;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 81]={
        -0.19943186E-03,-0.94956353E-01, 0.70004053E-02, 0.18818976E-02,
         0.30770853E-01, 0.25072621E-01,-0.66667576E-02,-0.73670633E-02,
        -0.26552579E-02, 0.65055750E-02, 0.35684884E-02,-0.32470161E-02,
         0.11213773E-01,-0.12469746E-01,-0.18769713E-02,-0.12980276E-02,
         0.46253074E-02,-0.48895217E-02, 0.47271373E-03, 0.93954196E-03,
         0.48848940E-03,-0.41380539E-03,-0.46671354E-03,-0.61491592E-03,
        -0.12102973E-02,-0.29639274E-03,-0.11250462E-02, 0.73100731E-03,
        -0.51156612E-03, 0.52140409E-03, 0.46942596E-03, 0.60408388E-03,
         0.11697524E-02,-0.52159361E-03,-0.11424850E-02,-0.57078942E-04,
         0.17822165E-04,-0.56660484E-03,-0.46262675E-03, 0.59224531E-03,
         0.87558076E-03, 0.27210667E-03,-0.61664946E-03, 0.22096341E-02,
         0.31120183E-02, 0.58596306E-04,-0.59677939E-04, 0.61067694E-04,
        -0.30899220E-03,-0.15568639E-03, 0.57468150E-03, 0.24872168E-03,
         0.50968450E-03,-0.51309547E-03,-0.47592839E-03, 0.84242946E-03,
         0.30530689E-03,-0.23598644E-03,-0.95267300E-04, 0.76619006E-04,
         0.11596172E-03,-0.39645022E-03, 0.41553951E-03,-0.14968478E-03,
         0.21308385E-03, 0.31722867E-03,-0.55748515E-03, 0.25799155E-03,
        -0.89943170E-03,-0.96333167E-03,-0.26523587E-03,-0.63187850E-03,
         0.34335564E-03,-0.14209244E-03,-0.92600042E-03,-0.44506247E-03,
         0.20939804E-03,-0.18669949E-03, 0.83229519E-04,-0.23611345E-04,
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
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_r5p65_ep25_dex                        =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]*x11                
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x21    *x41    
        +coeff[  6]    *x21        *x51
        +coeff[  7]    *x23            
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[  8]    *x23*x31        
        +coeff[  9]    *x21*x31        
        +coeff[ 10]*x11        *x41    
        +coeff[ 11]*x11*x22            
        +coeff[ 12]    *x21    *x42    
        +coeff[ 13]    *x23    *x41    
        +coeff[ 14]    *x22            
        +coeff[ 15]                *x52
        +coeff[ 16]    *x21*x31*x41    
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[ 17]*x11*x22    *x41    
        +coeff[ 18]        *x31        
        +coeff[ 19]*x11    *x31        
        +coeff[ 20]        *x31*x41    
        +coeff[ 21]            *x42    
        +coeff[ 22]*x11            *x51
        +coeff[ 23]            *x41*x51
        +coeff[ 24]            *x42*x51
        +coeff[ 25]*x11*x21            
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[ 26]    *x22    *x41    
        +coeff[ 27]    *x22        *x51
        +coeff[ 28]*x12*x21            
        +coeff[ 29]    *x21*x32        
        +coeff[ 30]    *x21    *x41*x51
        +coeff[ 31]    *x21        *x52
        +coeff[ 32]*x11        *x42    
        +coeff[ 33]            *x41*x52
        +coeff[ 34]*x11*x22*x31        
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[ 35]    *x24*x31        
        +coeff[ 36]*x13    *x31*x41    
        +coeff[ 37]    *x22*x31        
        +coeff[ 38]*x11*x21    *x41    
        +coeff[ 39]*x11    *x31*x41    
        +coeff[ 40]    *x22    *x41*x51
        +coeff[ 41]    *x22        *x52
        +coeff[ 42]*x12*x21    *x41    
        +coeff[ 43]    *x21*x31*x42    
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[ 44]    *x21    *x43    
        +coeff[ 45]*x11*x23*x31        
        +coeff[ 46]*x12                
        +coeff[ 47]        *x31    *x51
        +coeff[ 48]    *x23        *x51
        +coeff[ 49]*x11*x22        *x51
        +coeff[ 50]    *x21*x32*x41    
        +coeff[ 51]*x11*x21    *x41*x51
        +coeff[ 52]    *x21    *x42*x51
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[ 53]    *x24        *x51
        +coeff[ 54]            *x42*x52
        +coeff[ 55]    *x22    *x42*x51
        +coeff[ 56]    *x24*x31*x41    
        +coeff[ 57]*x11*x21*x31        
        +coeff[ 58]    *x21*x31    *x51
        +coeff[ 59]*x11    *x32        
        +coeff[ 60]        *x31*x41*x51
        +coeff[ 61]    *x22*x31*x41    
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[ 62]    *x22    *x42    
        +coeff[ 63]*x12*x21*x31        
        +coeff[ 64]    *x21    *x41*x52
        +coeff[ 65]*x11*x24            
        +coeff[ 66]    *x24    *x41    
        +coeff[ 67]*x11        *x43    
        +coeff[ 68]    *x23*x31*x41    
        +coeff[ 69]    *x23    *x42    
        +coeff[ 70]    *x23        *x52
    ;
    v_t_r5p65_ep25_dex                        =v_t_r5p65_ep25_dex                        
        +coeff[ 71]*x11*x22    *x42    
        +coeff[ 72]    *x22    *x41*x52
        +coeff[ 73]*x12        *x43    
        +coeff[ 74]    *x23    *x43    
        +coeff[ 75]*x11*x24*x31*x41    
        +coeff[ 76]*x13    *x31*x42    
        +coeff[ 77]            *x43*x53
        +coeff[ 78]    *x24            
        +coeff[ 79]*x12    *x31        
        ;

    return v_t_r5p65_ep25_dex                        ;
}
float y_r5p65_ep25_dex                        (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.4209658E-02;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
         0.10077215E-02, 0.68609692E-01,-0.11250547E-01, 0.21313868E-01,
        -0.52791037E-01,-0.35505150E-01,-0.46432256E-02, 0.33189717E-02,
         0.51286746E-01,-0.10262300E-01,-0.79622669E-02, 0.49575847E-02,
        -0.34150274E-02, 0.80788806E-02, 0.58398773E-02,-0.44455145E-01,
         0.28245660E-03, 0.12569465E-01,-0.19794848E-04,-0.11517981E-01,
        -0.36485938E-03, 0.55801950E-03, 0.13297167E-01, 0.62344274E-02,
        -0.36380803E-02,-0.35498370E-02,-0.39302288E-02, 0.12503776E-01,
         0.45435848E-02,-0.82136499E-03, 0.18757740E-02,-0.74706078E-02,
         0.29702694E-02,-0.13224955E-02, 0.63882354E-02,-0.18671316E-02,
         0.16946731E-04,-0.40107081E-02,-0.11398557E-02, 0.71904389E-03,
        -0.15684321E-03, 0.40604963E-03, 0.19105829E-02,-0.35555423E-02,
        -0.30899737E-02,-0.75220258E-03,-0.19996988E-02, 0.64610707E-04,
        -0.41913931E-03, 0.16754775E-03, 0.81835635E-03,-0.66632708E-03,
        -0.54429949E-03, 0.25199640E-02,-0.44991649E-03, 0.11380129E-02,
        -0.71649527E-03,-0.16505445E-02, 0.15707735E-02,-0.10636240E-02,
        -0.13718402E-02,-0.53510624E-02, 0.83058770E-03, 0.82752388E-03,
        -0.27622730E-02,-0.29232732E-02,-0.45562661E-02,-0.20548794E-02,
         0.26338937E-03, 0.61322912E-03, 0.63487684E-03,-0.14705838E-03,
         0.59455284E-03, 0.45273383E-03, 0.36694398E-02, 0.36682212E-03,
         0.75481221E-03,-0.21058301E-03,-0.29097070E-03, 0.36416706E-03,
        -0.29935325E-02, 0.95917337E-03, 0.61375013E-03,-0.40923603E-03,
        -0.14115724E-02,-0.83258573E-03,-0.12808206E-02,-0.21185793E-03,
         0.15752906E-02, 0.81520126E-03, 0.18241307E-03, 0.65426451E-04,
        -0.11785243E-03, 0.10134961E-03,-0.77521378E-04,-0.20637411E-03,
         0.32529497E-03,
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

    float v_y_r5p65_ep25_dex                        =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]        *x31        
        +coeff[  3]                *x51
        +coeff[  4]    *x21    *x41    
        +coeff[  5]    *x22            
        +coeff[  6]    *x21*x31        
        +coeff[  7]*x11        *x41    
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[  8]            *x41*x51
        +coeff[  9]*x11*x21            
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]            *x43    
        +coeff[ 13]    *x21    *x42    
        +coeff[ 14]                *x52
        +coeff[ 15]    *x22    *x41    
        +coeff[ 16]            *x42*x51
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 17]    *x23            
        +coeff[ 18]            *x44    
        +coeff[ 19]    *x22        *x51
        +coeff[ 20]        *x31*x43    
        +coeff[ 21]    *x24*x31        
        +coeff[ 22]            *x42    
        +coeff[ 23]        *x31*x41    
        +coeff[ 24]*x11*x21    *x41    
        +coeff[ 25]    *x21    *x41*x51
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 26]            *x41*x52
        +coeff[ 27]    *x24            
        +coeff[ 28]*x11*x23            
        +coeff[ 29]*x12                
        +coeff[ 30]    *x21*x31*x41    
        +coeff[ 31]    *x22*x31        
        +coeff[ 32]*x11*x22            
        +coeff[ 33]*x11*x21*x31        
        +coeff[ 34]    *x21    *x43    
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 35]*x11*x21        *x51
        +coeff[ 36]*x13                
        +coeff[ 37]    *x22    *x41*x51
        +coeff[ 38]*x11                
        +coeff[ 39]        *x32        
        +coeff[ 40]        *x31*x42    
        +coeff[ 41]*x11            *x51
        +coeff[ 42]*x11        *x42    
        +coeff[ 43]    *x22*x31*x41    
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 44]*x11*x21    *x42    
        +coeff[ 45]                *x53
        +coeff[ 46]*x11*x22    *x41    
        +coeff[ 47]*x11    *x33*x41    
        +coeff[ 48]    *x21            
        +coeff[ 49]*x11    *x31        
        +coeff[ 50]        *x31*x41*x51
        +coeff[ 51]*x12        *x41    
        +coeff[ 52]    *x22    *x42    
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 53]    *x21*x31*x42    
        +coeff[ 54]        *x31    *x52
        +coeff[ 55]    *x23*x31        
        +coeff[ 56]    *x22*x32        
        +coeff[ 57]*x11*x21*x31*x41    
        +coeff[ 58]    *x21    *x44    
        +coeff[ 59]    *x22*x31    *x51
        +coeff[ 60]*x11*x21    *x41*x51
        +coeff[ 61]    *x23    *x42    
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 62]*x12*x22            
        +coeff[ 63]    *x22        *x52
        +coeff[ 64]    *x24    *x41    
        +coeff[ 65]*x11*x22    *x42    
        +coeff[ 66]    *x22    *x42*x51
        +coeff[ 67]    *x23    *x43    
        +coeff[ 68]    *x22    *x45    
        +coeff[ 69]            *x43*x53
        +coeff[ 70]*x11    *x31*x41    
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 71]*x12    *x31        
        +coeff[ 72]*x11        *x43    
        +coeff[ 73]            *x43*x51
        +coeff[ 74]    *x23    *x41    
        +coeff[ 75]*x11    *x31*x42    
        +coeff[ 76]    *x21    *x42*x51
        +coeff[ 77]*x12            *x51
        +coeff[ 78]*x11*x21*x32        
        +coeff[ 79]    *x23        *x51
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 80]    *x22    *x43    
        +coeff[ 81]    *x21*x31*x43    
        +coeff[ 82]    *x21    *x41*x52
        +coeff[ 83]*x11*x21*x31    *x51
        +coeff[ 84]    *x23*x31*x41    
        +coeff[ 85]*x11*x22*x31*x41    
        +coeff[ 86]    *x22*x31*x41*x51
        +coeff[ 87]    *x21    *x42*x52
        +coeff[ 88]*x11*x24    *x41    
    ;
    v_y_r5p65_ep25_dex                        =v_y_r5p65_ep25_dex                        
        +coeff[ 89]*x11*x21    *x42*x52
        +coeff[ 90]        *x32*x41    
        +coeff[ 91]*x11    *x32        
        +coeff[ 92]    *x21*x31    *x51
        +coeff[ 93]        *x32    *x51
        +coeff[ 94]*x11        *x41*x51
        +coeff[ 95]    *x21        *x52
        +coeff[ 96]    *x21*x32*x41    
        ;

    return v_y_r5p65_ep25_dex                        ;
}
float p_r5p65_ep25_dex                        (float *x,int m){
    int ncoeff= 90;
    float avdat=  0.7786233E-03;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 91]={
         0.28769337E-03, 0.37090640E-03,-0.41250638E-02,-0.31384684E-02,
         0.22104739E-02,-0.37658652E-02,-0.81759243E-03,-0.11725784E-01,
         0.17808421E-02,-0.18900529E-02, 0.94354030E-03, 0.11257352E-01,
        -0.12044446E-02, 0.32464713E-02, 0.15105520E-02,-0.86301588E-03,
         0.81020140E-03,-0.85543292E-02, 0.22904049E-02,-0.32138284E-02,
        -0.19765934E-02, 0.22006347E-02,-0.19288680E-03, 0.18614321E-04,
         0.74079086E-03, 0.47764095E-03,-0.85756840E-03, 0.61283872E-03,
        -0.26372340E-03,-0.46860936E-03, 0.13118871E-02,-0.52902440E-03,
        -0.15023084E-02, 0.63055800E-03, 0.10079316E-03,-0.15980627E-03,
        -0.10018328E-03,-0.14128130E-03, 0.89681655E-03, 0.47915019E-03,
        -0.14072737E-03, 0.14699594E-02,-0.21959520E-02,-0.15157989E-03,
         0.66777196E-04, 0.22554911E-03, 0.24235668E-03, 0.15398644E-04,
         0.44646207E-03, 0.65746915E-03,-0.29971887E-03,-0.43442487E-05,
         0.30257666E-03,-0.10756797E-03,-0.70576416E-03, 0.43594762E-03,
        -0.34147405E-03, 0.38575559E-03, 0.21878853E-03,-0.24558906E-03,
        -0.68116514E-03, 0.17663611E-04, 0.16233481E-03, 0.78726101E-04,
         0.42925571E-04,-0.60410111E-03, 0.19713741E-03,-0.22094174E-03,
        -0.58634138E-04,-0.10921493E-03, 0.38111000E-03,-0.74091711E-03,
         0.11262017E-03,-0.20646369E-03, 0.97589439E-03,-0.60596358E-03,
        -0.54971135E-03, 0.18785492E-03, 0.52374153E-03,-0.26395914E-03,
         0.25412341E-03,-0.23489515E-03, 0.10701607E-03,-0.64162123E-04,
        -0.27828396E-03, 0.83704675E-04, 0.22271057E-04,-0.62407795E-04,
         0.50117596E-04, 0.49980776E-04,
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
    float x53 = x52*x5;

//                 function

    float v_p_r5p65_ep25_dex                        =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21*x31        
        +coeff[  7]    *x21    *x41    
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[  8]            *x42    
        +coeff[  9]    *x21        *x51
        +coeff[ 10]        *x31    *x51
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21            
        +coeff[ 13]    *x23            
        +coeff[ 14]                *x52
        +coeff[ 15]    *x22*x31        
        +coeff[ 16]*x11        *x41    
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 17]    *x22    *x41    
        +coeff[ 18]    *x21    *x42    
        +coeff[ 19]    *x22        *x51
        +coeff[ 20]    *x21    *x41*x51
        +coeff[ 21]    *x24            
        +coeff[ 22]    *x22*x31*x41    
        +coeff[ 23]        *x33*x41    
        +coeff[ 24]        *x31*x41    
        +coeff[ 25]    *x21*x31*x41    
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 26]            *x43    
        +coeff[ 27]*x11*x22            
        +coeff[ 28]    *x21        *x52
        +coeff[ 29]            *x41*x52
        +coeff[ 30]    *x22    *x42    
        +coeff[ 31]*x11*x21        *x51
        +coeff[ 32]    *x22    *x41*x51
        +coeff[ 33]*x11*x23            
        +coeff[ 34]*x11            *x51
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 35]        *x31*x42    
        +coeff[ 36]*x12                
        +coeff[ 37]*x11*x21*x31        
        +coeff[ 38]    *x23    *x41    
        +coeff[ 39]    *x23        *x51
        +coeff[ 40]                *x53
        +coeff[ 41]    *x21    *x43    
        +coeff[ 42]    *x23    *x42    
        +coeff[ 43]*x11                
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 44]*x11    *x31        
        +coeff[ 45]        *x31*x41*x51
        +coeff[ 46]            *x42*x51
        +coeff[ 47]*x11*x21    *x41    
        +coeff[ 48]*x11        *x42    
        +coeff[ 49]    *x21*x31*x42    
        +coeff[ 50]    *x22*x31    *x51
        +coeff[ 51]        *x34        
        +coeff[ 52]    *x21    *x42*x51
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 53]*x12        *x41    
        +coeff[ 54]*x11*x22    *x41    
        +coeff[ 55]    *x24        *x51
        +coeff[ 56]*x11*x21    *x41*x51
        +coeff[ 57]    *x25*x31        
        +coeff[ 58]    *x22    *x41*x52
        +coeff[ 59]*x11*x22*x31*x41    
        +coeff[ 60]    *x25*x31*x41    
        +coeff[ 61]        *x32    *x51
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 62]*x11    *x31*x41    
        +coeff[ 63]    *x21*x32*x41    
        +coeff[ 64]*x11        *x41*x51
        +coeff[ 65]    *x25            
        +coeff[ 66]            *x43*x51
        +coeff[ 67]*x11*x21    *x42    
        +coeff[ 68]*x12            *x51
        +coeff[ 69]*x11*x21*x31    *x51
        +coeff[ 70]    *x23    *x41*x51
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 71]    *x22    *x42*x51
        +coeff[ 72]*x12*x22            
        +coeff[ 73]*x11*x23    *x41    
        +coeff[ 74]    *x25    *x41    
        +coeff[ 75]*x11*x22    *x42    
        +coeff[ 76]    *x24    *x42    
        +coeff[ 77]*x11*x23        *x51
        +coeff[ 78]*x11*x24    *x41    
        +coeff[ 79]*x11*x23*x31*x41    
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 80]    *x24*x31*x42    
        +coeff[ 81]    *x21    *x43*x52
        +coeff[ 82]*x11        *x45    
        +coeff[ 83]*x12    *x31    *x52
        +coeff[ 84]    *x22*x33*x41*x51
        +coeff[ 85]        *x32        
        +coeff[ 86]    *x21*x31    *x51
        +coeff[ 87]    *x22*x32        
        +coeff[ 88]    *x21*x31*x41*x51
    ;
    v_p_r5p65_ep25_dex                        =v_p_r5p65_ep25_dex                        
        +coeff[ 89]    *x22        *x52
        ;

    return v_p_r5p65_ep25_dex                        ;
}
float l_r6_h90s_ep25_dex                      (float *x,int m){
    int ncoeff= 52;
    float avdat= -0.9308472E-02;
    float xmin[10]={
        -0.99915E-02,-0.49354E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 53]={
         0.88067399E-02,-0.43595850E+00, 0.48378124E-02,-0.98304115E-01,
         0.39512370E-01,-0.35672378E-01, 0.12610146E+00,-0.36474947E-01,
         0.18148681E-01,-0.37754003E-01,-0.13637834E-01, 0.32192204E-01,
         0.50860047E-01,-0.15567970E-01,-0.62779114E-01,-0.37580207E-02,
         0.47476934E-02, 0.24266064E-01,-0.24630014E-01, 0.13258898E-02,
         0.16745517E-02,-0.21371255E-02, 0.55904896E-02, 0.44525546E-03,
        -0.12264986E-02, 0.28099997E-02, 0.25762739E-02, 0.14107167E-02,
         0.27911428E-02,-0.25474126E-02, 0.12986844E-01,-0.57280078E-02,
        -0.18398060E-02,-0.16752514E-02, 0.72684913E-03,-0.11098804E-02,
         0.13433100E-02,-0.24199553E-02, 0.11642758E-01,-0.28395471E-02,
        -0.13817615E-02,-0.37381932E-03, 0.43904880E-03,-0.56383479E-03,
         0.43759480E-03,-0.18904897E-02, 0.33096322E-02,-0.16569305E-02,
         0.28990118E-02, 0.56273025E-02, 0.38262401E-02,-0.65088905E-02,
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
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_l_r6_h90s_ep25_dex                      =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]            *x41    
        +coeff[  3]                *x51
        +coeff[  4]*x11                
        +coeff[  5]    *x22            
        +coeff[  6]    *x21    *x41    
        +coeff[  7]    *x21        *x51
    ;
    v_l_r6_h90s_ep25_dex                      =v_l_r6_h90s_ep25_dex                      
        +coeff[  8]*x11        *x41    
        +coeff[  9]    *x23            
        +coeff[ 10]    *x23*x31        
        +coeff[ 11]    *x21*x31        
        +coeff[ 12]    *x21    *x42    
        +coeff[ 13]*x11*x22            
        +coeff[ 14]    *x23    *x41    
        +coeff[ 15]    *x23*x31*x41    
        +coeff[ 16]*x11    *x31        
    ;
    v_l_r6_h90s_ep25_dex                      =v_l_r6_h90s_ep25_dex                      
        +coeff[ 17]    *x21*x31*x41    
        +coeff[ 18]*x11*x22    *x41    
        +coeff[ 19]        *x31        
        +coeff[ 20]                *x52
        +coeff[ 21]*x11            *x51
        +coeff[ 22]*x11        *x42    
        +coeff[ 23]        *x31*x41    
        +coeff[ 24]            *x41*x51
        +coeff[ 25]    *x21*x32        
    ;
    v_l_r6_h90s_ep25_dex                      =v_l_r6_h90s_ep25_dex                      
        +coeff[ 26]    *x21    *x41*x51
        +coeff[ 27]    *x21        *x52
        +coeff[ 28]*x11    *x31*x41    
        +coeff[ 29]*x12*x21            
        +coeff[ 30]    *x21    *x43    
        +coeff[ 31]*x11*x22*x31        
        +coeff[ 32]    *x23*x31*x42    
        +coeff[ 33]            *x42    
        +coeff[ 34]*x11*x21            
    ;
    v_l_r6_h90s_ep25_dex                      =v_l_r6_h90s_ep25_dex                      
        +coeff[ 35]    *x22        *x51
        +coeff[ 36]            *x42*x51
        +coeff[ 37]*x11*x21    *x41    
        +coeff[ 38]    *x21*x31*x42    
        +coeff[ 39]*x12*x21    *x41    
        +coeff[ 40]    *x23*x32*x41    
        +coeff[ 41]*x12                
        +coeff[ 42]            *x41*x52
        +coeff[ 43]*x11*x21*x31        
    ;
    v_l_r6_h90s_ep25_dex                      =v_l_r6_h90s_ep25_dex                      
        +coeff[ 44]*x11    *x32        
        +coeff[ 45]    *x24            
        +coeff[ 46]    *x21*x32*x41    
        +coeff[ 47]    *x23        *x51
        +coeff[ 48]    *x22    *x43    
        +coeff[ 49]    *x21    *x44    
        +coeff[ 50]    *x24*x31*x42    
        +coeff[ 51]    *x23    *x44    
        ;

    return v_l_r6_h90s_ep25_dex                      ;
}
float x_r5p65_ep27_q3en                       (float *x,int m){
    int ncoeff= 56;
    float avdat= -0.4821754E-03;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 57]={
         0.17017144E-02,-0.21273466E-01, 0.13856858E+00, 0.18992445E+00,
        -0.59317835E-02, 0.29152339E-01, 0.20497121E-01,-0.59138052E-01,
        -0.22004616E-03, 0.17559230E-01, 0.23701756E-03,-0.49201269E-02,
        -0.85991593E-02,-0.15416239E-01, 0.69894576E-02,-0.27045673E-01,
         0.29784968E-01,-0.13668934E-02,-0.22755114E-02,-0.52418532E-02,
        -0.11019459E-01,-0.13564969E-02, 0.26349879E-02,-0.47627767E-02,
        -0.16618176E-04, 0.11610240E-01, 0.65027163E-02, 0.11889408E-02,
         0.19891190E-02,-0.24924353E-02,-0.13355545E-02, 0.26587655E-02,
         0.26788455E-02,-0.14852559E-03,-0.30850101E-03,-0.13524516E-02,
        -0.40340589E-03,-0.11259356E-02,-0.76945178E-03,-0.65641932E-03,
         0.13627097E-02, 0.22896128E-02,-0.42283717E-02,-0.52384981E-02,
         0.12186525E-03, 0.17593855E-03,-0.15221120E-03,-0.14551691E-03,
         0.36536576E-03, 0.26135886E-03, 0.36672098E-03,-0.78408368E-03,
         0.69885154E-03,-0.10869339E-02,-0.42343859E-03,-0.52960677E-03,
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

    float v_x_r5p65_ep27_q3en                       =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]                *x51
        +coeff[  3]    *x21            
        +coeff[  4]                *x52
        +coeff[  5]    *x21        *x51
        +coeff[  6]    *x22            
        +coeff[  7]    *x21    *x41    
    ;
    v_x_r5p65_ep27_q3en                       =v_x_r5p65_ep27_q3en                       
        +coeff[  8]            *x41*x52
        +coeff[  9]    *x23            
        +coeff[ 10]*x12*x21*x31        
        +coeff[ 11]            *x41    
        +coeff[ 12]*x11        *x41    
        +coeff[ 13]    *x21*x31        
        +coeff[ 14]*x11*x22            
        +coeff[ 15]    *x21    *x42    
        +coeff[ 16]    *x23    *x41    
    ;
    v_x_r5p65_ep27_q3en                       =v_x_r5p65_ep27_q3en                       
        +coeff[ 17]        *x31        
        +coeff[ 18]*x11    *x31        
        +coeff[ 19]            *x42    
        +coeff[ 20]    *x21*x31*x41    
        +coeff[ 21]        *x31*x41    
        +coeff[ 22]    *x22        *x51
        +coeff[ 23]    *x21    *x41*x51
        +coeff[ 24]            *x41*x53
        +coeff[ 25]*x11*x22    *x41    
    ;
    v_x_r5p65_ep27_q3en                       =v_x_r5p65_ep27_q3en                       
        +coeff[ 26]    *x23*x31        
        +coeff[ 27]*x12*x21            
        +coeff[ 28]*x11*x21    *x41    
        +coeff[ 29]*x11        *x42    
        +coeff[ 30]    *x21*x32        
        +coeff[ 31]    *x22    *x41    
        +coeff[ 32]*x11*x22*x31        
        +coeff[ 33]*x13    *x31*x41    
        +coeff[ 34]*x11*x21            
    ;
    v_x_r5p65_ep27_q3en                       =v_x_r5p65_ep27_q3en                       
        +coeff[ 35]            *x41*x51
        +coeff[ 36]    *x21        *x52
        +coeff[ 37]*x11    *x31*x41    
        +coeff[ 38]    *x21*x31    *x51
        +coeff[ 39]            *x42*x51
        +coeff[ 40]*x12*x21    *x41    
        +coeff[ 41]    *x23        *x51
        +coeff[ 42]    *x21*x31*x42    
        +coeff[ 43]    *x21    *x43    
    ;
    v_x_r5p65_ep27_q3en                       =v_x_r5p65_ep27_q3en                       
        +coeff[ 44]*x11*x22*x32    *x51
        +coeff[ 45]*x12                
        +coeff[ 46]        *x31    *x51
        +coeff[ 47]        *x32        
        +coeff[ 48]                *x53
        +coeff[ 49]*x12        *x41    
        +coeff[ 50]*x11*x21*x31        
        +coeff[ 51]            *x43    
        +coeff[ 52]*x11*x22        *x51
    ;
    v_x_r5p65_ep27_q3en                       =v_x_r5p65_ep27_q3en                       
        +coeff[ 53]    *x21*x32*x41    
        +coeff[ 54]*x13        *x41*x51
        +coeff[ 55]        *x33*x42    
        ;

    return v_x_r5p65_ep27_q3en                       ;
}
float t_r5p65_ep27_q3en                       (float *x,int m){
    int ncoeff= 80;
    float avdat=  0.5124495E-03;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 81]={
        -0.89949614E-03,-0.69620445E-01, 0.52552493E-02, 0.16819387E-02,
         0.22610147E-01,-0.52898075E-02, 0.18588293E-01,-0.39065690E-02,
        -0.54532783E-02,-0.19288927E-02, 0.47825738E-02, 0.26284559E-02,
        -0.13866472E-02,-0.24270522E-02, 0.71363337E-02, 0.44686237E-03,
         0.70066075E-03, 0.35241495E-02,-0.84451651E-02, 0.57909812E-03,
        -0.38270812E-03,-0.34128323E-02, 0.30540273E-03, 0.16496862E-03,
        -0.37583505E-03, 0.40676797E-03,-0.62070391E-03, 0.48955064E-03,
        -0.64825092E-03, 0.91588934E-03,-0.82733162E-03, 0.23260468E-02,
         0.30233563E-04,-0.34947411E-03,-0.80205238E-04, 0.30523996E-04,
        -0.24347789E-03, 0.22357295E-04,-0.13079520E-03, 0.39799791E-03,
         0.32581080E-03, 0.51829458E-03, 0.16821172E-02, 0.34335164E-04,
        -0.40078965E-04,-0.95709154E-04,-0.44339403E-03, 0.63668733E-04,
         0.27225699E-03, 0.11795543E-03, 0.12471018E-03,-0.44192775E-03,
         0.39640901E-03,-0.45662691E-03, 0.15175299E-03, 0.95844903E-03,
         0.17973396E-03,-0.16437624E-03,-0.23844764E-03, 0.24737323E-04,
        -0.12322883E-03, 0.99647499E-04,-0.41741715E-03, 0.12216525E-03,
        -0.16983057E-03,-0.92408758E-04,-0.22478495E-03,-0.14532209E-03,
         0.23561146E-03,-0.39719229E-03, 0.21009185E-03,-0.38148981E-03,
        -0.37582970E-03,-0.30781361E-03,-0.34588197E-03, 0.20889005E-03,
         0.10635910E-03,-0.92962396E-03, 0.16032183E-03,-0.19019768E-04,
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
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_r5p65_ep27_q3en                       =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]*x11                
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21    *x41    
        +coeff[  7]    *x21        *x51
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[  8]    *x23            
        +coeff[  9]    *x23*x31        
        +coeff[ 10]    *x21*x31        
        +coeff[ 11]*x11        *x41    
        +coeff[ 12]                *x52
        +coeff[ 13]*x11*x22            
        +coeff[ 14]    *x21    *x42    
        +coeff[ 15]        *x31        
        +coeff[ 16]*x11    *x31        
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[ 17]    *x21*x31*x41    
        +coeff[ 18]    *x23    *x41    
        +coeff[ 19]            *x42    
        +coeff[ 20]*x11            *x51
        +coeff[ 21]*x11*x22    *x41    
        +coeff[ 22]*x11*x21            
        +coeff[ 23]        *x31*x41    
        +coeff[ 24]*x12*x21            
        +coeff[ 25]    *x21*x32        
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[ 26]*x11*x21    *x41    
        +coeff[ 27]    *x21        *x52
        +coeff[ 28]    *x24            
        +coeff[ 29]*x11        *x42    
        +coeff[ 30]*x11*x22*x31        
        +coeff[ 31]    *x21    *x43    
        +coeff[ 32]*x13    *x31*x41    
        +coeff[ 33]    *x23*x31*x42    
        +coeff[ 34]*x12                
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[ 35]    *x22*x31        
        +coeff[ 36]    *x22    *x41    
        +coeff[ 37]    *x22        *x51
        +coeff[ 38]*x11*x21*x31        
        +coeff[ 39]*x11    *x31*x41    
        +coeff[ 40]            *x42*x51
        +coeff[ 41]    *x22    *x42    
        +coeff[ 42]    *x21*x31*x42    
        +coeff[ 43]        *x32        
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[ 44]        *x31    *x51
        +coeff[ 45]*x11*x21        *x51
        +coeff[ 46]    *x21    *x41*x51
        +coeff[ 47]*x11    *x32        
        +coeff[ 48]            *x43    
        +coeff[ 49]            *x41*x52
        +coeff[ 50]                *x53
        +coeff[ 51]*x12*x21    *x41    
        +coeff[ 52]    *x21*x32*x41    
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[ 53]    *x21    *x42*x51
        +coeff[ 54]    *x23*x31    *x51
        +coeff[ 55]    *x23    *x41*x51
        +coeff[ 56]    *x22*x31*x42    
        +coeff[ 57]*x12        *x43    
        +coeff[ 58]    *x23    *x42*x51
        +coeff[ 59]            *x41*x51
        +coeff[ 60]    *x21*x31    *x51
        +coeff[ 61]        *x31*x42    
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[ 62]    *x23        *x51
        +coeff[ 63]    *x22*x31*x41    
        +coeff[ 64]    *x22        *x52
        +coeff[ 65]*x12*x21*x31        
        +coeff[ 66]*x11*x21    *x42    
        +coeff[ 67]    *x21*x31*x41*x51
        +coeff[ 68]*x11*x24            
        +coeff[ 69]    *x24    *x41    
        +coeff[ 70]*x11        *x43    
    ;
    v_t_r5p65_ep27_q3en                       =v_t_r5p65_ep27_q3en                       
        +coeff[ 71]    *x24        *x51
        +coeff[ 72]    *x23*x31*x41    
        +coeff[ 73]    *x23        *x52
        +coeff[ 74]*x11*x22    *x42    
        +coeff[ 75]    *x22    *x41*x52
        +coeff[ 76]    *x24*x32        
        +coeff[ 77]    *x23    *x43    
        +coeff[ 78]*x13    *x31*x42    
        +coeff[ 79]*x13                
        ;

    return v_t_r5p65_ep27_q3en                       ;
}
float y_r5p65_ep27_q3en                       (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.3732981E-02;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
         0.13151467E-02, 0.67382812E-01,-0.15446213E-01, 0.23840152E-01,
        -0.65819845E-01,-0.38219169E-01,-0.54136431E-02, 0.39648917E-02,
         0.63369773E-01,-0.11496970E-01,-0.10644822E-01, 0.61538992E-02,
        -0.40749530E-02, 0.10896753E-01, 0.72544934E-02,-0.49492348E-01,
         0.50734438E-03, 0.15966650E-01,-0.58616125E-02, 0.11489292E-03,
        -0.13768747E-01,-0.20521651E-03, 0.59119362E-03, 0.15618074E-01,
         0.71706814E-02,-0.42266296E-02,-0.52319863E-02,-0.25284181E-02,
         0.13850410E-01,-0.41306438E-02,-0.14882846E-02, 0.81747054E-03,
        -0.85860753E-03, 0.23113042E-02,-0.79266801E-02,-0.14722147E-02,
        -0.37610685E-02,-0.33900891E-02, 0.50776708E-02, 0.50157221E-03,
        -0.14741353E-02, 0.51781256E-03, 0.23160800E-02, 0.79490169E-03,
         0.10875888E-02, 0.35389953E-02, 0.76170918E-02, 0.54537570E-02,
        -0.98675641E-03,-0.15027529E-02,-0.18085662E-02,-0.17739919E-02,
        -0.85000852E-02, 0.12879625E-02,-0.60885972E-02, 0.19293409E-04,
         0.12764987E-02, 0.14579014E-03,-0.14530504E-03,-0.74180373E-03,
        -0.56382193E-03,-0.53538562E-03, 0.13345978E-02, 0.48013584E-03,
        -0.81198267E-03, 0.10974103E-02,-0.11138470E-02, 0.81095955E-03,
        -0.32267657E-02,-0.20810347E-02, 0.20871116E-03, 0.34240356E-02,
        -0.17397966E-03, 0.75826648E-03, 0.44179862E-03, 0.13650233E-02,
         0.16316706E-02,-0.25052845E-03,-0.29779779E-03, 0.53300097E-03,
        -0.50537917E-03,-0.39021650E-02,-0.19868626E-02, 0.41019177E-03,
        -0.91255672E-03, 0.14599401E-03, 0.61801198E-03, 0.79520934E-04,
        -0.13024801E-03,-0.19151304E-03, 0.32305162E-03, 0.33968769E-03,
        -0.32942303E-03,-0.44667381E-02,-0.15911441E-03,-0.89110079E-03,
         0.64991554E-03,
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

    float v_y_r5p65_ep27_q3en                       =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]        *x31        
        +coeff[  3]                *x51
        +coeff[  4]    *x21    *x41    
        +coeff[  5]    *x22            
        +coeff[  6]    *x21*x31        
        +coeff[  7]*x11        *x41    
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[  8]            *x41*x51
        +coeff[  9]*x11*x21            
        +coeff[ 10]    *x21        *x51
        +coeff[ 11]        *x31    *x51
        +coeff[ 12]            *x43    
        +coeff[ 13]    *x21    *x42    
        +coeff[ 14]                *x52
        +coeff[ 15]    *x22    *x41    
        +coeff[ 16]            *x42*x51
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 17]    *x23            
        +coeff[ 18]    *x21    *x41*x51
        +coeff[ 19]            *x44    
        +coeff[ 20]    *x22        *x51
        +coeff[ 21]        *x31*x43    
        +coeff[ 22]    *x24*x31        
        +coeff[ 23]            *x42    
        +coeff[ 24]        *x31*x41    
        +coeff[ 25]*x11*x21    *x41    
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 26]            *x41*x52
        +coeff[ 27]*x11*x21        *x51
        +coeff[ 28]    *x24            
        +coeff[ 29]    *x22    *x41*x51
        +coeff[ 30]    *x21            
        +coeff[ 31]        *x32        
        +coeff[ 32]*x12                
        +coeff[ 33]    *x21*x31*x41    
        +coeff[ 34]    *x22*x31        
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 35]*x11*x21*x31        
        +coeff[ 36]    *x22*x31*x41    
        +coeff[ 37]*x11*x21    *x42    
        +coeff[ 38]*x11*x23            
        +coeff[ 39]*x11*x24            
        +coeff[ 40]*x11                
        +coeff[ 41]*x11            *x51
        +coeff[ 42]*x11        *x42    
        +coeff[ 43]*x11    *x31*x41    
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 44]        *x31*x41*x51
        +coeff[ 45]*x11*x22            
        +coeff[ 46]    *x21    *x43    
        +coeff[ 47]    *x23    *x41    
        +coeff[ 48]                *x53
        +coeff[ 49]*x11*x22    *x41    
        +coeff[ 50]*x11*x21*x31*x41    
        +coeff[ 51]*x11*x21    *x41*x51
        +coeff[ 52]    *x23    *x42    
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 53]    *x22        *x52
        +coeff[ 54]    *x22    *x42*x51
        +coeff[ 55]    *x23*x31*x42    
        +coeff[ 56]    *x24    *x43    
        +coeff[ 57]*x11    *x31        
        +coeff[ 58]        *x31*x42    
        +coeff[ 59]*x12        *x41    
        +coeff[ 60]    *x21        *x52
        +coeff[ 61]        *x31    *x52
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 62]            *x43*x51
        +coeff[ 63]        *x31*x42*x51
        +coeff[ 64]    *x22*x32        
        +coeff[ 65]    *x23        *x51
        +coeff[ 66]    *x22*x31    *x51
        +coeff[ 67]*x12*x22            
        +coeff[ 68]*x11*x22    *x42    
        +coeff[ 69]    *x22*x31*x41*x51
        +coeff[ 70]        *x32*x41    
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 71]    *x21*x31*x42    
        +coeff[ 72]*x12    *x31        
        +coeff[ 73]*x11        *x43    
        +coeff[ 74]*x11    *x31*x42    
        +coeff[ 75]    *x21    *x42*x51
        +coeff[ 76]    *x23*x31        
        +coeff[ 77]*x12            *x51
        +coeff[ 78]*x11*x21*x32        
        +coeff[ 79]    *x21    *x41*x52
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 80]*x11*x21*x31    *x51
        +coeff[ 81]    *x24    *x41    
        +coeff[ 82]    *x23*x31*x41    
        +coeff[ 83]*x11*x21        *x52
        +coeff[ 84]*x11*x22*x31*x41    
        +coeff[ 85]        *x34    *x51
        +coeff[ 86]*x11*x23        *x51
        +coeff[ 87]*x11    *x32        
        +coeff[ 88]    *x21*x31    *x51
    ;
    v_y_r5p65_ep27_q3en                       =v_y_r5p65_ep27_q3en                       
        +coeff[ 89]    *x22    *x42    
        +coeff[ 90]    *x21*x32*x41    
        +coeff[ 91]    *x21*x31*x41*x51
        +coeff[ 92]*x12        *x42    
        +coeff[ 93]    *x22    *x43    
        +coeff[ 94]*x12    *x31*x41    
        +coeff[ 95]    *x22*x31*x42    
        +coeff[ 96]*x11*x21    *x43    
        ;

    return v_y_r5p65_ep27_q3en                       ;
}
float p_r5p65_ep27_q3en                       (float *x,int m){
    int ncoeff= 90;
    float avdat=  0.5064601E-03;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 91]={
         0.27119872E-03, 0.15072488E-03,-0.41751061E-02, 0.44030853E-03,
         0.29025953E-02,-0.46768188E-02,-0.11492816E-02,-0.12682733E-01,
         0.22145249E-02,-0.19080935E-02, 0.10456054E-02, 0.12667133E-01,
        -0.15299680E-02, 0.32206795E-02, 0.16770373E-02,-0.10965647E-02,
         0.87133743E-03,-0.94812037E-02, 0.25462485E-02,-0.32731763E-02,
        -0.14356871E-02, 0.24551132E-02,-0.21779544E-03,-0.51408907E-03,
        -0.13581379E-02,-0.30278568E-05, 0.98751206E-03, 0.68552740E-03,
        -0.92321710E-03,-0.48716835E-03, 0.75821194E-03,-0.97568576E-04,
        -0.22164852E-03, 0.78423283E-04, 0.28880988E-03,-0.12449303E-03,
        -0.18198696E-03, 0.46062560E-03, 0.11416507E-02,-0.39666146E-03,
         0.25499139E-06,-0.76473906E-03, 0.97801076E-05, 0.59860467E-04,
        -0.10107296E-03, 0.21537967E-03, 0.71518874E-03,-0.65723827E-04,
         0.15000495E-03, 0.24887404E-03,-0.14667347E-03, 0.13897904E-02,
         0.10167228E-03,-0.12683171E-03,-0.76808530E-03,-0.32394126E-03,
        -0.19742106E-02,-0.33838485E-03,-0.95789396E-03, 0.82379917E-03,
        -0.58580365E-03,-0.66947898E-04, 0.10859541E-03,-0.17620590E-03,
         0.37152442E-03,-0.96365846E-04, 0.70526556E-03,-0.48658127E-03,
         0.13725953E-03,-0.71002032E-04, 0.21448339E-04,-0.20890388E-03,
        -0.47767838E-03,-0.62700827E-04,-0.33665248E-03,-0.63902189E-04,
         0.44819064E-03,-0.36319587E-03, 0.11176791E-03, 0.23989356E-03,
        -0.24766868E-03,-0.92476708E-04, 0.31541727E-03, 0.61809685E-03,
        -0.16885241E-03,-0.24048866E-03, 0.28934774E-04, 0.87402623E-04,
        -0.26340796E-04,-0.76754477E-04,
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
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;
    float x54 = x53*x5;

//                 function

    float v_p_r5p65_ep27_q3en                       =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21*x31        
        +coeff[  7]    *x21    *x41    
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[  8]            *x42    
        +coeff[  9]    *x21        *x51
        +coeff[ 10]        *x31    *x51
        +coeff[ 11]            *x41*x51
        +coeff[ 12]*x11*x21            
        +coeff[ 13]    *x23            
        +coeff[ 14]                *x52
        +coeff[ 15]    *x22*x31        
        +coeff[ 16]*x11        *x41    
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 17]    *x22    *x41    
        +coeff[ 18]    *x21    *x42    
        +coeff[ 19]    *x22        *x51
        +coeff[ 20]    *x21    *x41*x51
        +coeff[ 21]    *x24            
        +coeff[ 22]    *x22*x31*x41    
        +coeff[ 23]*x11*x21        *x51
        +coeff[ 24]    *x22    *x41*x51
        +coeff[ 25]        *x33*x41    
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 26]        *x31*x41    
        +coeff[ 27]    *x21*x31*x41    
        +coeff[ 28]            *x43    
        +coeff[ 29]            *x41*x52
        +coeff[ 30]*x11*x23            
        +coeff[ 31]*x11*x24            
        +coeff[ 32]*x11                
        +coeff[ 33]*x11    *x31        
        +coeff[ 34]        *x31*x41*x51
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 35]*x12                
        +coeff[ 36]*x11*x21*x31        
        +coeff[ 37]*x11        *x42    
        +coeff[ 38]    *x22    *x42    
        +coeff[ 39]    *x22*x31    *x51
        +coeff[ 40]        *x34        
        +coeff[ 41]*x11*x22    *x41    
        +coeff[ 42]    *x21        *x54
        +coeff[ 43]*x11            *x51
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 44]        *x31*x42    
        +coeff[ 45]            *x42*x51
        +coeff[ 46]*x11*x22            
        +coeff[ 47]*x11*x21    *x41    
        +coeff[ 48]*x11    *x31*x41    
        +coeff[ 49]    *x23        *x51
        +coeff[ 50]                *x53
        +coeff[ 51]    *x21    *x43    
        +coeff[ 52]    *x22        *x52
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 53]*x12        *x41    
        +coeff[ 54]    *x24    *x41    
        +coeff[ 55]*x11*x21    *x42    
        +coeff[ 56]    *x23    *x42    
        +coeff[ 57]*x11*x21    *x41*x51
        +coeff[ 58]    *x22    *x42*x51
        +coeff[ 59]    *x25    *x41    
        +coeff[ 60]*x11*x22    *x42    
        +coeff[ 61]    *x23*x31*x42    
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 62]        *x32        
        +coeff[ 63]    *x21        *x52
        +coeff[ 64]    *x23    *x41    
        +coeff[ 65]    *x22*x32        
        +coeff[ 66]    *x21*x31*x42    
        +coeff[ 67]    *x25            
        +coeff[ 68]            *x43*x51
        +coeff[ 69]    *x21*x31    *x52
        +coeff[ 70]    *x23*x32        
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 71]*x11*x21*x31*x41    
        +coeff[ 72]    *x23*x31*x41    
        +coeff[ 73]*x12            *x51
        +coeff[ 74]    *x22    *x43    
        +coeff[ 75]    *x23*x31    *x51
        +coeff[ 76]    *x23    *x41*x51
        +coeff[ 77]    *x22*x31*x41*x51
        +coeff[ 78]*x12*x22            
        +coeff[ 79]    *x25*x31        
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 80]*x11*x23    *x41    
        +coeff[ 81]    *x22*x31    *x52
        +coeff[ 82]    *x22    *x41*x52
        +coeff[ 83]*x11*x24    *x41    
        +coeff[ 84]*x11*x23*x31    *x51
        +coeff[ 85]*x11*x22*x33*x41    
        +coeff[ 86]        *x32    *x51
        +coeff[ 87]    *x21*x32*x41    
        +coeff[ 88]*x11        *x41*x51
    ;
    v_p_r5p65_ep27_q3en                       =v_p_r5p65_ep27_q3en                       
        +coeff[ 89]        *x31*x43    
        ;

    return v_p_r5p65_ep27_q3en                       ;
}
float l_r6_h90s_ep27_q3en                     (float *x,int m){
    int ncoeff= 59;
    float avdat= -0.5614310E-02;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99670E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 60]={
         0.64344076E-03,-0.27928358E+00,-0.32212082E-01, 0.24242336E-01,
        -0.30473128E-01, 0.79205185E-01,-0.15997348E-01, 0.11516576E-01,
        -0.23517901E-01,-0.82159424E-02, 0.30228083E-02, 0.20141939E-01,
         0.30944396E-01,-0.96670175E-02,-0.37057444E-01,-0.21888833E-02,
        -0.42580413E-02, 0.30070762E-02, 0.14877351E-01,-0.15173806E-01,
         0.79378701E-03,-0.18904225E-02, 0.15018963E-02,-0.18778697E-02,
         0.35587628E-02, 0.17589256E-02, 0.50875586E-02,-0.11286605E-02,
         0.87186473E-03,-0.16861580E-02, 0.17788789E-02,-0.16223585E-02,
        -0.23461925E-02, 0.94705410E-02,-0.35197607E-02,-0.55306125E-03,
        -0.34234705E-03, 0.37827043E-03, 0.96518412E-03,-0.42586465E-03,
         0.76195486E-02,-0.18507310E-02,-0.48606916E-03,-0.16855347E-03,
        -0.30700557E-03, 0.23316624E-03,-0.53872698E-03,-0.39874838E-03,
         0.30950637E-03, 0.23752832E-03,-0.25888628E-03,-0.29805035E-03,
         0.19841776E-02, 0.26637213E-02,-0.93676109E-03,-0.38474649E-02,
         0.16933086E-02, 0.28659394E-02,-0.43086126E-02,
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
    float x51 = x5;
    float x52 = x51*x5;

//                 function

    float v_l_r6_h90s_ep27_q3en                     =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]                *x51
        +coeff[  3]*x11                
        +coeff[  4]    *x22            
        +coeff[  5]    *x21    *x41    
        +coeff[  6]    *x21        *x51
        +coeff[  7]*x11        *x41    
    ;
    v_l_r6_h90s_ep27_q3en                     =v_l_r6_h90s_ep27_q3en                     
        +coeff[  8]    *x23            
        +coeff[  9]    *x23*x31        
        +coeff[ 10]            *x41    
        +coeff[ 11]    *x21*x31        
        +coeff[ 12]    *x21    *x42    
        +coeff[ 13]*x11*x22            
        +coeff[ 14]    *x23    *x41    
        +coeff[ 15]    *x23*x31*x41    
        +coeff[ 16]            *x42    
    ;
    v_l_r6_h90s_ep27_q3en                     =v_l_r6_h90s_ep27_q3en                     
        +coeff[ 17]*x11    *x31        
        +coeff[ 18]    *x21*x31*x41    
        +coeff[ 19]*x11*x22    *x41    
        +coeff[ 20]        *x31        
        +coeff[ 21]            *x41*x51
        +coeff[ 22]*x11*x21            
        +coeff[ 23]*x11            *x51
        +coeff[ 24]*x11        *x42    
        +coeff[ 25]    *x21*x32        
    ;
    v_l_r6_h90s_ep27_q3en                     =v_l_r6_h90s_ep27_q3en                     
        +coeff[ 26]    *x22    *x41    
        +coeff[ 27]    *x22        *x51
        +coeff[ 28]    *x21        *x52
        +coeff[ 29]*x11*x21    *x41    
        +coeff[ 30]*x11    *x31*x41    
        +coeff[ 31]*x12*x21            
        +coeff[ 32]    *x24            
        +coeff[ 33]    *x21    *x43    
        +coeff[ 34]*x11*x22*x31        
    ;
    v_l_r6_h90s_ep27_q3en                     =v_l_r6_h90s_ep27_q3en                     
        +coeff[ 35]    *x23*x31*x42    
        +coeff[ 36]*x12                
        +coeff[ 37]    *x22*x31        
        +coeff[ 38]            *x42*x51
        +coeff[ 39]*x11*x21*x31        
        +coeff[ 40]    *x21*x31*x42    
        +coeff[ 41]*x12*x21    *x41    
        +coeff[ 42]    *x23*x32*x41    
        +coeff[ 43]        *x31    *x51
    ;
    v_l_r6_h90s_ep27_q3en                     =v_l_r6_h90s_ep27_q3en                     
        +coeff[ 44]                *x52
        +coeff[ 45]            *x43    
        +coeff[ 46]    *x21*x31    *x51
        +coeff[ 47]    *x21    *x41*x51
        +coeff[ 48]            *x41*x52
        +coeff[ 49]*x11    *x32        
        +coeff[ 50]*x11        *x41*x51
        +coeff[ 51]*x12        *x41    
        +coeff[ 52]    *x21*x32*x41    
    ;
    v_l_r6_h90s_ep27_q3en                     =v_l_r6_h90s_ep27_q3en                     
        +coeff[ 53]    *x22    *x42    
        +coeff[ 54]*x11*x23            
        +coeff[ 55]    *x24    *x41    
        +coeff[ 56]    *x22*x31*x42    
        +coeff[ 57]    *x21    *x44    
        +coeff[ 58]    *x23    *x43    
        ;

    return v_l_r6_h90s_ep27_q3en                     ;
}
float x_r5p65_ep30_q3ex                       (float *x,int m){
    int ncoeff= 54;
    float avdat=  0.5975692E-02;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.46896E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99837E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49644E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 55]={
         0.38901684E-02, 0.32358199E+00, 0.63117698E-01,-0.21791870E-01,
         0.30947709E-01, 0.13465733E-01,-0.31705081E-01,-0.84152643E-05,
        -0.17130176E-01,-0.85458560E-02, 0.95687285E-02,-0.16913519E-02,
        -0.45343563E-02,-0.72826957E-02, 0.35903668E-02,-0.19998465E-01,
        -0.83773071E-02,-0.64781019E-02,-0.19524174E-03, 0.18764172E-01,
         0.15318341E-04,-0.89646236E-03,-0.76548464E-03,-0.78732392E-03,
        -0.19166162E-03,-0.12196316E-02,-0.19815657E-02,-0.16676235E-02,
         0.25166054E-02, 0.21361427E-02, 0.20291889E-02, 0.65113213E-02,
         0.35622129E-02, 0.71810547E-03, 0.83981227E-03,-0.11932317E-02,
        -0.80095109E-03,-0.12762116E-05, 0.16435661E-02, 0.22571194E-02,
        -0.15471209E-02, 0.37661619E-02,-0.35194459E-03,-0.65627124E-03,
        -0.87613479E-03,-0.12620106E-03, 0.46947977E-03,-0.73541101E-03,
         0.11360754E-02,-0.65925479E-03, 0.13535495E-02,-0.56318688E-03,
         0.23796114E-02, 0.18913638E-02,
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
    float x54 = x53*x5;

//                 function

    float v_x_r5p65_ep30_q3ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]                *x51
        +coeff[  2]    *x21            
        +coeff[  3]                *x52
        +coeff[  4]    *x21        *x51
        +coeff[  5]    *x22            
        +coeff[  6]    *x21    *x41    
        +coeff[  7]*x13                
    ;
    v_x_r5p65_ep30_q3ex                       =v_x_r5p65_ep30_q3ex                       
        +coeff[  8]*x11                
        +coeff[  9]    *x21*x31        
        +coeff[ 10]    *x23            
        +coeff[ 11]            *x41    
        +coeff[ 12]*x11        *x41    
        +coeff[ 13]            *x42    
        +coeff[ 14]*x11*x22            
        +coeff[ 15]    *x21    *x42    
        +coeff[ 16]    *x21    *x41*x51
    ;
    v_x_r5p65_ep30_q3ex                       =v_x_r5p65_ep30_q3ex                       
        +coeff[ 17]    *x21*x31*x41    
        +coeff[ 18]            *x41*x53
        +coeff[ 19]    *x23    *x41    
        +coeff[ 20]    *x22*x31*x41    
        +coeff[ 21]    *x21        *x54
        +coeff[ 22]        *x31        
        +coeff[ 23]*x11            *x51
        +coeff[ 24]*x11*x21            
        +coeff[ 25]*x11    *x31        
    ;
    v_x_r5p65_ep30_q3ex                       =v_x_r5p65_ep30_q3ex                       
        +coeff[ 26]            *x41*x51
        +coeff[ 27]        *x31*x41    
        +coeff[ 28]                *x53
        +coeff[ 29]    *x22        *x51
        +coeff[ 30]    *x22    *x41    
        +coeff[ 31]*x11*x22    *x41    
        +coeff[ 32]    *x23*x31        
        +coeff[ 33]*x12*x21            
        +coeff[ 34]*x11*x21    *x41    
    ;
    v_x_r5p65_ep30_q3ex                       =v_x_r5p65_ep30_q3ex                       
        +coeff[ 35]    *x21*x31    *x51
        +coeff[ 36]    *x21*x32        
        +coeff[ 37]*x12    *x31    *x51
        +coeff[ 38]*x11*x23            
        +coeff[ 39]    *x23        *x51
        +coeff[ 40]    *x21    *x42*x51
        +coeff[ 41]    *x22    *x42    
        +coeff[ 42]        *x31    *x51
        +coeff[ 43]*x11    *x31*x41    
    ;
    v_x_r5p65_ep30_q3ex                       =v_x_r5p65_ep30_q3ex                       
        +coeff[ 44]*x11        *x42    
        +coeff[ 45]            *x42*x51
        +coeff[ 46]    *x22*x31        
        +coeff[ 47]            *x43    
        +coeff[ 48]*x11*x22        *x51
        +coeff[ 49]    *x21    *x41*x52
        +coeff[ 50]*x11*x22*x31        
        +coeff[ 51]*x13        *x41*x51
        +coeff[ 52]    *x23    *x41*x51
    ;
    v_x_r5p65_ep30_q3ex                       =v_x_r5p65_ep30_q3ex                       
        +coeff[ 53]    *x22    *x42*x51
        ;

    return v_x_r5p65_ep30_q3ex                       ;
}
float t_r5p65_ep30_q3ex                       (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.2479329E-02;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.46896E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99837E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49644E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 51]={
         0.63250965E-03,-0.19915611E-01,-0.22842800E-02, 0.45461624E-03,
         0.11148836E+00, 0.59229392E-02,-0.10317651E-01, 0.17554996E-02,
         0.10522753E-02,-0.18974107E-02,-0.20552033E-02,-0.22142404E-02,
        -0.59869711E-03, 0.16187732E-03,-0.40223417E-03,-0.64957218E-03,
        -0.11446449E-02, 0.10942066E-02, 0.11743134E-02, 0.61968822E-04,
         0.15858719E-03, 0.15427097E-03,-0.20601245E-03,-0.20826903E-03,
         0.17976726E-03, 0.26147798E-03,-0.10703380E-03, 0.30437272E-03,
         0.74822264E-03,-0.22430508E-03, 0.21095257E-03,-0.14879792E-03,
         0.37206226E-03, 0.71508973E-03, 0.21692817E-04, 0.25848291E-03,
        -0.85052487E-03, 0.10593899E-02,-0.37588077E-04,-0.69601869E-04,
        -0.79951264E-04, 0.45775145E-03,-0.24638162E-03,-0.13710569E-03,
         0.33055386E-03,-0.24301175E-03, 0.13705554E-03, 0.28300460E-03,
        -0.89929177E-03,-0.56325807E-03,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_t_r5p65_ep30_q3ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]*x11                
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x21        *x51
        +coeff[  6]                *x52
        +coeff[  7]    *x22            
    ;
    v_t_r5p65_ep30_q3ex                       =v_t_r5p65_ep30_q3ex                       
        +coeff[  8]    *x21    *x41    
        +coeff[  9]            *x42    
        +coeff[ 10]    *x21    *x42    
        +coeff[ 11]    *x21    *x41*x51
        +coeff[ 12]    *x21        *x52
        +coeff[ 13]*x11*x21            
        +coeff[ 14]        *x31*x41    
        +coeff[ 15]            *x41*x51
        +coeff[ 16]    *x24            
    ;
    v_t_r5p65_ep30_q3ex                       =v_t_r5p65_ep30_q3ex                       
        +coeff[ 17]                *x53
        +coeff[ 18]    *x22    *x42    
        +coeff[ 19]        *x31        
        +coeff[ 20]    *x21*x31        
        +coeff[ 21]*x11        *x41    
        +coeff[ 22]*x11            *x51
        +coeff[ 23]*x11*x22            
        +coeff[ 24]            *x42*x51
        +coeff[ 25]    *x24*x31        
    ;
    v_t_r5p65_ep30_q3ex                       =v_t_r5p65_ep30_q3ex                       
        +coeff[ 26]        *x31    *x51
        +coeff[ 27]    *x22    *x41    
        +coeff[ 28]    *x22        *x51
        +coeff[ 29]    *x21*x31    *x51
        +coeff[ 30]*x11        *x42    
        +coeff[ 31]            *x43    
        +coeff[ 32]*x11*x23            
        +coeff[ 33]    *x23    *x41    
        +coeff[ 34]    *x23        *x51
    ;
    v_t_r5p65_ep30_q3ex                       =v_t_r5p65_ep30_q3ex                       
        +coeff[ 35]*x11*x22        *x51
        +coeff[ 36]    *x21    *x42*x51
        +coeff[ 37]    *x23    *x41*x51
        +coeff[ 38]*x12                
        +coeff[ 39]    *x23            
        +coeff[ 40]        *x31*x42    
        +coeff[ 41]    *x22    *x41*x51
        +coeff[ 42]    *x22        *x52
        +coeff[ 43]*x11*x21    *x42    
    ;
    v_t_r5p65_ep30_q3ex                       =v_t_r5p65_ep30_q3ex                       
        +coeff[ 44]    *x21    *x43    
        +coeff[ 45]    *x21*x31*x41*x51
        +coeff[ 46]    *x21        *x53
        +coeff[ 47]    *x24    *x41    
        +coeff[ 48]    *x24        *x51
        +coeff[ 49]    *x23        *x52
        ;

    return v_t_r5p65_ep30_q3ex                       ;
}
float y_r5p65_ep30_q3ex                       (float *x,int m){
    int ncoeff= 97;
    float avdat= -0.9482167E-04;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.46896E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99837E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49644E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
         0.12713786E-02, 0.26980892E-01,-0.19941127E-03,-0.14447949E-01,
        -0.10027437E-02, 0.14865476E-01, 0.10537067E-01,-0.51537436E-01,
         0.48216055E-02,-0.24351226E-01,-0.45447652E-02, 0.30887628E-02,
         0.50658319E-01,-0.75312941E-02,-0.80189072E-02, 0.36852094E-02,
        -0.33962331E-02, 0.65026805E-02,-0.38554627E-01, 0.11690071E-01,
        -0.54247985E-02,-0.71624937E-02,-0.13591558E-01, 0.10409278E-01,
        -0.64581288E-02, 0.46443156E-03, 0.89543751E-02, 0.22997980E-02,
         0.28976570E-02,-0.24462931E-02, 0.34545441E-02, 0.53016254E-03,
        -0.28232869E-03,-0.54911239E-03, 0.18516443E-02,-0.18620624E-02,
        -0.93988911E-03,-0.98363787E-03, 0.23224715E-02,-0.12931093E-02,
         0.19199615E-03, 0.31982193E-03, 0.61896228E-03, 0.11867899E-02,
         0.55131228E-02,-0.52685902E-03, 0.26490339E-02,-0.10342046E-02,
         0.42677517E-02,-0.20660353E-02, 0.11189656E-02,-0.61689113E-03,
         0.19843644E-02,-0.13242228E-02,-0.15554286E-02, 0.54115790E-04,
        -0.67368443E-02,-0.81689493E-03,-0.26020193E-02,-0.35641186E-02,
        -0.42846476E-03, 0.18346039E-03, 0.67514903E-03,-0.18107852E-02,
         0.11287912E-02,-0.51451527E-03,-0.25525867E-03,-0.11348826E-02,
         0.49046002E-03,-0.44643786E-03,-0.29901974E-02,-0.17318627E-02,
        -0.13979933E-02, 0.14843834E-02, 0.76304527E-03, 0.70667424E-03,
         0.75462117E-03, 0.14704125E-03, 0.10575323E-02,-0.13651897E-03,
         0.14603529E-03,-0.32830951E-03, 0.89398629E-04,-0.12143004E-03,
         0.46704168E-03, 0.28842853E-03, 0.34932126E-03, 0.29954544E-03,
        -0.24208972E-03,-0.19060659E-03,-0.25816453E-02, 0.33300702E-03,
         0.66885940E-03, 0.48148056E-03, 0.20858673E-03,-0.51401090E-03,
        -0.67694875E-03,
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
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_r5p65_ep30_q3ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]*x11                
        +coeff[  5]                *x51
        +coeff[  6]            *x42    
        +coeff[  7]    *x21    *x41    
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[  8]        *x31*x41    
        +coeff[  9]    *x22            
        +coeff[ 10]    *x21*x31        
        +coeff[ 11]*x11        *x41    
        +coeff[ 12]            *x41*x51
        +coeff[ 13]*x11*x21            
        +coeff[ 14]    *x21        *x51
        +coeff[ 15]        *x31    *x51
        +coeff[ 16]            *x43    
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 17]                *x52
        +coeff[ 18]    *x22    *x41    
        +coeff[ 19]    *x23            
        +coeff[ 20]    *x22*x31        
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]    *x22        *x51
        +coeff[ 23]    *x24            
        +coeff[ 24]    *x22    *x41*x51
        +coeff[ 25]    *x21    *x44    
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 26]    *x21    *x42    
        +coeff[ 27]    *x21*x31*x41    
        +coeff[ 28]*x11*x22            
        +coeff[ 29]*x11*x21        *x51
        +coeff[ 30]*x11*x23            
        +coeff[ 31]        *x32        
        +coeff[ 32]        *x31*x42    
        +coeff[ 33]*x12                
        +coeff[ 34]*x11        *x42    
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 35]*x11*x21    *x41    
        +coeff[ 36]*x11*x21*x31        
        +coeff[ 37]            *x41*x52
        +coeff[ 38]    *x22    *x42    
        +coeff[ 39]*x11*x22    *x41    
        +coeff[ 40]*x11    *x31        
        +coeff[ 41]*x11            *x51
        +coeff[ 42]*x11    *x31*x41    
        +coeff[ 43]        *x31*x41*x51
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 44]    *x21    *x43    
        +coeff[ 45]*x12        *x41    
        +coeff[ 46]    *x21*x31*x42    
        +coeff[ 47]    *x21        *x52
        +coeff[ 48]    *x23    *x41    
        +coeff[ 49]    *x22*x31*x41    
        +coeff[ 50]    *x23*x31        
        +coeff[ 51]                *x53
        +coeff[ 52]    *x23        *x51
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 53]    *x22*x31    *x51
        +coeff[ 54]*x11*x21    *x41*x51
        +coeff[ 55]            *x44*x51
        +coeff[ 56]    *x23    *x42    
        +coeff[ 57]            *x41*x53
        +coeff[ 58]*x11*x22    *x42    
        +coeff[ 59]    *x22    *x42*x51
        +coeff[ 60]*x11*x21    *x44    
        +coeff[ 61]        *x31    *x52
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 62]*x11        *x43    
        +coeff[ 63]*x11*x21    *x42    
        +coeff[ 64]    *x21    *x42*x51
        +coeff[ 65]    *x22*x32        
        +coeff[ 66]*x12            *x51
        +coeff[ 67]*x11*x21*x31*x41    
        +coeff[ 68]*x12*x22            
        +coeff[ 69]*x11*x21*x31    *x51
        +coeff[ 70]    *x24    *x41    
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 71]    *x23*x31*x41    
        +coeff[ 72]    *x22*x31*x41*x51
        +coeff[ 73]    *x24        *x51
        +coeff[ 74]*x11*x23        *x51
        +coeff[ 75]    *x24    *x43    
        +coeff[ 76]    *x24        *x52
        +coeff[ 77]        *x32*x41    
        +coeff[ 78]            *x42*x51
        +coeff[ 79]    *x21*x31    *x51
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 80]        *x32    *x51
        +coeff[ 81]        *x31*x43    
        +coeff[ 82]*x11        *x41*x51
        +coeff[ 83]*x12    *x31        
        +coeff[ 84]            *x43*x51
        +coeff[ 85]    *x21*x32*x41    
        +coeff[ 86]*x11    *x31*x42    
        +coeff[ 87]        *x31*x42*x51
        +coeff[ 88]*x12        *x42    
    ;
    v_y_r5p65_ep30_q3ex                       =v_y_r5p65_ep30_q3ex                       
        +coeff[ 89]*x11*x21*x32        
        +coeff[ 90]    *x22    *x43    
        +coeff[ 91]    *x21    *x41*x52
        +coeff[ 92]*x11*x21    *x43    
        +coeff[ 93]    *x24*x31        
        +coeff[ 94]*x11*x21        *x52
        +coeff[ 95]*x11*x23    *x41    
        +coeff[ 96]*x11*x22*x31*x41    
        ;

    return v_y_r5p65_ep30_q3ex                       ;
}
float p_r5p65_ep30_q3ex                       (float *x,int m){
    int ncoeff= 90;
    float avdat=  0.1811818E-02;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.46896E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99837E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49644E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 91]={
        -0.77757699E-03, 0.46565733E-02,-0.27665894E-01,-0.84101884E-02,
         0.13766836E-01, 0.17612525E-02, 0.21758491E-01,-0.25717770E-02,
         0.34831883E-02,-0.18947933E-01, 0.40861871E-02,-0.61320085E-02,
        -0.17217022E-02,-0.12971839E-02, 0.16932111E-01, 0.37189885E-02,
         0.80562587E-03, 0.25520576E-03, 0.58185156E-04, 0.68764799E-04,
        -0.10794124E-04,-0.19328375E-03,-0.10167353E-03,-0.55258591E-02,
        -0.23892296E-02,-0.36178017E-02,-0.45673097E-02, 0.15772204E-02,
         0.26107999E-02, 0.91914571E-03, 0.28739139E-02, 0.13578187E-02,
         0.33092088E-03,-0.12826709E-02, 0.53173525E-03, 0.55461260E-03,
         0.53364655E-03,-0.17282545E-02, 0.69250993E-04, 0.72255057E-06,
         0.49863308E-03,-0.29967690E-03,-0.69720822E-03,-0.17698412E-03,
         0.22980930E-03, 0.14787114E-02,-0.72404044E-03, 0.13701641E-03,
         0.18216336E-03,-0.24470547E-02, 0.31064203E-03, 0.24269437E-03,
         0.39369482E-03, 0.11363487E-02,-0.99606055E-03,-0.96841541E-04,
        -0.59430942E-04, 0.87883367E-04, 0.31844928E-03,-0.19850599E-03,
        -0.57782722E-03,-0.67352300E-03, 0.28760178E-03,-0.19639060E-03,
        -0.10728809E-02,-0.21790908E-03,-0.51012944E-03, 0.30801843E-02,
         0.45365287E-03, 0.15745625E-02,-0.30836026E-03, 0.94890303E-03,
        -0.63343839E-04, 0.12372366E-03, 0.19079969E-03,-0.18011991E-03,
         0.68271562E-03,-0.37634533E-03, 0.66386012E-03, 0.12790260E-03,
         0.69534732E-03, 0.62565587E-03,-0.25856169E-03,-0.22111551E-03,
         0.11664625E-02, 0.97870194E-04,-0.15588016E-03,-0.13076648E-02,
         0.35097840E-03,-0.11921878E-03,
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
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_p_r5p65_ep30_q3ex                       =avdat
        +coeff[  0]                    
        +coeff[  1]        *x31        
        +coeff[  2]            *x41    
        +coeff[  3]                *x51
        +coeff[  4]    *x22            
        +coeff[  5]    *x21*x31        
        +coeff[  6]    *x21    *x41    
        +coeff[  7]        *x31*x41    
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[  8]    *x21        *x51
        +coeff[  9]            *x41*x51
        +coeff[ 10]*x11*x21            
        +coeff[ 11]    *x23            
        +coeff[ 12]                *x52
        +coeff[ 13]*x11        *x41    
        +coeff[ 14]    *x22    *x41    
        +coeff[ 15]    *x22        *x51
        +coeff[ 16]    *x22    *x42    
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 17]    *x22*x31    *x51
        +coeff[ 18]        *x32*x42    
        +coeff[ 19]            *x44    
        +coeff[ 20]        *x33    *x51
        +coeff[ 21]    *x24*x31        
        +coeff[ 22]    *x22    *x44    
        +coeff[ 23]            *x42    
        +coeff[ 24]        *x31    *x51
        +coeff[ 25]    *x21    *x42    
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 26]    *x24            
        +coeff[ 27]*x11*x21    *x41    
        +coeff[ 28]            *x41*x52
        +coeff[ 29]    *x21            
        +coeff[ 30]    *x22*x31        
        +coeff[ 31]            *x43    
        +coeff[ 32]*x12                
        +coeff[ 33]*x11*x22            
        +coeff[ 34]*x11*x21*x31        
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 35]        *x31    *x52
        +coeff[ 36]*x11*x21        *x51
        +coeff[ 37]*x11*x23            
        +coeff[ 38]*x11            *x52
        +coeff[ 39]*x13                
        +coeff[ 40]*x11                
        +coeff[ 41]        *x32        
        +coeff[ 42]    *x21*x31*x41    
        +coeff[ 43]*x11            *x51
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 44]            *x42*x51
        +coeff[ 45]    *x22*x31*x41    
        +coeff[ 46]*x11        *x42    
        +coeff[ 47]    *x23        *x51
        +coeff[ 48]                *x53
        +coeff[ 49]    *x21    *x43    
        +coeff[ 50]    *x22    *x41*x51
        +coeff[ 51]*x12        *x41    
        +coeff[ 52]*x11*x22    *x41    
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 53]*x11*x21    *x42    
        +coeff[ 54]            *x41*x53
        +coeff[ 55]*x11*x21*x33*x41    
        +coeff[ 56]*x11    *x31        
        +coeff[ 57]        *x31*x42    
        +coeff[ 58]    *x21    *x41*x51
        +coeff[ 59]    *x21        *x52
        +coeff[ 60]    *x23*x31        
        +coeff[ 61]    *x23    *x41    
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 62]    *x22*x32        
        +coeff[ 63]*x11    *x31*x41    
        +coeff[ 64]    *x21*x31*x42    
        +coeff[ 65]    *x21    *x42*x51
        +coeff[ 66]            *x43*x51
        +coeff[ 67]    *x23    *x42    
        +coeff[ 68]*x11*x21    *x41*x51
        +coeff[ 69]    *x22    *x42*x51
        +coeff[ 70]*x12*x22            
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 71]*x11*x22    *x42    
        +coeff[ 72]        *x32*x41    
        +coeff[ 73]*x11        *x41*x51
        +coeff[ 74]        *x31*x43    
        +coeff[ 75]    *x21*x31*x41*x51
        +coeff[ 76]    *x25            
        +coeff[ 77]    *x22        *x52
        +coeff[ 78]    *x24    *x41    
        +coeff[ 79]*x11*x21*x32        
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 80]*x11*x21*x31*x41    
        +coeff[ 81]    *x23*x31*x41    
        +coeff[ 82]            *x42*x52
        +coeff[ 83]*x11        *x43    
        +coeff[ 84]    *x22    *x43    
        +coeff[ 85]*x11*x21*x31    *x51
        +coeff[ 86]        *x31    *x53
        +coeff[ 87]    *x25    *x41    
        +coeff[ 88]    *x22    *x41*x52
    ;
    v_p_r5p65_ep30_q3ex                       =v_p_r5p65_ep30_q3ex                       
        +coeff[ 89]        *x33*x41*x51
        ;

    return v_p_r5p65_ep30_q3ex                       ;
}
float l_r6_h90s_ep30_q3ex                     (float *x,int m){
    int ncoeff= 59;
    float avdat= -0.4454926E-02;
    float xmin[10]={
        -0.99915E-02,-0.47670E-01,-0.99616E-02,-0.29769E-01,-0.46896E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99837E-02, 0.49365E-01, 0.99935E-02, 0.25007E-01, 0.49644E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 60]={
        -0.13641941E-02,-0.27962425E+00,-0.31681750E-01, 0.24213351E-01,
        -0.33856858E-01, 0.80696344E-01, 0.11095006E-01,-0.23890873E-01,
        -0.85747056E-02, 0.31417292E-02, 0.20396680E-01,-0.89636967E-02,
        -0.82514640E-02, 0.31255119E-01,-0.93545942E-02,-0.39678350E-01,
        -0.18205348E-02, 0.79861860E-03,-0.46554371E-02, 0.30284186E-02,
         0.14601247E-01,-0.14662332E-01,-0.21851351E-02, 0.19086043E-02,
        -0.17419913E-02, 0.68765492E-02, 0.36380347E-02, 0.58312109E-03,
         0.17642663E-02, 0.14323267E-02, 0.17073822E-02,-0.16439519E-02,
        -0.28035676E-02, 0.72241067E-02,-0.35885957E-02, 0.30129824E-02,
        -0.37742846E-03,-0.81088510E-03,-0.82799449E-03,-0.11060389E-02,
        -0.15609120E-02, 0.56705028E-02,-0.12678871E-02,-0.17661980E-02,
         0.27628168E-02, 0.17850971E-03,-0.18674730E-03, 0.73162583E-03,
        -0.38605984E-03,-0.38327998E-03,-0.30472185E-03, 0.34297935E-02,
         0.75522502E-03,-0.47992449E-02, 0.16382027E-02, 0.30208430E-02,
         0.57412515E-03,-0.18284138E-02, 0.53797592E-03,
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

    float v_l_r6_h90s_ep30_q3ex                     =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]                *x51
        +coeff[  3]*x11                
        +coeff[  4]    *x22            
        +coeff[  5]    *x21    *x41    
        +coeff[  6]*x11        *x41    
        +coeff[  7]    *x23            
    ;
    v_l_r6_h90s_ep30_q3ex                     =v_l_r6_h90s_ep30_q3ex                     
        +coeff[  8]    *x23*x31        
        +coeff[  9]            *x41    
        +coeff[ 10]    *x21*x31        
        +coeff[ 11]    *x21        *x51
        +coeff[ 12]                *x52
        +coeff[ 13]    *x21    *x42    
        +coeff[ 14]*x11*x22            
        +coeff[ 15]    *x23    *x41    
        +coeff[ 16]    *x23*x31*x41    
    ;
    v_l_r6_h90s_ep30_q3ex                     =v_l_r6_h90s_ep30_q3ex                     
        +coeff[ 17]        *x31        
        +coeff[ 18]            *x42    
        +coeff[ 19]*x11    *x31        
        +coeff[ 20]    *x21*x31*x41    
        +coeff[ 21]*x11*x22    *x41    
        +coeff[ 22]            *x41*x51
        +coeff[ 23]*x11*x21            
        +coeff[ 24]*x11            *x51
        +coeff[ 25]    *x22    *x41    
    ;
    v_l_r6_h90s_ep30_q3ex                     =v_l_r6_h90s_ep30_q3ex                     
        +coeff[ 26]*x11        *x42    
        +coeff[ 27]    *x22*x31        
        +coeff[ 28]    *x21*x32        
        +coeff[ 29]                *x53
        +coeff[ 30]*x11    *x31*x41    
        +coeff[ 31]*x12*x21            
        +coeff[ 32]    *x24            
        +coeff[ 33]    *x21    *x43    
        +coeff[ 34]*x11*x22*x31        
    ;
    v_l_r6_h90s_ep30_q3ex                     =v_l_r6_h90s_ep30_q3ex                     
        +coeff[ 35]    *x23*x31*x42    
        +coeff[ 36]*x12                
        +coeff[ 37]    *x22        *x51
        +coeff[ 38]    *x21*x31    *x51
        +coeff[ 39]    *x21    *x41*x51
        +coeff[ 40]*x11*x21    *x41    
        +coeff[ 41]    *x21*x31*x42    
        +coeff[ 42]*x11*x23            
        +coeff[ 43]*x12*x21    *x41    
    ;
    v_l_r6_h90s_ep30_q3ex                     =v_l_r6_h90s_ep30_q3ex                     
        +coeff[ 44]    *x23*x32*x41    
        +coeff[ 45]        *x31*x41    
        +coeff[ 46]        *x31    *x51
        +coeff[ 47]            *x42*x51
        +coeff[ 48]*x11*x21*x31        
        +coeff[ 49]*x11        *x41*x51
        +coeff[ 50]*x12        *x41    
        +coeff[ 51]    *x22    *x42    
        +coeff[ 52]*x11        *x43    
    ;
    v_l_r6_h90s_ep30_q3ex                     =v_l_r6_h90s_ep30_q3ex                     
        +coeff[ 53]    *x24    *x41    
        +coeff[ 54]    *x22*x31*x42    
        +coeff[ 55]    *x21    *x44    
        +coeff[ 56]        *x33*x41*x51
        +coeff[ 57]*x11*x22    *x42    
        +coeff[ 58]*x11    *x32    *x52
        ;

    return v_l_r6_h90s_ep30_q3ex                     ;
}
float x_r5p65_ep5                             (float *x,int m){
    int ncoeff=  8;
    float avdat=  0.1135625E+00;
    float xmin[10]={
        -0.99915E-02,-0.54850E-01,-0.99424E-02,-0.30015E-01,-0.49704E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.62053E-01, 0.99768E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  9]={
        -0.24536755E-02, 0.18351061E-07,-0.48044781E-05, 0.67944609E-06,
        -0.10015689E-01,-0.30105026E-01, 0.78706129E-04, 0.25420710E-04,
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
    float x21 = x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x51 = x5;

//                 function

    float v_x_r5p65_ep5                             =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]                *x51
        +coeff[  3]    *x21            
        +coeff[  4]        *x31        
        +coeff[  5]            *x41    
        +coeff[  6]            *x42    
        +coeff[  7]        *x31*x41    
        ;

    return v_x_r5p65_ep5                             ;
}
float t_r5p65_ep5                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.1057822E+00;
    float xmin[10]={
        -0.99915E-02,-0.54850E-01,-0.99424E-02,-0.30015E-01,-0.49704E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.62053E-01, 0.99768E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 51]={
        -0.16734439E-02, 0.12542972E-04, 0.47982121E-04,-0.15780856E-03,
        -0.28254969E-01,-0.13073438E-03, 0.38821911E-03, 0.39539671E-04,
         0.13121557E-03, 0.13652297E-03, 0.49592207E-04,-0.93146700E-04,
        -0.62594314E-04, 0.15590165E-03,-0.11542413E-05, 0.17472259E-03,
        -0.54104348E-04, 0.21652648E-04,-0.35185491E-04, 0.50820934E-04,
         0.67594963E-04,-0.74246869E-04, 0.11653939E-03, 0.21797856E-04,
         0.69890506E-04,-0.16192475E-04, 0.22300468E-04, 0.53830663E-05,
         0.73984556E-05, 0.14974496E-04, 0.49343696E-04,-0.10309419E-03,
        -0.36512720E-04, 0.51159433E-04,-0.75126794E-04, 0.48208316E-04,
         0.34875006E-05,-0.13702978E-04,-0.11877309E-04,-0.19301349E-04,
        -0.26726164E-04, 0.15742857E-04, 0.10007546E-04, 0.15884887E-04,
        -0.12115683E-04,-0.19316214E-04,-0.20036920E-04,-0.39868646E-04,
        -0.18036639E-04,-0.42810439E-04,
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

    float v_t_r5p65_ep5                             =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]            *x41    
        +coeff[  5]                *x51
        +coeff[  6]    *x22    *x41    
        +coeff[  7]    *x22            
    ;
    v_t_r5p65_ep5                             =v_t_r5p65_ep5                             
        +coeff[  8]    *x21    *x41    
        +coeff[  9]    *x22*x31        
        +coeff[ 10]    *x21*x31        
        +coeff[ 11]        *x31*x41    
        +coeff[ 12]    *x23            
        +coeff[ 13]*x11*x21    *x41    
        +coeff[ 14]*x13*x21            
        +coeff[ 15]    *x22    *x42    
        +coeff[ 16]            *x42    
    ;
    v_t_r5p65_ep5                             =v_t_r5p65_ep5                             
        +coeff[ 17]            *x41*x51
        +coeff[ 18]*x11*x22            
        +coeff[ 19]*x11*x21*x31        
        +coeff[ 20]    *x21    *x42    
        +coeff[ 21]    *x23    *x41    
        +coeff[ 22]    *x22*x31*x41    
        +coeff[ 23]*x13*x23            
        +coeff[ 24]*x11*x21            
        +coeff[ 25]        *x32        
    ;
    v_t_r5p65_ep5                             =v_t_r5p65_ep5                             
        +coeff[ 26]*x11        *x41    
        +coeff[ 27]        *x31    *x51
        +coeff[ 28]                *x52
        +coeff[ 29]*x12        *x41    
        +coeff[ 30]    *x21*x31*x41    
        +coeff[ 31]*x11*x23            
        +coeff[ 32]*x11*x22    *x41    
        +coeff[ 33]*x11*x21    *x42    
        +coeff[ 34]*x11*x23    *x41    
    ;
    v_t_r5p65_ep5                             =v_t_r5p65_ep5                             
        +coeff[ 35]*x13*x21*x31*x41    
        +coeff[ 36]*x12                
        +coeff[ 37]        *x31*x42    
        +coeff[ 38]            *x43    
        +coeff[ 39]*x11*x22*x31        
        +coeff[ 40]    *x23*x31        
        +coeff[ 41]    *x22*x32        
        +coeff[ 42]*x11    *x33        
        +coeff[ 43]    *x21    *x43    
    ;
    v_t_r5p65_ep5                             =v_t_r5p65_ep5                             
        +coeff[ 44]    *x22    *x41*x51
        +coeff[ 45]*x11*x23*x31        
        +coeff[ 46]*x11*x22*x31*x41    
        +coeff[ 47]    *x23*x31*x41    
        +coeff[ 48]*x12    *x32*x41    
        +coeff[ 49]    *x23    *x42    
        ;

    return v_t_r5p65_ep5                             ;
}
float y_r5p65_ep5                             (float *x,int m){
    int ncoeff= 11;
    float avdat=  0.2715605E-02;
    float xmin[10]={
        -0.99915E-02,-0.54850E-01,-0.99424E-02,-0.30015E-01,-0.49704E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.62053E-01, 0.99768E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
         0.11991472E-02,-0.13028875E-04, 0.63618638E-01, 0.99837240E-02,
        -0.19210571E-03,-0.63569598E-04,-0.42676556E-05, 0.33191986E-05,
        -0.29401863E-05, 0.80371547E-05, 0.38741532E-05,
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

    float v_y_r5p65_ep5                             =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]*x11                
        +coeff[  4]    *x21    *x41    
        +coeff[  5]    *x21*x31        
        +coeff[  6]        *x31        
        +coeff[  7]    *x22            
    ;
    v_y_r5p65_ep5                             =v_y_r5p65_ep5                             
        +coeff[  8]*x11        *x41    
        +coeff[  9]    *x23            
        +coeff[ 10]*x11*x22            
        ;

    return v_y_r5p65_ep5                             ;
}
float p_r5p65_ep5                             (float *x,int m){
    int ncoeff= 11;
    float avdat=  0.2641549E-02;
    float xmin[10]={
        -0.99915E-02,-0.54850E-01,-0.99424E-02,-0.30015E-01,-0.49704E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.62053E-01, 0.99768E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
         0.86466532E-03, 0.57990693E-01,-0.65153441E-03,-0.13216129E-03,
        -0.15244633E-03,-0.66181885E-04, 0.10056368E-03, 0.24962384E-03,
        -0.71153692E-04,-0.12434259E-03, 0.11363328E-03,
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
    float x42 = x41*x4;

//                 function

    float v_p_r5p65_ep5                             =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]    *x21    *x41    
        +coeff[  3]*x11                
        +coeff[  4]    *x21*x31        
        +coeff[  5]            *x41    
        +coeff[  6]    *x22            
        +coeff[  7]    *x23            
    ;
    v_p_r5p65_ep5                             =v_p_r5p65_ep5                             
        +coeff[  8]*x11        *x41    
        +coeff[  9]    *x21    *x42    
        +coeff[ 10]*x11*x22            
        ;

    return v_p_r5p65_ep5                             ;
}
float l_r6_h90s_ep5                           (float *x,int m){
    int ncoeff=  8;
    float avdat= -0.1263345E-02;
    float xmin[10]={
        -0.99915E-02,-0.54850E-01,-0.99424E-02,-0.30015E-01,-0.49704E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.62053E-01, 0.99768E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  9]={
         0.98211435E-03,-0.22885765E-03, 0.98681112E-03, 0.30431938E-02,
        -0.18575359E-02,-0.41966679E-03,-0.28003626E-05, 0.44188009E-05,
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

    float v_l_r6_h90s_ep5                           =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]    *x22            
        +coeff[  5]            *x42    
        +coeff[  6]        *x31*x41    
        +coeff[  7]    *x22    *x41    
        ;

    return v_l_r6_h90s_ep5                           ;
}
float x_r5p65_ep7                             (float *x,int m){
    int ncoeff= 39;
    float avdat=  0.2719274E+00;
    float xmin[10]={
        -0.99915E-02,-0.50838E-01,-0.99424E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.51603E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 40]={
        -0.53151576E-02,-0.29583378E-02, 0.87299151E-03,-0.10946128E-01,
        -0.60755432E-01, 0.48339418E-02, 0.14901065E-02, 0.33637285E-02,
         0.51688496E-03,-0.11621349E-02,-0.70983142E-03, 0.98095927E-03,
         0.93227194E-03, 0.16052747E-02, 0.14030728E-03, 0.16337070E-03,
         0.16922961E-03, 0.27536249E-03,-0.24015066E-03, 0.79452549E-03,
         0.10470950E-03, 0.87701941E-04, 0.14471448E-03,-0.75617289E-04,
        -0.26390611E-03, 0.22981461E-03,-0.18192448E-03, 0.43645164E-03,
        -0.26682066E-04,-0.10515780E-03, 0.22723494E-04, 0.36576435E-04,
         0.81899976E-04,-0.54897799E-04,-0.97937911E-04, 0.14357502E-03,
        -0.10626337E-03,-0.21947408E-03, 0.30950306E-03,
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

    float v_x_r5p65_ep7                             =avdat
        +coeff[  0]                    
        +coeff[  1]                *x51
        +coeff[  2]    *x21            
        +coeff[  3]        *x31        
        +coeff[  4]            *x41    
        +coeff[  5]    *x22            
        +coeff[  6]*x11*x21            
        +coeff[  7]    *x22    *x41    
    ;
    v_x_r5p65_ep7                             =v_x_r5p65_ep7                             
        +coeff[  8]    *x21    *x41    
        +coeff[  9]            *x42    
        +coeff[ 10]        *x31*x41    
        +coeff[ 11]*x11*x21    *x41    
        +coeff[ 12]    *x22*x31        
        +coeff[ 13]    *x22    *x42    
        +coeff[ 14]*x11                
        +coeff[ 15]                *x52
        +coeff[ 16]            *x41*x51
    ;
    v_x_r5p65_ep7                             =v_x_r5p65_ep7                             
        +coeff[ 17]*x11*x21*x31        
        +coeff[ 18]    *x22        *x51
        +coeff[ 19]    *x22*x31*x41    
        +coeff[ 20]*x12                
        +coeff[ 21]*x11        *x41    
        +coeff[ 22]    *x21*x31        
        +coeff[ 23]        *x32        
        +coeff[ 24]    *x23            
        +coeff[ 25]    *x21    *x42    
    ;
    v_x_r5p65_ep7                             =v_x_r5p65_ep7                             
        +coeff[ 26]            *x43    
        +coeff[ 27]*x11*x21    *x42    
        +coeff[ 28]    *x23*x31*x41    
        +coeff[ 29]*x13*x21*x31*x41    
        +coeff[ 30]*x11    *x31        
        +coeff[ 31]        *x31    *x51
        +coeff[ 32]*x12        *x41    
        +coeff[ 33]*x11*x21        *x51
        +coeff[ 34]*x11*x22            
    ;
    v_x_r5p65_ep7                             =v_x_r5p65_ep7                             
        +coeff[ 35]    *x21*x31*x41    
        +coeff[ 36]        *x31*x42    
        +coeff[ 37]*x11*x23            
        +coeff[ 38]*x11*x21*x31*x41    
        ;

    return v_x_r5p65_ep7                             ;
}
float t_r5p65_ep7                             (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.2256751E+00;
    float xmin[10]={
        -0.99915E-02,-0.50838E-01,-0.99424E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.51603E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 51]={
        -0.42982278E-02, 0.16195700E-02,-0.12752294E-02,-0.33849847E-01,
        -0.61888625E-02, 0.95703145E-02,-0.24014779E-02, 0.47815945E-02,
         0.63466027E-05, 0.26661551E-02, 0.74699306E-03,-0.13133852E-02,
         0.15401125E-02, 0.27507450E-03, 0.26960319E-03, 0.34346920E-03,
         0.13291677E-02, 0.32707055E-02,-0.16155632E-03, 0.17283764E-03,
         0.25368668E-03, 0.13077777E-03,-0.58942568E-03, 0.55390177E-03,
         0.41384099E-03,-0.50580339E-03, 0.16510012E-02, 0.79101650E-03,
        -0.18595504E-03,-0.19558008E-03, 0.22059138E-03,-0.43710950E-03,
        -0.39999239E-03, 0.42804366E-03,-0.58975318E-04, 0.36266592E-04,
        -0.58460679E-04, 0.50995139E-04, 0.97548742E-04, 0.59575807E-04,
        -0.91653579E-04,-0.19504914E-03,-0.95806041E-04, 0.89474415E-04,
         0.21283522E-03,-0.35266770E-03,-0.29971203E-03, 0.25088555E-03,
        -0.18752608E-03,-0.56672317E-04,
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

    float v_t_r5p65_ep7                             =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]        *x31        
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x22            
        +coeff[  6]            *x42    
        +coeff[  7]    *x22    *x41    
    ;
    v_t_r5p65_ep7                             =v_t_r5p65_ep7                             
        +coeff[  8]*x13*x21            
        +coeff[  9]*x11*x21            
        +coeff[ 10]    *x21    *x41    
        +coeff[ 11]        *x31*x41    
        +coeff[ 12]*x11*x21    *x41    
        +coeff[ 13]*x11                
        +coeff[ 14]            *x41*x51
        +coeff[ 15]                *x52
        +coeff[ 16]    *x22*x31        
    ;
    v_t_r5p65_ep7                             =v_t_r5p65_ep7                             
        +coeff[ 17]    *x22    *x42    
        +coeff[ 18]    *x22*x33*x41    
        +coeff[ 19]*x12                
        +coeff[ 20]    *x21*x31        
        +coeff[ 21]*x11        *x41    
        +coeff[ 22]    *x23            
        +coeff[ 23]*x11*x21*x31        
        +coeff[ 24]    *x21    *x42    
        +coeff[ 25]    *x22        *x51
    ;
    v_t_r5p65_ep7                             =v_t_r5p65_ep7                             
        +coeff[ 26]    *x22*x31*x41    
        +coeff[ 27]*x11*x21    *x42    
        +coeff[ 28]        *x32        
        +coeff[ 29]*x11*x22            
        +coeff[ 30]    *x21*x31*x41    
        +coeff[ 31]            *x43    
        +coeff[ 32]*x11*x23            
        +coeff[ 33]*x11*x21*x31*x41    
        +coeff[ 34]        *x33*x42    
    ;
    v_t_r5p65_ep7                             =v_t_r5p65_ep7                             
        +coeff[ 35]*x11    *x31        
        +coeff[ 36]    *x21        *x51
        +coeff[ 37]        *x31    *x51
        +coeff[ 38]*x12        *x41    
        +coeff[ 39]*x11    *x31*x41    
        +coeff[ 40]        *x32*x41    
        +coeff[ 41]        *x31*x42    
        +coeff[ 42]*x11*x21        *x51
        +coeff[ 43]            *x42*x51
    ;
    v_t_r5p65_ep7                             =v_t_r5p65_ep7                             
        +coeff[ 44]    *x22*x32        
        +coeff[ 45]*x11*x23*x31        
        +coeff[ 46]*x11*x23    *x41    
        +coeff[ 47]    *x22    *x43    
        +coeff[ 48]    *x23*x33        
        +coeff[ 49]*x13                
        ;

    return v_t_r5p65_ep7                             ;
}
float y_r5p65_ep7                             (float *x,int m){
    int ncoeff= 81;
    float avdat= -0.7094019E-03;
    float xmin[10]={
        -0.99915E-02,-0.50838E-01,-0.99424E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.51603E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 82]={
         0.10550052E-02,-0.43983842E-03, 0.10051626E+00, 0.92178043E-02,
        -0.55397190E-02,-0.14493772E-02,-0.76256011E-03,-0.19930850E-02,
        -0.11011622E-02, 0.20652658E-02, 0.90169028E-03, 0.56010758E-03,
         0.19431761E-02,-0.11849164E-03,-0.15953172E-03, 0.16763224E-03,
        -0.21112851E-03, 0.19933526E-03, 0.52004645E-03,-0.30077543E-03,
         0.86643989E-03,-0.86843072E-04,-0.14141093E-03, 0.15749881E-03,
        -0.16262321E-03, 0.24187543E-03,-0.56532217E-03,-0.47772360E-03,
         0.14085000E-03, 0.54691691E-03, 0.18763556E-04, 0.35350629E-04,
         0.12346698E-03, 0.59142458E-04,-0.11210174E-03, 0.22789498E-03,
         0.11370079E-03, 0.42674289E-03, 0.12007695E-04,-0.38553986E-04,
         0.33542659E-04, 0.24113435E-04, 0.10940729E-03,-0.76555669E-04,
         0.57016863E-04, 0.68916408E-04,-0.93664283E-04, 0.21548725E-03,
         0.33137939E-03, 0.12802970E-03,-0.10047930E-04,-0.44075725E-04,
         0.11148157E-04,-0.19348197E-04,-0.17759394E-04, 0.16411645E-04,
         0.38833528E-04,-0.47258767E-04,-0.94754665E-04, 0.37526774E-04,
         0.16152735E-04,-0.52817602E-04, 0.35329016E-04,-0.38771163E-04,
        -0.95612741E-04,-0.59776667E-04,-0.84091385E-04,-0.48345311E-04,
        -0.20655460E-03, 0.65964421E-04, 0.56066465E-05, 0.60329285E-05,
        -0.58360229E-05,-0.15585383E-04, 0.22987995E-04, 0.17878914E-04,
         0.25678773E-04, 0.60226514E-04,-0.23796447E-04,-0.71225215E-04,
         0.37962811E-04,
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

    float v_y_r5p65_ep7                             =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]*x11                
        +coeff[  4]    *x21    *x41    
        +coeff[  5]    *x21*x31        
        +coeff[  6]*x11        *x41    
        +coeff[  7]    *x21    *x42    
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[  8]    *x21*x31*x41    
        +coeff[  9]    *x23            
        +coeff[ 10]*x11*x22            
        +coeff[ 11]    *x22            
        +coeff[ 12]    *x23    *x41    
        +coeff[ 13]        *x31        
        +coeff[ 14]            *x42    
        +coeff[ 15]*x11*x21            
        +coeff[ 16]*x11    *x31        
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 17]    *x21        *x51
        +coeff[ 18]    *x22    *x41    
        +coeff[ 19]*x11        *x42    
        +coeff[ 20]*x11*x22    *x41    
        +coeff[ 21]        *x31*x41    
        +coeff[ 22]    *x21*x32        
        +coeff[ 23]*x11*x21    *x41    
        +coeff[ 24]*x11    *x31*x41    
        +coeff[ 25]    *x21    *x41*x51
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 26]    *x21    *x43    
        +coeff[ 27]    *x21*x31*x42    
        +coeff[ 28]*x12*x21            
        +coeff[ 29]    *x23*x31        
        +coeff[ 30]                *x51
        +coeff[ 31]*x11            *x51
        +coeff[ 32]    *x22*x31        
        +coeff[ 33]    *x21*x31    *x51
        +coeff[ 34]    *x21*x32*x41    
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 35]*x11*x22*x31        
        +coeff[ 36]*x12*x21    *x41    
        +coeff[ 37]    *x23    *x42    
        +coeff[ 38]            *x41*x51
        +coeff[ 39]            *x43    
        +coeff[ 40]*x11*x21*x31        
        +coeff[ 41]*x11        *x41*x51
        +coeff[ 42]    *x22    *x42    
        +coeff[ 43]*x11        *x43    
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 44]    *x22*x31*x41    
        +coeff[ 45]    *x21    *x42*x51
        +coeff[ 46]    *x23        *x51
        +coeff[ 47]    *x23*x31*x41    
        +coeff[ 48]*x11*x22    *x42    
        +coeff[ 49]*x11*x22*x31*x41    
        +coeff[ 50]        *x32        
        +coeff[ 51]        *x31*x42    
        +coeff[ 52]*x12                
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 53]*x11    *x32        
        +coeff[ 54]    *x22        *x51
        +coeff[ 55]*x12        *x41    
        +coeff[ 56]*x11*x21    *x42    
        +coeff[ 57]*x11    *x31*x42    
        +coeff[ 58]    *x24            
        +coeff[ 59]*x11*x21*x31*x41    
        +coeff[ 60]    *x21*x31*x41*x51
        +coeff[ 61]*x11*x23            
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 62]*x12*x21*x31        
        +coeff[ 63]*x11*x22        *x51
        +coeff[ 64]    *x24    *x41    
        +coeff[ 65]*x11*x23    *x41    
        +coeff[ 66]    *x23    *x41*x51
        +coeff[ 67]    *x23*x33        
        +coeff[ 68]*x11*x24    *x42    
        +coeff[ 69]*x12*x23*x31*x41    
        +coeff[ 70]*x12    *x31        
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 71]*x11    *x31    *x51
        +coeff[ 72]    *x21        *x52
        +coeff[ 73]*x11    *x32*x41    
        +coeff[ 74]        *x31*x44    
        +coeff[ 75]    *x21    *x43*x51
        +coeff[ 76]*x11*x21*x31*x42    
        +coeff[ 77]    *x21*x31*x44    
        +coeff[ 78]*x12*x22    *x41    
        +coeff[ 79]*x11*x24    *x41    
    ;
    v_y_r5p65_ep7                             =v_y_r5p65_ep7                             
        +coeff[ 80]    *x21*x33*x41*x51
        ;

    return v_y_r5p65_ep7                             ;
}
float p_r5p65_ep7                             (float *x,int m){
    int ncoeff= 36;
    float avdat= -0.5841436E-03;
    float xmin[10]={
        -0.99915E-02,-0.50838E-01,-0.99424E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.51603E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 37]={
         0.29349717E-03, 0.43666326E-01,-0.12094827E-02,-0.89465296E-02,
         0.88374293E-03,-0.71152009E-03, 0.33059521E-02,-0.12331933E-02,
        -0.26431298E-02, 0.71922690E-03,-0.22172909E-02,-0.14181264E-02,
         0.13516437E-02, 0.38734970E-02,-0.17311731E-03, 0.26183203E-03,
        -0.32705395E-03, 0.10512479E-02,-0.11254995E-02, 0.14162805E-02,
        -0.17145592E-03, 0.19034174E-03, 0.32292103E-03, 0.12241564E-03,
        -0.35274526E-03,-0.74775255E-03, 0.17427781E-03, 0.37888947E-03,
        -0.37710823E-03, 0.23539727E-04,-0.85593063E-04, 0.16062906E-03,
         0.85885244E-04,-0.15353161E-03,-0.16094859E-03,-0.43896108E-03,
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

//                 function

    float v_p_r5p65_ep7                             =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]*x11                
        +coeff[  3]    *x21    *x41    
        +coeff[  4]    *x23*x31        
        +coeff[  5]            *x41    
        +coeff[  6]    *x23            
        +coeff[  7]*x11        *x41    
    ;
    v_p_r5p65_ep7                             =v_p_r5p65_ep7                             
        +coeff[  8]    *x21    *x42    
        +coeff[  9]    *x22            
        +coeff[ 10]    *x21*x31        
        +coeff[ 11]    *x21*x31*x41    
        +coeff[ 12]*x11*x22            
        +coeff[ 13]    *x23    *x41    
        +coeff[ 14]        *x31        
        +coeff[ 15]    *x21        *x51
        +coeff[ 16]*x11    *x31        
    ;
    v_p_r5p65_ep7                             =v_p_r5p65_ep7                             
        +coeff[ 17]    *x22    *x41    
        +coeff[ 18]    *x21    *x43    
        +coeff[ 19]*x11*x22    *x41    
        +coeff[ 20]            *x42    
        +coeff[ 21]*x11*x21            
        +coeff[ 22]    *x21    *x41*x51
        +coeff[ 23]*x11*x21    *x41    
        +coeff[ 24]*x11        *x42    
        +coeff[ 25]    *x21*x31*x42    
    ;
    v_p_r5p65_ep7                             =v_p_r5p65_ep7                             
        +coeff[ 26]*x12*x21            
        +coeff[ 27]*x11*x22*x31        
        +coeff[ 28]    *x23*x32        
        +coeff[ 29]                *x51
        +coeff[ 30]        *x31*x41    
        +coeff[ 31]    *x22*x31        
        +coeff[ 32]    *x21*x31    *x51
        +coeff[ 33]*x11    *x31*x41    
        +coeff[ 34]    *x21*x32*x41    
    ;
    v_p_r5p65_ep7                             =v_p_r5p65_ep7                             
        +coeff[ 35]    *x24    *x41    
        ;

    return v_p_r5p65_ep7                             ;
}
float l_r6_h90s_ep7                           (float *x,int m){
    int ncoeff= 45;
    float avdat= -0.2448206E-02;
    float xmin[10]={
        -0.99915E-02,-0.50838E-01,-0.99424E-02,-0.30015E-01,-0.49830E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99868E-02, 0.51603E-01, 0.99474E-02, 0.24967E-01, 0.49887E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 46]={
         0.17219717E-02,-0.17519855E-03, 0.79930900E-02, 0.53821463E-03,
        -0.34108991E-02,-0.66422898E-03,-0.10708183E-03, 0.10583933E-05,
        -0.10763309E-05, 0.63795540E-07, 0.32953230E-05, 0.11468471E-02,
        -0.24206698E-03, 0.11178716E-03,-0.12553712E-03,-0.24252562E-03,
        -0.24333145E-04,-0.58105732E-04,-0.34296470E-04,-0.10877138E-03,
        -0.22279592E-04,-0.18999714E-04, 0.60052589E-04,-0.39350598E-04,
        -0.22621790E-03, 0.46622801E-04, 0.13249766E-04, 0.87152957E-05,
        -0.88243187E-05,-0.10257411E-04,-0.25817741E-04, 0.16491120E-04,
        -0.12809383E-03,-0.42324045E-04, 0.33490887E-05, 0.11921316E-03,
        -0.47537837E-04, 0.22439439E-04,-0.18653393E-04, 0.71097002E-05,
         0.13120192E-04,-0.10008891E-04,-0.18318073E-04, 0.16678576E-04,
         0.28742736E-04,
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

    float v_l_r6_h90s_ep7                           =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]            *x41    
        +coeff[  3]                *x51
        +coeff[  4]    *x22            
        +coeff[  5]            *x42    
        +coeff[  6]    *x22*x31        
        +coeff[  7]        *x33        
    ;
    v_l_r6_h90s_ep7                           =v_l_r6_h90s_ep7                           
        +coeff[  8]        *x31    *x52
        +coeff[  9]    *x22*x33        
        +coeff[ 10]        *x33    *x52
        +coeff[ 11]        *x31        
        +coeff[ 12]*x11*x21            
        +coeff[ 13]        *x31*x41    
        +coeff[ 14]            *x41*x51
        +coeff[ 15]    *x22    *x41    
        +coeff[ 16]*x11                
    ;
    v_l_r6_h90s_ep7                           =v_l_r6_h90s_ep7                           
        +coeff[ 17]    *x21    *x41    
        +coeff[ 18]                *x52
        +coeff[ 19]*x11*x21    *x41    
        +coeff[ 20]    *x21*x31        
        +coeff[ 21]*x12                
        +coeff[ 22]    *x22        *x51
        +coeff[ 23]*x11*x21*x31        
        +coeff[ 24]    *x22    *x42    
        +coeff[ 25]    *x24*x31*x41    
    ;
    v_l_r6_h90s_ep7                           =v_l_r6_h90s_ep7                           
        +coeff[ 26]        *x32        
        +coeff[ 27]    *x21        *x51
        +coeff[ 28]        *x31    *x51
        +coeff[ 29]*x11        *x41    
        +coeff[ 30]    *x21    *x42    
        +coeff[ 31]*x11*x21        *x51
        +coeff[ 32]    *x22*x31*x41    
        +coeff[ 33]*x11*x21    *x42    
        +coeff[ 34]    *x23*x31*x41    
    ;
    v_l_r6_h90s_ep7                           =v_l_r6_h90s_ep7                           
        +coeff[ 35]    *x24    *x42    
        +coeff[ 36]*x11*x23*x31*x41    
        +coeff[ 37]    *x23            
        +coeff[ 38]    *x21*x31*x41    
        +coeff[ 39]            *x41*x52
        +coeff[ 40]*x11*x22            
        +coeff[ 41]*x12        *x41    
        +coeff[ 42]        *x31*x43    
        +coeff[ 43]    *x22    *x41*x51
    ;
    v_l_r6_h90s_ep7                           =v_l_r6_h90s_ep7                           
        +coeff[ 44]*x11*x23            
        ;

    return v_l_r6_h90s_ep7                           ;
}
float x_r5p65_fp                              (float *x,int m){
    int ncoeff= 45;
    float avdat=  0.2802399E-01;
    float xmin[10]={
        -0.99483E-02,-0.47003E-01,-0.99484E-02,-0.29097E-01,-0.46754E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99850E-02, 0.50561E-01, 0.99954E-02, 0.24794E-01, 0.49963E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 46]={
        -0.59442525E-02, 0.64600217E+00,-0.51850930E-01, 0.46876807E-01,
         0.53091255E-04,-0.23759795E-01, 0.16971052E-01,-0.28334511E-01,
         0.97197862E-02,-0.78712804E-02,-0.12761607E-01,-0.39359503E-02,
         0.56110858E-02,-0.24773827E-01,-0.85343920E-04, 0.63407626E-02,
        -0.66271576E-04,-0.37711903E-02,-0.27226191E-02,-0.34600163E-02,
         0.32030139E-02,-0.13300247E-01,-0.60335333E-02, 0.20369288E-01,
         0.81687337E-02,-0.16955660E-03,-0.74053626E-03,-0.15203432E-02,
        -0.13773559E-02,-0.85394556E-03, 0.34948757E-02,-0.18145323E-02,
         0.51325969E-02, 0.36937650E-02, 0.65314989E-02, 0.50207423E-02,
        -0.14698431E-02, 0.69349015E-03,-0.53871924E-03, 0.10336634E-02,
         0.11542466E-02, 0.17909295E-02, 0.26594901E-02, 0.97245368E-03,
        -0.39679562E-02,
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

    float v_x_r5p65_fp                              =avdat
        +coeff[  0]                    
        +coeff[  1]                *x51
        +coeff[  2]                *x52
        +coeff[  3]    *x21        *x51
        +coeff[  4]*x13                
        +coeff[  5]*x11                
        +coeff[  6]    *x22            
        +coeff[  7]    *x21    *x41    
    ;
    v_x_r5p65_fp                              =v_x_r5p65_fp                              
        +coeff[  8]    *x23            
        +coeff[  9]    *x21*x31        
        +coeff[ 10]            *x42    
        +coeff[ 11]*x11        *x41    
        +coeff[ 12]                *x53
        +coeff[ 13]    *x21    *x42    
        +coeff[ 14]*x12*x21    *x41*x51
        +coeff[ 15]    *x21            
        +coeff[ 16]*x11*x21            
    ;
    v_x_r5p65_fp                              =v_x_r5p65_fp                              
        +coeff[ 17]            *x41*x51
        +coeff[ 18]        *x31*x41    
        +coeff[ 19]    *x21        *x52
        +coeff[ 20]*x11*x22            
        +coeff[ 21]    *x21    *x41*x51
        +coeff[ 22]    *x21*x31*x41    
        +coeff[ 23]    *x23    *x41    
        +coeff[ 24]    *x22    *x42    
        +coeff[ 25]*x11*x22    *x43    
    ;
    v_x_r5p65_fp                              =v_x_r5p65_fp                              
        +coeff[ 26]        *x31        
        +coeff[ 27]            *x41    
        +coeff[ 28]*x11            *x51
        +coeff[ 29]*x11    *x31        
        +coeff[ 30]    *x22        *x51
        +coeff[ 31]    *x21*x31    *x51
        +coeff[ 32]    *x22    *x41    
        +coeff[ 33]*x11*x23            
        +coeff[ 34]*x11*x22    *x41    
    ;
    v_x_r5p65_fp                              =v_x_r5p65_fp                              
        +coeff[ 35]    *x23        *x51
        +coeff[ 36]    *x22*x31    *x51
        +coeff[ 37]*x12*x21            
        +coeff[ 38]*x11        *x41*x51
        +coeff[ 39]*x11*x21    *x41    
        +coeff[ 40]            *x42*x51
        +coeff[ 41]*x11*x22        *x51
        +coeff[ 42]    *x23*x31        
        +coeff[ 43]    *x22*x33        
    ;
    v_x_r5p65_fp                              =v_x_r5p65_fp                              
        +coeff[ 44]    *x22    *x43    
        ;

    return v_x_r5p65_fp                              ;
}
float t_r5p65_fp                              (float *x,int m){
    int ncoeff= 50;
    float avdat=  0.4979818E-02;
    float xmin[10]={
        -0.99483E-02,-0.47003E-01,-0.99484E-02,-0.29097E-01,-0.46754E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99850E-02, 0.50561E-01, 0.99954E-02, 0.24794E-01, 0.49963E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 51]={
        -0.16669293E-02,-0.22904230E-02,-0.19902101E-01, 0.40335843E-03,
         0.11170875E+00, 0.57479902E-02,-0.10345817E-01, 0.95167616E-03,
        -0.18813451E-02, 0.97721757E-03,-0.20227374E-02,-0.21169970E-02,
        -0.97489933E-03, 0.22652204E-04,-0.38467476E-03,-0.65281417E-03,
         0.10997613E-02, 0.16169017E-02, 0.37160124E-04, 0.16422603E-03,
         0.14530662E-03,-0.21370586E-03, 0.21952180E-03, 0.63891109E-03,
        -0.11477386E-03,-0.19712071E-03, 0.16808813E-03, 0.43605233E-03,
         0.20400625E-03, 0.23028678E-03,-0.23141128E-03, 0.74053043E-03,
         0.43068672E-03,-0.67794055E-03,-0.27917889E-04, 0.60764243E-04,
         0.41181193E-03, 0.26699429E-03, 0.36993268E-03,-0.26885435E-03,
        -0.55757409E-03,-0.29281634E-03, 0.99816549E-04,-0.16500732E-03,
         0.82177622E-03, 0.37268762E-03, 0.56361849E-03,-0.23480954E-03,
         0.24652818E-04,-0.95395131E-04,
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

    float v_t_r5p65_fp                              =avdat
        +coeff[  0]                    
        +coeff[  1]*x11                
        +coeff[  2]    *x21            
        +coeff[  3]            *x41    
        +coeff[  4]                *x51
        +coeff[  5]    *x21        *x51
        +coeff[  6]                *x52
        +coeff[  7]    *x22            
    ;
    v_t_r5p65_fp                              =v_t_r5p65_fp                              
        +coeff[  8]            *x42    
        +coeff[  9]    *x21    *x41    
        +coeff[ 10]    *x21    *x42    
        +coeff[ 11]    *x21    *x41*x51
        +coeff[ 12]    *x21        *x52
        +coeff[ 13]*x11*x21            
        +coeff[ 14]        *x31*x41    
        +coeff[ 15]            *x41*x51
        +coeff[ 16]                *x53
    ;
    v_t_r5p65_fp                              =v_t_r5p65_fp                              
        +coeff[ 17]    *x22    *x42    
        +coeff[ 18]        *x31        
        +coeff[ 19]    *x21*x31        
        +coeff[ 20]*x11        *x41    
        +coeff[ 21]*x11            *x51
        +coeff[ 22]            *x42*x51
        +coeff[ 23]*x11*x23            
        +coeff[ 24]        *x31    *x51
        +coeff[ 25]*x11*x22            
    ;
    v_t_r5p65_fp                              =v_t_r5p65_fp                              
        +coeff[ 26]    *x22*x31        
        +coeff[ 27]    *x22    *x41    
        +coeff[ 28]*x11        *x42    
        +coeff[ 29]    *x22        *x51
        +coeff[ 30]    *x21*x31    *x51
        +coeff[ 31]    *x23    *x41    
        +coeff[ 32]    *x23        *x51
        +coeff[ 33]    *x22    *x43    
        +coeff[ 34]*x12                
    ;
    v_t_r5p65_fp                              =v_t_r5p65_fp                              
        +coeff[ 35]*x11            *x52
        +coeff[ 36]    *x21    *x43    
        +coeff[ 37]*x11*x22        *x51
        +coeff[ 38]    *x22    *x41*x51
        +coeff[ 39]    *x21*x31*x41*x51
        +coeff[ 40]    *x21    *x42*x51
        +coeff[ 41]    *x22        *x52
        +coeff[ 42]    *x21        *x53
        +coeff[ 43]    *x23    *x42    
    ;
    v_t_r5p65_fp                              =v_t_r5p65_fp                              
        +coeff[ 44]    *x23    *x41*x51
        +coeff[ 45]    *x22    *x41*x52
        +coeff[ 46]    *x21    *x42*x52
        +coeff[ 47]*x12        *x42*x52
        +coeff[ 48]*x11    *x31        
        +coeff[ 49]    *x23            
        ;

    return v_t_r5p65_fp                              ;
}
float y_r5p65_fp                              (float *x,int m){
    int ncoeff= 97;
    float avdat=  0.5639369E-02;
    float xmin[10]={
        -0.99483E-02,-0.47003E-01,-0.99484E-02,-0.29097E-01,-0.46754E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99850E-02, 0.50561E-01, 0.99954E-02, 0.24794E-01, 0.49963E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 98]={
        -0.19296058E-02,-0.51570166E-01, 0.27344122E-02,-0.93162218E-02,
        -0.50216755E-02, 0.11279478E-01, 0.15175323E-01,-0.34888415E-02,
         0.42902725E-02, 0.15366526E-02,-0.31230224E-02, 0.93729245E-02,
        -0.65009743E-02, 0.60842647E-02, 0.55175251E-03,-0.24327685E-02,
        -0.64427039E-03, 0.15614829E-02,-0.41936492E-02, 0.27883588E-02,
         0.28632793E-02,-0.31744582E-02, 0.39615906E-02, 0.15259538E-02,
        -0.11263569E-02, 0.18399525E-02,-0.84832252E-03,-0.16185860E-02,
         0.18349258E-02, 0.19880584E-02,-0.25479298E-02, 0.25665967E-02,
        -0.42478740E-02, 0.58163196E-03,-0.27711299E-03, 0.37683902E-03,
         0.33312943E-03, 0.66483242E-03,-0.10914817E-02, 0.57368702E-03,
         0.44596259E-03,-0.58277731E-03, 0.14100304E-02,-0.62479503E-02,
        -0.15840332E-02, 0.25529452E-03, 0.17511895E-04, 0.10268333E-03,
        -0.64531423E-03, 0.39202051E-03,-0.11795257E-02, 0.71901706E-03,
         0.22619283E-03, 0.70423202E-03,-0.96763868E-03,-0.63421851E-03,
         0.11246884E-02, 0.31315238E-03, 0.75200910E-03, 0.11773120E-02,
         0.20276788E-02,-0.10136172E-04,-0.92420849E-03,-0.33276161E-03,
        -0.33472193E-03, 0.12474983E-03, 0.16295341E-03,-0.11427042E-03,
        -0.52064587E-03,-0.29842684E-03,-0.30736998E-03, 0.29963179E-03,
         0.60248951E-03, 0.14182120E-02,-0.47595846E-03,-0.46536210E-03,
         0.10116283E-02, 0.53102449E-04, 0.86998698E-04,-0.44701871E-03,
         0.11737989E-02,-0.35708610E-03,-0.57455822E-03,-0.99892577E-03,
        -0.34958200E-03, 0.88324212E-03,-0.10631260E-03,-0.45932076E-03,
         0.39041515E-04, 0.17719573E-03,-0.26224033E-03,-0.24606462E-03,
         0.19599324E-03, 0.10254631E-03, 0.10211332E-02,-0.12391292E-03,
         0.29535985E-03,
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
    float x43 = x42*x4;
    float x44 = x43*x4;
    float x51 = x5;
    float x52 = x51*x5;
    float x53 = x52*x5;

//                 function

    float v_y_r5p65_fp                              =avdat
        +coeff[  0]                    
        +coeff[  1]            *x41    
        +coeff[  2]    *x21            
        +coeff[  3]                *x51
        +coeff[  4]            *x42    
        +coeff[  5]    *x21    *x41    
        +coeff[  6]    *x22            
        +coeff[  7]            *x41*x51
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[  8]*x11*x21            
        +coeff[  9]    *x21        *x51
        +coeff[ 10]        *x31    *x51
        +coeff[ 11]    *x22    *x41    
        +coeff[ 12]    *x21    *x41*x51
        +coeff[ 13]            *x41*x52
        +coeff[ 14]*x11                
        +coeff[ 15]        *x31*x41    
        +coeff[ 16]*x11        *x41    
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 17]                *x52
        +coeff[ 18]    *x23            
        +coeff[ 19]    *x22*x31        
        +coeff[ 20]*x11*x21    *x41    
        +coeff[ 21]    *x22        *x51
        +coeff[ 22]    *x22    *x42    
        +coeff[ 23]    *x22    *x41*x53
        +coeff[ 24]        *x31        
        +coeff[ 25]            *x42*x51
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 26]    *x21    *x43    
        +coeff[ 27]    *x21        *x52
        +coeff[ 28]        *x31    *x52
        +coeff[ 29]    *x22*x31*x41    
        +coeff[ 30]    *x24            
        +coeff[ 31]    *x23        *x51
        +coeff[ 32]            *x41*x53
        +coeff[ 33]    *x21*x31        
        +coeff[ 34]        *x32        
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 35]            *x43    
        +coeff[ 36]*x12                
        +coeff[ 37]        *x31*x41*x51
        +coeff[ 38]*x11*x22            
        +coeff[ 39]*x11*x21*x31        
        +coeff[ 40]*x11        *x41*x51
        +coeff[ 41]*x11*x21        *x51
        +coeff[ 42]*x11*x21    *x42    
        +coeff[ 43]    *x22    *x41*x51
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 44]*x11*x23            
        +coeff[ 45]*x11*x22        *x51
        +coeff[ 46]*x13            *x51
        +coeff[ 47]*x11    *x31        
        +coeff[ 48]    *x21    *x42    
        +coeff[ 49]    *x21*x31*x41    
        +coeff[ 50]            *x43*x51
        +coeff[ 51]    *x21    *x42*x51
        +coeff[ 52]*x11            *x52
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 53]*x11*x21*x31*x41    
        +coeff[ 54]            *x42*x52
        +coeff[ 55]    *x22*x31    *x51
        +coeff[ 56]    *x21    *x41*x52
        +coeff[ 57]    *x21*x31    *x52
        +coeff[ 58]    *x21        *x53
        +coeff[ 59]    *x24        *x51
        +coeff[ 60]    *x22    *x41*x52
        +coeff[ 61]    *x23*x31*x42    
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 62]    *x21    *x41*x53
        +coeff[ 63]*x11            *x51
        +coeff[ 64]*x11        *x42    
        +coeff[ 65]    *x21*x32        
        +coeff[ 66]*x12        *x41    
        +coeff[ 67]*x11    *x31    *x51
        +coeff[ 68]    *x23*x31        
        +coeff[ 69]                *x53
        +coeff[ 70]*x12*x22            
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 71]*x11        *x41*x52
        +coeff[ 72]*x11*x22    *x42    
        +coeff[ 73]    *x22    *x42*x51
        +coeff[ 74]    *x21*x31*x42*x51
        +coeff[ 75]        *x31    *x53
        +coeff[ 76]    *x23    *x41*x51
        +coeff[ 77]    *x22    *x44    
        +coeff[ 78]    *x21    *x42*x52
        +coeff[ 79]            *x42*x53
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 80]    *x22        *x53
        +coeff[ 81]    *x23*x32    *x51
        +coeff[ 82]    *x24        *x52
        +coeff[ 83]*x11*x23    *x43    
        +coeff[ 84]*x12        *x41*x53
        +coeff[ 85]*x12*x21    *x42*x52
        +coeff[ 86]        *x32*x41    
        +coeff[ 87]    *x21*x31*x42    
        +coeff[ 88]*x12    *x31        
    ;
    v_y_r5p65_fp                              =v_y_r5p65_fp                              
        +coeff[ 89]    *x22*x32        
        +coeff[ 90]    *x21*x31*x41*x51
        +coeff[ 91]        *x31*x44    
        +coeff[ 92]*x11        *x42*x51
        +coeff[ 93]*x11*x21*x32        
        +coeff[ 94]    *x22    *x43    
        +coeff[ 95]        *x31*x41*x52
        +coeff[ 96]    *x22*x31*x42    
        ;

    return v_y_r5p65_fp                              ;
}
float p_r5p65_fp                              (float *x,int m){
    int ncoeff= 88;
    float avdat=  0.2297608E-02;
    float xmin[10]={
        -0.99483E-02,-0.47003E-01,-0.99484E-02,-0.29097E-01,-0.46754E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99850E-02, 0.50561E-01, 0.99954E-02, 0.24794E-01, 0.49963E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 89]={
        -0.14959753E-02, 0.46298648E-02,-0.26987864E-01,-0.85318657E-02,
         0.53787604E-03, 0.13722853E-01, 0.18886551E-02, 0.22069480E-01,
        -0.53961580E-02, 0.35984211E-02,-0.18753551E-01, 0.41305036E-02,
        -0.63949302E-02,-0.17759877E-02,-0.13606615E-02, 0.16882235E-01,
         0.37873187E-02, 0.20581431E-03,-0.32449287E-04,-0.76381853E-04,
        -0.22885614E-03,-0.31290925E-04,-0.24671957E-02,-0.23697692E-02,
        -0.34493215E-02,-0.45969114E-02, 0.16511740E-02, 0.48802304E-03,
        -0.23538223E-03,-0.41542589E-04, 0.16319572E-02, 0.29153801E-02,
        -0.61454240E-03,-0.16821142E-03, 0.26147766E-03, 0.32382293E-03,
         0.56857168E-03, 0.55961078E-03, 0.26093477E-02, 0.57431817E-03,
        -0.15089066E-04, 0.15181226E-03,-0.51479624E-05, 0.11337604E-02,
        -0.11280626E-02,-0.56517514E-03,-0.22123524E-02, 0.34761059E-03,
        -0.17754300E-02, 0.70493633E-03, 0.16550403E-03,-0.10261998E-02,
        -0.28096390E-03, 0.57421614E-04, 0.20149763E-03,-0.23439876E-03,
        -0.60254952E-03,-0.71841030E-03,-0.17004555E-03, 0.14620592E-02,
         0.94313943E-03, 0.21665540E-03,-0.10145664E-02,-0.30650580E-03,
         0.24816324E-03, 0.10238627E-02, 0.26821948E-02, 0.45195522E-03,
         0.14525824E-02,-0.28861628E-03,-0.28210404E-03,-0.40094671E-03,
         0.11944763E-03,-0.47110716E-04, 0.20725137E-03, 0.11435915E-03,
        -0.45387722E-04, 0.12042480E-03,-0.20470109E-03, 0.69184192E-04,
         0.65554312E-04, 0.64596150E-03, 0.67047827E-03, 0.45725741E-03,
         0.85069105E-03,-0.12466355E-03,-0.80779784E-04,-0.10617918E-02,
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
    float x53 = x52*x5;

//                 function

    float v_p_r5p65_fp                              =avdat
        +coeff[  0]                    
        +coeff[  1]        *x31        
        +coeff[  2]            *x41    
        +coeff[  3]                *x51
        +coeff[  4]*x11                
        +coeff[  5]    *x22            
        +coeff[  6]    *x21*x31        
        +coeff[  7]    *x21    *x41    
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[  8]            *x42    
        +coeff[  9]    *x21        *x51
        +coeff[ 10]            *x41*x51
        +coeff[ 11]*x11*x21            
        +coeff[ 12]    *x23            
        +coeff[ 13]                *x52
        +coeff[ 14]*x11        *x41    
        +coeff[ 15]    *x22    *x41    
        +coeff[ 16]    *x22        *x51
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 17]    *x22*x31    *x51
        +coeff[ 18]        *x33*x41    
        +coeff[ 19]        *x33    *x51
        +coeff[ 20]    *x24*x31        
        +coeff[ 21]            *x43*x52
        +coeff[ 22]        *x31*x41    
        +coeff[ 23]        *x31    *x51
        +coeff[ 24]    *x21    *x42    
        +coeff[ 25]    *x24            
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 26]*x11*x21    *x41    
        +coeff[ 27]    *x25            
        +coeff[ 28]*x11*x24            
        +coeff[ 29]*x11*x25            
        +coeff[ 30]    *x21            
        +coeff[ 31]    *x22*x31        
        +coeff[ 32]    *x21*x31*x41    
        +coeff[ 33]*x11            *x51
        +coeff[ 34]            *x42*x51
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 35]*x12                
        +coeff[ 36]*x11*x21*x31        
        +coeff[ 37]        *x31    *x52
        +coeff[ 38]            *x41*x52
        +coeff[ 39]*x11*x21        *x51
        +coeff[ 40]        *x34        
        +coeff[ 41]            *x45    
        +coeff[ 42]    *x22*x33*x41    
        +coeff[ 43]            *x43    
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 44]*x11*x22            
        +coeff[ 45]*x11        *x42    
        +coeff[ 46]    *x21    *x43    
        +coeff[ 47]    *x22    *x41*x51
        +coeff[ 48]*x11*x23            
        +coeff[ 49]*x11*x22    *x41    
        +coeff[ 50]    *x21    *x41*x52
        +coeff[ 51]            *x41*x53
        +coeff[ 52]        *x32        
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 53]        *x31*x42    
        +coeff[ 54]    *x21    *x41*x51
        +coeff[ 55]    *x21        *x52
        +coeff[ 56]    *x23*x31        
        +coeff[ 57]    *x23    *x41    
        +coeff[ 58]*x11    *x31*x41    
        +coeff[ 59]    *x22*x31*x41    
        +coeff[ 60]    *x22    *x42    
        +coeff[ 61]    *x23        *x51
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 62]    *x21*x31*x42    
        +coeff[ 63]    *x22        *x52
        +coeff[ 64]*x12        *x41    
        +coeff[ 65]*x11*x21    *x42    
        +coeff[ 66]    *x23    *x42    
        +coeff[ 67]*x11*x21    *x41*x51
        +coeff[ 68]    *x22    *x42*x51
        +coeff[ 69]*x12*x22            
        +coeff[ 70]*x11*x23*x31*x41    
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 71]            *x45*x51
        +coeff[ 72]*x12            *x53
        +coeff[ 73]*x11    *x31        
        +coeff[ 74]    *x22*x32        
        +coeff[ 75]                *x53
        +coeff[ 76]*x11    *x31    *x51
        +coeff[ 77]*x11        *x41*x51
        +coeff[ 78]    *x21*x31*x41*x51
        +coeff[ 79]*x11            *x52
    ;
    v_p_r5p65_fp                              =v_p_r5p65_fp                              
        +coeff[ 80]*x12    *x31        
        +coeff[ 81]    *x24    *x41    
        +coeff[ 82]*x11*x21*x31*x41    
        +coeff[ 83]    *x23*x31*x41    
        +coeff[ 84]    *x22    *x43    
        +coeff[ 85]        *x31    *x53
        +coeff[ 86]        *x34*x41    
        +coeff[ 87]    *x25    *x41    
        ;

    return v_p_r5p65_fp                              ;
}
float l_r6_h90s_fp                            (float *x,int m){
    int ncoeff= 60;
    float avdat= -0.1224601E-01;
    float xmin[10]={
        -0.99483E-02,-0.47003E-01,-0.99484E-02,-0.29097E-01,-0.46754E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.99850E-02, 0.50561E-01, 0.99954E-02, 0.24794E-01, 0.49963E-01,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 61]={
         0.97725238E-03,-0.28177884E+00,-0.33243656E-01, 0.24316324E-01,
        -0.36374688E-01, 0.79770543E-01,-0.26077900E-01,-0.24519954E-01,
        -0.84631955E-02, 0.42308006E-02, 0.20791048E-01, 0.11326477E-01,
         0.34268346E-01,-0.99035818E-02,-0.39993413E-01,-0.21219202E-02,
         0.12262885E-02,-0.49721370E-02,-0.28496159E-02, 0.29656182E-02,
         0.14563563E-01, 0.48898305E-02,-0.15408170E-01,-0.27373284E-02,
        -0.11038877E-02, 0.65900222E-03, 0.14654917E-02, 0.17425307E-02,
         0.65269405E-02,-0.23330017E-02, 0.32693150E-02,-0.16446104E-02,
        -0.27858040E-02, 0.63477177E-02, 0.81897145E-02,-0.15924384E-02,
        -0.35124230E-02,-0.35454781E-03,-0.43383561E-03,-0.81732951E-03,
        -0.19619593E-02, 0.18873460E-02, 0.37982611E-02,-0.18645512E-02,
        -0.62273988E-02, 0.84140331E-04,-0.88483974E-03,-0.39453633E-03,
        -0.11201642E-02, 0.56022446E-03,-0.69960125E-03, 0.22702177E-03,
         0.10768151E-02, 0.17525160E-02, 0.15261233E-02, 0.11015700E-02,
         0.55371865E-03,-0.74024283E-03,-0.55329630E-03,-0.26044340E-02,
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
    float x53 = x52*x5;
    float x54 = x53*x5;

//                 function

    float v_l_r6_h90s_fp                            =avdat
        +coeff[  0]                    
        +coeff[  1]    *x21            
        +coeff[  2]                *x51
        +coeff[  3]*x11                
        +coeff[  4]    *x22            
        +coeff[  5]    *x21    *x41    
        +coeff[  6]                *x52
        +coeff[  7]    *x23            
    ;
    v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
        +coeff[  8]    *x23*x31        
        +coeff[  9]            *x41    
        +coeff[ 10]    *x21*x31        
        +coeff[ 11]*x11        *x41    
        +coeff[ 12]    *x21    *x42    
        +coeff[ 13]*x11*x22            
        +coeff[ 14]    *x23    *x41    
        +coeff[ 15]    *x23*x31*x41    
        +coeff[ 16]        *x31        
    ;
    v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
        +coeff[ 17]            *x42    
        +coeff[ 18]            *x41*x51
        +coeff[ 19]*x11    *x31        
        +coeff[ 20]    *x21*x31*x41    
        +coeff[ 21]                *x53
        +coeff[ 22]*x11*x22    *x41    
        +coeff[ 23]    *x21        *x51
        +coeff[ 24]*x11            *x51
        +coeff[ 25]        *x31*x41    
    ;
    v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
        +coeff[ 26]*x11*x21            
        +coeff[ 27]    *x21*x32        
        +coeff[ 28]    *x22    *x41    
        +coeff[ 29]    *x21        *x52
        +coeff[ 30]*x11        *x42    
        +coeff[ 31]*x12*x21            
        +coeff[ 32]    *x24            
        +coeff[ 33]    *x21*x31*x42    
        +coeff[ 34]    *x21    *x43    
    ;
    v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
        +coeff[ 35]*x11*x23            
        +coeff[ 36]*x11*x22*x31        
        +coeff[ 37]*x11    *x33*x41    
        +coeff[ 38]*x12                
        +coeff[ 39]    *x21*x31    *x51
        +coeff[ 40]*x11*x21    *x41    
        +coeff[ 41]*x11    *x31*x41    
        +coeff[ 42]    *x22    *x42    
        +coeff[ 43]*x12*x21    *x41    
    ;
    v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
        +coeff[ 44]    *x24    *x41    
        +coeff[ 45]*x11*x23*x31        
        +coeff[ 46]*x11        *x41*x53
        +coeff[ 47]    *x23*x32*x41    
        +coeff[ 48]    *x21    *x41*x51
        +coeff[ 49]        *x31*x41*x51
        +coeff[ 50]*x11*x21*x31        
        +coeff[ 51]*x11    *x32        
        +coeff[ 52]    *x22*x31*x41    
    ;
    v_l_r6_h90s_fp                            =v_l_r6_h90s_fp                            
        +coeff[ 53]    *x21*x32*x41    
        +coeff[ 54]    *x21    *x42*x51
        +coeff[ 55]    *x21    *x41*x52
        +coeff[ 56]    *x21        *x53
        +coeff[ 57]                *x54
        +coeff[ 58]*x12*x21*x31        
        +coeff[ 59]    *x23    *x42    
        ;

    return v_l_r6_h90s_fp                            ;
}
