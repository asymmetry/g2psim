#ifndef G2P_PMODSQFS_H
#define G2P_PMODSQFS_H

class G2PPModsQFS
{
public:
    G2PPModsQFS();
    ~G2PPModsQFS();
    
    void operator() (int Z, int A, double Eb, double Ef, double theta, double *XS);
    void operator() (int Z, int A, double Eb, double Ef, double theta, double EPS, double EPSD, double FP, double *XS);
    void operator() (int Z, int A, double Eb, double Ef, double theta, double Tb, double Ta, double *XS);
    void operator() (int Z, int A, double Eb, double Ef, double theta, double EPS, double EPSD, double FP, double Tb, double Ta, double *XS);

private:
    double fEPS, fEPSD, fFP, fTb, fTa;
}; 

#endif
