#ifndef G2P_HALLBFIELD_H
#define G2P_HALLBFIELD_H

class G2PHallBField
{
public:
    G2PHallBField();
    ~G2PHallBField();
    
    void operator() (const double* pos, double* field);
};

#endif
