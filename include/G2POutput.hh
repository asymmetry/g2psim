#ifndef G2P_OUTPUT_H
#define G2P_OUTPUT_H

#include <vector>

#include "TObject.h"

using namespace std;

class TTree;
class G2PVar;

class G2POutput : public TObject
{
public:
    G2POutput();
    ~G2POutput();

    int Init();
    int Process();
    int End();

    bool TreeDefined() const { return fTree!=0; }
    TTree* GetTree() const { return fTree; }

protected:
    int Attach();

    int fNVar;
    double* fVar;
    vector<const char*> fVName;
    vector<G2PVar*> fVariables;

    TTree *fTree;

private:
    ClassDef(G2POutput, 1)
};

#endif
