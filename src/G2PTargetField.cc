#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"

#include "G2PHallBField.hh"

#include "G2PTargetField.hh"

//#define FIELD_DEBUG 1

using namespace std;

const double kDEG = TMath::Pi()/180.0;
const double kCM = 1.0e-2;

ClassImp(G2PTargetField);

G2PTargetField::G2PTargetField()
    :bIsInit(false), pMapFile(NULL), fZMin(0.0), fZMax(2.99),
     fRMin(2.99), fRMax(2.99), fStepZ(0.01), fStepR(0.01), nNumZ(0),
     nNumR(0), bRotation(false), fRatio(1.0)
{
    memset(fOrigin, 0, sizeof(fOrigin));
    memset(fEulerAngle, 0, sizeof(fEulerAngle));
    memset(fRotationMatrix, 0, sizeof(fRotationMatrix));
    fBField = NULL;
}

G2PTargetField::G2PTargetField(const char* name)
    :bIsInit(false), fZMin(0.0), fZMax(2.99), fRMin(0.0), fRMax(2.99),
     fStepZ(0.01), fStepR(0.01), nNumZ(0), nNumR(0), bRotation(false),
     fRatio(1.0)
{
    memset(fOrigin, 0, sizeof(fOrigin));
    memset(fEulerAngle, 0, sizeof(fEulerAngle));
    memset(fRotationMatrix, 0, sizeof(fRotationMatrix));
    fBField = NULL;

    map<string, int> name_map;
    name_map["uniform"] = 0;
    name_map["hallb"] = 1;

    switch (name_map[name]){
    case 0:
        pfFieldSelector = &G2PTargetField::CreateUniformMap;
        pMapFile = "uniformfield.map";
        break;
    case 1:
        pfFieldSelector = &G2PTargetField::CreateHallBMap;
        pMapFile = "hallbfield.map";
        break;
    default:
        pfFieldSelector = NULL;
        pMapFile = name;
        break;
    }
}

G2PTargetField::~G2PTargetField()
{
    // Nothing to do
}

void G2PTargetField::SetEulerAngle(double alpha, double beta, double gamma)
{
    // The Euler angle is defined using Z-X'-Z" convention
    fEulerAngle[0] = alpha*kDEG;
    fEulerAngle[1] = beta*kDEG;
    fEulerAngle[2] = gamma*kDEG;

    double s1 = sin(fEulerAngle[0]);
    double c1 = cos(fEulerAngle[0]);
    double s2 = sin(fEulerAngle[1]);
    double c2 = cos(fEulerAngle[1]);
    double s3 = sin(fEulerAngle[2]);
    double c3 = cos(fEulerAngle[2]);

    fRotationMatrix[0][0][0] =  c1*c3-c2*s1*s3;
    fRotationMatrix[0][1][0] = -c1*s3-c2*s1*c3;
    fRotationMatrix[0][2][0] =  s2*s1;
    fRotationMatrix[0][0][1] =  s1*c3+c2*c1*s3;
    fRotationMatrix[0][1][1] = -s1*s3+c2*c1*c3;
    fRotationMatrix[0][2][1] = -s2*c1;
    fRotationMatrix[0][0][2] =  s2*s3;
    fRotationMatrix[0][1][2] =  s2*c3;
    fRotationMatrix[0][2][2] =  c2;

    // inverse matrix
    fRotationMatrix[1][0][0] =  c1*c3-c2*s1*s3;
    fRotationMatrix[1][0][1] = -c1*s3-c2*c3*s1;
    fRotationMatrix[1][0][2] =  s1*s2;
    fRotationMatrix[1][1][0] =  c3*s1+c1*c2*s3;
    fRotationMatrix[1][1][1] =  c1*c2*c3-s1*s3;
    fRotationMatrix[1][1][2] = -c1*s2;
    fRotationMatrix[1][2][0] =  s2*s3;
    fRotationMatrix[1][2][1] =  c3*s2;
    fRotationMatrix[1][2][2] =  c2;

    bRotation = true;
}

void G2PTargetField::Init()
{
    bool isinit = false;
    
    nNumZ = int((fZMax-fZMin)/fStepZ+1e-8)+1;
    nNumR = int((fRMax-fRMin)/fStepR+1e-8)+1;

    fBField = new double**[nNumR];
    for (int i = 0; i<nNumR; i++) {
        fBField[i] = new double*[nNumZ];
        for (int j = 0; j<nNumZ; j++) {
            fBField[i][j] = new double[5];
            for (int k = 0; k<5; k++) fBField[i][j][k] = 0.0;
        }
    }

    if (ReadMap()) isinit = true;
    else if(CreateMap()) isinit = true;

#ifdef FIELD_DEBUG
    if (isinit) SaveRootFile();
#endif

    bIsInit = isinit;
}

void G2PTargetField::End()
{
    if (fBField!=NULL) {
        for (int i = 0; i<nNumR; i++) {
            for (int j = 0; j<nNumZ; j++) {
                delete[] fBField[i][j];
                fBField[i][j] = NULL;
            }
            delete[] fBField[i];
            fBField[i] = NULL;
        }
        delete[] fBField;
        fBField = NULL;
    }
}

void G2PTargetField::GetField(const double* x, double* b)
{
    double pos[3], field[3];
    pos[0] = x[0]; pos[1] = x[1]; pos[2] = x[2];
    
    TransLab2Field(pos);

#ifdef FIELD_DEBUG
    printf("G2PTargetField: %e\t%e\t%e\n", x[0], x[1], x[2]);
#endif
    
    Interpolate(pos, field, 2);

    TransField2Lab(field);

    b[0] = field[0]*fRatio;
    b[1] = field[1]*fRatio;
    b[2] = field[2]*fRatio;

#ifdef FIELD_DEBUG
    printf("G2PTargetField: %e\t%e\t%e\n", b[0], b[1], b[2]);
#endif
}

bool G2PTargetField::ReadMap()
{
    ifstream ifs;
    int count = 0;
    
    ifs.open(pMapFile);
    if (ifs.fail()) return false;

    const int LEN = 300;
    char buff[LEN];
    
    ifs.getline(buff, LEN); // eat the first line

    while (ifs.getline(buff, LEN)!=NULL) {
        TString tmpline(buff);

        if (tmpline.EndsWith("\n")) tmpline.Chop();

        istringstream ist(tmpline.Data());
        string tmpstr;
        vector<string> line_spl;
        while (ist>>tmpstr) {
            line_spl.push_back(tmpstr);
        }

        if (line_spl.empty()) continue; // ignore empty lines

        double tempZ = atof(line_spl[0].c_str())*kCM;
        double tempR = atof(line_spl[1].c_str())*kCM;

        if ((tempZ>=fZMin)&&(tempZ<=fZMax)&&(tempR>=fRMin)&&(tempR<=fRMax)) {
            // store the value
            double tempBz = atof(line_spl[2].c_str());
            double tempBr = atof(line_spl[3].c_str());
            double tempB  = atof(line_spl[4].c_str());

            int indexZ = int((tempZ-fZMin)/fStepZ+1e-8);
            int indexR = int((tempR-fRMin)/fStepR+1e-8);

            fBField[indexR][indexZ][0] = tempZ;
            fBField[indexR][indexZ][1] = tempR;
            fBField[indexR][indexZ][2] = tempBz;
            fBField[indexR][indexZ][3] = tempBr;
            fBField[indexR][indexZ][4] = tempB;

#ifdef FIELD_DEBUG
            printf("G2PTargetField: %e\t%e\t%e\t%e\t%e\n", tempZ, tempR, tempBz,tempBr, tempB);
#endif

            count++;
        }
    }

    ifs.close();

    if (count==0) return false;
    
    return true;  
}

bool G2PTargetField::CreateHallBMap()
{
    FILE* fp;

    if ((fp=fopen(pMapFile,"w"))==NULL) return false;

    fprintf(fp, "   z        r       Bz              Br              Btot\n");
    double x[3], b[3];

    x[1] = 0;
    
    for (int i = 0; i<nNumR; i++) {
        x[0] = i*fStepR;
        for (int j = 0; j<nNumZ; j++) {
            x[2] = j*fStepZ;

            GetHallBField(x, b);

            fBField[i][j][0] = j*fStepZ;
            fBField[i][j][1] = i*fStepR;
            fBField[i][j][2] = b[2];
            fBField[i][j][3] = sqrt(b[0]*b[0]+b[1]*b[1]);
            fBField[i][j][4] = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
            
            fprintf(fp, "%8.3f %8.3f\t%e\t%e\t%e\n", x[2]/kCM, x[0]/kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);
        }
    }

    fclose(fp);

    return true;
}

bool G2PTargetField::CreateUniformMap()
{
    FILE* fp;

    if ((fp=fopen(pMapFile,"w"))==NULL) return false;

    fprintf(fp, "   z        r       Bz              Br              Btot\n");

    for (int i = 0; i<nNumR; i++) {
        for (int j = 0; j<nNumZ; j++) {
            fBField[i][j][0] = j*fStepZ;
            fBField[i][j][1] = i*fStepR;
            fBField[i][j][2] = 1.0;
            fBField[i][j][3] = 0.0;
            fBField[i][j][4] = 1.0;
            
            fprintf(fp, "%8.3f %8.3f\t%e\t%e\t%e\n", fBField[i][j][0]/kCM, fBField[i][j][1]/kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);
        }
    }

    fclose(fp);

    return true;
}

bool G2PTargetField::Interpolate(const double* pos, double* b, const int order)
{

// Calculate the nth order Lagrange polynomial interpolation
// 1) Find out (B[R0][Z0],...,B[Rn][Zn]), which is a (n+1)by(n+1) matrix
// 2) Interpolate by Z first, get (B[R0][Z],...,B[Rn][Z]), which is a 1x(n+1)
//+matrix
// 3) Interpolate by R to get B

    // Only interpolate in the first octant 
    double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    double z = fabs(pos[2]);

    if ((r>fRMax)||(z>fZMax)||(r<fRMin)||(z<fZMin)) {
        b[0] = 0.0; b[1] = 0.0; b[2] = 0.0;
        return false;
    }
    
    int indexZ0 = int((z-fZMin)/fStepZ)-int((order-1)/2);
    int indexR0 = int((r-fRMin)/fStepR)-int((order-1)/2);

    if (indexZ0<0) indexZ0 = 0;
    if (indexR0<0) indexR0 = 0;
    if (indexZ0+order>nNumZ-1) indexZ0 = nNumZ-1-order;
    if (indexR0+order>nNumR-1) indexR0 = nNumR-1-order;

    // Lagrange polynomial interpolate on Z
    int indexZ, indexR, indexT;
    double tempBz[5], tempBr[5];
    for (int i = 0; i<=order; i++){
        indexR = indexR0+i;
        tempBz[i] = 0.0;
        tempBr[i] = 0.0;
        for (int j = 0; j<=order; j++) {
            double temp = 1.0;
            indexZ = indexZ0+j;
            for (int k = 0; k<=order; k++) {
                indexT = indexZ0+k;
                if (k!=j) temp*=(z-fBField[indexR][indexT][0])/(fBField[indexR][indexZ][0]-fBField[indexR][indexT][0]);
            }
            tempBz[i]+=temp*fBField[indexR][indexZ][2];
            tempBr[i]+=temp*fBField[indexR][indexZ][3];
        }
    }

    // R
    double Bz = 0, Br = 0;
    indexZ = indexZ0;
    for (int i = 0; i<=order; i++) {
        indexR = indexR0+i;
        double temp = 1.0;
        for (int j = 0; j<=order; j++) {
            indexT = indexR0+j;
            if (j!=i) temp*=(r-fBField[indexT][indexZ][1])/(fBField[indexR][indexZ][1]-fBField[indexT][indexZ][1]);
        }
        Bz+=temp*tempBz[i];
        Br+=temp*tempBr[i];
    }

    // Sign has already been considered
    b[0] = Br*(pos[0]/r)*(pos[2]/z);
    b[1] = Br*(pos[1]/r)*(pos[2]/z);
    b[2] = Bz;

    return  true;
}

void G2PTargetField::TransLab2Field(double* x)
{
    for (int i = 0; i<3; i++) x[i]-= fOrigin[i];
    if (bRotation) {
        double temp[3];
        for (int i = 0; i<3; i++) temp[i] = fRotationMatrix[0][i][0]*x[0]+fRotationMatrix[0][i][1]*x[1]+fRotationMatrix[0][i][2]*x[2];
        for (int i = 0; i<3; i++) x[i] = temp[i];
    }
}

void G2PTargetField::TransField2Lab(double* x)
{
    if (bRotation) {
        double temp[3];
        for (int i = 0; i<3; i++) temp[i] = fRotationMatrix[1][i][0]*x[0]+fRotationMatrix[1][i][1]*x[1]+fRotationMatrix[1][i][2]*x[2];
        for (int i = 0; i<3; i++) x[i] = temp[i];
    }
}

void G2PTargetField::SaveRootFile()
{
    TFile* file = new TFile("field.root","RECREATE");
	TTree* field = new TTree("field","field map");

    double x, y, z, r;
    double Bx, By, Bz, Br, Btot;
    
	field->Branch("x",&x,"x/D");
	field->Branch("y",&y,"y/D");
	field->Branch("r",&r,"r/D");
	field->Branch("z",&z,"z/D");
	field->Branch("Bx",&Bx,"Bx/D");
	field->Branch("By",&By,"By/D");
	field->Branch("Br",&Br,"Br/D");
	field->Branch("Bz",&Bz,"Bz/D");
	field->Branch("Btot",&Btot,"Btot/D");

    for (x = -1; x<=1; x+=0.02)
        for (y = -1; y<=1; y+=0.02)
            for (z = -1; z<=1;z+=0.02) {
                double xx[3];
                xx[0] = x;
                xx[1] = y;
                xx[2] = z;
                double b[3];
                GetField(xx, b);
                Bx = b[0];
                By = b[1];
                Bz = b[2];
                Btot = sqrt(Bx*Bx+By*By+Bz*Bz);
                Br = sqrt(Bx*Bx+By*By);
                field->Fill();
            }
    file->Write();
    file->Close();
}
