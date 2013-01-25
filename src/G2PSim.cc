#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "G2PGun.hh"
#include "G2PRand.hh"
#include "HRSRecUseDB.hh"
#include "HRSTransTCSNHCS.hh"

#include "../G2PXSection/G2PXSection.hh"
#include "../HRSTransport/HRSTransport.hh"

#include "G2PSim.hh"

//#define G2PSIM_DEBUG 1

using namespace Transform;

const double cDeg = TMath::Pi()/180.0;

ClassImp(G2PSim);

G2PSim::G2PSim()
    :bIsInit(false), pFile(NULL), pFileName(NULL), nIndex(1),
     nEvent(10000), bIsLeftArm(true), fHRSAngle(5.767*cDeg),
     fHRSMomentum(2.251), pGun(NULL), pHRS(NULL), pRecUseDB(NULL),
     pTree(NULL), pConfig(NULL), pRand(NULL), pfRunSelector(NULL)
{
    Clear();
}

G2PSim::~G2PSim()
{
    // Nothing to be done
}

void G2PSim::Init()
{
    bool noerror = true;
    
    pRand = new G2PRand();

    pGun->SetRand(pRand);
    pGun->SetHRSAngle(fHRSAngle);
    pGun->Init();

    pHRS->SetArm(bIsLeftArm);

    if (bIsLeftArm)
        pRecUseDB = new HRSRecUseDB("L","db_L.vdc.dat");
    else
        pRecUseDB = new HRSRecUseDB("R","db_R.vdc.dat");

    if (!pGun->IsInit()) noerror = false;

    if (noerror) {
        InitTree();
        if (pGun->IsUsingData())
            pfRunSelector = &G2PSim::RunData;
        else
            pfRunSelector = &G2PSim::RunSim;
    }

    bIsInit = noerror;
}

void G2PSim::Run()
{
    if ((pGun==NULL)||(pHRS==NULL)) {
        //
    }
    else {
        Init();
        if (bIsInit) {
            (this->*pfRunSelector)();
            End();
        }
    }
}

void G2PSim::End()
{
    pFile->Write();
    pFile->Close();
}

void G2PSim::Clear()
{
    memset(fV3bpm_lab, 0, sizeof(fV3bpm_lab));
    memset(fV5tg_tr, 0, sizeof(fV5tg_tr));
    memset(fV5tg_lab, 0, sizeof(fV5tg_lab));
    memset(fV5fpdata_tr, 0, sizeof(fV5fpdata_tr));
    memset(fV5fpdata_rot, 0, sizeof(fV5fpdata_rot));

    bIsGoodParticle = false;
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof(fV5fp_rot));
    memset(fV5rec_tr, 0, sizeof(fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof(fV5rec_lab));
    
    memset(fV5recdb_tr, 0, sizeof(fV5recdb_tr));
    memset(fV5recdb_lab, 0, sizeof(fV5recdb_lab));

    fXS = 0;
}

void G2PSim::InitTree()
{
    pFile = new TFile(pFileName, "recreate");
	pTree = new TTree("T", "sim result");
    pConfig = new TTree("config", "sim configure");

    pConfig->Branch("NEvent", &nEvent, "NEvent/I");

 	pConfig->Branch("IsLeftArm", &bIsLeftArm, "IsLeftArm/O");   
    pConfig->Branch("HRSAngle", &fHRSAngle, "HRSAngle/D");
    pConfig->Branch("HRSMomentum", &fHRSMomentum, "HRSMomentum/D");

    int hrssetting = pHRS->GetModelIndex();
    pConfig->Branch("HRSModel", &hrssetting, "HRSModel/I");
    
    int gunsetting = pGun->GetSetting();
    double posres = pGun->GetPosResolution();
    double angleres = pGun->GetAngleResolution();
    double deltares = pGun->GetDeltaResolution();
    pConfig->Branch("GunSetting", &gunsetting, "GunSetting/I");
	pConfig->Branch("GunPosRes", &posres, "GunPosRes/D");
	pConfig->Branch("GunAngleRes", &angleres, "GunAngleRes/D");
	pConfig->Branch("GunDeltaRes", &deltares, "GunDeltaRes/D");

    pConfig->Fill();

    pTree->Branch("Index", &nIndex,"Index/I");
    pTree->Branch("IsGood", &bIsGoodParticle, "IsGood/O");

    pTree->Branch("Xbpm_lab", &fV3bpm_lab[0],"Xbpm_lab/D");
	pTree->Branch("Ybpm_lab", &fV3bpm_lab[1],"Ybpm_lab/D");
	pTree->Branch("Zbpm_lab", &fV3bpm_lab[2],"Zbpm_lab/D");

    pTree->Branch("Xtg_tr",&fV5tg_tr[0],"Xtg_tr/D");
	pTree->Branch("Thetatg_tr",&fV5tg_tr[1],"Thetatg_tr/D");
	pTree->Branch("Ytg_tr",&fV5tg_tr[2],"Ytg_tr/D");
	pTree->Branch("Phitg_tr",&fV5tg_tr[3],"Phitg_tr/D");
	pTree->Branch("Delta",&fV5tg_tr[4],"Delta/D");

    pTree->Branch("Xtg_lab",&fV5tg_lab[0],"Xtg_lab/D");
	pTree->Branch("Thetatg_lab",&fV5tg_lab[1],"Thetatg_lab/D");
	pTree->Branch("Ytg_lab",&fV5tg_lab[2],"Ytg_lab/D");
	pTree->Branch("Phitg_lab",&fV5tg_lab[3],"Phitg_lab/D");
	pTree->Branch("Ztg_lab",&fV5tg_lab[4],"Ztg_lab/D");

    pTree->Branch("Xfpdata_tr",&fV5fpdata_tr[0],"Xfpdata_tr/D");
	pTree->Branch("Thetafpdata_tr",&fV5fpdata_tr[1],"Thetafpdata_tr/D");
	pTree->Branch("Yfpdata_tr",&fV5fpdata_tr[2],"Yfpdata_tr/D");
	pTree->Branch("Phifpdata_tr",&fV5fpdata_tr[3],"Phifpdata_tr/D");

    pTree->Branch("Xfpdata_rot",&fV5fpdata_rot[0],"Xfpdata_rot/D");
	pTree->Branch("Thetafpdata_rot",&fV5fpdata_rot[1],"Thetafpdata_rot/D");
	pTree->Branch("Yfpdata_rot",&fV5fpdata_rot[2],"Yfpdata_rot/D");
	pTree->Branch("Phifpdata_rot",&fV5fpdata_rot[3],"Phifpdata_rot/D");
    
	pTree->Branch("Xfp_tr",&fV5fp_tr[0],"Xfp_tr/D");
	pTree->Branch("Thetafp_tr",&fV5fp_tr[1],"Thetafp_tr/D");
	pTree->Branch("Yfp_tr",&fV5fp_tr[2],"Yfp_tr/D");
	pTree->Branch("Phifp_tr",&fV5fp_tr[3],"Phifp_tr/D");

    pTree->Branch("Xfp_rot",&fV5fp_rot[0],"Xfp_rot/D");
	pTree->Branch("Thetafp_rot",&fV5fp_rot[1],"Thetafp_rot/D");
	pTree->Branch("Yfp_rot",&fV5fp_rot[2],"Yfp_rot/D");
	pTree->Branch("Phifp_rot",&fV5fp_rot[3],"Phifp_rot/D");

    pTree->Branch("Xrec_tr",&fV5rec_tr[0],"Xrec_tr/D");
	pTree->Branch("Thetarec_tr",&fV5rec_tr[1],"Thetarec_tr/D");
	pTree->Branch("Yrec_tr",&fV5rec_tr[2],"Yrec_tr/D");
	pTree->Branch("Phirec_tr",&fV5rec_tr[3],"Phirec_tr/D");
	pTree->Branch("Deltarec",&fV5rec_tr[4],"Deltarec/D");

    pTree->Branch("Xrec_lab",&fV5rec_lab[0],"Xrec_lab/D");
	pTree->Branch("Thetarec_lab",&fV5rec_lab[1],"Thetarec_lab/D");
	pTree->Branch("Yrec_lab",&fV5rec_lab[2],"Yrec_lab/D");
	pTree->Branch("Phirec_lab",&fV5rec_lab[3],"Phirec_lab/D");
	pTree->Branch("Zrec_lab",&fV5tg_lab[4],"Zrec_lab/D");

    pTree->Branch("Xrecdb_tr",&fV5recdb_tr[0],"Xrecdb_tr/D");
	pTree->Branch("Thetarecdb_tr",&fV5recdb_tr[1],"Thetarecdb_tr/D");
	pTree->Branch("Yrecdb_tr",&fV5recdb_tr[2],"Yrecdb_tr/D");
	pTree->Branch("Phirecdb_tr",&fV5recdb_tr[3],"Phirecdb_tr/D");
	pTree->Branch("Deltarecdb",&fV5recdb_tr[4],"Deltarecdb/D");

    pTree->Branch("Xrecdb_lab",&fV5recdb_lab[0],"Xrecdb_lab/D");
	pTree->Branch("Thetarecdb_lab",&fV5recdb_lab[1],"Thetarecdb_lab/D");
	pTree->Branch("Yrecdb_lab",&fV5recdb_lab[2],"Yrecdb_lab/D");
	pTree->Branch("Phirecdb_lab",&fV5recdb_lab[3],"Phirecdb_lab/D");
    pTree->Branch("Zrecdb_lab",&fV5recdb_lab[3],"Zrecdb_lab/D");

    pTree->Branch("XS", &fXS, "XS/D");
}

void G2PSim::RunSim()
{
    while (nIndex<=nEvent) {
        if (!pGun->Shoot(fV3bpm_lab, fV5tg_tr)) break;

        fV5tg_tr[1] = tan(fV5tg_tr[1]);
        fV5tg_tr[3] = tan(fV5tg_tr[3]);

        bIsGoodParticle = pHRS->Forward(fV5tg_tr, fV5fp_tr);
        pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

#ifdef G2PSIM_DEBUG
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5fp_rot[0], fV5fp_rot[1], fV5fp_rot[2], fV5fp_rot[3],fV5fp_rot[4]);
#endif

        fV5fp_tr[4] = fV5tg_tr[0];
        bIsGoodParticle &= pHRS->Backward(fV5fp_tr, fV5rec_tr);
        pRecUseDB->CalcTargetCoords(fV5fp_rot, fV5recdb_tr);

#ifdef G2PSIM_DEBUG
        bIsGoodParticle = true;
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3],fV5rec_tr[4]);
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n\n", fV5recdb_tr[0], fV5recdb_tr[1], fV5recdb_tr[2], fV5recdb_tr[3],fV5recdb_tr[4]);
#endif

        pTree->Fill();

        if ((nIndex%10000)==0) printf("%d\n", nIndex);
        nIndex++;
    }
}

void G2PSim::RunData()
{
    while (nIndex<=nEvent) {
        if (!pGun->Shoot(fV3bpm_lab, fV5fpdata_tr)) break;

        pRecUseDB->TransTr2Rot(fV5fpdata_tr, fV5fpdata_rot);
        pRecUseDB->CalcTargetCoords(fV5fpdata_rot, fV5tg_tr);

        bIsGoodParticle = pHRS->Forward(fV5tg_tr, fV5fp_tr);
        bIsGoodParticle = true;
        pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

#ifdef G2PSIM_DEBUG
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5fp_rot[0], fV5fp_rot[1], fV5fp_rot[2], fV5fp_rot[3],fV5fp_rot[4]);
#endif

        double pV3[3];
        X_HCS2TCS(fV3bpm_lab[0], fV3bpm_lab[1], fV3bpm_lab[2], fHRSAngle, pV3[0], pV3[1], pV3[2]);
        fV5fpdata_tr[4] = pV3[0];
        
        bIsGoodParticle &= pHRS->Backward(fV5fpdata_tr, fV5rec_tr);
        pRecUseDB->CalcTargetCoords(fV5fpdata_rot, fV5recdb_tr);

        // debug first order matrix
        // double x = fV5fpdata_rot[0];
        // double t = fV5fpdata_rot[1];
        // double y = fV5fpdata_rot[2];
        // double p = fV5fpdata_rot[3];
        // double theta =
        //     +17.48692267240774*x*y*y*y  -4256.547756961809*t*y*y*y
        //     -53.08151549058515*p*y*y*y  -12.00336809935114*y*y*y
        //     -2.531635054380832*x*x*y*y  +407.2179196758938*t*x*y*y
        //     -35.40416955221132*p*x*y*y  +3.746203014043834*x*y*y
        //     -1979.415081762384*t*t*y*y  +3588.242501682196*p*t*y*y
        //     -216.2248844750452*t*y*y    +9.163630152410622*p*y*y
        //     +.3718881508359192*y*y      +.02694744865228241*x*x*x*y
        //     -10.28247234653471*t*x*x*y  +3.679685179351096*p*x*x*y
        //     -.06417739387337458*x*x*y   -163.9242516491097*t*t*x*y
        //     -487.6363651946646*p*t*x*y  +20.29805447384518*t*x*y
        //     +.5569414394411135*p*p*x*y  -3.159378888307759*p*x*y
        //     +.2608631481286963*x*y      -1.0931562345608137e-6*t*t*t*y
        //     -0.00161963081202737*p*t*t*y+30.95502791751757*t*t*y
        //     -.7998868733013537*p*p*t*y  +193.9365434624235*p*t*y
        //     -9.049510846030921*t*y      -131.6800111335095*p*p*p*y
        //     +3.580456906139325*p*p*y    -.5328508289514202*p*y
        //     -.002886285289093902*y      -.005059536705244698*x*x*x*x
        //     +.05062593746308184*t*x*x*x -.1151996633885437*p*x*x*x
        //     +0.008267251200365*x*x*x    -2.549815749792919*t*t*x*x
        //     +22.82053194803757*p*t*x*x  -.1299512338414789*t*x*x
        //     +.2900470084589791*p*p*x*x  +.2104480991760481*p*x*x
        //     +.01312850284081323*x*x     +100.9569632859276*t*t*t*x
        //     -.9439015875431279*p*t*t*x  +2.23112007154368*t*t*x
        //     -.07322848839775895*p*p*t*x -19.22564569062186*p*t*x
        //     +.6097347399877787*t*x      -.1411579432948645*p*p*p*x
        //     -.4187282396986695*p*p*x    -.2460694264573534*p*x
        //     +.01967314912792312*x       -5.705258709581791e-16*t*t*t*t
        //     +5.028505776493215*t*t*t    +113.9896671912491*p*t*t
        //     -7.44991666332372*t*t       +31.83284036371813*p*p*t
        //     -3.188034475496859*p*t      -2.330275574146642*t
        //     +26.37863589553575*p*p*p    -.6363653274747633*p*p
        //     -.04511461370913299*p       -.001384583770539381;
        // double phi =
        //     -57.49062074013359*y*y*y*y  -41.25822821622984*x*y*y*y
        //     +76.32579848133473*t*y*y*y  -33.26466741474866*p*y*y*y
        //     -47.92459307699169*y*y*y    +1.141695705176272*x*x*y*y
        //     +42.70549628115122*t*x*y*y  +31.60575876162463*p*x*y*y
        //     -3.88044439907499*x*y*y     -5254.876867731812*t*t*y*y
        //     -2626.63712948287*p*t*y*y   -209.4871104858254*t*y*y
        //     +30.49049843780868*p*p*y*y  +88.65456771072067*p*y*y
        //     -2.899251254522835*y*y      +.1134730363725824*x*x*x*y
        //     +5.118289105457366*t*x*x*y  +5.551407056874751*p*x*x*y
        //     +0.331644325169727*x*x*y    -243.389817134777*t*t*x*y
        //     -342.8276009352205*p*t*x*y  +7.239191616916298*t*x*y
        //     -38.10630757025967*p*p*x*y  +9.097384575461048*p*x*y
        //     +.4255909719467441*x*y      +13290.06866210993*t*t*t*y
        //     +5413.903772707059*p*t*t*y  -235.1404336622575*t*t*y
        //     -510.7604194064137*p*p*t*y  +384.2784288044106*p*t*y
        //     -5.577964573428799*t*y      +11.0402753330995*p*p*p*y
        //     -1.955059344396998*p*p*y    +13.35887790631663*p*y
        //     -.6644783223026642*y        -9.519481124758074e-4*x*x*x*x
        //     +0.0890421138965983*t*x*x*x +.02081427040640938*p*x*x*x
        //     -0.00422373506997289*x*x*x  +.2018503297213509*t*t*x*x
        //     +5.274915178643631*p*t*x*x  +.06449396249862054*t*x*x
        //     -11.80896161350661*p*p*x*x  -.06592709224574066*p*x*x
        //     +.004972958217155742*x*x    +11.94721053559095*t*t*t*x
        //     -428.2607333175339*p*t*t*x  +2.920232186756307*t*t*x
        //     +404.0818855951677*p*p*t*x  -4.023771997304622*p*t*x
        //     +.09873708516051571*t*x     -12.64400045655924*p*p*p*x
        //     -4.18426039204174*p*p*x     -1.001531316162395*p*x
        //     +.004947863711149666*x      -4.348047736641725*t*t*t*t
        //     -2147.338627556273*p*t*t*t  +8.857353357805392*t*t*t
        //     +23.10007439357093*p*p*t*t  -82.02000262264133*p*t*t
        //     -12.98138725729325*t*t      +3802.810315883664*p*p*p*t
        //     -120.2147702300636*p*p*t    +11.26691341836424*p*t
        //     -.05232590094653503*t       +10.25057892547751*p*p*p
        //     -11.94100712711096*p*p      +.2682033305056235*p
        //     +.004576526817886658;
        // double delta =
        //     +.2714583433872912*x*x*y*y  +.4453471337389856*x*y*y
        //     -45.35499421436646*t*y*y    +.5390967291291536*y*y
        //     +.002351872066926104*x*x*x*y+.4484017518424912*t*x*x*y
        //     -.006863686910071942*p*x*x*y+9.457055438135129e-4*x*x*y
        //     -43.08886263065754*t*t*x*y  -.08295451934571095*p*t*x*y
        //     +1.158243840431781*t*x*y    +.02590501284923457*x*y
        //     +1450.423216734052*t*t*t*y  -14.01673538204914*t*t*y
        //     -2.194379840060142*t*y      -.5012941564816951*p*y
        //     -.008752375165712667*y      -.001665233458527423*x*x*x*x
        //     +0.0173374034508381*t*x*x*x -.003388283973332728*p*x*x*x
        //     -4.5418632046300234e-4*x*x*x+1.394971853843119*t*t*x*x
        //     +.02836346993644056*p*t*x*x +.01394898811430735*t*x*x
        //     -0.00483757368419453*p*x*x  +.006429688182575492*x*x
        //     -132.7421239436652*t*t*t*x  +.2384221398095166*p*t*t*x
        //     +1.263362648305126*t*t*x    -8.915127387979794e-4*p*t*x
        //     +.3611426842675042*t*x      -.007611012040635297*p*x
        //     +.07454165297708666*x       +3390.507376300769*t*t*t*t
        //     +20.9659568389848*t*t*t     -3.412630778965652*t*t
        //     +.7203924897947827*p*t      +0.0134013423973891*t
        //     -.005387411438380335*p      -5.308984885022496e-4;
        // printf("%e\t%e\t%e\n", theta, phi, delta);
        
#ifdef G2PSIM_DEBUG
        bIsGoodParticle = true;
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3],fV5rec_tr[4]);
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n\n", fV5recdb_tr[0], fV5recdb_tr[1], fV5recdb_tr[2], fV5recdb_tr[3],fV5recdb_tr[4]);
#endif

        pTree->Fill();

        if ((nIndex%10000)==0) printf("%d\n", nIndex);
        nIndex++;
    }
}

        //throw away this event if delta_rec>=1.0
		// {
		// 	Transform::X_TCS2HCS(fV5tg_tr[0],fV5tg_tr[2],0.0,fHRSAngle,fV5tg_lab[0],fV5tg_lab[2],fV5tg_lab[4]);
		// 	Transform::P_TCS2HCS(fV5tg_tr[1],fV5tg_tr[3],fHRSAngle,fV5tg_lab[1],fV5tg_lab[3]);

		// 	Transform::X_TCS2HCS(fV5rec_tr[0],fV5rec_tr[2],0.0,fHRSAngle,fV5rec_lab[0],fV5rec_lab[2],fV5rec_lab[4]);
		// 	Transform::P_TCS2HCS(fV5rec_tr[1],fV5rec_tr[3],fHRSAngle,fV5rec_lab[1],fV5rec_lab[3]);
            
		// 	//Reconstruct use optics database
		// 	pRecUseDB->CalcTargetCoords(fV5fp_tr,fV5rec_db_tr);

		// 	pTree->Fill();
		// 	nIndex++;
		// }
