#ifndef TESTSNAKE_H
#define TESTSNAKE_H

// definition:
// iArm: Set Arm 0 means left arm, 1 means right arm
// iBPM: Set bpm resolution
// iDirection: Set sim direction: 0 means forward, 1 means backward, 2 means both
// iSource: Set data source: 1,2,3 means use delta,flat,gaussian distribution for input position and angle, 3 means use real data in file input_tg.dat and input_fp.dat
// iSetting: Set experiment setting: 10 means normal 484816 septa, 11 means 484816 septa with shim, 12 means 403216 septa with shim, 13 means 400016 septa with shim

void TestSNAKE(int iNEvent, int iArm, int iSetting, int iSource, int iDirection, double pHRSMomentum, double pBPMRes);

#endif
