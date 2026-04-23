#pragma once

#include <cmath>

struct KFParticleStruct {
  int id;
  int idMC;
  int pdg;
  float params[6];
  float covMat[36];
  bool isCovMatOK;
  float x;
  float y;
  float z;
  float pX;
  float pY;
  float pZ;
  float e;
  float s;
  int charge;
  float chi2;
  float ndf;
  float massHypo;

  bool isMother;
  bool isDaughter;

  float xErr;
  float yErr;
  float zErr;
  float pXErr;
  float pYErr;
  float pZErr;
  float p;
  float pT;
  float eta;
  float phi;
  float mass;
  float decayLength;
  float decayLengthXY;
  float lifeTime;
  float r;
  float rapidity;
  float theta;

  float eErr;
  float sErr;
  float pErr;
  float pTErr;
  float etaErr;
  float phiErr;
  float massErr;
  float decayLengthErr;
  float decayLengthXYErr;
  float lifeTimeErr;
  float rErr;
};