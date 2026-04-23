/*

Most lines of code are borrowed from
AliVParticle.h(cxx) and AliExternalTrackParam.h(cxx)
(https://alisw.github.io/)

Code to convert global tracks (x,y,z, px,py,pz, covmat[21])
to local ALICE tracks (y,z,sinp,tanl,q/pt, covmat[15]) and vice versa

*/

#pragma once

#include "TMath.h"

class ALICETrackParamsConverter {
public:
  ALICETrackParamsConverter(double xyz_[3], double pxpypz_[3], double cv_[21],
                            short sign_); // glob
  ALICETrackParamsConverter(double params_[5], double locCV_[15], double x_,
                            double alpha_, short sign_); // loc
  ~ALICETrackParamsConverter();                          //

  void GetXYZ(double xyz_[3]);
  void GetPxPyPz(double pxpypz_[3]);
  void GetGlobalCovMat(double cv_[21]);

  void GetLocalParams(double locPar_[5]);
  void GetLocalCovMat(double locC_[15]);

  void GetAlpha(double &alpha) { alpha = this->alpha; }
  void GetX(double &x) { x = this->x; }
  void GetTrackSign(short &sign) { sign = this->sign; }

  void PrintLocalParams();
  void PrintGlobalParams();

private:
  double xyz[3];
  double pxpypz[3];
  double globCV[21];

  double localParams[5];
  double locCV[15];

  double alpha;
  double x;
  short sign;

  double Alpha();
  const Double_t kAlmost1 = 1. - Double_t(FLT_EPSILON);
  const Double_t kAlmost0 = Double_t(FLT_MIN);

  template <typename T> Double_t GetSign(T x) { return (x > 0) ? 1 : -1; }

  void GlobalToLocal(Double_t xyz[3], Double_t pxpypz[3], Double_t cv[21],
                     Double_t fP[5], Double_t fC[15], Short_t sign);
  void GetCovarianceXYZPxPyPz(Double_t xyz[3], Double_t pxpypz[3],
                              Double_t cv[21], Double_t fP[5], Double_t fC[15],
                              Short_t charge);

  void CheckCovariance(double fC[15]);

  const Double_t kC0max = 100 * 100, // SigmaY<=100cm
      kC2max = 100 * 100,            // SigmaZ<=100cm
      kC5max = 1 * 1,                // SigmaSin<=1
      kC9max = 1 * 1,                // SigmaTan<=1
      kC14max = 100 * 100;           // Sigma1/Pt<=100 1/GeV
};

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//-=-=-=-=-=- Functions definitions -=-=-=-=-=-=-=-=-=-=-=-
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ALICETrackParamsConverter::ALICETrackParamsConverter(double xyz_[3],
                                                     double pxpypz_[3],
                                                     double cv_[21],
                                                     short sign_) {
  for (int i = 0; i < 3; i++) {
    xyz[i] = xyz_[i];
    pxpypz[i] = pxpypz_[i];
  }
  for (int i = 0; i < 21; i++)
    globCV[i] = cv_[i];

  x = xyz[0];
  alpha = Alpha();
  this->sign = sign_;

  GlobalToLocal(xyz_, pxpypz_, cv_, localParams, locCV, this->sign);
}

ALICETrackParamsConverter::ALICETrackParamsConverter(double params_[5],
                                                     double locCV_[15],
                                                     double x_, double alpha_,
                                                     short sign_) {
  for (int i = 0; i < 5; i++)
    localParams[i] = params_[i];
  for (int i = 0; i < 15; i++)
    locCV[i] = locCV_[i];

  this->x = x_;
  this->alpha = alpha_;
  this->sign = sign_;

  xyz[0] = x;
  // cout << "x = " << xyz[0] << endl;
  xyz[1] = params_[0];
  xyz[2] = params_[1];
  Double_t cs = TMath::Cos(alpha_), sn = TMath::Sin(alpha_), x = xyz[0];
  xyz[0] = x * cs - xyz[1] * sn;
  xyz[1] = x * sn + xyz[1] * cs;
  // LocalToGlobal(xyz, pxpypz, cv, params_, locCV_, sign_);
  // cout << "x = " << xyz[0] << endl;

  pxpypz[0] = params_[4];
  pxpypz[1] = params_[2];
  pxpypz[2] = params_[3];
  Double_t pt = 1. / TMath::Abs(pxpypz[0]);
  cs = TMath::Cos(alpha_);
  sn = TMath::Sin(alpha_);
  Double_t r = TMath::Sqrt((1. - pxpypz[1]) * (1. + pxpypz[1]));
  pxpypz[0] = pt * (r * cs - pxpypz[1] * sn);
  pxpypz[1] = pt * (pxpypz[1] * cs + r * sn);
  pxpypz[2] = pt * pxpypz[2];

  /// COVMAT CREATE /////
  GetCovarianceXYZPxPyPz(xyz, pxpypz, globCV, params_, locCV_, sign_);
}

ALICETrackParamsConverter::~ALICETrackParamsConverter() {}

double ALICETrackParamsConverter::Alpha() {
  double fAlpha = 0.;
  const double kSafe = 1e-5;
  Double_t radPos2 = xyz[0] * xyz[0] + xyz[1] * xyz[1];
  Double_t radMax = 45.; // approximately ITS outer radius

  if (radPos2 < radMax * radMax) { // inside the ITS
    fAlpha = TMath::ATan2(pxpypz[1], pxpypz[0]);
  } else { // outside the ITS
    Float_t phiPos = TMath::Pi() + TMath::ATan2(-xyz[1], -xyz[0]);
    fAlpha = TMath::DegToRad() *
             (20 * ((((Int_t)(phiPos * TMath::RadToDeg())) / 20)) + 10);
  }

  return fAlpha;
}

void ALICETrackParamsConverter::GlobalToLocal(Double_t xyz[3],
                                              Double_t pxpypz[3],
                                              Double_t cv[21], Double_t fP[5],
                                              Double_t fC[15], Short_t sign) {
  //
  // create external track parameters from the global parameters
  // x,y,z,px,py,pz and their 6x6 covariance matrix
  // A.Dainese 10.10.08

  // Calculate alpha: the rotation angle of the corresponding local system.
  //
  // For global radial position inside the beam pipe, alpha is the
  // azimuthal angle of the momentum projected on (x,y).
  //
  // For global radial position outside the ITS, alpha is the
  // azimuthal angle of the centre of the TPC sector in which the point
  // xyz lies
  //
  Double_t fAlpha = 0.;

  const double kSafe = 1e-5;
  Double_t radPos2 = xyz[0] * xyz[0] + xyz[1] * xyz[1];
  Double_t radMax = 45.;           // approximately ITS outer radius
  if (radPos2 < radMax * radMax) { // inside the ITS
    fAlpha = TMath::ATan2(pxpypz[1], pxpypz[0]);
  } else { // outside the ITS
    Float_t phiPos = TMath::Pi() + TMath::ATan2(-xyz[1], -xyz[0]);
    fAlpha = TMath::DegToRad() *
             (20 * ((((Int_t)(phiPos * TMath::RadToDeg())) / 20)) + 10);
  }
  //
  Double_t cs = TMath::Cos(fAlpha), sn = TMath::Sin(fAlpha);
  // protection:  avoid alpha being too close to 0 or +-pi/2
  if (TMath::Abs(sn) < 2 * kSafe) {
    if (fAlpha > 0)
      fAlpha += fAlpha < TMath::Pi() / 2. ? 2 * kSafe : -2 * kSafe;
    else
      fAlpha += fAlpha > -TMath::Pi() / 2. ? -2 * kSafe : 2 * kSafe;
    cs = TMath::Cos(fAlpha);
    sn = TMath::Sin(fAlpha);
  } else if (TMath::Abs(cs) < 2 * kSafe) {
    if (fAlpha > 0)
      fAlpha += fAlpha > TMath::Pi() / 2. ? 2 * kSafe : -2 * kSafe;
    else
      fAlpha += fAlpha > -TMath::Pi() / 2. ? 2 * kSafe : -2 * kSafe;
    cs = TMath::Cos(fAlpha);
    sn = TMath::Sin(fAlpha);
  }
  // Get the vertex of origin and the momentum
  TVector3 ver(xyz[0], xyz[1], xyz[2]);
  TVector3 mom(pxpypz[0], pxpypz[1], pxpypz[2]);
  //
  // Rotate to the local coordinate system
  ver.RotateZ(-fAlpha);
  mom.RotateZ(-fAlpha);

  //
  // x of the reference plane
  Double_t fX = ver.X();

  Double_t charge = (Double_t)sign;

  fP[0] = ver.Y();
  fP[1] = ver.Z();
  fP[2] = TMath::Sin(mom.Phi());
  fP[3] = mom.Pz() / mom.Pt();
  fP[4] = TMath::Sign(1 / mom.Pt(), charge);
  //
  if (TMath::Abs(1 - fP[2]) < kSafe)
    fP[2] = 1. - kSafe; // Protection
  else if (TMath::Abs(-1 - fP[2]) < kSafe)
    fP[2] = -1. + kSafe; // Protection
  //
  // Covariance matrix (formulas to be simplified)
  Double_t pt = 1. / TMath::Abs(fP[4]);
  Double_t r = TMath::Sqrt((1. - fP[2]) * (1. + fP[2]));
  //
  Double_t cv34 = TMath::Sqrt(cv[3] * cv[3] + cv[4] * cv[4]);
  //
  Int_t special = 0;
  double sgcheck = r * sn + fP[2] * cs;
  if (TMath::Abs(sgcheck) >= 1 - kSafe) { // special case: lab phi is +-pi/2
    special = 1;
    sgcheck = TMath::Sign(1.0, sgcheck);
  } else if (TMath::Abs(sgcheck) < kSafe) {
    sgcheck = TMath::Sign(1.0, cs);
    special = 2; // special case: lab phi is 0
  }
  //
  fC[0] = cv[0] + cv[2];
  fC[1] = TMath::Sign(cv34, -cv[3] * sn);
  fC[2] = cv[5];
  //
  if (special == 1) {
    double pti = 1. / pt;
    double pti2 = pti * pti;
    int q = GetSign(fP[4]);
    fC[3] = cv[6] * pti;
    fC[4] = -sgcheck * cv[8] * r * pti;
    fC[5] = TMath::Abs(cv[9] * r * r * pti2);
    fC[6] = (cv[10] * fP[3] - sgcheck * cv[15]) * pti / r;
    fC[7] = (cv[17] - sgcheck * cv[12] * fP[3]) * pti;
    fC[8] = (-sgcheck * cv[18] + cv[13] * fP[3]) * r * pti2;
    fC[9] = TMath::Abs(cv[20] - 2 * sgcheck * cv[19] * fP[3] +
                       cv[14] * fP[3] * fP[3]) *
            pti2;
    fC[10] = cv[10] * pti2 / r * q;
    fC[11] = -sgcheck * cv[12] * pti2 * q;
    fC[12] = cv[13] * r * pti * pti2 * q;
    fC[13] = (-sgcheck * cv[19] + cv[14] * fP[3]) * r * pti2 * pti;
    fC[14] = TMath::Abs(cv[14] * pti2 * pti2);
  } else if (special == 2) {
    double pti = 1. / pt;
    double pti2 = pti * pti;
    int q = GetSign(fP[4]);
    fC[3] = -cv[10] * pti * cs / sn;
    fC[4] = cv[12] * cs * pti;
    fC[5] = TMath::Abs(cv[14] * cs * cs * pti2);
    fC[6] = (sgcheck * cv[6] * fP[3] - cv[15]) * pti / sn;
    fC[7] = (cv[17] - sgcheck * cv[8] * fP[3]) * pti;
    fC[8] = (cv[19] - sgcheck * cv[13] * fP[3]) * cs * pti2;
    fC[9] = TMath::Abs(cv[20] - 2 * sgcheck * cv[18] * fP[3] +
                       cv[9] * fP[3] * fP[3]) *
            pti2;
    fC[10] = sgcheck * cv[6] * pti2 / sn * q;
    fC[11] = -sgcheck * cv[8] * pti2 * q;
    fC[12] = -sgcheck * cv[13] * cs * pti * pti2 * q;
    fC[13] = (-sgcheck * cv[18] + cv[9] * fP[3]) * pti2 * pti * q;
    fC[14] = TMath::Abs(cv[9] * pti2 * pti2);
  } else {
    Double_t m00 = -sn; // m10=cs;
    Double_t m23 = -pt * (sn + fP[2] * cs / r),
             m43 = -pt * pt * (r * cs - fP[2] * sn);
    Double_t m24 = pt * (cs - fP[2] * sn / r),
             m44 = -pt * pt * (r * sn + fP[2] * cs);
    Double_t m35 = pt, m45 = -pt * pt * fP[3];
    //
    m43 *= GetSign(fP[4]);
    m44 *= GetSign(fP[4]);
    m45 *= GetSign(fP[4]);
    //
    Double_t a1 = cv[13] - cv[9] * (m23 * m44 + m43 * m24) / m23 / m43;
    Double_t a2 = m23 * m24 - m23 * (m23 * m44 + m43 * m24) / m43;
    Double_t a3 = m43 * m44 - m43 * (m23 * m44 + m43 * m24) / m23;
    Double_t a4 = cv[14] + 2. * cv[9]; // cv[14]-2.*cv[9]*m24*m44/m23/m43;
    Double_t a5 = m24 * m24 - 2. * m24 * m44 * m23 / m43;
    Double_t a6 = m44 * m44 - 2. * m24 * m44 * m43 / m23;
    //
    fC[3] = (cv[10] * m43 - cv[6] * m44) / (m24 * m43 - m23 * m44) / m00;
    fC[10] = (cv[6] / m00 - fC[3] * m23) / m43;
    fC[6] = (cv[15] / m00 - fC[10] * m45) / m35;
    fC[4] = (cv[12] * m43 - cv[8] * m44) / (m24 * m43 - m23 * m44);
    fC[11] = (cv[8] - fC[4] * m23) / m43;
    fC[7] = cv[17] / m35 - fC[11] * m45 / m35;
    fC[5] = TMath::Abs((a4 * a3 - a6 * a1) / (a5 * a3 - a6 * a2));
    fC[14] = TMath::Abs((a1 - a2 * fC[5]) / a3);
    fC[12] = (cv[9] - fC[5] * m23 * m23 - fC[14] * m43 * m43) / m23 / m43;
    Double_t b1 = cv[18] - fC[12] * m23 * m45 - fC[14] * m43 * m45;
    Double_t b2 = m23 * m35;
    Double_t b3 = m43 * m35;
    Double_t b4 = cv[19] - fC[12] * m24 * m45 - fC[14] * m44 * m45;
    Double_t b5 = m24 * m35;
    Double_t b6 = m44 * m35;
    fC[8] = (b4 - b6 * b1 / b3) / (b5 - b6 * b2 / b3);
    fC[13] = b1 / b3 - b2 * fC[8] / b3;
    fC[9] =
        TMath::Abs((cv[20] - fC[14] * (m45 * m45) - fC[13] * 2. * m35 * m45) /
                   (m35 * m35));
  }

  CheckCovariance(fC);
  GetCovarianceXYZPxPyPz(xyz, pxpypz, globCV, fP, fC, sign);
  x = fX;

  return;
}

void ALICETrackParamsConverter::GetCovarianceXYZPxPyPz(
    Double_t xyz[3], Double_t pxpypz[3], Double_t cv[21], Double_t fP[5],
    Double_t fC[15], Short_t charge) {
  //---------------------------------------------------------------------
  // This function returns the global covariance matrix of the track params
  //
  // Cov(x,x) ... :   cv[0]
  // Cov(y,x) ... :   cv[1]  cv[2]
  // Cov(z,x) ... :   cv[3]  cv[4]  cv[5]
  // Cov(px,x)... :   cv[6]  cv[7]  cv[8]  cv[9]
  // Cov(py,x)... :   cv[10] cv[11] cv[12] cv[13] cv[14]
  // Cov(pz,x)... :   cv[15] cv[16] cv[17] cv[18] cv[19] cv[20]
  //
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------

  Double_t fAlpha = alpha;

  Double_t pt = 1. / TMath::Abs(fP[4]);
  Double_t cs = TMath::Cos(fAlpha), sn = TMath::Sin(fAlpha);
  Double_t r = TMath::Sqrt((1. - fP[2]) * (1. + fP[2]));

  Double_t m00 = -sn, m10 = cs;
  Double_t m23 = -pt * (sn + fP[2] * cs / r),
           m43 = -pt * pt * (r * cs - fP[2] * sn);
  Double_t m24 = pt * (cs - fP[2] * sn / r),
           m44 = -pt * pt * (r * sn + fP[2] * cs);
  Double_t m35 = pt, m45 = -pt * pt * fP[3];

  m43 *= GetSign(charge);
  m44 *= GetSign(charge);
  m45 *= GetSign(charge);

  cv[0] = fC[0] * m00 * m00;
  cv[1] = fC[0] * m00 * m10;
  cv[2] = fC[0] * m10 * m10;
  cv[3] = fC[1] * m00;
  cv[4] = fC[1] * m10;
  cv[5] = fC[2];
  cv[6] = m00 * (fC[3] * m23 + fC[10] * m43);
  cv[7] = m10 * (fC[3] * m23 + fC[10] * m43);
  cv[8] = fC[4] * m23 + fC[11] * m43;
  cv[9] =
      m23 * (fC[5] * m23 + fC[12] * m43) + m43 * (fC[12] * m23 + fC[14] * m43);
  cv[10] = m00 * (fC[3] * m24 + fC[10] * m44);
  cv[11] = m10 * (fC[3] * m24 + fC[10] * m44);
  cv[12] = fC[4] * m24 + fC[11] * m44;
  cv[13] =
      m23 * (fC[5] * m24 + fC[12] * m44) + m43 * (fC[12] * m24 + fC[14] * m44);
  cv[14] =
      m24 * (fC[5] * m24 + fC[12] * m44) + m44 * (fC[12] * m24 + fC[14] * m44);
  cv[15] = m00 * (fC[6] * m35 + fC[10] * m45);
  cv[16] = m10 * (fC[6] * m35 + fC[10] * m45);
  cv[17] = fC[7] * m35 + fC[11] * m45;
  cv[18] =
      m23 * (fC[8] * m35 + fC[12] * m45) + m43 * (fC[13] * m35 + fC[14] * m45);
  cv[19] =
      m24 * (fC[8] * m35 + fC[12] * m45) + m44 * (fC[13] * m35 + fC[14] * m45);
  cv[20] =
      m35 * (fC[9] * m35 + fC[13] * m45) + m45 * (fC[13] * m35 + fC[14] * m45);
}

void ALICETrackParamsConverter::PrintGlobalParams() {
  std::cout << "parameters:\t" << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2]
            << "\t" << pxpypz[0] << "\t" << pxpypz[1] << "\t" << pxpypz[2]
            << std::endl;
  std::cout << "sign:\t" << sign << std::endl;
  std::cout << "covmat:\t"
            << "\n"
            << globCV[0] << "\n"
            << globCV[1] << "\t" << globCV[2] << "\n"
            << globCV[3] << "\t" << globCV[4] << "\t" << globCV[5] << "\n"
            << globCV[6] << "\t" << globCV[7] << "\t" << globCV[8] << "\t"
            << globCV[9] << "\n"
            << globCV[10] << "\t" << globCV[11] << "\t" << globCV[12] << "\t"
            << globCV[13] << "\t" << globCV[14] << "\n"
            << globCV[15] << "\t" << globCV[16] << "\t" << globCV[17] << "\t"
            << globCV[18] << "\t" << globCV[19] << "\t" << globCV[20] << "\n";
  cout << endl;
}

void ALICETrackParamsConverter::PrintLocalParams() {
  std::cout << "parameters:\t" << localParams[0] << "\t" << localParams[1]
            << "\t" << localParams[2] << "\t" << localParams[3] << "\t"
            << localParams[4] << std::endl;
  std::cout << "sign:\t" << sign << std::endl;
  std::cout << "covmat:\t"
            << "\n"
            << locCV[0] << "\n"
            << locCV[1] << "\t" << locCV[2] << "\n"
            << locCV[3] << "\t" << locCV[4] << "\t" << locCV[5] << "\n"
            << locCV[6] << "\t" << locCV[7] << "\t" << locCV[8] << "\t"
            << locCV[9] << "\n"
            << locCV[10] << "\t" << locCV[11] << "\t" << locCV[12] << "\t"
            << locCV[13] << "\t" << locCV[14] << "\n";
  cout << endl;
}

//////////////////
///     Getters
//////////////////

void ALICETrackParamsConverter::GetXYZ(double xyz_[3]) {
  for (int i = 0; i < 3; i++)
    xyz_[i] = xyz[i];
}
void ALICETrackParamsConverter::GetPxPyPz(double pxpypz_[3]) {
  for (int i = 0; i < 3; i++)
    pxpypz_[i] = pxpypz[i];
}
void ALICETrackParamsConverter::GetGlobalCovMat(double cv_[21]) {
  for (int i = 0; i < 21; i++)
    cv_[i] = globCV[i];
}
void ALICETrackParamsConverter::GetLocalParams(double locPar_[5]) {
  for (int i = 0; i < 5; i++)
    locPar_[i] = localParams[i];
}
void ALICETrackParamsConverter::GetLocalCovMat(double locC_[15]) {
  for (int i = 0; i < 15; i++)
    locC_[i] = locCV[i];
}

/////// Check COvar matrix //////

void ALICETrackParamsConverter::CheckCovariance(double fC[15]) {

  // This function forces the diagonal elements of the covariance matrix to be
  // positive. In case the diagonal element is bigger than the maximal allowed
  // value, it is set to the limit and the off-diagonal elements that correspond
  // to it are set to zero.

  fC[0] = TMath::Abs(fC[0]);
  if (fC[0] > kC0max) {
    double scl = TMath::Sqrt(kC0max / fC[0]);
    fC[0] = kC0max;
    fC[1] *= scl;
    fC[3] *= scl;
    fC[6] *= scl;
    fC[10] *= scl;
  }
  fC[2] = TMath::Abs(fC[2]);
  if (fC[2] > kC2max) {
    double scl = TMath::Sqrt(kC2max / fC[2]);
    fC[2] = kC2max;
    fC[1] *= scl;
    fC[4] *= scl;
    fC[7] *= scl;
    fC[11] *= scl;
  }
  fC[5] = TMath::Abs(fC[5]);
  if (fC[5] > kC5max) {
    double scl = TMath::Sqrt(kC5max / fC[5]);
    fC[5] = kC5max;
    fC[3] *= scl;
    fC[4] *= scl;
    fC[8] *= scl;
    fC[12] *= scl;
  }
  fC[9] = TMath::Abs(fC[9]);
  if (fC[9] > kC9max) {
    double scl = TMath::Sqrt(kC9max / fC[9]);
    fC[9] = kC9max;
    fC[6] *= scl;
    fC[7] *= scl;
    fC[8] *= scl;
    fC[13] *= scl;
  }
  fC[14] = TMath::Abs(fC[14]);
  if (fC[14] > kC14max) {
    double scl = TMath::Sqrt(kC14max / fC[14]);
    fC[14] = kC14max;
    fC[10] *= scl;
    fC[11] *= scl;
    fC[12] *= scl;
    fC[13] *= scl;
  }

  // The part below is used for tests and normally is commented out
  //     TMatrixDSym m(5);
  //     TVectorD eig(5);

  //     m(0,0)=fC[0];
  //     m(1,0)=fC[1];  m(1,1)=fC[2];
  //     m(2,0)=fC[3];  m(2,1)=fC[4];  m(2,2)=fC[5];
  //     m(3,0)=fC[6];  m(3,1)=fC[7];  m(3,2)=fC[8];  m(3,3)=fC[9];
  //     m(4,0)=fC[10]; m(4,1)=fC[11]; m(4,2)=fC[12]; m(4,3)=fC[13];
  //     m(4,4)=fC[14];

  //     m(0,1)=m(1,0);
  //     m(0,2)=m(2,0); m(1,2)=m(2,1);
  //     m(0,3)=m(3,0); m(1,3)=m(3,1); m(2,3)=m(3,2);
  //     m(0,4)=m(4,0); m(1,4)=m(4,1); m(2,4)=m(4,2); m(3,4)=m(4,3);
  //     m.EigenVectors(eig);

  //     //    assert(eig(0)>=0 && eig(1)>=0 && eig(2)>=0 && eig(3)>=0 &&
  //     eig(4)>=0); if (!(eig(0)>=0 && eig(1)>=0 && eig(2)>=0 && eig(3)>=0 &&
  //     eig(4)>=0)) {
  //       AliWarning("Negative eigenvalues of the covariance matrix!");
  //       this->Print();
  //       eig.Print();
  //     }
}
