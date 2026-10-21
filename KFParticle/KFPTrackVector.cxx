/*
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "KFPTrackVector.h"

#include <iostream>

void KFPTrackVector::Resize(const int n)
{
  /** Resizes all vectors in the class to a given value.
   ** \param[in] n - new size of the vector
   **/
  for (int i = 0; i < 6; i++) { fP[i].resize(n); }
  for (int i = 0; i < 21; i++) { fC[i].resize(n); }
#ifdef NonhomogeneousField
  for (int i = 0; i < 10; i++) { fField[i].resize(n); }
#endif
  //     fChi2.resize(n);
  //     fNDF.resize(n);
  fId.resize(n);
  fPDG.resize(n);
  fQ.resize(n);
  fPVIndex.resize(n);
  fNPixelHits.resize(n);
}

void KFPTrackVector::Set(KFPTrackVector& v, int vSize, int offset)
{
  /** Copies "vSize" tracks from the KFPTrackVector "v" to the current object.
   ** Tracks are put starting from the "offset" position.
   ** \param[in] v - external KFPTrackVector with input tracks to be copied
   ** \param[in] vSize - number of tracks to be copied from "v"
   ** \param[in] offset - offset position in the current object, starting from which input tracks will be stored
   **/
  for (int iV = 0; iV < vSize; iV++) {
    for (int i = 0; i < 6; i++) { fP[i][offset + iV] = v.fP[i][iV]; }
    for (int i = 0; i < 21; i++) { fC[i][offset + iV] = v.fC[i][iV]; }
#ifdef NonhomogeneousField
    for (int i = 0; i < 10; i++) { fField[i][offset + iV] = v.fField[i][iV]; }
#endif
    //       fChi2[offset+iV] = v.fChi2[iV];
    //       fNDF[offset+iV] = v.fNDF[iV];
    fId[offset + iV]         = v.fId[iV];
    fPDG[offset + iV]        = v.fPDG[iV];
    fQ[offset + iV]          = v.fQ[iV];
    fPVIndex[offset + iV]    = v.fPVIndex[iV];
    fNPixelHits[offset + iV] = v.fNPixelHits[iV];
  }
}

#if 0
void KFPTrackVector::SetTracks(const KFPTrackVector& track, const PaddedVector<int_v>& trackIndex, const int nIndexes)
{
  /** The current object is resised to "nIndexes", tracks with indices "trackIndex" are copied to the current object.
   ** \param[in] track - input vector of tracks
   ** \param[in] trackIndex - indices of tracks in a vector "track", which should be stored to the current object
   ** \param[in] nIndexes - number of tracks to be copied, defines the new size of the current object
   **/
  if(nIndexes == 0) return;
  
  Resize(nIndexes);

  for(int iP=0; iP<6; iP++)
  {
    int iElement = 0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      float_v& vec = reinterpret_cast<float_v&>(fP[iP][iElement]);
      vec.gather(&(track.fP[iP][0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    float_v& vec = reinterpret_cast<float_v&>(fP[iP][iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fP[iP][0]), correctedIndices);
    
  }
  for(int iC=0; iC<21; iC++)
  {
    int iElement=0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      float_v& vec = reinterpret_cast<float_v&>(fC[iC][iElement]);
      vec.gather(&(track.fC[iC][0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    float_v& vec = reinterpret_cast<float_v&>(fC[iC][iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fC[iC][0]), correctedIndices);
  }
#ifdef NonhomogeneousField
  for(int iP=0; iP<10; iP++)
  {
    int iElement = 0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      float_v& vec = reinterpret_cast<float_v&>(fField[iP][iElement]);
      vec.gather(&(track.fField[iP][0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    float_v& vec = reinterpret_cast<float_v&>(fField[iP][iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fField[iP][0]), correctedIndices); 
  }
#endif
  {
    int iElement=0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      int_v& vec = reinterpret_cast<int_v&>(fId[iElement]);
      vec.gather(&(track.fId[0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    int_v& vec = reinterpret_cast<int_v&>(fId[iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fId[0]), correctedIndices);
  }
  {
    int iElement=0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      int_v& vec = reinterpret_cast<int_v&>(fPDG[iElement]);
      vec.gather(&(track.fPDG[0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    int_v& vec = reinterpret_cast<int_v&>(fPDG[iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fPDG[0]), correctedIndices);
  }
  {
    int iElement=0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      int_v& vec = reinterpret_cast<int_v&>(fQ[iElement]);
      vec.gather(&(track.fQ[0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    int_v& vec = reinterpret_cast<int_v&>(fQ[iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fQ[0]), correctedIndices);
  }
  {
    int iElement=0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      int_v& vec = reinterpret_cast<int_v&>(fPVIndex[iElement]);
      vec.gather(&(track.fPVIndex[0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    int_v& vec = reinterpret_cast<int_v&>(fPVIndex[iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fPVIndex[0]), correctedIndices);
  }
  {
    int iElement=0;
    for(iElement=0; iElement<nIndexes-SimdLen; iElement += SimdLen)
    {
      const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
      int_v& vec = reinterpret_cast<int_v&>(fNPixelHits[iElement]);
      vec.gather(&(track.fNPixelHits[0]), index);
    }
    const int_v& index = reinterpret_cast<const int_v&>(trackIndex[iElement]);
    int_v& vec = reinterpret_cast<int_v&>(fNPixelHits[iElement]);
    const int_v correctedIndices = select(int_v::indicesSequence(iElement)<nIndexes, index, 0);
    vec.gather(&(track.fNPixelHits[0]), correctedIndices);
  }
}
#else
void KFPTrackVector::SetTracks(const KFPTrackVector& track, const PaddedVector<int_v>& trackIndex, const int nIndexes)
{
  /** The current object is resised to "nIndexes", tracks with indices "trackIndex" are copied to the current object.
   ** \param[in] track - input vector of tracks
   ** \param[in] trackIndex - indices of tracks in a vector "track", which should be stored to the current object
   ** \param[in] nIndexes - number of tracks to be copied, defines the new size of the current object
   **/
  if (nIndexes == 0) { return; }

  Resize(nIndexes);

  for (int iP = 0; iP < 6; iP++) {
    for (int iElement = 0; iElement < nIndexes; iElement++) { fP[iP][iElement] = track.fP[iP][trackIndex[iElement]]; }
  }
  for (int iC = 0; iC < 21; iC++) {
    for (int iElement = 0; iElement < nIndexes; iElement++) { fC[iC][iElement] = track.fC[iC][trackIndex[iElement]]; }
  }
#ifdef NonhomogeneousField
  for (int iF = 0; iF < 10; iF++) {
    for (int iElement = 0; iElement < nIndexes; iElement++) {
      fField[iF][iElement] = track.fField[iF][trackIndex[iElement]];
    }
  }
#endif
  for (int iElement = 0; iElement < nIndexes; iElement++) { fId[iElement] = track.fId[trackIndex[iElement]]; }
  for (int iElement = 0; iElement < nIndexes; iElement++) { fPDG[iElement] = track.fPDG[trackIndex[iElement]]; }
  for (int iElement = 0; iElement < nIndexes; iElement++) { fQ[iElement] = track.fQ[trackIndex[iElement]]; }
  for (int iElement = 0; iElement < nIndexes; iElement++) { fPVIndex[iElement] = track.fPVIndex[trackIndex[iElement]]; }
  for (int iElement = 0; iElement < nIndexes; iElement++) {
    fNPixelHits[iElement] = track.fNPixelHits[trackIndex[iElement]];
  }
}
#endif

void KFPTrackVector::RotateXY(float_v alpha, int firstElement)
{
  /** Rotates SIMD vector of tracks starting from the position "firstElement" onto the angles "alpha" in the XY plane.
   ** Rotation matrix is:
   \verbatim
   {  cos(A), -sin(A),  0,        0,        0,   0 }
   {  sin(A),  cos(A),  0,        0,        0,   0 }
   {       0,       0,  1,        0,        0,   0 }
   {       0,       0,  0,   cos(A),  -sin(A),   0 }
   {       0,       0,  0,   sin(A),   cos(A),   0 }
   {       0,       0,  0,        0,        0,   1 }
   \endverbatim
   ** \param[in] alpha - rotation angles
   ** \param[in] firstElement - track index, starting from which SIMD vector of tracks will be rotated
   **/

  assert(firstElement >= 0 && firstElement < Size() && firstElement % float_v::size() == 0);
  int iblock = firstElement / float_v::size();

  float_v sA = sin(alpha);
  float_v cA = cos(alpha);

  float_v& x = fP[0].block(iblock);
  float_v& y = fP[1].block(iblock);

  const float_v xInit = x;
  const float_v yInit = y;

  x = -(xInit * sA + yInit * cA);
  y = xInit * cA - yInit * sA;

  float_v& px = fP[3].block(iblock);
  float_v& py = fP[4].block(iblock);

  const float_v pxInit = px;
  const float_v pyInit = py;

  px = -(pxInit * sA + pyInit * cA);
  py = pxInit * cA - pyInit * sA;

  float_v cov[21];
  for (int iC = 0; iC < 21; iC++) { cov[iC] = fC[iC].block(iblock); }

  fC[0].block(iblock) = cA * cA * cov[2] + 2 * cA * cov[1] * sA + cov[0] * sA * sA;

  fC[1].block(iblock) = -(cA * cA * cov[1]) + cA * (-cov[0] + cov[2]) * sA + cov[1] * sA * sA;
  fC[2].block(iblock) = cA * cA * cov[0] - 2 * cA * cov[1] * sA + cov[2] * sA * sA;

  fC[3].block(iblock) = -(cA * cov[4]) - cov[3] * sA;
  fC[4].block(iblock) = cA * cov[3] - cov[4] * sA;
  //  fC[5].block(iblock) = cov[5];

  fC[6].block(iblock) = cA * cA * cov[11] + cA * (cov[10] + cov[7]) * sA + cov[6] * sA * sA;
  fC[7].block(iblock) = -(cA * cA * cov[10]) + cA * (cov[11] - cov[6]) * sA + cov[7] * sA * sA;
  fC[8].block(iblock) = -(cA * cov[12]) - cov[8] * sA;
  fC[9].block(iblock) = cA * cA * cov[14] + 2 * cA * cov[13] * sA + cov[9] * sA * sA;

  fC[10].block(iblock) = -(cA * cA * cov[7]) + cA * (cov[11] - cov[6]) * sA + cov[10] * sA * sA;
  fC[11].block(iblock) = cA * cA * cov[6] - cA * (cov[10] + cov[7]) * sA + cov[11] * sA * sA;
  fC[12].block(iblock) = cA * cov[8] - cov[12] * sA;
  fC[13].block(iblock) = -(cA * cA * cov[13]) + cA * (cov[14] - cov[9]) * sA + cov[13] * sA * sA;
  fC[14].block(iblock) = cA * cA * cov[9] - 2 * cA * cov[13] * sA + cov[14] * sA * sA;

  fC[15].block(iblock) = -(cA * cov[16]) - cov[15] * sA;
  fC[16].block(iblock) = cA * cov[15] - cov[16] * sA;
  //  fC[17].block(iblock) = cov[17];
  fC[18].block(iblock) = -(cA * cov[19]) - cov[18] * sA;
  fC[19].block(iblock) = cA * cov[18] - cov[19] * sA;
  //  fC[20].block(iblock) = cov[20];
}  // RotateXY


void KFPTrackVector::PrintTrack(int n)
{
  /** Prints parameters of the track with index "n".
   ** \param[in] n - index of track to be printed
   **/
  for (int i = 0; i < 6; i++) { std::cout << fP[i][n] << " "; }
  std::cout << std::endl;

  for (int i = 0; i < 21; i++) { std::cout << fC[i][n] << " "; }
  std::cout << std::endl;

  std::cout << fId[n] << " " << fPDG[n] << " " << fQ[n] << " " << fPVIndex[n] << " " << fNPixelHits[n] << std::endl;
}

void KFPTrackVector::Print()
{
  /** Prints all field of the current object. **/

  std::cout << "NTracks " << Size() << std::endl;
  if (Size() == 0) { return; }

  std::cout << "Parameters: " << std::endl;
  for (int iP = 0; iP < 6; iP++) {
    std::cout << "  iP " << iP << ": ";
    for (int iTr = 0; iTr < Size(); iTr++) { std::cout << Parameter(iP)[iTr] << " "; }
    std::cout << std::endl;
  }

  std::cout << "Cov matrix: " << std::endl;
  for (int iC = 0; iC < 21; iC++) {
    std::cout << "  iC " << iC << ": ";
    for (int iTr = 0; iTr < Size(); iTr++) { std::cout << Covariance(iC)[iTr] << " "; }
    std::cout << std::endl;
  }

  std::cout << "Id: " << std::endl;
  for (int iTr = 0; iTr < Size(); iTr++) { std::cout << Id()[iTr] << " "; }
  std::cout << std::endl;

  std::cout << "Pdg: " << std::endl;
  for (int iTr = 0; iTr < Size(); iTr++) { std::cout << PDG()[iTr] << " "; }
  std::cout << std::endl;

  std::cout << "Q: " << std::endl;
  for (int iTr = 0; iTr < Size(); iTr++) { std::cout << Q()[iTr] << " "; }
  std::cout << std::endl;

  std::cout << "PV index: " << std::endl;
  for (int iTr = 0; iTr < Size(); iTr++) { std::cout << PVIndex()[iTr] << " "; }
  std::cout << std::endl;

  std::cout << "fNPixelHits: " << std::endl;
  for (int iTr = 0; iTr < Size(); iTr++) { std::cout << NPixelHits()[iTr] << " "; }
  std::cout << std::endl;

#ifdef NonhomogeneousField
  std::cout << "Field: " << std::endl;
  for (int iF = 0; iF < 6; iF++) {
    std::cout << "  iF " << iF << ": ";
    for (int iTr = 0; iTr < Size(); iTr++) { std::cout << FieldCoefficient(iF)[iTr] << " "; }
    std::cout << std::endl;
  }
#endif

  std::cout << "Last particle index: " << std::endl;
  std::cout << LastElectron() << " " << LastMuon() << " " << LastPion() << " " << LastKaon() << " " << LastProton()
            << " " << LastDeuteron() << " " << LastTritium() << " " << LastHe3() << " " << LastHe4() << " " << LastHe6()
            << " " << LastLi6() << " " << LastLi7() << " " << LastBe7() << std::endl;
}
