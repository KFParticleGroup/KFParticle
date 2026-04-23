#pragma once

#include "../ParticleStructures/KFParticleStructure.h"
#include "../ParticleStructures/MCParticleStructure.h"
#include "TFile.h"

class OutputObjInterface {
public:
  OutputObjInterface(TFile &outputFile) { this->outputFile = &outputFile; };
  virtual ~OutputObjInterface(){};

  void virtual Create() = 0;
  void virtual Fill() = 0;
  void virtual Fill(int iDecay, const MCParticleStruct &mcPart) = 0;
  void virtual Fill(int iDecay, const KFParticleStruct &kfPart) = 0;
  void virtual Fill(int iDecay, const MCParticleStruct &mcPart,
                    const KFParticleStruct &kfPartDecay,
                    const KFParticleStruct &kfPartProduction) = 0;
  void virtual Fill(int iDecay, const MCParticleStruct &mcPart,
                    const KFParticleStruct &kfPart) = 0;
  void virtual PreWriteProcess() = 0;
  void virtual Write() = 0;
  bool virtual Check() { return true; }

protected:
  TFile *outputFile;
};