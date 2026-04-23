/*
An example of user-defined class of output object.
The short instruction:
1) Name your class
2) Add output objects (histograms, plots, etc.) into private filed of the class
3) Realize methods Create(), Fill(..), PreWriteProcess() and Write().
If some of the methods are not need - leave them empty;
4) Add your class into "OutputMap.h" file in the way it is done for this class
*/

#pragma once

#include "OutputObjInterface.h"

class ExampleOutput : public OutputObjInterface {
public:
  ExampleOutput(TFile &outputFile) : OutputObjInterface(outputFile){};
  ~ExampleOutput(){};

  // methods that must be realized by user
  void Create() override;
  void Fill() override{};
  void Fill(const MCParticleStruct &mcPart) override;
  void Fill(const KFParticleStruct &kfPart) override;
  void Fill(const MCParticleStruct &mcPart,
            const KFParticleStruct &kfPart) override;
  void PreWriteProcess() override;
  void Write() override;

private:
  TH1F *fHistX; // user-defined histograms
  TH1F *fHistY;
  TH1F *fHistZ;
};

////---------- method implementations ----------

void ExampleOutput::Create() {
  fHistX = new TH1F("fHistX", "fHistX", 100, -1, 2);
  fHistY = new TH1F("fHistY", "fHistY", 100, -1, 2);
  fHistZ = new TH1F("fHistZ", "fHistZ", 100, -1, 2);
}

void ExampleOutput::Fill(const MCParticleStruct &mcPart) {
  fHistX->Fill(mcPart.initialX);
  fHistY->Fill(mcPart.initialX);
  fHistZ->Fill(mcPart.initialX);
}

void ExampleOutput::Fill(const KFParticleStruct &kfPart) {
  //....
}

void ExampleOutput::Fill(const MCParticleStruct &mcPart,
                         const KFParticleStruct &kfPart) {
  //....
}

void ExampleOutput::PreWriteProcess() {}

void ExampleOutput::Write() {
  TDirectory *ExampleOutput = outputFile->mkdir("ExampleOutput");
  ExampleOutput->cd();
  fHistX->Write();
  fHistY->Write();
  fHistZ->Write();
}