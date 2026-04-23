#pragma once

#include "OutputObjInterface.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <unordered_map>

// User-defined classes of output objects
#include "PullPtDepHistos.h"
#include "ResidPtDepHistos.h"
#include "TransportMethodCheck.h"

//#include "ExampleOutput.h"

struct OutputMap {
public:
  OutputMap(){};
  OutputMap(TFile &outFile);
  ~OutputMap(){};
  std::unordered_map<std::string, OutputObjInterface *> map;

  void CreateAll();
  void PreWriteProcessAll();
  void WriteAll();
  bool CheckAll();

private:
  // instances of user-defined classes of output objects
  ResidPtDepHistos *residPTOutput;
  PullPtDepHistos *pullPTOutput;
  TransportMethodCheck *transportCheckOutput;

  // ExampleOutput *exampleOutput;
};

//---------- add objects into map ------------
// here we add output objects in th map
OutputMap::OutputMap(TFile &outFile) {
  residPTOutput = new ResidPtDepHistos(outFile);
  map["residPTOutput"] = residPTOutput;

  pullPTOutput = new PullPtDepHistos(outFile);
  map["pullPTOutput"] = pullPTOutput;

  transportCheckOutput = new TransportMethodCheck(outFile);
  map["transportCheckOutput"] = transportCheckOutput;

  /*exampleOutput = new ExampleOutput(outFile);
  map["exampleOutput"] = exampleOutput;*/
}

//---------- methods realization ------------

void OutputMap::CreateAll() {
  for (auto &[name, obj] : map)
    obj->Create();
}

void OutputMap::PreWriteProcessAll() {
  for (auto &[name, obj] : map)
    obj->PreWriteProcess();
}

void OutputMap::WriteAll() {
  for (auto &[name, obj] : map)
    obj->Write();
}

bool OutputMap::CheckAll() {
  bool returnvalue = true;
  for (auto &[name, obj] : map) {
    bool value = obj->Check();
    if (!value)
      returnvalue = false;
  }
  return returnvalue;
}