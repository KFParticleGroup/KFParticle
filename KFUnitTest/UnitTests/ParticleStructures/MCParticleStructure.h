#pragma once

#include <cmath>

struct MCParticleStruct {
  MCParticleStruct(){};
  MCParticleStruct(const MCParticleStruct &part);

  int trackID;
  int parentID;
  int pdg;
  int charge;
  float mass;

  bool isMother;

  float initialX;
  float initialY;
  float initialZ;
  float finalX;
  float finalY;
  float finalZ;

  float trackLength;

  float initialPX;
  float initialPY;
  float initialPZ;
  float finalPX;
  float finalPY;
  float finalPZ;

  float ndf;
  float chi2;

  float initialPT() const {
    return sqrt(initialPX * initialPX + initialPY * initialPY);
  }
  float initialP() const {
    return sqrt(initialPX * initialPX + initialPY * initialPY +
                initialPZ * initialPZ);
  }
  float initialE() const {
    return sqrt(initialPX * initialPX + initialPY * initialPY +
                initialPZ * initialPZ + mass * mass);
  }
  float initialEta() const {
    return 0.5 * log((initialP() + initialPZ) / (initialP() - initialPZ));
  }
  float initialPhi() const { return atan2(initialPY, initialPX); }
  float initialTheta() const { return atan2(initialPT(), initialPZ); }

  float finalPT() const { return sqrt(finalPX * finalPX + finalPY * finalPY); }
  float finalP() const {
    return sqrt(finalPX * finalPX + finalPY * finalPY + finalPZ * finalPZ);
  }
  float finalE() const {
    return sqrt(finalPX * finalPX + finalPY * finalPY + finalPZ * finalPZ +
                mass * mass);
  }
  float finalEta() const {
    return 0.5 * log((finalP() + finalPZ) / (finalP() - finalPZ));
  }
  float finalPhi() const { return atan2(finalPY, finalPX); }
  float finalTheta() const { return atan2(finalPT(), finalPZ); }
};

MCParticleStruct::MCParticleStruct(const MCParticleStruct &part) {
  this->trackID = part.trackID;
  this->parentID = part.parentID;
  this->pdg = part.pdg;
  this->charge = part.charge;
  this->mass = part.mass;

  this->initialX = part.initialX;
  this->initialY = part.initialY;
  this->initialZ = part.initialZ;
  this->finalX = part.finalX;
  this->finalY = part.finalY;
  this->finalZ = part.finalZ;

  this->trackLength = part.trackLength;

  this->initialPX = part.initialPX;
  this->initialPY = part.initialPY;
  this->initialPZ = part.initialPZ;
  this->finalPX = part.finalPX;
  this->finalPY = part.finalPY;
  this->finalPZ = part.finalPZ;

  this->ndf = part.ndf;
  this->chi2 = part.chi2;
}