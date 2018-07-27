#include "TMath.h"
#include <math.h>  


float deltaPhi(float phi1, float phi2) {

  auto dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2) {

  float dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}
