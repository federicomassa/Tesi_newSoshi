#include "ParticleAnalysis/DeltaFunctions.h"
#include "TMath.h"

float deltaPhi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if (dphi >= TMath::Pi()) { dphi = 2.0*TMath::Pi()-dphi; }
  if (dphi < -TMath::Pi()) { dphi = 2.0*TMath::Pi()+dphi; }
  return dphi;
}

float deltaR(float phi1, float phi2, float eta1, float eta2) {
  float dphi = phi1 - phi2;
  if (dphi >= TMath::Pi()) { dphi = 2.0*TMath::Pi()-dphi; }
  if (dphi < -TMath::Pi()) { dphi = 2.0*TMath::Pi()+dphi; }
  float dR = TMath::Sqrt(dphi*dphi+(eta1-eta2)*(eta1-eta2));
  return dR;
}
