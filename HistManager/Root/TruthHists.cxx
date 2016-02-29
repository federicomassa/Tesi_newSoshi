/******************************************
 * Plot of tracks and tracking quantities
 *
 * G. Facini
 * Sep 10 10:03:25 CEST 2014
 *
 ******************************************/

#include <HistManager/TruthHists.h>
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"


TruthHists::TruthHists(TString name) {
   m_name = "TruthHist_" + name;
   m_label = "Truth Particle "; // don't forget the space
}

TruthHists::~TruthHists() {}

void TruthHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  /* particle properties */
  m_prodR       = declare1D(m_name, "prodR",       "production R [mm]", 110, 0.0,1100.0);
  m_prodRbarrel = declare1D(m_name, "prodRbarrel", "production R (|#eta(#gamma)|<1.5) [mm]", 110, 0.0,1100.0);
  m_prodZ       = declare1D(m_name, "prodZ",    "production Z [mm]", 310, 0.0,3100.0);
  m_decayR      = declare1D(m_name, "decayR",   "decay R [mm]",      100, 0.0, 1.0 );

  /* track parameterization */
  m_pt 			     = declare1D(m_name, "pt", 	 "p_{T} [GeV]", 	200,  0.0,     2000.0 ); 
  m_eta			     = declare1D(m_name, "eta",	 "#eta", 			  100,	 -5.0,     5.0    ); 
  m_abseta                   = declare1D(m_name, "abseta", "|#eta|",    25,   0.0,     5.0); 
  m_d0 			     = declare1D(m_name, "d0", 	 "d_{0}",		    120, -0.09,    0.09   ); 
  m_z0 			     = declare1D(m_name, "z0", 	 "z_{0}",		    120, -300.0,   300.0  ); 
  m_phi			     = declare1D(m_name, "phi",	 "#phi", 			  32,	 -3.2,     3.2    ); 

} // BookHists

void TruthHists::FillHists(const xAOD::TruthParticle* tpart, float weight) const {

  /* particle properties */
  if (tpart->hasProdVtx())  { 
    if (TMath::Abs(tpart->pdgId())==11) {

      bool isPhoton = false;
      float photonEta = 100.0;
      for (size_t iparent=0; iparent<tpart->prodVtx()->nIncomingParticles(); iparent++) {
        const xAOD::TruthParticle* parent = tpart->prodVtx()->incomingParticle(iparent);
        if (parent->pdgId()==22) { isPhoton=true; photonEta=parent->eta(); }
      }

      if (isPhoton) {
        m_prodR -> Fill(tpart->prodVtx()->perp(),weight); 
        m_prodZ -> Fill(TMath::Abs(tpart->prodVtx()->z()),weight); 
        if (TMath::Abs(photonEta)<1.5) {
          m_prodRbarrel -> Fill(tpart->prodVtx()->perp(),weight); 
        }
      }
    }
  }
  if (tpart->hasDecayVtx()) { m_decayR -> Fill(tpart->decayVtx()->perp(),weight); }

  /* track parameterization */
  m_pt 		 -> Fill(tpart->pt()* 1e-3,weight); 
  m_eta		 -> Fill(tpart->eta(),weight); 
  m_abseta -> Fill(TMath::Abs(tpart->eta()),weight); 
  if (tpart->isAvailable<float>("d0")) {
    m_d0 -> Fill(tpart->auxdata<float>("d0"),weight); 
  }
  if (tpart->isAvailable<float>("z0")) {
    m_z0 -> Fill(tpart->auxdata<float>("z0"),weight); 
  }
  m_phi -> Fill(tpart->phi(),weight); 

} // FillHists

