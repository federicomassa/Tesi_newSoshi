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

const double maxEta = 4.0;

TruthHists::TruthHists(TString name) {
  m_baseName = name;
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
  m_pt 			     = declare1D(m_name, "pt", 	 "p_{T} [GeV]", 	200,  0.0,     200.0 ); 
  m_ptVsEta 		     = declare2D(m_name, "ptVsEta", "p_{T} [GeV]", "|#eta|", 200,  0.0,     200.0, 40, 0, 4.0); 
  m_eta			     = declare1D(m_name, "eta",	 "#eta", 			  80,	 -maxEta,     maxEta   ); //da togliere per efficienza
  m_eta_wide		     = declare1D(m_name, "eta_wide",	 "#eta", 			  160,	 -8.0,     8.0    ); 
  m_abseta                   = declare1D(m_name, "abseta", "|#eta|",    20,   0.0,     maxEta); 
  m_d0 			     = declare1D(m_name, "d0", 	 "d_{0}",		    120, -0.09,    0.09   ); 
  m_z0 			     = declare1D(m_name, "z0", 	 "z_{0}",		    120, -300.0,   300.0  ); 
  m_phi			     = declare1D(m_name, "phi",	 "#phi", 			  32,	 -3.1415,     3.1415    ); 

  m_pdgId		     = declare1D(m_name, "pdgId",	 "pdgId", 			  5001,	 -2500,     2501    ); 

  m_prodVtx_X                = declare1D(m_name, "prodVtx_X", "Production vertex x [mm]", 100,-0.1, 0.1);
  m_prodVtx_Y                = declare1D(m_name, "prodVtx_Y", "Production vertex y [mm]", 100,-0.1, 0.1);
  m_prodVtx_Z                = declare1D(m_name, "prodVtx_Z", "Production vertex z [mm]", 4000,-200, 200);
  m_prodVtx_R                = declare1D(m_name, "prodVtx_R", "Production vertex r [mm]", 100,-0.1, 0.1);

  m_prodVtx_incoming         = declare1D(m_name, "prodVtx_incoming", "N of incoming particles in prod vertex", 100,-0.5, 99.5);
  m_prodVtx_outgoing         = declare1D(m_name, "prodVtx_outgoing", "N of outgoing particles in prod vertex", 100,-0.5, 99.5);

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
  m_ptVsEta      -> Fill(tpart->pt()* 1e-3,TMath::Abs(tpart->eta()), weight); 
  m_eta		 -> Fill(tpart->eta(),weight);
  m_eta_wide     -> Fill(tpart->eta(), weight);
  m_abseta -> Fill(TMath::Abs(tpart->eta()),weight); 
  if (tpart->isAvailable<float>("d0")) {
    m_d0 -> Fill(tpart->auxdata<float>("d0"),weight); 
  }
  // else if (tpart->hasProdVtx()) {
  //   const double x0 = tpart->prodVtx()->x();
  //   const double y0 = tpart->prodVtx()->y();
  //   m_d0 -> Fill(TMath::Sqrt(x0*x0 + y0*y0), weight);
  // }


  if (tpart->isAvailable<float>("z0")) {
    m_z0 -> Fill(tpart->auxdata<float>("z0"),weight); 
  }
  // else if (tpart->hasProdVtx()) {
  //   m_z0 -> Fill(tpart->prodVtx()->z(), weight);
  // }

  if (TMath::Abs(tpart->eta()) <= maxEta) //TODO: remember this. Only written if truth is detectable
    m_phi -> Fill(tpart->phi(),weight); 

  m_pdgId -> Fill(tpart->pdgId(), weight);

  //check prodVertex
  if (tpart->hasProdVtx()) {
    m_prodVtx_X -> Fill(tpart->prodVtx()->x(), weight);
    m_prodVtx_Y -> Fill(tpart->prodVtx()->y(), weight);
    m_prodVtx_Z -> Fill(tpart->prodVtx()->z(), weight);
    m_prodVtx_R -> Fill(TMath::Sqrt(TMath::Power(tpart->prodVtx()->x(),2) +
				    TMath::Power(tpart->prodVtx()->y(),2)),weight);
    m_prodVtx_incoming -> Fill(tpart->prodVtx()->nIncomingParticles(),weight);
    m_prodVtx_outgoing -> Fill(tpart->prodVtx()->nOutgoingParticles(),weight);
  }

		      
} // FillHists

