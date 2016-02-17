#include <HistManager/EventHists.h>

EventHists::EventHists(TString name) {
   m_name  = "EventHist_" + name;
   m_label = "Event "; // don't forget the space
}

EventHists::~EventHists() {}

void EventHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  m_pstrk              = declare1D(m_name, "nPstrk",	        "N pseudoTracks in Event",	      8, -0.5, 7.5);  
  m_primary_pstrk      = declare1D(m_name, "nPrimaryPstrk",   "N Primary pseudoTracks in Event" ,   4, -0.5, 3.5);
  m_secondary_pstrk    = declare1D(m_name, "nSecondaryPstrk", "N Secondary pseudoTracks in Event",  8, -0.5, 7.5);
  m_lost_pstrk         = declare1D(m_name, "nLostPstrk",      "N Lost pseudoTracks in Event",	      4, -0.5, 3.5);
  m_recotrk            = declare1D(m_name, "nRecotrk",	      "N recoTracks in Event",  	      4, -0.5, 3.5);  
  m_recotrkwide        = declare1D(m_name, "nRecotrkwide",    "N recoTracks in Event",  	   200, -0.5,  199.5);  
  m_sispseed           = declare1D(m_name, "nSispseed",	      "N siSpSeeded Tracks in Event",       15, -0.5, 14.5);  
  m_nPixelClusters     = declare1D(m_name, "nPixelClusters",  "N Pixel Clusters in Event",          41, -0.5, 40.5);
  m_nSCTClusters       = declare1D(m_name, "nSCTClusters",    "N SCT Clusters in Event",            50, 0.0, 500);
  m_nSiClusters        = declare1D(m_name, "nSiClusters",     "N Si Clusters in Event",             55, 0.0, 550);
  m_nPixelClusterswide = declare1D(m_name, "nPixelClusterswide", "N Pixel Clusters in Event",     200, 0.0, 4000);
  m_nSCTClusterswide   = declare1D(m_name, "nSCTClusterswide",   "N SCT Clusters in Event",       300, 0.0, 6000);
  m_nSiClusterswide    = declare1D(m_name, "nSiClusterswide",    "N Si Clusters in Event",        200, 0.0, 8000);

  m_trigWeight       = declare1D(m_name, "trigWeight",  "trigger weight",10,-0.5,9.5);
  m_muValue          = declare1D(m_name, "muValue",     "pilueup  #mu",30,0.0,30.0);
  m_muValueCorr      = declare1D(m_name, "muValueCorr", "pilueup  #mu",30,0.0,30.0);
  m_LumiBlock        = declare1D(m_name, "LumiBlock",   "Luminosity block",2000,0,2000);

  // For ITk
  m_recotrk1GeV      = declare1D(m_name, "nRecotrk1GeV",    "N recoTracks (p_{T}>1GeV) in Event", 200, -0.5, 199.5);  
  m_average2TrackEta = declare1D(m_name, "average2TrackEta","average #eta of 2 tracks",      30, -3.0, 3.0);
  m_photonEta        = declare1D(m_name, "photonEta",       "photon #eta",                   5, 0.0, 3.0);
  m_photonEta2Track  = declare1D(m_name, "photonEta2Track", "photon #eta",                   5, 0.0, 3.0);

  m_photonProdR        = declare1D(m_name, "photonProdR",       "photon production r[mm]",     50,  0.0, 10.0);
  m_photonProdR2Track  = declare1D(m_name, "photonProdR2Track", "photon production r[mm]",     50,  0.0, 10.0);
  m_photonDecayR       = declare1D(m_name, "photonDecayR",       "photon conversion r[mm]",   55,  0.0, 1100.0);
  m_photonDecayR2Track = declare1D(m_name, "photonDecayR2Track", "photon conversion r[mm]",   55,  0.0, 1100.0);

  m_nTruth            = declare1D(m_name, "nTruth",            "N Truth in Event",  100, -0.5, 99.5);
  m_nChargedTruth     = declare1D(m_name, "nChargedTruth",     "N ChargedTruth in Event",  100, -0.5, 99.5);
  m_nNeutralTruth     = declare1D(m_name, "nNeutralTruth",     "N NeutralTruth in Event",  100, -0.5, 99.5);
  m_nTruth1GeV        = declare1D(m_name, "nTruth1GeV",        "N Truth (p_{T}>1GeV) in Event",  100, -0.5, 99.5);
  m_nChargedTruth1GeV = declare1D(m_name, "nChargedTruth1GeV", "N ChargedTruth (p_{T}>1GeV) in Event",  100, -0.5, 99.5);
  m_nNeutralTruth1GeV = declare1D(m_name, "nNeutralTruth1GeV", "N NeutralTruth (p_{T}>1GeV) in Event",  100, -0.5, 99.5);
  m_nIsoTrack         = declare1D(m_name, "nIsoTrack",         "N Iso.Track",      200, -0.5, 199.5);
  m_nVertex           = declare1D(m_name, "nVertex",           "N vertex",         40, -0.5, 39.5);

  m_effPileupTrack1GeV = declare1D(m_name, "effPileupTrack1GeV", "Efficiency of pileup track (p_{T}>1GeV)",  100, 0.0, 1.0);

  m_effPileupTrack1GeV_vsTruth = declare2D(m_name, "effPileupTrack1GeV_vsTruth", "Truth particles (p_{T}>1GeV)", "Tracks (p_{T}>1GeV)",   10, -0.5, 9.5, 10, -0.5, 9.5); 


  m_nPixelClustersAtIBL           = declare1D(m_name, "nPixelClustersAtIBL",      "N Pixel clusters at IBL", 100, 0.0, 500);
  m_nPixelHitsAtIBL               = declare1D(m_name, "nPixelHitsAtIBL",          "N Pixel hits at IBL",     100, 0.0, 2000);
  m_nPixelClustersAtIBL_vsMu      = declare2D(m_name, "nPixelClustersAtIBL_vsMu", "N Pixel clusters at IBL", "pilueup #mu", 100, 0.0, 500, 30,0.0,30.0); 
  m_nPixelHitsAtIBL_vsMu          = declare2D(m_name, "nPixelHitsAtIBL_vsMu",     "N Pixel hits at IBL", "pilueup #mu", 100, 0.0, 2000, 30,0.0,30.0); 
  m_nPixelClustersAtIBL_wide      = declare1D(m_name, "nPixelClustersAtIBL_wide",      "N Pixel clusters at IBL", 100, 0.0, 5000);
  m_nPixelHitsAtIBL_wide          = declare1D(m_name, "nPixelHitsAtIBL_wide",          "N Pixel hits at IBL",     100, 0.0, 20000);
  m_nPixelClustersAtIBL_vsMu_wide = declare2D(m_name, "nPixelClustersAtIBL_vsMu_wide", "N Pixel clusters at IBL", "pilueup #mu", 100, 0.0, 5000, 30,0.0,30.0); 
  m_nPixelHitsAtIBL_vsMu_wide     = declare2D(m_name, "nPixelHitsAtIBL_vsMu_wide",     "N Pixel hits at IBL", "pilueup #mu", 100, 0.0, 20000, 30,0.0,30.0); 

  m_nPixelClustersAtBLayer           = declare1D(m_name, "nPixelClustersAtBLayer",      "N Pixel clusters at BLayer", 100, 0.0, 500);
  m_nPixelHitsAtBLayer               = declare1D(m_name, "nPixelHitsAtBLayer",          "N Pixel hits at BLayer",     100, 0.0, 2000);
  m_nPixelClustersAtBLayer_vsMu      = declare2D(m_name, "nPixelClustersAtBLayer_vsMu", "N Pixel clusters at BLayer", "pilueup #mu", 100, 0.0, 500, 30,0.0,30.0); 
  m_nPixelHitsAtBLayer_vsMu          = declare2D(m_name, "nPixelHitsAtBLayer_vsMu",     "N Pixel hits at BLayer", "pilueup #mu", 100, 0.0, 2000, 30,0.0,30.0); 
  m_nPixelClustersAtBLayer_wide      = declare1D(m_name, "nPixelClustersAtBLayer_wide",      "N Pixel clusters at BLayer", 100, 0.0, 5000);
  m_nPixelHitsAtBLayer_wide          = declare1D(m_name, "nPixelHitsAtBLayer_wide",          "N Pixel hits at BLayer",     100, 0.0, 20000);
  m_nPixelClustersAtBLayer_vsMu_wide = declare2D(m_name, "nPixelClustersAtBLayer_vsMu_wide", "N Pixel clusters at BLayer", "pilueup #mu", 100, 0.0, 5000, 30,0.0,30.0); 
  m_nPixelHitsAtBLayer_vsMu_wide     = declare2D(m_name, "nPixelHitsAtBLayer_vsMu_wide",     "N Pixel hits at BLayer", "pilueup #mu", 100, 0.0, 20000, 30,0.0,30.0); 


} // BookHists

void EventHists::FillHists(EventFeatures &evt, float weight) const {
  m_pstrk              -> Fill(evt.nPseudoTracks,weight);
  m_primary_pstrk      -> Fill(evt.nPrimaryPseudoTracks,weight);
  m_secondary_pstrk    -> Fill(evt.nSecondaryPseudoTracks,weight);
  m_lost_pstrk         -> Fill(evt.nLostPseudoTracks,weight);
  m_recotrk            -> Fill(evt.nRecoTracks,weight);   
  m_recotrkwide        -> Fill(evt.nRecoTracks,weight);   
  m_sispseed           -> Fill(evt.nSiSpSeeds,weight);   
  m_nPixelClusters     -> Fill(evt.nPixelClusters,weight); 
  m_nSCTClusters       -> Fill(evt.nSCTClusters,weight); 
  m_nSiClusters        -> Fill(evt.nPixelClusters+evt.nSCTClusters,weight); 
  m_nPixelClusterswide -> Fill(evt.nPixelClusters,weight); 
  m_nSCTClusterswide   -> Fill(evt.nSCTClusters,weight); 
  m_nSiClusterswide    -> Fill(evt.nPixelClusters+evt.nSCTClusters,weight); 
  m_trigWeight         -> Fill(evt.trigCounter,weight);
  m_muValue            -> Fill(evt.muValue,weight);
  m_muValueCorr        -> Fill(evt.muValueCorr,weight);
  m_LumiBlock          -> Fill(evt.LumiBlock,weight);

  m_recotrk1GeV  -> Fill(evt.nRecoTracks1GeV,weight);   
  m_nIsoTrack    -> Fill(1.0*evt.nIsoTrack,weight);
  m_nVertex      -> Fill(1.0*evt.nVertex,weight);

  if (evt.nRecoTracks==2) {
    m_average2TrackEta -> Fill( evt.average2TrackEta,weight);
  }

  if (evt.photonEta>-50.0) {
    m_photonEta -> Fill(TMath::Abs(evt.photonEta),weight);
    if (TMath::Abs(evt.photonEta)<0.5) {
      m_photonProdR  -> Fill(evt.photonProdR,weight);
      m_photonDecayR -> Fill(evt.photonDecayR,weight);
    }
    if (evt.nRecoTracks==2) {
      m_photonEta2Track -> Fill(TMath::Abs(evt.photonEta),weight);
      if (TMath::Abs(evt.photonEta)<0.5) {
        m_photonProdR2Track  -> Fill(evt.photonProdR,weight);
        m_photonDecayR2Track -> Fill(evt.photonDecayR,weight);
      }
    }
  }

  m_nTruth        -> Fill(1.0*evt.nTruth,weight);
  m_nChargedTruth -> Fill(1.0*evt.nChargedTruth,weight);
  m_nNeutralTruth -> Fill(1.0*evt.nNeutralTruth,weight);
  m_nTruth1GeV        -> Fill(1.0*evt.nTruth1GeV,weight);
  m_nChargedTruth1GeV -> Fill(1.0*evt.nChargedTruth1GeV,weight);
  m_nNeutralTruth1GeV -> Fill(1.0*evt.nNeutralTruth1GeV,weight);

  if (evt.nChargedTruth1GeV>0 && evt.nRecoTracks1GeV>1) {
    m_effPileupTrack1GeV -> Fill(1.0*(evt.nRecoTracks1GeV-1.0)/evt.nChargedTruth1GeV,weight);
    m_effPileupTrack1GeV_vsTruth -> Fill(1.0*evt.nChargedTruth1GeV,1.0*evt.nRecoTracks1GeV-1.0,weight);
  }

  m_nPixelClustersAtIBL           -> Fill(1.0*evt.nPixelClustersAtIBL,weight);
  m_nPixelHitsAtIBL               -> Fill(1.0*evt.nPixelHitsAtIBL,weight);
  m_nPixelClustersAtIBL_vsMu      -> Fill(1.0*evt.nPixelClustersAtIBL,evt.muValue,weight);
  m_nPixelHitsAtIBL_vsMu          -> Fill(1.0*evt.nPixelHitsAtIBL,evt.muValue,weight);
  m_nPixelClustersAtIBL_wide      -> Fill(1.0*evt.nPixelClustersAtIBL,weight);
  m_nPixelHitsAtIBL_wide          -> Fill(1.0*evt.nPixelHitsAtIBL,weight);
  m_nPixelClustersAtIBL_vsMu_wide -> Fill(1.0*evt.nPixelClustersAtIBL,evt.muValue,weight);
  m_nPixelHitsAtIBL_vsMu_wide     -> Fill(1.0*evt.nPixelHitsAtIBL,evt.muValue,weight);

  m_nPixelClustersAtBLayer           -> Fill(1.0*evt.nPixelClustersAtBLayer,weight);
  m_nPixelHitsAtBLayer               -> Fill(1.0*evt.nPixelHitsAtBLayer,weight);
  m_nPixelClustersAtBLayer_vsMu      -> Fill(1.0*evt.nPixelClustersAtBLayer,evt.muValue,weight);
  m_nPixelHitsAtBLayer_vsMu          -> Fill(1.0*evt.nPixelHitsAtBLayer,evt.muValue,weight);
  m_nPixelClustersAtBLayer_wide      -> Fill(1.0*evt.nPixelClustersAtBLayer,weight);
  m_nPixelHitsAtBLayer_wide          -> Fill(1.0*evt.nPixelHitsAtBLayer,weight);
  m_nPixelClustersAtBLayer_vsMu_wide -> Fill(1.0*evt.nPixelClustersAtBLayer,evt.muValue,weight);
  m_nPixelHitsAtBLayer_vsMu_wide     -> Fill(1.0*evt.nPixelHitsAtBLayer,evt.muValue,weight);

}

