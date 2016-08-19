#include <HistManager/EventHists.h>
#include <Utility/TrackHelper.h>
#include <Utility/InvariantMass.h>

#include <algorithm>

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
  m_recotrk            = declare1D(m_name, "nRecotrk",	      "N recoTracks in Event",  	      30, 0, 30);  
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
  m_recotrk1GeV      = declare1D(m_name, "nRecotrk1GeV",    "N recoTracks (p_{T}>1GeV) in Event", 40, 0, 40);  
  m_average2TrackEta = declare1D(m_name, "average2TrackEta","average #eta of 2 tracks",      30, -3.0, 3.0);
  m_photonEta        = declare1D(m_name, "photonEta",       "photon #eta",                   5, 0.0, 3.0);
  m_photonEta2Track  = declare1D(m_name, "photonEta2Track", "photon #eta",                   5, 0.0, 3.0);

  m_photonProdR        = declare1D(m_name, "photonProdR",       "photon production r[mm]",     50,  0.0, 10.0);
  m_photonProdR2Track  = declare1D(m_name, "photonProdR2Track", "photon production r[mm]",     50,  0.0, 10.0);
  m_photonDecayR       = declare1D(m_name, "photonDecayR",       "photon conversion r[mm]",   55,  0.0, 1100.0);
  m_photonDecayR2Track = declare1D(m_name, "photonDecayR2Track", "photon conversion r[mm]",   55,  0.0, 1100.0);

  m_nTruth            = declare1D(m_name, "nTruth",            "N Truth in Event",  100, -0.5, 99.5);
  m_nChargedTruth     = declare1D(m_name, "nChargedTruth",     "N ChargedTruth in Event",  200, 0, 200);
  m_nPrimaryChargedTruth     = declare1D(m_name, "nPrimaryChargedTruth",     "N Primary Charged Truth in Event",  100, 0, 100);
  m_nPrimaryChargedTruth1GeV = declare1D(m_name, "nPrimaryChargedTruth1GeV",     "N Primary Charged Truth in Event with p_{T} > 1 GeV",  100, -0, 100);
  m_nSecondaryChargedTruth     = declare1D(m_name, "nSecondaryChargedTruth",     "N Secondary Charged Truth in Event",  100, 0, 100);
  m_nSecondaryChargedTruth1GeV = declare1D(m_name, "nSecondaryChargedTruth1GeV",     "N Secondary Charged Truth in Event with p_{T} > 1 GeV",  100, -0, 100);

  m_nNeutralTruth     = declare1D(m_name, "nNeutralTruth",     "N NeutralTruth in Event",  100, -0.5, 99.5);
  m_nTruth1GeV        = declare1D(m_name, "nTruth1GeV",        "N Truth (p_{T}>1GeV) in Event",  100, -0.5, 99.5);
  m_nChargedTruth1GeV = declare1D(m_name, "nChargedTruth1GeV", "N ChargedTruth (p_{T}>1GeV) in Event",  100, 0, 100);
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

  m_truthMass                        = declare1D(m_name, "truthMass", "Truth parent mass [GeV]", 300,0.0,300.0);
  m_truthMassWithPhotons             = declare1D(m_name, "truthMassWithPhotons", "Truth parent mass [GeV]", 200,0,200);
  m_nPrimaryPhotons                  = declare1D(m_name, "nPrimaryPhotons", "n primary photons in event", 10, 0, 10);
  m_truthEta                         = declare1D(m_name, "truthEta", "Truth parent #eta", 100, 0,10.0);
  m_truthPhi                         = declare1D(m_name, "truthPhi", "Truth parent #eta", 100, -3.1416,3.1416);

  m_truthPtWithPhotons               = declare1D(m_name, "truthPtWithPhotons", "Truth parent #eta", 300, 0,300.0);
  m_truthEtaWithPhotons              = declare1D(m_name, "truthEtaWithPhotons", "Truth parent #eta", 100, 0,10.0);
  m_truthPhiWithPhotons              = declare1D(m_name, "truthPhiWithPhotons", "Truth parent #eta", 100, -3.1416,3.1416);


  m_visibleTruthMass                 = declare1D(m_name, "visibleTruthMass", "Truth parent mass [GeV]", 300,0.0,300.0);
  m_truthOutsideDetector             = declare1D(m_name, "truthOutsideDetector", "# hard truth outside detector", 10, 0.0, 10.0);
  m_truthNCenteredMuons              = declare1D(m_name, "truthNCenteredMuons", "N muons within 2.7", 5, 0.0, 5.0);
  m_truthPt                          = declare1D(m_name, "truthPt", "Truth parent pt [GeV]", 300,0.0,300.0);
  m_truthPz                          = declare1D(m_name, "truthPz", "Truth parent pz [GeV]", 400,-2000.0,2000.0);
  m_truthMaxPt                       = declare1D(m_name, "truthMaxPt", "product max Pt [GeV]", 100, 0.0, 300.0);
  m_truthMinPt                       = declare1D(m_name, "truthMinPt", "product min Pt [GeV]", 100, 0.0, 300.0);
  m_truthMaxEta                      = declare1D(m_name, "truthMaxEta", "product max |#eta|", 100, 0.0, 10.0);
  m_truthMinEta                      = declare1D(m_name, "truthMinEta", "product min |#eta|", 100, 0.0, 10.0);
  m_truthCharge                      = declare1D(m_name, "truthCharge", "Truth parent charge", 20,-10.0,10.0);
  m_truthSmallestDR                  = declare1D(m_name, "truthSmallestDR", "minimum #Delta R", 400, 0.0, 20.0);
  m_truthDZ                          = declare1D(m_name, "truthDZ", "#Delta z between products [mm]", 1000, -1000.0, 1000.0);

  //pt ranking, 1 is the largest
  m_truthPt1                         = declare1D(m_name, "truthPt1", "#mu p_{T, 1} [GeV]", 200, 0.0,200.0);   
  m_truthPt2                         = declare1D(m_name, "truthPt2", "#mu p_{T, 2} [GeV]", 200, 0.0,200.0);   
  m_truthPt3                         = declare1D(m_name, "truthPt3", "#mu p_{T, 3} [GeV]", 200, 0.0,200.0);   
  m_truthPt4                         = declare1D(m_name, "truthPt4", "#mu p_{T, 4} [GeV]", 200, 0.0,200.0);   

  m_isTruthPt1OnShell                = declare1D(m_name, "isTruthPt1OnShell", "is on shell", 2, 0.0, 2.0);
  m_isTruthPt2OnShell                = declare1D(m_name, "isTruthPt2OnShell", "is on shell", 2, 0.0, 2.0);
  m_isTruthPt3OnShell                = declare1D(m_name, "isTruthPt3OnShell", "is on shell", 2, 0.0, 2.0);
  m_isTruthPt4OnShell                = declare1D(m_name, "isTruthPt4OnShell", "is on shell", 2, 0.0, 2.0);
  m_isTruthPtRankingOrdered          = declare1D(m_name, "isTruthPtRankingOrdered", "is ranking ordered", 2, 0.0, 2.0);

  m_truthCutCode                     = declare1D(m_name, "truthCutCode", "truth cut code", 256, 0.0, 256);

  m_truthEta1VsEta2                  = declare2D(m_name, "truthEta1VsEta2", "#eta_1", "#eta_2", 200,-10.0,10.0,200,-10.0,10.0);

  m_visibleTruth_abseta              = declare1D(m_name, "visibleTruth_abseta", "truth |#eta|", 20, 0.0, 4.0);
  m_reco_abseta                      = declare1D(m_name, "reco_abseta", "reco |#eta|", 20,0.0,4.0);

  m_matchedMass                      = declare1D(m_name, "matchedMass", "Matched parent mass [GeV]", 300,0.0,300.0);
  m_matchedMaxPt                     = declare1D(m_name, "matchedMaxPt", "product max Pt [GeV]", 100, 0.0, 300.0);
  m_matchedMinPt                     = declare1D(m_name, "matchedMinPt", "product min Pt [GeV]", 100, 0.0, 300.0);
  m_matchedMaxEta                    = declare1D(m_name, "matchedMaxEta", "product max |#eta|", 100, 0.0, 10.0);
  m_matchedMinEta                    = declare1D(m_name, "matchedMinEta", "product min |#eta|", 100, 0.0, 10.0);
  m_matchedCharge                    = declare1D(m_name, "matchedCharge", "Matched parent charge", 20,-10.0,10.0);
  m_matchedSmallestDR                = declare1D(m_name, "matchedSmallestDR", "minimum #Delta R", 400, 0.0, 20.0);
  m_matchedTrueMaxEta                = declare1D(m_name, "matchedTrueMaxEta", "product max |#eta|", 100, 0.0,10.0);
  m_matchedDZ                        = declare1D(m_name, "matchedDZ", "#Delta z between products [mm]", 1000, -1000.0, 1000.0);

  m_isMatchingRight                  = declare1D(m_name, "isMatchingRight", "", 5, 0.0, 5.0);

  m_recoMass                         = declare1D(m_name, "recoMass", "Reco H mass [GeV]", 300,0.0,300.0);
  m_recoMass27                       = declare1D(m_name, "recoMass27", "Reco H mass [GeV]", 300,0.0,300.0);
  m_recoMass32                       = declare1D(m_name, "recoMass32", "Reco H mass [GeV]", 300,0.0,300.0);
  m_recoMass4                        = declare1D(m_name, "recoMass4", "Reco H mass [GeV]", 300,0.0,300.0);
  m_recoOnShellMass                  = declare1D(m_name, "recoOnShellMass", "Reco on-shell #mu#mu mass [GeV]", 200,0.0,400.0);
  m_recoOnShellEta1VsEta2            = declare2D(m_name, "recoOnShellEta1VsEta2", "Reco on-shell #eta_{1}", "Reco on-shell #eta_{2}", 40,-4.0,4.0, 40, -4.0, 4.0);
  m_recoOnVsOffShellMass             = declare2D(m_name, "recoOnVsOffShellMass", "Reco on-shell #mu#mu mass [GeV]","Reco off-shell #mu#mu mass [GeV]", 200,0.0,200.0, 200, 0.0, 200.0);
  m_recoOnShellMass27                = declare1D(m_name, "recoOnShellMass27", "Reco on-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoOnShellEta1VsEta2_27         = declare2D(m_name, "recoOnShellEta1VsEta2_27", "Reco on-shell #eta_{1}", "Reco on-shell #eta_{2}", 40,-4.0,4.0, 40, -4.0, 4.0);
  m_recoOnShellMass32                = declare1D(m_name, "recoOnShellMass32", "Reco on-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoOnShellEta1VsEta2_32         = declare2D(m_name, "recoOnShellEta1VsEta2_32", "Reco on-shell #eta_{1}", "Reco on-shell #eta_{2}", 40,-4.0,4.0, 40, -4.0, 4.0);
  m_recoOnShellMass4                 = declare1D(m_name, "recoOnShellMass4", "Reco on-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoOnShellEta1VsEta2_4          = declare2D(m_name, "recoOnShellEta1VsEta2_4", "Reco on-shell #eta_{1}", "Reco on-shell #eta_{2}", 40,-4.0,4.0, 40, -4.0, 4.0);
  m_recoOffShellMass                 = declare1D(m_name, "recoOffShellMass", "Reco off-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoOffShellMass27               = declare1D(m_name, "recoOffShellMass27", "Reco off-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoOffShellMass32               = declare1D(m_name, "recoOffShellMass32", "Reco off-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoOffShellMass4                = declare1D(m_name, "recoOffShellMass4", "Reco off-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoVisibleMass                  = declare1D(m_name, "recoVisibleMass", "Reco parent mass [GeV]", 300,0.0,300.0);
  m_recoVisibleOnShellMass           = declare1D(m_name, "recoVisibleOnShellMass", "Reco on-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoVisibleOffShellMass          = declare1D(m_name, "recoVisibleOffShellMass", "Reco off-shell #mu#mu mass [GeV]", 200,0.0,200.0);
  m_recoMassVsRecoOffShellMass      = declare2D(m_name, "recoMassVsRecoOffShellMass", "4#mu mass", "off-shell #mu#mu mass", 
						 300, 0, 300,
						 200, 0, 200);

  m_recoCharge                       = declare1D(m_name, "recoCharge", "Reco parent charge", 20,-10.0,10.0);
  m_recoTruthCharge                  = declare1D(m_name, "recoTruthCharge", "Reco parent truth charge", 20,-10.0,10.0);
  m_recoSmallestDR                   = declare1D(m_name, "recoSmallestDR", "minimum #Delta R", 400, 0.0, 20.0);

  m_recoPt1                          = declare1D(m_name, "recoPt1", "Pt [GeV]", 300, 0.0, 150.0);

  m_recoPt2                          = declare1D(m_name, "recoPt2", "Pt [GeV]", 300, 0.0, 150.0);

  m_recoPt3                          = declare1D(m_name, "recoPt3", "Pt [GeV]", 300, 0.0, 150.0);

  m_recoPt4                          = declare1D(m_name, "recoPt4", "Pt [GeV]", 300, 0.0, 150.0);

  m_recoTruthPt                      = declare1D(m_name, "recoTruthPt", "parent Pt [GeV]", 100, 0.0, 20.0);
  m_recoMaxPt                        = declare1D(m_name, "recoMaxPt", "product max Pt [GeV]", 100, 0.0, 300.0);
  m_recoMinPt                        = declare1D(m_name, "recoMinPt", "product min Pt [GeV]", 100, 0.0, 300.0);
  m_recoMaxEta                       = declare1D(m_name, "recoMaxEta", "product max |#eta|", 100, 0.0, 10.0);
  m_recoMinEta                       = declare1D(m_name, "recoMinEta", "product min |#eta|", 100, 0.0, 10.0);
  m_recoDZ                           = declare1D(m_name, "recoDZ", "#Delta z between products [mm]", 1000, -10.0, 10.0);
  m_recoMinZVsMaxZ                   = declare2D(m_name, "recoMinZVsMaxZ", "reco min z0 [mm]", "reco max z0 [mm]", 1000, -1000.0,1000.0,
						 1000, -1000.0,1000.0);
  m_recoNCenteredMuons               = declare1D(m_name, "recoNCenteredMuons", "N muons with |#eta| < 2.7", 5,0,5);
  m_recoNCenteredOffShellMuons       = declare1D(m_name, "recoNCenteredOffShellMuons", "N muons with |#eta| < 2.7", 5,0,5);
  
  //reco matching
  m_recoMatchedMaxDR                 = declare1D(m_name, "recoMatchedMaxDR", "max dR in reco matching", 2000, 0.0, 20.0);
  m_isRecoMatchingRight              = declare1D(m_name, "isRecoMatchingRight", "", 5, 0.0, 5.0);
  m_isRecoSetOutsideAcceptance       = declare1D(m_name, "isRecoSetOutsideAcceptance", "", 2, 0.0, 2.0);

  m_recoIsolation                    = declare1D(m_name, "recoMaxIsolation", "isolation", 100, 0, 5);
  m_recoIsolationVsPt                = declare2D(m_name, "recoIsolationVsPt", "isolation", "track p_{T}", 100, 0, 5, 200, 0, 200);
  m_recoSumPtVsPt                    = declare2D(m_name, "recoSumPtVsPt", "sum of p_{T}", "track p_{T}", 200, 0, 200, 200, 0, 200);



  m_recoCutCode                      = declare1D(m_name, "recoCutCode", "reco cut code", 51, -1.0, 50);
  m_recoCutCode27                    = declare1D(m_name, "recoCutCode27", "reco cut code", 51, -1.0, 50);
  m_recoCutCode32                    = declare1D(m_name, "recoCutCode32", "reco cut code", 51, -1.0, 50);
  m_recoCutCode4                     = declare1D(m_name, "recoCutCode4", "reco cut code", 51, -1.0, 50);

  m_recoPt                           = declare1D(m_name, "recoPt", "4#mu Pt [GeV]", 300, 0.0, 300.0);
  m_recoEta                          = declare1D(m_name, "recoEta", "4#mu |#eta|", 100, 0.0, 10);
  m_recoPhi                          = declare1D(m_name, "recoPhi", "4#mu #phi [rad]", 100, -3.1416, 3.1416);


  // m_biasMatchedMass                  = declare1D(m_name, "biasMatchedMass", "Bias parent mass [GeV]", 200,-100,100.0);
  // m_biasRecoMass                     = declare1D(m_name, "biasRecoMass", "Bias parent mass [GeV]", 200,-100,100.0);
  // m_biasMatchedVsRecoMass            = declare1D(m_name, "biasMatchedVsRecoMass", " - Matched parent mass [GeV]", 200,-100,100.0);


  //fake charge hists
  m_matchedChargeVsTruthMaxEta         = declare2D(m_name, "matchedChargeVsTruthMaxEta", "truth max |#eta|", "total matched charge", 50, 0.0, 10.0, 8, -4,4); 
  m_fakeChargeTruthMinEtaVsTruthMaxEta = declare2D(m_name, "fakeChargeTruthMinEtaVsTruthMaxEta", "truth min |#eta|", "truth max |#eta|", 50,0.0,10.0,50.0,0.0,10.0);
  m_fakeChargeTruthMinPtVsTruthMaxPt   = declare2D(m_name, "fakeChargeTruthMinPtVsTruthMaxPt", "truth min pt [GeV]", "truth max pt [GeV]", 30,0.0,300.0,30,0.0,300.0);
  m_fakeChargeBiasMatchedMass          = declare1D(m_name, "fakeChargeBiasMatchedMass", "bias Matched Mass [GeV]", 40,-200.0,200.0);
  m_fakeChargeBiasMatchedMaxEtaPtVsTruthMaxEta = declare2D(m_name,
						     "fakeChargeBiasMatchedMaxEtaPtVsTruthMaxEta",
						     "bias matched Pt [GeV]",						     "truth max |#eta|",
						     20,-100.0,100.0,
						     500,0.0,10.0);

  m_fakeChargeBiasMatchedMaxEtaPt = declare1D(m_name,
					      "fakeChargeBiasMatchedMaxEtaPt",
					      "bias matched Pt [GeV]",
					      20,-100.0,100.0);

  m_fakeChargePt                       = declare1D(m_name, "fakeChargePt", "p_{T} of fake charge track [GeV]", 30,0.0,300.0);
  m_fakeChargeEta                      = declare1D(m_name, "fakeChargeEta", "#eta of fake charge track", 40, -4.0,4.0);

  m_fakeChargeMaxEtaBiasCharge = declare1D(m_name,
					   "fakeChargeMaxEtaBiasCharge",
					   "bias charge",
					   8, -4.0,4.0);

  m_NFakeChargeMatchedTracks               = declare1D(m_name, "NFakeChargeMatchedTracks", "number of fake charge tracks", 5, 0.0, 5.0);
  m_NFakeChargeRecoTracks                  = declare1D(m_name, "NFakeChargeRecoTracks", "number of fake charge tracks", 5, 0.0, 5.0);
  m_NFakeChargeRecoVsDR                    = declare2D(m_name, "NFakeChargeRecoVsDR", "number of fake charge tracks", "#DeltaR", 5, 0.0, 5.0, 50, 0.0,5.0);

  //Higgs analysis

  m_truthOnShellMass                       = declare1D(m_name, "truthOnShellMass", "on-shell mass [GeV]", 150, 0.0,150.0);
  m_truthOffShellMass                      = declare1D(m_name, "truthOffShellMass", "off-shell mass [GeV]", 150, 0.0,150.0);
  m_visibleTruthOnShellMass                = declare1D(m_name, "truthVisibleOnShellMass", "on-shell mass [GeV]", 150, 0.0,150.0);
  m_visibleTruthOffShellMass               = declare1D(m_name, "truthVisibleOffShellMass", "off-shell mass [GeV]", 150, 0.0,150.0);

  m_truthMassVsTruthOffShellMass           = declare2D(m_name, "truthMassVsTruthOffShellMass", "4#mu mass", "off-shell #mu#mu mass", 
						       300, 0, 300,
						       200, 0, 200);


  // resolutions

  m_sigRecoMass = declare1D(m_name, "sigRecoMass", "reco - truth", 200, -50,50.0);
  m_sigRecoMass27 = declare1D(m_name, "sigRecoMass27", "reco - truth", 200, -50,50.0);
  m_sigRecoMass32 = declare1D(m_name, "sigRecoMass32", "reco - truth", 200, -50,50.0);
  m_sigRecoMass4 = declare1D(m_name, "sigRecoMass4", "reco - truth", 200, -50,50.0);

  m_sigRecoOnShellMass = declare1D(m_name, "sigRecoOnShellMass", "reco - truth", 200, -50,50.0);
  m_sigRecoOnShellMass27 = declare1D(m_name, "sigRecoOnShellMass27", "reco - truth", 200, -50,50.0);
  m_sigRecoOnShellMass32 = declare1D(m_name, "sigRecoOnShellMass32", "reco - truth", 200, -50,50.0);
  m_sigRecoOnShellMass4 = declare1D(m_name, "sigRecoOnShellMass4", "reco - truth", 200, -50,50.0);

  m_sigRecoOffShellMass = declare1D(m_name, "sigRecoOffShellMass", "reco - truth", 200, -50,50.0);
  m_sigRecoOffShellMass27 = declare1D(m_name, "sigRecoOffShellMass27", "reco - truth", 200, -50,50.0);
  m_sigRecoOffShellMass32 = declare1D(m_name, "sigRecoOffShellMass32", "reco - truth", 200, -50,50.0);
  m_sigRecoOffShellMass4 = declare1D(m_name, "sigRecoOffShellMass4", "reco - truth", 200, -50,50.0);

  m_sigPt = declare1D(m_name, "sigRecoPt", "reco - truth", 200, -50,50.0);
  m_sigPt27 = declare1D(m_name, "sigRecoPt27", "reco - truth", 200, -50,50.0);
  m_sigPt32 = declare1D(m_name, "sigRecoPt32", "reco - truth", 200, -50,50.0);
  m_sigPt4 = declare1D(m_name, "sigRecoPt4", "reco - truth", 200, -50,50.0);

  m_sigEta = declare1D(m_name, "sigRecoEta", "reco - truth", 200, -2,2.0);
  m_sigEta27 = declare1D(m_name, "sigRecoEta27", "reco - truth", 200, -2,2.0);
  m_sigEta32 = declare1D(m_name, "sigRecoEta32", "reco - truth", 200, -2,2.0);
  m_sigEta4 = declare1D(m_name, "sigRecoEta4", "reco - truth", 200, -2,2.0);

  m_sigPhi = declare1D(m_name, "sigRecoPhi", "reco - truth", 100, -0.2,0.2);
  m_sigPhi27 = declare1D(m_name, "sigRecoPhi27", "reco - truth", 100, -0.2,0.2);
  m_sigPhi32 = declare1D(m_name, "sigRecoPhi32", "reco - truth", 100, -0.2,0.2);
  m_sigPhi4 = declare1D(m_name, "sigRecoPhi4", "reco - truth", 100, -0.2,0.2);







} // BookHists

void EventHists::FillHists(EventFeatures &evt, float weight) const {

  double truthOnShellMass, truthOffShellMass;


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
  m_nPrimaryChargedTruth -> Fill(1.0*evt.nPrimaryChargedTruth,weight);
  m_nPrimaryChargedTruth1GeV -> Fill(1.0*evt.nPrimaryChargedTruth1GeV,weight);
  m_nSecondaryChargedTruth -> Fill(1.0*evt.nSecondaryChargedTruth,weight);
  m_nSecondaryChargedTruth1GeV -> Fill(1.0*evt.nSecondaryChargedTruth1GeV,weight);

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


  //physics hists
  if (evt.physicsEvent.isTruthSet()) {

    std::vector<xAOD::TruthParticleContainer::const_iterator> truth_v = evt.physicsEvent.GetTruth();      

    m_truthMass                     -> Fill(1.0*evt.physicsEvent.GetTruthMass()/1000.0, weight);


    m_recoCutCode           -> Fill(evt.physicsEvent.recoCutCode, weight);

    if (TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()) < 2.7)
      m_recoCutCode27            -> Fill(evt.physicsEvent.recoCutCode27, weight);
    else if (TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()) < 3.2)
      m_recoCutCode32            -> Fill(evt.physicsEvent.recoCutCode32, weight);
    else if (TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()) < 4.0)
      m_recoCutCode4            -> Fill(evt.physicsEvent.recoCutCode4, weight);
    

    if (evt.physicsEvent.arePhotonsSet()) {
      m_truthMassWithPhotons        -> Fill(1.0*evt.physicsEvent.GetTruthP4WithPhotons().M()/1000.0, weight);
      m_nPrimaryPhotons             -> Fill(evt.physicsEvent.GetPhotons().size());

    }

    if (evt.physicsEvent.isTruthVisible)
      m_visibleTruthMass                     -> Fill(1.0*evt.physicsEvent.GetTruthMass()/1000.0, weight);

    int NCenteredMuons = 0;
    for (auto itr = truth_v.begin();
	 itr != truth_v.end();
	 itr++) {
      if (TMath::Abs((**itr)->eta()) < 2.7)
	NCenteredMuons++;
    }

    
    m_truthNCenteredMuons -> Fill(NCenteredMuons, weight);
    
    
    m_truthOutsideDetector          -> Fill(1.0*evt.physicsEvent.GetTruthOutsideDetector(), weight);

    m_truthPt                       -> Fill(1.0*evt.physicsEvent.GetTruthP4().Pt()/1000.0, weight);
    m_truthEta                      -> Fill(evt.physicsEvent.GetTruthP4().Eta(), weight);
    m_truthPhi            -> Fill(1.0*evt.physicsEvent.GetTruthP4().Phi(), weight);


    if (evt.physicsEvent.isRecoSet()) {
      m_truthPtWithPhotons            -> Fill(1.0*evt.physicsEvent.GetTruthP4WithPhotons().Pt()/1000.0, weight);
      m_truthEtaWithPhotons            -> Fill(TMath::Abs(1.0*evt.physicsEvent.GetTruthP4WithPhotons().Eta()), weight);
      m_truthPhiWithPhotons  -> Fill(1.0*evt.physicsEvent.GetTruthP4WithPhotons().Phi());
    }

    m_truthPz                       -> Fill(1.0*evt.physicsEvent.GetTruthP4().Pz()/1000.0, weight);
    m_truthMaxPt                    -> Fill(1.0*evt.physicsEvent.GetTruthMaxPt()/1000.0, weight);
    m_truthMaxEta                   -> Fill(TMath::Abs(1.0*evt.physicsEvent.GetTruthWithMaxEta()->eta()), weight);
    m_truthCharge                   -> Fill(1.0*evt.physicsEvent.GetTruthCharge(), weight);
    m_truthSmallestDR               -> Fill(1.0*evt.physicsEvent.GetTruthSmallestDR(), weight);
    
    std::vector<double> rankedTruthPt;
    for (auto itr = truth_v.begin();
	 itr != truth_v.end();
	 itr++) {
      rankedTruthPt.push_back((**itr)->pt());
    }
    
    std::sort(rankedTruthPt.begin(), rankedTruthPt.end());
    
    if (truth_v.size() == 4) {
      m_truthPt1                      -> Fill(rankedTruthPt[0]/1000.0);
      m_truthPt2                      -> Fill(rankedTruthPt[1]/1000.0);
      m_truthPt3                      -> Fill(rankedTruthPt[2]/1000.0);
      m_truthPt4                      -> Fill(rankedTruthPt[3]/1000.0);

      bool isTruthPt1OnShell = (TMath::Abs((*(truth_v[0]))->pt() - rankedTruthPt[0]) < 1E-6 ||
				TMath::Abs((*(truth_v[1]))->pt() - rankedTruthPt[0]) < 1E-6);
      bool isTruthPt2OnShell = (TMath::Abs((*(truth_v[0]))->pt() - rankedTruthPt[1]) < 1E-6 ||
                                TMath::Abs((*(truth_v[1]))->pt() - rankedTruthPt[1]) < 1E-6);
      bool isTruthPt3OnShell = (TMath::Abs((*(truth_v[0]))->pt() - rankedTruthPt[2]) < 1E-6 ||
                                TMath::Abs((*(truth_v[1]))->pt() - rankedTruthPt[2]) < 1E-6);
      bool isTruthPt4OnShell = (TMath::Abs((*(truth_v[0]))->pt() - rankedTruthPt[3]) < 1E-6 ||
                                TMath::Abs((*(truth_v[1]))->pt() - rankedTruthPt[3]) < 1E-6);


      m_isTruthPt1OnShell             -> Fill(isTruthPt1OnShell, weight);
      m_isTruthPt2OnShell             -> Fill(isTruthPt2OnShell, weight);
      m_isTruthPt3OnShell             -> Fill(isTruthPt3OnShell, weight);
      m_isTruthPt4OnShell             -> Fill(isTruthPt4OnShell, weight);

      m_isTruthPtRankingOrdered       -> Fill(!isTruthPt1OnShell && !isTruthPt2OnShell && isTruthPt3OnShell && isTruthPt4OnShell, weight);

      m_truthCutCode                  -> Fill(evt.physicsEvent.truthCutCode, weight);


    }
    
    
    if (evt.physicsEvent.GetTruthOutsideDetector() == 0) {
      for (auto itr = truth_v.begin();
	   itr != truth_v.end();
	   itr++) {
	m_visibleTruth_abseta->Fill(TMath::Abs((**itr)->eta()), weight);
      }
    }
    
    m_truthEta1VsEta2               -> Fill((*(truth_v[0]))->eta(), (*(truth_v[1]))->eta(), weight);
    
    
    if (evt.physicsEvent.GetTruthMaxZ()->hasProdVtx() && evt.physicsEvent.GetTruthMinZ()->hasProdVtx())
      m_truthDZ                       -> Fill(evt.physicsEvent.GetTruthMaxZ()->prodVtx()->z() - evt.physicsEvent.GetTruthMinZ()->prodVtx()->z(), weight);
    
    //Higgs Analysis
    if (truth_v.size() == 4) {
      TLorentzVector tmp(0.0,0.0,0.0,0.0);
      for (auto itr = truth_v.begin();
	   itr != truth_v.begin() + 2;
	   itr++)
	tmp += (**itr)->p4();
      
      m_truthOnShellMass                    -> Fill(tmp.M()/1000.0);
      if (evt.physicsEvent.isTruthVisible)
	m_visibleTruthOnShellMass                    -> Fill(tmp.M()/1000.0);
    
      truthOnShellMass = tmp.M();
      // if (TMath::Abs(evt.physicsEvent.GetTruthP4WithPhotons().M()/1000.0 - 90) < 10)
      // 	std::cout << "NEL PICCO: ONVSOFF: " << tmp.M()/1000 << '\t';
	  
      
      tmp.SetPxPyPzE(0.0,0.0,0.0,0.0);
      for (auto itr = truth_v.begin() + 2;
	   itr != truth_v.end(); 
	   itr++) 
	tmp += (**itr)->p4();
      
      m_truthOffShellMass                   -> Fill(tmp.M()/1000.0);

      truthOffShellMass = tmp.M();
      // if (TMath::Abs(evt.physicsEvent.GetTruthP4WithPhotons().M()/1000.0 - 90) < 10)
      // 	std::cout << tmp.M()/1000 << std::endl;
      

      if (evt.physicsEvent.isTruthVisible)
	m_visibleTruthOffShellMass                   -> Fill(tmp.M()/1000.0);
      
      m_truthMassVsTruthOffShellMass        -> Fill(evt.physicsEvent.GetTruthMass()/1000.0, tmp.M()/1000.0, weight);
      
    }
    
  }
  
  if (evt.physicsEvent.isMatchedSet()) {

    m_matchedMass                      -> Fill(1.0*evt.physicsEvent.GetMatchedMass()/1000.0, weight);
    m_matchedMaxPt                     -> Fill(1.0*evt.physicsEvent.GetMatchedMaxPt()/1000.0, weight);
    m_matchedMaxEta                    -> Fill(TMath::Abs(1.0*evt.physicsEvent.GetMatchedWithMaxEta()->eta()), weight);
    m_matchedCharge                    -> Fill(1.0*evt.physicsEvent.GetMatchedCharge(), weight);
    m_matchedSmallestDR                -> Fill(1.0*evt.physicsEvent.GetMatchedSmallestDR(), weight);
     
    m_matchedDZ                        -> Fill(evt.physicsEvent.GetMatchedMaxZ()->z0() - evt.physicsEvent.GetMatchedMinZ()->z0(), weight);   
  }


  if (evt.physicsEvent.isRecoSet()) {

    std::vector<xAOD::TrackParticleContainer::const_iterator> reco_v = evt.physicsEvent.GetReco();
    std::vector<xAOD::TrackParticleContainer::const_iterator> onShell_reco_v = evt.physicsEvent.GetOnShellReco(/*reco_v, 91187.6, 105.6583715*/);
    std::vector<xAOD::TrackParticleContainer::const_iterator> offShell_reco_v = evt.physicsEvent.GetOffShellReco(/*reco_v,91187.6, 105.6583715*/);
    const std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> >& isolation = evt.physicsEvent.GetRecoIsolation();
    
    m_recoMass              -> Fill(1.0*evt.physicsEvent.GetRecoMass()/1000.0, weight);


    if (!evt.physicsEvent.isMatchedSet())
      m_isRecoSetOutsideAcceptance -> Fill(1);
    else
      m_isRecoSetOutsideAcceptance -> Fill(0);

    int NCenteredMuons = 0;
    for (auto itr = reco_v.begin();
	 itr != reco_v.end();
	 itr++) {
      if (TMath::Abs((**itr)->eta()) < 2.7)
	NCenteredMuons++;
    }

    
    double maxIsolation = -1;

    for (auto itr = reco_v.begin();
	 itr != reco_v.end();
	 itr++) {
      for (auto itr2 = isolation.begin();
	   itr2 != isolation.end();
	   itr2++) {
	if ((*itr2).first != (*itr))
	  continue;

	if ((*itr2).second > maxIsolation)
	  maxIsolation = (*itr2).second;

	m_recoIsolationVsPt -> Fill((*itr2).second, ((*(*itr2).first))->pt()/1000.0, weight);
	m_recoSumPtVsPt -> Fill((*itr2).second*((*(*itr2).first))->pt()/1000.0, (*((*itr2).first))->pt()/1000.0, weight);
      }
    }

    m_recoIsolation -> Fill(maxIsolation,weight);

    m_recoNCenteredMuons -> Fill(NCenteredMuons, weight);
    m_recoNCenteredOffShellMuons -> Fill(evt.physicsEvent.GetNCenteredOffShellMuons(/*reco_v, 91187.6, 105.6583715*/), weight);
    
    //ranking
    std::vector<double> rankedPt;
    for (auto itr = reco_v.begin();
	 itr != reco_v.end();
	 itr++) {
      rankedPt.push_back((**itr)->pt());
    }
    
    std::sort(rankedPt.begin(), rankedPt.end());
    
    if (reco_v.size() == 4) {
      m_recoPt1                      -> Fill(rankedPt[3]/1000.0);
      m_recoPt2                      -> Fill(rankedPt[2]/1000.0);
      m_recoPt3                      -> Fill(rankedPt[1]/1000.0);
      m_recoPt4                      -> Fill(rankedPt[0]/1000.0);
    }

    ///////////////////


    if (TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()) < 2.7) {
      m_recoMass27              -> Fill(1.0*evt.physicsEvent.GetRecoMass()/1000.0, weight);
      m_recoOnShellMass27       -> Fill(InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);
      m_recoOnShellEta1VsEta2_27 -> Fill((*(onShell_reco_v[0]))->eta(), (*(onShell_reco_v[1]))->eta(), weight);
      m_recoOffShellMass27      -> Fill(InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);

      m_sigRecoMass27 -> Fill((evt.physicsEvent.GetRecoMass() - evt.physicsEvent.GetTruthP4WithPhotons().M())/1000.0);
      m_sigRecoOnShellMass27 -> Fill((InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOnShellMass)/1000.0);
      m_sigRecoOffShellMass27 -> Fill((InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOffShellMass)/1000.0);
      m_sigPt27 -> Fill((evt.physicsEvent.GetRecoP4().Pt() - evt.physicsEvent.GetTruthP4WithPhotons().Pt())/1000.0);
      m_sigEta27 -> Fill(evt.physicsEvent.GetRecoP4().Eta() - evt.physicsEvent.GetTruthP4WithPhotons().Eta());
      m_sigPhi27 -> Fill(evt.physicsEvent.GetRecoP4().Phi() - evt.physicsEvent.GetTruthP4WithPhotons().Phi());

    }

    else if (TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()) < 3.2) {
      m_recoMass32              -> Fill(1.0*evt.physicsEvent.GetRecoMass()/1000.0, weight);
      m_recoOnShellMass32       -> Fill(InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);
      m_recoOnShellEta1VsEta2_32 -> Fill((*(onShell_reco_v[0]))->eta(), (*(onShell_reco_v[1]))->eta(), weight);
      m_recoOffShellMass32      -> Fill(InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);

      m_sigRecoMass32 -> Fill((evt.physicsEvent.GetRecoMass() - evt.physicsEvent.GetTruthP4WithPhotons().M())/1000.0);
      m_sigRecoOnShellMass32 -> Fill((InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOnShellMass)/1000.0);
      m_sigRecoOffShellMass32 -> Fill((InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOffShellMass)/1000.0);
      m_sigPt32 -> Fill((evt.physicsEvent.GetRecoP4().Pt() - evt.physicsEvent.GetTruthP4WithPhotons().Pt())/1000.0);
      m_sigEta32 -> Fill(evt.physicsEvent.GetRecoP4().Eta() - evt.physicsEvent.GetTruthP4WithPhotons().Eta());
      m_sigPhi32 -> Fill(evt.physicsEvent.GetRecoP4().Phi() - evt.physicsEvent.GetTruthP4WithPhotons().Phi());

    }
    else if (TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()) < 4.0) {
      m_recoMass4               -> Fill(1.0*evt.physicsEvent.GetRecoMass()/1000.0, weight);
      m_recoOnShellMass4        -> Fill(InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);
      m_recoOnShellEta1VsEta2_4 -> Fill((*(onShell_reco_v[0]))->eta(), (*(onShell_reco_v[1]))->eta(), weight);
      m_recoOffShellMass4       -> Fill(InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);

      m_sigRecoMass4 -> Fill((evt.physicsEvent.GetRecoMass() - evt.physicsEvent.GetTruthP4WithPhotons().M())/1000.0);
      m_sigRecoOnShellMass4 -> Fill((InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOnShellMass)/1000.0);
      m_sigRecoOffShellMass4 -> Fill((InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOffShellMass)/1000.0);
      m_sigPt4 -> Fill((evt.physicsEvent.GetRecoP4().Pt() - evt.physicsEvent.GetTruthP4WithPhotons().Pt())/1000.0);
      m_sigEta4 -> Fill(evt.physicsEvent.GetRecoP4().Eta() - evt.physicsEvent.GetTruthP4WithPhotons().Eta());
      m_sigPhi4 -> Fill(evt.physicsEvent.GetRecoP4().Phi() - evt.physicsEvent.GetTruthP4WithPhotons().Phi());

    }

    
    m_sigRecoMass -> Fill((evt.physicsEvent.GetRecoMass() - evt.physicsEvent.GetTruthP4WithPhotons().M())/1000.0);
    m_sigRecoOnShellMass -> Fill((InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOnShellMass)/1000.0);
    m_sigRecoOffShellMass -> Fill((InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass()) - truthOffShellMass)/1000.0);
    m_sigPt -> Fill((evt.physicsEvent.GetRecoP4().Pt() - evt.physicsEvent.GetTruthP4WithPhotons().Pt())/1000.0);
    m_sigEta -> Fill(evt.physicsEvent.GetRecoP4().Eta() - evt.physicsEvent.GetTruthP4WithPhotons().Eta());
    m_sigPhi -> Fill(evt.physicsEvent.GetRecoP4().Phi() - evt.physicsEvent.GetTruthP4WithPhotons().Phi());
    
    m_recoOnShellMass       -> Fill(InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);
    m_recoOnShellEta1VsEta2 -> Fill((*(onShell_reco_v[0]))->eta(), (*(onShell_reco_v[1]))->eta(), weight);
    m_recoOffShellMass      -> Fill(InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);
    m_recoMassVsRecoOffShellMass -> Fill(evt.physicsEvent.GetRecoMass()/1000.0, InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0, 1.0);
    m_recoOnVsOffShellMass       -> Fill(InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,
				    InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,
				    1.0);
    m_recoCharge            -> Fill(1.0*evt.physicsEvent.GetRecoCharge(), weight);
    //m_recoTruthCharge       -> Fill(1.0*evt.physicsEvent.GetRecoTruthCharge(), weight);
    m_recoSmallestDR        -> Fill(1.0*evt.physicsEvent.GetRecoSmallestDR(), weight);
    m_recoPt                -> Fill(1.0*evt.physicsEvent.GetRecoP4().Pt()/1000.0, weight);
    m_recoEta                -> Fill(TMath::Abs(1.0*evt.physicsEvent.GetRecoP4().Eta()), weight);
    m_recoPhi                -> Fill(1.0*evt.physicsEvent.GetRecoP4().Phi(), weight);
    //m_recoTruthPt           -> Fill(1.0*evt.physicsEvent.GetRecoTruthP4().Pt()/1000.0, weight);
    m_recoMaxPt             -> Fill(1.0*evt.physicsEvent.GetRecoMaxPt()/1000.0, weight);
    if (xAOD::TrackHelper::truthParticle(evt.physicsEvent.GetRecoWithMaxEta()))
      m_recoMaxEta                    -> Fill(TMath::Abs(1.0*xAOD::TrackHelper::truthParticle(evt.physicsEvent.GetRecoWithMaxEta())->eta()), weight);

    m_recoDZ                -> Fill((*(reco_v[1]))->z0() - (*(reco_v[0]))->z0(), weight);    

    m_recoMinZVsMaxZ        -> Fill(evt.physicsEvent.GetRecoMinZ()->z0(), evt.physicsEvent.GetRecoMaxZ()->z0(), weight);

    //reco matching
    m_recoMatchedMaxDR      -> Fill(1.0*evt.physicsEvent.GetRecoMatchedMaxDR(), weight);


    for (auto itr = reco_v.begin();
	 itr != reco_v.end();
	 itr++) {
      m_reco_abseta->Fill(TMath::Abs((**itr)->eta()), weight);
    }

    if (evt.physicsEvent.isTruthVisible) {
      m_recoVisibleMass              -> Fill(1.0*evt.physicsEvent.GetRecoMass()/1000.0, weight);
      m_recoVisibleOnShellMass       -> Fill(InvariantMass::TrackInvariantMass(onShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);
      m_recoVisibleOffShellMass      -> Fill(InvariantMass::TrackInvariantMass(offShell_reco_v, evt.physicsEvent.GetDecayMass())/1000.0,1.0);
    }
    

  }


  // if (evt.physicsEvent.isMatchedSet()) {
  //   m_biasMatchedMass                      -> Fill(1.0*(evt.physicsEvent.GetMatchedMass()-evt.physicsEvent.GetTruthMass())/1000.0, weight);
  // }

  // if (evt.physicsEvent.isRecoSet()) {
  //   m_biasRecoMass              -> Fill(1.0*(evt.physicsEvent.GetRecoMass()-evt.physicsEvent.GetTruthMass())/1000.0, weight);
  // }

  // if (evt.physicsEvent.isMatchedSet() && evt.physicsEvent.isRecoSet()) {
  //   m_biasMatchedVsRecoMass     -> Fill(1.0*(evt.physicsEvent.GetRecoMass()-evt.physicsEvent.GetMatchedMass())/1000.0, weight);
  // }


    


  //fake charge in matched tracks
  if (evt.physicsEvent.isMatchedSet()) {
    if (TMath::Abs(evt.physicsEvent.GetMatchedCharge()) > 1E-6) {

      m_matchedChargeVsTruthMaxEta     -> Fill(TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()), evt.physicsEvent.GetMatchedCharge(), weight);        
      m_fakeChargeTruthMinEtaVsTruthMaxEta -> Fill(evt.physicsEvent.GetTruthMinEta(), TMath::Abs(evt.physicsEvent.GetTruthWithMaxEta()->eta()), weight);
      m_fakeChargeTruthMinPtVsTruthMaxPt   -> Fill(evt.physicsEvent.GetTruthMinPt()/1000.0, evt.physicsEvent.GetTruthMaxPt()/1000.0, weight);
      m_fakeChargeBiasMatchedMass          -> Fill((evt.physicsEvent.GetMatchedMass() - evt.physicsEvent.GetTruthMass())/1000.0, weight);
      m_fakeChargeBiasMatchedMaxEtaPtVsTruthMaxEta             -> Fill((evt.physicsEvent.GetMatchedWithMaxEta()->pt() - evt.physicsEvent.GetTruthWithMaxEta()->pt())/1000.0, TMath::Abs(evt.physicsEvent.GetMatchedWithMaxEta()->eta()), weight);
      m_fakeChargePt             -> Fill(evt.physicsEvent.GetMatchedFake()[0]->pt()/1000.0, weight);
      m_fakeChargeEta            -> Fill(evt.physicsEvent.GetMatchedFake()[0]->eta(), weight);
      m_fakeChargeBiasMatchedMaxEtaPt -> Fill(evt.physicsEvent.GetMatchedWithMaxEta()->pt()/1000.0, weight);
      m_fakeChargeMaxEtaBiasCharge -> Fill(evt.physicsEvent.GetMatchedWithMaxEta()->charge() - evt.physicsEvent.GetTruthWithMaxEta()->charge(), weight);
    }
    else {
      m_matchedTrueMaxEta         -> Fill(TMath::Abs(xAOD::TrackHelper::truthParticle(evt.physicsEvent.GetMatchedWithMaxEta())->eta()), weight);
    }
  }
  


  if (evt.physicsEvent.isPairSet()) {
      
    unsigned int nFakeChargeTracks = 0;
    int isMatchingRight = 0;
    bool isTruthParticle = true;

    std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> > truthMatchedPair = evt.physicsEvent.GetPair();
    for (std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >::iterator itr = truthMatchedPair.begin();
	 itr != truthMatchedPair.end();
	 itr++) {
      if (TMath::Abs((*((*itr).first))->charge() - (*((*itr).second))->charge()) > 1E-6)
	nFakeChargeTracks++;

      if (xAOD::TrackHelper::truthParticle(*((*itr).second))) {
	if (TMath::Abs((*((*itr).first))->eta() - xAOD::TrackHelper::truthParticle(*((*itr).second))->eta()) < 1E-6)
	  isMatchingRight++;
      }
      else {isTruthParticle = false; break;}
	

    }
    
    m_NFakeChargeMatchedTracks -> Fill(nFakeChargeTracks, weight);
    
    if (isTruthParticle)
      m_isMatchingRight -> Fill(isMatchingRight, weight);
    else
      m_isMatchingRight -> Fill(-1.0, weight);
    
  }

  
  if (evt.physicsEvent.isRecoPairSet()) {
      
      unsigned int nFakeChargeTracks = 0;
      int isRecoMatchingRight = 0;
      bool isTruthParticle = true;
      
      auto truthMatchedPair = evt.physicsEvent.GetRecoPair();
      for (auto itr = truthMatchedPair.begin();
	   itr != truthMatchedPair.end();
	   itr++) {
	if (TMath::Abs((*((*itr).first))->charge() - (*((*itr).second))->charge()) > 1E-6)
	  nFakeChargeTracks++;

	if (xAOD::TrackHelper::truthParticle(*((*itr).second))) {
	  if (TMath::Abs((*((*itr).first))->eta() - xAOD::TrackHelper::truthParticle(*((*itr).second))->eta()) < 1E-10)
	  
	    isRecoMatchingRight++;
	  else {
	    std::cout << "reco link not corresponding!" << std::endl;
	    std::cout << "Pdg id: " << xAOD::TrackHelper::truthParticle(*((*itr).second))->pdgId() << std::endl;
	    std::cout << "Eta: " << (*((*itr).first))->eta() << '\t' << xAOD::TrackHelper::truthParticle(*((*itr).second))->eta() << std::endl;
	    std::cout << "Phi: " << (*((*itr).first))->phi() << '\t' << xAOD::TrackHelper::truthParticle(*((*itr).second))->phi() << std::endl;
	    std::cout << "Prod vtx: " << xAOD::TrackHelper::truthParticle(*((*itr).second))->prodVtx()->z() << '\t' << xAOD::TrackHelper::truthParticle(*((*itr).second))->prodVtx()->perp() << std::endl;
	    std::cout << "Pt: " << xAOD::TrackHelper::truthParticle(*((*itr).second))->pt() << std::endl;
	  }
	}
	else {isTruthParticle = false; break;}
	  


      }
      
  

      m_NFakeChargeRecoTracks -> Fill(nFakeChargeTracks,weight);
      m_NFakeChargeRecoVsDR   -> Fill(nFakeChargeTracks, evt.physicsEvent.GetRecoMatchedMaxDR(), weight);

      if (isTruthParticle)
	m_isRecoMatchingRight -> Fill(isRecoMatchingRight, weight);
      else
	m_isRecoMatchingRight -> Fill(-1.0, weight);

  }

}


   
    
  
 
      
    

      




