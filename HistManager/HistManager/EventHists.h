#ifndef HistManager_EventHists_H
#define HistManager_EventHists_H

#include <HistManager/EventFeaturesInterface.h>
#include <HistManager/ManageHists.h>

class EventHists : public ManageHists {

  public:

    EventHists(TString name);
    ~EventHists();

    void BookHists();
    void FillHists(EventFeatures &evt, float weight) const; 

  private:

    TString m_name;
    TString m_label;

  public:

    TH1F* m_pstrk; //!
    TH1F* m_primary_pstrk; //!
    TH1F* m_secondary_pstrk; //!
    TH1F* m_lost_pstrk; //!
    TH1F* m_recotrk; //!    
    TH1F* m_recotrkwide; //!    
    TH1F* m_sispseed; //!       
    TH1F* m_nPixelClusters; //!
    TH1F* m_nSCTClusters; //!
    TH1F* m_nSiClusters; //!
    TH1F* m_nPixelClusterswide; //!
    TH1F* m_nSCTClusterswide; //!
    TH1F* m_nSiClusterswide; //!
    TH1F* m_trigWeight; //!
    TH1F* m_muValue; //!
    TH1F* m_muValueCorr; //!
    TH1F* m_LumiBlock; //!

    // For ITk
    TH1F* m_recotrk1GeV; //!
    TH1F* m_nIsoTrack; //!
    TH1F* m_nVertex; //!
    TH1F* m_average2TrackEta; //!
    TH1F* m_photonEta; //!
    TH1F* m_photonEta2Track; //!
    TH1F* m_photonProdR; //!
    TH1F* m_photonProdR2Track; //!
    TH1F* m_photonDecayR; //!
    TH1F* m_photonDecayR2Track; //!
    TH1F* m_nTruth; //!
    TH1F* m_nChargedTruth; //!
    TH1F* m_nPrimaryChargedTruth; //!
    TH1F* m_nNeutralTruth; //!
    TH1F* m_nTruth1GeV; //!
    TH1F* m_nChargedTruth1GeV; //!
    TH1F* m_nNeutralTruth1GeV; //!
    TH1F* m_effPileupTrack1GeV; //!
    TH2F* m_effPileupTrack1GeV_vsTruth; //!

    TH1F* m_nPixelClustersAtIBL; //!
    TH1F* m_nPixelHitsAtIBL; //!
    TH2F* m_nPixelClustersAtIBL_vsMu; //!
    TH2F* m_nPixelHitsAtIBL_vsMu; //!
    TH1F* m_nPixelClustersAtIBL_wide; //!
    TH1F* m_nPixelHitsAtIBL_wide; //!
    TH2F* m_nPixelClustersAtIBL_vsMu_wide; //!
    TH2F* m_nPixelHitsAtIBL_vsMu_wide; //!

    TH1F* m_nPixelClustersAtBLayer; //!
    TH1F* m_nPixelHitsAtBLayer; //!
    TH2F* m_nPixelClustersAtBLayer_vsMu; //!
    TH2F* m_nPixelHitsAtBLayer_vsMu; //!
    TH1F* m_nPixelClustersAtBLayer_wide; //!
    TH1F* m_nPixelHitsAtBLayer_wide; //!
    TH2F* m_nPixelClustersAtBLayer_vsMu_wide; //!
    TH2F* m_nPixelHitsAtBLayer_vsMu_wide; //!

    TH1F* m_truthMass; //!
    TH1F* m_truthOutsideDetector; //!
    TH1F* m_truthNCenteredMuons; //!
    TH1F* m_truthPt; //!
    TH1F* m_truthEta; //!
    TH1F* m_truthPz; //!
    TH1F* m_truthMaxPt; //!
    TH1F* m_truthMinPt; //!
    TH1F* m_truthMaxEta; //!
    TH1F* m_truthMinEta; //!
    TH1F* m_truthCharge; //!
    TH1F* m_truthSmallestDR; //!
    TH1F* m_truthDZ; //!

    TH2F* m_truthEta1VsEta2; //!
    
    TH1F* m_truthPt1; //!
    TH1F* m_truthPt2; //!
    TH1F* m_truthPt3; //!
    TH1F* m_truthPt4; //!

    TH1F* m_isTruthPt1OnShell; //!
    TH1F* m_isTruthPt2OnShell; //!
    TH1F* m_isTruthPt3OnShell; //!
    TH1F* m_isTruthPt4OnShell; //!
    TH1F* m_isTruthPtRankingOrdered; //!
    TH1F* m_truthCutCode; //!

    TH1F* m_visibleTruth_abseta; //!
    TH1F* m_reco_abseta; //!

    TH1F* m_matchedMass; //!
    TH1F* m_matchedMaxPt; //!
    TH1F* m_matchedMinPt; //!
    TH1F* m_matchedMaxEta; //!
    TH1F* m_matchedTrueMaxEta; //!
    TH1F* m_matchedMinEta; //!
    TH1F* m_matchedCharge; //!
    TH1F* m_matchedSmallestDR; //!
    TH1F* m_matchedDZ; //!
    TH1F* m_isMatchingRight; //!

    TH1F* m_recoMass; //!
    TH1F* m_recoMass27; //!
    TH1F* m_recoMass32; //!
    TH1F* m_recoMass4; //!
    TH1F* m_recoCharge; //!
    TH1F* m_recoTruthCharge; //!
    TH1F* m_recoMinDR; //!
    TH1F* m_recoPt; //!
    TH1F* m_recoPt1; //!
    TH1F* m_recoPt2; //!
    TH1F* m_recoPt3; //!
    TH1F* m_recoPt4; //!
    TH1F* m_recoTruthPt; //!
    TH1F* m_recoMaxPt; //!
    TH1F* m_recoMinPt; //!
    TH1F* m_recoMaxEta; //!
    TH1F* m_recoMinEta; //!
    TH1F* m_recoSmallestDR; //!
    TH1F* m_recoDZ; //!
    TH2F* m_recoMinZVsMaxZ; //!
    TH1F* m_recoNCenteredMuons; //!
    TH1F* m_recoNCenteredOffShellMuons; //!
    //reco matching
    TH1F* m_recoMatchedMaxDR; //!
    TH1F* m_isRecoMatchingRight; //!

    TH1F* m_recoCutCode; //!
    TH1F* m_recoCutCode27; //!
    TH1F* m_recoCutCode32; //!
    TH1F* m_recoCutCode4;  //!

    TH1F* m_biasMatchedMass; //!
    TH1F* m_biasRecoMass; //!
    TH1F* m_biasMatchedVsRecoMass; //!

    //fake charge in matched tracks
    TH2F* m_matchedChargeVsTruthMaxEta; //!
    TH2F* m_fakeChargeTruthMinEtaVsTruthMaxEta; //!
    TH2F* m_fakeChargeTruthMinPtVsTruthMaxPt; //!
    TH1F* m_fakeChargePt; //!
    TH1F* m_fakeChargeEta; //!
    TH1F* m_fakeChargeBiasMatchedMass; //!
    TH2F* m_fakeChargeBiasMatchedMaxEtaPtVsTruthMaxEta; //!
    TH1F* m_fakeChargeBiasMatchedMaxEtaPt;
    TH1F* m_fakeChargeMaxEtaBiasCharge; //!
    TH1F* m_NFakeChargeMatchedTracks; //!
    TH2F* m_NFakeChargeRecoVsDR; //!
    TH1F* m_NFakeChargeRecoTracks; //!


    //Higgs Analysis
    TH1F* m_truthOnShellMass; //!
    TH1F* m_truthOffShellMass; //!
    TH2F* m_truthMassVsTruthOffShellMass; //!
    TH1F* m_visibleTruthMass; //!
    TH1F* m_visibleTruthOnShellMass; //!
    TH1F* m_visibleTruthOffShellMass; //!
    TH1F* m_recoVisibleMass; //!
    TH1F* m_recoOnShellMass; //!
    TH2F* m_recoOnShellEta1VsEta2; //!
    TH2F* m_recoOnVsOffShellMass; //!
    TH1F* m_recoOnShellMass27; //!
    TH2F* m_recoOnShellEta1VsEta2_27; //!
    TH1F* m_recoOnShellMass32; //!
    TH2F* m_recoOnShellEta1VsEta2_32; //!
    TH1F* m_recoOnShellMass4; //!
    TH2F* m_recoOnShellEta1VsEta2_4; //!
    TH1F* m_recoOffShellMass; //!
    TH1F* m_recoOffShellMass27; //!
    TH1F* m_recoOffShellMass32; //!
    TH1F* m_recoOffShellMass4; //!

    TH2F* m_recoMassVsRecoOffShellMass; //!

    TH1F* m_recoVisibleOnShellMass; //!
    TH1F* m_recoVisibleOffShellMass; //!
    
    TH1F* m_recoIsolation; //!
    TH2F* m_recoIsolationVsPt; //!
    TH2F* m_recoSumPtVsPt; //!

};

#endif //HistManager_EventHists_H
