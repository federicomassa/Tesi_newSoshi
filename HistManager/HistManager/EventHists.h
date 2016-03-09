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

};

#endif //HistManager_EventHists_H
