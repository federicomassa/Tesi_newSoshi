#ifndef HistManager_TrackHists_H
#define HistManager_TrackHists_H

#include <TMath.h>
#include "HistManager/ManageHists.h"
#include "HistManager/TrackHelper.h"
#include <set>
#include <vector>

#ifndef __MAKECINT__
#include "xAODTracking/TrackParticle.h"
#endif // not __MAKECINT__

class TrackHists : public ManageHists {

  public:

    TrackHists(TString name);
    ~TrackHists();

    void BookHists();
#ifndef __MAKECINT__
    void FillHists(const xAOD::TrackParticle* trk, float weight) const;


    //* Determine if the TruthParticle has been used before */
    // CHANGED BY FEDERICO MASSA: now using trackHelper to retrieve truthParticle: tracks do not have barcodes!!!
    bool hasBeenUsed( const xAOD::TrackParticle* trk) {
      if (m_usedBarcodes.find(xAOD::TrackHelper::truthParticle(trk)->auxdata<int>("barcode"))==m_usedBarcodes.end()) { return false; } 
      else { return true; }
    }

    //* Determine if the TruthParticle has been used before */
    void resetBarcodes() {
      m_usedBarcodes.clear();
    }

#endif // not __MAKECINT__

  private:
    const double m_etaMax = 5.0;
    const int m_etaVectorSize = TMath::Nint(m_etaMax)*5;
    TString m_name;
    TString m_label;
    bool m_isPseudo;
    mutable std::set<int> m_usedBarcodes; //!

    TH1F* m_muValue;
    TH1F* m_lumibk;
    TH1F* m_pixClus;

    /* particle properties */
    TH1F* m_prodR; //!
    TH1F* m_prodZ; //!
    TH1F* m_truthPt; //!
    TH1F* m_truthPtnarrow; //!
    TH1F* m_truthEta; //!
    TH1F* m_truthD0; //!
    TH1F* m_truthD0_wide; //!
    TH1F* m_truthZ0; //!
    TH1F* m_truthPhi; //!
    TH1F* m_biasPt; //!
    TH1F* m_biasQPt; //!
    TH1F* m_biasEta; //!
    TH1F* m_biasPhi; //!
    TH1F* m_biasD0; //!
    TH1F* m_biasZ0; //!

    std::vector<TH1F*> m_biasPt_abseta; //!
    std::vector<TH1F*> m_biasQPt_abseta; //!
    std::vector<TH1F*> m_biasPhi_abseta; //!
    std::vector<TH1F*> m_biasD0_abseta; //!
    std::vector<TH1F*> m_biasZ0_abseta; //!

    /* track parameterization */
    TH1F* m_qoverp; //!     
    TH1F* m_pt; //!       
    TH1F* m_ptnarrow; //!       
    TH1F* m_eta; //!        
    TH1F* m_phi; //!        
    TH1F* m_d0; //!       
    TH1F* m_z0; //!       
    TH1F* m_z0Corr1; //!       
    TH1F* m_z0Corr2; //!       

    /* track parameterization errors */
    TH1F* m_qoverp_err; //!       
    TH1F* m_eta_err; //!        
    TH1F* m_d0_err; //!       
    TH1F* m_d0signif; //!       
    TH1F* m_z0_err; //!       
    TH1F* m_phi_err; //!        

    TH1F* m_d0sign; //!       
    TH1F* m_d0signsig; //!       
    TH1F* m_d0signsigP; //!       
    TH1F* m_d0signsigN; //!       

    /* track properties */
    TH1F* m_nPixHits; //! 
    TH1F* m_nSCTHits; //! 
    TH1F* m_nSiHits; //!

    TH2F* m_eta_nPixHits; //!
    TH2F* m_eta_nSCTHits; //!
    TH2F* m_eta_nSiHits; //!

    TH1F* m_nGangedPix; //! 
    TH1F* m_nGangedPixFF; //! 
    TH1F* m_nPixLay; //! 
    TH1F* m_nPixSharedHits; //!        
    TH1F* m_nPixSplitHits; //!        
    TH1F* m_nPixSpoiltHits; //!       
    TH1F* m_nPixOutliers; //!       
    TH1F* m_nPixHoles; //!        
    TH1F* m_nPixelDeadSensors; //! 
    TH1F* m_nSCTSharedHits; //!       
    TH1F* m_nSCTSpoiltHits; //!       
    TH1F* m_nSCTOutliers; //!       
    TH1F* m_nSCTHoles; //!   
    TH1F* m_nSCTDoubleHoles; //!        
    TH1F* m_nSCTDeadSensors; //!
    TH1F* m_nSiDeadSensors; //!

    TH1F* m_nTRTHits; //! 
    TH1F* m_nTRTOutliers; //!       
    TH1F* m_nTRTHoles; //!        
    TH1F* m_nTRTHTHits; //! 

    TH1F* m_chiSqPerDof; //!
    TH1F* m_nOutliers; //!
    TH1F* m_stdDevChi2OS; //!
    TH1F* m_eProbComb; //!
    TH1F* m_eProbHT; //!
    TH1F* m_eProbToT; //!
    TH1F* m_eProbBrem; //!

    TH1F* m_truthMatchProb; //!
    TH1F* m_matchingDR; //!
    TH1F* m_abseta; //!
    TH1F* m_abseta_dr001; //!
    TH1F* m_abseta_dr002; //!
    TH1F* m_abseta_dr003; //!
    TH1F* m_abseta_dr004; //!
    TH1F* m_abseta_dr005; //!
    TH1F* m_abseta_dr006; //!
    TH1F* m_abseta_dr007; //!
    TH1F* m_abseta_dr008; //!
    TH1F* m_abseta_dr009; //!

    // Cluster study
    TH1F* m_IBLCharge; //!
    TH1F* m_IBLHitSize; //!
    TH1F* m_IBLHitSizeZ; //!
    TH1F* m_IBLHitSizePhi; //!
    TH1F* m_IBLlocalPhi; //!
    TH1F* m_IBLlocalTheta; //!
    TH2F* m_IBLlocalPhi_sizePhi; //!
    TH2F* m_IBLlocalTheta_sizeZ; //!
    TH1F* m_IBLbiasedResidualX; //!
    TH1F* m_IBLbiasedResidualY; //!
    TH1F* m_IBLunbiasedResidualX; //!
    TH1F* m_IBLunbiasedResidualY; //!
    TH1F* m_IBLResidualX; //!
    TH1F* m_IBLResidualY; //!

    TH1F* m_BLayerCharge; //!
    TH1F* m_BLayerHitSize; //!
    TH1F* m_BLayerHitSizeZ; //!
    TH1F* m_BLayerHitSizePhi; //!
    TH1F* m_BLayerlocalPhi; //!
    TH1F* m_BLayerlocalTheta; //!
    TH2F* m_BLayerlocalPhi_sizePhi; //!
    TH2F* m_BLayerlocalTheta_sizeZ; //!
    TH1F* m_BLayerbiasedResidualX; //!
    TH1F* m_BLayerbiasedResidualY; //!
    TH1F* m_BLayerunbiasedResidualX; //!
    TH1F* m_BLayerunbiasedResidualY; //!
    TH1F* m_BLayerResidualX; //!
    TH1F* m_BLayerResidualY; //!

};

#endif //HistManager_TrackHists_H
