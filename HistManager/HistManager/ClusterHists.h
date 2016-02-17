#ifndef HistManager_ClusterHists_H
#define HistManager_ClusterHists_H

#include "HistManager/ManageHists.h"
#ifndef __MAKECINT__
#include "xAODTracking/TrackMeasurementValidation.h"
#endif // not __MAKECINT__

class ClusterHists : public ManageHists {

 public:

  ClusterHists(TString name);
  ~ClusterHists();

  void BookHists();
#ifndef __MAKECINT__
  void FillHists(const xAOD::TrackMeasurementValidation* prd, float weight) const;
#endif // not __MAKECINT__


  private:

    TString m_name;
    TString m_label;

  public:

    /* layout */
    TH2F* m_SiliconLayout_full_xy; //!
    TH2F* m_SiliconLayout_full_zr; //!
    TH2F* m_SiliconLayout_quarter_xy; //!
    TH2F* m_SiliconLayout_quarter_zr; //!
    TH2F* m_PixelLayout_full_xy; //!
    TH2F* m_PixelLayout_full_zr; //!
    TH2F* m_PixelLayout_quarter_xy; //!
    TH2F* m_PixelLayout_quarter_zr; //!

    /* positions */
    TH1F* m_localX; //!
    TH1F* m_localY; //!
    TH1F* m_localXError; //!
    TH1F* m_localYError; //!
    TH1F* m_globalX; //!
    TH1F* m_globalY; //!
    TH1F* m_globalZ; //!

 
    /* basic properties */
    TH1F* m_layer; //!
    TH1F* m_nPart; //!
    TH1F* m_nPriPart; //!
    TH1F* m_size; //!
    TH1F* m_sizePhi; //!
    TH1F* m_sizeZ; //!
    TH1F* m_charge; //!
    TH1F* m_ToT; //!
    TH1F* m_isFake; //!
    TH1F* m_gangedPixel; //!
    TH1F* m_isSplit; //!
    TH1F* m_splitProb1; //!
    TH1F* m_splitProb2; //!
    TH1F* m_NN_sizeX; //!
    TH1F* m_NN_sizeY; //!
    TH1F* m_NN_phiBS; //!
    TH1F* m_NN_thetaBS; //!

    /* fun with 2D */
    TH2F* m_size_VS_charge; //!
    TH2F* m_sizePhi_VS_charge; //!
    TH2F* m_sizeZ_VS_charge; //!
    TH2F* m_size_VS_ToT; //!
    TH2F* m_sizePhi_VS_ToT; //!
    TH2F* m_sizeZ_VS_ToT; //!

    // hits in cluster
    TH1F* m_nHitCluster; //!
    TH1F* m_hitCharge; //!
    TH1F* m_hitTot; //!
    TH1F* m_hitMxTot; //!
    TH1F* m_hitMnTot; //!
    TH1F* m_hitAveTot; //!
    TH1F* m_hitMxCharge; //!
    TH1F* m_hitMnCharge; //!
    TH1F* m_hitAveCharge; //!
    TH1F* m_hit1Charge; //!
    TH1F* m_hit1Tot; //!

};

#endif //HistManager_ClusterHists_H
