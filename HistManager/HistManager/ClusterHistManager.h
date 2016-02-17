#ifndef HistManager_ClusterHistManager_H
#define HistManager_ClusterHistManager_H

#include "HistManager/ClusterHists.h"
#ifndef __MAKECINT__
#include "xAODTracking/TrackMeasurementValidation.h"
#endif // not __MAKECINT__

#include <iostream>

class ClusterHistManager {

 public:

  ClusterHistManager(TString name, int layer = -1, bool more = false);
  ~ClusterHistManager();

  void Init(EL::Worker* wk);
#ifndef __MAKECINT__
  void FillHists(const xAOD::TrackMeasurementValidation* prd, float weight) const;
#endif // not __MAKECINT__


  private:

    TString m_name;
    int m_layer;       // -1 = all
    bool m_more;       // false = default

  public:

    ClusterHists* m_all;          //!
    ClusterHists* m_single;       //!
    ClusterHists* m_single_split; //!
    ClusterHists* m_merged;       //!
    ClusterHists* m_merged_split; //!
    ClusterHists* m_merged2;       //!
    ClusterHists* m_merged2_split; //!
    ClusterHists* m_merged3;       //!
    ClusterHists* m_merged3_split; //!
    ClusterHists* m_merged4;       //!
    ClusterHists* m_merged4_split; //!

};

#endif //HistManager_ClusterHistManager_H
