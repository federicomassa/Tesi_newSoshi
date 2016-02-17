#ifndef HistManager_TrackHistManager_H
#define HistManager_TrackHistManager_H

#include "HistManager/TrackHists.h"
#ifndef __MAKECINT__
#include "xAODTracking/TrackMeasurementValidation.h"
#include <HistManager/TrackHelper.h>
#endif // not __MAKECINT__

#include <iostream>

class TrackHistManager {

 public:

  TrackHistManager(TString name, bool doPrimary, bool doSecondary, bool doFake, bool doDups, bool splitOnHits);
  ~TrackHistManager();

  void Init(EL::Worker* wk);
  void resetBarcodes();
#ifndef __MAKECINT__
  void FillHists(const xAOD::TrackParticle* trk, float weight) const;
#endif // not __MAKECINT__


  private:

    TString m_name;
    bool m_doPrimary;
    bool m_doSecondary;
    bool m_doFake;
    bool m_doDups;
    bool m_splitOnHits;

  public:

    TrackHists* m_all;          //!
    TrackHists* m_primary;      //!
    TrackHists* m_primary_dup;  //!
    TrackHists* m_secondary;    //!
    TrackHists* m_secondary_dup;//!
    TrackHists* m_fake;         //!
    TrackHists* m_true;
    TrackHists* m_primary_bhit;      //!
    TrackHists* m_primary_bhitShare; //!
    TrackHists* m_primary_bhitSplit; //!
    TrackHists* m_primary_bhitMerge; //!

};

#endif //HistManager_TrackHistManager_H
