#ifndef HistManager_TruthHistManager_H
#define HistManager_TruthHistManager_H

#include "HistManager/TruthHists.h"

#include <iostream>
#include <vector>

class TruthHistManager {

 public:

  TruthHistManager(TString name, bool doAll, bool doGun, bool doPileup);
  ~TruthHistManager();

  void Init(EL::Worker* wk);
  void FillHists(const xAOD::TruthParticle* tp, float weight, bool isHardTruth, bool isPileupTruth) const;


  private:

    TString m_name;
    bool m_doAll;
    bool m_doGun;
    bool m_doPileup;

    TruthHists* m_all;          //!
    TruthHists* m_gun;          //!
    TruthHists* m_pileup;       //!

    std::vector<TruthHists*> m_histList;

 public:
    const std::vector<TruthHists*>& GetList() const {
      return m_histList;
    }
    TString GetName() const {return m_name;}
};

#endif //HistManager_TruthHistManager_H
