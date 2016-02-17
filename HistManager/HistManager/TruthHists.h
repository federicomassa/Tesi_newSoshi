#ifndef HistManager_TruthHists_H
#define HistManager_TruthHists_H

#include "HistManager/ManageHists.h"
#ifndef __MAKECINT__
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#endif // not __MAKECINT__

class TruthHists : public ManageHists {

 public:

  TruthHists(TString name);
  ~TruthHists();

  void BookHists();
#ifndef __MAKECINT__
  void FillHists(const xAOD::TruthParticle* trk, float weight) const;
#endif // not __MAKECINT__


  private:

    TString m_name;
    TString m_label;

  public:

    /* particle/mother properties */
    TH1F* m_prodR; //!
    TH1F* m_prodRbarrel; //!
    TH1F* m_prodZ; //!
    TH1F* m_decayR; //!

    /* track parameterization */
    TH1F* m_pt; //!			  
    TH1F* m_eta; //!			  
    TH1F* m_abseta; //!			  
    TH1F* m_d0; //!			  
    TH1F* m_z0; //!			  
    TH1F* m_phi; //!			  

};

#endif //HistManager_TruthHists_H
