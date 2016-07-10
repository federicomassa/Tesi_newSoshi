#ifndef HistManager_TruthHists_H
#define HistManager_TruthHists_H

#include "HistManager/ManageHists.h"
#ifndef __MAKECINT__
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#endif // not __MAKECINT__

class TruthHists : public ManageHists {

  friend class RunHists;
  
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
    TString GetName() const {return m_baseName;}
    /* particle/mother properties */
    TH1F* m_prodR; //!
    TH1F* m_prodRbarrel; //!
    TH1F* m_prodZ; //!
    TH1F* m_decayR; //!

    /* track parameterization */
    TH1F* m_pt; //!			  
    TH2F* m_ptVsEta; //!			  
    TH1F* m_eta; //!
    TH1F* m_pdgId; //!
    TH1F* m_eta_wide; //!
    TH1F* m_abseta; //!			  
    TH1F* m_d0; //!			  
    TH1F* m_z0; //!			  
    TH1F* m_phi; //!			  

    TH1F* m_prodVtx_X; //!
    TH1F* m_prodVtx_Y; //!
    TH1F* m_prodVtx_Z; //!
    TH1F* m_prodVtx_R; //!

    TH1F* m_prodVtx_incoming; //!
    TH1F* m_prodVtx_outgoing; //!

    TString m_baseName;

};

#endif //HistManager_TruthHists_H
