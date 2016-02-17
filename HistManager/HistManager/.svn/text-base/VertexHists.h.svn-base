#ifndef HistManager_VertexHists_H
#define HistManager_VertexHists_H

#include "HistManager/ManageHists.h"
#ifndef __MAKECINT__
#include "xAODTracking/Vertex.h"
#endif // not __MAKECINT__

class VertexHists : public ManageHists {

 public:

  VertexHists(TString name);
  ~VertexHists();

  void BookHists();
#ifndef __MAKECINT__
  void FillHists(const xAOD::Vertex* vtx, float weight) const;
#endif // not __MAKECINT__


  private:

    TString m_name;
    TString m_label;

  public:

    /* positions */
    TH1F* m_globalX; //!
    TH1F* m_globalY; //!
    TH1F* m_globalZ; //!
    TH1F* m_globalXError; //!
    TH1F* m_globalYError; //!
    TH1F* m_globalZError; //!

    TH1F* m_chi2; //!
    TH1F* m_ndof; //!
    TH1F* m_chi2dof; //!
    TH1F* m_numTrk; //!

    TH1F* m_Lxy; //!
    TH1F* m_Lxysig; //!


};

#endif //HistManager_VertexHists_H
