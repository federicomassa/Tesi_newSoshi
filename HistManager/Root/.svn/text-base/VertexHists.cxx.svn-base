#include <HistManager/VertexHists.h>


VertexHists::VertexHists(TString name) {
   m_name = "VertexHist_" + name;
   m_label = "Vertex "; // don't forget the space
}

VertexHists::~VertexHists() {}

void VertexHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  /* positions */
  m_globalX = declare1D(m_name, "globalX", m_label + "global X [mm]",150,-150.0,150.0);
  m_globalY = declare1D(m_name, "globalY", m_label + "global Y [mm]",150,-150.0,150.0);
  m_globalZ = declare1D(m_name, "globalZ", m_label + "global Z [mm]",100,-200.0,200.0);

  m_globalXError = declare1D(m_name, "globalXError", m_label + "global X error [mm]",100,0.0, 5.0);
  m_globalYError = declare1D(m_name, "globalYError", m_label + "global Y error [mm]",100,0.0, 5.0);
  m_globalZError = declare1D(m_name, "globalZError", m_label + "global Z error [mm]",100,0.0, 5.0);

  m_chi2    = declare1D(m_name, "chi2",    m_label + "#chi^{2}",        100, 0.0,50.0);
  m_ndof    = declare1D(m_name, "ndof",    m_label + "# DoF",            50,-0.5,49.5);
  m_chi2dof = declare1D(m_name, "chi2dof", m_label + "#chi^{2}/dof",    100, 0.0,10.0);
  m_numTrk  = declare1D(m_name, "numTrk",  m_label + "Number of tracks", 50,-0.5,49.5);

  m_Lxy    = declare1D(m_name, "Lxy", m_label + "Decay length [mm]",100,0,100);
  m_Lxysig = declare1D(m_name, "Lxysig", m_label + "Decay length significance",100,0,100);


} // BookHists

void VertexHists::FillHists(const xAOD::Vertex* vtx, float weight) const {

  /* positions */
  m_globalX -> Fill(vtx->x(),weight);
  m_globalY -> Fill(vtx->y(),weight);
  m_globalZ -> Fill(vtx->z(),weight);

  const AmgSymMatrix(3) covVtx = vtx->covariancePosition();

  float xerr = sqrt(covVtx(0,0));
  float yerr = sqrt(covVtx(1,1));
  float zerr = sqrt(covVtx(2,2));
  m_globalXError -> Fill(xerr,weight);
  m_globalYError -> Fill(yerr,weight);
  m_globalZError -> Fill(zerr,weight);

  m_chi2    -> Fill(vtx->chiSquared(),weight);
  m_ndof    -> Fill(vtx->numberDoF(),weight);
  m_chi2dof -> Fill(vtx->chiSquared()/vtx->numberDoF(),weight);
  m_numTrk  -> Fill(1.0*vtx->nTrackParticles(),weight);

  float Lxy    = sqrt(vtx->x()*vtx->x()+vtx->y()*vtx->y());  // this is unsigend version
  float Lxysig = Lxy/sqrt(xerr*xerr+yerr*yerr);

  m_Lxy    -> Fill(Lxy,weight);
  m_Lxysig -> Fill(Lxysig,weight);



} // FillHists

