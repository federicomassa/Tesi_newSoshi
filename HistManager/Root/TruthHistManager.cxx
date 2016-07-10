/******************************************
 * Plot of truth quantities
 *
 * G. Facini + F.Massa
 * Apr 5 15:53:25 CEST 2016
 *
 ******************************************/

#include <HistManager/TruthHistManager.h>

TruthHistManager::TruthHistManager(TString name) {
   m_name = name;
   m_all               = new TruthHists(name + "_all");   
}

TruthHistManager::~TruthHistManager() {}

void TruthHistManager::Init(EL::Worker* wk) {

  std::cout << "Init TruthHistManager " << m_name << std::endl;

  //Now initialize the worker and add TruthHists objects to vector, needed for RunHistManager
  m_all       -> Init( wk );
  m_histList       . push_back(m_all);

} // Init

void TruthHistManager::FillHists(const xAOD::TruthParticle* tp, float weight) const {
  //if isHardTruth == false and isPileupTruth == false just fill the common part
  m_all->FillHists( tp,weight);
  return;

} // FillHists
