/******************************************
 * Plot of truth quantities
 *
 * G. Facini + F.Massa
 * Apr 5 15:53:25 CEST 2016
 *
 ******************************************/

#include <HistManager/TruthHistManager.h>

TruthHistManager::TruthHistManager(TString name, bool doAll, bool doGun, bool doPileup) {
   m_name = name;
   
   m_doAll   = doAll;
   m_doGun = doGun;
   m_doPileup  = doPileup;

   m_all               = new TruthHists(name+"_all");
   if(m_doGun) {
      m_gun        = new TruthHists(name+"_gun");
   }
   if(m_doPileup) { 
      m_pileup      = new TruthHists(name+"_pileup");
   }
   
}

TruthHistManager::~TruthHistManager() {}

void TruthHistManager::Init(EL::Worker* wk) {

  std::cout << "Init TruthHistManager " << m_name << std::endl;

  //Now initialize the worker and add TruthHists objects to vector, needed for RunHistManager

  if(m_doAll) {
    m_all       -> Init( wk );
    m_histList       . push_back(m_all);
  }
  if(m_doGun) { 
    m_gun     -> Init( wk );
    m_histList       . push_back(m_gun);
  }
  if(m_doPileup) {   
    m_pileup          -> Init( wk );
    m_histList       . push_back(m_pileup);
  }

} // Init

void TruthHistManager::FillHists(const xAOD::TruthParticle* tp, float weight, bool isHardTruth = false, bool isPileupTruth = false) const {
  //if isHardTruth == false and isPileupTruth == false just fill the common part
  if (isHardTruth && isPileupTruth) {
    std::cout << "ERROR\t in TruthHistManager::FillHists(...): truth particle cannot be both gun and pileup" << std::endl;
    exit(1);
  }

  if( m_doAll)
    m_all->FillHists( tp,weight);
  if (m_doGun && isHardTruth)
    m_gun->FillHists( tp, weight);
  if (m_doPileup && isPileupTruth) 
    m_pileup->FillHists(tp, weight);
   

  return;

} // FillHists
