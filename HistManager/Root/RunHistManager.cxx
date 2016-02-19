/******************************************
 * Plot of run info extracted from track plots
 *
 * F. Massa
 * Feb 19 15:50:00 CEST 2016
 *
 ******************************************/
#include <HistManager/TrackHists.h>
#include <HistManager/TrackHistManager.h>
#include <HistManager/RunHistManager.h>
#include <HistManager/RunHists.h>

//name identifies the type of tracks reconstructed, for example reco (matched) or all (matched+pileup).
//Every RunHists object contained in the m_runHistList vector matches the correspondent TrackHists object contained in m_trackHistList
RunHistManager::RunHistManager(TrackHistManager* trackHistManager):m_trackHistList(trackHistManager->GetList()) {
  m_name = trackHistManager->GetName();
}

RunHistManager::~RunHistManager() {}

void RunHistManager::Init(EL::Worker* wk) {

  std::cout << "Init RunHistManager " << m_name << std::endl;

  //Now initialize the worker and add RunHists objects to vector, needed for RunHistManager
  for (std::vector<TrackHists*>::const_iterator trk_itr = m_trackHistList.begin(); trk_itr != m_trackHistList.end(); trk_itr++) {
    m_runHistList.push_back(new RunHists(*trk_itr));
    m_runHistList.back()->Init(wk);
  }
} //Init

void RunHistManager::FillHists(float weight) const {
  for (std::vector<RunHists*>::const_iterator run_itr = m_runHistList.begin(); run_itr != m_runHistList.end(); run_itr++) {
    (*run_itr)->FillHists(weight);
  }
} // FillHists

