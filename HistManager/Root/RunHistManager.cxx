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
#include <HistManager/assert.h>

//name identifies the type of tracks reconstructed, for example reco (matched) or all (matched+pileup).
//Every RunHists object contained in the m_runHistList vector matches the correspondent TrackHists object contained in m_trackHistList
RunHistManager::RunHistManager(const TrackHistManager* const trackHistManager, const TruthHists* const truthHist):
  m_trackHistList(trackHistManager->GetList()),
  m_truthHist(truthHist) 
{
  m_name = trackHistManager->GetName();
}

RunHistManager::~RunHistManager() {}

void RunHistManager::Init(EL::Worker* wk) {

  std::cout << "Init RunHistManager " << m_name << std::endl;

  m_wk = wk;
  //Now initialize the worker and add RunHists objects to vector, needed for RunHistManager
  for (std::vector<TrackHists*>::const_iterator trk_itr = m_trackHistList.begin(); trk_itr != m_trackHistList.end(); trk_itr++) {
    m_runHistList.push_back(new RunHists(*trk_itr, m_truthHist));
    m_runHistList.back()->Init(wk);
  }

  BookHists();

} //Init

TGraphErrors* RunHistManager::declareGraphErrors(TString cat, TString name, TString xaxis, TString yaxis) {
  TGraphErrors* tmp = new TGraphErrors();
  tmp->SetNameTitle(TString(cat+"__"+name),name);
  tmp->GetXaxis()->SetTitle(xaxis);
  tmp->GetYaxis()->SetTitle(yaxis);
  m_wk->addOutput(tmp);
  return tmp;
}; // declareGraphErrors


void RunHistManager::BookHists() {
  m_fakeProb_abseta = declareGraphErrors("RunHistManager", "fakeProb_abseta", "|#eta|", "Fake probability");
}

void RunHistManager::FillHists(float weight) const {
  for (std::vector<RunHists*>::const_iterator run_itr = m_runHistList.begin(); run_itr != m_runHistList.end(); run_itr++) {
    (*run_itr)->FillHists(weight);
  }

  // ========================================= FAKES ================================ //
  // fake probability(eta in I) = P(reco track is fake | truth eta of reco track is in I)

  //take address of TrackHists object correspondent to all reco tracks
  TrackHists* trackHist_reco_all = m_trackHistList[0];
  Assert("RunHistManager::FillHists()\t trackHist_reco_all mismatched", trackHist_reco_all->GetName().Contains("reco_all"));
  //take address of TrackHists object correspondent to fake reco tracks
  TrackHists* trackHist_reco_fake = m_trackHistList[3];
  Assert("RunHistManager::FillHists()\t trackHist_reco_fake mismatched", trackHist_reco_all->GetName().Contains("reco_all"));

  //we are using true eta
  TH1F* reco_all_abseta = trackHist_reco_all->m_truthAbseta;
  TH1F* reco_fake_abseta = trackHist_reco_fake->m_truthAbseta;

  Assert("RunHistManager::FillHists()\t All and Fake hists incompatible", 
	 reco_all_abseta->GetXaxis()->GetXmin() == reco_fake_abseta->GetXaxis()->GetXmin() &&
	 reco_all_abseta->GetXaxis()->GetXmin() == 0 &&
	 reco_all_abseta->GetXaxis()->GetXmax() == reco_fake_abseta->GetXaxis()->GetXmax() &&
	 reco_all_abseta->GetXaxis()->GetBinWidth(0) == reco_fake_abseta->GetXaxis()->GetBinWidth(0)
);
  
  const double etaInterval = reco_all_abseta->GetXaxis()->GetBinWidth(0);
  int nBin = TMath::Nint((reco_all_abseta->GetXaxis()->GetXmax() - reco_all_abseta->GetXaxis()->GetXmin())/etaInterval);
  
  Assert("RunHistManager::FillHists()\t nBin in fake probability is zero", nBin != 0);

  for (int i = 0; i < nBin; i++) {
    double num = reco_fake_abseta->GetBinContent(i+1);
    double den = reco_all_abseta->GetBinContent(i+1);
    
    if (den == 0) {
      num = 0.0;
      den = 1.0;
    }
   
    double fakeProb = num/den;
    
    m_fakeProb_abseta->SetPoint(i, etaInterval + double(i)*etaInterval, fakeProb);
    m_fakeProb_abseta->SetPointError(i, etaInterval/2.0, TMath::Sqrt(fakeProb*(1-fakeProb)/den));
  }

} // FillHists

