#include "HistManager/RunHists.h"
#include "HistManager/TrackHistManager.h"
#include "HistManager/TrackHists.h"

RunHists::RunHists(const TrackHists* const trackHist):m_trackHist(trackHist) {
  m_name = "RunHist_" + trackHist->GetName();
  m_label = "Run ";
}

RunHists::~RunHists() {
  
}

void RunHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;
  m_testHist = declare1D(m_name, "testHist", "x axis", 100,0,100); 

}

void RunHists::FillHists(float weight) const {
  m_testHist->Fill(m_trackHist->m_pt->GetMean());
}
