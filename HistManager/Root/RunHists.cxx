#include "HistManager/RunHists.h"
#include "HistManager/TrackHistManager.h"
#include "HistManager/TrackHists.h"
#include "HistManager/assert.h"

RunHists::RunHists(const TrackHists* const trackHist, const TruthHists* const truthHist) : 
  m_trackHist(trackHist),
  m_truthHist(truthHist) 
{
  m_name = "RunHist_" + trackHist->GetName();
  m_label = "Run ";
}

RunHists::~RunHists() {
  
}

void RunHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;
  m_sysPt_abseta   = declareGraphErrors(m_name, "sysPt_abseta", "|#eta|", "#Delta(p_{T}) bias [MeV]"); 
  m_sysQPt_abseta  = declareGraphErrors(m_name, "sysQPt_abseta", "|#eta|", "p_{T}^{truth} x #Delta(q/p_{T}) bias"); //check if q/pt x pt
  m_sysPhi_abseta  = declareGraphErrors(m_name, "sysPhi_abseta", "|#eta|", "#Delta(#phi) bias"); 
  m_sysD0_abseta   = declareGraphErrors(m_name, "sysD0_abseta", "|#eta|", "#Delta(D0) bias [mm]"); 
  m_sysZ0_abseta   = declareGraphErrors(m_name, "sysZ0_abseta", "|#eta|", "#Delta(Z0) bias [mm]"); 

  m_sigPt_abseta   = declareGraphErrors(m_name, "sigPt_abseta", "|#eta|", "#sigma(p_{T}) [MeV]"); 
  m_sigQPt_abseta  = declareGraphErrors(m_name, "sigQPt_abseta", "|#eta|", "p_{T}^{truth} x #sigma(q/p_{T})"); //check if q/pt x pt
  m_sigPhi_abseta  = declareGraphErrors(m_name, "sigPhi_abseta", "|#eta|", "#sigma(#phi)"); 
  m_sigD0_abseta   = declareGraphErrors(m_name, "sigD0_abseta", "|#eta|", "#sigma(D0) [mm]"); 
  m_sigZ0_abseta   = declareGraphErrors(m_name, "sigZ0_abseta", "|#eta|", "#sigma(Z0) [mm]"); 

  m_eff_abseta     = declareGraphErrors(m_name, "eff_abseta", "|#eta|", "Efficiency"); 
  m_eff_phi        = declareGraphErrors(m_name, "eff_phi", "#phi", "Efficiency"); 

}

void RunHists::FillHists(float weight) const {

  // ================================ Systematic bias and Sigma vs eta ============================ //

  const double m_etaInterval = m_trackHist->m_etaInterval;

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPt_abseta.begin(); itr != m_trackHist->m_biasPt_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPt_abseta.begin();
    m_sysPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_sysPt_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

    m_sigPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetRMS());
    m_sigPt_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetRMSError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasQPt_abseta.begin(); itr != m_trackHist->m_biasQPt_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasQPt_abseta.begin();
    m_sysQPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_sysQPt_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

    m_sigQPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetRMS());
    m_sigQPt_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetRMSError());

  }
  
  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPhi_abseta.begin(); itr != m_trackHist->m_biasPhi_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPhi_abseta.begin();
    m_sysPhi_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_sysPhi_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());
    
    m_sigPhi_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetRMS());
    m_sigPhi_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetRMSError());
    
  }
  
  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasD0_abseta.begin(); itr != m_trackHist->m_biasD0_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasD0_abseta.begin();
    m_sysD0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_sysD0_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

    m_sigD0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetRMS());
    m_sigD0_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetRMSError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasZ0_abseta.begin(); itr != m_trackHist->m_biasZ0_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasZ0_abseta.begin();
    m_sysZ0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_sysZ0_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

    m_sigZ0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetRMS());
    m_sigZ0_abseta->SetPointError(index,m_etaInterval/2.0,(*itr)->GetRMSError());

  }

  // ==================================== Efficiency ========================================== //
  // NOTE: bins doesn't follow usual ROOT convention, they are centered, so bin(0) does not contain any info

  //check if hist of abseta in truth and track hists are compatible
  Assert("Track and Truth abseta histograms are incompatible", 
	 m_trackHist->m_abseta->GetXaxis()->GetBinWidth(1) == m_truthHist->m_abseta->GetXaxis()->GetBinWidth(1) &&
	 m_trackHist->m_abseta->GetXaxis()->GetBinWidth(1) == m_etaInterval &&
	 m_trackHist->m_abseta->GetXaxis()->GetXmin() == m_truthHist->m_abseta->GetXaxis()->GetXmin() &&
	 TMath::Abs(m_trackHist->m_abseta->GetXaxis()->GetXmin()) < 1e-9 &&
	 m_trackHist->m_abseta->GetXaxis()->GetXmax() == m_truthHist->m_abseta->GetXaxis()->GetXmax()
	 );

  int nBins = TMath::Nint((m_trackHist->m_abseta->GetXaxis()->GetXmax() - m_trackHist->m_abseta->GetXaxis()->GetXmin())/m_etaInterval);

  for (int i = 0; i < nBins; i++) {
    Double_t num = m_trackHist->m_abseta->GetBinContent(i+1);
    Double_t den = m_truthHist->m_abseta->GetBinContent(i+1);
    Double_t eff = 0.0;

    if (TMath::Abs(den) > 1e-6) 
      eff = num/den;
    else {
      eff = 0.0;
      den = 1.0;
	}

    m_eff_abseta->SetPoint      (i, m_etaInterval/2.0 + m_etaInterval*double(i), eff);
    m_eff_abseta->SetPointError (i, m_etaInterval/2.0, TMath::Sqrt(eff*(1-eff)/den));
    
  }


//check if hist of phi in truth and track hists are compatible
  Assert("Track and Truth phi histograms are incompatible", 
	 m_trackHist->m_phi->GetXaxis()->GetBinWidth(1) == m_truthHist->m_phi->GetXaxis()->GetBinWidth(1) &&
	 m_trackHist->m_phi->GetXaxis()->GetXmin() == m_truthHist->m_phi->GetXaxis()->GetXmin() &&
	 m_trackHist->m_phi->GetXaxis()->GetXmax() == m_truthHist->m_phi->GetXaxis()->GetXmax()
	 );
  const double m_phiInterval = m_trackHist->m_phi->GetXaxis()->GetBinWidth(0);

  nBins = TMath::Nint((m_trackHist->m_phi->GetXaxis()->GetXmax() - m_trackHist->m_phi->GetXaxis()->GetXmin())/m_phiInterval);

  for (int i = 0; i < nBins; i++) {
    Double_t num = m_trackHist->m_phi->GetBinContent(i+1);
    Double_t den = m_truthHist->m_phi->GetBinContent(i+1);
    Double_t eff = 0.0;

    if (TMath::Abs(den) > 1e-6) 
      eff = num/den;
    else {
      eff = 0.0;
      den = 1.0;
	}

    m_eff_phi->SetPoint      (i, m_trackHist->m_phi->GetXaxis()->GetXmin() + m_phiInterval*double(i), eff);
    m_eff_phi->SetPointError (i, m_phiInterval/2.0, TMath::Sqrt(eff*(1-eff)/den));
    
  }

}
