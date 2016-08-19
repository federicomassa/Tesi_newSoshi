#include "HistManager/GaussFitter.h"
#include "HistManager/RunHists.h"
#include "HistManager/TrackHistManager.h"
#include "HistManager/TrackHists.h"
#include "Utility/assert.h"
#include "HistManager/PrintMessage.h"

const Verbosity verbosity_level = Verbosity::ERROR;

RunHists::RunHists(const TrackHists* const trackHist, const TruthHists* const truthHist) : 
  m_trackHist(trackHist),
  m_truthHist(truthHist) 
{
  m_name = "RunHist_" + trackHist->GetName() + "_" + truthHist->GetName();
  m_label = "Run ";
}

RunHists::~RunHists() {
  
}

void RunHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  m_sysPt_abseta   = declareGraphErrors(m_name, "sysPt_abseta", "|#eta|", "#Delta(p_{T}) bias [MeV]"); 
  m_sysQPt_abseta  = declareGraphErrors(m_name, "sysQPt_abseta", "|#eta|", "p_{T}^{truth} x #Delta(q/p_{T}) bias");
  m_sysPhi_abseta  = declareGraphErrors(m_name, "sysPhi_abseta", "|#eta|", "#Delta(#phi) bias"); 
  m_sysD0_abseta   = declareGraphErrors(m_name, "sysD0_abseta", "|#eta|", "#Delta(D0) bias [mm]"); 
  m_sysZ0_abseta   = declareGraphErrors(m_name, "sysZ0_abseta", "|#eta|", "#Delta(Z0) bias [mm]"); 
  m_sysEta_abseta  = declareGraphErrors(m_name, "sysEta_abseta", "|#eta|", "#Delta(#eta)/|#eta_{truth}| bias");

  m_sigPt_abseta   = declareGraphErrors(m_name, "sigPt_abseta", "|#eta|", "#sigma(p_{T}) [MeV]"); 
  m_sigQPt_abseta  = declareGraphErrors(m_name, "sigQPt_abseta", "|#eta|", "p_{T}^{truth} x #sigma(q/p_{T})");
  m_sigEta_abseta  = declareGraphErrors(m_name, "sigEta_abseta", "|#eta|", "#sigma(#eta)/|#eta_{truth}|");
  m_sigPhi_abseta  = declareGraphErrors(m_name, "sigPhi_abseta", "|#eta|", "#sigma(#phi)"); 
  m_sigD0_abseta   = declareGraphErrors(m_name, "sigD0_abseta", "|#eta|", "#sigma(D0) [mm]"); 
  m_sigZ0_abseta   = declareGraphErrors(m_name, "sigZ0_abseta", "|#eta|", "#sigma(Z0) [mm]"); 

  m_sysPt_poseta   = declareGraphErrors(m_name, "sysPt_poseta", "|#eta|", "#Delta(p_{T}) bias [MeV]: positive #eta"); 
  m_sysQPt_poseta  = declareGraphErrors(m_name, "sysQPt_poseta", "|#eta|", "p_{T}^{truth} x #Delta(q/p_{T}) bias: positive #eta");
  m_sysPhi_poseta  = declareGraphErrors(m_name, "sysPhi_poseta", "|#eta|", "#Delta(#phi) bias: positive #eta"); 
  m_sysD0_poseta   = declareGraphErrors(m_name, "sysD0_poseta", "|#eta|", "#Delta(D0) bias [mm]: positive #eta"); 
  m_sysZ0_poseta   = declareGraphErrors(m_name, "sysZ0_poseta", "|#eta|", "#Delta(Z0) bias [mm]: positive #eta"); 
  m_sysEta_poseta  = declareGraphErrors(m_name, "sysEta_poseta", "|#eta|", "#Delta(#eta)/|#eta_{truth}| bias: positive #eta");

  m_sigPt_poseta   = declareGraphErrors(m_name, "sigPt_poseta", "|#eta|", "#sigma(p_{T}) [MeV]: positive #eta"); 
  m_sigQPt_poseta  = declareGraphErrors(m_name, "sigQPt_poseta", "|#eta|", "p_{T}^{truth} x #sigma(q/p_{T}): positive #eta");
  m_sigEta_poseta  = declareGraphErrors(m_name, "sigEta_poseta", "|#eta|", "#sigma(#eta)/|#eta_{truth}|: positive #eta");
  m_sigPhi_poseta  = declareGraphErrors(m_name, "sigPhi_poseta", "|#eta|", "#sigma(#phi): positive #eta"); 
  m_sigD0_poseta   = declareGraphErrors(m_name, "sigD0_poseta", "|#eta|", "#sigma(D0) [mm]: positive #eta"); 
  m_sigZ0_poseta   = declareGraphErrors(m_name, "sigZ0_poseta", "|#eta|", "#sigma(Z0) [mm]: positive #eta"); 

  m_sysPt_negeta   = declareGraphErrors(m_name, "sysPt_negeta", "|#eta|", "#Delta(p_{T}) bias [MeV]: negative #eta"); 
  m_sysQPt_negeta  = declareGraphErrors(m_name, "sysQPt_negeta", "|#eta|", "p_{T}^{truth} x #Delta(q/p_{T}) bias: negative #eta");
  m_sysPhi_negeta  = declareGraphErrors(m_name, "sysPhi_negeta", "|#eta|", "#Delta(#phi) bias: negative #eta"); 
  m_sysD0_negeta   = declareGraphErrors(m_name, "sysD0_negeta", "|#eta|", "#Delta(D0) bias [mm]: negative #eta"); 
  m_sysZ0_negeta   = declareGraphErrors(m_name, "sysZ0_negeta", "|#eta|", "#Delta(Z0) bias [mm]: negative #eta"); 
  m_sysEta_negeta  = declareGraphErrors(m_name, "sysEta_negeta", "|#eta|", "#Delta(#eta)/|#eta_{truth}| bias: negative #eta");

  m_sigPt_negeta   = declareGraphErrors(m_name, "sigPt_negeta", "|#eta|", "#sigma(p_{T}) [MeV]: negative #eta"); 
  m_sigQPt_negeta  = declareGraphErrors(m_name, "sigQPt_negeta", "|#eta|", "p_{T}^{truth} x #sigma(q/p_{T}): negative #eta");
  m_sigEta_negeta  = declareGraphErrors(m_name, "sigEta_negeta", "|#eta|", "#sigma(#eta)/|#eta_{truth}|: negative #eta");
  m_sigPhi_negeta  = declareGraphErrors(m_name, "sigPhi_negeta", "|#eta|", "#sigma(#phi): negative #eta"); 
  m_sigD0_negeta   = declareGraphErrors(m_name, "sigD0_negeta", "|#eta|", "#sigma(D0) [mm]: negative #eta"); 
  m_sigZ0_negeta   = declareGraphErrors(m_name, "sigZ0_negeta", "|#eta|", "#sigma(Z0) [mm]: negative #eta"); 

  m_PixHits        = declareGraphErrors(m_name, "PixHits", "|#eta|", "<PixHits>");
  m_SCTHits        = declareGraphErrors(m_name, "SCTHits", "|#eta|", "<SCTHits>");
  m_SiHits        = declareGraphErrors(m_name, "SiHits", "|#eta|", "<SiHits>");

  m_chiSqPerDof   = declareGraphErrors(m_name, "chiSqPerDof", "|#eta|", "<#chi^{2}/DOF>");


  m_eff_abseta     = declareGraphErrors(m_name, "eff_abseta", "|#eta|", "Efficiency"); 
  m_eff_phi        = declareGraphErrors(m_name, "eff_phi", "#phi", "Efficiency"); 

}

void RunHists::FillHists(float weight) const {

  // ================================ Systematic bias and Sigma vs eta ============================ //

  const double m_etaInterval = m_trackHist->m_etaInterval;

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPt_abseta.begin(); itr != m_trackHist->m_biasPt_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPt_abseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysPt_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigPt_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasQPt_abseta.begin(); itr != m_trackHist->m_biasQPt_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasQPt_abseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysQPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysQPt_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigQPt_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigQPt_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasEta_abseta.begin(); itr != m_trackHist->m_biasEta_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasEta_abseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysEta_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysEta_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigEta_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigEta_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPhi_abseta.begin(); itr != m_trackHist->m_biasPhi_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPhi_abseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysPhi_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysPhi_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigPhi_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigPhi_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasD0_abseta.begin(); itr != m_trackHist->m_biasD0_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasD0_abseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysD0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysD0_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigD0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigD0_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasZ0_abseta.begin(); itr != m_trackHist->m_biasZ0_abseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasZ0_abseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysZ0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysZ0_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigZ0_abseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigZ0_abseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  // ================================ Systematic bias and Sigma vs eta for positive eta ============================ //

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPt_poseta.begin(); itr != m_trackHist->m_biasPt_poseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPt_poseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysPt_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysPt_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigPt_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigPt_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasQPt_poseta.begin(); itr != m_trackHist->m_biasQPt_poseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasQPt_poseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysQPt_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysQPt_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigQPt_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigQPt_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasEta_poseta.begin(); itr != m_trackHist->m_biasEta_poseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasEta_poseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysEta_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysEta_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigEta_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigEta_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPhi_poseta.begin(); itr != m_trackHist->m_biasPhi_poseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPhi_poseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysPhi_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysPhi_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigPhi_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigPhi_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasD0_poseta.begin(); itr != m_trackHist->m_biasD0_poseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasD0_poseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysD0_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysD0_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigD0_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigD0_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasZ0_poseta.begin(); itr != m_trackHist->m_biasZ0_poseta.end(); itr++) {
    int index = itr - m_trackHist->m_biasZ0_poseta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysZ0_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysZ0_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigZ0_poseta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigZ0_poseta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  // ================================ Systematic bias and Sigma vs eta for negative eta ============================ //

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPt_negeta.begin(); itr != m_trackHist->m_biasPt_negeta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPt_negeta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysPt_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysPt_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigPt_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigPt_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasQPt_negeta.begin(); itr != m_trackHist->m_biasQPt_negeta.end(); itr++) {
    int index = itr - m_trackHist->m_biasQPt_negeta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysQPt_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysQPt_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigQPt_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigQPt_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasEta_negeta.begin(); itr != m_trackHist->m_biasEta_negeta.end(); itr++) {
    int index = itr - m_trackHist->m_biasEta_negeta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysEta_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysEta_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigEta_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigEta_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasPhi_negeta.begin(); itr != m_trackHist->m_biasPhi_negeta.end(); itr++) {
    int index = itr - m_trackHist->m_biasPhi_negeta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysPhi_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysPhi_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigPhi_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigPhi_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasD0_negeta.begin(); itr != m_trackHist->m_biasD0_negeta.end(); itr++) {
    int index = itr - m_trackHist->m_biasD0_negeta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysD0_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysD0_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigD0_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigD0_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_biasZ0_negeta.begin(); itr != m_trackHist->m_biasZ0_negeta.end(); itr++) {
    int index = itr - m_trackHist->m_biasZ0_negeta.begin();

    GaussFitter fitter(*itr, verbosity_level);
    fitter.Fit();
    m_sysZ0_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetMean());
    m_sysZ0_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetMeanError());

    m_sigZ0_negeta->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), fitter.GetSigma());
    m_sigZ0_negeta->SetPointError(index,m_etaInterval/2.0,fitter.GetSigmaError());

  }


  // ==================================== Efficiency ========================================== //
  // NOTE: bins doesn't follow usual ROOT convention, they are centered, so bin(0) does not contain any info
  // Efficiency(x in I) = P(track is reconstructed| truth has x in I).. criteria for reconstruction is specified by which
  //                                                                    TrackHist_reco_? we choose..
  //                                                                    classical efficiency would be
  //                                                                    reco_true/truth_stable
  //                                                                    chosen when plotting

  //check if hist of abseta in truth and track hists are compatible
  Assert("Track and Truth abseta histograms are incompatible", 
	 m_trackHist->m_truthAbseta->GetXaxis()->GetBinWidth(1) == m_truthHist->m_abseta->GetXaxis()->GetBinWidth(1) &&
	 m_trackHist->m_truthAbseta->GetXaxis()->GetBinWidth(1) == m_etaInterval &&
	 m_trackHist->m_truthAbseta->GetXaxis()->GetXmin() == m_truthHist->m_abseta->GetXaxis()->GetXmin() &&
	 TMath::Abs(m_trackHist->m_truthAbseta->GetXaxis()->GetXmin()) < 1e-9 &&
	 m_trackHist->m_truthAbseta->GetXaxis()->GetXmax() == m_truthHist->m_abseta->GetXaxis()->GetXmax()
	 );

  int nBins = TMath::Nint((m_trackHist->m_truthAbseta->GetXaxis()->GetXmax() - m_trackHist->m_truthAbseta->GetXaxis()->GetXmin())/m_etaInterval);

  for (int i = 0; i < nBins; i++) {
    Double_t num = m_trackHist->m_truthAbseta->GetBinContent(i+1);
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
	 m_trackHist->m_truthPhi->GetXaxis()->GetBinWidth(1) == m_truthHist->m_phi->GetXaxis()->GetBinWidth(1) &&
	 m_trackHist->m_truthPhi->GetXaxis()->GetXmin() == m_truthHist->m_phi->GetXaxis()->GetXmin() &&
	 m_trackHist->m_truthPhi->GetXaxis()->GetXmax() == m_truthHist->m_phi->GetXaxis()->GetXmax()
	 );
  const double m_phiInterval = m_trackHist->m_truthPhi->GetXaxis()->GetBinWidth(0);

  nBins = TMath::Nint((m_trackHist->m_truthPhi->GetXaxis()->GetXmax() - m_trackHist->m_truthPhi->GetXaxis()->GetXmin())/m_phiInterval);

  for (int i = 0; i < nBins; i++) {
    Double_t num = m_trackHist->m_truthPhi->GetBinContent(i+1);
    Double_t den = m_truthHist->m_phi->GetBinContent(i+1);
    Double_t eff = 0.0;

    if (TMath::Abs(den) > 1e-6) 
      eff = num/den;
    else {
      eff = 0.0;
      den = 1.0;
	}

    m_eff_phi->SetPoint      (i, m_trackHist->m_truthPhi->GetXaxis()->GetXmin() + m_phiInterval*double(i), eff);
    m_eff_phi->SetPointError (i, m_phiInterval/2.0, TMath::Sqrt(eff*(1-eff)/den));
    
  }


  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_PixHits.begin(); itr != m_trackHist->m_PixHits.end(); itr++) {
    int index = itr - m_trackHist->m_PixHits.begin();

    m_PixHits->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_PixHits->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_SCTHits.begin(); itr != m_trackHist->m_SCTHits.end(); itr++) {
    int index = itr - m_trackHist->m_SCTHits.begin();

    m_SCTHits->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_SCTHits->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_SiHits.begin(); itr != m_trackHist->m_SiHits.end(); itr++) {
    int index = itr - m_trackHist->m_SiHits.begin();

    m_SiHits->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_SiHits->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

  }

  for (std::vector<TH1F*>::const_iterator itr = m_trackHist->m_chiSqPerDof2D.begin(); itr != m_trackHist->m_chiSqPerDof2D.end(); itr++) {
    int index = itr - m_trackHist->m_chiSqPerDof2D.begin();

    m_chiSqPerDof->SetPoint(index, m_etaInterval/2.0 + m_etaInterval*double(index), (*itr)->GetMean());
    m_chiSqPerDof->SetPointError(index,m_etaInterval/2.0,(*itr)->GetMeanError());

  }





}
