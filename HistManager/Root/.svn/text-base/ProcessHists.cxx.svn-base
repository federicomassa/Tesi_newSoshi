#include <HistManager/ProcessHists.h>

ProcessHists::ProcessHists(TString name) {
   m_name  = "ProcessHist_" + name;
   m_label = "Process "; // don't forget the space
}

ProcessHists::~ProcessHists() {}

void ProcessHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  m_numLep  = declare1D(m_name, "numLep", "Number of leptons",                4, -0.5, 3.5);
  m_ptLep1  = declare1D(m_name, "ptLep1", "Leading lepton p_{T} [GeV]",     100, 25.0, 225.0);
  m_ptLep2  = declare1D(m_name, "ptLep2", "Sub-leading lepton p_{T} [GeV]", 100, 25.0, 225.0);
  m_etaLep1 = declare1D(m_name, "etaLep1", "Leading lepton #eta",            50, -2.5, 2.5);
  m_etaLep2 = declare1D(m_name, "etaLep2", "Sub-leading lepton #eta",        50, -2.5, 2.5);
  m_phiLep1 = declare1D(m_name, "phiLep1", "Leading lepton #phi",            64, -3.2, 3.2);
  m_phiLep2 = declare1D(m_name, "phiLep2", "Sub-leading lepton #phi",        64, -3.2, 3.2);

  m_massLep12 = declare1D(m_name, "massLep12", "di-lepton mass [GeV]",  100, 25.0, 225.0);
  m_ptLep12   = declare1D(m_name, "ptLep12",   "di-lepton p_{T} [GeV]", 100, 25.0, 225.0);
  m_dRLep12   = declare1D(m_name, "dRLep12",   "dR(lep1,lep2)",          80,  0.0,   8.0);
  m_dEtaLep12 = declare1D(m_name, "dEtaLep12", "d#eta(lep1,lep2)",       40,  0.0,   4.0);

  m_numJet  = declare1D(m_name, "numJet", "Number of jets",           12, -0.5, 11.5);
  m_ptJet1  = declare1D(m_name, "ptJet1", "Leading jet p_{T} [GeV]", 100, 25.0, 525.0);
  m_etaJet1 = declare1D(m_name, "etaJet1", "Leading jet #eta",        50, -5.0, 5.0);
  m_phiJet1 = declare1D(m_name, "phiJet1", "Leading jet #phi",        64, -3.2, 3.2);

} // BookHists

void ProcessHists::FillHists(ProcessFeatures &proc, float weight) const {

  m_numLep  -> Fill(proc.nLeptons,weight);
  if (proc.nLeptons>1) {
    m_ptLep1  -> Fill(proc.Lepton[0].Pt()/1000.0,weight);
    m_ptLep2  -> Fill(proc.Lepton[1].Pt()/1000.0,weight);
    m_etaLep1 -> Fill(proc.Lepton[0].Eta(),weight);
    m_etaLep2 -> Fill(proc.Lepton[1].Eta(),weight);
    m_phiLep1 -> Fill(proc.Lepton[0].Phi(),weight);
    m_phiLep2 -> Fill(proc.Lepton[1].Phi(),weight);

    m_massLep12 -> Fill((proc.Lepton[0]+proc.Lepton[1]).M()/1000.0,weight);
    m_ptLep12   -> Fill((proc.Lepton[0]+proc.Lepton[1]).Pt()/1000.0,weight);
    m_dRLep12   -> Fill(proc.Lepton[0].DeltaR(proc.Lepton[1]),weight);
    m_dEtaLep12 -> Fill(TMath::Abs(proc.Lepton[0].Eta()-proc.Lepton[1].Eta()),weight);
  }

  m_numJet  -> Fill(proc.nJets,weight);
  if (proc.nJets>0) {
    m_ptJet1  -> Fill(proc.Jet[0].Pt()/1000.0,weight);
    m_etaJet1 -> Fill(proc.Jet[0].Eta(),weight);
    m_phiJet1 -> Fill(proc.Jet[0].Phi(),weight);
  }

}

