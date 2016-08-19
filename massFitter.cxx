#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "plotStyle.h"
#include "GaussFitter.h"

using namespace std;

void massFitter(const char* fname1, const char* fname2, const char* fname3, const char* latexOut) {

  ofstream latex(latexOut);

  std::vector<TFile*> f;
  f.push_back(new TFile(fname1));
  f.push_back(new TFile(fname2));
  f.push_back(new TFile(fname3));

  std::vector<TH1F*> onShell, offShell, tot, pt, eta, phi,
    onShell27, offShell27, tot27,
    onShell32, offShell32, tot32,
    onShell4, offShell4, tot4,
    pt27, pt32, pt4,
    eta27, eta32, eta4,
    phi27, phi32, phi4;
  

  for (unsigned int i = 0; i < 3; i++) {

    onShell.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass"));
    offShell.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass"));
    tot.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass"));
    pt.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPt"));
    eta.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoEta"));
    phi.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPhi"));
    

    onShell27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass27"));
    offShell27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass27"));
    tot27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass27"));
    pt27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPt27"));
    eta27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoEta27"));
    phi27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPhi27"));


    onShell32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass32"));
    offShell32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass32"));
    tot32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass32"));
    pt32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPt32"));
    eta32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoEta32"));
    phi32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPhi32"));


    onShell4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass4"));
    offShell4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass4"));
    tot4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass4"));
    pt4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPt4"));
    eta4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoEta4"));
    phi4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoPhi4"));


  }

  std::vector<std::pair<std::vector<TH1F*>, const string > > v;
  v.push_back(std::make_pair(onShell, "\\textbf{On-shell mass (GeV)}"));
  v.push_back(std::make_pair(offShell, "\\textbf{Off-shell mass  (GeV)}"));
  v.push_back(std::make_pair(tot, "\\textbf{$4\\mu$ mass (GeV)}"));
  v.push_back(std::make_pair(pt, "\\textbf{$4\\mu$ $p_{T}$ (GeV)}"));
  v.push_back(std::make_pair(eta, "\\textbf{$4\\mu$ $\\eta$}"));
  v.push_back(std::make_pair(phi, "\\textbf{$4\\mu$ $\\phi$}"));


  v.push_back(std::make_pair(onShell27, "\\textbf{On-shell  - 2.7 (GeV)}"));
  v.push_back(std::make_pair(offShell27, "\\textbf{Off-shell  - 2.7 (GeV)}"));
  v.push_back(std::make_pair(tot27, "\\textbf{$4\\mu$  - 2.7 (GeV)}"));
  v.push_back(std::make_pair(pt27, "\\textbf{$4\\mu$ $p_{T}$ - 2.7 (GeV)}"));
  v.push_back(std::make_pair(eta27, "\\textbf{$4\\mu$ $\\eta$} - 2.7"));
  v.push_back(std::make_pair(phi27, "\\textbf{$4\\mu$ $\\phi$} - 2.7"));


  v.push_back(std::make_pair(onShell32, "\\textbf{On-shell  - 3.2 (GeV)}"));
  v.push_back(std::make_pair(offShell32, "\\textbf{Off-shell  - 3.2 (GeV)}"));
  v.push_back(std::make_pair(tot32, "\\textbf{$4\\mu$  - 3.2 (GeV)}"));
  v.push_back(std::make_pair(pt32, "\\textbf{$4\\mu$ $p_{T}$ - 3.2 (GeV)}"));
  v.push_back(std::make_pair(eta32, "\\textbf{$4\\mu$ $\\eta$} - 3.2"));
  v.push_back(std::make_pair(phi32, "\\textbf{$4\\mu$ $\\phi$} - 3.2"));


  v.push_back(std::make_pair(onShell4, "\\textbf{On-shell  - 4.0 (GeV)}"));
  v.push_back(std::make_pair(offShell4, "\\textbf{Off-shell  - 4.0 (GeV)}"));
  v.push_back(std::make_pair(tot4, "\\textbf{$4\\mu$  - 4.0 (GeV)}"));
  v.push_back(std::make_pair(pt4, "\\textbf{$4\\mu$ $p_{T}$ - 4.0 (GeV)}"));
  v.push_back(std::make_pair(eta4, "\\textbf{$4\\mu$ $\\eta$} - 4.0"));
  v.push_back(std::make_pair(phi4, "\\textbf{$4\\mu$ $\\phi$} - 4.0"));


  for (unsigned int i = 0; i < v.size(); i++) {

    latex << v[i].second << " & ";

    for (int j = 0; j < 3; j++) {

      GaussFitter fitter((v[i].first)[j], Verbosity::ERROR, 3);

      fitter.Fit();

      if (i == 0 && j == 1) {

	TCanvas* c = new TCanvas();
	fitter.DrawHistWithFit(c);

      }
      
      std::cout << fitter.GetReducedChiSquare() << std::endl;

      int precision = -TMath::FloorNint(TMath::Log10(fitter.GetSigmaError()));

      latex << std::fixed << std::setprecision(precision) << fitter.GetSigma() <<
	" $\\pm$ "  << fitter.GetSigmaError();

      if (j == 2)
	latex << " \\\\ \\hline";
      else
	latex << " & ";
      
    if (j == 2 && (i+1)%3 == 0)
      latex << " \\hline ";

    latex << std::endl;      


    }


  }
  
  latex.close();

}
