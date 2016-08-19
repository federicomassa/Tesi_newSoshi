#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TMath.h>
#include <THStack.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <iomanip>

#include "GaussFitter.h"

const double cross_H = 50000; //fb

const double BR_HZZ = 0.03;
const double BR_ZZ4mu = 0.03*0.03; //including photons??

const double eff_gamma_H = 0.812;
const double eff_gamma_err_H = 0.001;

//const double totCross_H = cross_H*BR_HZZ*BR_ZZ4mu*eff_gamma_H;
const double totCross_H = 1.64*eff_gamma_H;

const double eff_gen_ZZ = 0.25;
const double eff_gamma_ZZ = 0.808;
const double eff_gamma_err_ZZ = 0.002;

//const double totCross_ZZ = 25.4*eff_gamma_ZZ*eff_gen_ZZ; //fb, including BR ZZ->4mu
//const double totCross_ZZ = 6.9*eff_gamma_ZZ;
const double totCross_ZZ = 24*eff_gamma_ZZ*0.25;


const double intLum = 3000; //fb^-1

double power_err(double s, double b) {

  return (s*(TMath::Sqrt(s)+TMath::Sqrt(b)) + 2*b*TMath::Sqrt(s))/(2*(s+b)*TMath::Sqrt(s+b));

}

void SB(const std::string& signalFile, const std::string& bkgFile, const std::string& option = "", const std::string& name = "", ofstream* f = NULL) {

  double coverage;
  if (name.find("27") != std::string::npos)
    coverage = 2.7;
  else if (name.find("32") != std::string::npos)
    coverage = 3.2;
  else
    coverage = 4.0;
    

  TFile* f_H  = new TFile(signalFile.c_str());
  TFile* f_ZZ = new TFile(bkgFile.c_str());

  TH1F* truthMass_H  = (TH1F*) f_H->Get("EventHist_all__truthMass");
  TH1F* truthMass_ZZ = (TH1F*) f_ZZ->Get("EventHist_all__truthMass");

  TH1F* tmp_H = (TH1F*) f_H->Get(("EventHist_all__recoMass" + option).c_str());
  //  std::cout << tmp_H->GetNbinsX() << '\t' << tmp_H->GetBinLowEdge(1) << '\t' << tmp_H->GetBinLowEdge(tmp_H->GetNbinsX()+1) << std::endl;

  TH1F* null_H  = new TH1F("null_H", "null_H", tmp_H->GetNbinsX(), tmp_H->GetBinLowEdge(1), tmp_H->GetBinLowEdge(tmp_H->GetNbinsX()+1));

  TH1F* mass_H  = new TH1F("mass_H", "mass_H", tmp_H->GetNbinsX(), tmp_H->GetBinLowEdge(1), tmp_H->GetBinLowEdge(tmp_H->GetNbinsX()+1));

  for (int i = 0; i < tmp_H->GetNbinsX(); i++) {
    mass_H -> SetBinContent(i, tmp_H->GetBinContent(i));
    null_H -> SetBinContent(i, 0.0);
    null_H -> SetBinError(i, tmp_H -> GetBinError(i));
  }

  TH1F* tmp_ZZ = (TH1F*) f_ZZ->Get(("EventHist_all__recoMass" + option).c_str());

  TH1F* mass_ZZ = new TH1F("mass_ZZ", "mass_ZZ", tmp_ZZ->GetNbinsX(), tmp_ZZ->GetBinLowEdge(1), tmp_ZZ->GetBinLowEdge(tmp_ZZ->GetNbinsX()+1));

  TH1F* null_ZZ  = new TH1F("null_ZZ", "null_ZZ", tmp_ZZ->GetNbinsX(), tmp_ZZ->GetBinLowEdge(1), tmp_ZZ->GetBinLowEdge(tmp_ZZ->GetNbinsX()+1));

  for (int i = 0; i < tmp_ZZ->GetNbinsX(); i++) {
    mass_ZZ -> SetBinContent(i, tmp_ZZ->GetBinContent(i));
    null_ZZ -> SetBinContent(i, 0.0);
    null_ZZ -> SetBinError(i, tmp_ZZ -> GetBinError(i));
  }


  const double eff_analysis_H = tmp_H->GetEntries()/truthMass_H->GetEntries();
  const double eff_analysis_err_H = TMath::Sqrt(eff_analysis_H*(1-eff_analysis_H)/truthMass_H->GetEntries());

  const double eff_analysis_ZZ = tmp_ZZ->GetEntries()/truthMass_ZZ->GetEntries();
  const double eff_analysis_err_ZZ = TMath::Sqrt(eff_analysis_ZZ*(1-eff_analysis_ZZ)/truthMass_ZZ->GetEntries());

  THStack* massStack = new THStack("massStack", "");


  double scale_H = totCross_H*eff_analysis_H*intLum/tmp_H->GetEntries();
  double scale_ZZ = totCross_ZZ*eff_analysis_ZZ*intLum/tmp_ZZ->GetEntries();

  double scale_err_H = scale_H*TMath::Sqrt(TMath::Power(eff_gamma_err_H/eff_gamma_H,2) + TMath::Power(eff_analysis_err_H/eff_analysis_H,2));
  double scale_err_ZZ = scale_ZZ*TMath::Sqrt(TMath::Power(eff_gamma_err_ZZ/eff_gamma_ZZ,2) + TMath::Power(eff_analysis_err_ZZ/eff_analysis_ZZ,2));

  //  std::cout << tmp_H->GetEntries() << '\t' << tmp_ZZ->GetEntries() << std::endl;

  mass_ZZ -> Scale(scale_ZZ);
  mass_ZZ -> SetFillColor(kBlue);
  mass_ZZ -> SetMarkerColor(kBlue);
  mass_ZZ -> SetMarkerStyle(22);

  massStack -> Add(mass_ZZ);

  mass_H -> SetFillColor(kRed);
  mass_H -> SetFillStyle(3001);
  mass_H -> SetMarkerColor(kRed);
  mass_H -> SetMarkerStyle(21);
  mass_H -> SetMarkerSize(0.5);
  mass_H->Scale(scale_H);

  TCanvas* c = new TCanvas();

  GaussFitter fitter(mass_H);
  fitter.Fit();
  

  for (int i = 0; i < tmp_ZZ->GetNbinsX(); i++)
    mass_H -> SetBinError(i, scale_err_H); //TODO: errore sulla scala?

  massStack -> Add(mass_H, "e2lp");
  //  massStack -> Add(null_hist, "elp");

  TPad* p1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  p1->SetBottomMargin(0.02);
  p1->Draw();
  p1->cd();

  massStack->Draw();
  massStack -> GetXaxis() -> SetRangeUser(100,150);
  massStack -> GetYaxis() -> SetTitle("Events at 3000 fb^{-1} / 1 GeV");
  massStack -> GetXaxis() -> SetLabelSize(0.);

  massStack->Draw();

  TLegend* leg = new TLegend(0.15,0.85,0.35,0.70);

  leg->AddEntry(mass_H, "Higgs + ZZ", "lp");
  leg->AddEntry(mass_ZZ, "ZZ", "lp");

  leg->Draw();

  
  TPaveText* text = new TPaveText(135, 500, 148, 440);
  std::stringstream ss;
  ss << coverage;
  std::string coverage_str;
  ss >> coverage_str;
  
  text->AddText(("|#eta^{#mu}_{max}| < " + coverage_str).c_str());
  text->Draw();
  text->SetTextSize(0.8*text->GetTextSize());
  text->Draw();

  c->cd();

  TPad* p2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  

  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.2);
  //p2->SetLeftMargin(0.4);
  p2->Draw();
  p2->cd();
  

  TF1* fitFunc = fitter.DrawFitFunction();

  for (int i = 0; i < tmp_H->GetNbinsX(); i++)
    mass_H->SetBinError(i, 0.0);
  
  mass_H->SetTitle("");
  mass_H->SetStats(0);
  mass_H->GetXaxis()->SetLabelFont(43);
  mass_H->GetXaxis()->SetLabelSize(15);
  mass_H->GetYaxis()->SetLabelSize(0);
  mass_H->DrawCopy("AH");
  mass_H->GetXaxis()->SetRangeUser(100,150);
  mass_H->DrawCopy("AH");

  TGaxis* axis = new TGaxis(100,0, 100, 540, 0, 540, 502, "R");
  axis->SetLabelFont(43);
  axis->SetLabelSize(12);
  axis->SetTitle("");
  axis->SetTitleSize(0.05);
  axis->SetTitleOffset(1.1);
  p2->cd();
  axis->Draw();

  TGaxis* axis2 = new TGaxis(100,0, 150, 0, 100, 150, 10+500, "R");
  axis2->SetLabelFont(43);
  axis2->SetLabelSize(12);
  axis2->SetTitle("Four-muon invariant mass [GeV]");
  axis2->SetTitleSize(0.09);
  axis2->SetTitleOffset(1.2);
  axis2->SetTitleFont(40);
  //axis2->SetTitleOffset(1.01);
  //  axis2->SetTickSize(0.03);
  p2->cd();
  axis2->Draw();


  double lowEdge = fitter.GetMean() - 1.5*fitter.GetSigma();
  double upperEdge = fitter.GetMean() + 1.5*fitter.GetSigma();

  double integral_H = mass_H->Integral(mass_H->FindBin(lowEdge),
				       mass_H->FindBin(upperEdge));

  double integral_ZZ = mass_ZZ->Integral(mass_ZZ->FindBin(lowEdge),
				       mass_ZZ->FindBin(upperEdge));

  std::cout << "==================== FIT RESULT ====================" << std::endl;

  std::cout << "Amplitude: " << fitter.GetAmplitude() << " +/- " << fitter.GetAmplitudeError() << std::endl;
  std::cout << "Mean:      " << fitter.GetMean() << " +/- " << fitter.GetMeanError() << std::endl;
  std::cout << "Sigma:     " << fitter.GetSigma() << " +/- " << fitter.GetSigmaError() << std::endl;

  std::cout << "====================================================" << std::endl;


  std::cout << "Scale H:  " << scale_H << std::endl;
  std::cout << "Scale ZZ: " << scale_ZZ << std::endl;

  std::cout << "S: " << integral_H << " +/- " << TMath::Sqrt(integral_H) << " +/- " << integral_H*scale_err_H/scale_H << " (sys) +/- " <<
    TMath::Sqrt(integral_H*scale_H) << " (extr. stat)" << std::endl;

  if (f != NULL) {

    *f << TMath::Nint(integral_H) << " $\\pm$ " << TMath::Nint(TMath::Sqrt(integral_H)) << " & ";

  }

  std::cout << "Fit S: " << fitFunc->Integral(lowEdge, upperEdge) << std::endl;

  std::cout << "B: " << integral_ZZ << " +/- " << TMath::Sqrt(integral_ZZ) << " +/- " << integral_ZZ*scale_err_ZZ/scale_ZZ << " (sys) +/- " <<
    TMath::Sqrt(integral_ZZ*scale_ZZ) << " (extr. stat)" << std::endl;

  if (f != NULL)
    *f << TMath::Nint(integral_ZZ) << " $\\pm$ " << TMath::Nint(TMath::Sqrt(integral_ZZ)) << " & ";

  std::cout << "S+B: " << integral_H + integral_ZZ << " +/- " << TMath::Sqrt(integral_H+integral_ZZ) << " (stat 3000 fb^-1) +/- " <<
    TMath::Sqrt(TMath::Power(integral_H*scale_err_H/scale_H,2) + TMath::Power(integral_ZZ*scale_err_ZZ/scale_ZZ,2)) << " (sys) " << 
    TMath::Sqrt(integral_H*scale_H + integral_ZZ*scale_ZZ) << " (extr. stat) " << std::endl;

  std::cout << "S/err(S+B): " << integral_H/TMath::Sqrt(TMath::Power(integral_H*scale_err_H/scale_H,2) + TMath::Power(integral_ZZ*scale_err_ZZ/scale_ZZ,2)) << std::endl;


  std::cout << "S/B: " << integral_H/integral_ZZ << std::endl;

  if (f != NULL)
    *f << std::fixed << std::setprecision(1) << integral_H/integral_ZZ << " & ";

  std::cout << "S/sqrt(S+B): " << integral_H/TMath::Sqrt(integral_H + integral_ZZ) << std::endl;
  
  if (f != NULL) {

    int precision = -TMath::FloorNint(TMath::Log10(power_err(integral_H, integral_ZZ)));
    *f << std::fixed << std::setprecision(precision) << integral_H/TMath::Sqrt(integral_H + integral_ZZ) << " $\\pm$ " << power_err(integral_H,integral_ZZ) << " \\\\ \\hline" << std::endl;
  }

  fitFunc->Draw("same");


  c->Print(("SB"+name+".png").c_str(), "png");


}
