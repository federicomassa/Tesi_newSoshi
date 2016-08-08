#include <TH1F.h>
#include <TH2F.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>

#include <string>


void TestPlots(const char* rootFileName, const string& histogramTag) {

  const string extension = ".png";
  const char* format = "png";

  TFile* f = new TFile(rootFileName);
  
  TH1F* truthMass = (TH1F*) f->Get("EventHist_all__truthMass");
  TH1F* truthOnShellMass = (TH1F*) f->Get("EventHist_all__truthOnShellMass");
  TH1F* truthOffShellMass = (TH1F*) f->Get("EventHist_all__truthOffShellMass");
  TH2F* truthMassVsTruthOffShellMass = (TH2F*) f->Get("EventHist_all__truthMassVsTruthOffShellMass");

  TH1F* recoMass = (TH1F*) f->Get("EventHist_all__recoMass");
  recoMass->GetXaxis()->SetRangeUser(50,160);

  TH1F* recoOnShellMass = (TH1F*) f->Get("EventHist_all__recoOnShellMass");
  TH1F* recoOffShellMass = (TH1F*) f->Get("EventHist_all__recoOffShellMass");

  TH1F* recoCutCode = (TH1F*) f->Get("EventHist_all__recoCutCode");
  recoCutCode->GetXaxis()->SetRangeUser(-1,21);

  TCanvas* c = new TCanvas();

  ////////////////

  truthMass -> Draw();
  c -> Print((histogramTag + "_" + "truthMass" + extension).c_str(), format);

  truthOnShellMass -> Draw();
  c -> Print((histogramTag + "_" + "truthOnShellMass" + extension).c_str(), format);

  truthOffShellMass -> Draw();
  c -> Print((histogramTag + "_" + "truthOffShellMass" + extension).c_str(), format);

  truthMassVsTruthOffShellMass -> Draw("COL0");
  c -> Print((histogramTag + "_" + "truthMassVsTruthOffShellMass" + extension).c_str(), format);

  recoMass -> Draw();
  c -> Print((histogramTag + "_" + "recoMass" + extension).c_str(), format);

  recoOnShellMass -> Draw();
  c -> Print((histogramTag + "_" + "recoOnShellMass" + extension).c_str(), format);

  recoOffShellMass -> Draw();
  c -> Print((histogramTag + "_" + "recoOffShellMass" + extension).c_str(), format);

  recoCutCode -> Draw();
  c -> Print((histogramTag + "_" + "recoCutCode" + extension).c_str(), format);

  ///////////////

  delete c;

}
