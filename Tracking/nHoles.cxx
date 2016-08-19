#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TPaveText.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "../plotStyle.h"
#include "../GaussFitter.h"


void nHoles(const std::string& pt) {

  std::vector<TFile*> file;

  file.push_back(new TFile(("/afs/cern.ch/user/f/fmassa/newSoshi/TrackingResults/ExtBrl4/datasetPi_" + pt + "GeV/hist-pi" + pt + "pu200.root").c_str()));
  file.push_back(new TFile(("/afs/cern.ch/user/f/fmassa/newSoshi/TrackingResults/ExtBrl4/datasetPi_" + pt + "GeV/hist-pi" + pt + "pu200.root").c_str()));
  file.push_back(new TFile(("/afs/cern.ch/user/f/fmassa/newSoshi/TrackingResults/ExtBrl4/datasetPi_" + pt + "GeV/hist-pi" + pt + "pu200.root").c_str()));
  
  std::vector<std::string> legend;
  legend.push_back("ExtBrl4");
  legend.push_back("IExtBrl4");
  legend.push_back("InclBrl4");


  std::vector<int> color;
  color.push_back(kRed);
  color.push_back(kGreen);
  color.push_back(kBlue);
  
  std::vector<TH1F*> holes;
  holes.push_back((TH1F*) file[0]->Get("TrackHist_reco_all__nPixHoles"));
  holes.push_back((TH1F*) file[1]->Get("TrackHist_reco_all__nPixHoles"));
  holes.push_back((TH1F*) file[2]->Get("TrackHist_reco_all__nPixHoles"));

  std::vector<TGraphErrors*> chiSq;
  chiSq.push_back((TGraphErrors*) file[0]->Get("RunHist_reco_all_hard_all__chiSqPerDof"));
  chiSq.push_back((TGraphErrors*) file[1]->Get("RunHist_reco_all_hard_all__chiSqPerDof"));
  chiSq.push_back((TGraphErrors*) file[2]->Get("RunHist_reco_all_hard_all__chiSqPerDof"));

  std::vector<TH1F*> pull;
  pull.push_back((TH1F*) file[0]->Get("TrackHist_reco_all__qoverpsignif"));
  pull.push_back((TH1F*) file[0]->Get("TrackHist_reco_all__etasignif"));
  pull.push_back((TH1F*) file[0]->Get("TrackHist_reco_all__phisignif"));
  pull.push_back((TH1F*) file[0]->Get("TrackHist_reco_all__d0signif"));
  pull.push_back((TH1F*) file[0]->Get("TrackHist_reco_all__z0signif"));
  
  std::vector<std::string> pullTitle;
  pullTitle.push_back("#frac{q}{p} significance");
  pullTitle.push_back("#eta significance");
  pullTitle.push_back("#phi significance");
  pullTitle.push_back("d_{0} significance");
  pullTitle.push_back("z_{0} significance");

  std::vector<std::string> pullName;
  pullName.push_back("qoverp_pull.png");
  pullName.push_back("eta_pull.png");
  pullName.push_back("phi_pull.png");
  pullName.push_back("d0_pull.png");
  pullName.push_back("z0_pull.png");

  std::cout << "qui" << std::endl;

  ////////////////////

  TCanvas* c = new TCanvas();
  TLegend* leg = new TLegend(0.60,0.85,0.88,0.65);

  TMultiGraph* mg = new TMultiGraph;

  std::cout << "qui" << std::endl;

  for (std::vector<TGraphErrors*>::iterator itr = chiSq.begin();
       itr != chiSq.end();
       itr++) {
    (*itr)->SetMarkerColor(color[itr-chiSq.begin()]);
    (*itr)->SetLineColor(color[itr-chiSq.begin()]);
    (*itr)->SetMarkerSize(2);
    mg->Add(*itr);
    leg->AddEntry(*itr, legend[itr-chiSq.begin()].c_str(), "lp");
  }

  std::cout << "qui" << std::endl;

  mg->Draw("APE");

  mg->SetTitle("");
  mg->GetXaxis()->SetTitle("|\\eta_{truth}|");
  mg->GetYaxis()->SetTitle("\\chi^{2}/DOF");

  leg->Draw();

  c->Print("chiSquarePerDof_eta.png", "png");

  std::cout << "qui" << std::endl;
  delete c;
  delete leg;
  delete mg;

  std::cout << "qui" << std::endl;
  /////////////////////////

  std::cout << "qui" << std::endl;  

  c = new TCanvas();
  
  leg = new TLegend(0.60,0.85,0.88,0.65);

  for (std::vector<TH1F*>::iterator itr = holes.begin();
       itr != holes.end();
       itr++) {
    (*itr)->SetTitle("");
    (*itr)->SetLineWidth(2);
    (*itr)->SetLineColor(color[itr-holes.begin()]);
    (*itr)->Scale(1.0/(*itr)->GetEntries());
    (*itr)->GetYaxis()->SetTitle("A.U.");
    (*itr)->GetXaxis()->SetTitle("Number of pixel holes");
    (*itr)->GetXaxis()->SetRangeUser(0.0,3.0);
    (*itr)->SetStats(0);

    if (itr == holes.begin())
      (*itr)->Draw();
    else
      (*itr)->Draw("same");
    
    leg->AddEntry(*itr, legend[itr-holes.begin()].c_str(), "lp");

  }
  //leg->AddEntry(holes, legend[itr - tracks.begin()].c_str(), "lp");
 
  leg->Draw();
  //leg->Draw();
  (c->cd())->SetLogy();

  c->Print("nPixHoles.png", "png");
       
  delete leg;
  delete c;

  std::cout << "qui" << std::endl;
  ////////////////////////



  for (std::vector<TH1F*>::iterator itr = pull.begin();
       itr != pull.end();
       itr++) {

    c = new TCanvas();

    (*itr)->SetTitle("");
    (*itr)->SetLineWidth(2);
    (*itr)->SetLineColor(kBlue);
    (*itr)->Scale(1.0/(*itr)->GetEntries());
    (*itr)->GetYaxis()->SetTitle("A.U.");
    (*itr)->GetXaxis()->SetTitle(pullTitle[itr-pull.begin()].c_str());
    (*itr)->GetXaxis()->SetRangeUser(-5,5.0);
    (*itr)->SetStats(0);

    GaussFitter fitter((*itr), Verbosity::ERROR, 2);
    fitter.Fit();

    fitter.DrawHistWithFit(c);    

    string mean, meanError, sigma, sigmaError;

    double meanPrecision = -TMath::FloorNint(TMath::Log10(fitter.GetMeanError()));
    double sigmaPrecision = -TMath::FloorNint(TMath::Log10(fitter.GetSigmaError()));

    

    std::stringstream ss;
    ss << TMath::Nint(fitter.GetMean()*TMath::Power(10,meanPrecision))*TMath::Power(10,-meanPrecision);
    ss >> mean;
    ss.clear();

    ss << TMath::Nint(fitter.GetMeanError()*TMath::Power(10,meanPrecision))*TMath::Power(10,-meanPrecision);;
    ss >> meanError;
    ss.clear();

    ss << TMath::Nint(fitter.GetSigma()*TMath::Power(10,sigmaPrecision))*TMath::Power(10,-sigmaPrecision);;
    ss >> sigma;
    ss.clear();

    ss << TMath::Nint(fitter.GetSigmaError()*TMath::Power(10,sigmaPrecision))*TMath::Power(10,-sigmaPrecision);
    ss >> sigmaError;
    ss.clear();

    TPaveText* text;

    if (itr - pull.begin() != 3)
      text = new TPaveText(1.6,0.07,4.5,0.06);
    else
      text = new TPaveText(1.6,0.05,4.5,0.04);

    text->AddText(("Fit mean: " + mean + " #pm " + meanError).c_str());
    text->AddText(("Fit #sigma: " + sigma + " #pm " + sigmaError).c_str());

    text->Draw();

    c->Print(pullName[itr-pull.begin()].c_str(),"png");
    
    delete c;

  }
  //leg->AddEntry(pull, legend[itr - tracks.begin()].c_str(), "lp");
 
  //leg->Draw();


  std::cout << "qui" << std::endl;
}
