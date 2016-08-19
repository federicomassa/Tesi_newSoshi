#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "../plotStyle.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

void nHits(const std::string& layout, double* pixMean = 0, double* pixMeanError = 0, double* SCTMean = 0, double* SCTMeanError = 0, 
	   double* siMean = 0, double* siMeanError = 0) {

  std::vector<TFile*> file;

  file.push_back(new TFile(("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_" + layout + "/hist-pi15pu0.root").c_str()));
  file.push_back(new TFile(("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_" + layout + "/hist-pi15pu50.root").c_str()));
  file.push_back(new TFile(("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_" + layout + "/hist-pi15pu140.root").c_str()));
  file.push_back(new TFile(("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_" + layout + "/hist-pi15pu200.root").c_str()));

  std::vector<std::string> legend;
  legend.push_back("<#mu> = 0");
  legend.push_back("<#mu> = 50");
  legend.push_back("<#mu> = 140");
  legend.push_back("<#mu> = 200");


  std::vector<int> color;
  color.push_back(kRed);
  color.push_back(kBlue);
  color.push_back(kGreen);
  color.push_back(kBlack);
  
  TH2F* pixel = (TH2F*) file[3]->Get("TrackHist_reco_all__eta_nPixHits");

  TH2F* SCT = (TH2F*) file[3]->Get("TrackHist_reco_all__eta_nSCTHits");

  TH2F* total = (TH2F*) file[3]->Get("TrackHist_reco_all__eta_nSiHits");

  TH2F* holes = (TH2F*) file[3]->Get("TrackHist_reco_all__nPixHoles");


  /////////////////////////

  TCanvas *c = new TCanvas();
  
  TLegend* leg = new TLegend(0.60,0.85,0.88,0.65);

  if (pixMean && pixMeanError) {
    *pixMean = pixel->GetMean(2);
    *pixMeanError = pixel->GetMeanError(2);
  }    

  for (int i = 0; i < pixel->GetNbinsX(); i++) {
    
    double sum = 0;
    
    for (int j = 0; j < pixel->GetNbinsY(); j++)
      sum += pixel->GetBinContent(i+1, j+1);

    for (int j = 0; j < pixel->GetNbinsY(); j++) {
      if (sum > 1E-6)
	pixel->SetBinContent(i+1,j+1,pixel->GetBinContent(i+1,j+1)/sum);
      else
	pixel->SetBinContent(i+1,j+1,0);
    }
    
  }
  
  pixel->SetTitle("");
  pixel->GetXaxis()->SetTitle("Track |#eta|");
  pixel->GetYaxis()->SetTitle("Number of pixel hits");
  pixel->SetStats(0);
  pixel->Draw("COLZ");
  
  //leg->AddEntry(pixel, legend[itr - tracks.begin()].c_str(), "lp");
  

    

  //leg->Draw();

  c->Print(("nPixHits_eta_" + layout + ".png").c_str(), "png");

       
  delete leg;
  delete c;

  /////////////////////////

  c = new TCanvas();
  
  leg = new TLegend(0.60,0.85,0.88,0.65);

  if (SCTMean && SCTMeanError) {
    *SCTMean = SCT->GetMean(2);
    *SCTMeanError = SCT->GetMeanError(2);
  }

  for (int i = 0; i < SCT->GetNbinsX(); i++) {
    
    double sum = 0;
    
    for (int j = 0; j < SCT->GetNbinsY(); j++)
      sum += SCT->GetBinContent(i+1, j+1);

    for (int j = 0; j < SCT->GetNbinsY(); j++) {
      if (sum > 1E-6)
	SCT->SetBinContent(i+1,j+1,SCT->GetBinContent(i+1,j+1)/sum);
      else
    	SCT->SetBinContent(i+1,j+1,0);
    }
    
  }
  
  SCT->SetTitle("");
  SCT->GetXaxis()->SetTitle("Track |#eta|");
  SCT->GetYaxis()->SetTitle("Number of SCT hits");
  SCT->SetStats(0);
  SCT->Draw("COLZ");
  
  //leg->AddEntry(SCT, legend[itr - tracks.begin()].c_str(), "lp");
  
  //leg->Draw();
  

  c->Print(("nSCTHits_eta_" + layout + ".png").c_str(), "png");
  
  
  delete leg;
  delete c;
  
  /////////////////////////
  
  c = new TCanvas();
  
  leg = new TLegend(0.60,0.85,0.88,0.65);
  
  if (siMean && siMeanError) {
    *siMean = total->GetMean(2);
    *siMeanError = total->GetMeanError(2);
  }
  
  for (int i = 0; i < total->GetNbinsX(); i++) {
    
    double sum = 0;
    
    for (int j = 0; j < total->GetNbinsY(); j++)
      sum += total->GetBinContent(i+1, j+1);
    
    for (int j = 0; j < total->GetNbinsY(); j++) {
      if (sum > 1E-6)
	total->SetBinContent(i+1,j+1,total->GetBinContent(i+1,j+1)/sum);
      else 
	total->SetBinContent(i+1,j+1,0);
    }
    
  }
  
  total->SetTitle("");
  total->GetXaxis()->SetTitle("Track |#eta|");
  total->GetYaxis()->SetTitle("Number of total hits");
  total->SetStats(0);
  total->Draw("COLZ");
  
  //leg->AddEntry(total, legend[itr - tracks.begin()].c_str(), "lp");
 

  //leg->Draw();

  c->Print(("nSiHits_eta_" + layout + ".png").c_str(), "png");


       
  delete leg;
  delete c;

  /////////////////////////
  
  c = new TCanvas();
  
  leg = new TLegend(0.60,0.85,0.88,0.65);

  holes->SetTitle("");
  holes->SetLineWidth(2);
  holes->SetLineColor(kBlue);
  holes->Scale(1.0/holes->GetEntries());
  holes->GetYaxis()->SetTitle("A.U.");
  holes->GetXaxis()->SetTitle("Number of pixel holes");
  holes->SetStats(0);
  holes->Draw("COLZ");
  
  //leg->AddEntry(holes, legend[itr - tracks.begin()].c_str(), "lp");
 

  //leg->Draw();

  c->Print(("nSiHits_eta_" + layout + ".png").c_str(), "png");


       
  delete leg;
  delete c;


}
