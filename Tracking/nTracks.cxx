#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "../plotStyle.h"

#include <string>
#include <vector>

void nTracks() {

  std::vector<TFile*> file;

  file.push_back(new TFile("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_ExtBrl4/hist-pi15pu0.root"));
  file.push_back(new TFile("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_ExtBrl4/hist-pi15pu50.root"));
  file.push_back(new TFile("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_ExtBrl4/hist-pi15pu140.root"));
  file.push_back(new TFile("/afs/cern.ch/user/f/fmassa/newSoshi/datasetPi_15GeV_ExtBrl4/hist-pi15pu200.root"));

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
  
  std::vector<TH1F*> tracks;
  std::vector<TH1F*> primary;
  std::vector<TH1F*> secondary;
  TH1F* secondary1GeV = (TH1F*) file[0]->Get("EventHist_all__nSecondaryChargedTruth1GeV");

  TH1F* pixel_matched = (TH1F*) file[3]->Get("TrackHist_reco_all__nPixHits");
  TH1F* pixel_unmatched = (TH1F*) file[3]->Get("TrackHist_pileup_all__nPixHits");
  TH1F* SCT_matched = (TH1F*) file[3]->Get("TrackHist_reco_all__nSCTHits");
  TH1F* SCT_unmatched = (TH1F*) file[3]->Get("TrackHist_pileup_all__nSCTHits");
  TH1F* total_matched = (TH1F*) file[3]->Get("TrackHist_reco_all__nSiHits");
  TH1F* total_unmatched = (TH1F*) file[3]->Get("TrackHist_pileup_all__nSiHits");
  TH1F* chisquare_matched = (TH1F*) file[3]->Get("TrackHist_reco_all__chiSqPerDof");
  TH1F* chisquare_unmatched = (TH1F*) file[3]->Get("TrackHist_pileup_all__chiSqPerDof");


  for (std::vector<TFile*>::iterator itr = file.begin();
       itr != file.end();
       itr++) {
    tracks.push_back((TH1F*) (*itr)->Get("EventHist_all__nRecotrk"));
    primary.push_back((TH1F*) (*itr)->Get("EventHist_all__nPrimaryChargedTruth")); 
    secondary.push_back((TH1F*) (*itr)->Get("EventHist_all__nSecondaryChargedTruth")); 
  }

  TCanvas *c = new TCanvas();
  
  TLegend* leg = new TLegend(0.60,0.85,0.88,0.65);

  for (std::vector<TH1F*>::iterator itr = tracks.begin();
       itr != tracks.end();
       itr++) {
    
    plotStyle((*itr), "", 0, "Number of reconstructed tracks in event","A.U.", color[itr - tracks.begin()], 0, 7);
    (*itr)->Scale(1.0/(*itr)->GetEntries());
    
    if (itr == tracks.begin())
      (*itr)->Draw();
    else
      (*itr)->Draw("same");

    leg->AddEntry((*itr), legend[itr - tracks.begin()].c_str(), "lp");

  }

  leg->Draw();

  c->Print("nRecoTracks.png", "png");

       
  delete leg;
  delete c;


  /////////////////////////////////

  c = new TCanvas();
  leg = new TLegend(0.60,0.85,0.88,0.65);

  for (std::vector<TH1F*>::iterator itr = primary.begin();
       itr != primary.end();
       itr++) {
    
    plotStyle((*itr), "", 0, "Number of primary charged particles in event","A.U.", color[itr - primary.begin()], 0, 18);
    (*itr)->Scale(1.0/(*itr)->GetEntries());

    if (itr == primary.begin())
      (*itr)->Draw();
    else
      (*itr)->Draw("same");

    leg->AddEntry((*itr), legend[itr - primary.begin()].c_str(), "lp");

  }

  leg->Draw();

  c->Print("nPrimaryChargedTruth.png", "png");

       
  delete leg;
  delete c;
  
  //////////////////////

  c = new TCanvas();
  leg = new TLegend(0.60,0.85,0.88,0.65);

  for (std::vector<TH1F*>::iterator itr = secondary.begin();
       itr != secondary.end();
       itr++) {
    
    plotStyle((*itr), "", 0, "Number of secondary charged particles in event","A.U.", color[itr - secondary.begin()], 0, 100);
    (*itr)->Scale(1.0/(*itr)->GetEntries());

    if (itr == secondary.begin())
      (*itr)->Draw();
    else
      (*itr)->Draw("same");

    (c->cd())->SetLogy();

    leg->AddEntry((*itr), legend[itr - secondary.begin()].c_str(), "lp");

  }

  leg->Draw();

  c->Print("nSecondaryChargedTruth.png", "png");

       
  delete leg;
  delete c;

  /////////////////////////////////

  c = new TCanvas();

  plotStyle(secondary1GeV, "", 0, "Number of secondary charged particles with p_{T} > 1 GeV","A.U.", kBlue, 0, 12);
  secondary1GeV->Scale(1.0/secondary1GeV->GetEntries());

  secondary1GeV->Draw();
  
  (c->cd())->SetLogy();
  

  c->Print("nSecondaryChargedTruth1GeV.png", "png");
  
  delete c;
  
  //////////////////// 

  c = new TCanvas();
  leg = new TLegend(0.60,0.85,0.88,0.65);
    
  plotStyle(pixel_matched, "", 0, "Number of pixel hits per track","A.U.", kBlue, 0, 25);
  plotStyle(pixel_unmatched, "", 0, "Number of pixel hits per track","A.U.", kRed, 0, 25);
  
  pixel_matched->Scale(1.0/pixel_matched->GetEntries());
  pixel_unmatched->Scale(1.0/pixel_unmatched->GetEntries());

  if (pixel_matched->GetMaximum() > pixel_unmatched->GetMaximum()) {
    pixel_matched->Draw();
    pixel_unmatched->Draw("same");
  }
  else {
    pixel_unmatched->Draw();
    pixel_matched->Draw("same");
  }

  
  leg->AddEntry(pixel_matched, "Matched tracks", "lp");
  leg->AddEntry(pixel_unmatched, "Unmatched tracks", "lp");
 

  leg->Draw();
  
  c->Print("nPixHits.png", "png");
  
  
  delete leg;
  delete c;

  //////////////////// 

  c = new TCanvas();
  leg = new TLegend(0.60,0.85,0.88,0.65);
    
  plotStyle(SCT_matched, "", 0, "Number of SCT hits per track","A.U.", kBlue, 0, 18);
  plotStyle(SCT_unmatched, "", 0, "Number of SCT hits per track","A.U.", kRed, 0, 18);
  
  SCT_matched->Scale(1.0/SCT_matched->GetEntries());
  SCT_unmatched->Scale(1.0/SCT_unmatched->GetEntries());


  if (SCT_matched->GetMaximum() > SCT_unmatched->GetMaximum()) {
    SCT_matched->Draw();
    SCT_unmatched->Draw("same");
  }
  else {
    SCT_unmatched->Draw();
    SCT_matched->Draw("same");
  }
  
  leg->AddEntry(SCT_matched, "Matched tracks", "lp");
  leg->AddEntry(SCT_unmatched, "Unmatched tracks", "lp");
 

  leg->Draw();
  
  c->Print("nSCTHits.png", "png");
  
  
  delete leg;
  delete c;

  //////////////////// 

  c = new TCanvas();
  leg = new TLegend(0.60,0.85,0.88,0.65);
    
  plotStyle(total_matched, "", 0, "Number of total hits per track","A.U.", kBlue, 0, 30);
  plotStyle(total_unmatched, "", 0, "Number of total hits per track","A.U.", kRed, 0, 30);
  
  total_matched->Scale(1.0/total_matched->GetEntries());
  total_unmatched->Scale(1.0/total_unmatched->GetEntries());

  if (total_matched->GetMaximum() > total_unmatched->GetMaximum()) {
    total_matched->Draw();
    total_unmatched->Draw("same");
  }
  else {
    total_unmatched->Draw();
    total_matched->Draw("same");
  }
  
  leg->AddEntry(total_matched, "Matched tracks", "lp");
  leg->AddEntry(total_unmatched, "Unmatched tracks", "lp");
 

  leg->Draw();
  
  c->Print("nSiHits.png", "png");
  
  
  delete leg;
  delete c;

  ///////////////

  c = new TCanvas();
  leg = new TLegend(0.60,0.85,0.88,0.65);
    
  plotStyle(chisquare_matched, "", 0, "Track #chi^{2}/DOF","A.U.", kBlue, 0, 4);
  plotStyle(chisquare_unmatched, "", 0, "Track #chi^{2}/DOF","A.U.", kRed, 0, 4);
  
  chisquare_matched->Scale(1.0/chisquare_matched->GetEntries());
  chisquare_unmatched->Scale(1.0/chisquare_unmatched->GetEntries());

  if (chisquare_matched->GetMaximum() > chisquare_unmatched->GetMaximum()) {
    chisquare_matched->Draw();
    chisquare_unmatched->Draw("same");
  }
  else {
    chisquare_unmatched->Draw();
    chisquare_matched->Draw("same");
  }

  
  leg->AddEntry(chisquare_matched, "Matched tracks", "lp");
  leg->AddEntry(chisquare_unmatched, "Unmatched tracks", "lp");
 

  leg->Draw();
  
  c->Print("chiSquarePerDof.png", "png");
  
  
  delete leg;
  delete c;
  

}
