#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

void CutEfficiency(const string& fnameExt, const string& fnameIExt = "", const string& fnameIncl = "", const string& latex = "") {

  bool writeLatex = latex.size() > 0 && fnameExt.size() > 0
    && fnameIExt.size() > 0 && fnameIncl.size() > 0;

  ofstream latexOut, latexOut27, latexOut32, latexOut4;

  if (writeLatex) {
    latexOut.open((latex + ".txt").c_str());
    latexOut27.open((latex+"27.txt").c_str());
    latexOut32.open((latex+"32.txt").c_str());
    latexOut4.open((latex+"4.txt").c_str());
  }

  std::vector<TFile*> f;
  f.push_back(new TFile(fnameExt.c_str()));
  if (fnameIExt.size() > 0)
    f.push_back(new TFile(fnameIExt.c_str()));
  if (fnameIncl.size() > 0)
    f.push_back(new TFile(fnameIncl.c_str()));

  std::vector<TH1F*> recoCutCode;
  std::vector<TH1F*> recoCutCode27;
  std::vector<TH1F*> recoCutCode32;
  std::vector<TH1F*> recoCutCode4;

  for (unsigned int i = 0; i < f.size(); i++) {
    recoCutCode.push_back((TH1F*) f[i]->Get("EventHist_all__recoCutCode"));
    recoCutCode27.push_back((TH1F*) f[i]->Get("EventHist_all__recoCutCode27"));
    recoCutCode32.push_back((TH1F*) f[i]->Get("EventHist_all__recoCutCode32"));
    recoCutCode4.push_back((TH1F*) f[i]->Get("EventHist_all__recoCutCode4"));
     }

  std::vector<std::string> cutName;
  cutName.push_back("\\textbf{Number of tracks}");
  cutName.push_back("\\textbf{Hits}            ");
  cutName.push_back("\\boldmath$p^T_4$         ");
  cutName.push_back("\\boldmath$\\Delta z$      ");
  cutName.push_back("\\textbf{Isolation}       ");
  cutName.push_back("\\textbf{Charge}          ");
  cutName.push_back("\\boldmath$p^T_1$         ");
  cutName.push_back("\\boldmath$p^T_2$         ");
  cutName.push_back("\\boldmath$p^T_3$         ");
  cutName.push_back("\\boldmath$p^T_4$         ");
  cutName.push_back("\\textbf{DR}              ");
  cutName.push_back("\\textbf{Off-shell muons} ");
  cutName.push_back("\\textbf{On-shell mass}   ");
  cutName.push_back("\\textbf{Off-shell mass}  ");

  if (true) {

    ///////// NO ETA SELECTION ////////

    std::vector<double> entriesOld,entries;
    // for (unsigned int i = 0; i < f.size(); i++) {
    //   entriesOld.push_back(recoCutCode[i]->GetEntries() - recoCutCode[i]->GetBinContent(1));   
    //   entries.push_back(entriesOld[i] - recoCutCode[i]->GetBinContent(3));
    // }

    for (unsigned int i = 0; i < f.size(); i++) {
      entriesOld.push_back(recoCutCode[i]->GetEntries());   
      entries.push_back(entriesOld[i] - recoCutCode[i]->GetBinContent(1));
    }


    std::cout <<    "Initial:          " << entriesOld[0] << std::endl;
    
    for (std::vector<std::string>::iterator itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries[0] << '\t' << entries[0]/entriesOld[0] << std::endl;
      
      if (writeLatex)
	latexOut << std::fixed << std::setprecision(3) << (*itr) << " & " << entries[0]/entriesOld[0] << " & " << entries[1]/entriesOld[1] << " & " << entries[2]/entriesOld[2] << " \\\\ \\hline" << std::endl;
      	       
					       
      entriesOld = entries;

      for (unsigned int i = 0; i < f.size(); i++)
	entries[i] -= recoCutCode[i]->GetBinContent(3 + itr - cutName.begin());
      
      
      
    }
    
    std::cout << std::endl;
    
    
  }

  if (writeLatex) 
    latexOut << "\\textbf{Total}" << " & " << recoCutCode[0]->GetBinContent(2)/(recoCutCode[0]->GetEntries()) << " & " << 
      recoCutCode[1]->GetBinContent(2)/(recoCutCode[1]->GetEntries()) << " & " <<
      recoCutCode[2]->GetBinContent(2)/(recoCutCode[2]->GetEntries()) << " \\\\ \\hline" << std::endl;
  
  



  if (true) {

    ///////// ETA 27 ////////

    std::vector<double> entriesOld,entries;
    // for (unsigned int i = 0; i < f.size(); i++) {
    //   entriesOld.push_back(recoCutCode27[i]->GetEntries() - recoCutCode27[i]->GetBinContent(1));   
    //   entries.push_back(entriesOld[i] - recoCutCode27[i]->GetBinContent(3));
    // }

    for (unsigned int i = 0; i < f.size(); i++) {
      entriesOld.push_back(recoCutCode27[i]->GetEntries());   
      entries.push_back(entriesOld[i] - recoCutCode27[i]->GetBinContent(1));
    }


    std::cout <<    "Initial:          " << entriesOld[0] << std::endl;
    
    for (std::vector<std::string>::iterator itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries[0] << '\t' << entries[0]/entriesOld[0] << std::endl;
      
      if (writeLatex)
	latexOut27 << std::fixed << std::setprecision(3) << (*itr) << " & " << entries[0]/entriesOld[0] << " & " << entries[1]/entriesOld[1] << " & " << entries[2]/entriesOld[2] << " \\\\ \\hline" << std::endl;
      	       
					       
      entriesOld = entries;

      for (unsigned int i = 0; i < f.size(); i++)
	entries[i] -= recoCutCode27[i]->GetBinContent(3 + itr - cutName.begin());
      
      
      
    }
    
    std::cout << std::endl;
    //////////////////////////////////
    
  }

  if (writeLatex) 
    latexOut27 << "\\textbf{Total}" << " & " << recoCutCode27[0]->GetBinContent(2)/(recoCutCode27[0]->GetEntries()) << " & " << 
      recoCutCode27[1]->GetBinContent(2)/(recoCutCode27[1]->GetEntries()) << " & " <<
      recoCutCode27[2]->GetBinContent(2)/(recoCutCode27[2]->GetEntries()) << " \\\\ \\hline" << std::endl;
  

  // ==============

  if (true) {

    ///////// ETA 32 ////////

    std::vector<double> entriesOld,entries;
    // for (unsigned int i = 0; i < f.size(); i++) {
    //   entriesOld.push_back(recoCutCode32[i]->GetEntries() - recoCutCode32[i]->GetBinContent(1));   
    //   entries.push_back(entriesOld[i] - recoCutCode32[i]->GetBinContent(3));
    // }

    for (unsigned int i = 0; i < f.size(); i++) {
      entriesOld.push_back(recoCutCode32[i]->GetEntries());   
      entries.push_back(entriesOld[i] - recoCutCode32[i]->GetBinContent(1));
    }


    std::cout <<    "Initial:          " << entriesOld[0] << std::endl;
    
    for (std::vector<std::string>::iterator itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries[0] << '\t' << entries[0]/entriesOld[0] << std::endl;
      
      if (writeLatex)
	latexOut32 << std::fixed << std::setprecision(3) << (*itr) << " & " << entries[0]/entriesOld[0] << " & " << entries[1]/entriesOld[1] << " & " << entries[2]/entriesOld[2] << " \\\\ \\hline" << std::endl;
      	       
					       
      entriesOld = entries;

      for (unsigned int i = 0; i < f.size(); i++)
	entries[i] -= recoCutCode32[i]->GetBinContent(3 + itr - cutName.begin());
      
      
      
    }
    
    std::cout << std::endl;
    //////////////////////////////////
    
  }

  if (writeLatex) 
    latexOut32 << "\\textbf{Total}" << " & " << recoCutCode32[0]->GetBinContent(2)/(recoCutCode32[0]->GetEntries()) << " & " << 
      recoCutCode32[1]->GetBinContent(2)/(recoCutCode32[1]->GetEntries()) << " & " <<
      recoCutCode32[2]->GetBinContent(2)/(recoCutCode32[2]->GetEntries()) << " \\\\ \\hline" << std::endl;
  

  // ===============

  if (true) {

    ///////// ETA 4 ////////

    std::vector<double> entriesOld,entries;
    // for (unsigned int i = 0; i < f.size(); i++) {
    //   entriesOld.push_back(recoCutCode4[i]->GetEntries() - recoCutCode4[i]->GetBinContent(1));   
    //   entries.push_back(entriesOld[i] - recoCutCode4[i]->GetBinContent(3));
    // }

    for (unsigned int i = 0; i < f.size(); i++) {
      entriesOld.push_back(recoCutCode4[i]->GetEntries());   
      entries.push_back(entriesOld[i] - recoCutCode4[i]->GetBinContent(1));
    }


    std::cout <<    "Initial:          " << entriesOld[0] << std::endl;
    
    for (std::vector<std::string>::iterator itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries[0] << '\t' << entries[0]/entriesOld[0] << std::endl;
      
      if (writeLatex)
	latexOut4 <<  std::fixed << std::setprecision(3) << (*itr) << " & " << entries[0]/entriesOld[0] << " & " << entries[1]/entriesOld[1] << " & " << entries[2]/entriesOld[2] << " \\\\ \\hline" << std::endl;
      	       
					       
      entriesOld = entries;

      for (unsigned int i = 0; i < f.size(); i++)
	entries[i] -= recoCutCode4[i]->GetBinContent(3 + itr - cutName.begin());
      
      
      
    }
    
    std::cout << std::endl;
    //////////////////////////////////
    
  }

  if (writeLatex) 
    latexOut4 << "\\textbf{Total}" << " & " << recoCutCode4[0]->GetBinContent(2)/(recoCutCode4[0]->GetEntries()) << " & " << 
      recoCutCode4[1]->GetBinContent(2)/(recoCutCode4[1]->GetEntries()) << " & " <<
      recoCutCode4[2]->GetBinContent(2)/(recoCutCode4[2]->GetEntries()) << " \\\\ \\hline" << std::endl;
  
   




}

