#include "HistManager/GaussFitter.h"
#include "HistManager/PrintMessage.h"
#include "HistManager/assert.h"
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>

const unsigned int GaussFitter::m_maxIterationNumber = 1000;
const double GaussFitter::m_minThreshold = 1e-6;
const double GaussFitter::m_defaultThreshold = 0.0001;

GaussFitter::GaussFitter(TH1F* h, const Verbosity& verbosity_level = Verbosity::ERROR) :
  m_fitHist(h),
  m_verbosityLevel(verbosity_level)
{
  Assert("Histogram pointer is null", h != NULL);
  Assert("Histogram range not properly set", h->GetXaxis()->GetXmin() < h->GetXaxis()->GetXmax());
  m_iterationNumber = 0;
  m_fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  m_fitFunc->SetParName(0, "Amplitude");
  m_fitFunc->SetParName(1, "Mean");
  m_fitFunc->SetParName(2, "Sigma");
  m_threshold = m_defaultThreshold;
  m_fitStatus = -1;
  m_mean = -1E9;
  m_amplitude = -1E9;
  m_sigma = -1E9;
}

TF1* GaussFitter::Fit() {

  std::stringstream ss;
  std::string iterationNumber;

  ss << ++m_iterationNumber;
  ss >> iterationNumber;
  PrintMessage("GaussFitter::Fit()\tFitting histogram... iteration " + iterationNumber, Verbosity::DEBUG, m_verbosityLevel);

  if (m_fitHist->GetEntries() == 0) 
    {
      PrintMessage("GaussFitter::Fit()\tHistogram is empty. Returning null pointer", Verbosity::WARNING, m_verbosityLevel);
      return NULL;
	}

  if (m_iterationNumber == 1) {
    //Initial parameters are taken directly from histogram
    m_mean = m_fitHist->GetMean();
    m_amplitude = m_fitHist->GetBinContent(m_fitHist->FindBin(m_mean));
    m_sigma = m_fitHist->GetRMS();
  }
  else {
    if (m_fitFunc->GetParameter("Amplitude") != TMath::QuietNaN() &&
	m_fitFunc->GetParameter("Mean") != TMath::QuietNaN() &&
	m_fitFunc->GetParameter("Sigma") != TMath::QuietNaN()) {
      
      m_amplitude = m_fitFunc->GetParameter("Amplitude");
      m_mean = m_fitFunc->GetParameter("Mean");
      m_sigma = m_fitFunc->GetParameter("Sigma");
    }
    else {
      PrintMessage("GaussFitter::Fit()\tWARNING\tFitting did not end successfully", Verbosity::WARNING, m_verbosityLevel);
    }
  }
  
  m_fitFunc -> SetRange(m_mean - 4*m_sigma, m_mean + 4*m_sigma);
  m_fitFunc -> SetParameter(0, m_amplitude);
  m_fitFunc -> SetParameter(1, m_mean);
  m_fitFunc -> SetParameter(2, m_sigma);

  m_fitResult = m_fitHist -> Fit(m_fitFunc, "RSNQ");
  m_fitStatus = m_fitResult;

  if (m_fitStatus != 0) {
    PrintMessage("GaussFitter::Fit()\tWARNING\tFit not successful. Returning null", Verbosity::WARNING, m_verbosityLevel);
    return NULL;
  }

  if (TMath::Abs((m_fitFunc->GetParameter("Amplitude") - m_amplitude)/m_amplitude) < m_threshold &&
      TMath::Abs((m_fitFunc->GetParameter("Mean") - m_mean)/m_mean) < m_threshold &&
      TMath::Abs((m_fitFunc->GetParameter("Sigma") - m_sigma)/m_sigma) < m_threshold) {
    return m_fitFunc;
  }
  else {
    if (m_iterationNumber < m_maxIterationNumber)
      return Fit();
    else {
      PrintMessage("GaussFitter::Fit()\tReached maximum number of iterations", Verbosity::INFO, m_verbosityLevel);
      return m_fitFunc;
    }
  }

  return m_fitFunc;

}

void GaussFitter::SetThreshold(double threshold = m_defaultThreshold) {
  if (threshold > m_minThreshold)
    m_threshold = threshold;
  else
    m_threshold = m_defaultThreshold;
}

void GaussFitter::DrawFitFunction(TCanvas* c) const {
  if (m_fitFunc != NULL) {
    if (m_fitFunc->GetParameter("Amplitude") != TMath::QuietNaN() &&
	m_fitFunc->GetParameter("Mean") != TMath::QuietNaN() &&
	m_fitFunc->GetParameter("Sigma") != TMath::QuietNaN()) {
      
      c->cd();
      m_fitFunc->DrawClone();
    }
    
    else {
      PrintMessage("DrawFitFunction():\t Nothing to draw", Verbosity::ERROR, m_verbosityLevel);
    }
  }
  else {
    PrintMessage("GaussFitter::DrawFitFunction()\tWARNING\tFitting function is null", Verbosity::ERROR, m_verbosityLevel);
  }
}

void GaussFitter::DrawHistogram(TCanvas* c) const {
  if (m_fitHist != NULL) {
    c->cd();
    m_fitHist->DrawCopy("same");
  }
  else {
    PrintMessage("GaussFitter::DrawHistogram()\tERROR\tHistogram is null", Verbosity::ERROR, m_verbosityLevel);
  }
}

void GaussFitter::DrawHistWithFit(TCanvas* c) const {
  if (m_fitHist != NULL) {
    if (m_fitFunc != NULL) {
      c->cd();
      m_fitHist->DrawCopy();
      m_fitFunc->Draw("same");
    }
    else {
      PrintMessage("GaussFitter::DrawHistWithFit()\tERROR\tFitting function is null", Verbosity::ERROR, m_verbosityLevel);
    }
  }
  else {
    PrintMessage("GaussFitter::DrawHistWithFit()\tERROR\tHistogram is null", Verbosity::ERROR, m_verbosityLevel);
  }
}

Double_t GaussFitter::GetAmplitude() const {
  if (m_fitFunc != NULL) {
    if (m_fitFunc->GetParameter("Amplitude") != TMath::QuietNaN() && m_fitStatus == 0) 
      return m_fitFunc->GetParameter("Amplitude"); 
    else 
      return 0.0;
      }
 
  
  else {
    PrintMessage("GaussFitter::GetAmplitude()\tWARNING\tFitting function is null", Verbosity::WARNING, m_verbosityLevel);
    return -1E9;
  }
}

Double_t GaussFitter::GetMean() const {
  if (m_fitFunc != NULL) {
    if (m_fitFunc->GetParameter("Mean") != TMath::QuietNaN() && m_fitStatus == 0) 
      return m_fitFunc->GetParameter("Mean"); 
    else 
      return 0.0;
      }
 
  
  else {
    PrintMessage("GaussFitter::GetMean()\tWARNING\tFitting function is null", Verbosity::WARNING, m_verbosityLevel);
    return -1E9;
  }
}


Double_t GaussFitter::GetSigma() const {
  if (m_fitFunc != NULL) {
    if (m_fitFunc->GetParameter("Sigma") != TMath::QuietNaN() && m_fitStatus == 0) 
      return TMath::Abs(m_fitFunc->GetParameter("Sigma")); 
    else 
      return 0.0;
      }
 
  
  else {
    PrintMessage("GaussFitter::GetSigma()\tWARNING\tFitting function is null", Verbosity::WARNING, m_verbosityLevel);
    return -1E9;
  }
}

Double_t GaussFitter::GetMeanError() const {
  if (m_fitFunc != NULL) {
    if (m_fitFunc->GetParameter("Mean") != TMath::QuietNaN() && m_fitStatus == 0) 
      return TMath::Sqrt((*m_fitResult).GetCovarianceMatrix()(1,1)); 
    else 
      return 0.0;
      }
 
  
  else {
    PrintMessage("GaussFitter::GetMeanError()\tWARNING\tFitting function is null", Verbosity::WARNING, m_verbosityLevel);
    return -1E9;
  }
}

Double_t GaussFitter::GetAmplitudeError() const {
  if (m_fitFunc != NULL) {
    if (m_fitFunc->GetParameter("Amplitude") != TMath::QuietNaN() && m_fitStatus == 0) 
      return TMath::Sqrt((*m_fitResult).GetCovarianceMatrix()(0,0)); 
    else 
      return 0.0;
      }
 
  
  else {
    PrintMessage("GaussFitter::GetAmplitudeError()\tWARNING\tFitting function is null", Verbosity::WARNING, m_verbosityLevel);
    return -1E9;
  }
}

Double_t GaussFitter::GetSigmaError() const {
  if (m_fitFunc != NULL) {
    if (m_fitFunc->GetParameter("Sigma") != TMath::QuietNaN() && m_fitStatus == 0) 
      return TMath::Sqrt((*m_fitResult).GetCovarianceMatrix()(2,2)); 
    else 
      return 0.0;
      }
 
  
  else {
    PrintMessage("GaussFitter::GetSigmaError()\tWARNING\tFitting function is null", Verbosity::WARNING, m_verbosityLevel);
    return -1E9;
  }
}


