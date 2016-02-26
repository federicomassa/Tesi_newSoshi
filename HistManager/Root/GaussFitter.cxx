#include "HistManager/GaussFitter.h"
#include "HistManager/assert.h"
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>

const unsigned int GaussFitter::m_maxIterationNumber = 1000;
const double GaussFitter::m_minThreshold = 1e-6;
const double GaussFitter::m_defaultThreshold = 0.0001;

GaussFitter::GaussFitter(TH1F* h) :
  m_fitHist(h)
{
  Assert("Histogram pointer is null", h != NULL);
  Assert("Histogram range not properly set", h->GetXaxis()->GetXmin() < h->GetXaxis()->GetXmax());
  m_iterationNumber = 0;
  m_fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  m_fitFunc->SetParName(0, "Amplitude");
  m_fitFunc->SetParName(1, "Mean");
  m_fitFunc->SetParName(2, "Sigma");
  m_threshold = m_defaultThreshold;
}

TF1* GaussFitter::Fit() {
  std::cout << "GaussFitter::Fit()\tFitting histogram... iteration " << ++m_iterationNumber << std::endl;

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
      std::cout << "WARNING: Fitting did not end successfully" << std::endl;
    }
  }
  
  m_fitFunc -> SetRange(m_mean - 4*m_sigma, m_mean + 4*m_sigma);
  m_fitFunc -> SetParameter(0, m_amplitude);
  m_fitFunc -> SetParameter(1, m_mean);
  m_fitFunc -> SetParameter(2, m_sigma);

  m_fitResult = m_fitHist -> Fit(m_fitFunc, "RSN");
  m_fitStatus = m_fitResult;

  Assert("ERROR: fit ended with nonzero status", m_fitStatus == 0);
  if (TMath::Abs((m_fitFunc->GetParameter("Amplitude") - m_amplitude)/m_amplitude) < m_threshold &&
      TMath::Abs((m_fitFunc->GetParameter("Mean") - m_mean)/m_mean) < m_threshold &&
      TMath::Abs((m_fitFunc->GetParameter("Sigma") - m_sigma)/m_sigma) < m_threshold) {
    return m_fitFunc;
  }
  else {
    if (m_iterationNumber < m_maxIterationNumber)
      return Fit();
    else {
      std::cout << "GaussFitter::Fit()\tReached maximum number of iterations" << std::endl;
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

void GaussFitter::DrawFitFunction(TCanvas* c) {
  if (m_fitFunc->GetParameter("Amplitude") != TMath::QuietNaN() &&
      m_fitFunc->GetParameter("Mean") != TMath::QuietNaN() &&
      m_fitFunc->GetParameter("Sigma") != TMath::QuietNaN()) {

    c->cd();
    m_fitFunc->DrawClone();
  }

  else {
    std::cout << "In DrawFitFunction():\t Nothing to draw" << std::endl;
  }

}

void GaussFitter::DrawHistogram(TCanvas* c) {
  
  c->cd();
  m_fitHist->DrawCopy("same");

}

void GaussFitter::DrawHistWithFit(TCanvas* c) {

  c->cd();
  m_fitHist->DrawCopy();
  m_fitFunc->Draw("same");

}
