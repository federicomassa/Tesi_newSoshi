#ifndef GAUS_FITTER_H
#define GAUS_FITTER_H

#include <TFitResultPtr.h>
#include "HistManager/PrintMessage.h"

class TF1;
class TH1F;
class TCanvas;

class GaussFitter {

  TF1* m_fitFunc;
  TH1F* m_fitHist;
  
  unsigned int m_iterationNumber;
  static const unsigned int m_maxIterationNumber;

  const Verbosity& m_verbosityLevel;
  
  double m_mean;
  double m_amplitude;
  double m_sigma;

  TFitResultPtr m_fitResult;
  int m_fitStatus;

  double m_threshold;
  static const double m_minThreshold;
  static const double m_defaultThreshold;

 public:
  GaussFitter(TH1F*, const Verbosity&);
  TF1* Fit();
  void SetThreshold(double);
  void DrawFitFunction(TCanvas*) const;
  void DrawHistogram(TCanvas*) const;
  void DrawHistWithFit(TCanvas*) const;
  Double_t GetMean() const;
  Double_t GetAmplitude() const;
  Double_t GetSigma() const;
  Double_t GetMeanError() const;
  Double_t GetAmplitudeError() const;
  Double_t GetSigmaError() const;
  
};

#endif //GAUS_FITTER_H
