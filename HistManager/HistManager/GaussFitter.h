#ifndef GAUS_FITTER_H
#define GAUS_FITTER_H

#include <TFitResultPtr.h>

class TF1;
class TH1F;
class TCanvas;

class GaussFitter {

  TF1* m_fitFunc;
  TH1F* m_fitHist;
  
  unsigned int m_iterationNumber;
  static const unsigned int m_maxIterationNumber;

  double m_mean;
  double m_amplitude;
  double m_sigma;

  TFitResultPtr m_fitResult;
  int m_fitStatus;

  double m_threshold;
  static const double m_minThreshold;
  static const double m_defaultThreshold;

 public:
  GaussFitter(TH1F*);
  TF1* Fit();
  void SetThreshold(double);
  void DrawFitFunction(TCanvas*);
  void DrawHistogram(TCanvas*);
  void DrawHistWithFit(TCanvas*);


};

#endif //GAUS_FITTER_H
