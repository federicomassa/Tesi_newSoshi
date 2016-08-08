#ifndef ParticleAnalysis_HiggsAnalysis_H
#define ParticleAnalysis_HiggsAnalysis_H

// c++ include(s):
#include <iterator>
#include <vector>
#include <iostream>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include <EventLoop/Algorithm.h>
#include "AthLinks/ElementLink.h"

// package include(s):
#include <HistManager/TrackHistManager.h>
#include <HistManager/TruthHistManager.h>
#include <HistManager/EventHists.h>
#include <HistManager/VertexHists.h>
#include <HistManager/ClusterHistManager.h>
#include <HistManager/RunHistManager.h>

#include <Utility/assert.h>

#ifndef __MAKECINT__
#include "xAODTruth/TruthParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackStateValidationContainer.h"
#include "xAODTracking/TrackMeasurementValidationContainer.h"
#endif // not __MAKECINT__


class HiggsAnalysis : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  xAOD::TEvent *m_event;  //!
  int m_eventCounter;     //!

  TrackHists*       trkHist_all; //!
  TruthHistManager* truHist_all; //!
  TruthHistManager* truHist_hard; //!
  TruthHistManager* truHist_pileup; //!
  TruthHistManager* truHist_photons; //!
  TruthHistManager* truHist_visibleHard; //!
  TrackHistManager* trkHist_matched; //!
  TrackHistManager* trkHist_fakeCharge; //!
  TrackHistManager* trkHist_reco; //!
  ClusterHists*     clusHist_all; //!
  VertexHists*      vtxHist_secondary; //!
  EventHists*       eventHist_all; //!
  RunHistManager*   runHist_reco; //!
// Everything in the header file  that refers to the xAOD edm needs 
#ifndef __MAKECINT__
	
#endif // not __MAKECINT__


  // this is a standard constructor
  HiggsAnalysis ();

  void SetTargetParticleID(int idTarget) { m_idTarget = idTarget; }
  void SetTruthPtCut(std::vector<double> truthPtCut) {
    if (truthPtCut.size() != m_nTarget) {
      std::cout << "ERROR\t in SetTruthPtCut: wrong vector size" << std::endl;
      exit(1);
    }

    m_truthPtCut = truthPtCut;

  }
  void SetPtCut(double ptCut) { m_ptCut = ptCut;}
  void SetHitCut(int hitCut) {m_hitCut = hitCut;}
  void SetParentMass(double parentMass) {m_parentMass = parentMass;}
  void SetIntermediateMass(double intermediateMass) {m_intermediateMass = intermediateMass;}
  void SetDecayMass(double decayMass) {m_decayMass = decayMass;}
  void SetParentMassTolerance(double tolerance) {m_parentMassTolerance = tolerance;}
  void SetIntermediateMassTolerance(double intermediateMassTolerance) {m_intermediateMassTolerance = intermediateMassTolerance;}
  void SetNTarget(int nTarget) {m_nTarget = nTarget;}
  void SetDRCut(double dRCut) {m_dRCut = dRCut;}
  void SetEtaMax(double etaMax) {m_etaMax = etaMax;}
  void SetOnShellMassCut(const double& low, const double& high) {
    Assert("In ParticleAnalysis::SetOnShellMassCut(...): low > high", high - low > 1E-6);

    m_onShellMassCutLow = low;
    m_onShellMassCutHigh = high;

  }

  void SetOffShellMassCut(const double& low, const double& high) {
    Assert("In ParticleAnalysis::SetOffShellMassCut(...): low > high", high - low > 1E-6);

    m_offShellMassCutLow = low;
    m_offShellMassCutHigh = high;

  }

  void SetMassCut(const double& low, const double& high) {
    Assert("In ParticleAnalysis::SetMassCut(...): low > high", high - low > 1E-6);

    m_massCutLow = low;
    m_massCutHigh = high;

  }



  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

private:
  unsigned int m_idTarget, m_nTarget;
  std::vector<double> m_truthPtCut;
  double m_ptCut, m_parentMass,m_intermediateMass, m_decayMass, m_parentMassTolerance,m_intermediateMassTolerance, m_dRCut, m_etaMax;

  double m_massCutLow, m_massCutHigh;
  double m_onShellMassCutLow, m_onShellMassCutHigh;
  double m_offShellMassCutLow, m_offShellMassCutHigh;

  int m_hitCut;
  //Reconstruction cuts
  double m_trkHist_reco_ptCut;
  int m_trkHist_reco_hitsCut;

  // this is needed to distribute the algorithm to the workers
  ClassDef(HiggsAnalysis, 1);
};

#endif
