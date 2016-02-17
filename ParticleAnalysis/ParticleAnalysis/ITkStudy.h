#ifndef ParticleAnalysis_ITkStudy_H
#define ParticleAnalysis_ITkStudy_H

// c++ include(s):
#include <iterator>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include <EventLoop/Algorithm.h>
#include "AthLinks/ElementLink.h"

// package include(s):
#include <HistManager/TrackHistManager.h>
#include <HistManager/TruthHists.h>
#include <HistManager/EventHists.h>
#include <HistManager/VertexHists.h>
#include <HistManager/ClusterHistManager.h>


#ifndef __MAKECINT__
#include "xAODTruth/TruthParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackStateValidationContainer.h"
#include "xAODTracking/TrackMeasurementValidationContainer.h"
#endif // not __MAKECINT__


class ITkStudy : public EL::Algorithm
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
  TruthHists*       truHist_stable; //!
  TrackHistManager* trkHist_reco; //!
  ClusterHists*     clusHist_all; //!
  VertexHists*      vtxHist_secondary; //!
  EventHists*       eventHist_all; //!

// Everything in the header file  that refers to the xAOD edm needs 
#ifndef __MAKECINT__
	
#endif // not __MAKECINT__


  // this is a standard constructor
  ITkStudy ();

  void SetTargetParticleID(int idTarget) { m_idTarget = idTarget; }
  void SetTargetParticleE(int EnergyTarget) { m_EnergyTarget = EnergyTarget; }

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
  int m_idTarget,m_EnergyTarget;

  // this is needed to distribute the algorithm to the workers
  ClassDef(ITkStudy, 1);
};

#endif
