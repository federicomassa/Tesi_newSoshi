/******************************************
 *
 * Code to analysize simple samples for 
 * study tracking efficiency
 *
 * G. Facini
 * Wed Sep 10 07:21:58 CEST 2014
 *
 ******************************************/

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ParticleAnalysis/ITkStudy.h>
#include <HistManager/RunHists.h>
#include <HistManager/EventFeaturesInterface.h>
#include <HistManager/TrackHelper.h>

// EDM includes: - if move to header file will not compile!
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

#include "xAODTracking/VertexContainer.h"

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include <map>
#include <limits>
#include <string>
#include <sstream>
#include <iostream>

// this is needed to distribute the algorithm to the workers
ClassImp(ITkStudy)

float deltaPhi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if (dphi >= TMath::Pi()) { dphi = 2.0*TMath::Pi()-dphi; }
  if (dphi < -TMath::Pi()) { dphi = 2.0*TMath::Pi()+dphi; }
  return dphi;
}

float deltaR(float phi1, float phi2, float eta1, float eta2) {
  float dphi = phi1 - phi2;
  if (dphi >= TMath::Pi()) { dphi = 2.0*TMath::Pi()-dphi; }
  if (dphi < -TMath::Pi()) { dphi = 2.0*TMath::Pi()+dphi; }
  float dR = TMath::Sqrt(dphi*dphi+(eta1-eta2)*(eta1-eta2));
  return dR;
}

ITkStudy::ITkStudy() {

}

EL::StatusCode ITkStudy::setupJob(EL::Job& job) {

  job.useXAOD ();
  xAOD::Init( "PartAnaxAOD" ).ignore(); // call before opening first file

  // count number of events
  m_eventCounter = 0;
  std::cout << "Set counter to 0 " << m_eventCounter << std::endl;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkStudy::histInitialize() {

  // In TrackHistManager constructor:
  // 1st bool is to fill plots for tracks matched to primary particles
  // 2nd bool is to fill plots for tracks matched to secondary particles
  // 3rd bool is to fill plots for fake tracks 
  // 4th bool is if to fill duplicates (need 1st and 2nd be true)
  // 5th bool is if to fill plots based on hit content on track (need 1st be true)
  
  // recoTracks histograms
  trkHist_all       = new TrackHists("TrackAll"); //every track, no cuts here
  truHist_stable    = new TruthHists("stable"); //Filled for every truth hs particle, no cuts
  trkHist_reco      = new TrackHistManager("reco",true,true,true,true,true); //every track that matches the hs
  clusHist_all      = new ClusterHists("SiliconAll");
  vtxHist_secondary = new VertexHists("Vertex");
  eventHist_all     = new EventHists("all");
  runHist_reco      = new RunHistManager(trkHist_reco);

  trkHist_all       -> Init(wk());
  truHist_stable    -> Init(wk());
  trkHist_reco      -> Init(wk());
  clusHist_all      -> Init(wk());
  vtxHist_secondary -> Init(wk());
  eventHist_all     -> Init(wk());
  runHist_reco      -> Init(wk());

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkStudy::fileExecute () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkStudy::changeInput(bool /* firstFile */) {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkStudy::initialize() {

  m_event = wk()->xaodEvent();

  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int
	m_eventCounter=0;
	
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkStudy::execute() {
  // print every 100 events, so we know where we are:
  if(m_eventCounter % 1 == 0) {
    Info("execute()", "Event number = %i", m_eventCounter );
  }

  if (m_eventCounter == 11090 || m_eventCounter == 40465) 
    {
      m_eventCounter++;
      return EL::StatusCode::SUCCESS;
    }
  //std::cout << "Event count " << m_eventCounter << std::endl;
  m_eventCounter++;
  if(m_eventCounter > 1e5) { 
    if(m_eventCounter == 1e5) { 
      std::cout << "STOPPING after event " << m_eventCounter << std::endl;
      std::cout << "Look in here to change -> ITkStudy :: execute()" << std::endl;
    }
    return EL::StatusCode::SUCCESS; 
  }

  EventFeatures eventFeatures;

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;

  //---------------------------
  // Reset barcode map in histogram class
  //---------------------------
  trkHist_reco->resetBarcodes();

  //---------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  xAOD::TStore store; // needed for shallow copies

  //---------------------
  // TRUTH PARTICLES
  //---------------------
  // get track container of interest
  const xAOD::TruthParticleContainer* truthPartsOrig = 0;
  if ( !m_event->retrieve( truthPartsOrig, "TruthParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Truth Particle container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //---------------------
  // RECONSTRUCTED TRACKS
  //---------------------
  // get track container of interest
  const xAOD::TrackParticleContainer* recoTracksOrig = 0;
  if ( !m_event->retrieve( recoTracksOrig, "InDetTrackParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Reconstructed Track container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  // make a shallow copy
  std::pair< xAOD::TrackParticleContainer*, xAOD::ShallowAuxContainer* > recoTracksShallowCopyPair = xAOD::shallowCopyContainer( *recoTracksOrig );
  xAOD::TrackParticleContainer* recoTracks = recoTracksShallowCopyPair.first;

  //---------------------
  // Pixel Clusters
  //---------------------
  // get cluster container of interest
  const xAOD::TrackMeasurementValidationContainer* pixClustersOrig = 0;
  if ( !m_event->retrieve( pixClustersOrig, "PixelClusters" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Pixel Cluster container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  // make a shallow copy
  std::pair< xAOD::TrackMeasurementValidationContainer*, xAOD::ShallowAuxContainer* > pixClustersPair = xAOD::shallowCopyContainer( *pixClustersOrig );
  xAOD::TrackMeasurementValidationContainer* pixClusters =  pixClustersPair.first;

  //---------------------
  // SCT Clusters
  //---------------------
  // get cluster container of interest
  const xAOD::TrackMeasurementValidationContainer* sctClustersOrig = 0;
  if ( !m_event->retrieve( sctClustersOrig, "SCT_Clusters" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve SCT Cluster container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  // make a shallow copy
  std::pair< xAOD::TrackMeasurementValidationContainer*, xAOD::ShallowAuxContainer* > sctClustersPair = xAOD::shallowCopyContainer( *sctClustersOrig );
  xAOD::TrackMeasurementValidationContainer* sctClusters =  sctClustersPair.first;

  //==================
  // Vertex container
  //==================
  const xAOD::VertexContainer* primaryVertex = 0;
  if (!m_event->retrieve(primaryVertex,"PrimaryVertices").isSuccess()) {
    Error("execute()", "Failed to retrieve PrimaryVertices container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //---------------------
  // FIND THE MOTHER
  //---------------------
  int numTruth = 0;
  int numChargedTruth = 0;
  int numNeutralTruth = 0;
  int numTruth1GeV = 0;
  int numChargedTruth1GeV = 0;
  int numNeutralTruth1GeV = 0;

  //Create the new container and its auxiliary store for only primary truth particles 
  //need this since link to mother and children is not yet in the AR...
  double pxQuark = 0.0;
  double pyQuark = 0.0;
  double pzQuark = 0.0;
  xAOD::TruthParticleContainer* truthPriParts = new xAOD::TruthParticleContainer();
  xAOD::TruthParticleAuxContainer* truthPriPartsAux = new xAOD::TruthParticleAuxContainer();
  truthPriParts->setStore( truthPriPartsAux ); //< Connect the two


  // ---------------------- ITERATION OVER TRUTH PARTICLES ---------------------------- //
  for (xAOD::TruthParticleContainer::const_iterator truthPart_itr=truthPartsOrig->begin(); truthPart_itr!=truthPartsOrig->end(); truthPart_itr++) {

    if (TMath::Abs((*truthPart_itr)->pdgId())<6 || (*truthPart_itr)->pdgId()==21) {
      pxQuark = (*truthPart_itr)->px();
      pyQuark = (*truthPart_itr)->py();
      pzQuark = (*truthPart_itr)->pz();
    }

    if ((*truthPart_itr)->status()!=1) { continue; }

    if (TMath::Abs((*truthPart_itr)->pdgId())!=13) {
      numTruth++;
      if (TMath::Abs((*truthPart_itr)->charge())==1) { numChargedTruth++; }
      if (TMath::Abs((*truthPart_itr)->charge())==0) { numNeutralTruth++; }

      if ((*truthPart_itr)->pt()>1000.0) {
        numTruth1GeV++;
        if (TMath::Abs((*truthPart_itr)->charge())==1) { numChargedTruth1GeV++; }
        if (TMath::Abs((*truthPart_itr)->charge())==0) { numNeutralTruth1GeV++; }
      }
    }
    if ((*truthPart_itr)->charge()==0) { continue; }

    // Copy this particle to the new container:
    xAOD::TruthParticle* tp = new xAOD::TruthParticle();
    tp->makePrivateStore( **truthPart_itr );
    truthPriParts->push_back( tp );
  }
  // ---------------------- END ITERATION OVER TRUTH PARTICLES ---------------------------- //


  //==========================
  // Count number of vertices
  //==========================
  int numVtx = 0;
  for (xAOD::VertexContainer::const_iterator vtx_itr=primaryVertex->begin(); vtx_itr!=primaryVertex->end(); vtx_itr++) {
    if ((*vtx_itr)->vertexType()!=0) { numVtx++; }
  }


  //================================================
  // Create new track container with truth matching
  //================================================
  xAOD::TrackParticleContainer* itkTrack = new xAOD::TrackParticleContainer();
  xAOD::TrackParticleAuxContainer* itkTrackAux = new xAOD::TrackParticleAuxContainer();
  itkTrack->setStore( itkTrackAux ); //< Connect the two

  int nIsoTrk = 0;
  int numRecoTracks = 0;
  int numRecoTracks1GeV = 0;
  for (xAOD::TrackParticleContainer::const_iterator recoTrk_itr=recoTracks->begin(); recoTrk_itr!=recoTracks->end(); recoTrk_itr++) {

    numRecoTracks++;
    if ((*recoTrk_itr)->pt()>1000.0) { numRecoTracks1GeV++; }

    xAOD::TruthParticleContainer::const_iterator genMatched_itr;
    xAOD::TruthParticleContainer::const_iterator genTargetMatched_itr;

    float mindR = 100.0;
    float matchedR = 100.0;
    for (xAOD::TruthParticleContainer::const_iterator genPrt_itr=truthPartsOrig->begin(); genPrt_itr!=truthPartsOrig->end(); genPrt_itr++) {
      if((*genPrt_itr)->barcode() <= 10e3 || (*genPrt_itr)->barcode() >= 200e3) { continue; }
      if ((*genPrt_itr)->status()!=1) { continue; }
      if ((*genPrt_itr)->pdgId()==22) { continue; }
      float dR = deltaR((*genPrt_itr)->phi(),(*recoTrk_itr)->phi(),(*genPrt_itr)->eta(),(*recoTrk_itr)->eta());
      if (dR<mindR) { genMatched_itr = genPrt_itr; mindR = dR; }
      if (TMath::Abs((*genPrt_itr)->pdgId())==m_idTarget && TMath::Abs((*genPrt_itr)->pt()-m_EnergyTarget)<0.1) {
        genTargetMatched_itr = genPrt_itr; matchedR = dR;
      }
    }

    if (matchedR<0.05) {
      nIsoTrk++;
    }

    xAOD::TrackParticle* tp = new xAOD::TrackParticle();
    tp->makePrivateStore( **recoTrk_itr);

    // set new variables
    static SG::AuxElement::Decorator<float> mindRmatching("mindRmatching");
    mindRmatching(*tp) = mindR;

    static SG::AuxElement::Decorator<float> matchedDR("matchedDR");
    matchedDR(*tp) = matchedR;

    static SG::AuxElement::Decorator<float> truthMatchProbX("truthMatchProbX");
    truthMatchProbX(*tp) = (*recoTrk_itr)->auxdata<float>( "truthMatchProbability");

    static SG::AuxElement::Decorator<xAOD::TruthParticleContainer::const_iterator> genMatchedIter("genMatchedIter");
    genMatchedIter(*tp) = genMatched_itr;

    static SG::AuxElement::Decorator<xAOD::TruthParticleContainer::const_iterator> genTargetMatchedIter("genTargetMatchedIter");
    genTargetMatchedIter(*tp) = genTargetMatched_itr;

    float alpha = 166.7820476;
    float rho = (*recoTrk_itr)->pt()/2.0/alpha;
    float d0  = (*recoTrk_itr)->d0();
    float phi = (*recoTrk_itr)->phi();
    float x0 = -(rho+d0)*TMath::Sin(phi);
    float y0 =  (rho+d0)*TMath::Cos(phi);

    int aIsign = 0;
    float aProd = pxQuark*x0+pyQuark*y0;
    if (aProd>0) { aIsign =  1; }
    if (aProd<0) { aIsign = -1; }

    static SG::AuxElement::Decorator<int> signToQuark("signToQuark");
    signToQuark(*tp) = aIsign;

    itkTrack->push_back( tp );
  }

  //---------------------
  // FILL PLOTS!!!
  //---------------------
  eventFeatures.nRecoTracks = numRecoTracks;
  eventFeatures.nRecoTracks1GeV = numRecoTracks1GeV;
  // eventFeatures.nPixelClusters = pixClustersOrig->size();
  // eventFeatures.nSCTClusters = sctClustersOrig->size();
  // eventFeatures.nSiClusters = pixClustersOrig->size() + sctClustersOrig->size();

  eventFeatures.nSiSpSeeds = 0;
  eventFeatures.nPseudoTracks = 0;
  eventFeatures.nPrimaryPseudoTracks   = 0;
  eventFeatures.nSecondaryPseudoTracks = 0;

  eventFeatures.average2TrackEta = 0;
  if ((itkTrack)->size()==2) {
    float aveEta = 0.0;
    for (xAOD::TrackParticleContainer::const_iterator itkTrk_itr=itkTrack->begin(); itkTrk_itr!=itkTrack->end(); itkTrk_itr++) {
      aveEta += (*itkTrk_itr)->eta();
    }
    eventFeatures.average2TrackEta = 0.5*aveEta;
  }

  // Search for conversion photon
  bool isPhoton = false;
  float phoEta = 100.0;
  float phoProdR  = -100.0;
  float phoDecayR = -100.0;
  for (xAOD::TruthParticleContainer::const_iterator truItr=truthPriParts->begin(); truItr!=truthPriParts->end(); truItr++) {
    if (TMath::Abs((*truItr)->pdgId())==11) {
      if ((*truItr)->hasProdVtx()) { 
        for (size_t iparent=0; iparent<(*truItr)->prodVtx()->nIncomingParticles(); iparent++) {
          const xAOD::TruthParticle* parent = (*truItr)->prodVtx()->incomingParticle(iparent);
          if (parent->pdgId()==22 && parent->pt()>10000.0) { isPhoton=true; phoEta=parent->eta(); phoProdR=parent->prodVtx()->perp(); phoDecayR=parent->decayVtx()->perp(); }
        }
      }
    }
  }
  if (isPhoton) {
    eventFeatures.photonEta = phoEta;
    eventFeatures.photonProdR = phoProdR;
    eventFeatures.photonDecayR = phoDecayR;
  }

  eventFeatures.nTruth = numTruth;
  eventFeatures.nChargedTruth = numChargedTruth;
  eventFeatures.nNeutralTruth = numNeutralTruth;

  eventFeatures.nTruth1GeV = numTruth1GeV;
  eventFeatures.nChargedTruth1GeV = numChargedTruth1GeV;
  eventFeatures.nNeutralTruth1GeV = numNeutralTruth1GeV;

  eventFeatures.nIsoTrack = nIsoTrk;
  eventFeatures.nVertex = numVtx;
  

  eventHist_all->FillHists( eventFeatures, 1.0 ); 

  for (xAOD::TruthParticleContainer::const_iterator truItr=truthPriParts->begin(); truItr!=truthPriParts->end(); truItr++) {
    if (TMath::Abs((*truItr)->pdgId())==m_idTarget && TMath::Abs((*truItr)->pt()-m_EnergyTarget)<0.1) {
      truHist_stable->FillHists( (*truItr), 1.0 );
    }
  }

  // Fill secondary vertex information
  for (xAOD::VertexContainer::const_iterator vtx_itr=primaryVertex->begin(); vtx_itr!=primaryVertex->end(); vtx_itr++) {
    if ((*vtx_itr)->vertexType()==3) {  // not a primary vertex
      vtxHist_secondary -> FillHists((*vtx_itr),1.0);
    }
  }

  // Fill selected tracks in the histogram
  xAOD::TrackParticleContainer::const_iterator itkTrk_itr_matched;
  float mindRMatched = 100.0;
  for (xAOD::TrackParticleContainer::const_iterator itkTrk_itr=itkTrack->begin(); itkTrk_itr!=itkTrack->end(); itkTrk_itr++) {
    trkHist_all -> FillHists((*itkTrk_itr),1.0);

    float mindR = (*itkTrk_itr)->auxdata<float>("matchedDR");
    if (mindR<mindRMatched) { itkTrk_itr_matched = itkTrk_itr; mindRMatched = mindR; }
  }
//  if (mindRMatched<0.02) {

  if (mindRMatched<0.1) {
    trkHist_reco->FillHists( (*itkTrk_itr_matched), 1.0 );
  }

  // Fill all pixel clusters
  xAOD::TrackMeasurementValidationContainer::const_iterator pixCluster_itr = pixClustersOrig->begin();
  xAOD::TrackMeasurementValidationContainer::const_iterator pixCluster_end = pixClustersOrig->end();
  for( ; pixCluster_itr != pixCluster_end; ++pixCluster_itr ) {
    clusHist_all -> FillHists( (*pixCluster_itr), 1.0 );
  }

  // Fill all strip clusters
  xAOD::TrackMeasurementValidationContainer::const_iterator sctCluster_itr = sctClustersOrig->begin();
  xAOD::TrackMeasurementValidationContainer::const_iterator sctCluster_end = sctClustersOrig->end();
  for( ; sctCluster_itr != sctCluster_end; ++sctCluster_itr ) {
    clusHist_all -> FillHists( (*sctCluster_itr), 1.0 );
  }

  //---------------------
  // CLEAN UP
  //---------------------
  delete truthPriParts;
  delete truthPriPartsAux;
  delete recoTracksShallowCopyPair.first;
  delete recoTracksShallowCopyPair.second;
  delete pixClustersPair.first;
  delete pixClustersPair.second;
  delete sctClustersPair.first;
  delete sctClustersPair.second;

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ITkStudy::postExecute() {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkStudy::finalize() {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkStudy::histFinalize() {
  runHist_reco->FillHists(1.0);
  return EL::StatusCode::SUCCESS;
}

