/******************************************
 *
 * Code to analysize physics samples of 
 * neutral particle decay into equal particles
 * like Z0 -> 2 mu or H -> 4 mu 
 *
 * G. Facini + modifications F.Massa
 * Wed Sep 10 07:21:58 CEST 2014
 *
 ******************************************/

//TO SET PARAMETERS SEE util/ITkPhysics.cxx

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ParticleAnalysis/ITkMultiParticle.h>
#include <Utility/DeltaFunctions.h>
#include <Utility/combinations.h>
#include <Utility/PhysicsEvent.h>
#include <HistManager/RunHists.h>
#include <HistManager/EventFeaturesInterface.h>
#include <Utility/TrackHelper.h>

// EDM includes: - if move to header file will not compile!
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

#include "xAODTracking/VertexContainer.h"

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include <TLorentzVector.h>

#include <map>
#include <limits>
#include <string>
#include <sstream>
#include <iostream>

// this is needed to distribute the algorithm to the workers
ClassImp(ITkMultiParticle)

ITkMultiParticle::ITkMultiParticle() {

}

EL::StatusCode ITkMultiParticle::setupJob(EL::Job& job) {

  job.useXAOD ();
  xAOD::Init( "PartAnaxAOD" ).ignore(); // call before opening first file

  // count number of events
  m_eventCounter = 0;
  std::cout << "Set counter to 0 " << m_eventCounter << std::endl;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkMultiParticle::histInitialize() {

  // In TruthHistManager constructor:
  // 1st bool is to fill all truth particles
  // 2nd bool is to fill gun truth particles
  // 3rd bool is to fill pileup truth particles

  // In TrackHistManager constructor:
  // 1st bool is to fill plots for tracks matched to primary particles
  // 2nd bool is to fill plots for tracks matched to secondary particles
  // 3rd bool is to fill plots for fake tracks 
  // 4th bool is if to fill duplicates (need 1st and 2nd be true)
  // 5th bool is if to fill plots based on hit content on track (need 1st be true)
  
  // recoTracks histograms

  trkHist_all       = new TrackHists("TrackAll"); //every track, no cuts here
  truHist_all       = new TruthHistManager("all"); //every truth status, no cuts
  truHist_hard      = new TruthHistManager("hard"); //
  truHist_pileup    = new TruthHistManager("pileup"); //
  truHist_all       -> Init(wk());
  truHist_hard      -> Init(wk());
  truHist_pileup    -> Init(wk());
  trkHist_matched      = new TrackHistManager("matched",true,true,true,true,true); //every track that matches the hs
  trkHist_reco = new TrackHistManager("reco",true,true,true,true,true); //every track that is recognized as part of signal
  trkHist_fakeCharge = new TrackHistManager("fakeCharge", true, true, true, true, true);
  clusHist_all      = new ClusterHists("SiliconAll");
  vtxHist_secondary = new VertexHists("Vertex");
  eventHist_all     = new EventHists("all");
  runHist_reco      = new RunHistManager(trkHist_reco, truHist_hard); //runHist reco is between target track and gun truth

  trkHist_all       -> Init(wk());
  trkHist_matched      -> Init(wk());
  trkHist_fakeCharge   -> Init(wk());
  trkHist_reco -> Init(wk());
  clusHist_all      -> Init(wk());
  vtxHist_secondary -> Init(wk());
  eventHist_all     -> Init(wk());
  runHist_reco      -> Init(wk());

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkMultiParticle::fileExecute () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkMultiParticle::changeInput(bool /* firstFile */) {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkMultiParticle::initialize() {

  m_event = wk()->xaodEvent();

  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int
  m_eventCounter=0;
	
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkMultiParticle::execute() {
  // print every 100 events, so we know where we are:
  if(m_eventCounter % 200 == 0) {
    Info("execute()", "Event number = %i", m_eventCounter );
  }

  m_eventCounter++;

  if (/*m_eventCounter == 11090 || m_eventCounter == 40465*/ false) 
    {
      return EL::StatusCode::SUCCESS;
    }
  //std::cout << "Event count " << m_eventCounter << std::endl;
  if(m_eventCounter > 1e5) { 
    if(m_eventCounter == 1e5) { 
      std::cout << "STOPPING after event " << m_eventCounter << std::endl;
      std::cout << "Look in here to change -> ITkMultiParticle :: execute()" << std::endl;
    }
    return EL::StatusCode::SUCCESS; 
  }

  EventFeatures eventFeatures;
  PhysicsEvent eventReco(m_nTarget, m_etaMax, m_decayMass, m_parentMass);

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;

  //---------------------------
  // Reset barcode map in histogram class
  //---------------------------
  trkHist_matched->resetBarcodes();
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
  int numPrimaryChargedTruth = 0;
  int numNeutralTruth = 0;
  int numTruth1GeV = 0;
  int numChargedTruth1GeV = 0;
  int numNeutralTruth1GeV = 0;

  double largeNumber = 10000000;

  //Create the new container and its auxiliary store for only primary truth particles 
  //need this since link to mother and children is not yet in the AR...
  double pxQuark = 0.0;
  double pyQuark = 0.0;
  double pzQuark = 0.0;
  xAOD::TruthParticleContainer* truthPriParts = new xAOD::TruthParticleContainer();
  xAOD::TruthParticleAuxContainer* truthPriPartsAux = new xAOD::TruthParticleAuxContainer();
  truthPriParts->setStore( truthPriPartsAux ); //< Connect the two

  std::vector<xAOD::TruthParticleContainer::const_iterator> genTargetMatched_itr;

  // ---------------------- ITERATION OVER TRUTH PARTICLES ---------------------------- //
  for (xAOD::TruthParticleContainer::const_iterator truthPart_itr=truthPartsOrig->begin(); truthPart_itr!=truthPartsOrig->end(); truthPart_itr++) {

    if (TMath::Abs((*truthPart_itr)->pdgId())<6 || (*truthPart_itr)->pdgId()==21) {
      pxQuark = (*truthPart_itr)->px();
      pyQuark = (*truthPart_itr)->py();
      pzQuark = (*truthPart_itr)->pz();
    }

    if ((*truthPart_itr)->status()!=1) { continue; }

    if (/*TMath::Abs((*truthPart_itr)->pdgId())!=13*/true) {
      numTruth++;
      if (TMath::Abs((*truthPart_itr)->charge())==1) { 
	numChargedTruth++;
	if ((*truthPart_itr)->barcode() < 200000)
	  numPrimaryChargedTruth++; 
      } //only primary (not interaction)
      if (TMath::Abs((*truthPart_itr)->charge())==0) { numNeutralTruth++; }

      if ((*truthPart_itr)->pt()>1000.0) {
        numTruth1GeV++;
        if (TMath::Abs((*truthPart_itr)->charge())==1) {numChargedTruth1GeV++;}
        if (TMath::Abs((*truthPart_itr)->charge())==0) { numNeutralTruth1GeV++; }
      }
    }
    if ((*truthPart_itr)->charge()==0) { continue; }

    //WARNING::truthPriParts will not contain truth outside detector eta range
    //if (TMath::Abs((*truthPart_itr)->eta()) > m_etaMax) {continue; }

    // Copy this particle to the new container:
    xAOD::TruthParticle* tp = new xAOD::TruthParticle();
    tp->makePrivateStore( **truthPart_itr );

    const double p = TMath::Sqrt(TMath::Power(tp->px(),2) +
				 TMath::Power(tp->py(),2) +
				 TMath::Power(tp->pz(),2));

    if (tp->hasProdVtx() && p > 1E-10) {
      static SG::AuxElement::Decorator<float> x0("x0");
      x0(*tp) = tp->prodVtx()->x();
      
      static SG::AuxElement::Decorator<float> y0("y0");
      y0(*tp) = tp->prodVtx()->y();
      
      static SG::AuxElement::Decorator<float> d0("d0");
      //d0(*tp) = TMath::Sign(tp->prodVtx()->perp(),
      //			    tp->prodVtx()->x()*tp->py() - tp->prodVtx()->y()*tp->px());
      d0(*tp) = 0.0;
      
      static SG::AuxElement::Decorator<float> z0("z0");
      z0(*tp) = tp->prodVtx()->z();
    }


    truthPriParts->push_back( tp );    
  }

  // ---------------------- END ITERATION OVER TRUTH PARTICLES ---------------------------- //

  // Find decay products candidates of neutral parent

  bool isValidMatching = true; //if in the event a truth target goes over the geometrical acceptance the tracks should not be matched

  for (xAOD::TruthParticleContainer::const_iterator truItr=truthPriParts->begin(); truItr!=truthPriParts->end(); truItr++) {
    if (TMath::Abs((*truItr)->pdgId()) == m_idTarget
	/* && (*truItr)->pt() > m_ptCut */)
      genTargetMatched_itr.push_back(truItr);
  }

  ///////////// Now find the right decay products if not already /////////////

  if (genTargetMatched_itr.size() == m_nTarget) {
    //it's the expected number of particles, check if the invariant mass is right
    TLorentzVector momentum(0.0,0.0,0.0,0.0);
    double totalCharge = 0.0;
    for (unsigned int i = 0; i < m_nTarget; i++) {
      momentum += (*(genTargetMatched_itr[i]))->p4();
      totalCharge += (*(genTargetMatched_itr[i]))->charge();
    }
    //parent mass and total charge check
    if (TMath::Abs(totalCharge) > 1E-6) {
      std::cout << "ERROR\t found " << m_nTarget << " particles but not corresponding to parent" << std::endl;
      exit(1);
    }
  }
  
  else if (genTargetMatched_itr.size() < m_nTarget) {
    //there are less particles than expected. Error: parent not found (cut too tight?)
    std::cout << "ERROR\t found less particles than expected. Cut too tight? Exiting." << std::endl;
    exit(1);
  }

  else if (genTargetMatched_itr.size() > m_nTarget) {
    //there are more particles than expected. Cut too loose? Look for the right subset

    //store all combinations
    std::vector<std::vector<xAOD::TruthParticleContainer::const_iterator> > genTargetCombinations = 
      combinations<xAOD::TruthParticleContainer::const_iterator>(genTargetMatched_itr, m_nTarget);

    std::vector<std::vector<xAOD::TruthParticleContainer::const_iterator> >::iterator bestCombination_itr;

    //now choose the one nearest to parent mass
    double minResidual = largeNumber;
    for (std::vector<std::vector<xAOD::TruthParticleContainer::const_iterator> >::iterator combItr = genTargetCombinations.begin();
	 combItr != genTargetCombinations.end(); combItr++) {
      TLorentzVector momentum(0.0,0.0,0.0,0.0);
      double totalCharge = 0.0;
      for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator genItr = (*combItr).begin(); genItr != (*combItr).end(); genItr++) {
	momentum += (**genItr)->p4();
	totalCharge += (**genItr)->charge();
      } 
      
      if (TMath::Abs(totalCharge) < 1E-6) {
	double residual = TMath::Abs(momentum.M() - m_parentMass);
	if (residual < minResidual) {
	  minResidual = residual;
	  bestCombination_itr = combItr;
	}
      }

    }

    if (TMath::Abs(minResidual - largeNumber) > 1E-6) {
    //now that we have the best combination we overwrite the vector containing the matched particles (may reconsider this for fake study)
      genTargetMatched_itr.clear();
      for (unsigned int i = 0; i < m_nTarget; i++) 
	genTargetMatched_itr.push_back((*bestCombination_itr).at(i));
      
    }
    else {
      std::cout << "ERROR in ITkMultiParticle::execute(): No truth combination valid. Exiting" << std::endl;
      exit(1);
    }

  }

  eventReco.SetTruth(genTargetMatched_itr);

  for (unsigned int i = 0; i < m_nTarget; i++) {
    isValidMatching = isValidMatching && (TMath::Abs((*(genTargetMatched_itr[i]))->eta()) < m_etaMax);
  }

  //here put every event info about truth particles
  /*  eventFeatures.truthMass = eventReco.GetTruthMass();
  eventFeatures.truthPt   = eventReco.GetTruthP4().Pt();
  eventFeatures.truthMaxPt = eventReco.GetTruthMaxPt();
  eventFeatures.truthMinPt = eventReco.GetTruthMinPt();
  eventFeatures.truthMaxEtaPt = eventReco.GetTruthMaxEtaPt();
  eventFeatures.truthMaxEta = eventReco.GetTruthMaxEta();
  eventFeatures.truthMinEta = eventReco.GetTruthMinEta();
  eventFeatures.truthCharge = eventReco.GetTruthCharge();
  eventFeatures.truthSmallestDR = eventReco.GetTruthSmallestDR();
  eventFeatures.truthOutsideDetector = eventReco.GetTruthOutsideDetector();
  */

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
    xAOD::TruthParticleContainer::const_iterator genTargetMatched_itr_element;

    float mindR = largeNumber;
    float matchedR = largeNumber;
    for (xAOD::TruthParticleContainer::const_iterator genPrt_itr=truthPartsOrig->begin(); genPrt_itr!=truthPartsOrig->end(); genPrt_itr++) {
      if((*genPrt_itr)->barcode() <= 10e3 || (*genPrt_itr)->barcode() >= 200e3) { continue; }
      if ((*genPrt_itr)->status()!=1) { continue; }
      if ((*genPrt_itr)->pdgId()==22) { continue; }
      float dR = deltaR((*genPrt_itr)->phi(),(*recoTrk_itr)->phi(),(*genPrt_itr)->eta(),(*recoTrk_itr)->eta());
      if (dR<mindR) { genMatched_itr = genPrt_itr; mindR = dR; }
      //      if (TMath::Abs((*genPrt_itr)->pdgId())==m_idTarget && TMath::Abs((*genPrt_itr)->pt()-m_EnergyTarget)<0.1) {
      //  genTargetMatched_itr = genPrt_itr; matchedR = dR;
      //  }
    }

    //now find the R distance between the track and the closest target truth
    for (unsigned int i = 0; i < genTargetMatched_itr.size(); i++) {
      float dR = deltaR((*(genTargetMatched_itr[i]))->phi(), (*recoTrk_itr)->phi(), (*(genTargetMatched_itr[i]))->eta(), (*recoTrk_itr)->eta());
      if (dR < matchedR) {
	matchedR = dR;
	genTargetMatched_itr_element = genTargetMatched_itr[i];
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
    genTargetMatchedIter(*tp) = genTargetMatched_itr_element;

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
  } //track for

  ///////////////////////////////////////////////////////////////////////////
  //matched = tracks associated with truth minimizing dR distance
  // (candidates) = true simulation of analysis, look for combination that more approximates parentMass

  /////////////// Now find the tracks corresponding to every truth matched and fill a pair vector///
  std::vector<xAOD::TrackParticleContainer::const_iterator> itkTrk_matched;
  std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> > truth_trackMatched;
  std::vector<xAOD::TrackParticleContainer::const_iterator> itkTrk_itr_candidates;

  if (isValidMatching) {
    for (unsigned int i = 0; i < genTargetMatched_itr.size(); i++) {
      float mindR = largeNumber;
      xAOD::TrackParticleContainer::const_iterator itkTrk_itr_matched;
      for (xAOD::TrackParticleContainer::const_iterator itkTrk_itr=itkTrack->begin(); itkTrk_itr!=itkTrack->end(); itkTrk_itr++) {
	float dR = deltaR((*(genTargetMatched_itr[i]))->phi(), (*itkTrk_itr)->phi(), (*(genTargetMatched_itr[i]))->eta(), (*itkTrk_itr)->eta());
	if (dR < mindR) { 
	  itkTrk_itr_matched = itkTrk_itr;
	  mindR = dR;
	}  
      }
      itkTrk_matched.push_back(itkTrk_itr_matched);
      truth_trackMatched.push_back(std::make_pair(genTargetMatched_itr[i],itkTrk_itr_matched));
    }

    eventReco.SetPair(truth_trackMatched);
    eventReco.SetMatched(itkTrk_matched);
  }

 

  for (xAOD::TrackParticleContainer::const_iterator itkTrk_itr=itkTrack->begin(); itkTrk_itr!=itkTrack->end(); itkTrk_itr++) {
    //parent reconstruction without truth matching
    //const xAOD::TruthParticle* truthParticle = xAOD::TrackHelper::truthParticle(*itkTrk_itr);
    //insert all reconstruction cuts over single track here
    if ((*itkTrk_itr)->pt() > m_ptCut ) {
      uint8_t getInt = 0;
      (*itkTrk_itr)->summaryValue(getInt, xAOD::numberOfPixelHits);
      int nPixHits = getInt;
      (*itkTrk_itr)->summaryValue(getInt, xAOD::numberOfSCTHits);
      int nSCTHits = getInt;

      if ((nPixHits + nSCTHits) >= m_hitCut) {
	itkTrk_itr_candidates.push_back(itkTrk_itr);
      }
    }
  }

  bool isValidCandidate = true; //isValidCandidate == true when there are a number of tracks >= m_nTarget 


  if (itkTrk_itr_candidates.size() == m_nTarget) {
    double totalCharge = 0.0;
    TLorentzVector candidateMomentum(0.0,0.0,0.0,0.0);
    for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator candItr = itkTrk_itr_candidates.begin();
	 candItr != itkTrk_itr_candidates.end();
	 candItr++) {
      totalCharge += (**candItr)->charge();
      double px = (**candItr)->pt()*TMath::Cos((**candItr)->phi());
      double py = (**candItr)->pt()*TMath::Sin((**candItr)->phi());
      double pz = (**candItr)->pt()*TMath::SinH((**candItr)->eta());
      double e = TMath::Sqrt(TMath::Power(px,2)+TMath::Power(py,2)+TMath::Power(pz,2) + TMath::Power(m_decayMass,2));
      TLorentzVector tmp(px,py,pz,e);
      candidateMomentum += tmp;
    }

    //apply here cuts on set of candidates
    if (TMath::Abs(totalCharge) > 1E-6 ||
	TMath::Abs(candidateMomentum.M() - m_parentMass) > m_parentMassTolerance)
      isValidCandidate = false;
  }

  else if (itkTrk_itr_candidates.size() < m_nTarget) {
    isValidCandidate = false;
  }

  else if (itkTrk_itr_candidates.size() > m_nTarget) {
    std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> > combVector = 
      combinations<xAOD::TrackParticleContainer::const_iterator> (itkTrk_itr_candidates, m_nTarget);

    std::vector<xAOD::TrackParticleContainer::const_iterator> best_recoTrk_itr;

    double minResidual = largeNumber;
    for (std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator combItr = combVector.begin();
	 combItr != combVector.end(); combItr++) {
      TLorentzVector combMomentum(0.0,0.0,0.0,0.0);
      double totalCharge = 0.0;

      for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator combElement = (*combItr).begin();
	   combElement != (*combItr).end(); combElement++) {
	double px = (**combElement)->pt()*TMath::Cos((**combElement)->phi());
	double py = (**combElement)->pt()*TMath::Sin((**combElement)->phi());
	double pz = (**combElement)->pt()*TMath::SinH((**combElement)->eta());
	double e = TMath::Sqrt(TMath::Power(px,2) + TMath::Power(py,2) + TMath::Power(pz,2) + TMath::Power(m_decayMass,2));
	totalCharge += (**combElement)->charge();
	TLorentzVector tmp(px,py,pz,e);
	combMomentum += tmp;
      }

      if (TMath::Abs(totalCharge) < 1E-6) {
	double candidateParentMass = combMomentum.M();
	double residual = TMath::Abs(candidateParentMass - m_parentMass);
	//apply here cuts on candidate set
	if (residual < minResidual &&
	    residual < m_parentMassTolerance) {
	  minResidual = residual;
	  best_recoTrk_itr = (*combItr);
	}
      }

    }

    if (TMath::Abs(minResidual - largeNumber) > 1E-6)
      itkTrk_itr_candidates = best_recoTrk_itr;
    else
      isValidCandidate = false;
  }
  
  
  ////////////////////////////////////////////////////////////////////
  // now look into candidates for parent decay and zero charge
  //  int nCount = 0;

  //only for Higgs compatibility
  std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> > dummyPair;

  if (isValidCandidate) {
    eventReco.SetReco(itkTrk_itr_candidates, dummyPair);
    for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator candidate_itr = itkTrk_itr_candidates.begin();
	 candidate_itr != itkTrk_itr_candidates.end();
	 candidate_itr++)
      trkHist_reco -> FillHists(**candidate_itr, 1.0);
  
  
    //here put every event info on reco tracks
    /*    eventFeatures.recoMass = eventReco.GetRecoMass();
    eventFeatures.recoCharge = eventReco.GetRecoCharge();
    eventFeatures.recoTruthCharge = eventReco.GetRecoTruthCharge();
    eventFeatures.recoPt = eventReco.GetRecoP4().Pt();
    eventFeatures.recoTruthPt = eventReco.GetRecoTruthP4().Pt();
    eventFeatures.recoMaxPt = eventReco.GetRecoMaxPt();
    eventFeatures.recoMinPt = eventReco.GetRecoMinPt();
    eventFeatures.recoMaxEta = eventReco.GetRecoMaxEta();
    eventFeatures.recoMinEta = eventReco.GetRecoMinEta();
    eventFeatures.recoSmallestDR = eventReco.GetRecoSmallestDR();
    */

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
  eventFeatures.nPrimaryChargedTruth = numPrimaryChargedTruth;
  eventFeatures.nNeutralTruth = numNeutralTruth;

  eventFeatures.nTruth1GeV = numTruth1GeV;
  eventFeatures.nChargedTruth1GeV = numChargedTruth1GeV;
  eventFeatures.nNeutralTruth1GeV = numNeutralTruth1GeV;

  eventFeatures.nIsoTrack = nIsoTrk;
  eventFeatures.nVertex = numVtx;



  for (xAOD::TruthParticleContainer::const_iterator truItr=truthPriParts->begin(); truItr!=truthPriParts->end(); truItr++) {
    bool isSignalTruth = false;
    for (unsigned int i = 0; i < genTargetMatched_itr.size(); i++) {
      isSignalTruth = isSignalTruth || ((truItr) == (genTargetMatched_itr[i]));
    }
    truHist_all->FillHists((*truItr),1.0);
    if (isSignalTruth) 
      truHist_hard->FillHists((*truItr),1.0);
    else
      truHist_pileup->FillHists((*truItr),1.0);
  }

  // Fill secondary vertex information
  for (xAOD::VertexContainer::const_iterator vtx_itr=primaryVertex->begin(); vtx_itr!=primaryVertex->end(); vtx_itr++) {
    if ((*vtx_itr)->vertexType()==3) {  // not a primary vertex
      vtxHist_secondary -> FillHists((*vtx_itr),1.0);
    }
  }

  // Fill selected tracks in the histogram
  for (xAOD::TrackParticleContainer::const_iterator itkTrk_itr=itkTrack->begin(); itkTrk_itr!=itkTrack->end(); itkTrk_itr++) {
    trkHist_all -> FillHists((*itkTrk_itr),1.0);
  }
//  if (mindRMatched<0.02) {

  uint8_t getInt(0);

  for (std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator,xAOD::TrackParticleContainer::const_iterator> >::iterator itr = 
	 truth_trackMatched.begin(); itr != truth_trackMatched.end(); itr++) {
      trkHist_matched->FillHists(*((*itr).second), 1.0 );
      if (eventReco.isMatchedSet())
	if (TMath::Abs(eventReco.GetMatchedCharge()) > 1E-6)
	  trkHist_fakeCharge -> FillHists(*((*itr).second), 1.0);
      
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

  /*
  //here put every event info about matched tracks
  if (isValidMatching) {
    eventFeatures.matchedMass = eventReco.GetMatchedMass();
    eventFeatures.matchedMaxPt = eventReco.GetMatchedMaxPt();
    eventFeatures.matchedMinPt = eventReco.GetMatchedMinPt();
    eventFeatures.matchedMaxEta = eventReco.GetMatchedMaxEta();
    eventFeatures.matchedMaxEtaPt = eventReco.GetMatchedMaxEtaPt();
    eventFeatures.matchedMinEta = eventReco.GetMatchedMinEta();
    eventFeatures.matchedCharge = eventReco.GetMatchedCharge();
    eventFeatures.matchedSmallestDR = eventReco.GetMatchedSmallestDR();
  }
  */

  eventFeatures.physicsEvent = eventReco;
  eventHist_all->FillHists( eventFeatures, 1.0 ); 

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


EL::StatusCode ITkMultiParticle::postExecute() {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkMultiParticle::finalize() {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ITkMultiParticle::histFinalize() {
  runHist_reco->FillHists(1.0);
  return EL::StatusCode::SUCCESS;
}

