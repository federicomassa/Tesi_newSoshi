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

//TO SET PARAMETERS SEE util/ITkHiggs.cxx

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ParticleAnalysis/HiggsAnalysis.h>
#include <Utility/DeltaFunctions.h>
#include <Utility/combinations.h>
#include <Utility/PhysicsEvent.h>
#include <Utility/InvariantMass.h>
#include <Utility/TrackHelper.h>
#include <Utility/HiggsEventSelector.h>

#include <HistManager/RunHists.h>
#include <HistManager/EventFeaturesInterface.h>

// EDM includes: - if move to header file will not compile!
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

#include "xAODTracking/VertexContainer.h"

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include <TLorentzVector.h>
#include <set>
#include <map>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <iostream>

// this is needed to distribute the algorithm to the workers
ClassImp(HiggsAnalysis)

void SortIterate(std::vector<xAOD::TrackParticleContainer::const_iterator>& v, std::vector<xAOD::TrackParticleContainer::const_iterator>& copy) {

  if (v.size() == 0)
    return;

  double maxP = -1;
  std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator maxP_itr;

  for (auto itr = v.begin();
       itr != v.end();
       itr++) {
    if ((**itr)->pt() > maxP) {
      maxP = ((**itr)->pt());
      maxP_itr = itr;
    }
  }

  copy.push_back(*maxP_itr);
  v.erase(maxP_itr);

  SortIterate(v,copy);
  
}


void SortVector(std::vector<xAOD::TrackParticleContainer::const_iterator>& v) {
  
  std::vector<xAOD::TrackParticleContainer::const_iterator> copy;

  SortIterate(v,copy);

  v = copy;

}



HiggsAnalysis::HiggsAnalysis() {

}

EL::StatusCode HiggsAnalysis::setupJob(EL::Job& job) {

  job.useXAOD ();
  xAOD::Init( "PartAnaxAOD" ).ignore(); // call before opening first file

  // count number of events
  m_eventCounter = 0;
  m_moreTargets = 0;
  std::cout << "Set counter to 0 " << m_eventCounter << std::endl;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsAnalysis::histInitialize() {

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
  truHist_photons    = new TruthHistManager("photons"); //
  truHist_visibleHard = new TruthHistManager("visibleHard"); // eta + kinematics cut

  truHist_all       -> Init(wk());
  truHist_hard      -> Init(wk());
  truHist_pileup    -> Init(wk());
  truHist_photons    -> Init(wk());
  truHist_visibleHard -> Init(wk());
  trkHist_matched      = new TrackHistManager("matched",true,true,true,false,false); //every track that matches the hs
  trkHist_reco = new TrackHistManager("reco",true,true,true,false,false); //every track that is recognized as part of signal
  trkHist_pileup = new TrackHistManager("pileup",true,true,true,false,false); //every track that is recognized as part of signal
  trkHist_fakeCharge = new TrackHistManager("fakeCharge", true, true, true, false,false);
  clusHist_all      = new ClusterHists("SiliconAll");
  vtxHist_secondary = new VertexHists("Vertex");
  eventHist_all     = new EventHists("all");
  runHist_reco      = new RunHistManager(trkHist_reco, truHist_hard); //runHist reco is between target track and gun truth

  trkHist_all       -> Init(wk());
  trkHist_matched      -> Init(wk());
  trkHist_fakeCharge   -> Init(wk());
  trkHist_reco -> Init(wk());
  trkHist_pileup -> Init(wk());
  clusHist_all      -> Init(wk());
  vtxHist_secondary -> Init(wk());
  eventHist_all     -> Init(wk());
  runHist_reco      -> Init(wk());

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsAnalysis::fileExecute () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsAnalysis::changeInput(bool /* firstFile */) {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsAnalysis::initialize() {

  m_event = wk()->xaodEvent();

  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int
  m_eventCounter=0;
	
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsAnalysis::execute() {
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
  if(m_eventCounter > 1E5) { 
    if(m_eventCounter == 1E5) { 
      std::cout << "STOPPING after event " << m_eventCounter << std::endl;
      std::cout << "Look in here to change -> HiggsAnalysis :: execute()" << std::endl;
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
  trkHist_pileup->resetBarcodes();

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
  int numPrimaryChargedTruth1GeV = 0;
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
	if ((*truthPart_itr)->barcode() < 200000) {
	  numPrimaryChargedTruth++; 
	  if ((*truthPart_itr)->pt() > 1000.0 && TMath::Abs((*truthPart_itr)->eta()) < 4.0)
	    numPrimaryChargedTruth1GeV++;
	}
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

  double vtxZ = -1000;

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

    std::vector<const xAOD::TruthVertex*> vtx;
  
    for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator genItr = (genTargetMatched_itr).begin(); genItr != (genTargetMatched_itr).end(); genItr++) {
      vtx.push_back((**genItr)->prodVtx());
      if (vtx.back() == 0) {
	std::cout << "Found 4 muons but without vertex" << std::endl;
	return EL::StatusCode::SUCCESS;
      }
      
      vtxZ = vtx[0]->z();
      for (auto itr = vtx.begin() + 1;
	   itr != vtx.end();
	   itr++) {
	if (TMath::Abs((*itr)->z() - vtxZ) > 1E-6) {
	  std::cout << "Found 4 muons but not from the same vertex" << std::endl;
	  return EL::StatusCode::SUCCESS;
	}
      }

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
      
      //same z constaint
      std::vector<const xAOD::TruthVertex*> vtx;
      bool isCombInvalid = false;

      for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator genItr = (*combItr).begin(); genItr != (*combItr).end(); genItr++) {
	vtx.push_back((**genItr)->prodVtx());
	if (vtx.back() == 0) {
	  isCombInvalid = true;
	  break;
	}
	momentum += (**genItr)->p4();
	totalCharge += (**genItr)->charge();
      }
      
      if (isCombInvalid)
	continue;

      vtxZ = vtx[0]->z();
      for (auto itr = vtx.begin() + 1;
	   itr != vtx.end();
	   itr++) {
	if (TMath::Abs((*itr)->z() - vtxZ) > 1E-6)
	  isCombInvalid = true;
      }
      

      if (isCombInvalid)
	continue;


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
      std::cout << "ERROR in HiggsAnalysis::execute(): No truth combination valid. Exiting" << std::endl;
      exit(1);
    }

  }


  /*  if (TMath::Abs(InvariantMass::TruthInvariantMass(genTargetMatched_itr) - (50000+160000)/2.0) > (160000-50000)/2.0)
    return EL::StatusCode::SUCCESS;
  */

  //photons?
  std::vector<const xAOD::TruthParticle*> primaryPhotons;
  const xAOD::TruthVertex* primaryVtx = (*(genTargetMatched_itr[0]))->prodVtx(); 
  double primaryVtxZ = primaryVtx->z();

  if (primaryVtx == NULL) {
    std::cout << "In HiggsAnalysis::execute(): primary vertex not existent?" << std::endl;
    exit(1);
  }

  if (primaryVtx->nOutgoingParticles() > m_nTarget) {
   
    for (int i = 0; 
	 i < primaryVtx->nOutgoingParticles();
	 i++) {

      const xAOD::TruthParticle* truthPart = primaryVtx->outgoingParticle(i);

      if (truthPart->pdgId() == 22) {
	primaryPhotons.push_back(truthPart);
      }
    }
  }



  eventReco.SetPhotons(primaryPhotons);
  eventReco.SetTruth(genTargetMatched_itr);

  // if (TMath::Abs(eventReco.GetTruthP4WithPhotons().M() - 125000) > 1000) {

  //   std::cout << "N photons:  " << '\t' << primaryPhotons.size() << std::endl;
  //   std::cout << "N particles:" << '\t' << primaryVtx->nOutgoingParticles() << std::endl;
  // }

  // else {
  //   return EL::StatusCode::SUCCESS;
  // }

  for (auto itr = primaryPhotons.begin();
       itr != primaryPhotons.end();
       itr++) 
    truHist_photons -> FillHists(*itr, 1.0);

  /*
       //     on shell cuts
  isValidMatching = isValidMatching && ((*(genTargetMatched_itr[0]))->pt() > m_truthPtCut[0]);
  isValidMatching = isValidMatching && ((*(genTargetMatched_itr[1]))->pt() > m_truthPtCut[1]);
       //     off shell cuts (both > [3] and at least one > [2]) 
  isValidMatching = isValidMatching && ((*(genTargetMatched_itr[2]))->pt() > m_truthPtCut[3]);
  isValidMatching = isValidMatching && ((*(genTargetMatched_itr[3]))->pt() > m_truthPtCut[3]);
  
  bool atLeastOneOverCut = false;
  atLeastOneOverCut = atLeastOneOverCut || ((*(genTargetMatched_itr[2]))->pt() > m_truthPtCut[2]);
  atLeastOneOverCut = atLeastOneOverCut || ((*(genTargetMatched_itr[3]))->pt() > m_truthPtCut[2]);

  isValidMatching = isValidMatching && atLeastOneOverCut;

  */

  std::vector<double> rankedPt;
  for (auto itr = genTargetMatched_itr.begin();
       itr != genTargetMatched_itr.end();
       itr++) {
    rankedPt.push_back((**itr)->pt());
  }

  std::sort(rankedPt.begin(), rankedPt.end());

  std::vector<xAOD::TruthParticleContainer::const_iterator> rankedTruth;

  for (auto rankedItr = rankedPt.begin();
       rankedItr != rankedPt.end();
       rankedItr++) {
    for (auto truthItr = genTargetMatched_itr.begin();
	 truthItr != genTargetMatched_itr.end();
	 truthItr++) {
     
      if (TMath::Abs((**truthItr)->pt() - *rankedItr) < 1E-6) {
	rankedTruth.push_back(*truthItr);
	continue;
      }

    }
  }

  //cut code: eta1 eta2 eta3 eta4 pt1 pt2 pt3 pt4, where 1,2,3,4 are pt in ascending order: 00000000 all cuts are passed
  bool eta1Cut = TMath::Abs((*(rankedTruth[0]))->eta()) > m_etaMax;
  bool pt1Cut  = (*(rankedTruth[0]))->pt() < m_truthPtCut[0];
  bool eta2Cut = TMath::Abs((*(rankedTruth[1]))->eta()) > m_etaMax;
  bool pt2Cut  = (*(rankedTruth[1]))->pt() < m_truthPtCut[1];
  bool eta3Cut = TMath::Abs((*(rankedTruth[2]))->eta()) > m_etaMax;
  bool pt3Cut  = (*(rankedTruth[2]))->pt() < m_truthPtCut[2];
  bool eta4Cut = TMath::Abs((*(rankedTruth[3]))->eta()) > m_etaMax;
  bool pt4Cut  = (*(rankedTruth[3]))->pt() < m_truthPtCut[3];

  int nCentered = 0;

  for (auto itr = genTargetMatched_itr.begin();
       itr != genTargetMatched_itr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) < 2.5 &&
	TMath::Abs((**itr)->eta()) > 0.1)
      nCentered++;
      }

  int truthCutCode = 
    pt4Cut*TMath::Power(2,0) +
    pt3Cut*TMath::Power(2,1) +
    pt2Cut*TMath::Power(2,2) +
    pt1Cut*TMath::Power(2,3) +
    eta4Cut*TMath::Power(2,4) +
    eta3Cut*TMath::Power(2,5) +
    eta2Cut*TMath::Power(2,6) +
    eta1Cut*TMath::Power(2,7);

  //no truth cut for the time being

  if (truthCutCode < 16)
    isValidMatching = true;
  else
    isValidMatching = false;

  eventReco.truthCutCode = truthCutCode;

  isValidMatching = isValidMatching && (nCentered >= 0);

  if (isValidMatching) {

    for (auto itr = genTargetMatched_itr.begin();
	 itr != genTargetMatched_itr.end();
	 itr++)
      truHist_visibleHard -> FillHists(**itr, 1.0);
    eventReco.isTruthVisible = true;
  }

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

 
  //INITIALLY WITHOUT CUTS
  for (xAOD::TrackParticleContainer::const_iterator itkTrk_itr=itkTrack->begin(); itkTrk_itr!=itkTrack->end(); itkTrk_itr++) {
    //parent reconstruction without truth matching
    //const xAOD::TruthParticle* truthParticle = xAOD::TrackHelper::truthParticle(*itkTrk_itr);
    //insert all reconstruction cuts over single track here
    if ((*itkTrk_itr)->pt() > 0 ) {
      uint8_t getInt = 0;
      (*itkTrk_itr)->summaryValue(getInt, xAOD::numberOfPixelHits);
      int nPixHits = getInt;
      (*itkTrk_itr)->summaryValue(getInt, xAOD::numberOfSCTHits);
      int nSCTHits = getInt;

      
      //change here the analysis eta coverage
      if ((nPixHits + nSCTHits) >= 0 && TMath::Abs((*itkTrk_itr)->eta()) < 4.0) {
	itkTrk_itr_candidates.push_back(itkTrk_itr);
      }
    }
  }

  //in this version of the program the recoCutCode is controlled here, it is 
  //only updated when it does not pass a cut

  if (itkTrk_itr_candidates.size() >= 4)
    eventReco.recoCutCode = 0;
  else {
    eventReco.recoCutCode = -1;
  }

  SortVector(itkTrk_itr_candidates);  

  // NOW ADD HIT CUT
  for (auto itr = itkTrk_itr_candidates.begin();
       itr != itkTrk_itr_candidates.end();
       itr++) {
    //parent reconstruction without truth matching
    //const xAOD::TruthParticle* truthParticle = xAOD::TrackHelper::truthParticle(*itkTrk_itr);
    //insert all reconstruction cuts over single track here
      uint8_t getInt = 0;
      (**itr)->summaryValue(getInt, xAOD::numberOfPixelHits);
      int nPixHits = getInt;
      (**itr)->summaryValue(getInt, xAOD::numberOfSCTHits);
      int nSCTHits = getInt;

      if ((nPixHits + nSCTHits) < m_hitCut) {
	itkTrk_itr_candidates.erase(itr);
	itr--;
      }
    
  }

  if (itkTrk_itr_candidates.size() < 4 && eventReco.recoCutCode == 0) {
    eventReco.recoCutCode = 1; 
  }

  //NOW MIN PT CUT
  for (auto itr = itkTrk_itr_candidates.begin();
       itr != itkTrk_itr_candidates.end();
       itr++) {
    if ((**itr)->pt() < m_truthPtCut[0]) {
	itkTrk_itr_candidates.erase(itr);
	itr--;
      }
    
  }

  if (itkTrk_itr_candidates.size() < 4 && eventReco.recoCutCode == 0) {
    eventReco.recoCutCode = 2; 
  }


  bool isValidCandidate = true; //selection on reco tracks

  //NOW WE TRY TO DECREASE THE RATE OF EVENTS IN WHICH MORE THAN 4 TRACKS PASS THE CUTS
  bool moreThan4 = false;

  if (itkTrk_itr_candidates.size() > 4) {
    std::cout << "More than 4" << std::endl;
    moreThan4 = true;
  }

  ///////////// 1. DELTA Z //////////////////
  //we select only the tracks that are inside a region of 3*sigma_z around the primaryVtx
  //As the sigma_z resolution depend heavily on the eta of the track, we divide the resolution
  //into eta regions... sigma_z = 50 um for eta < 2.7; = 150 um for 2.7 < eta < 3.2; = 500 um for 3.2 < eta < 4.0

  double deltaZCut27 = 5*0.05; //mm
  double deltaZCut32 = 5*0.15; //mm
  double deltaZCut40 = 5*0.5; //mm

  for (auto itr = itkTrk_itr_candidates.begin();
       itr != itkTrk_itr_candidates.end();
       itr++) {

    if (TMath::Abs((**itr)->eta()) < 2.7) {
      if (TMath::Abs((**itr)->z0() - primaryVtxZ) > deltaZCut27) {
	itkTrk_itr_candidates.erase(itr);
	itr--;
      }
    }

    else if (TMath::Abs((**itr)->eta()) < 3.2) {
      if (TMath::Abs((**itr)->z0() - primaryVtxZ) > deltaZCut32) {
	itkTrk_itr_candidates.erase(itr);
	itr--;
      }
    }

    else if (TMath::Abs((**itr)->eta()) < 4.0) {
      if (TMath::Abs((**itr)->z0() - primaryVtxZ) > deltaZCut40) {
	itkTrk_itr_candidates.erase(itr);
	itr--;
      }
    }

    else {
      std::cout << "HiggsAnalysis::Execute(): TRACK BEYOND ETA 4?" << std::endl;
      exit(1);
    }

  }

  if (itkTrk_itr_candidates.size() < 4 && eventReco.recoCutCode == 0)
    eventReco.recoCutCode = 3;

  //   if (itkTrk_itr_candidates.size() > 4) 
  //    std::cout << "Still more than 4" << std::endl;

  
  // else if (moreThan4 && itkTrk_itr_candidates.size() <= 4 && itkTrk_itr_candidates.size() != 0) 
  //   std::cout << "Now: " << itkTrk_itr_candidates.size() << '\t' << "primary is: " << primaryVtxZ << std::endl;

  // else if (moreThan4 && itkTrk_itr_candidates.size() == 0) {
  //   std::cout << "NOOOW: " << itkTrk_itr_candidates.size() << '\t' << "primary is: " << primaryVtxZ << std::endl;
  // }

  
  ////// 2. PT ISOLATION ///////////
  //isolation cut: sum(pt in dR < 0.1) / pt(hard track)
  std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> > recoIsolation;
  for (auto candidateItr = itkTrk_itr_candidates.begin();
       candidateItr != itkTrk_itr_candidates.end();
       candidateItr++) {
    double ptSum = 0;
  
    for (auto tracks = itkTrack->begin();
	 tracks != itkTrack->end(); 
	 tracks++) {
      
      //if it's the muon candidate do not consider it
      if (*candidateItr == tracks)
	continue;

      if (deltaR((**candidateItr)->phi(), (*tracks)->phi(), (**candidateItr)->eta(), (*tracks)->eta()) < 0.1)
	ptSum += (*tracks)->pt();
      
    }
    
    if (ptSum/(**candidateItr)->pt() > 1) {
      itkTrk_itr_candidates.erase(candidateItr);
      candidateItr--;
    }

    else
      recoIsolation.push_back(std::make_pair(*candidateItr,ptSum/(**candidateItr)->pt()));

  }

  if (itkTrk_itr_candidates.size() < 4 && eventReco.recoCutCode == 0)
    eventReco.recoCutCode = 4;

  //NOW THERE IS NOTHING ELSE WE CAN DO, IF STILL > 4 TARGETS THEN COMBINATORY PROBLEM...
  
  if (itkTrk_itr_candidates.size() > 4) 
    std::cout << "STILL MORE THAN 4!!!!!" << std::endl;
 
  if (moreThan4 && itkTrk_itr_candidates.size() <= 4) 
    std::cout << "\t ----> " << itkTrk_itr_candidates.size() << std::endl;


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
        if (TMath::Abs(totalCharge) > 1E-6 || false
	    /*	TMath::Abs(candidateMomentum.M() - m_parentMass) > m_parentMassTolerance*/)
	  isValidCandidate = false;
	
  }

  else if (itkTrk_itr_candidates.size() < m_nTarget) {
    // if (isValidMatching)
    //   std::cout << "< Target" << std::endl;
    isValidCandidate = false;
  }

  else if (itkTrk_itr_candidates.size() > m_nTarget) {

    m_moreTargets++;
    // std::cout << "More than 4 candidates: " << m_moreTargets << std::endl;

    
    

    isValidCandidate = false; //initially

    std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> > combVector = 
      combinations<xAOD::TrackParticleContainer::const_iterator> (itkTrk_itr_candidates, m_nTarget);

    std::vector<xAOD::TrackParticleContainer::const_iterator> best_recoTrk_itr;

    //    double minResidual = largeNumber;
    for (std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator combItr = combVector.begin();
	 combItr != combVector.end(); combItr++) {
      // TLorentzVector combMomentum(0.0,0.0,0.0,0.0);
      double totalCharge = 0.0;

      for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator combElement = (*combItr).begin();
	   combElement != (*combItr).end(); combElement++) {
	/*
	double px = (**combElement)->pt()*TMath::Cos((**combElement)->phi());
	double py = (**combElement)->pt()*TMath::Sin((**combElement)->phi());
	double pz = (**combElement)->pt()*TMath::SinH((**combElement)->eta());
	double e = TMath::Sqrt(TMath::Power(px,2) + TMath::Power(py,2) + TMath::Power(pz,2) + TMath::Power(m_decayMass,2));
	*/

	totalCharge += (**combElement)->charge();

	/*	TLorentzVector tmp(px,py,pz,e);
	combMomentum += tmp;
	*/

      }

      if (TMath::Abs(totalCharge) < 1E-6) {
	

	//only consider combinations that passes ranked ptCut

	/////////////////////

	HiggsEventSelector selector(*combItr);
	selector.SetRankedPtCut(m_truthPtCut);
	
	if (selector.Check() != 0) {
	  combVector.erase(combItr);
	  combItr--; //to account for the shift emerging from the erase
	  continue;
	}
	////////////////////

	


	/*
	double candidateParentMass = combMomentum.M();
	//Take the combination with mass closest to the Higgs --> BIAS IN ZZ->4mu
	double residual = TMath::Abs(candidateParentMass - m_parentMass);
	//apply here cuts on candidate set
	if (residual < minResidual && true
	    /*	    residual < m_parentMassTolerance / *) {
	  minResidual = residual;
	  best_recoTrk_itr = (*combItr);
	}
      */

      }

      else {
	combVector.erase(combItr);
	combItr--; //to account for the shift emerging from the erase
      }

    }

    //now we have only neutral combinations that passes the ranked pt cut
    // Now look for the combination with pair closest to the Z, NB the combinations have been ordered during Check() (1-2 on shell, 3-4 off shell)

    double minOnShellResidual = largeNumber;
    std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator bestOnShellCombination;

    for (std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator combItr = combVector.begin();
	 combItr != combVector.end(); combItr++) {
      std::vector<xAOD::TrackParticleContainer::const_iterator> onShellReco;
      onShellReco.push_back((*combItr)[0]);
      onShellReco.push_back((*combItr)[1]);

      double residual = TMath::Abs(InvariantMass::TrackInvariantMass(onShellReco, m_decayMass) - m_intermediateMass);
     
      if (residual < minOnShellResidual) {
	minOnShellResidual = residual;
	bestOnShellCombination = combItr;
      }
    }

    double maxMass = -1;

    //now, among the combinations that have the best on shell combination, choose the one with the highest 
    // 4mu mass
    for (std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator combItr = combVector.begin();
	 combItr != combVector.end(); combItr++) {

      //contains the best combination?
      if (std::find((*combItr).begin(),(*combItr).end(),(*bestOnShellCombination)[0]) == (*combItr).end() )
	continue;
      if (std::find((*combItr).begin(),(*combItr).end(),(*bestOnShellCombination)[1]) == (*combItr).end() )
	continue;

      //it contains the best combination
      double combMass = InvariantMass::TrackInvariantMass(*combItr, m_decayMass);
      if (combMass > maxMass) {
	maxMass = combMass;
	best_recoTrk_itr = *combItr;
      }

    }
    
    
    
    if (combVector.size() != 0) {
      isValidCandidate = true;
      itkTrk_itr_candidates = best_recoTrk_itr;
    }
    
  }
  


  
 

  ////////////////////////////////////////////////////////////////////
  // now look into candidates for parent decay and zero charge
  //  int nCount = 0;
  double maxEta = -1;  
  //cut on intermediate particles
  if (isValidCandidate) {  //!!!!!!!!!!!
    HiggsEventSelector selector(itkTrk_itr_candidates);
    int result;
    selector.SetRankedPtCut(m_truthPtCut);
    result = selector.Check();
    if (result != 0 && eventReco.recoCutCode == 0)
      eventReco.recoCutCode = result + 3;
    //selector.SetMassCut(m_massCutLow, m_massCutHigh);


    selector.SetDRCut(0.1);
    result = selector.Check();
    if (result != 0 && eventReco.recoCutCode == 0)
      eventReco.recoCutCode = 9;

    selector.SetCenteredOffShellMuonCut(0.0, 2.7);
    result = selector.Check();
    if (result != 0 && eventReco.recoCutCode == 0)
      eventReco.recoCutCode = 10;


    selector.SetOnShellMassCut(m_onShellMassCutLow, m_onShellMassCutHigh);
    result = selector.Check();
    if (result != 0 && eventReco.recoCutCode == 0)
      eventReco.recoCutCode = 11;


    selector.SetOffShellMassCut(m_offShellMassCutLow, m_offShellMassCutHigh);
    result = selector.Check();
    if (result != 0 && eventReco.recoCutCode == 0)
      eventReco.recoCutCode = 12;

    //selector.SetCenteredMuonCut(0.1, 2.7, 4);

    //selector.SetIsolationCut(recoIsolation, 100);
    
    if (maxEta < 2.7)
      eventReco.recoCutCode27 = eventReco.recoCutCode;
    else if (maxEta < 3.2)
      eventReco.recoCutCode32 = eventReco.recoCutCode;
    else if (maxEta < 4)
      eventReco.recoCutCode4  = eventReco.recoCutCode; 
    

    for (auto itr = itkTrk_itr_candidates.begin();
	 itr != itkTrk_itr_candidates.end();
	 itr++) {
      if (TMath::Abs((**itr)->eta()) > maxEta)
	maxEta = TMath::Abs((**itr)->eta());
    }

    if (eventReco.recoCutCode == 0)
      isValidCandidate = true;
    else
      isValidCandidate = false;

    /*
    std::vector<xAOD::TrackParticleContainer::const_iterator> onShell_candidates = PhysicsEvent::GetOnShellReco(itkTrk_itr_candidates, m_intermediateMass, m_decayMass);
    isValidCandidate = isValidCandidate && (TMath::Abs(InvariantMass::TrackInvariantMass(onShell_candidates, m_decayMass) - m_intermediateMass) < m_intermediateMassTolerance); 
    */
  }
    

  if (isValidCandidate /*&& isValidMatching*/) {
    eventReco.SetReco(itkTrk_itr_candidates, recoIsolation);

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

  if (isValidMatching) { //only for thesis plots
    eventFeatures.nRecoTracks = numRecoTracks;
    eventFeatures.nRecoTracks1GeV = numRecoTracks1GeV;
  }
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


  if (isValidMatching) { //only to make the thesis plots

  eventFeatures.nTruth = numTruth;
  eventFeatures.nChargedTruth = numChargedTruth;
  eventFeatures.nPrimaryChargedTruth = numPrimaryChargedTruth;
  eventFeatures.nPrimaryChargedTruth1GeV = numPrimaryChargedTruth1GeV;
  eventFeatures.nNeutralTruth = numNeutralTruth;

  eventFeatures.nTruth1GeV = numTruth1GeV;
  eventFeatures.nChargedTruth1GeV = numChargedTruth1GeV;
  eventFeatures.nNeutralTruth1GeV = numNeutralTruth1GeV;

  }

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
    else {
      truHist_pileup->FillHists((*truItr),1.0);
      // if ((*truItr)->prodVtx()) {
      // 	std::cout << "Pileup d0: " << (*truItr)->prodVtx()->perp() << std::endl;
      // }
    }
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

    if (isValidMatching) {
      if (itkTrk_itr != (*(truth_trackMatched.begin())).second &&
	  itkTrk_itr != (*(truth_trackMatched.begin()+1)).second &&
	  itkTrk_itr != (*(truth_trackMatched.begin()+2)).second &&
	  itkTrk_itr != (*(truth_trackMatched.begin()+3)).second) {
	trkHist_pileup -> FillHists((*itkTrk_itr),1.0);
      }
    }

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


EL::StatusCode HiggsAnalysis::postExecute() {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsAnalysis::finalize() {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HiggsAnalysis::histFinalize() {
  runHist_reco->FillHists(1.0);
  return EL::StatusCode::SUCCESS;
}

