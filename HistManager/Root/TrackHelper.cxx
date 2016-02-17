#include "HistManager/TrackHelper.h"
#include <TMath.h>

// utility function
const double deltaR(const double& eta1, const double& phi1, const double& eta2, const double& phi2) {
  return TMath::Sqrt(TMath::Power(eta2-eta1,2) + TMath::Power(phi2-phi1,2));
}

  // question : should be warning not info?  does "warning" exist?
int xAOD::TrackHelper::numberOfSiHits(const xAOD::TrackParticle *trkPart) {
  if(!trkPart) { return 0; }

  uint8_t tmp(0);
  int nSiHits(0);

  if( trkPart->summaryValue(tmp, xAOD::numberOfPixelHits) ) { nSiHits += tmp; }
  else { 
    Info("xAOD::TrackHelper::numberOfSiHits", "numberOfPixelHits is not stored in TrackParticle object"); 
    return -1;
  }
  if( trkPart->summaryValue(tmp, xAOD::numberOfSCTHits   ) ) { nSiHits += tmp; }
  else { 
    Info("xAOD::TrackHelper::numberOfSiHits", "numberOfPixelHits is not stored in TrackParticle object"); 
    return -1;
  }

  return nSiHits;
} // numberOfSiHits

int xAOD::TrackHelper::numberOfSiDeadSensors(const xAOD::TrackParticle *trkPart) {
  if(!trkPart) { return 0; }

  uint8_t tmp(0);
  int nSiDeadSensors(0);

  if( trkPart->summaryValue(tmp, xAOD::numberOfPixelDeadSensors) ) { nSiDeadSensors += tmp; }
  else { 
    std::cout << "xAOD::TrackHelper::numberOfSiDeadSensors: numberOfPixelDeadSensors is not stored in TrackParticle object" << std::endl;
    return -1;
  }
  if( trkPart->summaryValue(tmp, xAOD::numberOfSCTDeadSensors) ) { nSiDeadSensors += tmp; }
  else { 
    Info("xAOD::TrackHelper::numberOfSiDeadSensors", "numberOfSCTDeadSensors is not stored in TrackParticle object"); 
    return -1;
  }

  return nSiDeadSensors;
} // numberOfSiDeadSensors

float xAOD::TrackHelper::chiSqPerDoF(const xAOD::TrackParticle *trkPart) {
  float value( trkPart->numberDoF() );
  if( value == 0 ) {
    Info("xAOD::TrackHelper::chi2PerDoF", "Degress of Freedom is 0");
    return -1;
  }
  value = ( trkPart->chiSquared()/value );
  return value;
}

float xAOD::TrackHelper::truthMatchProb(const xAOD::TrackParticle *trkPart) {
  if( ! trkPart->isAvailable< float >( "truthMatchProbability" ) ) { return -999.; }
  return trkPart->auxdata< float >( "truthMatchProbability");
}

const xAOD::TruthParticle* xAOD::TrackHelper::truthParticle(const xAOD::TrackParticle *trkPart) {
  typedef ElementLink< xAOD::TruthParticleContainer > Link_t;
  static const char* NAME = "truthParticleLink";
  if( ! trkPart->isAvailable< Link_t >( NAME ) ) {
    return 0;
  }
  const Link_t& link = trkPart->auxdata< Link_t >( NAME );
  if( ! link.isValid() ) {
    return 0;
  }

  return *link;
} // truthParticle

const std::vector<ElementLink< xAOD::TrackStateValidationContainer > > xAOD::TrackHelper::measureOnTrack(const xAOD::TrackParticle *trkPart) {
  static const char* NAME("IDDET1_msosLink");
  typedef std::vector< ElementLink< xAOD::TrackStateValidationContainer > > vecLink_t;
  if( ! trkPart->isAvailable< vecLink_t >( NAME ) ) {
    const vecLink_t vecLink;
    return vecLink;
  }
  const vecLink_t& vecLink = trkPart->auxdata< vecLink_t >( NAME );
  return vecLink;
} // measureOnTrack

void xAOD::TrackHelper::countMissingAndWrongHitsOnTrack(const xAOD::TrackParticle *trk, const xAOD::TrackParticle *refTrk ) {

  // firstly, make sure that we are checking a track against its matching pseudotrack
  int trkBarcode = trk->auxdata< int >( "barcode" ); 
  //std::cout << "\t track barcode "  << trkBarcode << std::endl;  
  if( trkBarcode != refTrk->auxdata< int >( "barcode" ) ) 
  {
     std::cout << "xAOD::TrackHelper::countMissingAndWrongHitsOnTrack"          << std::endl;
     std::cout << "Error - not checking against right pseudoTrack. Returning. " << std::endl; 
     return;
  }

  static SG::AuxElement::Decorator< int > nCorrectPixel("nCorrectPixel"); 
  static SG::AuxElement::Decorator< int > nIncorrectPixel("nIncorrectPixel");   
  static SG::AuxElement::Decorator< float > fracCorrectPixel("fracCorrectPixel"); 
  static SG::AuxElement::Decorator< float > fracIncorrectPixel("fracIncorrectPixel");  
  static SG::AuxElement::Decorator< int > nMissingPixel("nMissingPixel");   
  
  static SG::AuxElement::Decorator< int > nCorrectSCT("nCorrectSCT"); 
  static SG::AuxElement::Decorator< int > nIncorrectSCT("nIncorrectSCT");   
  static SG::AuxElement::Decorator< float > fracCorrectSCT("fracCorrectSCT"); 
  static SG::AuxElement::Decorator< float > fracIncorrectSCT("fracIncorrectSCT");  
  static SG::AuxElement::Decorator< int > nMissingSCT("nMissingSCT");     

  uint8_t getInt(0); 

  // initialise dressing
  nCorrectPixel(*trk)      = 0;
  nIncorrectPixel(*trk)    = 0;
  fracCorrectPixel(*trk)   = 0.0;
  fracIncorrectPixel(*trk) = 0.0;
  refTrk->summaryValue( getInt, xAOD::numberOfPixelHits );
  nMissingPixel(*trk)      =  getInt;
    
  nCorrectSCT(*trk)      = 0;
  nIncorrectSCT(*trk)    = 0;
  fracCorrectSCT(*trk)   = 0.0;
  fracIncorrectSCT(*trk) = 0.0;
  refTrk->summaryValue( getInt, xAOD::numberOfSCTHits );   
  nMissingSCT(*trk)      =  getInt;
  
  // ****************************************************************************************************************************
  //  check the number of incorrectly assigned hits on track: 
  //  1) retrieve ElementLinks to MSOSs on track
  //  2) loop on the trk MSOSs and get link to corresponding PRD (no cluster barcode info on MSOS, but yes on PRD...)
  //  3) check PRD barcode against trk barcode
  
  typedef std::vector< ElementLink< xAOD::TrackStateValidationContainer > > vecLink_t;
  vecLink_t  linksToMSOS = measureOnTrack( trk );

  if(linksToMSOS.size() == 0){ 
     std::cout << "\t Error - not retrieved vector of ElementLink to pixel MSOS of this track. Return..."  << std::endl; 
     return;  
  }

  uint8_t nPix(0), nCorrectPix(0), nIncorrectPix(0); 
  uint8_t nSct(0), nCorrectSct(0), nIncorrectSct(0);   
  for(auto msos_it = linksToMSOS.begin(); msos_it != linksToMSOS.end(); ++msos_it)
  {
    //Check if the element link is valid
    if( ! (*msos_it).isValid() ) { continue; }
    const xAOD::TrackStateValidation* msos = *(*msos_it); 
    if( msos->type() ==  6 ) { continue; }  // must not be a hole
    if( !(msos->trackMeasurementValidationLink().isValid() && *(msos->trackMeasurementValidationLink())) )
    {
       std::cout << "\t Invalid link to PRD for this msos - Skipping to next msos" << std::endl;    
    }

    // ********* This seems ineffective to prevent nMissingHits from being < 0 in few cases ****************
    // check if already considered this MSOS's PRD's sibling (if any)
      // siblings are copies of the parent with a new local position
      // the identifier and positions can be different but the layer, size, and charge information will be the same
    /*
    bool skip(false);
    for( auto  msos_newit = linksToMSOS.begin(); msos_newit != msos_it; ++msos_newit) {
       const  xAOD::MeasurementStateOnSurface *this_msos = *(*msos_newit);
       if ( xAOD::TrackHelper::areSiblings( msos, this_msos) ) { skip = true; break; }
    }
    if(skip) { continue; }
    */
    //**************************************************************************************************
    static SG::AuxElement::ConstAccessor< std::vector<int> > acc_sihit_barcode("sihit_barcode");
    
    //Get pixel clusters
    if( msos->detType() == 1){ 
      const xAOD::TrackMeasurementValidation* pixelCluster =  *(msos->trackMeasurementValidationLink());        
      // count pixel hits on track 
      nPix++;  
      const std::vector<int>   &  sihit_barcode   = acc_sihit_barcode(*pixelCluster);
      bool isMeasurementCorrect(false); 
      for( unsigned int i(0); i < sihit_barcode.size(); ++i ){
        //std::cout << "\t\t sihit_barcode[" << i << "] = " << sihit_barcode[i] << std::endl;
        if( sihit_barcode[i] == trkBarcode ){
          isMeasurementCorrect = true; break;
        }
      }
      if( isMeasurementCorrect ){ nCorrectPix++;} 
      else { nIncorrectPix++;}
    }
    //Get SCT clusters
    else if( msos->detType() == 2){   
      const xAOD::TrackMeasurementValidation* sctCluster =  *(msos->trackMeasurementValidationLink());        
      // count SCT hits on track 
      nSct++;  
      const std::vector<int>   &  sihit_barcode   = acc_sihit_barcode(*sctCluster);
      bool isMeasurementCorrect(false); 
      for( unsigned int i(0); i < sihit_barcode.size(); ++i ){
        //std::cout << "\t\t sihit_barcode[" << i << "] = " << sihit_barcode[i] << std::endl;
        if( sihit_barcode[i] == trkBarcode ){
          isMeasurementCorrect = true; break;
        }
      }
      if( isMeasurementCorrect ){ nCorrectSct++;} 
      else { nIncorrectSct++;}
    }
  
  }
  
  //dress the TrackParticle
  nCorrectPixel(*trk)      = nCorrectPix;
  nIncorrectPixel(*trk)    = nIncorrectPix;
  fracCorrectPixel(*trk)   = static_cast<float>(nCorrectPix)/nPix;
  fracIncorrectPixel(*trk) = static_cast<float>(nIncorrectPix)/nPix;
  // count missing hits on track 
  refTrk->summaryValue( getInt, xAOD::numberOfPixelHits );
  nMissingPixel(*trk) =  getInt - nCorrectPix;  

  nCorrectSCT(*trk)      = nCorrectSct;
  nIncorrectSCT(*trk)    = nIncorrectSct;
  fracCorrectSCT(*trk)   = static_cast<float>(nCorrectSct)/nSct;
  fracIncorrectSCT(*trk) = static_cast<float>(nIncorrectSct)/nSct;
  // count missing hits on track 
  refTrk->summaryValue( getInt, xAOD::numberOfSCTHits );
  nMissingSCT(*trk) =  getInt - nCorrectSct;  

/*
uint8_t getInt2(0);
trk->summaryValue( getInt2, xAOD::numberOfPixelHits ); 
std::cout << "************" <<  std::endl; 
std::cout << "  " << std::endl;
std::cout << "this track" << std::endl;
std::cout << "\t nPixHits  " <<  static_cast<int>(getInt2)     << std::endl;
std::cout << "\t nCorrectPixel  "     << static_cast<int>(nCorrectPix)   << std::endl;
std::cout << "\t nIncorrectPixel  "   << static_cast<int>(nIncorrectPix) << std::endl;
std::cout << "\t nMissingPixel "      << getInt - nCorrectPix << std::endl;
//std::cout << "\t nMissingPixel "      << nMissingPix  << std::endl;
std::cout << "\t fracCorrectPixel "   << static_cast<float>(nCorrectPix)/nPix   << std::endl;
std::cout << "\t fracIncorrectPixel " << static_cast<float>(nIncorrectPix)/nPix << std::endl;
std::cout << "  ----- " << std::endl;
std::cout << "reference track  " << std::endl;
std::cout << "\t nPixHits  " << static_cast<int>(getInt)  << std::endl;
std::cout << "************" <<  std::endl;   
*/
}

bool xAOD::TrackHelper::areSiblings( const xAOD::TrackStateValidation* iMSOS, const xAOD::TrackStateValidation* jMSOS ) {
  if( !iMSOS->trackMeasurementValidationLink().isValid() ) { return false; }
  if( !jMSOS->trackMeasurementValidationLink().isValid() ) { return false; }
  const xAOD::TrackMeasurementValidation* iclus = *(iMSOS->trackMeasurementValidationLink());
  const xAOD::TrackMeasurementValidation* jclus = *(jMSOS->trackMeasurementValidationLink());
  if( !iclus || !jclus) { return false; }
  return areSiblings( iclus, jclus );
} // areSiblings

bool xAOD::TrackHelper::areSiblings( const xAOD::TrackMeasurementValidation* iPRD, const xAOD::TrackMeasurementValidation* jPRD) { 
  // only pixel clusters are split
  if( iPRD->isAvailable< int >( "side" ) ) { return false; } // only SCT clusters have a side
  // siblings are copies of the parent with a new local position
  // the identifier and positions can be different but the layer, size, and charge information will be the same
  if( !iPRD->isAvailable< int >( "layer" ) ) { return false; }
  if( !jPRD->isAvailable< int >( "layer" ) ) { return false; }
  if( iPRD->auxdata< int >( "layer" ) != jPRD->auxdata< int >( "layer" ) ) { return false; }
  if( !iPRD->isAvailable< int >( "sizePhi" ) ) { return false; }
  if( !jPRD->isAvailable< int >( "sizePhi" ) ) { return false; }
  if( iPRD->auxdata< int >( "sizePhi" ) != jPRD->auxdata< int >( "sizePhi" ) ) { return false; }
  if( !iPRD->isAvailable< int >( "sizeZ" ) ) { return false; }
  if( !jPRD->isAvailable< int >( "sizeZ" ) ) { return false; }
  if( iPRD->auxdata< int >( "sizeZ" ) != jPRD->auxdata< int >( "sizeZ" ) ) { return false; }
  if( !iPRD->isAvailable< float >( "charge" ) ) { return false; }
  if( !jPRD->isAvailable< float >( "charge" ) ) { return false; }
  if( iPRD->auxdata< float >( "charge" ) != jPRD->auxdata< float >( "charge" ) ) { return false; }
  return true;
} // areSiblings

std::multimap< xAOD::TrackMeasurementValidation*, xAOD::TrackMeasurementValidation* > xAOD::TrackHelper::getSiblingMap( const xAOD::TrackMeasurementValidationContainer* clusters ) {
  std::multimap< xAOD::TrackMeasurementValidation*, xAOD::TrackMeasurementValidation* > family;
  xAOD::TrackMeasurementValidationContainer::const_iterator clus_itr = clusters->begin();
  xAOD::TrackMeasurementValidationContainer::const_iterator clus_jtr = clusters->begin();
  xAOD::TrackMeasurementValidationContainer::const_iterator clus_end = clusters->end();
  for( ; clus_itr != clus_end; ++clus_itr ) {
    clus_jtr = clus_itr;
    clus_jtr++;
    for( ; clus_jtr != clus_end; ++clus_itr ) {
      if( areSiblings( *clus_itr, *clus_jtr ) ) {
        family.insert( std::pair<xAOD::TrackMeasurementValidation*,xAOD::TrackMeasurementValidation*>(*clus_itr,*clus_jtr) );
        family.insert( std::pair<xAOD::TrackMeasurementValidation*,xAOD::TrackMeasurementValidation*>(*clus_jtr,*clus_itr) );
      } // if are siblings
    } // loop over rest of clusters
  } // full loop over clusters
  return family;
} // getSiblingMap

bool xAOD::TrackHelper :: isPrimary( const xAOD::TruthParticle* tp )
{
  if(!tp) 
    return false;
    
  if( tp->status() != 1 )     //Check to see if its a stable particle
    return false;
  
  //Note this const cast will not  be required in the future -- there is bug in the truth edm
  if( const_cast<xAOD::TruthParticle*>(tp)->isNeutral() )         //Clearly for tracking we don't care about Neutrals 
    return false;

  if( tp->barcode() == 0 ||   //Barcode of zero indicates there was no truth paticle found for this track 
      tp->barcode() >= 200e3) 
    return false;
 
  if( tp->pt() < 500.   ||  fabs(tp->eta()) > 2.5 ) // Check the particle is within acceptance 
    return false;
            
  // make sure particle decays before the last pixel layer
  if( !tp->hasProdVtx() || tp->prodVtx()->perp() > 110) 
    return false;

  return true;
}


bool xAOD::TrackHelper :: isSecondary( const xAOD::TruthParticle* tp )
{
  if(!tp) 
    return false;
  
  if( tp->status() != 1 )     //Check to see if its a stable particle
    return false;
  //Note this const cast will not  be required in the future -- there is bug in the truth edm
  if( const_cast<xAOD::TruthParticle*>(tp)->isNeutral() )         //Clearly for tracking we don't care about Neutrals 
    return false;

  if( tp->barcode() < 200e3 )    //larger than this indicates secondary  particles
    return false;
 
  if( tp->pt() < 500.   || fabs(tp->eta()) >  2.5 ) // Check the particle is within acceptance 
    return false;

  // make sure particle decays before the last SCT layer
  if( !tp->hasProdVtx() || tp->prodVtx()->perp() > 500)
    return false;

  return true;
}

bool xAOD::TrackHelper :: isFake_TruthMatchProb( const xAOD::TrackParticle* tp )
{

  static const char* NAME = "truthMatchProbability";
  if( ! tp->isAvailable< float >( NAME ) ) {
    Warning("isFake(..)", "Track Particle has no MatchProb is this data?" );
    return true;
  }

  const float truthProb = tp->auxdata< float >( NAME );
  if( truthProb < 0.7 ) {
    return true;
  }
  const xAOD::TruthParticle* truth = xAOD::TrackHelper::truthParticle( tp );

  if (truth && truth->barcode()!=0)
      return false;
  return true;    
} 

bool xAOD::TrackHelper :: isFake( const xAOD::TrackParticle* tp )
{
  const xAOD::TruthParticle* truth = xAOD::TrackHelper::truthParticle( tp );

  if (truth != NULL)
    {
      const double dR = deltaR(truth->eta(), truth->phi(), tp->eta(), tp->phi());
      if (dR < 0.02 && truth->barcode() != 0) return false;
      else return true;
  }
  else //truthParticle not retrieved
    {
      std::cerr << "In xAOD::TrackHelper::isFake(const xAOD::TrackParticle*) - Null TruthParticle pointer" << std::endl;
      return true;
  }
  
  std::cerr << "In xAOD::TrackHelper::isFake(const xAOD::TrackParticle*) - Should have never showed this message" << std::endl;
  return true;    
}
 
bool xAOD::TrackHelper :: isPrimary( const xAOD::TrackParticle* tp )
{
  if(!tp) 
    return false;
    
  if( isFake(tp) )
    return false;
  
  const xAOD::TruthParticle* truth =  xAOD::TrackHelper::truthParticle( tp );
  if( !truth  ){
    return false;
  }

  return isPrimary( truth );
}


bool xAOD::TrackHelper :: isSecondary( const xAOD::TrackParticle* tp)
{
  if(!tp) 
    return false;

  if( isFake(tp) )
    return false;

  
  const xAOD::TruthParticle* truth =  xAOD::TrackHelper::truthParticle( tp );
  if( !truth  ){
    return false;
  }

  return isSecondary( truth );
}

