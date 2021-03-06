#include <HistManager/ClusterHists.h>


ClusterHists::ClusterHists(TString name) {
   m_name = "ClusterHist_" + name;
   m_label = "Cluster "; // don't forget the space
}

ClusterHists::~ClusterHists() {}

void ClusterHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  /* layout */
  m_SiliconLayout_full_xy    = declare2D(m_name,"SiliconLayout_full_xy","x[mm]","y[mm]",200,-1100,1100,200,-1100,1100);
  m_SiliconLayout_full_zr    = declare2D(m_name,"SiliconLayout_full_zr","z[mm]","r[mm]",200,-3100,3100,200,   0,1100);
  m_SiliconLayout_quarter_xy = declare2D(m_name,"SiliconLayout_quarter_xy","x[mm]","y[mm]",200,0,1100,200,0,1100);
  m_SiliconLayout_quarter_zr = declare2D(m_name,"SiliconLayout_quarter_zr","z[mm]","r[mm]",200,0,3100,200,0,1100);

  m_PixelLayout_full_xy    = declare2D(m_name,"PixelLayout_full_xy","x[mm]","y[mm]",200,-350,350,200,-350,350);
  m_PixelLayout_full_zr    = declare2D(m_name,"PixelLayout_full_zr","z[mm]","r[mm]",200,-2000,2000,200,   0,350);
  m_PixelLayout_quarter_xy = declare2D(m_name,"PixelLayout_quarter_xy","x[mm]","y[mm]",200,0,350,200,0,350);
//  m_PixelLayout_quarter_zr = declare2D(m_name,"PixelLayout_quarter_zr","z[mm]","r[mm]",200,0,2000,200,0,350);
  m_PixelLayout_quarter_zr = declare2D(m_name,"PixelLayout_quarter_zr","z[mm]","r[mm]",800,0,200,200,30,50);

  /* positions */
  m_localX = declare1D(m_name, "localX", m_label + "local X [mm]",    120, -9., 9. );
  m_localY = declare1D(m_name, "localY", m_label + "local Y [mm]",    120, -31., 31. );
  m_localXError = declare1D(m_name, "localXError", m_label + "local X error [mm]",    60, 0., 5. );
  m_localYError = declare1D(m_name, "localYError", m_label + "local Y error [mm]",    60, 0., 20. );
  m_globalX = declare1D(m_name, "globalX", m_label + "global X [mm]",    150, -150., 150. );
  m_globalY = declare1D(m_name, "globalY", m_label + "global Y [mm]",    150, -150., 150. );
  m_globalZ = declare1D(m_name, "globalZ", m_label + "global Z [mm]",    100, -200., 200. );


  /* basic properties */
  m_layer   = declare1D(m_name, "layer",   m_label + "layer",     4,  -0.5, 3.5 );
  m_nPart   = declare1D(m_name, "nPart",   m_label + "N Particles",  10, -0.5, 9.5);
  m_nPriPart= declare1D(m_name, "nPriPart",m_label + "N Primary Particles",  10, -0.5, 9.5);
  m_size    = declare1D(m_name, "size",    m_label + "Size",      20, -0.5, 19.5);
  m_sizePhi = declare1D(m_name, "sizePhi", m_label + "Size #phi", 10, -0.5, 9.5);
  m_sizeZ   = declare1D(m_name, "sizeZ",   m_label + "Size z",    10, -0.5, 9.5);
  m_charge  = declare1D(m_name, "charge",  m_label + "charge [e]",120, 0., 240e3);
  m_ToT     = declare1D(m_name, "ToT",     m_label + "ToT [bc]",  120, 0., 240.);
  m_isFake  = declare1D(m_name, "isFake",  m_label + "Is Fake",   2, -0.5, 1.5);
  m_gangedPixel = declare1D(m_name, "gangedPixel", m_label + "Ganged Pixel",   2, -0.5, 1.5);
  m_isSplit     = declare1D(m_name, "isSplit",     m_label + "Is Split",   2, -0.5, 1.5);
  m_splitProb1  = declare1D(m_name, "splitProb1",  m_label + "Split Prob: 2 Trk Hypo",   100, -1.0, 1.0);
  m_splitProb2  = declare1D(m_name, "splitProb2",  m_label + "Split Prob: 3 Trk Hypo",   100, -1.0, 1.0);
  m_NN_sizeX    = declare1D(m_name, "sizeX",       m_label + "Size x", 10, -0.5, 9.5);
  m_NN_sizeY    = declare1D(m_name, "sizeY",       m_label + "Size y", 10, -0.5, 9.5);
  m_NN_phiBS    = declare1D(m_name, "NN_phiBS",    m_label + "NN phi BS", 120, 0.0, 0.40);
  m_NN_thetaBS  = declare1D(m_name, "NN_thetaBS",  m_label + "NN theta BS", 100, -1.0, 1.0);

  /* fun with 2D */
  m_size_VS_charge    = declare2D( m_name, "size_VS_charge",    m_label + "Size",      m_label + "charge [e]", 20, -0.5, 19.5, 60, 0., 240e3 );
  m_sizePhi_VS_charge = declare2D( m_name, "sizePhi_VS_charge", m_label + "Size #phi", m_label + "charge [e]", 10, -0.5, 9.5, 60, 0., 240e3 );
  m_sizeZ_VS_charge   = declare2D( m_name, "sizeZ_VS_charge",   m_label + "Size z",    m_label + "charge [e]", 10, -0.5, 9.5, 60, 0., 240e3 );
  m_size_VS_ToT       = declare2D( m_name, "size_VS_ToT",       m_label + "Size",      m_label + "ToT [bc]", 20, -0.5, 19.5, 60, 0., 240 );
  m_sizePhi_VS_ToT    = declare2D( m_name, "sizePhi_VS_ToT",    m_label + "Size #phi", m_label + "ToT [bc]", 10, -0.5, 9.5, 60, 0., 240 );
  m_sizeZ_VS_ToT      = declare2D( m_name, "sizeZ_VS_ToT",      m_label + "Size z",    m_label + "ToT [bc]", 10, -0.5, 9.5, 60, 0., 240 );

  // hits in cluster
  m_nHitCluster  = declare1D(m_name, "nHitCluster",  "Number of hits in cluster", 20, -0.5, 19.5);
  m_hitCharge    = declare1D(m_name, "hitCharge",    "Charge [e]",     200,    0, 40000);
  m_hitTot       = declare1D(m_name, "hitTot",       "ToT [bc]",        18, -0.5, 17.5);
  m_hitMxTot     = declare1D(m_name, "hitMxTot",     "max.ToT [bc]",    18, -0.5, 17.5);
  m_hitMnTot     = declare1D(m_name, "hitMnTot",     "min.ToT [bc]",    18, -0.5, 17.5);
  m_hitAveTot    = declare1D(m_name, "hitAveTot",    "ave.ToT [bc]",    18, -0.5, 17.5);
  m_hitMxCharge  = declare1D(m_name, "hitMxCharge",  "max.Charge [e]", 200,    0, 40000);
  m_hitMnCharge  = declare1D(m_name, "hitMnCharge",  "min.Charge [e]", 200,    0, 40000);
  m_hitAveCharge = declare1D(m_name, "hitAveCharge", "ave.Charge [e]", 200,    0, 40000);
  m_hit1Charge   = declare1D(m_name, "hit1Charge",   "Charge [e]",     200,    0, 40000);
  m_hit1Tot      = declare1D(m_name, "hit1Tot",      "ToT [bc]",        18, -0.5, 17.5);

} // BookHists

void ClusterHists::FillHists(const xAOD::TrackMeasurementValidation* prd, float weight) const {

  /* layout */
  double rrr = sqrt(prd->globalX()*prd->globalX()+prd->globalY()*prd->globalY());

  if (rrr<380.0) {
    m_SiliconLayout_full_zr    -> Fill(prd->globalZ(),rrr);
    m_SiliconLayout_quarter_zr -> Fill(prd->globalZ(),rrr);
    m_PixelLayout_full_zr      -> Fill(prd->globalZ(),rrr);
    m_PixelLayout_quarter_zr   -> Fill(prd->globalZ(),rrr);
    if (fabs(prd->globalZ())<750.0) {
      m_SiliconLayout_full_xy    -> Fill(prd->globalX(),prd->globalY());
      m_SiliconLayout_quarter_xy -> Fill(prd->globalX(),prd->globalY());
      m_PixelLayout_full_xy      -> Fill(prd->globalX(),prd->globalY());
      m_PixelLayout_quarter_xy   -> Fill(prd->globalX(),prd->globalY());
    }
  }
  if (rrr>380.0) {
    m_SiliconLayout_full_zr    -> Fill(prd->globalZ(),rrr);
    m_SiliconLayout_quarter_zr -> Fill(prd->globalZ(),rrr);
    if (fabs(prd->globalZ())<1350.0) {
      m_SiliconLayout_full_xy    -> Fill(prd->globalX(),prd->globalY());
      m_SiliconLayout_quarter_xy -> Fill(prd->globalX(),prd->globalY());
    }
  }

  /* positions */
  m_localX -> Fill(prd->localX(),weight);
  m_localY -> Fill(prd->localY(),weight);
  m_localXError -> Fill(prd->localXError(),weight);
  m_localYError -> Fill(prd->localYError(),weight);
  m_globalX -> Fill(prd->globalX(),weight);
  m_globalY -> Fill(prd->globalY(),weight);
  m_globalZ -> Fill(prd->globalZ(),weight);

  /* basic properties */
  m_layer -> Fill(prd->auxdata<int>("layer"),weight);
  if (prd->isAvailable<int>("nParticles")) {
    m_nPart    -> Fill(prd->auxdata<int>("nParticles"),weight);
    m_nPriPart -> Fill(prd->auxdata<int>("nPrimaryParticles"),weight);
  }

  // only store for pixel
  if (rrr<380.0) {
    //m_size    -> Fill(prd->auxdata<int>("size"),weight);
    m_sizePhi -> Fill(prd->auxdata<int>("sizePhi"),weight);
    m_sizeZ   -> Fill(prd->auxdata<int>("sizeZ"),weight);
    m_charge  -> Fill(prd->auxdata<float>("charge"),weight);
    m_ToT     -> Fill(prd->auxdata<int>("ToT"),weight);
    m_isFake  -> Fill(prd->auxdata<char>("isFake"),weight);
    m_gangedPixel -> Fill(prd->auxdata<char>("gangedPixel"),weight);
    m_isSplit     -> Fill(prd->auxdata<int>("isSplit"),weight);
    m_splitProb1  -> Fill(prd->auxdata<float>("splitProbability1"),weight);
    m_splitProb2  -> Fill(prd->auxdata<float>("splitProbability2"),weight);
    //  m_NN_sizeX    -> Fill( prd->auxdata< int >( "NN_sizeX" ) );
    //  m_NN_sizeY    -> Fill( prd->auxdata< int >( "NN_sizeY" ) );
    //  m_NN_phiBS    -> Fill( prd->auxdata< float >( "NN_phiBS" ) );
    //  m_NN_thetaBS  -> Fill( prd->auxdata< float >( "NN_thetaBS" ) );

    /* fun with 2D */
    //    m_size_VS_charge    -> Fill(prd->auxdata<int>("size"),    prd->auxdata<float>("charge"),weight);
    m_sizePhi_VS_charge -> Fill(prd->auxdata<int>("sizePhi"), prd->auxdata<float>("charge"),weight);
    m_sizeZ_VS_charge   -> Fill(prd->auxdata<int>("sizeZ"),   prd->auxdata<float>("charge"),weight);
    //    m_size_VS_ToT       -> Fill(prd->auxdata<int>("size"),    prd->auxdata<int>("ToT"),weight);
    m_sizePhi_VS_ToT    -> Fill(prd->auxdata<int>("sizePhi"), prd->auxdata<int>("ToT"),weight);
    m_sizeZ_VS_ToT      -> Fill(prd->auxdata<int>("sizeZ"),   prd->auxdata<int>("ToT"),weight);
  }

  std::vector<int>  totList;
  std::vector<float> chList;
  std::vector<float> CTerm,ATerm,ETerm;
  if (prd->isAvailable<std::vector<int>>("rdo_tot"))      { totList = prd->auxdata<std::vector<int>>("rdo_tot"); }
  if (prd->isAvailable<std::vector<float>>("rdo_charge")) { chList  = prd->auxdata<std::vector<float>>("rdo_charge"); }
  if (prd->isAvailable<std::vector<float>>("rdo_Cterm"))  { CTerm = prd->auxdata<std::vector<float>>("rdo_Cterm"); }
  if (prd->isAvailable<std::vector<float>>("rdo_Aterm"))  { ATerm = prd->auxdata<std::vector<float>>("rdo_Aterm"); }
  if (prd->isAvailable<std::vector<float>>("rdo_Eterm"))  { ETerm = prd->auxdata<std::vector<float>>("rdo_Eterm"); }
 
  // Occupancy
  if (prd->isAvailable<std::vector<int>>("rdo_phi_pixel_index")) {
    int mxTot = 0;
    int mnTot = 10000;
    int avTot = 0;
    float mxCharge = 0;
    float mnCharge = 100000000;
    float avCharge = 0;

    int nHits = prd->auxdata<std::vector<int>>("rdo_phi_pixel_index").size();
    m_nHitCluster -> Fill(1.0*nHits,weight);
    int eta_module = prd->auxdataConst< int >("eta_module"); // -20,20
    int phi_module = prd->auxdataConst< int >("phi_module"); // 0-50

    for (int i=0; i<nHits; i++) {

      int eta = prd->auxdata<std::vector<int>>("rdo_eta_pixel_index")[i];
      int phi = prd->auxdata<std::vector<int>>("rdo_phi_pixel_index")[i];

      int tot = 0;
      if (totList.size()==nHits) { tot = totList[i]; }

      float charge = 0.0;
      if (chList.size()==nHits) { charge = chList[i]; }

      // re-calculate tot
      if (tot==0) {
        if (CTerm.size()==nHits) {
          if (charge+CTerm[i]!=0.0) { tot = ATerm[i]*(charge+ETerm[i])/(charge+CTerm[i]); }
        }
      }

      // charge calibration
      float newCharge = 0.0;
      if (ATerm.size()==nHits) {
        if (ATerm[i]>0.0 && tot/ATerm[i]<1.0) {
          newCharge = (CTerm[i]*tot/ATerm[i]-ETerm[i])/(1.0-tot/ATerm[i]);
        } 
        else {
          newCharge=0.0;
        }
      }

      // compare nominal and new charge
//      if (TMath::Abs(newCharge-charge)>0.1) {
//        std::cout << "STSTST " << charge << " " << newCharge << " " << tot << " " << totList.size() << " " << ATerm[i] << " " << CTerm[i] << " " << ETerm[i] << std::endl;
//      }


      m_hitCharge    -> Fill(charge,weight);
      m_hitTot       -> Fill(1.0*tot,weight);
      avTot += tot;
      if (tot>mxTot) { mxTot=tot; }
      if (tot<mnTot) { mnTot=tot; }

      avCharge += charge;
      if (charge>mxCharge) { mxCharge=charge; }
      if (charge<mnCharge) { mnCharge=charge; }
    }

    m_hitMxTot  -> Fill(1.0*mxTot,weight);
    m_hitMnTot  -> Fill(1.0*mnTot,weight);
    m_hitAveTot -> Fill(1.0*avTot/prd->auxdata<std::vector<int>>("rdo_phi_pixel_index").size(),weight);
    m_hitMxCharge  -> Fill(mxCharge,weight);
    m_hitMnCharge  -> Fill(mnCharge,weight);
    m_hitAveCharge -> Fill(avCharge/prd->auxdata<std::vector<int>>("rdo_phi_pixel_index").size(),weight);

    if (nHits==1) {
      int tot = 0;
      if (totList.size()==nHits) { tot = totList[0]; }

      float charge = 0.0;
      if (chList.size()==nHits) { charge = chList[0]; }

      // re-calculate tot
      if (tot==0) {
        if (CTerm.size()==nHits) {
          if (charge+CTerm[0]!=0.0) { tot = ATerm[0]*(charge+ETerm[0])/(charge+CTerm[0]); }
        }
      }

      // charge calibration
      float newCharge = 0.0;
      if (ATerm.size()==nHits) {
        if (ATerm[0]>0.0 && tot/ATerm[0]<1.0) {
          newCharge = (CTerm[0]*tot/ATerm[0]-ETerm[0])/(1.0-tot/ATerm[0]);
        } 
        else {
          newCharge=0.0;
        }
      }

      m_hit1Charge -> Fill(charge,weight);
      m_hit1Tot    -> Fill(1.0*tot,weight);
    }

  }


} // FillHists

