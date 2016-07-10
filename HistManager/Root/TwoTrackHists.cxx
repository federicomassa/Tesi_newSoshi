#include <HistManager/TwoTrackHists.h>
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include <Utility/TrackHelper.h>

#include "xAODTracking/TrackParticleContainer.h"


TwoTrackHists::TwoTrackHists(TString name) {
   m_name = "TwoTrackHist_" + name;
}

TwoTrackHists::~TwoTrackHists() {}

void TwoTrackHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  m_TargetLayerOfInterest        = declare1D(m_name, "TargetLayerOfInterest", "target layer of interest", 6,-0.5,5.5);
  m_ClusterPatternIndex          = declare1D(m_name, "ClusterPatternIndex", "cluster pattern index", 70,-0.5,69.5);

  m_minMSOSLocalDistanceAtLayer  = declare1D(m_name, "minMSOSLocalDistanceAtLayer", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX = declare1D(m_name, "minMSOSLocalDistanceAtLayerX", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY = declare1D(m_name, "minMSOSLocalDistanceAtLayerY", "minimum local distance of MSOS [mm]",100,0,10.0);

  m_minMSOSLocalDistanceAtLayer0   = declare1D(m_name, "minMSOSLocalDistanceAtLayer0", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer1   = declare1D(m_name, "minMSOSLocalDistanceAtLayer1", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer2   = declare1D(m_name, "minMSOSLocalDistanceAtLayer2", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer3   = declare1D(m_name, "minMSOSLocalDistanceAtLayer3", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer4   = declare1D(m_name, "minMSOSLocalDistanceAtLayer4", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer5   = declare1D(m_name, "minMSOSLocalDistanceAtLayer5", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer6   = declare1D(m_name, "minMSOSLocalDistanceAtLayer6", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer7   = declare1D(m_name, "minMSOSLocalDistanceAtLayer7", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer8   = declare1D(m_name, "minMSOSLocalDistanceAtLayer8", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer9   = declare1D(m_name, "minMSOSLocalDistanceAtLayer9", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer10  = declare1D(m_name, "minMSOSLocalDistanceAtLayer10", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer11  = declare1D(m_name, "minMSOSLocalDistanceAtLayer11", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer12  = declare1D(m_name, "minMSOSLocalDistanceAtLayer12", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer13  = declare1D(m_name, "minMSOSLocalDistanceAtLayer13", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer14  = declare1D(m_name, "minMSOSLocalDistanceAtLayer14", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer15  = declare1D(m_name, "minMSOSLocalDistanceAtLayer15", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer16  = declare1D(m_name, "minMSOSLocalDistanceAtLayer16", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer17  = declare1D(m_name, "minMSOSLocalDistanceAtLayer17", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer18  = declare1D(m_name, "minMSOSLocalDistanceAtLayer18", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer19  = declare1D(m_name, "minMSOSLocalDistanceAtLayer19", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer20  = declare1D(m_name, "minMSOSLocalDistanceAtLayer20", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer21  = declare1D(m_name, "minMSOSLocalDistanceAtLayer21", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer22  = declare1D(m_name, "minMSOSLocalDistanceAtLayer22", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer23  = declare1D(m_name, "minMSOSLocalDistanceAtLayer23", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer24  = declare1D(m_name, "minMSOSLocalDistanceAtLayer24", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer25  = declare1D(m_name, "minMSOSLocalDistanceAtLayer25", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer26  = declare1D(m_name, "minMSOSLocalDistanceAtLayer26", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer27  = declare1D(m_name, "minMSOSLocalDistanceAtLayer27", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer28  = declare1D(m_name, "minMSOSLocalDistanceAtLayer28", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer29  = declare1D(m_name, "minMSOSLocalDistanceAtLayer29", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer30  = declare1D(m_name, "minMSOSLocalDistanceAtLayer30", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer31  = declare1D(m_name, "minMSOSLocalDistanceAtLayer31", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer32  = declare1D(m_name, "minMSOSLocalDistanceAtLayer32", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer33  = declare1D(m_name, "minMSOSLocalDistanceAtLayer33", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer34  = declare1D(m_name, "minMSOSLocalDistanceAtLayer34", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer35  = declare1D(m_name, "minMSOSLocalDistanceAtLayer35", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer36  = declare1D(m_name, "minMSOSLocalDistanceAtLayer36", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer37  = declare1D(m_name, "minMSOSLocalDistanceAtLayer37", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer38  = declare1D(m_name, "minMSOSLocalDistanceAtLayer38", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer39  = declare1D(m_name, "minMSOSLocalDistanceAtLayer39", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer40  = declare1D(m_name, "minMSOSLocalDistanceAtLayer40", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer41  = declare1D(m_name, "minMSOSLocalDistanceAtLayer41", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer42  = declare1D(m_name, "minMSOSLocalDistanceAtLayer42", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer43  = declare1D(m_name, "minMSOSLocalDistanceAtLayer43", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer44  = declare1D(m_name, "minMSOSLocalDistanceAtLayer44", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer45  = declare1D(m_name, "minMSOSLocalDistanceAtLayer45", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer46  = declare1D(m_name, "minMSOSLocalDistanceAtLayer46", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer47  = declare1D(m_name, "minMSOSLocalDistanceAtLayer47", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer48  = declare1D(m_name, "minMSOSLocalDistanceAtLayer48", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer49  = declare1D(m_name, "minMSOSLocalDistanceAtLayer49", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer50  = declare1D(m_name, "minMSOSLocalDistanceAtLayer50", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer51  = declare1D(m_name, "minMSOSLocalDistanceAtLayer51", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer52  = declare1D(m_name, "minMSOSLocalDistanceAtLayer52", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer53  = declare1D(m_name, "minMSOSLocalDistanceAtLayer53", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer54  = declare1D(m_name, "minMSOSLocalDistanceAtLayer54", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer55  = declare1D(m_name, "minMSOSLocalDistanceAtLayer55", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer56  = declare1D(m_name, "minMSOSLocalDistanceAtLayer56", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer57  = declare1D(m_name, "minMSOSLocalDistanceAtLayer57", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer58  = declare1D(m_name, "minMSOSLocalDistanceAtLayer58", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer59  = declare1D(m_name, "minMSOSLocalDistanceAtLayer59", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer60  = declare1D(m_name, "minMSOSLocalDistanceAtLayer60", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer61  = declare1D(m_name, "minMSOSLocalDistanceAtLayer61", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer62  = declare1D(m_name, "minMSOSLocalDistanceAtLayer62", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer63  = declare1D(m_name, "minMSOSLocalDistanceAtLayer63", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer64  = declare1D(m_name, "minMSOSLocalDistanceAtLayer64", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer65  = declare1D(m_name, "minMSOSLocalDistanceAtLayer65", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer66  = declare1D(m_name, "minMSOSLocalDistanceAtLayer66", "minimum local distance of MSOS [mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayer67  = declare1D(m_name, "minMSOSLocalDistanceAtLayer67", "minimum local distance of MSOS [mm]",100,0,10.0);

  m_minMSOSLocalDistanceAtLayerX0   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX0", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX1   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX1", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX2   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX2", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX3   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX3", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX4   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX4", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX5   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX5", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX6   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX6", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX7   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX7", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX8   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX8", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX9   = declare1D(m_name, "minMSOSLocalDistanceAtLayerX9", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX10  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX10", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX11  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX11", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX12  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX12", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX13  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX13", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX14  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX14", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX15  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX15", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX16  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX16", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX17  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX17", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX18  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX18", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX19  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX19", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX20  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX20", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX21  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX21", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX22  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX22", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX23  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX23", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX24  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX24", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX25  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX25", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX26  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX26", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX27  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX27", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX28  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX28", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX29  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX29", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX30  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX30", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX31  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX31", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX32  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX32", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX33  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX33", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX34  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX34", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX35  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX35", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX36  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX36", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX37  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX37", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX38  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX38", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX39  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX39", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX40  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX40", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX41  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX41", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX42  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX42", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX43  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX43", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX44  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX44", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX45  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX45", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX46  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX46", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX47  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX47", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX48  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX48", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX49  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX49", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX50  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX50", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX51  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX51", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX52  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX52", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX53  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX53", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX54  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX54", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX55  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX55", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX56  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX56", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX57  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX57", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX58  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX58", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX59  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX59", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX60  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX60", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX61  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX61", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX62  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX62", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX63  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX63", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX64  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX64", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX65  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX65", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX66  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX66", "minimum local distance of MSOS X[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerX67  = declare1D(m_name, "minMSOSLocalDistanceAtLayerX67", "minimum local distance of MSOS X[mm]",100,0,10.0);

  m_minMSOSLocalDistanceAtLayerY0   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY0", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY1   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY1", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY2   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY2", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY3   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY3", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY4   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY4", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY5   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY5", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY6   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY6", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY7   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY7", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY8   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY8", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY9   = declare1D(m_name, "minMSOSLocalDistanceAtLayerY9", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY10  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY10", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY11  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY11", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY12  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY12", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY13  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY13", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY14  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY14", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY15  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY15", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY16  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY16", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY17  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY17", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY18  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY18", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY19  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY19", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY20  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY20", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY21  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY21", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY22  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY22", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY23  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY23", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY24  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY24", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY25  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY25", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY26  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY26", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY27  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY27", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY28  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY28", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY29  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY29", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY30  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY30", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY31  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY31", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY32  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY32", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY33  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY33", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY34  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY34", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY35  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY35", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY36  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY36", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY37  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY37", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY38  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY38", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY39  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY39", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY40  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY40", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY41  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY41", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY42  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY42", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY43  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY43", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY44  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY44", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY45  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY45", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY46  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY46", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY47  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY47", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY48  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY48", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY49  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY49", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY50  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY50", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY51  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY51", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY52  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY52", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY53  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY53", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY54  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY54", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY55  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY55", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY56  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY56", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY57  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY57", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY58  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY58", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY59  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY59", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY60  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY60", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY61  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY61", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY62  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY62", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY63  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY63", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY64  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY64", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY65  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY65", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY66  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY66", "minimum local distance of MSOS Y[mm]",100,0,10.0);
  m_minMSOSLocalDistanceAtLayerY67  = declare1D(m_name, "minMSOSLocalDistanceAtLayerY67", "minimum local distance of MSOS Y[mm]",100,0,10.0);

  m_minClusLocalDistanceAtLayer  = declare1D(m_name, "minClusLocalDistanceAtLayer", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayerX = declare1D(m_name, "minClusLocalDistanceAtLayerX", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayerY = declare1D(m_name, "minClusLocalDistanceAtLayerY", "minimum local distance of clusters [mm]",100,0,10.0);

  m_minClusLocalDistanceAtLayer0   = declare1D(m_name, "minClusLocalDistanceAtLayer0", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer1   = declare1D(m_name, "minClusLocalDistanceAtLayer1", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer2   = declare1D(m_name, "minClusLocalDistanceAtLayer2", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer3   = declare1D(m_name, "minClusLocalDistanceAtLayer3", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer4   = declare1D(m_name, "minClusLocalDistanceAtLayer4", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer5   = declare1D(m_name, "minClusLocalDistanceAtLayer5", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer6   = declare1D(m_name, "minClusLocalDistanceAtLayer6", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer7   = declare1D(m_name, "minClusLocalDistanceAtLayer7", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer8   = declare1D(m_name, "minClusLocalDistanceAtLayer8", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer9   = declare1D(m_name, "minClusLocalDistanceAtLayer9", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer10  = declare1D(m_name, "minClusLocalDistanceAtLayer10", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer11  = declare1D(m_name, "minClusLocalDistanceAtLayer11", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer12  = declare1D(m_name, "minClusLocalDistanceAtLayer12", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer13  = declare1D(m_name, "minClusLocalDistanceAtLayer13", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer14  = declare1D(m_name, "minClusLocalDistanceAtLayer14", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer15  = declare1D(m_name, "minClusLocalDistanceAtLayer15", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer16  = declare1D(m_name, "minClusLocalDistanceAtLayer16", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer17  = declare1D(m_name, "minClusLocalDistanceAtLayer17", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer18  = declare1D(m_name, "minClusLocalDistanceAtLayer18", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer19  = declare1D(m_name, "minClusLocalDistanceAtLayer19", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer20  = declare1D(m_name, "minClusLocalDistanceAtLayer20", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer21  = declare1D(m_name, "minClusLocalDistanceAtLayer21", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer22  = declare1D(m_name, "minClusLocalDistanceAtLayer22", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer23  = declare1D(m_name, "minClusLocalDistanceAtLayer23", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer24  = declare1D(m_name, "minClusLocalDistanceAtLayer24", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer25  = declare1D(m_name, "minClusLocalDistanceAtLayer25", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer26  = declare1D(m_name, "minClusLocalDistanceAtLayer26", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer27  = declare1D(m_name, "minClusLocalDistanceAtLayer27", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer28  = declare1D(m_name, "minClusLocalDistanceAtLayer28", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer29  = declare1D(m_name, "minClusLocalDistanceAtLayer29", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer30  = declare1D(m_name, "minClusLocalDistanceAtLayer30", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer31  = declare1D(m_name, "minClusLocalDistanceAtLayer31", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer32  = declare1D(m_name, "minClusLocalDistanceAtLayer32", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer33  = declare1D(m_name, "minClusLocalDistanceAtLayer33", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer34  = declare1D(m_name, "minClusLocalDistanceAtLayer34", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer35  = declare1D(m_name, "minClusLocalDistanceAtLayer35", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer36  = declare1D(m_name, "minClusLocalDistanceAtLayer36", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer37  = declare1D(m_name, "minClusLocalDistanceAtLayer37", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer38  = declare1D(m_name, "minClusLocalDistanceAtLayer38", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer39  = declare1D(m_name, "minClusLocalDistanceAtLayer39", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer40  = declare1D(m_name, "minClusLocalDistanceAtLayer40", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer41  = declare1D(m_name, "minClusLocalDistanceAtLayer41", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer42  = declare1D(m_name, "minClusLocalDistanceAtLayer42", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer43  = declare1D(m_name, "minClusLocalDistanceAtLayer43", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer44  = declare1D(m_name, "minClusLocalDistanceAtLayer44", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer45  = declare1D(m_name, "minClusLocalDistanceAtLayer45", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer46  = declare1D(m_name, "minClusLocalDistanceAtLayer46", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer47  = declare1D(m_name, "minClusLocalDistanceAtLayer47", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer48  = declare1D(m_name, "minClusLocalDistanceAtLayer48", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer49  = declare1D(m_name, "minClusLocalDistanceAtLayer49", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer50  = declare1D(m_name, "minClusLocalDistanceAtLayer50", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer51  = declare1D(m_name, "minClusLocalDistanceAtLayer51", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer52  = declare1D(m_name, "minClusLocalDistanceAtLayer52", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer53  = declare1D(m_name, "minClusLocalDistanceAtLayer53", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer54  = declare1D(m_name, "minClusLocalDistanceAtLayer54", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer55  = declare1D(m_name, "minClusLocalDistanceAtLayer55", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer56  = declare1D(m_name, "minClusLocalDistanceAtLayer56", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer57  = declare1D(m_name, "minClusLocalDistanceAtLayer57", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer58  = declare1D(m_name, "minClusLocalDistanceAtLayer58", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer59  = declare1D(m_name, "minClusLocalDistanceAtLayer59", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer60  = declare1D(m_name, "minClusLocalDistanceAtLayer60", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer61  = declare1D(m_name, "minClusLocalDistanceAtLayer61", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer62  = declare1D(m_name, "minClusLocalDistanceAtLayer62", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer63  = declare1D(m_name, "minClusLocalDistanceAtLayer63", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer64  = declare1D(m_name, "minClusLocalDistanceAtLayer64", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer65  = declare1D(m_name, "minClusLocalDistanceAtLayer65", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer66  = declare1D(m_name, "minClusLocalDistanceAtLayer66", "minimum local distance of clusters [mm]",100,0,10.0);
  m_minClusLocalDistanceAtLayer67  = declare1D(m_name, "minClusLocalDistanceAtLayer67", "minimum local distance of clusters [mm]",100,0,10.0);

  m_minClusGlobalDistanceAtLayer  = declare1D(m_name, "minClusGlobalDistanceAtLayer", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayerX = declare1D(m_name, "minClusGlobalDistanceAtLayerX", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayerY = declare1D(m_name, "minClusGlobalDistanceAtLayerY", "minimum distance of clusters [mm]",100,0,10.0);

  m_minClusGlobalDistanceAtLayer0 = declare1D(m_name, "minClusGlobalDistanceAtLayer0", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer1 = declare1D(m_name, "minClusGlobalDistanceAtLayer1", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer2 = declare1D(m_name, "minClusGlobalDistanceAtLayer2", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer3 = declare1D(m_name, "minClusGlobalDistanceAtLayer3", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer4 = declare1D(m_name, "minClusGlobalDistanceAtLayer4", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer5 = declare1D(m_name, "minClusGlobalDistanceAtLayer5", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer6 = declare1D(m_name, "minClusGlobalDistanceAtLayer6", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer7 = declare1D(m_name, "minClusGlobalDistanceAtLayer7", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer8 = declare1D(m_name, "minClusGlobalDistanceAtLayer8", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer9 = declare1D(m_name, "minClusGlobalDistanceAtLayer9", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer10 = declare1D(m_name, "minClusGlobalDistanceAtLayer10", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer11 = declare1D(m_name, "minClusGlobalDistanceAtLayer11", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer12 = declare1D(m_name, "minClusGlobalDistanceAtLayer12", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer13 = declare1D(m_name, "minClusGlobalDistanceAtLayer13", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer14 = declare1D(m_name, "minClusGlobalDistanceAtLayer14", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer15 = declare1D(m_name, "minClusGlobalDistanceAtLayer15", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer16 = declare1D(m_name, "minClusGlobalDistanceAtLayer16", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer17 = declare1D(m_name, "minClusGlobalDistanceAtLayer17", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer18 = declare1D(m_name, "minClusGlobalDistanceAtLayer18", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer19 = declare1D(m_name, "minClusGlobalDistanceAtLayer19", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer20 = declare1D(m_name, "minClusGlobalDistanceAtLayer20", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer21 = declare1D(m_name, "minClusGlobalDistanceAtLayer21", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer22 = declare1D(m_name, "minClusGlobalDistanceAtLayer22", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer23 = declare1D(m_name, "minClusGlobalDistanceAtLayer23", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer24 = declare1D(m_name, "minClusGlobalDistanceAtLayer24", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer25 = declare1D(m_name, "minClusGlobalDistanceAtLayer25", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer26 = declare1D(m_name, "minClusGlobalDistanceAtLayer26", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer27 = declare1D(m_name, "minClusGlobalDistanceAtLayer27", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer28 = declare1D(m_name, "minClusGlobalDistanceAtLayer28", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer29 = declare1D(m_name, "minClusGlobalDistanceAtLayer29", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer30 = declare1D(m_name, "minClusGlobalDistanceAtLayer30", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer31 = declare1D(m_name, "minClusGlobalDistanceAtLayer31", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer32 = declare1D(m_name, "minClusGlobalDistanceAtLayer32", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer33 = declare1D(m_name, "minClusGlobalDistanceAtLayer33", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer34 = declare1D(m_name, "minClusGlobalDistanceAtLayer34", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer35 = declare1D(m_name, "minClusGlobalDistanceAtLayer35", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer36 = declare1D(m_name, "minClusGlobalDistanceAtLayer36", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer37 = declare1D(m_name, "minClusGlobalDistanceAtLayer37", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer38 = declare1D(m_name, "minClusGlobalDistanceAtLayer38", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer39 = declare1D(m_name, "minClusGlobalDistanceAtLayer39", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer40 = declare1D(m_name, "minClusGlobalDistanceAtLayer40", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer41 = declare1D(m_name, "minClusGlobalDistanceAtLayer41", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer42 = declare1D(m_name, "minClusGlobalDistanceAtLayer42", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer43 = declare1D(m_name, "minClusGlobalDistanceAtLayer43", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer44 = declare1D(m_name, "minClusGlobalDistanceAtLayer44", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer45 = declare1D(m_name, "minClusGlobalDistanceAtLayer45", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer46 = declare1D(m_name, "minClusGlobalDistanceAtLayer46", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer47 = declare1D(m_name, "minClusGlobalDistanceAtLayer47", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer48 = declare1D(m_name, "minClusGlobalDistanceAtLayer48", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer49 = declare1D(m_name, "minClusGlobalDistanceAtLayer49", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer50 = declare1D(m_name, "minClusGlobalDistanceAtLayer50", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer51 = declare1D(m_name, "minClusGlobalDistanceAtLayer51", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer52 = declare1D(m_name, "minClusGlobalDistanceAtLayer52", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer53 = declare1D(m_name, "minClusGlobalDistanceAtLayer53", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer54 = declare1D(m_name, "minClusGlobalDistanceAtLayer54", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer55 = declare1D(m_name, "minClusGlobalDistanceAtLayer55", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer56 = declare1D(m_name, "minClusGlobalDistanceAtLayer56", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer57 = declare1D(m_name, "minClusGlobalDistanceAtLayer57", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer58 = declare1D(m_name, "minClusGlobalDistanceAtLayer58", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer59 = declare1D(m_name, "minClusGlobalDistanceAtLayer59", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer60 = declare1D(m_name, "minClusGlobalDistanceAtLayer60", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer61 = declare1D(m_name, "minClusGlobalDistanceAtLayer61", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer62 = declare1D(m_name, "minClusGlobalDistanceAtLayer62", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer63 = declare1D(m_name, "minClusGlobalDistanceAtLayer63", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer64 = declare1D(m_name, "minClusGlobalDistanceAtLayer64", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer65 = declare1D(m_name, "minClusGlobalDistanceAtLayer65", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer66 = declare1D(m_name, "minClusGlobalDistanceAtLayer66", "minimum distance of clusters [mm]",100,0,10.0);
  m_minClusGlobalDistanceAtLayer67 = declare1D(m_name, "minClusGlobalDistanceAtLayer67", "minimum distance of clusters [mm]",100,0,10.0);

  m_pixLocal_shiftX = declare1D(m_name, "pixLocal_shiftX", "cluster shift from MSOS x[mm]",100,-10,10);
  m_pixLocal_shiftY = declare1D(m_name, "pixLocal_shiftY", "cluster shift from MSOS y[mm]",100,-10,10);

  m_PixelLayout_xy = declare2D(m_name,"PixelLayout_xy","x[mm]","y[mm]",150,-150,150,150,-150,150);
  m_PixelLayout_zr = declare2D(m_name,"PixelLayout_zr","z[mm]","r[mm]",200,-800,800,150,0,150);
  m_PixelLayout2Hits_xy = declare2D(m_name,"PixelLayout2Hits_xy","x[mm]","y[mm]",150,-150,150,150,-150,150);
  m_PixelLayout2Hits_zr = declare2D(m_name,"PixelLayout2Hits_zr","z[mm]","r[mm]",200,-800,800,150,0,150);

  m_PixelLayoutZoom_xy = declare2D(m_name,"PixelLayoutZoom_xy","x[mm]","y[mm]",240,-60,60,240,-60,60);
  m_PixelLayoutZoom_zr = declare2D(m_name,"PixelLayoutZoom_zr","z[mm]","r[mm]",200,-400,400,160,20,60);
  m_PixelLayout2HitsZoom_xy = declare2D(m_name,"PixelLayout2HitsZoom_xy","x[mm]","y[mm]",240,-60,60,240,-60,60);
  m_PixelLayout2HitsZoom_zr = declare2D(m_name,"PixelLayout2HitsZoom_zr","z[mm]","r[mm]",200,-400,400,160,20,60);

  m_etaWeight = declare2D(m_name,"etaWeight","#eta(leading track)","#eta(sub-leading track)",60,-3.0,3.0,60,-3.0,3.0);
  m_ptWeight  = declare2D(m_name,"ptWeight","p_{T}(leading track) [GeV]","p_{T}(sub-leading track) [GeV]",200,0.0,50.0,200,0.0,50.0);

  m_mindRJet1    = declare1D(m_name, "mindRJet1",    "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P0  = declare1D(m_name, "mindRJet1P0",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P1  = declare1D(m_name, "mindRJet1P1",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P2  = declare1D(m_name, "mindRJet1P2",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P3  = declare1D(m_name, "mindRJet1P3",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P4  = declare1D(m_name, "mindRJet1P4",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P5  = declare1D(m_name, "mindRJet1P5",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P6  = declare1D(m_name, "mindRJet1P6",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P7  = declare1D(m_name, "mindRJet1P7",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P8  = declare1D(m_name, "mindRJet1P8",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P9  = declare1D(m_name, "mindRJet1P9",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P10 = declare1D(m_name, "mindRJet1P10", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P11 = declare1D(m_name, "mindRJet1P11", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P12 = declare1D(m_name, "mindRJet1P12", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P13 = declare1D(m_name, "mindRJet1P13", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P14 = declare1D(m_name, "mindRJet1P14", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P15 = declare1D(m_name, "mindRJet1P15", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P16 = declare1D(m_name, "mindRJet1P16", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P17 = declare1D(m_name, "mindRJet1P17", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P18 = declare1D(m_name, "mindRJet1P18", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P19 = declare1D(m_name, "mindRJet1P19", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P20 = declare1D(m_name, "mindRJet1P20", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P21 = declare1D(m_name, "mindRJet1P21", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P22 = declare1D(m_name, "mindRJet1P22", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P23 = declare1D(m_name, "mindRJet1P23", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P24 = declare1D(m_name, "mindRJet1P24", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P25 = declare1D(m_name, "mindRJet1P25", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P26 = declare1D(m_name, "mindRJet1P26", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P27 = declare1D(m_name, "mindRJet1P27", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P28 = declare1D(m_name, "mindRJet1P28", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P29 = declare1D(m_name, "mindRJet1P29", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P30 = declare1D(m_name, "mindRJet1P30", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P31 = declare1D(m_name, "mindRJet1P31", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P32 = declare1D(m_name, "mindRJet1P32", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P33 = declare1D(m_name, "mindRJet1P33", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P34 = declare1D(m_name, "mindRJet1P34", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P35 = declare1D(m_name, "mindRJet1P35", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P36 = declare1D(m_name, "mindRJet1P36", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P37 = declare1D(m_name, "mindRJet1P37", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P38 = declare1D(m_name, "mindRJet1P38", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P39 = declare1D(m_name, "mindRJet1P39", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P40 = declare1D(m_name, "mindRJet1P40", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P41 = declare1D(m_name, "mindRJet1P41", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P42 = declare1D(m_name, "mindRJet1P42", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P43 = declare1D(m_name, "mindRJet1P43", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P44 = declare1D(m_name, "mindRJet1P44", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P45 = declare1D(m_name, "mindRJet1P45", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P46 = declare1D(m_name, "mindRJet1P46", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P47 = declare1D(m_name, "mindRJet1P47", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P48 = declare1D(m_name, "mindRJet1P48", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P49 = declare1D(m_name, "mindRJet1P49", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P50 = declare1D(m_name, "mindRJet1P50", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P51 = declare1D(m_name, "mindRJet1P51", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P52 = declare1D(m_name, "mindRJet1P52", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P53 = declare1D(m_name, "mindRJet1P53", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P54 = declare1D(m_name, "mindRJet1P54", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P55 = declare1D(m_name, "mindRJet1P55", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P56 = declare1D(m_name, "mindRJet1P56", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P57 = declare1D(m_name, "mindRJet1P57", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P58 = declare1D(m_name, "mindRJet1P58", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P59 = declare1D(m_name, "mindRJet1P59", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P60 = declare1D(m_name, "mindRJet1P60", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P61 = declare1D(m_name, "mindRJet1P61", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P62 = declare1D(m_name, "mindRJet1P62", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P63 = declare1D(m_name, "mindRJet1P63", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P64 = declare1D(m_name, "mindRJet1P64", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P65 = declare1D(m_name, "mindRJet1P65", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P66 = declare1D(m_name, "mindRJet1P66", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet1P67 = declare1D(m_name, "mindRJet1P67", "#Delta R(track,jet)",20,0,0.2);

  m_mindRJet2    = declare1D(m_name, "mindRJet2",    "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P0  = declare1D(m_name, "mindRJet2P0",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P1  = declare1D(m_name, "mindRJet2P1",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P2  = declare1D(m_name, "mindRJet2P2",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P3  = declare1D(m_name, "mindRJet2P3",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P4  = declare1D(m_name, "mindRJet2P4",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P5  = declare1D(m_name, "mindRJet2P5",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P6  = declare1D(m_name, "mindRJet2P6",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P7  = declare1D(m_name, "mindRJet2P7",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P8  = declare1D(m_name, "mindRJet2P8",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P9  = declare1D(m_name, "mindRJet2P9",  "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P10 = declare1D(m_name, "mindRJet2P10", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P11 = declare1D(m_name, "mindRJet2P11", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P12 = declare1D(m_name, "mindRJet2P12", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P13 = declare1D(m_name, "mindRJet2P13", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P14 = declare1D(m_name, "mindRJet2P14", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P15 = declare1D(m_name, "mindRJet2P15", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P16 = declare1D(m_name, "mindRJet2P16", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P17 = declare1D(m_name, "mindRJet2P17", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P18 = declare1D(m_name, "mindRJet2P18", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P19 = declare1D(m_name, "mindRJet2P19", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P20 = declare1D(m_name, "mindRJet2P20", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P21 = declare1D(m_name, "mindRJet2P21", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P22 = declare1D(m_name, "mindRJet2P22", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P23 = declare1D(m_name, "mindRJet2P23", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P24 = declare1D(m_name, "mindRJet2P24", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P25 = declare1D(m_name, "mindRJet2P25", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P26 = declare1D(m_name, "mindRJet2P26", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P27 = declare1D(m_name, "mindRJet2P27", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P28 = declare1D(m_name, "mindRJet2P28", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P29 = declare1D(m_name, "mindRJet2P29", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P30 = declare1D(m_name, "mindRJet2P30", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P31 = declare1D(m_name, "mindRJet2P31", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P32 = declare1D(m_name, "mindRJet2P32", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P33 = declare1D(m_name, "mindRJet2P33", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P34 = declare1D(m_name, "mindRJet2P34", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P35 = declare1D(m_name, "mindRJet2P35", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P36 = declare1D(m_name, "mindRJet2P36", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P37 = declare1D(m_name, "mindRJet2P37", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P38 = declare1D(m_name, "mindRJet2P38", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P39 = declare1D(m_name, "mindRJet2P39", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P40 = declare1D(m_name, "mindRJet2P40", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P41 = declare1D(m_name, "mindRJet2P41", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P42 = declare1D(m_name, "mindRJet2P42", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P43 = declare1D(m_name, "mindRJet2P43", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P44 = declare1D(m_name, "mindRJet2P44", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P45 = declare1D(m_name, "mindRJet2P45", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P46 = declare1D(m_name, "mindRJet2P46", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P47 = declare1D(m_name, "mindRJet2P47", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P48 = declare1D(m_name, "mindRJet2P48", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P49 = declare1D(m_name, "mindRJet2P49", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P50 = declare1D(m_name, "mindRJet2P50", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P51 = declare1D(m_name, "mindRJet2P51", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P52 = declare1D(m_name, "mindRJet2P52", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P53 = declare1D(m_name, "mindRJet2P53", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P54 = declare1D(m_name, "mindRJet2P54", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P55 = declare1D(m_name, "mindRJet2P55", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P56 = declare1D(m_name, "mindRJet2P56", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P57 = declare1D(m_name, "mindRJet2P57", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P58 = declare1D(m_name, "mindRJet2P58", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P59 = declare1D(m_name, "mindRJet2P59", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P60 = declare1D(m_name, "mindRJet2P60", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P61 = declare1D(m_name, "mindRJet2P61", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P62 = declare1D(m_name, "mindRJet2P62", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P63 = declare1D(m_name, "mindRJet2P63", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P64 = declare1D(m_name, "mindRJet2P64", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P65 = declare1D(m_name, "mindRJet2P65", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P66 = declare1D(m_name, "mindRJet2P66", "#Delta R(track,jet)",20,0,0.2);
  m_mindRJet2P67 = declare1D(m_name, "mindRJet2P67", "#Delta R(track,jet)",20,0,0.2);

  const int nbins=21;
  double xbins[22] = {10,14,20,30,40,50,65,80,100,140,170,200,240,300,340,400,500,600,700,800,900,1000};

  m_jetpT_locT = declare2D(m_name,"jetpT_locT","Jet p_{T} [GeV]","#delta_{T}^{min}",nbins, xbins,1000,0.0,10.0);
  m_jetpT_locX = declare2D(m_name,"jetpT_locX","Jet p_{T} [GeV]","#delta_{X}^{min}",nbins, xbins,1000,0.0,10.0);
  m_jetpT_locY = declare2D(m_name,"jetpT_locY","Jet p_{T} [GeV]","#delta_{Y}^{min}",nbins, xbins,1000,0.0,10.0);

  m_jetpT_locTNP = declare2D(m_name,"jetpT_locTNP","Jet p_{T} [GeV]","#delta_{T}^{min}",nbins, xbins,1000,0.0,10.0);
  m_jetpT_locXNP = declare2D(m_name,"jetpT_locXNP","Jet p_{T} [GeV]","#delta_{X}^{min}",nbins, xbins,1000,0.0,10.0);
  m_jetpT_locYNP = declare2D(m_name,"jetpT_locYNP","Jet p_{T} [GeV]","#delta_{Y}^{min}",nbins, xbins,1000,0.0,10.0);

  m_jetpT_locTSP = declare2D(m_name,"jetpT_locTSP","Jet p_{T} [GeV]","#delta_{T}^{min}",nbins, xbins,1000,0.0,10.0);
  m_jetpT_locXSP = declare2D(m_name,"jetpT_locXSP","Jet p_{T} [GeV]","#delta_{X}^{min}",nbins, xbins,1000,0.0,10.0);
  m_jetpT_locYSP = declare2D(m_name,"jetpT_locYSP","Jet p_{T} [GeV]","#delta_{Y}^{min}",nbins, xbins,1000,0.0,10.0);

} // BookHists


void TwoTrackHists::FillHists(const xAOD::TrackParticle* trk, float weight) const {

  uint8_t getInt(0);   // for accessing summary information
  float   getFlt(0.0); // for accessing summary information

  if (trk->isAvailable<int>("ClusterPatternIndex")) {
    int clusPatternIndex = trk->auxdata<int>("ClusterPatternIndex");
    m_TargetLayerOfInterest -> Fill(1.0*trk->auxdata<int>("TargetLayerOfInterest"),weight);
    m_ClusterPatternIndex -> Fill(1.0*clusPatternIndex,weight);

    float minMSOSLocalDistance  = trk->auxdata<float>("minMSOSLocalDistanceAtLayer");
    float minMSOSLocalDistanceX = trk->auxdata<float>("minMSOSLocalDistanceAtLayerX");
    float minMSOSLocalDistanceY = trk->auxdata<float>("minMSOSLocalDistanceAtLayerY");
    if (minMSOSLocalDistance==0) { minMSOSLocalDistance=100.0; }

    m_minMSOSLocalDistanceAtLayer  -> Fill(minMSOSLocalDistance,weight);
    m_minMSOSLocalDistanceAtLayerX -> Fill(minMSOSLocalDistanceX,weight);
    m_minMSOSLocalDistanceAtLayerY -> Fill(minMSOSLocalDistanceY,weight);
    if (clusPatternIndex== 0) { 
      m_minMSOSLocalDistanceAtLayer0  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX0  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY0  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 1) { 
      m_minMSOSLocalDistanceAtLayer1  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX1  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY1  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 2) { 
      m_minMSOSLocalDistanceAtLayer2  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX2  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY2  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 3) { 
      m_minMSOSLocalDistanceAtLayer3  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX3  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY3  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 4) { 
      m_minMSOSLocalDistanceAtLayer4  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX4  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY4  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 5) { 
      m_minMSOSLocalDistanceAtLayer5  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX5  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY5  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 6) { 
      m_minMSOSLocalDistanceAtLayer6  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX6  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY6  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 7) { 
      m_minMSOSLocalDistanceAtLayer7  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX7  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY7  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 8) { 
      m_minMSOSLocalDistanceAtLayer8  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX8  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY8  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex== 9) { 
      m_minMSOSLocalDistanceAtLayer9  -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX9  -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY9  -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==10) { 
      m_minMSOSLocalDistanceAtLayer10 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX10 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY10 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==11) { 
      m_minMSOSLocalDistanceAtLayer11 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX11 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY11 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==12) { 
      m_minMSOSLocalDistanceAtLayer12 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX12 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY12 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==13) { 
      m_minMSOSLocalDistanceAtLayer13 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX13 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY13 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==14) { 
      m_minMSOSLocalDistanceAtLayer14 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX14 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY14 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==15) { 
      m_minMSOSLocalDistanceAtLayer15 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX15 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY15 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==16) { 
      m_minMSOSLocalDistanceAtLayer16 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX16 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY16 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==17) { 
      m_minMSOSLocalDistanceAtLayer17 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX17 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY17 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==18) { 
      m_minMSOSLocalDistanceAtLayer18 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX18 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY18 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==19) { 
      m_minMSOSLocalDistanceAtLayer19 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX19 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY19 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==20) { 
      m_minMSOSLocalDistanceAtLayer20 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX20 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY20 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==21) { 
      m_minMSOSLocalDistanceAtLayer21 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX21 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY21 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==22) { 
      m_minMSOSLocalDistanceAtLayer22 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX22 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY22 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==23) { 
      m_minMSOSLocalDistanceAtLayer23 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX23 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY23 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==24) { 
      m_minMSOSLocalDistanceAtLayer24 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX24 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY24 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==25) { 
      m_minMSOSLocalDistanceAtLayer25 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX25 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY25 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==26) { 
      m_minMSOSLocalDistanceAtLayer26 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX26 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY26 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==27) { 
      m_minMSOSLocalDistanceAtLayer27 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX27 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY27 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==28) { 
      m_minMSOSLocalDistanceAtLayer28 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX28 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY28 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==29) { 
      m_minMSOSLocalDistanceAtLayer29 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX29 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY29 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==30) { 
      m_minMSOSLocalDistanceAtLayer30 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX30 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY30 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==31) { 
      m_minMSOSLocalDistanceAtLayer31 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX31 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY31 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==32) { 
      m_minMSOSLocalDistanceAtLayer32 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX32 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY32 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==33) { 
      m_minMSOSLocalDistanceAtLayer33 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX33 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY33 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==34) { 
      m_minMSOSLocalDistanceAtLayer34 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX34 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY34 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==35) { 
      m_minMSOSLocalDistanceAtLayer35 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX35 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY35 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==36) { 
      m_minMSOSLocalDistanceAtLayer36 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX36 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY36 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==37) { 
      m_minMSOSLocalDistanceAtLayer37 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX37 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY37 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==38) { 
      m_minMSOSLocalDistanceAtLayer38 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX38 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY38 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==39) { 
      m_minMSOSLocalDistanceAtLayer39 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX39 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY39 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==40) { 
      m_minMSOSLocalDistanceAtLayer40 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX40 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY40 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==41) { 
      m_minMSOSLocalDistanceAtLayer41 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX41 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY41 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==42) { 
      m_minMSOSLocalDistanceAtLayer42 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX42 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY42 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==43) { 
      m_minMSOSLocalDistanceAtLayer43 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX43 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY43 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==44) { 
      m_minMSOSLocalDistanceAtLayer44 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX44 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY44 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==45) { 
      m_minMSOSLocalDistanceAtLayer45 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX45 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY45 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==46) { 
      m_minMSOSLocalDistanceAtLayer46 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX46 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY46 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==47) { 
      m_minMSOSLocalDistanceAtLayer47 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX47 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY47 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==48) { 
      m_minMSOSLocalDistanceAtLayer48 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX48 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY48 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==49) { 
      m_minMSOSLocalDistanceAtLayer49 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX49 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY49 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==50) { 
      m_minMSOSLocalDistanceAtLayer50 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX50 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY50 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==51) { 
      m_minMSOSLocalDistanceAtLayer51 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX51 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY51 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==52) { 
      m_minMSOSLocalDistanceAtLayer52 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX52 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY52 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==53) { 
      m_minMSOSLocalDistanceAtLayer53 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX53 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY53 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==54) { 
      m_minMSOSLocalDistanceAtLayer54 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX54 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY54 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==55) { 
      m_minMSOSLocalDistanceAtLayer55 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX55 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY55 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==56) { 
      m_minMSOSLocalDistanceAtLayer56 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX56 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY56 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==57) { 
      m_minMSOSLocalDistanceAtLayer57 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX57 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY57 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==58) { 
      m_minMSOSLocalDistanceAtLayer58 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX58 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY58 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==59) { 
      m_minMSOSLocalDistanceAtLayer59 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX59 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY59 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==60) { 
      m_minMSOSLocalDistanceAtLayer60 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX60 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY60 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==61) { 
      m_minMSOSLocalDistanceAtLayer61 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX61 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY61 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==62) { 
      m_minMSOSLocalDistanceAtLayer62 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX62 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY62 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==63) { 
      m_minMSOSLocalDistanceAtLayer63 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX63 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY63 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==64) { 
      m_minMSOSLocalDistanceAtLayer64 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX64 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY64 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==65) { 
      m_minMSOSLocalDistanceAtLayer65 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX65 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY65 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==66) { 
      m_minMSOSLocalDistanceAtLayer66 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX66 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY66 -> Fill(minMSOSLocalDistanceY,weight); 
    }
    if (clusPatternIndex==67) { 
      m_minMSOSLocalDistanceAtLayer67 -> Fill(minMSOSLocalDistance,weight); 
      m_minMSOSLocalDistanceAtLayerX67 -> Fill(minMSOSLocalDistanceX,weight); 
      m_minMSOSLocalDistanceAtLayerY67 -> Fill(minMSOSLocalDistanceY,weight); 
    }

    float minClusLocalDistance = trk->auxdata<float>("minClusLocalDistanceAtLayer");
    if (minClusLocalDistance==0) { minClusLocalDistance=100.0; }

    m_minClusLocalDistanceAtLayer  -> Fill(minClusLocalDistance,weight);
    m_minClusLocalDistanceAtLayerX -> Fill(trk->auxdata<float>("minClusLocalDistanceAtLayerX"),weight);
    m_minClusLocalDistanceAtLayerY -> Fill(trk->auxdata<float>("minClusLocalDistanceAtLayerY"),weight);
    if (clusPatternIndex== 0) { m_minClusLocalDistanceAtLayer0  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 1) { m_minClusLocalDistanceAtLayer1  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 2) { m_minClusLocalDistanceAtLayer2  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 3) { m_minClusLocalDistanceAtLayer3  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 4) { m_minClusLocalDistanceAtLayer4  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 5) { m_minClusLocalDistanceAtLayer5  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 6) { m_minClusLocalDistanceAtLayer6  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 7) { m_minClusLocalDistanceAtLayer7  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 8) { m_minClusLocalDistanceAtLayer8  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 9) { m_minClusLocalDistanceAtLayer9  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==10) { m_minClusLocalDistanceAtLayer10 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==11) { m_minClusLocalDistanceAtLayer11 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==12) { m_minClusLocalDistanceAtLayer12 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==13) { m_minClusLocalDistanceAtLayer13 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==14) { m_minClusLocalDistanceAtLayer14 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==15) { m_minClusLocalDistanceAtLayer15 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==16) { m_minClusLocalDistanceAtLayer16 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==17) { m_minClusLocalDistanceAtLayer17 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==18) { m_minClusLocalDistanceAtLayer18 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==19) { m_minClusLocalDistanceAtLayer19 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==20) { m_minClusLocalDistanceAtLayer20 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==21) { m_minClusLocalDistanceAtLayer21 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==22) { m_minClusLocalDistanceAtLayer22 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==23) { m_minClusLocalDistanceAtLayer23 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==24) { m_minClusLocalDistanceAtLayer24 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==25) { m_minClusLocalDistanceAtLayer25 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==26) { m_minClusLocalDistanceAtLayer26 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==27) { m_minClusLocalDistanceAtLayer27 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==28) { m_minClusLocalDistanceAtLayer28 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==29) { m_minClusLocalDistanceAtLayer29 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==30) { m_minClusLocalDistanceAtLayer30 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==31) { m_minClusLocalDistanceAtLayer31 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==32) { m_minClusLocalDistanceAtLayer32 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==33) { m_minClusLocalDistanceAtLayer33 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==34) { m_minClusLocalDistanceAtLayer34 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==35) { m_minClusLocalDistanceAtLayer35 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==36) { m_minClusLocalDistanceAtLayer36 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==37) { m_minClusLocalDistanceAtLayer37 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==38) { m_minClusLocalDistanceAtLayer38 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==39) { m_minClusLocalDistanceAtLayer39 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==40) { m_minClusLocalDistanceAtLayer40 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==41) { m_minClusLocalDistanceAtLayer41 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==42) { m_minClusLocalDistanceAtLayer42 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==43) { m_minClusLocalDistanceAtLayer43 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==44) { m_minClusLocalDistanceAtLayer44 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==45) { m_minClusLocalDistanceAtLayer45 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==46) { m_minClusLocalDistanceAtLayer46 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==47) { m_minClusLocalDistanceAtLayer47 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==48) { m_minClusLocalDistanceAtLayer48 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==49) { m_minClusLocalDistanceAtLayer49 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==50) { m_minClusLocalDistanceAtLayer50 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==51) { m_minClusLocalDistanceAtLayer51 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==52) { m_minClusLocalDistanceAtLayer52 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==53) { m_minClusLocalDistanceAtLayer53 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==54) { m_minClusLocalDistanceAtLayer54 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==55) { m_minClusLocalDistanceAtLayer55 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==56) { m_minClusLocalDistanceAtLayer56 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==57) { m_minClusLocalDistanceAtLayer57 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==58) { m_minClusLocalDistanceAtLayer58 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==59) { m_minClusLocalDistanceAtLayer59 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==60) { m_minClusLocalDistanceAtLayer60 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==61) { m_minClusLocalDistanceAtLayer61 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==62) { m_minClusLocalDistanceAtLayer62 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==63) { m_minClusLocalDistanceAtLayer63 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==64) { m_minClusLocalDistanceAtLayer64 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==65) { m_minClusLocalDistanceAtLayer65 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==66) { m_minClusLocalDistanceAtLayer66 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==67) { m_minClusLocalDistanceAtLayer67 -> Fill(minClusLocalDistance,weight); }

    float minClusGlobalDistance = trk->auxdata<float>("minClusGlobalDistanceAtLayer");
    if (minClusGlobalDistance==0) { minClusGlobalDistance=100.0; }

    m_minClusGlobalDistanceAtLayer -> Fill(minClusGlobalDistance,weight);
    m_minClusGlobalDistanceAtLayer -> Fill(trk->auxdata<float>("minClusGlobalDistanceAtLayerX"),weight);
    m_minClusGlobalDistanceAtLayer -> Fill(trk->auxdata<float>("minClusGlobalDistanceAtLayerY"),weight);
    if (clusPatternIndex== 0) { m_minClusGlobalDistanceAtLayer0  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 1) { m_minClusGlobalDistanceAtLayer1  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 2) { m_minClusGlobalDistanceAtLayer2  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 3) { m_minClusGlobalDistanceAtLayer3  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 4) { m_minClusGlobalDistanceAtLayer4  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 5) { m_minClusGlobalDistanceAtLayer5  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 6) { m_minClusGlobalDistanceAtLayer6  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 7) { m_minClusGlobalDistanceAtLayer7  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 8) { m_minClusGlobalDistanceAtLayer8  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex== 9) { m_minClusGlobalDistanceAtLayer9  -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==10) { m_minClusGlobalDistanceAtLayer10 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==11) { m_minClusGlobalDistanceAtLayer11 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==12) { m_minClusGlobalDistanceAtLayer12 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==13) { m_minClusGlobalDistanceAtLayer13 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==14) { m_minClusGlobalDistanceAtLayer14 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==15) { m_minClusGlobalDistanceAtLayer15 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==16) { m_minClusGlobalDistanceAtLayer16 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==17) { m_minClusGlobalDistanceAtLayer17 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==18) { m_minClusGlobalDistanceAtLayer18 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==19) { m_minClusGlobalDistanceAtLayer19 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==20) { m_minClusGlobalDistanceAtLayer20 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==21) { m_minClusGlobalDistanceAtLayer21 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==22) { m_minClusGlobalDistanceAtLayer22 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==23) { m_minClusGlobalDistanceAtLayer23 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==24) { m_minClusGlobalDistanceAtLayer24 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==25) { m_minClusGlobalDistanceAtLayer25 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==26) { m_minClusGlobalDistanceAtLayer26 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==27) { m_minClusGlobalDistanceAtLayer27 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==28) { m_minClusGlobalDistanceAtLayer28 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==29) { m_minClusGlobalDistanceAtLayer29 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==30) { m_minClusGlobalDistanceAtLayer30 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==31) { m_minClusGlobalDistanceAtLayer31 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==32) { m_minClusGlobalDistanceAtLayer32 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==33) { m_minClusGlobalDistanceAtLayer33 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==34) { m_minClusGlobalDistanceAtLayer34 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==35) { m_minClusGlobalDistanceAtLayer35 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==36) { m_minClusGlobalDistanceAtLayer36 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==37) { m_minClusGlobalDistanceAtLayer37 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==38) { m_minClusGlobalDistanceAtLayer38 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==39) { m_minClusGlobalDistanceAtLayer39 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==40) { m_minClusGlobalDistanceAtLayer40 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==41) { m_minClusGlobalDistanceAtLayer41 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==42) { m_minClusGlobalDistanceAtLayer42 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==43) { m_minClusGlobalDistanceAtLayer43 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==44) { m_minClusGlobalDistanceAtLayer44 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==45) { m_minClusGlobalDistanceAtLayer45 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==46) { m_minClusGlobalDistanceAtLayer46 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==47) { m_minClusGlobalDistanceAtLayer47 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==48) { m_minClusGlobalDistanceAtLayer48 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==49) { m_minClusGlobalDistanceAtLayer49 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==50) { m_minClusGlobalDistanceAtLayer50 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==51) { m_minClusGlobalDistanceAtLayer51 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==52) { m_minClusGlobalDistanceAtLayer52 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==53) { m_minClusGlobalDistanceAtLayer53 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==54) { m_minClusGlobalDistanceAtLayer54 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==55) { m_minClusGlobalDistanceAtLayer55 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==56) { m_minClusGlobalDistanceAtLayer56 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==57) { m_minClusGlobalDistanceAtLayer57 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==58) { m_minClusGlobalDistanceAtLayer58 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==59) { m_minClusGlobalDistanceAtLayer59 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==60) { m_minClusGlobalDistanceAtLayer60 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==61) { m_minClusGlobalDistanceAtLayer61 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==62) { m_minClusGlobalDistanceAtLayer62 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==63) { m_minClusGlobalDistanceAtLayer63 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==64) { m_minClusGlobalDistanceAtLayer64 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==65) { m_minClusGlobalDistanceAtLayer65 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==66) { m_minClusGlobalDistanceAtLayer66 -> Fill(minClusLocalDistance,weight); }
    if (clusPatternIndex==67) { m_minClusGlobalDistanceAtLayer67 -> Fill(minClusLocalDistance,weight); }

    int nRefLayer = trk->auxdata<int>("nRefHitAtLayerOfInterest");
    const xAOD::TrackStateValidation *msosRef1 = trk->auxdata<const xAOD::TrackStateValidation *>("RefMSOS1");
    const xAOD::TrackStateValidation *msosRef2 = trk->auxdata<const xAOD::TrackStateValidation *>("RefMSOS2");
    const xAOD::TrackMeasurementValidation *pixRefClus1 = trk->auxdata<const xAOD::TrackMeasurementValidation *>("PixelRefClus1");
    const xAOD::TrackMeasurementValidation *pixRefClus2 = trk->auxdata<const xAOD::TrackMeasurementValidation *>("PixelRefClus2");

    int nTarLayer = trk->auxdata<int>("nTarHitAtLayerOfInterest");
    const xAOD::TrackStateValidation *msosTar1 = trk->auxdata<const xAOD::TrackStateValidation *>("TarMSOS1");
    const xAOD::TrackStateValidation *msosTar2 = trk->auxdata<const xAOD::TrackStateValidation *>("TarMSOS2");
    const xAOD::TrackMeasurementValidation *pixTarClus1 = trk->auxdata<const xAOD::TrackMeasurementValidation *>("PixelTarClus1");
    const xAOD::TrackMeasurementValidation *pixTarClus2 = trk->auxdata<const xAOD::TrackMeasurementValidation *>("PixelTarClus2");

    if (nRefLayer>0) {
      m_pixLocal_shiftX -> Fill(pixRefClus1->localX()-msosRef1->localX(),weight);
      m_pixLocal_shiftY -> Fill(pixRefClus1->localY()-msosRef1->localY(),weight);
      if (fabs(pixRefClus1->globalZ())<450.0) {
        m_PixelLayout_xy     -> Fill(pixRefClus1->globalX(),pixRefClus1->globalY(),weight);
        m_PixelLayoutZoom_xy -> Fill(pixRefClus1->globalX(),pixRefClus1->globalY(),weight);
      }
      m_PixelLayout_zr     -> Fill(pixRefClus1->globalZ(),sqrt(pixRefClus1->globalX()*pixRefClus1->globalX()+pixRefClus1->globalY()*pixRefClus1->globalY()),weight);
      m_PixelLayoutZoom_zr -> Fill(pixRefClus1->globalZ(),sqrt(pixRefClus1->globalX()*pixRefClus1->globalX()+pixRefClus1->globalY()*pixRefClus1->globalY()),weight);
    }
    if (nRefLayer>1) {
      m_pixLocal_shiftX -> Fill(pixRefClus2->localX()-msosRef2->localX(),weight);
      m_pixLocal_shiftY -> Fill(pixRefClus2->localY()-msosRef2->localY(),weight);
      if (fabs(pixRefClus2->globalZ())<450.0) {
        m_PixelLayout_xy          -> Fill(pixRefClus2->globalX(),pixRefClus2->globalY(),weight);
        m_PixelLayout2Hits_xy     -> Fill(pixRefClus1->globalX(),pixRefClus1->globalY(),weight);
        m_PixelLayout2Hits_xy     -> Fill(pixRefClus2->globalX(),pixRefClus2->globalY(),weight);
        m_PixelLayoutZoom_xy      -> Fill(pixRefClus2->globalX(),pixRefClus2->globalY(),weight);
        m_PixelLayout2HitsZoom_xy -> Fill(pixRefClus1->globalX(),pixRefClus1->globalY(),weight);
        m_PixelLayout2HitsZoom_xy -> Fill(pixRefClus2->globalX(),pixRefClus2->globalY(),weight);
      }
      m_PixelLayout_zr          -> Fill(pixRefClus2->globalZ(),sqrt(pixRefClus2->globalX()*pixRefClus2->globalX()+pixRefClus2->globalY()*pixRefClus2->globalY()),weight);
      m_PixelLayout2Hits_zr     -> Fill(pixRefClus1->globalZ(),sqrt(pixRefClus1->globalX()*pixRefClus1->globalX()+pixRefClus1->globalY()*pixRefClus1->globalY()),weight);
      m_PixelLayout2Hits_zr     -> Fill(pixRefClus2->globalZ(),sqrt(pixRefClus2->globalX()*pixRefClus2->globalX()+pixRefClus2->globalY()*pixRefClus2->globalY()),weight);
      m_PixelLayoutZoom_zr      -> Fill(pixRefClus2->globalZ(),sqrt(pixRefClus2->globalX()*pixRefClus2->globalX()+pixRefClus2->globalY()*pixRefClus2->globalY()),weight);
      m_PixelLayout2HitsZoom_zr -> Fill(pixRefClus1->globalZ(),sqrt(pixRefClus1->globalX()*pixRefClus1->globalX()+pixRefClus1->globalY()*pixRefClus1->globalY()),weight);
      m_PixelLayout2HitsZoom_zr -> Fill(pixRefClus2->globalZ(),sqrt(pixRefClus2->globalX()*pixRefClus2->globalX()+pixRefClus2->globalY()*pixRefClus2->globalY()),weight);
    }

    if (nTarLayer>0) {
      m_pixLocal_shiftX -> Fill(pixTarClus1->localX()-msosTar1->localX(),weight);
      m_pixLocal_shiftY -> Fill(pixTarClus1->localY()-msosTar1->localY(),weight);
      if (fabs(pixTarClus1->globalZ())<450.0) {
        m_PixelLayout_xy     -> Fill(pixTarClus1->globalX(),pixTarClus1->globalY(),weight);
        m_PixelLayoutZoom_xy -> Fill(pixTarClus1->globalX(),pixTarClus1->globalY(),weight);
      }
      m_PixelLayout_zr     -> Fill(pixTarClus1->globalZ(),sqrt(pixTarClus1->globalX()*pixTarClus1->globalX()+pixTarClus1->globalY()*pixTarClus1->globalY()),weight);
      m_PixelLayoutZoom_zr -> Fill(pixTarClus1->globalZ(),sqrt(pixTarClus1->globalX()*pixTarClus1->globalX()+pixTarClus1->globalY()*pixTarClus1->globalY()),weight);
    }
    if (nTarLayer>1) {
      m_pixLocal_shiftX -> Fill(pixTarClus2->localX()-msosTar2->localX(),weight);
      m_pixLocal_shiftY -> Fill(pixTarClus2->localY()-msosTar2->localY(),weight);
      if (fabs(pixTarClus2->globalZ())<450.0) {
        m_PixelLayout_xy          -> Fill(pixTarClus2->globalX(),pixTarClus2->globalY(),weight);
        m_PixelLayout2Hits_xy     -> Fill(pixTarClus1->globalX(),pixTarClus1->globalY(),weight);
        m_PixelLayout2Hits_xy     -> Fill(pixTarClus2->globalX(),pixTarClus2->globalY(),weight);
        m_PixelLayoutZoom_xy      -> Fill(pixTarClus2->globalX(),pixTarClus2->globalY(),weight);
        m_PixelLayout2HitsZoom_xy -> Fill(pixTarClus1->globalX(),pixTarClus1->globalY(),weight);
        m_PixelLayout2HitsZoom_xy -> Fill(pixTarClus2->globalX(),pixTarClus2->globalY(),weight);
      }
      m_PixelLayout_zr          -> Fill(pixTarClus2->globalZ(),sqrt(pixTarClus2->globalX()*pixTarClus2->globalX()+pixTarClus2->globalY()*pixTarClus2->globalY()),weight);
      m_PixelLayout2Hits_zr     -> Fill(pixTarClus1->globalZ(),sqrt(pixTarClus1->globalX()*pixTarClus1->globalX()+pixTarClus1->globalY()*pixTarClus1->globalY()),weight);
      m_PixelLayout2Hits_zr     -> Fill(pixTarClus2->globalZ(),sqrt(pixTarClus2->globalX()*pixTarClus2->globalX()+pixTarClus2->globalY()*pixTarClus2->globalY()),weight);
      m_PixelLayoutZoom_zr      -> Fill(pixTarClus2->globalZ(),sqrt(pixTarClus2->globalX()*pixTarClus2->globalX()+pixTarClus2->globalY()*pixTarClus2->globalY()),weight);
      m_PixelLayout2HitsZoom_zr -> Fill(pixTarClus1->globalZ(),sqrt(pixTarClus1->globalX()*pixTarClus1->globalX()+pixTarClus1->globalY()*pixTarClus1->globalY()),weight);
      m_PixelLayout2HitsZoom_zr -> Fill(pixTarClus2->globalZ(),sqrt(pixTarClus2->globalX()*pixTarClus2->globalX()+pixTarClus2->globalY()*pixTarClus2->globalY()),weight);
    }


    xAOD::TrackParticleContainer::const_iterator subtrk = trk->auxdata<xAOD::TrackParticleContainer::const_iterator>("twoRecoSiblingTrkIter");

    m_etaWeight -> Fill(trk->eta(),(*subtrk)->eta(),weight);
    m_ptWeight  -> Fill(trk->pt()/1000.0,(*subtrk)->pt()/1000.0,weight);


    // Correlation with jet
    float mnDRJ1 = trk->auxdata<float>("mindRJet1");
    float mnDRJ2 = trk->auxdata<float>("mindRJet2");
    float ptJ1   = trk->auxdata<float>("pTJet1");
    float ptJ2   = trk->auxdata<float>("pTJet2");


    if (TMath::Abs(ptJ1-ptJ2)<0.01 && mnDRJ1<0.4 && mnDRJ2<0.4) {
      m_mindRJet1 -> Fill(mnDRJ1,weight);
      if (clusPatternIndex== 0) { m_mindRJet1P0  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 1) { m_mindRJet1P1  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 2) { m_mindRJet1P2  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 3) { m_mindRJet1P3  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 4) { m_mindRJet1P4  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 5) { m_mindRJet1P5  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 6) { m_mindRJet1P6  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 7) { m_mindRJet1P7  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 8) { m_mindRJet1P8  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex== 9) { m_mindRJet1P9  -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==10) { m_mindRJet1P10 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==11) { m_mindRJet1P11 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==12) { m_mindRJet1P12 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==13) { m_mindRJet1P13 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==14) { m_mindRJet1P14 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==15) { m_mindRJet1P15 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==16) { m_mindRJet1P16 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==17) { m_mindRJet1P17 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==18) { m_mindRJet1P18 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==19) { m_mindRJet1P19 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==20) { m_mindRJet1P20 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==21) { m_mindRJet1P21 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==22) { m_mindRJet1P22 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==23) { m_mindRJet1P23 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==24) { m_mindRJet1P24 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==25) { m_mindRJet1P25 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==26) { m_mindRJet1P26 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==27) { m_mindRJet1P27 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==28) { m_mindRJet1P28 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==29) { m_mindRJet1P29 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==30) { m_mindRJet1P30 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==31) { m_mindRJet1P31 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==32) { m_mindRJet1P32 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==33) { m_mindRJet1P33 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==34) { m_mindRJet1P34 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==35) { m_mindRJet1P35 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==36) { m_mindRJet1P36 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==37) { m_mindRJet1P37 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==38) { m_mindRJet1P38 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==39) { m_mindRJet1P39 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==40) { m_mindRJet1P40 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==41) { m_mindRJet1P41 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==42) { m_mindRJet1P42 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==43) { m_mindRJet1P43 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==44) { m_mindRJet1P44 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==45) { m_mindRJet1P45 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==46) { m_mindRJet1P46 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==47) { m_mindRJet1P47 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==48) { m_mindRJet1P48 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==49) { m_mindRJet1P49 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==50) { m_mindRJet1P50 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==51) { m_mindRJet1P51 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==52) { m_mindRJet1P52 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==53) { m_mindRJet1P53 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==54) { m_mindRJet1P54 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==55) { m_mindRJet1P55 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==56) { m_mindRJet1P56 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==57) { m_mindRJet1P57 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==58) { m_mindRJet1P58 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==59) { m_mindRJet1P59 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==60) { m_mindRJet1P60 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==61) { m_mindRJet1P61 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==62) { m_mindRJet1P62 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==63) { m_mindRJet1P63 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==64) { m_mindRJet1P64 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==65) { m_mindRJet1P65 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==66) { m_mindRJet1P66 -> Fill(mnDRJ1,weight); }
      if (clusPatternIndex==67) { m_mindRJet1P67 -> Fill(mnDRJ1,weight); }

      m_mindRJet2 -> Fill(mnDRJ1,weight);
      if (clusPatternIndex== 0) { m_mindRJet2P0  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 1) { m_mindRJet2P1  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 2) { m_mindRJet2P2  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 3) { m_mindRJet2P3  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 4) { m_mindRJet2P4  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 5) { m_mindRJet2P5  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 6) { m_mindRJet2P6  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 7) { m_mindRJet2P7  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 8) { m_mindRJet2P8  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex== 9) { m_mindRJet2P9  -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==10) { m_mindRJet2P10 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==11) { m_mindRJet2P11 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==12) { m_mindRJet2P12 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==13) { m_mindRJet2P13 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==14) { m_mindRJet2P14 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==15) { m_mindRJet2P15 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==16) { m_mindRJet2P16 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==17) { m_mindRJet2P17 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==18) { m_mindRJet2P18 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==19) { m_mindRJet2P19 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==20) { m_mindRJet2P20 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==21) { m_mindRJet2P21 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==22) { m_mindRJet2P22 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==23) { m_mindRJet2P23 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==24) { m_mindRJet2P24 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==25) { m_mindRJet2P25 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==26) { m_mindRJet2P26 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==27) { m_mindRJet2P27 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==28) { m_mindRJet2P28 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==29) { m_mindRJet2P29 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==30) { m_mindRJet2P30 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==31) { m_mindRJet2P31 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==32) { m_mindRJet2P32 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==33) { m_mindRJet2P33 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==34) { m_mindRJet2P34 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==35) { m_mindRJet2P35 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==36) { m_mindRJet2P36 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==37) { m_mindRJet2P37 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==38) { m_mindRJet2P38 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==39) { m_mindRJet2P39 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==40) { m_mindRJet2P40 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==41) { m_mindRJet2P41 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==42) { m_mindRJet2P42 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==43) { m_mindRJet2P43 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==44) { m_mindRJet2P44 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==45) { m_mindRJet2P45 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==46) { m_mindRJet2P46 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==47) { m_mindRJet2P47 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==48) { m_mindRJet2P48 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==49) { m_mindRJet2P49 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==50) { m_mindRJet2P50 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==51) { m_mindRJet2P51 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==52) { m_mindRJet2P52 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==53) { m_mindRJet2P53 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==54) { m_mindRJet2P54 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==55) { m_mindRJet2P55 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==56) { m_mindRJet2P56 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==57) { m_mindRJet2P57 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==58) { m_mindRJet2P58 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==59) { m_mindRJet2P59 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==60) { m_mindRJet2P60 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==61) { m_mindRJet2P61 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==62) { m_mindRJet2P62 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==63) { m_mindRJet2P63 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==64) { m_mindRJet2P64 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==65) { m_mindRJet2P65 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==66) { m_mindRJet2P66 -> Fill(mnDRJ2,weight); }
      if (clusPatternIndex==67) { m_mindRJet2P67 -> Fill(mnDRJ2,weight); }

      m_jetpT_locT -> Fill(ptJ1/1000.0,minMSOSLocalDistance,weight);
      m_jetpT_locX -> Fill(ptJ1/1000.0,trk->auxdata<float>("minClusLocalDistanceAtLayerX"),weight);
      m_jetpT_locY -> Fill(ptJ1/1000.0,trk->auxdata<float>("minClusLocalDistanceAtLayerY"),weight);

      if (clusPatternIndex==15 || clusPatternIndex==42 || clusPatternIndex==58) {
        m_jetpT_locTNP -> Fill(ptJ1/1000.0,minMSOSLocalDistance,weight);
        m_jetpT_locXNP -> Fill(ptJ1/1000.0,trk->auxdata<float>("minClusLocalDistanceAtLayerX"),weight);
        m_jetpT_locYNP -> Fill(ptJ1/1000.0,trk->auxdata<float>("minClusLocalDistanceAtLayerY"),weight);
      }
      if (clusPatternIndex==39 || clusPatternIndex==54 || clusPatternIndex==56 || clusPatternIndex==64 || clusPatternIndex==65 || clusPatternIndex==66 || clusPatternIndex==67) {
        m_jetpT_locTSP -> Fill(ptJ1/1000.0,minMSOSLocalDistance,weight);
        m_jetpT_locXSP -> Fill(ptJ1/1000.0,trk->auxdata<float>("minClusLocalDistanceAtLayerX"),weight);
        m_jetpT_locYSP -> Fill(ptJ1/1000.0,trk->auxdata<float>("minClusLocalDistanceAtLayerY"),weight);
      }

    }

  }

}

