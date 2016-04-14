#!/bin/sh


echo "Are you sure you want to run on all the samples? It could take some time..."
select yn in "Yes" "No"; do
    case $yn in
        Yes ) 
	    ITkRun AnalysisResults/ExtBrl4/datasetPi_15GeV /afs/cern.ch/work/f/fmassa/public/newGeo/out/ExtBrl4/datasetPi_15GeV 211 15000;
	    ITkRun AnalysisResults/IExtBrl4/datasetPi_15GeV /afs/cern.ch/work/f/fmassa/public/newGeo/out/IExtBrl4/datasetPi_15GeV 211 15000;
	    ITkRun AnalysisResults/InclBrl4/datasetPi_15GeV /afs/cern.ch/work/f/fmassa/public/newGeo/out/InclBrl4/datasetPi_15GeV 211 15000;
	    break;;
        No ) exit;;
    esac
done