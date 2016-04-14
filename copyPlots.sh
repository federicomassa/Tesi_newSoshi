#!/bin/bash
cd /afs/cern.ch/user/f/fmassa/public/
rm -rf Plots_ptCut3GeV_withPu200.tar
cp -r /afs/cern.ch/user/f/fmassa/newSoshi/Plots/* Plots_ptCut3GeV_withPu200/
tar -cvf Plots_ptCut3GeV_withPu200.tar Plots_ptCut3GeV_withPu200/
cd /afs/cern.ch/user/f/fmassa/newSoshi