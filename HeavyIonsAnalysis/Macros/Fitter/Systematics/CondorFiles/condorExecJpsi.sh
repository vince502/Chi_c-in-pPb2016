#!/bin/bash
cd /eos/cms/store/group/phys_heavyions/okukral/Chi_c/Workspace/CMSSW_8_0_30/src/HeavyIonsAnalysis/Macros/Fitter/Systematics/
cmsenv
root -l -b -q 'SystFits_Jpsi.C('$1', '$2')' 
