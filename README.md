#Ponia-OniaPhoton

This package is mean to be run after the BPH CompactSkim. 

* Setup: (it has being tested on 7_4x, 7_6x, 8_0x but should run in any of the recent cmssw releases)

```
export SCRAM_ARCH=slc6_amd64_gcc530
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_8_0_24_patch1
cd CMSSW_8_0_24_patch1/src/
cmsenv
git clone https://github.com/alberto-sanchez/Ponia-OniaPhoton.git Ponia/OniaPhoton
scram b
```

* Run: (use your favorite input sample)

```
cmsRun Ponia/OniaPhoton/test/runChiConAODSIM.py (for chic reconstruction using AODSIM)
```

In test directory you can find other examples to run over data in AOD, and BPHSKIM, as well as for MINIAOD,
however this latest format is not available at the moment. 

