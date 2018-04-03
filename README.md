# EGTagAndProbe
Set of tools to evaluate trigger performance on T&amp;P

Based on TauTagAndProbe package developed by L. Cadamuro & O. Davignon

### Install instructions
```
cmsrel CMSSW_10_0_0_pre3
cd CMSSW_10_0_0_pre3/src
cmsenv
# MVA EleID Fall 2017
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
git clone https://github.com/tstreble/EGTagAndProbe
cd EGTagAndProbe
git checkout HLT_SingleTau
cd -
scram b -j4
```

### Producing TagAndProbe ntuples
Set flag isMC and isMINIAOD according to sample in test/test.py

HLT path used specified in python/MCAnalysis_cff.py (MC) or python/tagAndProbe_cff.py (data)

Launch test.py

Standard Z->ee TagAndProbe selections are applied

Information related to reconstructed taus matched to probe electrons are also saved

Basic trigger information is saved as well:
- eleProbeTriggerBits encodes the trigger decision (including trigger filter matching) stored in a bitwise fashion following the order specified in HLTLISTPROBED in python/XXX_cff.py: the decision of the 2nd algorithm can for instance be checked with (eleProbeTriggerBits>>1)&1
- basic kinematic quantities of the matched filter object are stored as well (pt, eta, phi) as vectors in the same order


### Submit job on the Grid
Modify crab3_config.py: change requestName, inputDataSet, outLFNDirBase, outputDatasetTag, storageSite
```
cd CMSSW_10_0_0_pre3/src/EGTagAndProbe/EGTagAndProbe/test
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms
crab submit -c crab3_config.py
```
