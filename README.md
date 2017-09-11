# MuonPerformance
```
scram p -n isolation CMSSW CMSSW_9_3_0_pre4
cd isolation/src
cmsenv
git clone git@github.com:jshlee/MuonPerformance.git
cmsenv
git cms-init -q
git cms-merge-topic jshlee:muonPuppiIso
git clone git@github.com:jshlee/MuonPerformance.git
scram b -j 20
```
