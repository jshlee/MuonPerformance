# MuonPerformance
```
scram p -n muonPerf CMSSW CMSSW_10_0_0_pre1
cd muonPerf/src
cmsenv
git cms-init -q
git cms-merge-topic jshlee:muonIso-10x
git clone git@github.com:jshlee/MuonPerformance.git
scram b -j 20
```
