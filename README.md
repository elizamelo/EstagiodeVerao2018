# EstagiodeVerao2018
# ggNtuplesAnalyzer
ggNtuplesAnalyzer @ Farm/Fermis

### Setup
```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git clone https://github.com/elizamelo/EstagiodeVerao2018.git
cd EstagiodeVerao2018/ggNtuplesAnalyzer/
```

### Analysis Workflow
```
./run_ana_ggNtuplesAnalyzer.sh
```

This will select dimuon pairs:
- event pass **HLT_IsoMu24** HLT Path
- leading **(tight)** muon: **pt > 27.0** and **|eta| < 2.4**
- trailing **(tight)** muon: **pt > 2.0** and **|eta| < 2.4**
- muons with **opposite charge**

Histograms files can be found at:
```
ls -lha outputFiles
```

### Plotting (i.e. dimuon mass distribution):
Load root:
```
root -l outputFiles/histos_ggNtuplesAnalyzer.root
```

Inside root:
```
((TH1D*)_file0->Get("h_DiMuon_Mass"))->Draw()
```
