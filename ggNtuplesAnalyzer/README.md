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



## ggnTuples Content
### Data
```
******************************************************************************
*Tree    :EventTree : Event data (tag V08_00_26_06)                          *
*Entries :    12429 : Total =        80705538 bytes  File  Size =   26592078 *
*        :          : Tree compression factor =   3.03                       *
******************************************************************************
*Br    0 :run       : run/I                                                  *
*Entries :    12429 : Total  Size=      50332 bytes  File Size  =        478 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression= 104.32     *
*............................................................................*
*Br    1 :event     : event/L                                                *
*Entries :    12429 : Total  Size=     100222 bytes  File Size  =      49699 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.01     *
*............................................................................*
*Br    2 :lumis     : lumis/I                                                *
*Entries :    12429 : Total  Size=      50344 bytes  File Size  =        653 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=  76.37     *
*............................................................................*
*Br    3 :isData    : isData/O                                               *
*Entries :    12429 : Total  Size=      12979 bytes  File Size  =        166 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  75.34     *
*............................................................................*
*Br    4 :nVtx      : nVtx/I                                                 *
*Entries :    12429 : Total  Size=      50338 bytes  File Size  =      11660 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.28     *
*............................................................................*
*Br    5 :nGoodVtx  : nGoodVtx/I                                             *
*Entries :    12429 : Total  Size=      50362 bytes  File Size  =      11559 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.31     *
*............................................................................*
*Br    6 :nTrksPV   : nTrksPV/I                                              *
*Entries :    12429 : Total  Size=      50356 bytes  File Size  =        439 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression= 113.61     *
*............................................................................*
*Br    7 :isPVGood  : isPVGood/O                                             *
*Entries :    12429 : Total  Size=      12989 bytes  File Size  =        229 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  54.62     *
*............................................................................*
*Br    8 :vtx       : vtx/F                                                  *
*Entries :    12429 : Total  Size=      50332 bytes  File Size  =      40488 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.23     *
*............................................................................*
*Br    9 :vty       : vty/F                                                  *
*Entries :    12429 : Total  Size=      50332 bytes  File Size  =      40152 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.24     *
*............................................................................*
*Br   10 :vtz       : vtz/F                                                  *
*Entries :    12429 : Total  Size=      50332 bytes  File Size  =      46315 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   11 :rho       : rho/F                                                  *
*Entries :    12429 : Total  Size=      50332 bytes  File Size  =      44283 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   12 :rhoCentral : rhoCentral/F                                          *
*Entries :    12429 : Total  Size=      50374 bytes  File Size  =      43104 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.16     *
*............................................................................*
*Br   13 :HLTEleMuX : HLTEleMuX/l                                            *
*Entries :    12429 : Total  Size=     100254 bytes  File Size  =      16247 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   6.14     *
*............................................................................*
*Br   14 :HLTPho    : HLTPho/l                                               *
*Entries :    12429 : Total  Size=     100230 bytes  File Size  =       1790 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=  55.72     *
*............................................................................*
*Br   15 :HLTJet    : HLTJet/l                                               *
*Entries :    12429 : Total  Size=     100230 bytes  File Size  =       2279 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=  43.77     *
*............................................................................*
*Br   16 :HLTEleMuXIsPrescaled : HLTEleMuXIsPrescaled/l                      *
*Entries :    12429 : Total  Size=     100342 bytes  File Size  =       1066 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=  93.62     *
*............................................................................*
*Br   17 :HLTPhoIsPrescaled : HLTPhoIsPrescaled/l                            *
*Entries :    12429 : Total  Size=     100318 bytes  File Size  =       1054 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=  94.68     *
*............................................................................*
*Br   18 :HLTJetIsPrescaled : HLTJetIsPrescaled/l                            *
*Entries :    12429 : Total  Size=     100318 bytes  File Size  =       1062 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=  93.96     *
*............................................................................*
*Br   19 :phoPrescale : vector<int>                                          *
*Entries :    12429 : Total  Size=     624214 bytes  File Size  =      30211 *
*Baskets :       22 : Basket Size=      32000 bytes  Compression=  20.64     *
*............................................................................*
*Br   20 :metFilters : metFilters/I                                          *
*Entries :    12429 : Total  Size=      50374 bytes  File Size  =       1099 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=  45.39     *
*............................................................................*
*Br   21 :pfMET     : pfMET/F                                                *
*Entries :    12429 : Total  Size=      50344 bytes  File Size  =      44803 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   22 :pfMETPhi  : pfMETPhi/F                                             *
*Entries :    12429 : Total  Size=      50362 bytes  File Size  =      46354 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   23 :pfMETsumEt : pfMETsumEt/F                                          *
*Entries :    12429 : Total  Size=      50374 bytes  File Size  =      43568 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.14     *
*............................................................................*
*Br   24 :pfMETmEtSig : pfMETmEtSig/F                                        *
*Entries :    12429 : Total  Size=      50380 bytes  File Size  =      44616 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   25 :pfMETSig  : pfMETSig/F                                             *
*Entries :    12429 : Total  Size=      50362 bytes  File Size  =      45722 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.09     *
*............................................................................*
*Br   26 :pfMET_T1JERUp : pfMET_T1JERUp/F                                    *
*Entries :    12429 : Total  Size=      50392 bytes  File Size  =      44818 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   27 :pfMET_T1JERDo : pfMET_T1JERDo/F                                    *
*Entries :    12429 : Total  Size=      50392 bytes  File Size  =      44818 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   28 :pfMET_T1JESUp : pfMET_T1JESUp/F                                    *
*Entries :    12429 : Total  Size=      50392 bytes  File Size  =      44846 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   29 :pfMET_T1JESDo : pfMET_T1JESDo/F                                    *
*Entries :    12429 : Total  Size=      50392 bytes  File Size  =      44822 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   30 :pfMET_T1UESUp : pfMET_T1UESUp/F                                    *
*Entries :    12429 : Total  Size=      50392 bytes  File Size  =      44800 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   31 :pfMET_T1UESDo : pfMET_T1UESDo/F                                    *
*Entries :    12429 : Total  Size=      50392 bytes  File Size  =      44797 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   32 :pfMETPhi_T1JESUp : pfMETPhi_T1JESUp/F                              *
*Entries :    12429 : Total  Size=      50410 bytes  File Size  =      46353 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   33 :pfMETPhi_T1JESDo : pfMETPhi_T1JESDo/F                              *
*Entries :    12429 : Total  Size=      50410 bytes  File Size  =      46379 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   34 :pfMETPhi_T1UESUp : pfMETPhi_T1UESUp/F                              *
*Entries :    12429 : Total  Size=      50410 bytes  File Size  =      46353 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   35 :pfMETPhi_T1UESDo : pfMETPhi_T1UESDo/F                              *
*Entries :    12429 : Total  Size=      50410 bytes  File Size  =      46308 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   36 :nPho      : nPho/I                                                 *
*Entries :    12429 : Total  Size=      50338 bytes  File Size  =       6322 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   7.89     *
*............................................................................*
*Br   37 :phoE      : vector<float>                                          *
*Entries :    12429 : Total  Size=     241752 bytes  File Size  =     111304 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.17     *
*............................................................................*
*Br   38 :phoEt     : vector<float>                                          *
*Entries :    12429 : Total  Size=     241766 bytes  File Size  =     110672 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.18     *
*............................................................................*
*Br   39 :phoEta    : vector<float>                                          *
*Entries :    12429 : Total  Size=     241780 bytes  File Size  =     114353 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.11     *
*............................................................................*
*Br   40 :phoPhi    : vector<float>                                          *
*Entries :    12429 : Total  Size=     241780 bytes  File Size  =     114585 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.11     *
*............................................................................*
*Br   41 :phoCalibE : vector<float>                                          *
*Entries :    12429 : Total  Size=     241822 bytes  File Size  =     111541 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.16     *
*............................................................................*
*Br   42 :phoCalibEt : vector<float>                                         *
*Entries :    12429 : Total  Size=     241836 bytes  File Size  =     110628 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.18     *
*............................................................................*
*Br   43 :phoSCE    : vector<float>                                          *
*Entries :    12429 : Total  Size=     241780 bytes  File Size  =     111256 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.17     *
*............................................................................*
*Br   44 :phoSCRawE : vector<float>                                          *
*Entries :    12429 : Total  Size=     241822 bytes  File Size  =     111443 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.17     *
*............................................................................*
*Br   45 :phoESEnP1 : vector<float>                                          *
*Entries :    12429 : Total  Size=     241822 bytes  File Size  =      79792 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.02     *
*............................................................................*
*Br   46 :phoESEnP2 : vector<float>                                          *
*Entries :    12429 : Total  Size=     241822 bytes  File Size  =      79902 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.02     *
*............................................................................*
*Br   47 :phoSCEta  : vector<float>                                          *
*Entries :    12429 : Total  Size=     241808 bytes  File Size  =     114480 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.11     *
*............................................................................*
*Br   48 :phoSCPhi  : vector<float>                                          *
*Entries :    12429 : Total  Size=     241808 bytes  File Size  =     114607 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.11     *
*............................................................................*
*Br   49 :phoSCEtaWidth : vector<float>                                      *
*Entries :    12429 : Total  Size=     241878 bytes  File Size  =     109318 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.21     *
*............................................................................*
*Br   50 :phoSCPhiWidth : vector<float>                                      *
*Entries :    12429 : Total  Size=     241878 bytes  File Size  =     110452 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br   51 :phoSCBrem : vector<float>                                          *
*Entries :    12429 : Total  Size=     241822 bytes  File Size  =     111117 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.17     *
*............................................................................*
*Br   52 :phohasPixelSeed : vector<int>                                      *
*Entries :    12429 : Total  Size=     241906 bytes  File Size  =      49830 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.84     *
*............................................................................*
*Br   53 :phoEleVeto : vector<int>                                           *
*Entries :    12429 : Total  Size=     241836 bytes  File Size  =      48880 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.94     *
*............................................................................*
*Br   54 :phoR9     : vector<float>                                          *
*Entries :    12429 : Total  Size=     241766 bytes  File Size  =     108450 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.22     *
*............................................................................*
*Br   55 :phoHoverE : vector<float>                                          *
*Entries :    12429 : Total  Size=     241822 bytes  File Size  =      71859 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.36     *
*............................................................................*
*Br   56 :phoE1x3   : vector<float>                                          *
*Entries :    12429 : Total  Size=     241794 bytes  File Size  =     110206 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br   57 :phoE1x5   : vector<float>                                          *
*Entries :    12429 : Total  Size=     241794 bytes  File Size  =     112102 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.15     *
*............................................................................*
*Br   58 :phoE2x2   : vector<float>                                          *
*Entries :    12429 : Total  Size=     241794 bytes  File Size  =     110160 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br   59 :phoE2x5Max : vector<float>                                         *
*Entries :    12429 : Total  Size=     241836 bytes  File Size  =     111746 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.16     *
*............................................................................*
*Br   60 :phoE5x5   : vector<float>                                          *
*Entries :    12429 : Total  Size=     241794 bytes  File Size  =     111938 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.16     *
*............................................................................*
*Br   61 :phoESEffSigmaRR : vector<float>                                    *
*Entries :    12429 : Total  Size=     241906 bytes  File Size  =      78497 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.07     *
*............................................................................*
*Br   62 :phoSigmaIEtaIEtaFull5x5 : vector<float>                            *
*Entries :    12429 : Total  Size=     242018 bytes  File Size  =     109178 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.21     *
*............................................................................*
*Br   63 :phoSigmaIEtaIPhiFull5x5 : vector<float>                            *
*Entries :    12429 : Total  Size=     242018 bytes  File Size  =     113617 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.13     *
*............................................................................*
*Br   64 :phoSigmaIPhiIPhiFull5x5 : vector<float>                            *
*Entries :    12429 : Total  Size=     242018 bytes  File Size  =     108306 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.23     *
*............................................................................*
*Br   65 :phoE1x3Full5x5 : vector<float>                                     *
*Entries :    12429 : Total  Size=     241892 bytes  File Size  =     110135 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br   66 :phoE1x5Full5x5 : vector<float>                                     *
*Entries :    12429 : Total  Size=     241892 bytes  File Size  =     111845 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.16     *
*............................................................................*
*Br   67 :phoE2x2Full5x5 : vector<float>                                     *
*Entries :    12429 : Total  Size=     241892 bytes  File Size  =     109845 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.20     *
*............................................................................*
*Br   68 :phoE2x5MaxFull5x5 : vector<float>                                  *
*Entries :    12429 : Total  Size=     241934 bytes  File Size  =     111820 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.16     *
*............................................................................*
*Br   69 :phoE5x5Full5x5 : vector<float>                                     *
*Entries :    12429 : Total  Size=     241892 bytes  File Size  =     111808 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.16     *
*............................................................................*
*Br   70 :phoR9Full5x5 : vector<float>                                       *
*Entries :    12429 : Total  Size=     241864 bytes  File Size  =     108011 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.23     *
*............................................................................*
*Br   71 :phoPFChIso : vector<float>                                         *
*Entries :    12429 : Total  Size=     241836 bytes  File Size  =      93479 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.58     *
*............................................................................*
*Br   72 :phoPFPhoIso : vector<float>                                        *
*Entries :    12429 : Total  Size=     241850 bytes  File Size  =      95753 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.52     *
*............................................................................*
*Br   73 :phoPFNeuIso : vector<float>                                        *
*Entries :    12429 : Total  Size=     241850 bytes  File Size  =      72127 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.35     *
*............................................................................*
*Br   74 :phoPFChWorstIso : vector<float>                                    *
*Entries :    12429 : Total  Size=     241906 bytes  File Size  =      98224 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br   75 :phoCITKChIso : vector<float>                                       *
*Entries :    12429 : Total  Size=     241864 bytes  File Size  =      94266 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.56     *
*............................................................................*
*Br   76 :phoCITKPhoIso : vector<float>                                      *
*Entries :    12429 : Total  Size=     241878 bytes  File Size  =      95763 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.52     *
*............................................................................*
*Br   77 :phoCITKNeuIso : vector<float>                                      *
*Entries :    12429 : Total  Size=     241878 bytes  File Size  =      72190 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.34     *
*............................................................................*
*Br   78 :phoIDMVA  : vector<float>                                          *
*Entries :    12429 : Total  Size=     241808 bytes  File Size  =     109911 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.20     *
*............................................................................*
*Br   79 :phoFiredSingleTrgs : vector<unsigned int>                          *
*Entries :    12429 : Total  Size=     241948 bytes  File Size  =      41201 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   5.86     *
*............................................................................*
*Br   80 :phoFiredDoubleTrgs : vector<unsigned int>                          *
*Entries :    12429 : Total  Size=     241948 bytes  File Size  =      60677 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.98     *
*............................................................................*
*Br   81 :phoFiredL1Trgs : vector<unsigned int>                              *
*Entries :    12429 : Total  Size=     241892 bytes  File Size  =      78238 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   3.08     *
*............................................................................*
*Br   82 :phoSeedTime : vector<float>                                        *
*Entries :    12429 : Total  Size=     241850 bytes  File Size  =     111468 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.16     *
*............................................................................*
*Br   83 :phoSeedEnergy : vector<float>                                      *
*Entries :    12429 : Total  Size=     241878 bytes  File Size  =     110068 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br   84 :phoxtalBits : vector<unsigned short>                               *
*Entries :    12429 : Total  Size=     208533 bytes  File Size  =      37754 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   5.51     *
*............................................................................*
*Br   85 :phoIDbit  : vector<unsigned short>                                 *
*Entries :    12429 : Total  Size=     208494 bytes  File Size  =      44817 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.64     *
*............................................................................*
*Br   86 :phoScale_stat_up : vector<float>                                   *
*Entries :    12429 : Total  Size=     241920 bytes  File Size  =      49624 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.86     *
*............................................................................*
*Br   87 :phoScale_stat_dn : vector<float>                                   *
*Entries :    12429 : Total  Size=     241920 bytes  File Size  =      51285 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.71     *
*............................................................................*
*Br   88 :phoScale_syst_up : vector<float>                                   *
*Entries :    12429 : Total  Size=     241920 bytes  File Size  =      59104 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.08     *
*............................................................................*
*Br   89 :phoScale_syst_dn : vector<float>                                   *
*Entries :    12429 : Total  Size=     241920 bytes  File Size  =      59918 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.03     *
*............................................................................*
*Br   90 :phoScale_gain_up : vector<float>                                   *
*Entries :    12429 : Total  Size=     241920 bytes  File Size  =      38194 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   6.32     *
*............................................................................*
*Br   91 :phoScale_gain_dn : vector<float>                                   *
*Entries :    12429 : Total  Size=     241920 bytes  File Size  =      38224 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   6.31     *
*............................................................................*
*Br   92 :phoResol_rho_up : vector<float>                                    *
*Entries :    12429 : Total  Size=     241906 bytes  File Size  =      58869 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.10     *
*............................................................................*
*Br   93 :phoResol_rho_dn : vector<float>                                    *
*Entries :    12429 : Total  Size=     241906 bytes  File Size  =      58876 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   4.10     *
*............................................................................*
*Br   94 :phoResol_phi_up : vector<float>                                    *
*Entries :    12429 : Total  Size=     241906 bytes  File Size  =     110291 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br   95 :phoResol_phi_dn : vector<float>                                    *
*Entries :    12429 : Total  Size=     241906 bytes  File Size  =     110291 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br   96 :npfPho    : npfPho/I                                               *
*Entries :    12429 : Total  Size=      50350 bytes  File Size  =        437 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression= 114.12     *
*............................................................................*
*Br   97 :pfphoEt   : vector<float>                                          *
*Entries :    12429 : Total  Size=     175188 bytes  File Size  =      19543 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br   98 :pfphoEta  : vector<float>                                          *
*Entries :    12429 : Total  Size=     175200 bytes  File Size  =      19516 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br   99 :pfphoPhi  : vector<float>                                          *
*Entries :    12429 : Total  Size=     175200 bytes  File Size  =      19516 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  100 :nEle      : nEle/I                                                 *
*Entries :    12429 : Total  Size=      50338 bytes  File Size  =       7869 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   6.34     *
*............................................................................*
*Br  101 :eleCharge : vector<int>                                            *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =      51080 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.50     *
*............................................................................*
*Br  102 :eleChargeConsistent : vector<int>                                  *
*Entries :    12429 : Total  Size=     230647 bytes  File Size  =      50520 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.55     *
*............................................................................*
*Br  103 :eleEn     : vector<float>                                          *
*Entries :    12429 : Total  Size=     230465 bytes  File Size  =     104143 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.21     *
*............................................................................*
*Br  104 :eleSCEn   : vector<float>                                          *
*Entries :    12429 : Total  Size=     230491 bytes  File Size  =     104453 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.20     *
*............................................................................*
*Br  105 :eleEcalEn : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =     102576 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.24     *
*............................................................................*
*Br  106 :eleESEnP1 : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =      74823 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.07     *
*............................................................................*
*Br  107 :eleESEnP2 : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =      75027 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.07     *
*............................................................................*
*Br  108 :eleD0     : vector<float>                                          *
*Entries :    12429 : Total  Size=     230465 bytes  File Size  =     107491 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.14     *
*............................................................................*
*Br  109 :eleDz     : vector<float>                                          *
*Entries :    12429 : Total  Size=     230465 bytes  File Size  =     108019 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.13     *
*............................................................................*
*Br  110 :eleSIP    : vector<float>                                          *
*Entries :    12429 : Total  Size=     230478 bytes  File Size  =     104921 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br  111 :elePt     : vector<float>                                          *
*Entries :    12429 : Total  Size=     230465 bytes  File Size  =     103436 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.22     *
*............................................................................*
*Br  112 :eleEta    : vector<float>                                          *
*Entries :    12429 : Total  Size=     230478 bytes  File Size  =     105344 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.18     *
*............................................................................*
*Br  113 :elePhi    : vector<float>                                          *
*Entries :    12429 : Total  Size=     230478 bytes  File Size  =     105622 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.18     *
*............................................................................*
*Br  114 :eleR9     : vector<float>                                          *
*Entries :    12429 : Total  Size=     230465 bytes  File Size  =      91700 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.51     *
*............................................................................*
*Br  115 :eleCalibPt : vector<float>                                         *
*Entries :    12429 : Total  Size=     230530 bytes  File Size  =     103316 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.23     *
*............................................................................*
*Br  116 :eleCalibEn : vector<float>                                         *
*Entries :    12429 : Total  Size=     230530 bytes  File Size  =     104097 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.21     *
*............................................................................*
*Br  117 :eleSCEta  : vector<float>                                          *
*Entries :    12429 : Total  Size=     230504 bytes  File Size  =     105475 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.18     *
*............................................................................*
*Br  118 :eleSCPhi  : vector<float>                                          *
*Entries :    12429 : Total  Size=     230504 bytes  File Size  =     105543 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.18     *
*............................................................................*
*Br  119 :eleSCRawEn : vector<float>                                         *
*Entries :    12429 : Total  Size=     230530 bytes  File Size  =     104456 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.20     *
*............................................................................*
*Br  120 :eleSCEtaWidth : vector<float>                                      *
*Entries :    12429 : Total  Size=     230569 bytes  File Size  =     102497 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.24     *
*............................................................................*
*Br  121 :eleSCPhiWidth : vector<float>                                      *
*Entries :    12429 : Total  Size=     230569 bytes  File Size  =     103166 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.23     *
*............................................................................*
*Br  122 :eleHoverE : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =      81868 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.81     *
*............................................................................*
*Br  123 :eleEoverP : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =      96964 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.37     *
*............................................................................*
*Br  124 :eleEoverPout : vector<float>                                       *
*Entries :    12429 : Total  Size=     230556 bytes  File Size  =     103921 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.21     *
*............................................................................*
*Br  125 :eleEoverPInv : vector<float>                                       *
*Entries :    12429 : Total  Size=     230556 bytes  File Size  =      96082 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.39     *
*............................................................................*
*Br  126 :eleBrem   : vector<float>                                          *
*Entries :    12429 : Total  Size=     230491 bytes  File Size  =      95336 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.41     *
*............................................................................*
*Br  127 :eledEtaAtVtx : vector<float>                                       *
*Entries :    12429 : Total  Size=     230556 bytes  File Size  =      98781 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.33     *
*............................................................................*
*Br  128 :eledPhiAtVtx : vector<float>                                       *
*Entries :    12429 : Total  Size=     230556 bytes  File Size  =     101938 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.26     *
*............................................................................*
*Br  129 :eledEtaAtCalo : vector<float>                                      *
*Entries :    12429 : Total  Size=     230569 bytes  File Size  =      99695 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.31     *
*............................................................................*
*Br  130 :eleSigmaIEtaIEtaFull5x5 : vector<float>                            *
*Entries :    12429 : Total  Size=     230699 bytes  File Size  =     100919 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.28     *
*............................................................................*
*Br  131 :eleSigmaIPhiIPhiFull5x5 : vector<float>                            *
*Entries :    12429 : Total  Size=     230699 bytes  File Size  =     101739 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.26     *
*............................................................................*
*Br  132 :eleConvVeto : vector<int>                                          *
*Entries :    12429 : Total  Size=     230543 bytes  File Size  =      46513 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.95     *
*............................................................................*
*Br  133 :eleMissHits : vector<int>                                          *
*Entries :    12429 : Total  Size=     230543 bytes  File Size  =      50594 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.55     *
*............................................................................*
*Br  134 :eleESEffSigmaRR : vector<float>                                    *
*Entries :    12429 : Total  Size=     230595 bytes  File Size  =      69766 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.30     *
*............................................................................*
*Br  135 :elePFChIso : vector<float>                                         *
*Entries :    12429 : Total  Size=     230530 bytes  File Size  =      90987 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.53     *
*............................................................................*
*Br  136 :elePFPhoIso : vector<float>                                        *
*Entries :    12429 : Total  Size=     230543 bytes  File Size  =      93176 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.47     *
*............................................................................*
*Br  137 :elePFNeuIso : vector<float>                                        *
*Entries :    12429 : Total  Size=     230543 bytes  File Size  =      72321 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.18     *
*............................................................................*
*Br  138 :elePFPUIso : vector<float>                                         *
*Entries :    12429 : Total  Size=     230530 bytes  File Size  =      91938 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.50     *
*............................................................................*
*Br  139 :elePFClusEcalIso : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      93915 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.45     *
*............................................................................*
*Br  140 :elePFClusHcalIso : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      88527 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.60     *
*............................................................................*
*Br  141 :elePFMiniIso : vector<float>                                       *
*Entries :    12429 : Total  Size=     230556 bytes  File Size  =      91927 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.50     *
*............................................................................*
*Br  142 :eleIDMVA  : vector<float>                                          *
*Entries :    12429 : Total  Size=     230504 bytes  File Size  =     101500 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.27     *
*............................................................................*
*Br  143 :eleIDMVAHZZ : vector<float>                                        *
*Entries :    12429 : Total  Size=     230543 bytes  File Size  =     102024 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.25     *
*............................................................................*
*Br  144 :eledEtaseedAtVtx : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =     105084 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.19     *
*............................................................................*
*Br  145 :eleE1x5   : vector<float>                                          *
*Entries :    12429 : Total  Size=     230491 bytes  File Size  =      96124 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.39     *
*............................................................................*
*Br  146 :eleE2x5   : vector<float>                                          *
*Entries :    12429 : Total  Size=     230491 bytes  File Size  =      96103 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.39     *
*............................................................................*
*Br  147 :eleE5x5   : vector<float>                                          *
*Entries :    12429 : Total  Size=     230491 bytes  File Size  =      96209 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.39     *
*............................................................................*
*Br  148 :eleE1x5Full5x5 : vector<float>                                     *
*Entries :    12429 : Total  Size=     230582 bytes  File Size  =     104373 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.20     *
*............................................................................*
*Br  149 :eleE2x5Full5x5 : vector<float>                                     *
*Entries :    12429 : Total  Size=     230582 bytes  File Size  =     104291 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.21     *
*............................................................................*
*Br  150 :eleE5x5Full5x5 : vector<float>                                     *
*Entries :    12429 : Total  Size=     230582 bytes  File Size  =     104549 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.20     *
*............................................................................*
*Br  151 :eleR9Full5x5 : vector<float>                                       *
*Entries :    12429 : Total  Size=     230556 bytes  File Size  =      98751 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.33     *
*............................................................................*
*Br  152 :eleEcalDrivenSeed : vector<int>                                    *
*Entries :    12429 : Total  Size=     230621 bytes  File Size  =      49552 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.64     *
*............................................................................*
*Br  153 :eleDr03EcalRecHitSumEt : vector<float>                             *
*Entries :    12429 : Total  Size=     230686 bytes  File Size  =      96095 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.39     *
*............................................................................*
*Br  154 :eleDr03HcalDepth1TowerSumEt : vector<float>                        *
*Entries :    12429 : Total  Size=     230751 bytes  File Size  =      89881 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.56     *
*............................................................................*
*Br  155 :eleDr03HcalDepth2TowerSumEt : vector<float>                        *
*Entries :    12429 : Total  Size=     230751 bytes  File Size  =      61232 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.76     *
*............................................................................*
*Br  156 :eleDr03HcalTowerSumEt : vector<float>                              *
*Entries :    12429 : Total  Size=     230673 bytes  File Size  =      90300 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.55     *
*............................................................................*
*Br  157 :eleDr03TkSumPt : vector<float>                                     *
*Entries :    12429 : Total  Size=     230582 bytes  File Size  =      87249 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.64     *
*............................................................................*
*Br  158 :elecaloEnergy : vector<float>                                      *
*Entries :    12429 : Total  Size=     230569 bytes  File Size  =     102612 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.24     *
*............................................................................*
*Br  159 :eleTrkdxy : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =     107505 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.14     *
*............................................................................*
*Br  160 :eleKFHits : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =      66891 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.44     *
*............................................................................*
*Br  161 :eleKFChi2 : vector<float>                                          *
*Entries :    12429 : Total  Size=     230517 bytes  File Size  =      98705 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.33     *
*............................................................................*
*Br  162 :eleGSFChi2 : vector<float>                                         *
*Entries :    12429 : Total  Size=     230530 bytes  File Size  =     103131 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.23     *
*............................................................................*
*Br  163 :eleGSFPt  : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175200 bytes  File Size  =      19516 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  164 :eleGSFEta : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175212 bytes  File Size  =      19573 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.93     *
*............................................................................*
*Br  165 :eleGSFPhi : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175212 bytes  File Size  =      19573 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.93     *
*............................................................................*
*Br  166 :eleGSFCharge : vector<vector<float> >                              *
*Entries :    12429 : Total  Size=     175248 bytes  File Size  =      19667 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.89     *
*............................................................................*
*Br  167 :eleGSFHits : vector<vector<int> >                                  *
*Entries :    12429 : Total  Size=     175224 bytes  File Size  =      19539 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  168 :eleGSFMissHits : vector<vector<int> >                              *
*Entries :    12429 : Total  Size=     175272 bytes  File Size  =      19605 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.91     *
*............................................................................*
*Br  169 :eleGSFNHitsMax : vector<vector<int> >                              *
*Entries :    12429 : Total  Size=     175272 bytes  File Size  =      19605 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.91     *
*............................................................................*
*Br  170 :eleGSFVtxProb : vector<vector<float> >                             *
*Entries :    12429 : Total  Size=     175260 bytes  File Size  =      19619 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.91     *
*............................................................................*
*Br  171 :eleGSFlxyPV : vector<vector<float> >                               *
*Entries :    12429 : Total  Size=     175236 bytes  File Size  =      19575 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.93     *
*............................................................................*
*Br  172 :eleGSFlxyBS : vector<vector<float> >                               *
*Entries :    12429 : Total  Size=     175236 bytes  File Size  =      19575 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.93     *
*............................................................................*
*Br  173 :eleBCEn   : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175188 bytes  File Size  =      19543 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  174 :eleBCEta  : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175200 bytes  File Size  =      19516 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  175 :eleBCPhi  : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175200 bytes  File Size  =      19516 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  176 :eleBCS25  : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175200 bytes  File Size  =      19516 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  177 :eleBCS15  : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     175200 bytes  File Size  =      19516 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  178 :eleBCSieie : vector<vector<float> >                                *
*Entries :    12429 : Total  Size=     175224 bytes  File Size  =      19539 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  179 :eleBCSieip : vector<vector<float> >                                *
*Entries :    12429 : Total  Size=     175224 bytes  File Size  =      19539 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  180 :eleBCSipip : vector<vector<float> >                                *
*Entries :    12429 : Total  Size=     175224 bytes  File Size  =      19539 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  181 :eleFiredSingleTrgs : vector<unsigned int>                          *
*Entries :    12429 : Total  Size=     230634 bytes  File Size  =      60184 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.82     *
*............................................................................*
*Br  182 :eleFiredDoubleTrgs : vector<unsigned int>                          *
*Entries :    12429 : Total  Size=     230634 bytes  File Size  =      43028 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   5.35     *
*............................................................................*
*Br  183 :eleFiredL1Trgs : vector<unsigned int>                              *
*Entries :    12429 : Total  Size=     230582 bytes  File Size  =      72642 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.17     *
*............................................................................*
*Br  184 :eleIDbit  : vector<unsigned short>                                 *
*Entries :    12429 : Total  Size=     202808 bytes  File Size  =      46593 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.34     *
*............................................................................*
*Br  185 :eleScale_stat_up : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      51991 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.43     *
*............................................................................*
*Br  186 :eleScale_stat_dn : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      53597 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.29     *
*............................................................................*
*Br  187 :eleScale_syst_up : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      59465 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.87     *
*............................................................................*
*Br  188 :eleScale_syst_dn : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      59802 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.85     *
*............................................................................*
*Br  189 :eleScale_gain_up : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      42379 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   5.43     *
*............................................................................*
*Br  190 :eleScale_gain_dn : vector<float>                                   *
*Entries :    12429 : Total  Size=     230608 bytes  File Size  =      42386 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   5.43     *
*............................................................................*
*Br  191 :eleResol_rho_up : vector<float>                                    *
*Entries :    12429 : Total  Size=     230595 bytes  File Size  =      61186 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.76     *
*............................................................................*
*Br  192 :eleResol_rho_dn : vector<float>                                    *
*Entries :    12429 : Total  Size=     230595 bytes  File Size  =      61196 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.76     *
*............................................................................*
*Br  193 :eleResol_phi_up : vector<float>                                    *
*Entries :    12429 : Total  Size=     230595 bytes  File Size  =     100225 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.30     *
*............................................................................*
*Br  194 :eleResol_phi_dn : vector<float>                                    *
*Entries :    12429 : Total  Size=     230595 bytes  File Size  =     100225 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.30     *
*............................................................................*
*Br  195 :npfHF     : npfHF/I                                                *
*Entries :    12429 : Total  Size=      50344 bytes  File Size  =       2104 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=  23.70     *
*............................................................................*
*Br  196 :pfHFEn    : vector<float>                                          *
*Entries :    12429 : Total  Size=     177592 bytes  File Size  =      29234 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.06     *
*............................................................................*
*Br  197 :pfHFECALEn : vector<float>                                         *
*Entries :    12429 : Total  Size=     177640 bytes  File Size  =      29252 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.06     *
*............................................................................*
*Br  198 :pfHFHCALEn : vector<float>                                         *
*Entries :    12429 : Total  Size=     177640 bytes  File Size  =      29015 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.11     *
*............................................................................*
*Br  199 :pfHFPt    : vector<float>                                          *
*Entries :    12429 : Total  Size=     177592 bytes  File Size  =      28541 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.21     *
*............................................................................*
*Br  200 :pfHFEta   : vector<float>                                          *
*Entries :    12429 : Total  Size=     177604 bytes  File Size  =      29239 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.06     *
*............................................................................*
*Br  201 :pfHFPhi   : vector<float>                                          *
*Entries :    12429 : Total  Size=     177604 bytes  File Size  =      29604 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.98     *
*............................................................................*
*Br  202 :pfHFIso   : vector<float>                                          *
*Entries :    12429 : Total  Size=     177604 bytes  File Size  =      29522 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.00     *
*............................................................................*
*Br  203 :nMu       : nMu/I                                                  *
*Entries :    12429 : Total  Size=      50332 bytes  File Size  =       5395 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   9.24     *
*............................................................................*
*Br  204 :muPt      : vector<float>                                          *
*Entries :    12429 : Total  Size=     238072 bytes  File Size  =     103840 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.29     *
*............................................................................*
*Br  205 :muEn      : vector<float>                                          *
*Entries :    12429 : Total  Size=     238072 bytes  File Size  =     105893 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.24     *
*............................................................................*
*Br  206 :muEta     : vector<float>                                          *
*Entries :    12429 : Total  Size=     238085 bytes  File Size  =     109331 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.17     *
*............................................................................*
*Br  207 :muPhi     : vector<float>                                          *
*Entries :    12429 : Total  Size=     238085 bytes  File Size  =     109288 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.17     *
*............................................................................*
*Br  208 :muCharge  : vector<int>                                            *
*Entries :    12429 : Total  Size=     238124 bytes  File Size  =      46364 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   5.13     *
*............................................................................*
*Br  209 :muType    : vector<int>                                            *
*Entries :    12429 : Total  Size=     238098 bytes  File Size  =      50262 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.73     *
*............................................................................*
*Br  210 :muIDbit   : vector<unsigned short>                                 *
*Entries :    12429 : Total  Size=     206606 bytes  File Size  =      46979 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.39     *
*............................................................................*
*Br  211 :muD0      : vector<float>                                          *
*Entries :    12429 : Total  Size=     238072 bytes  File Size  =     111360 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.13     *
*............................................................................*
*Br  212 :muDz      : vector<float>                                          *
*Entries :    12429 : Total  Size=     238072 bytes  File Size  =     111611 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.13     *
*............................................................................*
*Br  213 :muSIP     : vector<float>                                          *
*Entries :    12429 : Total  Size=     238085 bytes  File Size  =     109109 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.18     *
*............................................................................*
*Br  214 :muChi2NDF : vector<float>                                          *
*Entries :    12429 : Total  Size=     238137 bytes  File Size  =     100580 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.36     *
*............................................................................*
*Br  215 :muInnerD0 : vector<float>                                          *
*Entries :    12429 : Total  Size=     238137 bytes  File Size  =     111483 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.13     *
*............................................................................*
*Br  216 :muInnerDz : vector<float>                                          *
*Entries :    12429 : Total  Size=     238137 bytes  File Size  =     111771 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.13     *
*............................................................................*
*Br  217 :muTrkLayers : vector<int>                                          *
*Entries :    12429 : Total  Size=     238163 bytes  File Size  =      55481 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.28     *
*............................................................................*
*Br  218 :muPixelLayers : vector<int>                                        *
*Entries :    12429 : Total  Size=     238189 bytes  File Size  =      47105 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   5.05     *
*............................................................................*
*Br  219 :muPixelHits : vector<int>                                          *
*Entries :    12429 : Total  Size=     238163 bytes  File Size  =      50336 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.72     *
*............................................................................*
*Br  220 :muMuonHits : vector<int>                                           *
*Entries :    12429 : Total  Size=     238150 bytes  File Size  =      61806 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.85     *
*............................................................................*
*Br  221 :muStations : vector<int>                                           *
*Entries :    12429 : Total  Size=     238150 bytes  File Size  =      52143 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.56     *
*............................................................................*
*Br  222 :muMatches : vector<int>                                            *
*Entries :    12429 : Total  Size=     238137 bytes  File Size  =      53347 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.45     *
*............................................................................*
*Br  223 :muTrkQuality : vector<int>                                         *
*Entries :    12429 : Total  Size=     238176 bytes  File Size  =      38768 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   6.13     *
*............................................................................*
*Br  224 :muIsoTrk  : vector<float>                                          *
*Entries :    12429 : Total  Size=     238124 bytes  File Size  =     102279 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.32     *
*............................................................................*
*Br  225 :muPFChIso : vector<float>                                          *
*Entries :    12429 : Total  Size=     238137 bytes  File Size  =     104804 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.27     *
*............................................................................*
*Br  226 :muPFPhoIso : vector<float>                                         *
*Entries :    12429 : Total  Size=     238150 bytes  File Size  =     102984 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.31     *
*............................................................................*
*Br  227 :muPFNeuIso : vector<float>                                         *
*Entries :    12429 : Total  Size=     238150 bytes  File Size  =      84068 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.83     *
*............................................................................*
*Br  228 :muPFPUIso : vector<float>                                          *
*Entries :    12429 : Total  Size=     238137 bytes  File Size  =     103169 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.30     *
*............................................................................*
*Br  229 :muPFChIso03 : vector<float>                                        *
*Entries :    12429 : Total  Size=     238163 bytes  File Size  =     102655 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.32     *
*............................................................................*
*Br  230 :muPFPhoIso03 : vector<float>                                       *
*Entries :    12429 : Total  Size=     238176 bytes  File Size  =     100832 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.36     *
*............................................................................*
*Br  231 :muPFNeuIso03 : vector<float>                                       *
*Entries :    12429 : Total  Size=     238176 bytes  File Size  =      74531 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.19     *
*............................................................................*
*Br  232 :muPFPUIso03 : vector<float>                                        *
*Entries :    12429 : Total  Size=     238163 bytes  File Size  =      98664 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.41     *
*............................................................................*
*Br  233 :muPFMiniIso : vector<float>                                        *
*Entries :    12429 : Total  Size=     238163 bytes  File Size  =     101148 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.35     *
*............................................................................*
*Br  234 :muFiredTrgs : vector<unsigned int>                                 *
*Entries :    12429 : Total  Size=     238163 bytes  File Size  =      61623 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.86     *
*............................................................................*
*Br  235 :muFiredL1Trgs : vector<unsigned int>                               *
*Entries :    12429 : Total  Size=     238189 bytes  File Size  =      70789 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.36     *
*............................................................................*
*Br  236 :muInnervalidFraction : vector<float>                               *
*Entries :    12429 : Total  Size=     238280 bytes  File Size  =      54371 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.37     *
*............................................................................*
*Br  237 :musegmentCompatibility : vector<float>                             *
*Entries :    12429 : Total  Size=     238306 bytes  File Size  =      96534 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br  238 :muchi2LocalPosition : vector<float>                                *
*Entries :    12429 : Total  Size=     238267 bytes  File Size  =     105897 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.24     *
*............................................................................*
*Br  239 :mutrkKink : vector<float>                                          *
*Entries :    12429 : Total  Size=     238137 bytes  File Size  =     105879 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.24     *
*............................................................................*
*Br  240 :muBestTrkPtError : vector<float>                                   *
*Entries :    12429 : Total  Size=     238228 bytes  File Size  =     106020 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.24     *
*............................................................................*
*Br  241 :muBestTrkPt : vector<float>                                        *
*Entries :    12429 : Total  Size=     238163 bytes  File Size  =     104031 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   2.28     *
*............................................................................*
*Br  242 :muBestTrkType : vector<int>                                        *
*Entries :    12429 : Total  Size=     238189 bytes  File Size  =      37987 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   6.26     *
*............................................................................*
*Br  243 :nJet      : nJet/I                                                 *
*Entries :    12429 : Total  Size=      50338 bytes  File Size  =       9036 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   5.52     *
*............................................................................*
*Br  244 :jetPt     : vector<float>                                          *
*Entries :    12429 : Total  Size=     347281 bytes  File Size  =     217524 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.59     *
*............................................................................*
*Br  245 :jetEn     : vector<float>                                          *
*Entries :    12429 : Total  Size=     347281 bytes  File Size  =     220138 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.58     *
*............................................................................*
*Br  246 :jetEta    : vector<float>                                          *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     225616 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.54     *
*............................................................................*
*Br  247 :jetPhi    : vector<float>                                          *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     225145 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.54     *
*............................................................................*
*Br  248 :jetRawPt  : vector<float>                                          *
*Entries :    12429 : Total  Size=     347332 bytes  File Size  =     217924 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.59     *
*............................................................................*
*Br  249 :jetRawEn  : vector<float>                                          *
*Entries :    12429 : Total  Size=     347332 bytes  File Size  =     220049 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.58     *
*............................................................................*
*Br  250 :jetMt     : vector<float>                                          *
*Entries :    12429 : Total  Size=     347281 bytes  File Size  =     217635 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.59     *
*............................................................................*
*Br  251 :jetArea   : vector<float>                                          *
*Entries :    12429 : Total  Size=     347315 bytes  File Size  =     112902 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.07     *
*............................................................................*
*Br  252 :jetLeadTrackPt : vector<float>                                     *
*Entries :    12429 : Total  Size=     347434 bytes  File Size  =     163719 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   2.12     *
*............................................................................*
*Br  253 :jetLeadTrackEta : vector<float>                                    *
*Entries :    12429 : Total  Size=     347451 bytes  File Size  =     208936 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.66     *
*............................................................................*
*Br  254 :jetLeadTrackPhi : vector<float>                                    *
*Entries :    12429 : Total  Size=     347451 bytes  File Size  =     215316 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.61     *
*............................................................................*
*Br  255 :jetLepTrackPID : vector<int>                                       *
*Entries :    12429 : Total  Size=     347434 bytes  File Size  =      89379 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.88     *
*............................................................................*
*Br  256 :jetLepTrackPt : vector<float>                                      *
*Entries :    12429 : Total  Size=     347417 bytes  File Size  =     118721 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   2.92     *
*............................................................................*
*Br  257 :jetLepTrackEta : vector<float>                                     *
*Entries :    12429 : Total  Size=     347434 bytes  File Size  =     141716 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   2.45     *
*............................................................................*
*Br  258 :jetLepTrackPhi : vector<float>                                     *
*Entries :    12429 : Total  Size=     347434 bytes  File Size  =     143136 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   2.42     *
*............................................................................*
*Br  259 :jetCSV2BJetTags : vector<float>                                    *
*Entries :    12429 : Total  Size=     347451 bytes  File Size  =     205880 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.68     *
*............................................................................*
*Br  260 :jetJetProbabilityBJetTags : vector<float>                          *
*Entries :    12429 : Total  Size=     347621 bytes  File Size  =     195250 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.78     *
*............................................................................*
*Br  261 :jetpfCombinedMVAV2BJetTags : vector<float>                         *
*Entries :    12429 : Total  Size=     347638 bytes  File Size  =     190756 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.82     *
*............................................................................*
*Br  262 :jetPFLooseId : vector<bool>                                        *
*Entries :    12429 : Total  Size=     218243 bytes  File Size  =      48813 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   4.46     *
*............................................................................*
*Br  263 :jetID     : vector<int>                                            *
*Entries :    12429 : Total  Size=     347281 bytes  File Size  =      56944 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   6.09     *
*............................................................................*
*Br  264 :jetPUID   : vector<float>                                          *
*Entries :    12429 : Total  Size=     347315 bytes  File Size  =     214418 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.62     *
*............................................................................*
*Br  265 :jetPUFullID : vector<int>                                          *
*Entries :    12429 : Total  Size=     347383 bytes  File Size  =      77828 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   4.46     *
*............................................................................*
*Br  266 :jetJECUnc : vector<float>                                          *
*Entries :    12429 : Total  Size=     347349 bytes  File Size  =     216125 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.60     *
*............................................................................*
*Br  267 :jetFiredTrgs : vector<unsigned int>                                *
*Entries :    12429 : Total  Size=     347400 bytes  File Size  =      50106 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   6.92     *
*............................................................................*
*Br  268 :jetCHF    : vector<float>                                          *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     212281 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.63     *
*............................................................................*
*Br  269 :jetNHF    : vector<float>                                          *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     190650 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.82     *
*............................................................................*
*Br  270 :jetCEF    : vector<float>                                          *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =      96377 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.60     *
*............................................................................*
*Br  271 :jetNEF    : vector<float>                                          *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     220552 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   1.57     *
*............................................................................*
*Br  272 :jetNCH    : vector<int>                                            *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     107926 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.21     *
*............................................................................*
*Br  273 :jetNNP    : vector<int>                                            *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     108749 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.19     *
*............................................................................*
*Br  274 :jetMUF    : vector<float>                                          *
*Entries :    12429 : Total  Size=     347298 bytes  File Size  =     128727 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   2.69     *
*............................................................................*
*Br  275 :jetVtxPt  : vector<float>                                          *
*Entries :    12429 : Total  Size=     347332 bytes  File Size  =      97391 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.56     *
*............................................................................*
*Br  276 :jetVtxMass : vector<float>                                         *
*Entries :    12429 : Total  Size=     347366 bytes  File Size  =      96957 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.58     *
*............................................................................*
*Br  277 :jetVtxNtrks : vector<float>                                        *
*Entries :    12429 : Total  Size=     347383 bytes  File Size  =      76849 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   4.51     *
*............................................................................*
*Br  278 :jetVtx3DVal : vector<float>                                        *
*Entries :    12429 : Total  Size=     347383 bytes  File Size  =      97758 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.55     *
*............................................................................*
*Br  279 :jetVtx3DSig : vector<float>                                        *
*Entries :    12429 : Total  Size=     347383 bytes  File Size  =      97291 *
*Baskets :       13 : Basket Size=      32000 bytes  Compression=   3.56     *
*............................................................................*
*Br  280 :nAK8Jet   : nAK8Jet/I                                              *
*Entries :    12429 : Total  Size=      50356 bytes  File Size  =       3499 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=  14.25     *
*............................................................................*
*Br  281 :AK8JetPt  : vector<float>                                          *
*Entries :    12429 : Total  Size=     182472 bytes  File Size  =      38418 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.74     *
*............................................................................*
*Br  282 :AK8JetEn  : vector<float>                                          *
*Entries :    12429 : Total  Size=     182472 bytes  File Size  =      38793 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.69     *
*............................................................................*
*Br  283 :AK8JetRawPt : vector<float>                                        *
*Entries :    12429 : Total  Size=     182508 bytes  File Size  =      38558 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.72     *
*............................................................................*
*Br  284 :AK8JetRawEn : vector<float>                                        *
*Entries :    12429 : Total  Size=     182508 bytes  File Size  =      38972 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.67     *
*............................................................................*
*Br  285 :AK8JetEta : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      39482 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.61     *
*............................................................................*
*Br  286 :AK8JetPhi : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      39487 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.61     *
*............................................................................*
*Br  287 :AK8JetMass : vector<float>                                         *
*Entries :    12429 : Total  Size=     182496 bytes  File Size  =      38546 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.72     *
*............................................................................*
*Br  288 :AK8Jet_tau1 : vector<float>                                        *
*Entries :    12429 : Total  Size=     182508 bytes  File Size  =      38938 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.67     *
*............................................................................*
*Br  289 :AK8Jet_tau2 : vector<float>                                        *
*Entries :    12429 : Total  Size=     182508 bytes  File Size  =      39008 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.67     *
*............................................................................*
*Br  290 :AK8Jet_tau3 : vector<float>                                        *
*Entries :    12429 : Total  Size=     182508 bytes  File Size  =      38899 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.68     *
*............................................................................*
*Br  291 :AK8JetCHF : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      39033 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.66     *
*............................................................................*
*Br  292 :AK8JetNHF : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      39172 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.65     *
*............................................................................*
*Br  293 :AK8JetCEF : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      31985 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.69     *
*............................................................................*
*Br  294 :AK8JetNEF : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      39027 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.66     *
*............................................................................*
*Br  295 :AK8JetNCH : vector<int>                                            *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      34125 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.33     *
*............................................................................*
*Br  296 :AK8JetNNP : vector<int>                                            *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      34080 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.34     *
*............................................................................*
*Br  297 :AK8JetMUF : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      34438 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.28     *
*............................................................................*
*Br  298 :AK8Jetnconstituents : vector<int>                                  *
*Entries :    12429 : Total  Size=     182604 bytes  File Size  =      34384 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.30     *
*............................................................................*
*Br  299 :AK8JetPFLooseId : vector<bool>                                     *
*Entries :    12429 : Total  Size=     177102 bytes  File Size  =      30287 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.83     *
*............................................................................*
*Br  300 :AK8JetPFTightLepVetoId : vector<bool>                              *
*Entries :    12429 : Total  Size=     177186 bytes  File Size  =      30441 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.80     *
*............................................................................*
*Br  301 :AK8JetSoftDropMass : vector<float>                                 *
*Entries :    12429 : Total  Size=     182592 bytes  File Size  =      39211 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.64     *
*............................................................................*
*Br  302 :AK8JetSoftDropMassCorr : vector<float>                             *
*Entries :    12429 : Total  Size=     182640 bytes  File Size  =      39302 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.63     *
*............................................................................*
*Br  303 :AK8JetPrunedMass : vector<float>                                   *
*Entries :    12429 : Total  Size=     182568 bytes  File Size  =      39026 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.66     *
*............................................................................*
*Br  304 :AK8JetPrunedMassCorr : vector<float>                               *
*Entries :    12429 : Total  Size=     182616 bytes  File Size  =      39140 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.65     *
*............................................................................*
*Br  305 :AK8JetpfBoostedDSVBTag : vector<float>                             *
*Entries :    12429 : Total  Size=     182640 bytes  File Size  =      39479 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.61     *
*............................................................................*
*Br  306 :AK8JetDSVnewV4 : vector<float>                                     *
*Entries :    12429 : Total  Size=     182544 bytes  File Size  =      38962 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.67     *
*............................................................................*
*Br  307 :AK8JetCSV : vector<float>                                          *
*Entries :    12429 : Total  Size=     182484 bytes  File Size  =      39088 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.66     *
*............................................................................*
*Br  308 :AK8JetJECUnc : vector<float>                                       *
*Entries :    12429 : Total  Size=     182520 bytes  File Size  =      38815 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.69     *
*............................................................................*
*Br  309 :AK8JetL2L3corr : vector<float>                                     *
*Entries :    12429 : Total  Size=     182544 bytes  File Size  =      38192 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.77     *
*............................................................................*
*Br  310 :AK8puppiPt : vector<float>                                         *
*Entries :    12429 : Total  Size=     182496 bytes  File Size  =      38474 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.73     *
*............................................................................*
*Br  311 :AK8puppiMass : vector<float>                                       *
*Entries :    12429 : Total  Size=     182520 bytes  File Size  =      38836 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.69     *
*............................................................................*
*Br  312 :AK8puppiEta : vector<float>                                        *
*Entries :    12429 : Total  Size=     182508 bytes  File Size  =      39464 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.61     *
*............................................................................*
*Br  313 :AK8puppiPhi : vector<float>                                        *
*Entries :    12429 : Total  Size=     182508 bytes  File Size  =      39501 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.61     *
*............................................................................*
*Br  314 :AK8puppiTau1 : vector<float>                                       *
*Entries :    12429 : Total  Size=     182520 bytes  File Size  =      38964 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.67     *
*............................................................................*
*Br  315 :AK8puppiTau2 : vector<float>                                       *
*Entries :    12429 : Total  Size=     182520 bytes  File Size  =      38922 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.68     *
*............................................................................*
*Br  316 :AK8puppiTau3 : vector<float>                                       *
*Entries :    12429 : Total  Size=     182520 bytes  File Size  =      38920 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.68     *
*............................................................................*
*Br  317 :AK8puppiSDL2L3corr : vector<float>                                 *
*Entries :    12429 : Total  Size=     182592 bytes  File Size  =      36822 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.94     *
*............................................................................*
*Br  318 :AK8puppiSDMass : vector<float>                                     *
*Entries :    12429 : Total  Size=     182544 bytes  File Size  =      37400 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.87     *
*............................................................................*
*Br  319 :AK8puppiSDMassL2L3Corr : vector<float>                             *
*Entries :    12429 : Total  Size=     182640 bytes  File Size  =      37406 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.87     *
*............................................................................*
*Br  320 :nAK8SDSJ  : vector<int>                                            *
*Entries :    12429 : Total  Size=     182472 bytes  File Size  =      30899 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.89     *
*............................................................................*
*Br  321 :AK8SDSJPt : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     193240 bytes  File Size  =      44287 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.35     *
*............................................................................*
*Br  322 :AK8SDSJEta : vector<vector<float> >                                *
*Entries :    12429 : Total  Size=     193252 bytes  File Size  =      44812 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.30     *
*............................................................................*
*Br  323 :AK8SDSJPhi : vector<vector<float> >                                *
*Entries :    12429 : Total  Size=     193252 bytes  File Size  =      44732 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.31     *
*............................................................................*
*Br  324 :AK8SDSJMass : vector<vector<float> >                               *
*Entries :    12429 : Total  Size=     193264 bytes  File Size  =      44703 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.31     *
*............................................................................*
*Br  325 :AK8SDSJE  : vector<vector<float> >                                 *
*Entries :    12429 : Total  Size=     193228 bytes  File Size  =      44184 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.36     *
*............................................................................*
*Br  326 :AK8SDSJCharge : vector<vector<int> >                               *
*Entries :    12429 : Total  Size=     193288 bytes  File Size  =      37356 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.16     *
*............................................................................*
*Br  327 :AK8SDSJFlavour : vector<vector<int> >                              *
*Entries :    12429 : Total  Size=     193300 bytes  File Size  =      30371 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.35     *
*............................................................................*
*Br  328 :AK8SDSJCSV : vector<vector<float> >                                *
*Entries :    12429 : Total  Size=     193252 bytes  File Size  =      43603 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.42     *
*............................................................................*
*Br  329 :nAK8puppiSDSJ : vector<int>                                        *
*Entries :    12429 : Total  Size=     182532 bytes  File Size  =      31211 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.83     *
*............................................................................*
*Br  330 :AK8puppiSDSJPt : vector<vector<float> >                            *
*Entries :    12429 : Total  Size=     192980 bytes  File Size  =      43735 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.40     *
*............................................................................*
*Br  331 :AK8puppiSDSJEta : vector<vector<float> >                           *
*Entries :    12429 : Total  Size=     192992 bytes  File Size  =      44624 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.31     *
*............................................................................*
*Br  332 :AK8puppiSDSJPhi : vector<vector<float> >                           *
*Entries :    12429 : Total  Size=     192992 bytes  File Size  =      44595 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.32     *
*............................................................................*
*Br  333 :AK8puppiSDSJMass : vector<vector<float> >                          *
*Entries :    12429 : Total  Size=     193004 bytes  File Size  =      44099 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.36     *
*............................................................................*
*Br  334 :AK8puppiSDSJE : vector<vector<float> >                             *
*Entries :    12429 : Total  Size=     192968 bytes  File Size  =      44018 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.37     *
*............................................................................*
*Br  335 :AK8puppiSDSJCharge : vector<vector<int> >                          *
*Entries :    12429 : Total  Size=     193028 bytes  File Size  =      36700 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.25     *
*............................................................................*
*Br  336 :AK8puppiSDSJFlavour : vector<vector<int> >                         *
*Entries :    12429 : Total  Size=     193040 bytes  File Size  =      30547 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.30     *
*............................................................................*
*Br  337 :AK8puppiSDSJCSV : vector<vector<float> >                           *
*Entries :    12429 : Total  Size=     192992 bytes  File Size  =      43502 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.42     *
*............................................................................*
*Br  338 :ndiMu     : ndiMu/I                                                *
*Entries :    12429 : Total  Size=      50344 bytes  File Size  =       5227 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   9.54     *
*............................................................................*
*Br  339 :diMuIndex1 : vector<int>                                           *
*Entries :    12429 : Total  Size=     191728 bytes  File Size  =      34597 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.53     *
*............................................................................*
*Br  340 :diMuIndex2 : vector<int>                                           *
*Entries :    12429 : Total  Size=     191728 bytes  File Size  =      34643 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.52     *
*............................................................................*
*Br  341 :diMuVtxIsValid : vector<int>                                       *
*Entries :    12429 : Total  Size=     191776 bytes  File Size  =      34395 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.56     *
*............................................................................*
*Br  342 :diMuVtx   : vector<float>                                          *
*Entries :    12429 : Total  Size=     191692 bytes  File Size  =      54557 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.50     *
*............................................................................*
*Br  343 :diMuVty   : vector<float>                                          *
*Entries :    12429 : Total  Size=     191692 bytes  File Size  =      54345 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.52     *
*............................................................................*
*Br  344 :diMuVtz   : vector<float>                                          *
*Entries :    12429 : Total  Size=     191692 bytes  File Size  =      54834 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.49     *
*............................................................................*
*Br  345 :diMuChi2  : vector<float>                                          *
*Entries :    12429 : Total  Size=     191704 bytes  File Size  =      55038 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.47     *
*............................................................................*
*Br  346 :diMuNDF   : vector<int>                                            *
*Entries :    12429 : Total  Size=     191692 bytes  File Size  =      34517 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.54     *
*............................................................................*
*Br  347 :diMuVtxProb : vector<float>                                        *
*Entries :    12429 : Total  Size=     191740 bytes  File Size  =      50996 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.75     *
*............................................................................*
*Br  348 :diMu_CosAlpha : vector<float>                                      *
*Entries :    12429 : Total  Size=     191764 bytes  File Size  =      53521 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.57     *
*............................................................................*
*Br  349 :diMu_Lxy  : vector<float>                                          *
*Entries :    12429 : Total  Size=     191704 bytes  File Size  =      55340 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.46     *
*............................................................................*
*Br  350 :diMu_Rxy  : vector<float>                                          *
*Entries :    12429 : Total  Size=     191704 bytes  File Size  =      54709 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.50     *
*............................................................................*
*Br  351 :diMu_eLxy : vector<float>                                          *
*Entries :    12429 : Total  Size=     191716 bytes  File Size  =      54236 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.53     *
*............................................................................*
*Br  352 :diMu_SLxy : vector<float>                                          *
*Entries :    12429 : Total  Size=     191716 bytes  File Size  =      54715 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.49     *
*............................................................................*
*Br  353 :diMu_ctau : vector<float>                                          *
*Entries :    12429 : Total  Size=     191716 bytes  File Size  =      55450 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.45     *
*............................................................................*
*Br  354 :diMu_ctauErr : vector<float>                                       *
*Entries :    12429 : Total  Size=     191752 bytes  File Size  =      54067 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   3.54     *
*............................................................................*
*Br  355 :diMuVtxIsValid_KinFit : vector<int>                                *
*Entries :    12429 : Total  Size=     191860 bytes  File Size  =      35258 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.43     *
*............................................................................*
*Br  356 :diMuMass_KinFit : vector<float>                                    *
*Entries :    12429 : Total  Size=     176408 bytes  File Size  =      23780 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   7.40     *
*............................................................................*
*Br  357 :diMuPx_KinFit : vector<float>                                      *
*Entries :    12429 : Total  Size=     191764 bytes  File Size  =      36361 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  358 :diMuPy_KinFit : vector<float>                                      *
*Entries :    12429 : Total  Size=     191764 bytes  File Size  =      36372 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  359 :diMuPz_KinFit : vector<float>                                      *
*Entries :    12429 : Total  Size=     191764 bytes  File Size  =      36359 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  360 :diMuVtx_KinFit : vector<float>                                     *
*Entries :    12429 : Total  Size=     191776 bytes  File Size  =      36145 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.29     *
*............................................................................*
*Br  361 :diMuVty_KinFit : vector<float>                                     *
*Entries :    12429 : Total  Size=     191776 bytes  File Size  =      36186 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.29     *
*............................................................................*
*Br  362 :diMuVtz_KinFit : vector<float>                                     *
*Entries :    12429 : Total  Size=     191776 bytes  File Size  =      36180 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.29     *
*............................................................................*
*Br  363 :diMuChi2_KinFit : vector<float>                                    *
*Entries :    12429 : Total  Size=     191788 bytes  File Size  =      36374 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  364 :diMuNDF_KinFit : vector<int>                                       *
*Entries :    12429 : Total  Size=     191776 bytes  File Size  =      34955 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.47     *
*............................................................................*
*Br  365 :diMuVtxProb_KinFit : vector<float>                                 *
*Entries :    12429 : Total  Size=     191824 bytes  File Size  =      35819 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.34     *
*............................................................................*
*Br  366 :diMuCosAlpha_KinFit : vector<float>                                *
*Entries :    12429 : Total  Size=     191836 bytes  File Size  =      36285 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.27     *
*............................................................................*
*Br  367 :diMu_Lxy_KinFit : vector<float>                                    *
*Entries :    12429 : Total  Size=     191788 bytes  File Size  =      36397 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  368 :diMu_Rxy_KinFit : vector<float>                                    *
*Entries :    12429 : Total  Size=     191788 bytes  File Size  =      36345 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  369 :diMu_eLxy_KinFit : vector<float>                                   *
*Entries :    12429 : Total  Size=     191800 bytes  File Size  =      36164 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.29     *
*............................................................................*
*Br  370 :diMu_SLxy_KinFit : vector<float>                                   *
*Entries :    12429 : Total  Size=     191800 bytes  File Size  =      36171 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.29     *
*............................................................................*
*Br  371 :diMu_ctau_KinFit : vector<float>                                   *
*Entries :    12429 : Total  Size=     191800 bytes  File Size  =      36215 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.28     *
*............................................................................*
*Br  372 :diMu_ctauErr_KinFit : vector<float>                                *
*Entries :    12429 : Total  Size=     191836 bytes  File Size  =      36362 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  373 :diMu_Mmu1_KinFit : vector<float>                                   *
*Entries :    12429 : Total  Size=     191800 bytes  File Size  =      34917 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.48     *
*............................................................................*
*Br  374 :diMu_Mmu2_KinFit : vector<float>                                   *
*Entries :    12429 : Total  Size=     191800 bytes  File Size  =      34917 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.48     *
*............................................................................*
*Br  375 :diMu_mu1Px_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36401 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.25     *
*............................................................................*
*Br  376 :diMu_mu1Py_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36396 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  377 :diMu_mu1Pz_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36409 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.25     *
*............................................................................*
*Br  378 :diMu_mu2Px_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36422 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.25     *
*............................................................................*
*Br  379 :diMu_mu2Py_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36399 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  380 :diMu_mu2Pz_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36410 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.25     *
*............................................................................*
*Br  381 :diMu_mu1En_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36354 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  382 :diMu_mu2En_KinFit : vector<float>                                  *
*Entries :    12429 : Total  Size=     191812 bytes  File Size  =      36351 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.26     *
*............................................................................*
*Br  383 :diMu_mu1Charge_KinFit : vector<int>                                *
*Entries :    12429 : Total  Size=     191860 bytes  File Size  =      35382 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.41     *
*............................................................................*
*Br  384 :diMu_mu2Charge_KinFit : vector<int>                                *
*Entries :    12429 : Total  Size=     191860 bytes  File Size  =      35399 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   5.40     *
*............................................................................*
```

### MC
```
******************************************************************************
*Tree    :EventTree : Event data (tag V08_00_26_06)                          *
*Entries :     3277 : Total =        33610240 bytes  File  Size =   10278747 *
*        :          : Tree compression factor =   3.26                       *
******************************************************************************
*Br    0 :run       : run/I                                                  *
*Entries :     3277 : Total  Size=      13649 bytes  File Size  =        174 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  75.76     *
*............................................................................*
*Br    1 :event     : event/L                                                *
*Entries :     3277 : Total  Size=      26775 bytes  File Size  =       5855 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   4.49     *
*............................................................................*
*Br    2 :lumis     : lumis/I                                                *
*Entries :     3277 : Total  Size=      13659 bytes  File Size  =        307 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  42.95     *
*............................................................................*
*Br    3 :isData    : isData/O                                               *
*Entries :     3277 : Total  Size=       3827 bytes  File Size  =        125 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  26.84     *
*............................................................................*
*Br    4 :nVtx      : nVtx/I                                                 *
*Entries :     3277 : Total  Size=      13654 bytes  File Size  =       3797 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   3.47     *
*............................................................................*
*Br    5 :nGoodVtx  : nGoodVtx/I                                             *
*Entries :     3277 : Total  Size=      13674 bytes  File Size  =       3743 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   3.52     *
*............................................................................*
*Br    6 :nTrksPV   : nTrksPV/I                                              *
*Entries :     3277 : Total  Size=      13669 bytes  File Size  =        168 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  78.49     *
*............................................................................*
*Br    7 :isPVGood  : isPVGood/O                                             *
*Entries :     3277 : Total  Size=       3837 bytes  File Size  =        170 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  19.75     *
*............................................................................*
*Br    8 :vtx       : vtx/F                                                  *
*Entries :     3277 : Total  Size=      13649 bytes  File Size  =      10772 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.22     *
*............................................................................*
*Br    9 :vty       : vty/F                                                  *
*Entries :     3277 : Total  Size=      13649 bytes  File Size  =      10466 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.26     *
*............................................................................*
*Br   10 :vtz       : vtz/F                                                  *
*Entries :     3277 : Total  Size=      13649 bytes  File Size  =      12263 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   11 :rho       : rho/F                                                  *
*Entries :     3277 : Total  Size=      13649 bytes  File Size  =      11602 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.14     *
*............................................................................*
*Br   12 :rhoCentral : rhoCentral/F                                          *
*Entries :     3277 : Total  Size=      13684 bytes  File Size  =      11582 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.14     *
*............................................................................*
*Br   13 :HLTEleMuX : HLTEleMuX/l                                            *
*Entries :     3277 : Total  Size=      26795 bytes  File Size  =       7702 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   3.41     *
*............................................................................*
*Br   14 :HLTPho    : HLTPho/l                                               *
*Entries :     3277 : Total  Size=      26780 bytes  File Size  =       4340 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   6.06     *
*............................................................................*
*Br   15 :HLTJet    : HLTJet/l                                               *
*Entries :     3277 : Total  Size=      26780 bytes  File Size  =       2376 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  11.07     *
*............................................................................*
*Br   16 :HLTEleMuXIsPrescaled : HLTEleMuXIsPrescaled/l                      *
*Entries :     3277 : Total  Size=      26850 bytes  File Size  =        239 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression= 110.08     *
*............................................................................*
*Br   17 :HLTPhoIsPrescaled : HLTPhoIsPrescaled/l                            *
*Entries :     3277 : Total  Size=      26835 bytes  File Size  =        236 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression= 111.46     *
*............................................................................*
*Br   18 :HLTJetIsPrescaled : HLTJetIsPrescaled/l                            *
*Entries :     3277 : Total  Size=      26835 bytes  File Size  =        236 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression= 111.46     *
*............................................................................*
*Br   19 :pdf       : vector<float>                                          *
*Entries :     3277 : Total  Size=     138519 bytes  File Size  =      80117 *
*Baskets :        5 : Basket Size=      32000 bytes  Compression=   1.72     *
*............................................................................*
*Br   20 :pthat     : pthat/F                                                *
*Entries :     3277 : Total  Size=      13659 bytes  File Size  =        166 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  79.43     *
*............................................................................*
*Br   21 :processID : processID/F                                            *
*Entries :     3277 : Total  Size=      13679 bytes  File Size  =        185 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  71.29     *
*............................................................................*
*Br   22 :genWeight : genWeight/F                                            *
*Entries :     3277 : Total  Size=      13679 bytes  File Size  =        182 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  72.47     *
*............................................................................*
*Br   23 :genHT     : genHT/F                                                *
*Entries :     3277 : Total  Size=      13659 bytes  File Size  =        166 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  79.43     *
*............................................................................*
*Br   24 :genPho1   : genPho1/F                                              *
*Entries :     3277 : Total  Size=      13669 bytes  File Size  =        168 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  78.49     *
*............................................................................*
*Br   25 :genPho2   : genPho2/F                                              *
*Entries :     3277 : Total  Size=      13669 bytes  File Size  =        168 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  78.49     *
*............................................................................*
*Br   26 :EventTag  : TString                                                *
*Entries :     3277 : Total  Size=      16963 bytes  File Size  =       4741 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   3.47     *
*............................................................................*
*Br   27 :nPUInfo   : nPUInfo/I                                              *
*Entries :     3277 : Total  Size=      13669 bytes  File Size  =        179 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  73.67     *
*............................................................................*
*Br   28 :nPU       : vector<int>                                            *
*Entries :     3277 : Total  Size=     256823 bytes  File Size  =      80314 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   3.19     *
*............................................................................*
*Br   29 :puBX      : vector<int>                                            *
*Entries :     3277 : Total  Size=     256836 bytes  File Size  =      10341 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=  24.79     *
*............................................................................*
*Br   30 :puTrue    : vector<float>                                          *
*Entries :     3277 : Total  Size=     256862 bytes  File Size  =      28998 *
*Baskets :        9 : Basket Size=      32000 bytes  Compression=   8.84     *
*............................................................................*
*Br   31 :nMC       : nMC/I                                                  *
*Entries :     3277 : Total  Size=      13649 bytes  File Size  =       4235 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   3.11     *
*............................................................................*
*Br   32 :mcPID     : vector<int>                                            *
*Entries :     3277 : Total  Size=     518521 bytes  File Size  =     127739 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   4.05     *
*............................................................................*
*Br   33 :mcVtx     : vector<float>                                          *
*Entries :     3277 : Total  Size=     518521 bytes  File Size  =      76098 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   6.81     *
*............................................................................*
*Br   34 :mcVty     : vector<float>                                          *
*Entries :     3277 : Total  Size=     518521 bytes  File Size  =      75460 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   6.86     *
*............................................................................*
*Br   35 :mcVtz     : vector<float>                                          *
*Entries :     3277 : Total  Size=     518521 bytes  File Size  =      75531 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   6.86     *
*............................................................................*
*Br   36 :mcPt      : vector<float>                                          *
*Entries :     3277 : Total  Size=     518500 bytes  File Size  =     408921 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   1.27     *
*............................................................................*
*Br   37 :mcMass    : vector<float>                                          *
*Entries :     3277 : Total  Size=     518542 bytes  File Size  =     158514 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   3.27     *
*............................................................................*
*Br   38 :mcEta     : vector<float>                                          *
*Entries :     3277 : Total  Size=     518521 bytes  File Size  =     446509 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   1.16     *
*............................................................................*
*Br   39 :mcPhi     : vector<float>                                          *
*Entries :     3277 : Total  Size=     518521 bytes  File Size  =     407731 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   1.27     *
*............................................................................*
*Br   40 :mcE       : vector<float>                                          *
*Entries :     3277 : Total  Size=     518479 bytes  File Size  =     436916 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   1.19     *
*............................................................................*
*Br   41 :mcEt      : vector<float>                                          *
*Entries :     3277 : Total  Size=     518500 bytes  File Size  =     408988 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   1.27     *
*............................................................................*
*Br   42 :mcGMomPID : vector<int>                                            *
*Entries :     3277 : Total  Size=     518605 bytes  File Size  =      69764 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   7.42     *
*............................................................................*
*Br   43 :mcMomPID  : vector<int>                                            *
*Entries :     3277 : Total  Size=     518584 bytes  File Size  =      81824 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   6.33     *
*............................................................................*
*Br   44 :mcMomPt   : vector<float>                                          *
*Entries :     3277 : Total  Size=     518563 bytes  File Size  =     115393 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   4.49     *
*............................................................................*
*Br   45 :mcMomMass : vector<float>                                          *
*Entries :     3277 : Total  Size=     518605 bytes  File Size  =      96010 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   5.39     *
*............................................................................*
*Br   46 :mcMomEta  : vector<float>                                          *
*Entries :     3277 : Total  Size=     518584 bytes  File Size  =     131390 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   3.94     *
*............................................................................*
*Br   47 :mcMomPhi  : vector<float>                                          *
*Entries :     3277 : Total  Size=     518584 bytes  File Size  =     117015 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   4.43     *
*............................................................................*
*Br   48 :mcStatusFlag : vector<unsigned short>                              *
*Entries :     3277 : Total  Size=     282604 bytes  File Size  =      38680 *
*Baskets :       10 : Basket Size=      32000 bytes  Compression=   7.29     *
*............................................................................*
*Br   49 :mcParentage : vector<int>                                          *
*Entries :     3277 : Total  Size=     518647 bytes  File Size  =      63878 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   8.11     *
*............................................................................*
*Br   50 :mcStatus  : vector<int>                                            *
*Entries :     3277 : Total  Size=     518584 bytes  File Size  =      88788 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   5.83     *
*............................................................................*
*Br   51 :mcCalIsoDR03 : vector<float>                                       *
*Entries :     3277 : Total  Size=     518668 bytes  File Size  =      86232 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   6.01     *
*............................................................................*
*Br   52 :mcTrkIsoDR03 : vector<float>                                       *
*Entries :     3277 : Total  Size=     518668 bytes  File Size  =      81728 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   6.34     *
*............................................................................*
*Br   53 :mcCalIsoDR04 : vector<float>                                       *
*Entries :     3277 : Total  Size=     518668 bytes  File Size  =      93971 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   5.51     *
*............................................................................*
*Br   54 :mcTrkIsoDR04 : vector<float>                                       *
*Entries :     3277 : Total  Size=     518668 bytes  File Size  =      88022 *
*Baskets :       17 : Basket Size=      32000 bytes  Compression=   5.88     *
*............................................................................*
*Br   55 :genMET    : genMET/F                                               *
*Entries :     3277 : Total  Size=      13664 bytes  File Size  =      12546 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.05     *
*............................................................................*
*Br   56 :genMETPhi : genMETPhi/F                                            *
*Entries :     3277 : Total  Size=      13679 bytes  File Size  =      12231 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   57 :metFilters : metFilters/I                                          *
*Entries :     3277 : Total  Size=      13684 bytes  File Size  =        277 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  47.62     *
*............................................................................*
*Br   58 :pfMET     : pfMET/F                                                *
*Entries :     3277 : Total  Size=      13659 bytes  File Size  =      11839 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   59 :pfMETPhi  : pfMETPhi/F                                             *
*Entries :     3277 : Total  Size=      13674 bytes  File Size  =      12189 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   60 :pfMETsumEt : pfMETsumEt/F                                          *
*Entries :     3277 : Total  Size=      13684 bytes  File Size  =      11490 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.15     *
*............................................................................*
*Br   61 :pfMETmEtSig : pfMETmEtSig/F                                        *
*Entries :     3277 : Total  Size=      13689 bytes  File Size  =      11735 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   62 :pfMETSig  : pfMETSig/F                                             *
*Entries :     3277 : Total  Size=      13674 bytes  File Size  =      12059 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.09     *
*............................................................................*
*Br   63 :pfMET_T1JERUp : pfMET_T1JERUp/F                                    *
*Entries :     3277 : Total  Size=      13699 bytes  File Size  =      11838 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   64 :pfMET_T1JERDo : pfMET_T1JERDo/F                                    *
*Entries :     3277 : Total  Size=      13699 bytes  File Size  =      11830 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   65 :pfMET_T1JESUp : pfMET_T1JESUp/F                                    *
*Entries :     3277 : Total  Size=      13699 bytes  File Size  =      11848 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   66 :pfMET_T1JESDo : pfMET_T1JESDo/F                                    *
*Entries :     3277 : Total  Size=      13699 bytes  File Size  =      11830 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   67 :pfMET_T1UESUp : pfMET_T1UESUp/F                                    *
*Entries :     3277 : Total  Size=      13699 bytes  File Size  =      11826 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   68 :pfMET_T1UESDo : pfMET_T1UESDo/F                                    *
*Entries :     3277 : Total  Size=      13699 bytes  File Size  =      11829 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   69 :pfMETPhi_T1JESUp : pfMETPhi_T1JESUp/F                              *
*Entries :     3277 : Total  Size=      13714 bytes  File Size  =      12199 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   70 :pfMETPhi_T1JESDo : pfMETPhi_T1JESDo/F                              *
*Entries :     3277 : Total  Size=      13714 bytes  File Size  =      12196 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   71 :pfMETPhi_T1UESUp : pfMETPhi_T1UESUp/F                              *
*Entries :     3277 : Total  Size=      13714 bytes  File Size  =      12189 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   72 :pfMETPhi_T1UESDo : pfMETPhi_T1UESDo/F                              *
*Entries :     3277 : Total  Size=      13714 bytes  File Size  =      12207 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.08     *
*............................................................................*
*Br   73 :nPho      : nPho/I                                                 *
*Entries :     3277 : Total  Size=      13654 bytes  File Size  =       1320 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   9.99     *
*............................................................................*
*Br   74 :phoE      : vector<float>                                          *
*Entries :     3277 : Total  Size=      59252 bytes  File Size  =      23749 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.47     *
*............................................................................*
*Br   75 :phoEt     : vector<float>                                          *
*Entries :     3277 : Total  Size=      59259 bytes  File Size  =      23569 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.49     *
*............................................................................*
*Br   76 :phoEta    : vector<float>                                          *
*Entries :     3277 : Total  Size=      59266 bytes  File Size  =      24702 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.38     *
*............................................................................*
*Br   77 :phoPhi    : vector<float>                                          *
*Entries :     3277 : Total  Size=      59266 bytes  File Size  =      24629 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.39     *
*............................................................................*
*Br   78 :phoCalibE : vector<float>                                          *
*Entries :     3277 : Total  Size=      59287 bytes  File Size  =      23710 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.48     *
*............................................................................*
*Br   79 :phoCalibEt : vector<float>                                         *
*Entries :     3277 : Total  Size=      59294 bytes  File Size  =      23549 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.50     *
*............................................................................*
*Br   80 :phoSCE    : vector<float>                                          *
*Entries :     3277 : Total  Size=      59266 bytes  File Size  =      23748 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.48     *
*............................................................................*
*Br   81 :phoSCRawE : vector<float>                                          *
*Entries :     3277 : Total  Size=      59287 bytes  File Size  =      23749 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.48     *
*............................................................................*
*Br   82 :phoESEnP1 : vector<float>                                          *
*Entries :     3277 : Total  Size=      59287 bytes  File Size  =      14328 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.10     *
*............................................................................*
*Br   83 :phoESEnP2 : vector<float>                                          *
*Entries :     3277 : Total  Size=      59287 bytes  File Size  =      14418 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.08     *
*............................................................................*
*Br   84 :phoSCEta  : vector<float>                                          *
*Entries :     3277 : Total  Size=      59280 bytes  File Size  =      24662 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.38     *
*............................................................................*
*Br   85 :phoSCPhi  : vector<float>                                          *
*Entries :     3277 : Total  Size=      59280 bytes  File Size  =      24609 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.39     *
*............................................................................*
*Br   86 :phoSCEtaWidth : vector<float>                                      *
*Entries :     3277 : Total  Size=      59315 bytes  File Size  =      23767 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.47     *
*............................................................................*
*Br   87 :phoSCPhiWidth : vector<float>                                      *
*Entries :     3277 : Total  Size=      59315 bytes  File Size  =      23892 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br   88 :phoSCBrem : vector<float>                                          *
*Entries :     3277 : Total  Size=      59287 bytes  File Size  =      23893 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br   89 :phohasPixelSeed : vector<int>                                      *
*Entries :     3277 : Total  Size=      59329 bytes  File Size  =      10285 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.72     *
*............................................................................*
*Br   90 :phoEleVeto : vector<int>                                           *
*Entries :     3277 : Total  Size=      59294 bytes  File Size  =       9358 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   6.28     *
*............................................................................*
*Br   91 :phoR9     : vector<float>                                          *
*Entries :     3277 : Total  Size=      59259 bytes  File Size  =      22650 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.60     *
*............................................................................*
*Br   92 :phoHoverE : vector<float>                                          *
*Entries :     3277 : Total  Size=      59287 bytes  File Size  =      12555 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.68     *
*............................................................................*
*Br   93 :phoE1x3   : vector<float>                                          *
*Entries :     3277 : Total  Size=      59273 bytes  File Size  =      24011 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.45     *
*............................................................................*
*Br   94 :phoE1x5   : vector<float>                                          *
*Entries :     3277 : Total  Size=      59273 bytes  File Size  =      23983 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.45     *
*............................................................................*
*Br   95 :phoE2x2   : vector<float>                                          *
*Entries :     3277 : Total  Size=      59273 bytes  File Size  =      23928 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br   96 :phoE2x5Max : vector<float>                                         *
*Entries :     3277 : Total  Size=      59294 bytes  File Size  =      23819 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.47     *
*............................................................................*
*Br   97 :phoE5x5   : vector<float>                                          *
*Entries :     3277 : Total  Size=      59273 bytes  File Size  =      23936 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br   98 :phoESEffSigmaRR : vector<float>                                    *
*Entries :     3277 : Total  Size=      59329 bytes  File Size  =      14670 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.01     *
*............................................................................*
*Br   99 :phoSigmaIEtaIEtaFull5x5 : vector<float>                            *
*Entries :     3277 : Total  Size=      59385 bytes  File Size  =      23006 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.56     *
*............................................................................*
*Br  100 :phoSigmaIEtaIPhiFull5x5 : vector<float>                            *
*Entries :     3277 : Total  Size=      59385 bytes  File Size  =      25040 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.35     *
*............................................................................*
*Br  101 :phoSigmaIPhiIPhiFull5x5 : vector<float>                            *
*Entries :     3277 : Total  Size=      59385 bytes  File Size  =      23307 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.52     *
*............................................................................*
*Br  102 :phoE1x3Full5x5 : vector<float>                                     *
*Entries :     3277 : Total  Size=      59322 bytes  File Size  =      23951 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br  103 :phoE1x5Full5x5 : vector<float>                                     *
*Entries :     3277 : Total  Size=      59322 bytes  File Size  =      23945 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br  104 :phoE2x2Full5x5 : vector<float>                                     *
*Entries :     3277 : Total  Size=      59322 bytes  File Size  =      23923 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br  105 :phoE2x5MaxFull5x5 : vector<float>                                  *
*Entries :     3277 : Total  Size=      59343 bytes  File Size  =      23868 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br  106 :phoE5x5Full5x5 : vector<float>                                     *
*Entries :     3277 : Total  Size=      59322 bytes  File Size  =      23778 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.47     *
*............................................................................*
*Br  107 :phoR9Full5x5 : vector<float>                                       *
*Entries :     3277 : Total  Size=      59308 bytes  File Size  =      22583 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.60     *
*............................................................................*
*Br  108 :phoPFChIso : vector<float>                                         *
*Entries :     3277 : Total  Size=      59294 bytes  File Size  =      14522 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.05     *
*............................................................................*
*Br  109 :phoPFPhoIso : vector<float>                                        *
*Entries :     3277 : Total  Size=      59301 bytes  File Size  =      18692 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.15     *
*............................................................................*
*Br  110 :phoPFNeuIso : vector<float>                                        *
*Entries :     3277 : Total  Size=      59301 bytes  File Size  =      14684 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.00     *
*............................................................................*
*Br  111 :phoPFChWorstIso : vector<float>                                    *
*Entries :     3277 : Total  Size=      59329 bytes  File Size  =      19777 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.97     *
*............................................................................*
*Br  112 :phoPFRandConeChIso : vector<float>                                 *
*Entries :     3277 : Total  Size=      59350 bytes  File Size  =      14086 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.18     *
*............................................................................*
*Br  113 :phoCITKChIso : vector<float>                                       *
*Entries :     3277 : Total  Size=      59308 bytes  File Size  =      14468 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.06     *
*............................................................................*
*Br  114 :phoCITKPhoIso : vector<float>                                      *
*Entries :     3277 : Total  Size=      59315 bytes  File Size  =      18710 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.14     *
*............................................................................*
*Br  115 :phoCITKNeuIso : vector<float>                                      *
*Entries :     3277 : Total  Size=      59315 bytes  File Size  =      14704 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.00     *
*............................................................................*
*Br  116 :phoIDMVA  : vector<float>                                          *
*Entries :     3277 : Total  Size=      59280 bytes  File Size  =      23253 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.53     *
*............................................................................*
*Br  117 :phoFiredSingleTrgs : vector<ULong64_t>                             *
*Entries :     3277 : Total  Size=      71998 bytes  File Size  =      12960 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.51     *
*............................................................................*
*Br  118 :phoFiredDoubleTrgs : vector<ULong64_t>                             *
*Entries :     3277 : Total  Size=      71998 bytes  File Size  =      14786 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.83     *
*............................................................................*
*Br  119 :phoFiredL1Trgs : vector<ULong64_t>                                 *
*Entries :     3277 : Total  Size=      71970 bytes  File Size  =      18123 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.94     *
*............................................................................*
*Br  120 :phoSeedTime : vector<float>                                        *
*Entries :     3277 : Total  Size=      59301 bytes  File Size  =      24310 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.42     *
*............................................................................*
*Br  121 :phoSeedEnergy : vector<float>                                      *
*Entries :     3277 : Total  Size=      59315 bytes  File Size  =      24127 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.44     *
*............................................................................*
*Br  122 :phoxtalBits : vector<unsigned short>                               *
*Entries :     3277 : Total  Size=      52977 bytes  File Size  =       8663 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   6.06     *
*............................................................................*
*Br  123 :phoIDbit  : vector<unsigned short>                                 *
*Entries :     3277 : Total  Size=      52956 bytes  File Size  =      11504 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.56     *
*............................................................................*
*Br  124 :phoScale_stat_up : vector<float>                                   *
*Entries :     3277 : Total  Size=      59336 bytes  File Size  =      11296 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.21     *
*............................................................................*
*Br  125 :phoScale_stat_dn : vector<float>                                   *
*Entries :     3277 : Total  Size=      59336 bytes  File Size  =      11355 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.18     *
*............................................................................*
*Br  126 :phoScale_syst_up : vector<float>                                   *
*Entries :     3277 : Total  Size=      59336 bytes  File Size  =      13527 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.35     *
*............................................................................*
*Br  127 :phoScale_syst_dn : vector<float>                                   *
*Entries :     3277 : Total  Size=      59336 bytes  File Size  =      13634 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.31     *
*............................................................................*
*Br  128 :phoScale_gain_up : vector<float>                                   *
*Entries :     3277 : Total  Size=      59336 bytes  File Size  =       8799 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   6.68     *
*............................................................................*
*Br  129 :phoScale_gain_dn : vector<float>                                   *
*Entries :     3277 : Total  Size=      59336 bytes  File Size  =       8799 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   6.68     *
*............................................................................*
*Br  130 :phoResol_rho_up : vector<float>                                    *
*Entries :     3277 : Total  Size=      59329 bytes  File Size  =      13072 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.50     *
*............................................................................*
*Br  131 :phoResol_rho_dn : vector<float>                                    *
*Entries :     3277 : Total  Size=      59329 bytes  File Size  =      13066 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.50     *
*............................................................................*
*Br  132 :phoResol_phi_up : vector<float>                                    *
*Entries :     3277 : Total  Size=      59329 bytes  File Size  =      23169 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.54     *
*............................................................................*
*Br  133 :phoResol_phi_dn : vector<float>                                    *
*Entries :     3277 : Total  Size=      59329 bytes  File Size  =      23169 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.54     *
*............................................................................*
*Br  134 :nEle      : nEle/I                                                 *
*Entries :     3277 : Total  Size=      13654 bytes  File Size  =       1819 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   7.25     *
*............................................................................*
*Br  135 :eleCharge : vector<int>                                            *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      11318 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.70     *
*............................................................................*
*Br  136 :eleChargeConsistent : vector<int>                                  *
*Entries :     3277 : Total  Size=      53709 bytes  File Size  =      11467 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.64     *
*............................................................................*
*Br  137 :eleEn     : vector<float>                                          *
*Entries :     3277 : Total  Size=      53611 bytes  File Size  =      18630 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.85     *
*............................................................................*
*Br  138 :eleSCEn   : vector<float>                                          *
*Entries :     3277 : Total  Size=      53625 bytes  File Size  =      18713 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.84     *
*............................................................................*
*Br  139 :eleEcalEn : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      18037 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.95     *
*............................................................................*
*Br  140 :eleESEnP1 : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      14808 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.59     *
*............................................................................*
*Br  141 :eleESEnP2 : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      14896 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.57     *
*............................................................................*
*Br  142 :eleD0     : vector<float>                                          *
*Entries :     3277 : Total  Size=      53611 bytes  File Size  =      19139 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.78     *
*............................................................................*
*Br  143 :eleDz     : vector<float>                                          *
*Entries :     3277 : Total  Size=      53611 bytes  File Size  =      19138 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.78     *
*............................................................................*
*Br  144 :eleSIP    : vector<float>                                          *
*Entries :     3277 : Total  Size=      53618 bytes  File Size  =      18852 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.82     *
*............................................................................*
*Br  145 :elePt     : vector<float>                                          *
*Entries :     3277 : Total  Size=      53611 bytes  File Size  =      18470 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.88     *
*............................................................................*
*Br  146 :eleEta    : vector<float>                                          *
*Entries :     3277 : Total  Size=      53618 bytes  File Size  =      18784 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.83     *
*............................................................................*
*Br  147 :elePhi    : vector<float>                                          *
*Entries :     3277 : Total  Size=      53618 bytes  File Size  =      18783 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.83     *
*............................................................................*
*Br  148 :eleR9     : vector<float>                                          *
*Entries :     3277 : Total  Size=      53611 bytes  File Size  =      15034 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.53     *
*............................................................................*
*Br  149 :eleCalibPt : vector<float>                                         *
*Entries :     3277 : Total  Size=      53646 bytes  File Size  =      18477 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.88     *
*............................................................................*
*Br  150 :eleCalibEn : vector<float>                                         *
*Entries :     3277 : Total  Size=      53646 bytes  File Size  =      18598 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.86     *
*............................................................................*
*Br  151 :eleSCEta  : vector<float>                                          *
*Entries :     3277 : Total  Size=      53632 bytes  File Size  =      18775 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.83     *
*............................................................................*
*Br  152 :eleSCPhi  : vector<float>                                          *
*Entries :     3277 : Total  Size=      53632 bytes  File Size  =      18773 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.83     *
*............................................................................*
*Br  153 :eleSCRawEn : vector<float>                                         *
*Entries :     3277 : Total  Size=      53646 bytes  File Size  =      18718 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.84     *
*............................................................................*
*Br  154 :eleSCEtaWidth : vector<float>                                      *
*Entries :     3277 : Total  Size=      53667 bytes  File Size  =      18559 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.86     *
*............................................................................*
*Br  155 :eleSCPhiWidth : vector<float>                                      *
*Entries :     3277 : Total  Size=      53667 bytes  File Size  =      18550 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.87     *
*............................................................................*
*Br  156 :eleHoverE : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      13421 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.96     *
*............................................................................*
*Br  157 :eleEoverP : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      16397 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.24     *
*............................................................................*
*Br  158 :eleEoverPout : vector<float>                                       *
*Entries :     3277 : Total  Size=      53660 bytes  File Size  =      18483 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.88     *
*............................................................................*
*Br  159 :eleEoverPInv : vector<float>                                       *
*Entries :     3277 : Total  Size=      53660 bytes  File Size  =      15691 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.39     *
*............................................................................*
*Br  160 :eleBrem   : vector<float>                                          *
*Entries :     3277 : Total  Size=      53625 bytes  File Size  =      15424 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.45     *
*............................................................................*
*Br  161 :eledEtaAtVtx : vector<float>                                       *
*Entries :     3277 : Total  Size=      53660 bytes  File Size  =      17918 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.97     *
*............................................................................*
*Br  162 :eledPhiAtVtx : vector<float>                                       *
*Entries :     3277 : Total  Size=      53660 bytes  File Size  =      18353 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.90     *
*............................................................................*
*Br  163 :eledEtaAtCalo : vector<float>                                      *
*Entries :     3277 : Total  Size=      53667 bytes  File Size  =      18042 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.95     *
*............................................................................*
*Br  164 :eleSigmaIEtaIEtaFull5x5 : vector<float>                            *
*Entries :     3277 : Total  Size=      53737 bytes  File Size  =      18127 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.93     *
*............................................................................*
*Br  165 :eleSigmaIPhiIPhiFull5x5 : vector<float>                            *
*Entries :     3277 : Total  Size=      53737 bytes  File Size  =      18180 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.93     *
*............................................................................*
*Br  166 :eleConvVeto : vector<int>                                          *
*Entries :     3277 : Total  Size=      53653 bytes  File Size  =      11062 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.80     *
*............................................................................*
*Br  167 :eleMissHits : vector<int>                                          *
*Entries :     3277 : Total  Size=      53653 bytes  File Size  =      11936 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.45     *
*............................................................................*
*Br  168 :eleESEffSigmaRR : vector<float>                                    *
*Entries :     3277 : Total  Size=      53681 bytes  File Size  =      13102 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.06     *
*............................................................................*
*Br  169 :elePFChIso : vector<float>                                         *
*Entries :     3277 : Total  Size=      53646 bytes  File Size  =      13286 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.00     *
*............................................................................*
*Br  170 :elePFPhoIso : vector<float>                                        *
*Entries :     3277 : Total  Size=      53653 bytes  File Size  =      14956 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.55     *
*............................................................................*
*Br  171 :elePFNeuIso : vector<float>                                        *
*Entries :     3277 : Total  Size=      53653 bytes  File Size  =      12554 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.23     *
*............................................................................*
*Br  172 :elePFPUIso : vector<float>                                         *
*Entries :     3277 : Total  Size=      53646 bytes  File Size  =      15420 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.45     *
*............................................................................*
*Br  173 :elePFClusEcalIso : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      14883 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.57     *
*............................................................................*
*Br  174 :elePFClusHcalIso : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      14342 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.71     *
*............................................................................*
*Br  175 :elePFMiniIso : vector<float>                                       *
*Entries :     3277 : Total  Size=      53660 bytes  File Size  =      14254 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.73     *
*............................................................................*
*Br  176 :eleIDMVA  : vector<float>                                          *
*Entries :     3277 : Total  Size=      53632 bytes  File Size  =      18424 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.88     *
*............................................................................*
*Br  177 :eleIDMVAHZZ : vector<float>                                        *
*Entries :     3277 : Total  Size=      53653 bytes  File Size  =      18469 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.88     *
*............................................................................*
*Br  178 :eledEtaseedAtVtx : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      18617 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.86     *
*............................................................................*
*Br  179 :eleE1x5   : vector<float>                                          *
*Entries :     3277 : Total  Size=      53625 bytes  File Size  =      15695 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.39     *
*............................................................................*
*Br  180 :eleE2x5   : vector<float>                                          *
*Entries :     3277 : Total  Size=      53625 bytes  File Size  =      15760 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.37     *
*............................................................................*
*Br  181 :eleE5x5   : vector<float>                                          *
*Entries :     3277 : Total  Size=      53625 bytes  File Size  =      15665 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.39     *
*............................................................................*
*Br  182 :eleE1x5Full5x5 : vector<float>                                     *
*Entries :     3277 : Total  Size=      53674 bytes  File Size  =      18709 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.84     *
*............................................................................*
*Br  183 :eleE2x5Full5x5 : vector<float>                                     *
*Entries :     3277 : Total  Size=      53674 bytes  File Size  =      18671 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.85     *
*............................................................................*
*Br  184 :eleE5x5Full5x5 : vector<float>                                     *
*Entries :     3277 : Total  Size=      53674 bytes  File Size  =      18624 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.85     *
*............................................................................*
*Br  185 :eleR9Full5x5 : vector<float>                                       *
*Entries :     3277 : Total  Size=      53660 bytes  File Size  =      17723 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.00     *
*............................................................................*
*Br  186 :eleEcalDrivenSeed : vector<int>                                    *
*Entries :     3277 : Total  Size=      53695 bytes  File Size  =      11575 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.59     *
*............................................................................*
*Br  187 :eleDr03EcalRecHitSumEt : vector<float>                             *
*Entries :     3277 : Total  Size=      53730 bytes  File Size  =      15614 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.41     *
*............................................................................*
*Br  188 :eleDr03HcalDepth1TowerSumEt : vector<float>                        *
*Entries :     3277 : Total  Size=      53765 bytes  File Size  =      14689 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.62     *
*............................................................................*
*Br  189 :eleDr03HcalDepth2TowerSumEt : vector<float>                        *
*Entries :     3277 : Total  Size=      53765 bytes  File Size  =      12084 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.40     *
*............................................................................*
*Br  190 :eleDr03HcalTowerSumEt : vector<float>                              *
*Entries :     3277 : Total  Size=      53723 bytes  File Size  =      14725 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.61     *
*............................................................................*
*Br  191 :eleDr03TkSumPt : vector<float>                                     *
*Entries :     3277 : Total  Size=      53674 bytes  File Size  =      12793 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.16     *
*............................................................................*
*Br  192 :elecaloEnergy : vector<float>                                      *
*Entries :     3277 : Total  Size=      53667 bytes  File Size  =      18049 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.95     *
*............................................................................*
*Br  193 :eleTrkdxy : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      19145 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.78     *
*............................................................................*
*Br  194 :eleKFHits : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      13646 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.89     *
*............................................................................*
*Br  195 :eleKFChi2 : vector<float>                                          *
*Entries :     3277 : Total  Size=      53639 bytes  File Size  =      17370 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.06     *
*............................................................................*
*Br  196 :eleGSFChi2 : vector<float>                                         *
*Entries :     3277 : Total  Size=      53646 bytes  File Size  =      18506 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.87     *
*............................................................................*
*Br  197 :eleGSFPt  : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46544 bytes  File Size  =       5142 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.96     *
*............................................................................*
*Br  198 :eleGSFEta : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46550 bytes  File Size  =       5155 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.93     *
*............................................................................*
*Br  199 :eleGSFPhi : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46550 bytes  File Size  =       5155 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.93     *
*............................................................................*
*Br  200 :eleGSFCharge : vector<vector<float> >                              *
*Entries :     3277 : Total  Size=      46568 bytes  File Size  =       5178 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.90     *
*............................................................................*
*Br  201 :eleGSFHits : vector<vector<int> >                                  *
*Entries :     3277 : Total  Size=      46556 bytes  File Size  =       5151 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  202 :eleGSFMissHits : vector<vector<int> >                              *
*Entries :     3277 : Total  Size=      46580 bytes  File Size  =       5165 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.92     *
*............................................................................*
*Br  203 :eleGSFNHitsMax : vector<vector<int> >                              *
*Entries :     3277 : Total  Size=      46580 bytes  File Size  =       5165 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.92     *
*............................................................................*
*Br  204 :eleGSFVtxProb : vector<vector<float> >                             *
*Entries :     3277 : Total  Size=      46574 bytes  File Size  =       5168 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.91     *
*............................................................................*
*Br  205 :eleGSFlxyPV : vector<vector<float> >                               *
*Entries :     3277 : Total  Size=      46562 bytes  File Size  =       5153 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  206 :eleGSFlxyBS : vector<vector<float> >                               *
*Entries :     3277 : Total  Size=      46562 bytes  File Size  =       5153 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  207 :eleBCEn   : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46538 bytes  File Size  =       5145 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  208 :eleBCEta  : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46544 bytes  File Size  =       5142 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.96     *
*............................................................................*
*Br  209 :eleBCPhi  : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46544 bytes  File Size  =       5142 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.96     *
*............................................................................*
*Br  210 :eleBCS25  : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46544 bytes  File Size  =       5142 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.96     *
*............................................................................*
*Br  211 :eleBCS15  : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46544 bytes  File Size  =       5142 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.96     *
*............................................................................*
*Br  212 :eleBCSieie : vector<vector<float> >                                *
*Entries :     3277 : Total  Size=      46556 bytes  File Size  =       5151 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  213 :eleBCSieip : vector<vector<float> >                                *
*Entries :     3277 : Total  Size=      46556 bytes  File Size  =       5151 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  214 :eleBCSipip : vector<vector<float> >                                *
*Entries :     3277 : Total  Size=      46556 bytes  File Size  =       5151 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.94     *
*............................................................................*
*Br  215 :eleFiredSingleTrgs : vector<ULong64_t>                             *
*Entries :     3277 : Total  Size=      60702 bytes  File Size  =      12324 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.88     *
*............................................................................*
*Br  216 :eleFiredDoubleTrgs : vector<ULong64_t>                             *
*Entries :     3277 : Total  Size=      60702 bytes  File Size  =      10493 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.73     *
*............................................................................*
*Br  217 :eleFiredL1Trgs : vector<ULong64_t>                                 *
*Entries :     3277 : Total  Size=      60674 bytes  File Size  =      13757 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.37     *
*............................................................................*
*Br  218 :eleIDbit  : vector<unsigned short>                                 *
*Entries :     3277 : Total  Size=      50044 bytes  File Size  =      11192 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.43     *
*............................................................................*
*Br  219 :eleScale_stat_up : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      11685 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.55     *
*............................................................................*
*Br  220 :eleScale_stat_dn : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      11847 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.49     *
*............................................................................*
*Br  221 :eleScale_syst_up : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      12755 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.17     *
*............................................................................*
*Br  222 :eleScale_syst_dn : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      12790 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.16     *
*............................................................................*
*Br  223 :eleScale_gain_up : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      10247 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.19     *
*............................................................................*
*Br  224 :eleScale_gain_dn : vector<float>                                   *
*Entries :     3277 : Total  Size=      53688 bytes  File Size  =      10247 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.19     *
*............................................................................*
*Br  225 :eleResol_rho_up : vector<float>                                    *
*Entries :     3277 : Total  Size=      53681 bytes  File Size  =      13000 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.09     *
*............................................................................*
*Br  226 :eleResol_rho_dn : vector<float>                                    *
*Entries :     3277 : Total  Size=      53681 bytes  File Size  =      13021 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.08     *
*............................................................................*
*Br  227 :eleResol_phi_up : vector<float>                                    *
*Entries :     3277 : Total  Size=      53681 bytes  File Size  =      17608 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.02     *
*............................................................................*
*Br  228 :eleResol_phi_dn : vector<float>                                    *
*Entries :     3277 : Total  Size=      53681 bytes  File Size  =      17608 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.02     *
*............................................................................*
*Br  229 :nMu       : nMu/I                                                  *
*Entries :     3277 : Total  Size=      13649 bytes  File Size  =       1844 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   7.15     *
*............................................................................*
*Br  230 :muPt      : vector<float>                                          *
*Entries :     3277 : Total  Size=      76428 bytes  File Size  =      41549 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.83     *
*............................................................................*
*Br  231 :muEn      : vector<float>                                          *
*Entries :     3277 : Total  Size=      76428 bytes  File Size  =      41703 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.82     *
*............................................................................*
*Br  232 :muEta     : vector<float>                                          *
*Entries :     3277 : Total  Size=      76435 bytes  File Size  =      42861 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.77     *
*............................................................................*
*Br  233 :muPhi     : vector<float>                                          *
*Entries :     3277 : Total  Size=      76435 bytes  File Size  =      42879 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.77     *
*............................................................................*
*Br  234 :muCharge  : vector<int>                                            *
*Entries :     3277 : Total  Size=      76456 bytes  File Size  =      15114 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   5.03     *
*............................................................................*
*Br  235 :muType    : vector<int>                                            *
*Entries :     3277 : Total  Size=      76442 bytes  File Size  =      16539 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.59     *
*............................................................................*
*Br  236 :muIDbit   : vector<unsigned short>                                 *
*Entries :     3277 : Total  Size=      61537 bytes  File Size  =      14001 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.36     *
*............................................................................*
*Br  237 :muD0      : vector<float>                                          *
*Entries :     3277 : Total  Size=      76428 bytes  File Size  =      43651 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.74     *
*............................................................................*
*Br  238 :muDz      : vector<float>                                          *
*Entries :     3277 : Total  Size=      76428 bytes  File Size  =      43742 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.74     *
*............................................................................*
*Br  239 :muSIP     : vector<float>                                          *
*Entries :     3277 : Total  Size=      76435 bytes  File Size  =      42491 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.79     *
*............................................................................*
*Br  240 :muChi2NDF : vector<float>                                          *
*Entries :     3277 : Total  Size=      76463 bytes  File Size  =      38312 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.98     *
*............................................................................*
*Br  241 :muInnerD0 : vector<float>                                          *
*Entries :     3277 : Total  Size=      76463 bytes  File Size  =      43685 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.74     *
*............................................................................*
*Br  242 :muInnerDz : vector<float>                                          *
*Entries :     3277 : Total  Size=      76463 bytes  File Size  =      43784 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.74     *
*............................................................................*
*Br  243 :muTrkLayers : vector<int>                                          *
*Entries :     3277 : Total  Size=      76477 bytes  File Size  =      20363 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.73     *
*............................................................................*
*Br  244 :muPixelLayers : vector<int>                                        *
*Entries :     3277 : Total  Size=      76491 bytes  File Size  =      16540 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.59     *
*............................................................................*
*Br  245 :muPixelHits : vector<int>                                          *
*Entries :     3277 : Total  Size=      76477 bytes  File Size  =      17747 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.28     *
*............................................................................*
*Br  246 :muMuonHits : vector<int>                                           *
*Entries :     3277 : Total  Size=      76470 bytes  File Size  =      23453 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.24     *
*............................................................................*
*Br  247 :muStations : vector<int>                                           *
*Entries :     3277 : Total  Size=      76470 bytes  File Size  =      18310 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.15     *
*............................................................................*
*Br  248 :muMatches : vector<int>                                            *
*Entries :     3277 : Total  Size=      76463 bytes  File Size  =      18831 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.03     *
*............................................................................*
*Br  249 :muTrkQuality : vector<int>                                         *
*Entries :     3277 : Total  Size=      76484 bytes  File Size  =      12656 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   6.00     *
*............................................................................*
*Br  250 :muIsoTrk  : vector<float>                                          *
*Entries :     3277 : Total  Size=      76456 bytes  File Size  =      28237 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.69     *
*............................................................................*
*Br  251 :muPFChIso : vector<float>                                          *
*Entries :     3277 : Total  Size=      76463 bytes  File Size  =      30372 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.50     *
*............................................................................*
*Br  252 :muPFPhoIso : vector<float>                                         *
*Entries :     3277 : Total  Size=      76470 bytes  File Size  =      32805 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.32     *
*............................................................................*
*Br  253 :muPFNeuIso : vector<float>                                         *
*Entries :     3277 : Total  Size=      76470 bytes  File Size  =      28507 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.67     *
*............................................................................*
*Br  254 :muPFPUIso : vector<float>                                          *
*Entries :     3277 : Total  Size=      76463 bytes  File Size  =      41178 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.84     *
*............................................................................*
*Br  255 :muPFChIso03 : vector<float>                                        *
*Entries :     3277 : Total  Size=      76477 bytes  File Size  =      26661 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.85     *
*............................................................................*
*Br  256 :muPFPhoIso03 : vector<float>                                       *
*Entries :     3277 : Total  Size=      76484 bytes  File Size  =      28148 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.70     *
*............................................................................*
*Br  257 :muPFNeuIso03 : vector<float>                                       *
*Entries :     3277 : Total  Size=      76484 bytes  File Size  =      23266 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.27     *
*............................................................................*
*Br  258 :muPFPUIso03 : vector<float>                                        *
*Entries :     3277 : Total  Size=      76477 bytes  File Size  =      39889 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.90     *
*............................................................................*
*Br  259 :muPFMiniIso : vector<float>                                        *
*Entries :     3277 : Total  Size=      76477 bytes  File Size  =      26435 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.87     *
*............................................................................*
*Br  260 :muFiredTrgs : vector<ULong64_t>                                    *
*Entries :     3277 : Total  Size=     106392 bytes  File Size  =      26682 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   3.97     *
*............................................................................*
*Br  261 :muFiredL1Trgs : vector<ULong64_t>                                  *
*Entries :     3277 : Total  Size=     106408 bytes  File Size  =      24588 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   4.31     *
*............................................................................*
*Br  262 :muInnervalidFraction : vector<float>                               *
*Entries :     3277 : Total  Size=      76540 bytes  File Size  =      19241 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   3.95     *
*............................................................................*
*Br  263 :musegmentCompatibility : vector<float>                             *
*Entries :     3277 : Total  Size=      76554 bytes  File Size  =      36795 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   2.07     *
*............................................................................*
*Br  264 :muchi2LocalPosition : vector<float>                                *
*Entries :     3277 : Total  Size=      76533 bytes  File Size  =      40467 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.88     *
*............................................................................*
*Br  265 :mutrkKink : vector<float>                                          *
*Entries :     3277 : Total  Size=      76463 bytes  File Size  =      41192 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.84     *
*............................................................................*
*Br  266 :muBestTrkPtError : vector<float>                                   *
*Entries :     3277 : Total  Size=      76512 bytes  File Size  =      41833 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.82     *
*............................................................................*
*Br  267 :muBestTrkPt : vector<float>                                        *
*Entries :     3277 : Total  Size=      76477 bytes  File Size  =      41633 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.82     *
*............................................................................*
*Br  268 :muBestTrkType : vector<int>                                        *
*Entries :     3277 : Total  Size=      76491 bytes  File Size  =      11350 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   6.69     *
*............................................................................*
*Br  269 :npfPho    : npfPho/I                                               *
*Entries :     3277 : Total  Size=      13664 bytes  File Size  =        167 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  78.96     *
*............................................................................*
*Br  270 :pfphoEt   : vector<float>                                          *
*Entries :     3277 : Total  Size=      46538 bytes  File Size  =       5145 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.95     *
*............................................................................*
*Br  271 :pfphoEta  : vector<float>                                          *
*Entries :     3277 : Total  Size=      46544 bytes  File Size  =       5142 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.96     *
*............................................................................*
*Br  272 :pfphoPhi  : vector<float>                                          *
*Entries :     3277 : Total  Size=      46544 bytes  File Size  =       5142 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.96     *
*............................................................................*
*Br  273 :npfHF     : npfHF/I                                                *
*Entries :     3277 : Total  Size=      13659 bytes  File Size  =        934 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  14.12     *
*............................................................................*
*Br  274 :pfHFEn    : vector<float>                                          *
*Entries :     3277 : Total  Size=      47860 bytes  File Size  =       9504 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.99     *
*............................................................................*
*Br  275 :pfHFECALEn : vector<float>                                         *
*Entries :     3277 : Total  Size=      47884 bytes  File Size  =       9513 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.98     *
*............................................................................*
*Br  276 :pfHFHCALEn : vector<float>                                         *
*Entries :     3277 : Total  Size=      47884 bytes  File Size  =       9470 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   5.00     *
*............................................................................*
*Br  277 :pfHFPt    : vector<float>                                          *
*Entries :     3277 : Total  Size=      47860 bytes  File Size  =       9059 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   5.23     *
*............................................................................*
*Br  278 :pfHFEta   : vector<float>                                          *
*Entries :     3277 : Total  Size=      47866 bytes  File Size  =       9420 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   5.03     *
*............................................................................*
*Br  279 :pfHFPhi   : vector<float>                                          *
*Entries :     3277 : Total  Size=      47866 bytes  File Size  =       9611 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.93     *
*............................................................................*
*Br  280 :pfHFIso   : vector<float>                                          *
*Entries :     3277 : Total  Size=      47866 bytes  File Size  =       9531 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.97     *
*............................................................................*
*Br  281 :nJet      : nJet/I                                                 *
*Entries :     3277 : Total  Size=      13654 bytes  File Size  =       2263 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   5.83     *
*............................................................................*
*Br  282 :jetPt     : vector<float>                                          *
*Entries :     3277 : Total  Size=      84900 bytes  File Size  =      49820 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.69     *
*............................................................................*
*Br  283 :jetEn     : vector<float>                                          *
*Entries :     3277 : Total  Size=      84900 bytes  File Size  =      50931 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.66     *
*............................................................................*
*Br  284 :jetEta    : vector<float>                                          *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      52386 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.61     *
*............................................................................*
*Br  285 :jetPhi    : vector<float>                                          *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      52203 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.62     *
*............................................................................*
*Br  286 :jetRawPt  : vector<float>                                          *
*Entries :     3277 : Total  Size=      84924 bytes  File Size  =      49945 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.69     *
*............................................................................*
*Br  287 :jetRawEn  : vector<float>                                          *
*Entries :     3277 : Total  Size=      84924 bytes  File Size  =      50946 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.66     *
*............................................................................*
*Br  288 :jetMt     : vector<float>                                          *
*Entries :     3277 : Total  Size=      84900 bytes  File Size  =      49866 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.69     *
*............................................................................*
*Br  289 :jetArea   : vector<float>                                          *
*Entries :     3277 : Total  Size=      84916 bytes  File Size  =      27420 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   3.08     *
*............................................................................*
*Br  290 :jetLeadTrackPt : vector<float>                                     *
*Entries :     3277 : Total  Size=      84972 bytes  File Size  =      37801 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.23     *
*............................................................................*
*Br  291 :jetLeadTrackEta : vector<float>                                    *
*Entries :     3277 : Total  Size=      84980 bytes  File Size  =      46693 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.81     *
*............................................................................*
*Br  292 :jetLeadTrackPhi : vector<float>                                    *
*Entries :     3277 : Total  Size=      84980 bytes  File Size  =      47983 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.76     *
*............................................................................*
*Br  293 :jetLepTrackPID : vector<int>                                       *
*Entries :     3277 : Total  Size=      84972 bytes  File Size  =      20931 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   4.04     *
*............................................................................*
*Br  294 :jetLepTrackPt : vector<float>                                      *
*Entries :     3277 : Total  Size=      84964 bytes  File Size  =      28945 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.92     *
*............................................................................*
*Br  295 :jetLepTrackEta : vector<float>                                     *
*Entries :     3277 : Total  Size=      84972 bytes  File Size  =      34837 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.42     *
*............................................................................*
*Br  296 :jetLepTrackPhi : vector<float>                                     *
*Entries :     3277 : Total  Size=      84972 bytes  File Size  =      35170 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.40     *
*............................................................................*
*Br  297 :jetCSV2BJetTags : vector<float>                                    *
*Entries :     3277 : Total  Size=      84980 bytes  File Size  =      40392 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.09     *
*............................................................................*
*Br  298 :jetJetProbabilityBJetTags : vector<float>                          *
*Entries :     3277 : Total  Size=      85060 bytes  File Size  =      32172 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.63     *
*............................................................................*
*Br  299 :jetpfCombinedMVAV2BJetTags : vector<float>                         *
*Entries :     3277 : Total  Size=      85068 bytes  File Size  =      35579 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.38     *
*............................................................................*
*Br  300 :jetDeepCSVTags_b : vector<float>                                   *
*Entries :     3277 : Total  Size=      84988 bytes  File Size  =      34857 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.42     *
*............................................................................*
*Br  301 :jetDeepCSVTags_bb : vector<float>                                  *
*Entries :     3277 : Total  Size=      84996 bytes  File Size  =      35186 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.40     *
*............................................................................*
*Br  302 :jetDeepCSVTags_c : vector<float>                                   *
*Entries :     3277 : Total  Size=      84988 bytes  File Size  =      34686 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.44     *
*............................................................................*
*Br  303 :jetDeepCSVTags_cc : vector<float>                                  *
*Entries :     3277 : Total  Size=      84996 bytes  File Size  =      34803 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.43     *
*............................................................................*
*Br  304 :jetDeepCSVTags_udsg : vector<float>                                *
*Entries :     3277 : Total  Size=      85012 bytes  File Size  =      34310 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.46     *
*............................................................................*
*Br  305 :jetPartonID : vector<int>                                          *
*Entries :     3277 : Total  Size=      84948 bytes  File Size  =      19433 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   4.35     *
*............................................................................*
*Br  306 :jetHadFlvr : vector<int>                                           *
*Entries :     3277 : Total  Size=      84940 bytes  File Size  =      19113 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   4.42     *
*............................................................................*
*Br  307 :jetGenJetEn : vector<float>                                        *
*Entries :     3277 : Total  Size=      84948 bytes  File Size  =      46689 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.81     *
*............................................................................*
*Br  308 :jetGenJetPt : vector<float>                                        *
*Entries :     3277 : Total  Size=      84948 bytes  File Size  =      46069 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.83     *
*............................................................................*
*Br  309 :jetGenJetEta : vector<float>                                       *
*Entries :     3277 : Total  Size=      84956 bytes  File Size  =      48084 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.76     *
*............................................................................*
*Br  310 :jetGenJetPhi : vector<float>                                       *
*Entries :     3277 : Total  Size=      84956 bytes  File Size  =      47968 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.76     *
*............................................................................*
*Br  311 :jetGenPartonID : vector<int>                                       *
*Entries :     3277 : Total  Size=      84972 bytes  File Size  =      12428 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   6.80     *
*............................................................................*
*Br  312 :jetGenEn  : vector<float>                                          *
*Entries :     3277 : Total  Size=      84924 bytes  File Size  =      12445 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   6.78     *
*............................................................................*
*Br  313 :jetGenPt  : vector<float>                                          *
*Entries :     3277 : Total  Size=      84924 bytes  File Size  =      12445 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   6.78     *
*............................................................................*
*Br  314 :jetGenEta : vector<float>                                          *
*Entries :     3277 : Total  Size=      84932 bytes  File Size  =      12472 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   6.77     *
*............................................................................*
*Br  315 :jetGenPhi : vector<float>                                          *
*Entries :     3277 : Total  Size=      84932 bytes  File Size  =      12472 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   6.77     *
*............................................................................*
*Br  316 :jetGenPartonMomID : vector<int>                                    *
*Entries :     3277 : Total  Size=      84996 bytes  File Size  =      12480 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   6.77     *
*............................................................................*
*Br  317 :jetP4Smear : vector<float>                                         *
*Entries :     3277 : Total  Size=      84940 bytes  File Size  =      46537 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.81     *
*............................................................................*
*Br  318 :jetP4SmearUp : vector<float>                                       *
*Entries :     3277 : Total  Size=      84956 bytes  File Size  =      47235 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.79     *
*............................................................................*
*Br  319 :jetP4SmearDo : vector<float>                                       *
*Entries :     3277 : Total  Size=      84956 bytes  File Size  =      46090 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.83     *
*............................................................................*
*Br  320 :jetPFLooseId : vector<bool>                                        *
*Entries :     3277 : Total  Size=      56211 bytes  File Size  =      13315 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   4.18     *
*............................................................................*
*Br  321 :jetID     : vector<int>                                            *
*Entries :     3277 : Total  Size=      84900 bytes  File Size  =      18678 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   4.52     *
*............................................................................*
*Br  322 :jetPUID   : vector<float>                                          *
*Entries :     3277 : Total  Size=      84916 bytes  File Size  =      51486 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.64     *
*............................................................................*
*Br  323 :jetPUFullID : vector<int>                                          *
*Entries :     3277 : Total  Size=      84948 bytes  File Size  =      20939 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   4.03     *
*............................................................................*
*Br  324 :jetJECUnc : vector<float>                                          *
*Entries :     3277 : Total  Size=      84932 bytes  File Size  =      49943 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.69     *
*............................................................................*
*Br  325 :jetFiredTrgs : vector<ULong64_t>                                   *
*Entries :     3277 : Total  Size=     123252 bytes  File Size  =      22241 *
*Baskets :        5 : Basket Size=      32000 bytes  Compression=   5.52     *
*............................................................................*
*Br  326 :jetCHF    : vector<float>                                          *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      47482 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.78     *
*............................................................................*
*Br  327 :jetNHF    : vector<float>                                          *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      41934 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.01     *
*............................................................................*
*Br  328 :jetCEF    : vector<float>                                          *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      16776 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   5.03     *
*............................................................................*
*Br  329 :jetNEF    : vector<float>                                          *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      49688 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.70     *
*............................................................................*
*Br  330 :jetNCH    : vector<int>                                            *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      25275 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   3.34     *
*............................................................................*
*Br  331 :jetNNP    : vector<int>                                            *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      25845 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   3.27     *
*............................................................................*
*Br  332 :jetMUF    : vector<float>                                          *
*Entries :     3277 : Total  Size=      84908 bytes  File Size  =      33767 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   2.50     *
*............................................................................*
*Br  333 :jetVtxPt  : vector<float>                                          *
*Entries :     3277 : Total  Size=      84924 bytes  File Size  =      16303 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   5.18     *
*............................................................................*
*Br  334 :jetVtxMass : vector<float>                                         *
*Entries :     3277 : Total  Size=      84940 bytes  File Size  =      16281 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   5.19     *
*............................................................................*
*Br  335 :jetVtxNtrks : vector<float>                                        *
*Entries :     3277 : Total  Size=      84948 bytes  File Size  =      15012 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   5.63     *
*............................................................................*
*Br  336 :jetVtx3DVal : vector<float>                                        *
*Entries :     3277 : Total  Size=      84948 bytes  File Size  =      16377 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   5.16     *
*............................................................................*
*Br  337 :jetVtx3DSig : vector<float>                                        *
*Entries :     3277 : Total  Size=      84948 bytes  File Size  =      16336 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   5.17     *
*............................................................................*
*Br  338 :nAK8Jet   : nAK8Jet/I                                              *
*Entries :     3277 : Total  Size=      13669 bytes  File Size  =        209 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  63.10     *
*............................................................................*
*Br  339 :AK8JetPt  : vector<float>                                          *
*Entries :     3277 : Total  Size=      46592 bytes  File Size  =       5373 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  340 :AK8JetEn  : vector<float>                                          *
*Entries :     3277 : Total  Size=      46592 bytes  File Size  =       5375 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  341 :AK8JetRawPt : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5373 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  342 :AK8JetRawEn : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5377 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  343 :AK8JetEta : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5374 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  344 :AK8JetPhi : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5376 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  345 :AK8JetMass : vector<float>                                         *
*Entries :     3277 : Total  Size=      46604 bytes  File Size  =       5374 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  346 :AK8Jet_tau1 : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5381 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  347 :AK8Jet_tau2 : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5380 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  348 :AK8Jet_tau3 : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5380 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  349 :AK8JetCHF : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5374 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  350 :AK8JetNHF : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5376 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  351 :AK8JetCEF : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5351 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.62     *
*............................................................................*
*Br  352 :AK8JetNEF : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5376 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  353 :AK8JetNCH : vector<int>                                            *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5348 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.62     *
*............................................................................*
*Br  354 :AK8JetNNP : vector<int>                                            *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5343 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.63     *
*............................................................................*
*Br  355 :AK8JetMUF : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5348 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.62     *
*............................................................................*
*Br  356 :AK8Jetnconstituents : vector<int>                                  *
*Entries :     3277 : Total  Size=      46658 bytes  File Size  =       5366 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.60     *
*............................................................................*
*Br  357 :AK8JetPFLooseId : vector<bool>                                     *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5690 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.10     *
*............................................................................*
*Br  358 :AK8JetPFTightLepVetoId : vector<bool>                              *
*Entries :     3277 : Total  Size=      46640 bytes  File Size  =       5788 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   7.96     *
*............................................................................*
*Br  359 :AK8JetSoftDropMass : vector<float>                                 *
*Entries :     3277 : Total  Size=      46652 bytes  File Size  =       5401 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.54     *
*............................................................................*
*Br  360 :AK8JetSoftDropMassCorr : vector<float>                             *
*Entries :     3277 : Total  Size=      46676 bytes  File Size  =       5394 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.55     *
*............................................................................*
*Br  361 :AK8JetPrunedMass : vector<float>                                   *
*Entries :     3277 : Total  Size=      46640 bytes  File Size  =       5379 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  362 :AK8JetPrunedMassCorr : vector<float>                               *
*Entries :     3277 : Total  Size=      46664 bytes  File Size  =       5395 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.55     *
*............................................................................*
*Br  363 :AK8JetpfBoostedDSVBTag : vector<float>                             *
*Entries :     3277 : Total  Size=      46676 bytes  File Size  =       5387 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.56     *
*............................................................................*
*Br  364 :AK8JetDSVnewV4 : vector<float>                                     *
*Entries :     3277 : Total  Size=      46628 bytes  File Size  =       5376 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  365 :AK8JetCSV : vector<float>                                          *
*Entries :     3277 : Total  Size=      46598 bytes  File Size  =       5361 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.60     *
*............................................................................*
*Br  366 :AK8JetJECUnc : vector<float>                                       *
*Entries :     3277 : Total  Size=      46616 bytes  File Size  =       5371 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  367 :AK8JetL2L3corr : vector<float>                                     *
*Entries :     3277 : Total  Size=      46628 bytes  File Size  =       5387 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.56     *
*............................................................................*
*Br  368 :AK8puppiPt : vector<float>                                         *
*Entries :     3277 : Total  Size=      46604 bytes  File Size  =       5383 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  369 :AK8puppiMass : vector<float>                                       *
*Entries :     3277 : Total  Size=      46616 bytes  File Size  =       5373 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  370 :AK8puppiEta : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5380 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  371 :AK8puppiPhi : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5378 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  372 :AK8puppiTau1 : vector<float>                                       *
*Entries :     3277 : Total  Size=      46616 bytes  File Size  =       5376 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  373 :AK8puppiTau2 : vector<float>                                       *
*Entries :     3277 : Total  Size=      46616 bytes  File Size  =       5379 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  374 :AK8puppiTau3 : vector<float>                                       *
*Entries :     3277 : Total  Size=      46616 bytes  File Size  =       5377 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.58     *
*............................................................................*
*Br  375 :AK8puppiSDL2L3corr : vector<float>                                 *
*Entries :     3277 : Total  Size=      46652 bytes  File Size  =       5367 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.59     *
*............................................................................*
*Br  376 :AK8puppiSDMass : vector<float>                                     *
*Entries :     3277 : Total  Size=      46628 bytes  File Size  =       5361 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.60     *
*............................................................................*
*Br  377 :AK8puppiSDMassL2L3Corr : vector<float>                             *
*Entries :     3277 : Total  Size=      46676 bytes  File Size  =       5360 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.61     *
*............................................................................*
*Br  378 :AK8JetPartonID : vector<int>                                       *
*Entries :     3277 : Total  Size=      46628 bytes  File Size  =       5366 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.59     *
*............................................................................*
*Br  379 :AK8JetHadFlvr : vector<int>                                        *
*Entries :     3277 : Total  Size=      46622 bytes  File Size  =       5332 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.65     *
*............................................................................*
*Br  380 :AK8JetGenJetIndex : vector<int>                                    *
*Entries :     3277 : Total  Size=      46646 bytes  File Size  =       5331 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.65     *
*............................................................................*
*Br  381 :AK8JetGenJetEn : vector<float>                                     *
*Entries :     3277 : Total  Size=      46628 bytes  File Size  =       5385 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.56     *
*............................................................................*
*Br  382 :AK8JetGenJetPt : vector<float>                                     *
*Entries :     3277 : Total  Size=      46628 bytes  File Size  =       5388 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.56     *
*............................................................................*
*Br  383 :AK8JetGenJetEta : vector<float>                                    *
*Entries :     3277 : Total  Size=      46634 bytes  File Size  =       5390 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.56     *
*............................................................................*
*Br  384 :AK8JetGenJetPhi : vector<float>                                    *
*Entries :     3277 : Total  Size=      46634 bytes  File Size  =       5387 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.56     *
*............................................................................*
*Br  385 :AK8JetGenPartonID : vector<int>                                    *
*Entries :     3277 : Total  Size=      46646 bytes  File Size  =       5329 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.65     *
*............................................................................*
*Br  386 :AK8JetGenEn : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5334 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.64     *
*............................................................................*
*Br  387 :AK8JetGenPt : vector<float>                                        *
*Entries :     3277 : Total  Size=      46610 bytes  File Size  =       5334 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.64     *
*............................................................................*
*Br  388 :AK8JetGenEta : vector<float>                                       *
*Entries :     3277 : Total  Size=      46616 bytes  File Size  =       5328 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.65     *
*............................................................................*
*Br  389 :AK8JetGenPhi : vector<float>                                       *
*Entries :     3277 : Total  Size=      46616 bytes  File Size  =       5328 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.65     *
*............................................................................*
*Br  390 :AK8JetGenPartonMomID : vector<int>                                 *
*Entries :     3277 : Total  Size=      46664 bytes  File Size  =       5346 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.63     *
*............................................................................*
*Br  391 :AK8JetP4Smear : vector<float>                                      *
*Entries :     3277 : Total  Size=      46622 bytes  File Size  =       5381 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  392 :AK8JetP4SmearUp : vector<float>                                    *
*Entries :     3277 : Total  Size=      46634 bytes  File Size  =       5385 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.56     *
*............................................................................*
*Br  393 :AK8JetP4SmearDo : vector<float>                                    *
*Entries :     3277 : Total  Size=      46634 bytes  File Size  =       5380 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.57     *
*............................................................................*
*Br  394 :nAK8SDSJ  : vector<int>                                            *
*Entries :     3277 : Total  Size=      46592 bytes  File Size  =       5339 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.63     *
*............................................................................*
*Br  395 :AK8SDSJPt : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46662 bytes  File Size  =       5514 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.37     *
*............................................................................*
*Br  396 :AK8SDSJEta : vector<vector<float> >                                *
*Entries :     3277 : Total  Size=      46668 bytes  File Size  =       5515 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.37     *
*............................................................................*
*Br  397 :AK8SDSJPhi : vector<vector<float> >                                *
*Entries :     3277 : Total  Size=      46668 bytes  File Size  =       5519 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.37     *
*............................................................................*
*Br  398 :AK8SDSJMass : vector<vector<float> >                               *
*Entries :     3277 : Total  Size=      46674 bytes  File Size  =       5534 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.34     *
*............................................................................*
*Br  399 :AK8SDSJE  : vector<vector<float> >                                 *
*Entries :     3277 : Total  Size=      46656 bytes  File Size  =       5517 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.37     *
*............................................................................*
*Br  400 :AK8SDSJCharge : vector<vector<int> >                               *
*Entries :     3277 : Total  Size=      46686 bytes  File Size  =       5457 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.46     *
*............................................................................*
*Br  401 :AK8SDSJFlavour : vector<vector<int> >                              *
*Entries :     3277 : Total  Size=      46692 bytes  File Size  =       5482 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.42     *
*............................................................................*
*Br  402 :AK8SDSJCSV : vector<vector<float> >                                *
*Entries :     3277 : Total  Size=      46668 bytes  File Size  =       5486 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.42     *
*............................................................................*
*Br  403 :nAK8puppiSDSJ : vector<int>                                        *
*Entries :     3277 : Total  Size=      46622 bytes  File Size  =       5348 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.62     *
*............................................................................*
*Br  404 :AK8puppiSDSJPt : vector<vector<float> >                            *
*Entries :     3277 : Total  Size=      46672 bytes  File Size  =       5505 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.38     *
*............................................................................*
*Br  405 :AK8puppiSDSJEta : vector<vector<float> >                           *
*Entries :     3277 : Total  Size=      46678 bytes  File Size  =       5512 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.37     *
*............................................................................*
*Br  406 :AK8puppiSDSJPhi : vector<vector<float> >                           *
*Entries :     3277 : Total  Size=      46678 bytes  File Size  =       5515 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.37     *
*............................................................................*
*Br  407 :AK8puppiSDSJMass : vector<vector<float> >                          *
*Entries :     3277 : Total  Size=      46684 bytes  File Size  =       5500 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.39     *
*............................................................................*
*Br  408 :AK8puppiSDSJE : vector<vector<float> >                             *
*Entries :     3277 : Total  Size=      46666 bytes  File Size  =       5507 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.38     *
*............................................................................*
*Br  409 :AK8puppiSDSJCharge : vector<vector<int> >                          *
*Entries :     3277 : Total  Size=      46696 bytes  File Size  =       5478 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.43     *
*............................................................................*
*Br  410 :AK8puppiSDSJFlavour : vector<vector<int> >                         *
*Entries :     3277 : Total  Size=      46702 bytes  File Size  =       5489 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.41     *
*............................................................................*
*Br  411 :AK8puppiSDSJCSV : vector<vector<float> >                           *
*Entries :     3277 : Total  Size=      46678 bytes  File Size  =       5496 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   8.40     *
*............................................................................*
```
