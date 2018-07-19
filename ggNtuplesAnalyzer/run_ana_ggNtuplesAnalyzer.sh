#!/bin/bash

rm -rf outputFiles
mkdir outputFiles


############################################################
# Data - 35.86/fb - 2016 Data 

#time root -b -q -l "run_Ana.C+(\"filesLists/localFiles/data.txt\", -1, \"Run2016B\")"  &
############################################################
## MC
# time root -b -q -l "run_Ana.C+(\"filesLists/filesFromEOS/ZToJPsiGamma_RunIISummer16MiniAODv2.txt\", -1, \"Run2016B\")"  &
time root -b -q -l "run_Ana.C+(\"filesLists/localFiles/mc.txt\", -1, \"MCTEST\")" &

wait

ls -lha outputFiles

hadd -f outputFiles/histos_ggNtuplesAnalyzer.root outputFiles/*.root

rm run_Ana_C.d run_Ana_C.so run_Ana_C_ACLiC_dict_rdict.pcm

echo "Done!"
