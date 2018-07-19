#include <iostream>
using namespace std;

#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"

#if defined(CINT) && !defined(MAKECINT)
#include "plugins/getFilesList.C+"
#include "ana_ggNtuplesAnalyzer.C+"
#else
#include "plugins/getFilesList.C"
#include "ana_ggNtuplesAnalyzer.C"
#endif


Int_t run_Ana(string filesList, int nFiles, string outFileAppend) {
	gSystem->AddIncludePath("-Iexternal");
	gStyle->SetOptStat(0);

    ana_ggNtuplesAnalyzer(getFilesList(filesList) , nFiles, outFileAppend);
 	
 	return 0;
}