#ifndef STD_H
#define STD_H

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TFormula.h"
#include "TPRegexp.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TExec.h"
#include "TColor.h"
#include "THStack.h"
#include "TDirectory.h"
#include "TList.h"
#include "TF1.h"
#include "TLine.h"
#include "TRolke.h"
#include "Riostream.h"
#include "TEventList.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TVectorT.h"
#include "TBox.h"
#include "TPaletteAxis.h"
#include "TArrow.h"
#include "TCut.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TEllipse.h"
#include "TIterator.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>      // for the sprintf call
#include <string>
#include <sstream>      // for the int to string operation (stringstream call)
#include <cstring>      // for the string functions
#include <math.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include <time.h>
#include <iomanip>
#include <assert.h>
#include <cstdlib>

using namespace std;

#define _FAT(x)                    { cout << "FATAL: " << __FILE__ << " " << __LINE__ << ": " << x << endl; exit(-1); }
#define _ERR(enable,x) if(enable==1) cout << "ERROR: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _WRN(enable,x) if(enable==1) cout << "WARNING: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _DBG(enable,x) if(enable==1) cout << "DEBUG: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _INF(enable,x) if(enable==1) cout << "INFO: "  << __FILE__ << " " << __LINE__ << ": " << x << endl;

#endif
