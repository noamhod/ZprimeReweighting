#include "Pythia8/Pythia.h"

#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
using namespace Pythia8;

#define _FAT(x)                    { cout << "FATAL: " << __FILE__ << " " << __LINE__ << ": " << x << endl; exit(-1); }
#define _ERR(enable,x) if(enable==1) cout << "ERROR: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _WRN(enable,x) if(enable==1) cout << "WARNING: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _DBG(enable,x) if(enable==1) cout << "DEBUG: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _INF(enable,x) if(enable==1) cout << "INFO: "  << __FILE__ << " " << __LINE__ << ": " << x << endl;

int main(int argc, char* argv[])
{
	if(argc!=2) _FAT("need to specify one arg");
	TString name = argv[1];
	_INF(1,"arg="<<name);

	gROOT->ProcessLine(".L Loader.C+"); // for the vector branches...
	
	// Generator.
	Pythia* pythia = new Pythia();

	// Pick new random number seed for each run, based on clock.
	pythia->readString("Random:setSeed = on");
	pythia->readString("Random:seed = 0");

	// Process selection.
	if(name=="DY") pythia->readString("WeakSingleBoson:ffbar2gmZ = on");
	else
	{
		pythia->readString("NewGaugeBoson:ffbar2gmZZprime = on");
		pythia->readString("Zprime:gmZmode = 0"); // 0 = full gamma^*/Z^0/Z'^0 structure, with interference included.
		if(name.Contains("MinZP"))
		{
			/// chi model
			pythia->readString("Zprime:vd = 0.783156008298");
			pythia->readString("Zprime:ad = -0.391578004149");
			pythia->readString("Zprime:vu = 0.0");
			pythia->readString("Zprime:au = 0.391578004149");
			pythia->readString("Zprime:ve = -0.783156008298");
			pythia->readString("Zprime:ae = -0.391578004149");
			pythia->readString("Zprime:vnue = -0.587367006224");
			pythia->readString("Zprime:anue = -0.587367006224");
		}
	}

	if(name.Contains("ZP"))
	{
		pythia->readString("32:m0 = 3000.");
		pythia->readString("32:mMin = 2000.");
		pythia->readString("32:mMax = 4000.");
	}
	pythia->readString("PhaseSpace:mHatMin = 2000.");
	pythia->readString("PhaseSpace:mHatMax = 4000.");


	int    idOut  = 13;
	string sidOut = "13";
	if(name=="DY")
	{
		pythia->readString("23:onMode = off"); // switch off all of the Z0 decay modes
		pythia->readString("23:onIfAny = "+sidOut); // switch on the Z0->mu-mu+ decay mode only
	}
	else
	{
		pythia->readString("32:onMode = off");
		pythia->readString("32:onIfAny = "+sidOut);
	}

	// Initialize.
	pythia->readString("Beams:eCM = 13000.");
	pythia->init();
	
	// list changes before the run
	pythia->settings.listChanged();
	pythia->particleData.listChanged();

	//ROOT
	TFile* file = new TFile(name+".tree.root","recreate");
	TTree* tree = new TTree(name,name);
	vector<TLorentzVector>* p4     = new vector<TLorentzVector>;
	vector<int>*            id     = new vector<int>;
	vector<int>*            status = new vector<int>;
	float alphaS;
	float alphaE;
	tree->Branch("p4",     &p4);
	tree->Branch("id",     &id);
	tree->Branch("status", &status);
	tree->Branch("alphaS", &alphaS);
	tree->Branch("alphaE", &alphaE);
	
	_INF(1,"Starting run");
	
	// Begin event loop. Generate event. Skip if error. List first one.
	int resonance = (name.Contains("ZP")) ? 32 : 23;
	for(int iEvent=0 ; iEvent<200000 ; ++iEvent)
	{
		if(!pythia->next()) continue;
		
		if(iEvent<1) { pythia->info.list(); pythia->process.list(); pythia->event.list(); }
		
		// branches
		p4->clear();
		id->clear();
		status->clear();
		alphaS = pythia->info.alphaS();
		alphaE = pythia->info.alphaEM();
		
		int iZ = 0;
		for(int i=0 ; i<pythia->process.size() ; i++)
		{
			if(pythia->process[i].id()==resonance)
			{
				iZ = i;
				int momA = pythia->process[i].mother1();
				int momB = pythia->process[i].mother2();
				int dauA = pythia->process[i].daughter1();
				int dauB = pythia->process[i].daughter2();

				TLorentzVector pA,pB,pX,pC,pD;
				pA.SetPxPyPzE(pythia->process[momA].px(),pythia->process[momA].py(),pythia->process[momA].pz(),pythia->process[momA].e());
				pB.SetPxPyPzE(pythia->process[momB].px(),pythia->process[momB].py(),pythia->process[momB].pz(),pythia->process[momB].e());
				pX.SetPxPyPzE(pythia->process[i].px(),   pythia->process[i].py(),   pythia->process[i].pz(),   pythia->process[i].e());
				pC.SetPxPyPzE(pythia->process[dauA].px(),pythia->process[dauA].py(),pythia->process[dauA].pz(),pythia->process[dauA].e());
				pD.SetPxPyPzE(pythia->process[dauB].px(),pythia->process[dauB].py(),pythia->process[dauB].pz(),pythia->process[dauB].e());
				p4->push_back(pA);
				p4->push_back(pB);
				p4->push_back(pX);
				p4->push_back(pC);
				p4->push_back(pD);
				id->push_back(pythia->process[momA].id());
				id->push_back(pythia->process[momB].id());
				id->push_back(pythia->process[i].id());
				id->push_back(pythia->process[dauA].id());
				id->push_back(pythia->process[dauB].id());
				status->push_back(pythia->process[momA].status());
				status->push_back(pythia->process[momB].status());
				status->push_back(pythia->process[i].status());
				status->push_back(pythia->process[dauA].status());
				status->push_back(pythia->process[dauB].status());
				//////////
				break; ///
				//////////
			}
		}
		
		///////////////////
		tree->Fill(); /////
		///////////////////
		
		if(iEvent%10000==0 && iEvent>0) cout << "Event: " << iEvent << endl;
	} // end for iEvent<prm.nEvents


	// List changes.
	pythia->settings.listChanged();
	pythia->particleData.listChanged();
	pythia->stat();

	file = tree->GetCurrentFile();
	file->cd();
	tree->Write("", TObject::kOverwrite);
	file->Write();
	file->Close();
	
	ofstream* f = new ofstream(name+".sigma.dat", ios_base::app);
	(*f) << name << "\t\t"
	<< pythia->info.nTried() << "\t\t"
	<< pythia->info.nSelected() << "\t\t"
	<< pythia->info.nAccepted() << "\t\t"
	<< pythia->info.sigmaGen() << "\t\t"
	<< pythia->info.sigmaErr() << endl;
	f->close();
	delete f;

	return 0;
}