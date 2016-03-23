/////////////////////////////////////////////
//// root -b -l -q ZprimeReweighting.C++ ////
/////////////////////////////////////////////

#include "ZprimeReweighting.h"

void plot(TH1* h1, TH1* h2, TH1* h3, TString model, TString sL, TString fname)
{
	TCanvas* cnv = new TCanvas("cnv","",600,600);
	cnv->Divide(1,2);
	TVirtualPad* tvp_hists = cnv->cd(1);
	TVirtualPad* tvp_ratio = cnv->cd(2);
	cnv->Draw();
	tvp_hists->SetPad(0.00, 0.20, 1.00, 1.00);
	tvp_ratio->SetPad(0.00, 0.02, 1.00, 0.20);
	tvp_hists->SetBottomMargin(0.012);
	tvp_ratio->SetBottomMargin(0.20);
	tvp_ratio->SetTopMargin(0.012);
	tvp_hists->SetTicks(1,1);
	tvp_ratio->SetTicks(1,1);	
	 
	TString sXtitle = h1->GetXaxis()->GetTitle();
	TH1F* hr = (TH1F*)h2->Clone("ratio");
	hr->Sumw2();
	hr->SetTitle(";"+sXtitle+";Gen/Rwt");
	hr->Divide(h3);
	
	float rmin=+0.8;
	float rmax=+1.2;
	hr->SetMarkerStyle(20);
	hr->SetMarkerSize(0.5);
	hr->SetMarkerColor(kBlack);
	hr->SetLineColor(kBlack);
	hr->SetLineStyle(1);
	hr->SetLineWidth(2);
	float xLabelSize = hr->GetXaxis()->GetLabelSize()*1.85;
	float yLabelSize = hr->GetYaxis()->GetLabelSize()*1.85;
	float xTitleSize = hr->GetXaxis()->GetTitleSize()*1.85;
	float yTitleSize = hr->GetYaxis()->GetTitleSize()*1.85;
	float titleSize  = hr->GetTitleSize()            *1.85;
	hr->GetXaxis()->SetLabelSize(xLabelSize);
	hr->GetYaxis()->SetLabelSize(yLabelSize);
	hr->GetXaxis()->SetTitleSize(xTitleSize);
	hr->GetYaxis()->SetTitleSize(yTitleSize);
	hr->SetTitleSize(titleSize);
	hr->GetYaxis()->SetTitleOffset(0.55);
	hr->GetXaxis()->SetTitleOffset(0.83);
	hr->SetMinimum(rmin);
	hr->SetMaximum(rmax);
	TLine* lineR = new TLine(hr->GetXaxis()->GetXmin(),1.,hr->GetXaxis()->GetXmax(),1.);
	
	tvp_hists->SetBottomMargin(0);
	tvp_ratio->SetTopMargin(0);
	tvp_ratio->SetBottomMargin(0.20);
	
	TLegend* leg = new TLegend(0.57,0.67,0.88,0.92,"","brNDC");
	leg->SetFillStyle(4000); // will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->AddEntry((TObject*)NULL, "#it{q#bar{q}}#rightarrow#it{#gamma}/#it{Z}/#it{Z}'_{"+model+"}#rightarrow#it{#mu}#it{#mu}  (#it{m}_{#it{Z'}} = 3 TeV)", "");
	leg->AddEntry((TObject*)NULL, "Normalised to "+sL+" (500k events)", "");
	leg->AddEntry(h2,"|DY+#it{Z}'_{"+model+"}|^{2} generated","ple");
	leg->AddEntry(h3,"|DY+#it{Z}'_{"+model+"}|^{2} reweighted","ple");
	leg->AddEntry(h1,"|DY|^{2}","ple");
	
	tvp_hists->cd();
	tvp_hists->SetLogy();
	h2->SetMinimum( h1->GetMinimum()*0.5 );
	h2->SetMaximum( (h2->GetMaximum()>h1->GetMaximum()) ? h2->GetMaximum()*3 : h1->GetMaximum()*3 );
	h2->Draw();
	h3->Draw("same");
	h1->Draw("same");
	leg->Draw("same");
	tvp_hists->Update();
	tvp_hists->RedrawAxis();
	
	tvp_ratio->cd();
	tvp_ratio->SetGridy();
	hr->Draw("e1p");
	lineR->Draw("same");
	tvp_ratio->Update();
	tvp_ratio->RedrawAxis();
	
	cnv->Update();
	cnv->SaveAs(fname);
}



void ZprimeReweighting()
{
	gROOT->ProcessLine(".L Loader.C+"); // for the vector branches...
	
	
	Int_t icol=0; // WHITE
	Int_t font=42; // Helvetica
	Double_t tsize=0.05;
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetFrameFillColor(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetStatColor(icol);
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.08);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetTitleXOffset(1.05);
	gStyle->SetTitleYOffset(0.95);
	gStyle->SetTextFont(font);
	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font,"x");
	gStyle->SetTitleFont(font,"x");
	gStyle->SetLabelFont(font,"y");
	gStyle->SetTitleFont(font,"y");
	gStyle->SetLabelFont(font,"z");
	gStyle->SetTitleFont(font,"z");
	gStyle->SetLabelSize(tsize*0.85,"x");
	gStyle->SetTitleSize(tsize*1.10,"x");
	gStyle->SetLabelSize(tsize*0.85,"y");
	gStyle->SetTitleSize(tsize*1.10,"y");
	gStyle->SetLabelSize(tsize*0.85,"z");
	gStyle->SetTitleSize(tsize*1.10,"z");
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(0.6);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	gStyle->SetEndErrorSize(0.);
	// gStyle.SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPaintTextFormat("4.1f");
	
	TFile* fDY = new TFile("DY.tree.500k.root","READ");
	TTree* tDY = (TTree*)fDY->Get("DY");
	
	TFile* fZP = new TFile("ZP.tree.500k.root","READ");
	TTree* tZP = (TTree*)fZP->Get("ZP");
	
	TFile* fMinZP = new TFile("MinZP.tree.200k.root","READ");
	// TFile* fMinZP = new TFile("ZP.tree.500k.root","READ"); // !!!!!!!!!!!!!!!
	TTree* tMinZP = (TTree*)fMinZP->Get("MinZP"); 
	// TTree* tMinZP = (TTree*)fMinZP->Get("ZP"); // !!!!!!!!!!!!!!!!!!!!
	
	float mb2fb = 1.e12;
	float nDY    = tDY->GetEntries();
	float nZP    = tZP->GetEntries();
	float nMinZP = tMinZP->GetEntries();
	float sigmaDY    = 1.56068e-13*mb2fb; // mb->fb
	float sigmaZP    = 1.32885e-12*mb2fb; // mb->fb
	float sigmaMinZP = 9.1092e-13*mb2fb; // mb->fb  !!!!!!!!!!!!!!!!!!!!
	float L = 1.; // fb-1
	TString sL = "1 fb^{-1}";
	float LDY    = nDY/sigmaDY;
	float LZP    = nZP/sigmaZP;
	float LMinZP = nMinZP/sigmaMinZP;
	float scaleDY    = L/LDY;
	float scaleZP    = L/LZP;
	float scaleMinZP = L/LMinZP;
	
	int nBins  = 50;
	float mMin = 2000;
	float mMax = 4000;
	TString binwidth = tstr((mMax-mMin)/nBins,1)+" [1/GeV]";
	
	TH1F* hDY       = new TH1F("hDY",       ";#it{m}_{#it{#mu#mu}} [GeV];Events / "+binwidth,nBins,mMin,mMax);
	TH1F* hZPgen    = new TH1F("hZPgen",    ";#it{m}_{#it{#mu#mu}} [GeV];Events / "+binwidth,nBins,mMin,mMax);
	TH1F* hZPrwt    = new TH1F("hZPrwt",    ";#it{m}_{#it{#mu#mu}} [GeV];Events / "+binwidth,nBins,mMin,mMax);
	TH1F* hMinZPgen = new TH1F("hMinZPgen", ";#it{m}_{#it{#mu#mu}} [GeV];Events / "+binwidth,nBins,mMin,mMax);
	TH1F* hMinZPrwt = new TH1F("hMinZPrwt", ";#it{m}_{#it{#mu#mu}} [GeV];Events / "+binwidth,nBins,mMin,mMax);
	
	hDY->Sumw2();
	hDY->SetLineColor(kBlack);
	hDY->SetLineWidth(2);
	hDY->SetMarkerColor(kBlack);
	hDY->SetMarkerStyle(20);
	
	hZPgen->Sumw2();
	hZPgen->SetLineColor(kRed);
	hZPgen->SetLineWidth(2);
	hZPgen->SetMarkerColor(kRed);
	hZPgen->SetMarkerStyle(20);
	
	hZPrwt->Sumw2();
	hZPrwt->SetLineColor(kAzure+9);
	hZPrwt->SetLineWidth(2);
	hZPrwt->SetMarkerColor(kAzure+9);
	hZPrwt->SetMarkerStyle(24);
	
	hMinZPgen->Sumw2();
	hMinZPgen->SetLineColor(kRed);
	hMinZPgen->SetLineWidth(2);
	hMinZPgen->SetMarkerColor(kRed);
	hMinZPgen->SetMarkerStyle(20);
	
	hMinZPrwt->Sumw2();
	hMinZPrwt->SetLineColor(kAzure+9);
	hMinZPrwt->SetLineWidth(2);
	hMinZPrwt->SetMarkerColor(kAzure+9);
	hMinZPrwt->SetMarkerStyle(24);
	
	
	
	////////////////////////////////////////
	// theoretical stuff... ////////////////
	setFermions(); /////////////////////////
	setFixedWidth(false); //////////////////
	setScaleWidth(true); ///////////////////
	resetZPmass(); /////////////////////////
	resetfgZP(); ///////////////////////////
	////////////////////////////////////////
	setZPmass(mZ0); ////////////////////////
	cout << "SM Z width is " << wTotZP() << "\n"<< endl;
	setZPmass(3000); ///////////////////////
	cout << "Z' SSM width is " << wTotZP() << "\n" << endl;
	setZPmass(3000); ///////////////////////
	/////// Z'_{#chi} parameters ///////////
	double gamma = sqrt(41./24.)*sw; ///////
	double theta = asin(-sqrt(16./41.)); ///
	setModelPrime(theta,gamma); ////////////
	cout << "Z'_chi width is " << wTotMinZP() << " (in pythia: 36.83 GeV) \n" << endl;
	for(ui2fermion::iterator it=ui2f.begin() ; it!=ui2f.end() ; ++it)
	{
		cout << "gV(" << namef(it->first) << ")=" << gMinZPV(it->first) << endl;
		cout << "gA(" << namef(it->first) << ")=" << gMinZPA(it->first) << endl;
	}
	
	////////////////////////////////////////
	// return;
	
	
	
	vector<TLorentzVector>* p4     = 0;
	vector<int>*            id     = 0;
	vector<int>*            status = 0;
	float alphaS;
	float alphaE;
	tDY->SetBranchAddress("p4",     &p4);
	tDY->SetBranchAddress("id",     &id);
	tDY->SetBranchAddress("status", &status);
	tDY->SetBranchAddress("alphaS", &alphaS);
	tDY->SetBranchAddress("alphaE", &alphaE);
	for(int i=0 ; i<nDY ; ++i)
	{
		tDY->GetEntry(i);
		TLorentzVector p = p4->at(3)+p4->at(4);
		hDY->Fill(p.M(),scaleDY);
		
		if(id->size()!=5) _ERR(1,"id->size()!=5");
		if(p4->size()!=5) _ERR(1,"p4->size()!=5");
		float m12 = (p4->at(0)+p4->at(1)).M();
		float m34 = (p4->at(3)+p4->at(4)).M();
		float m   = p4->at(2).M();
		if(fabs(m-m12)/m>0.1)      _ERR(1,"fabs(m-m12)/m>0.1: (m="+str(m)+", m12="+str(m12)+")");
		if(fabs(m-m34)/m>0.1)      _ERR(1,"fabs(m-m34)/m>0.1: (m="+str(m)+", m34="+str(m34)+")");
		if(fabs(m12-m34)/m12>0.1)  _ERR(1,"fabs(m12-m34)/m12>0.1: (m12="+str(m12)+", m34="+str(m34)+")");
		if(id->at(0)+id->at(1)!=0) _ERR(1,"id->at(0)+id->at(1)!=0: ("+str(id->at(0))+","+str(id->at(1))+")");
		if(id->at(3)+id->at(4)!=0) _ERR(1,"id->at(3)+id->at(4)!=0: ("+str(id->at(3))+","+str(id->at(4))+")");
		
		// (1) take the true quark and the true lepton and boost it to the cmf.
		// (2) calculate cos(theta*) between the tree-level lepton and the quark
		// (3) get the quark id
		// (4) get the mass of the intermediate state (Z)
		int ilepton = (id->at(3)>0) ? 3 : 4; // outgoing leptons are written in enties 3 and 4 of the vectors
		int iquark  = (id->at(0)>0) ? 0 : 1; // incoming partons are written in enties 0 and 1 of the vectors
		int iaquark = (iquark==0) ? 1 : 0;
		TLorentzVector tlvlepton, tlvqa, tlvqb;
		tlvlepton.SetPtEtaPhiE(p4->at(ilepton).Pt(), p4->at(ilepton).Eta(),p4->at(ilepton).Phi(),p4->at(ilepton).E());
		tlvqa.SetPxPyPzE(p4->at(iquark).Px(),p4->at(iquark).Py(),p4->at(iquark).Pz(),p4->at(iquark).E());
		tlvqb.SetPxPyPzE(p4->at(iaquark).Px(),p4->at(iaquark).Py(),p4->at(iaquark).Pz(),p4->at(iaquark).E());
		TLorentzVector* tlvleptonBoosted = boost(tlvqa,tlvqb,tlvlepton);
		TLorentzVector* tlvqaBoosted     = boost(tlvqa,tlvqb,tlvqa);
		TLorentzVector* tlvqbBoosted     = boost(tlvqa,tlvqb,tlvqb);
		TVector3 tv3q, tv3lep;
		tv3q.SetXYZ(tlvqaBoosted->Px(),tlvqaBoosted->Py(),tlvqaBoosted->Pz());
		tv3lep.SetXYZ(tlvleptonBoosted->Px(),tlvleptonBoosted->Py(),tlvleptonBoosted->Pz());
		float cost = tv3q.Dot(tv3lep)/(tv3q.Mag()*tv3lep.Mag()); // truth cos(theta*)
		int idIn   = (int)fabs(id->at(0));
		int idOut  = (int)fabs(id->at(3));
	
		setAlphaEM(alphaE); // !!!
		setAlphaST(alphaS); // !!!
		float sHat = m34*m34;
		minZprimeFormat = false;
		float wZPSSM = weightZP(cost,sHat,idIn,idOut);
		minZprimeFormat = true;
		float wMinZP = weightMinZP(cost,sHat,idIn,idOut);	
		hZPrwt->Fill(p.M(),scaleDY*wZPSSM);
		hMinZPrwt->Fill(p.M(),scaleDY*wMinZP);
	}
	
	p4     = 0;
	id     = 0;
	status = 0;
	tZP->SetBranchAddress("p4",     &p4);
	tZP->SetBranchAddress("id",     &id);
	tZP->SetBranchAddress("status", &status);
	tZP->SetBranchAddress("alphaS", &alphaS);
	tZP->SetBranchAddress("alphaE", &alphaE);
	for(int i=0 ; i<nZP ; ++i)
	{
		tZP->GetEntry(i);
		TLorentzVector p = p4->at(3)+p4->at(4);
		hZPgen->Fill(p.M(),scaleZP);
	}
	
	p4     = 0;
	id     = 0;
	status = 0;
	tMinZP->SetBranchAddress("p4",     &p4);
	tMinZP->SetBranchAddress("id",     &id);
	tMinZP->SetBranchAddress("status", &status);
	tMinZP->SetBranchAddress("alphaS", &alphaS);
	tMinZP->SetBranchAddress("alphaE", &alphaE);
	for(int i=0 ; i<nMinZP ; ++i)
	{
		tMinZP->GetEntry(i);
		TLorentzVector p = p4->at(3)+p4->at(4);
		hMinZPgen->Fill(p.M(),scaleMinZP);
	}
	
	plot(hDY,hZPgen,hZPrwt,"SSM",sL,"ZprimeReweighting.pdf(");
	plot(hDY,hMinZPgen,hMinZPrwt,"#chi",sL,"ZprimeReweighting.pdf)");
}