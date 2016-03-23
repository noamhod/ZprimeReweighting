#include "std.h"

//// Some enums
enum pdgFermions
{
	dQ=1, uQ=2, sQ=3, cQ=4, bQ=5, tQ=6,
	eL=11, veL=12, mL=13, vmL=14, tL=15, vtL=16,
};
enum pdgBosons
{
	GLU=21, GAMMA=22, Z=23, WPLUS=24,
	ZPRIME0=32,GRV=5000039,KK=5000023,
};



//// Some unit conversions, constants and flags
static const float MeV2TeV = 1.e-6;
static const float MeV2GeV = 1.e-3;
static const float GeV2TeV = 1.e-3;
static const float GeV2MeV = 1.e+3;
static const float TeV2MeV = 1.e+6;
static const float TeV2GeV = 1.e+3;

static const double mb2fb = 1.e12;
static const double nb2fb = 1.e+6;
static const double pb2fb = 1.e+3;
static const double nb2mb = 1.e-6;

static const float muonMass = 0.105658367; // GeV
static const float elecMass = 0.00051099891; // GeV

static const double pi   = 3.14159265;
static const double sq2 = sqrt(2.);
static const double f12 = 1./2.;
static const double f13 = 1./3.;
static const double f23 = 2./3.;
static const double f16 = 1./6.;

static const double mZ0 = 91.18760; // 91.1876;   // GeV ->PDG
static const double wZ0 = 2.50419; // 2.4952;    // GeV ->PDG
static const double mZ2 = mZ0*mZ0;   // GeV^2
static const double vev = 246.; // GeV (the SM Higgs vacuum expectation value)

static const double sw2     = 0.2312015; // ->PDG
static const double cw2     = 0.7687985; // ->PDG
static const double sw      = sqrt(sw2);
static const double cw      = sqrt(cw2);
static const double cwsw    = sqrt(sw2*cw2);
static       double alphaEM = 0.008138; //<--at 3 TeV // 0.00729735; //<-PDG // IT MAY FLOAT !!!
static       double alphaST = 0.086; //<--at 3 TeV // 0.11856; // IT MAY FLOAT !!!
static       double Qe      = sqrt(4.*pi*alphaEM);
static       void setAlphaEM(double alpha) { alphaEM = alpha; Qe = sqrt(4.*pi*alphaEM); }
static       void setAlphaST(double alpha) { alphaST = alpha; }
static       bool minZprimeFormat = false;


static const double GeV2mb  = 0.38937930419;      // 1/GeV^2 to mb (mili-barn) conversion constant
static const double GeV2nb  = 1e6*0.38937930419;  // 1/GeV^2 to nb (nano-barn) conversion constant
static const double GeV2pb  = 1e9*0.38937930419;  // 1/GeV^2 to pb (pico-barn) conversion constant
static const double GeV2fb  = 1e12*0.38937930419;  // 1/GeV^2 to fb (femto-barn) conversion constant

static bool doScale      = false;
static bool doScaleWidth = true; // turn on/off the g^2 scale in the BW width (default = true)
void setCouplingsScale(bool doscale) { doScale = doscale; }
void setScaleWidth(bool doscale)     { doScaleWidth = doscale; }

static const double mZPinit = 2000.; // GeV
static double mZP           = 2000.; // GeV
void setZPmass(double m) { mZP = m; }
void resetZPmass() { mZP = mZPinit; }

static double thetaP = 0.;
static double gammaP = 1.;
void setThetaPrime(double theta)
{
	// if(theta<0 || theta>pi) _FAT("thetaP<0 || thetaP>pi: "<<theta);
	thetaP = theta;
}
void setGammaPrime(double gamma)
{
	if(gamma<0) _FAT("gammaP<0: "<<gamma);
	gammaP = gamma;
}
void setModelPrime(double theta, double gamma)
{
	setThetaPrime(theta);
	setGammaPrime(gamma);
}

static const double       min_weight   = 1.e-30;
static const double       max_weight   = 1.e+10;
static bool               dokFactors   = false;
static bool               doFixedWidth = false;
void setkFactors(bool dokF)      { dokFactors = dokF; }
void setFixedWidth(bool doFixed) { doFixedWidth = doFixed; }




//// Some classes
class fermion
{
	public:
		fermion(string Name, unsigned int PDGid, double Mass, double Charge, double WeakI3, double y, double bl)
		{
			this->name = Name;
			this->id = PDGid;
			this->M  = Mass;
			this->Q  = Charge;
			this->I3 = WeakI3;
			this->Y  = y;
			this->BL = bl;
			if(fabs(this->id)<10) this->Nc = 3;
			else                  this->Nc = 1;
		}
		// ~fermion();

	public:
		string name;
		unsigned int id;
		unsigned int Nc;
		double M; // MeV
		double Q; // in units of the proton's charge
		double I3;
		double Y;
		double BL;
};


/// Some typedefs
typedef complex<double>             dcomplex;
typedef complex<int>                icomplex;
typedef map<unsigned int, fermion*> ui2fermion;
typedef map<string, fermion*>       s2fermion;
ui2fermion ui2f;
s2fermion  s2f;

static dcomplex Im = dcomplex(0.,1.);


//// Some general functions
string str(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return str;
}
TString tstr(float x, int prcn=-1)
{
	stringstream strm;
	string str;
	if(prcn!=-1) strm << setprecision(prcn) << fixed << x;
	else         strm                       << fixed << x;
	strm >> str;
	return (TString)str;
}
inline bool isnaninf(double x)
{
	if(std::isinf(x)) { _WRN(1,"value is infinity");     return true; }
	if(std::isnan(x)) { _WRN(1,"value is not a number"); return true; }
	return false;
}
inline TLorentzVector* boost( TLorentzVector* pBoost, TLorentzVector* p )
{
	TLorentzVector* pBoosted = (TLorentzVector*)p->Clone("");
	pBoosted->Boost(-1.*pBoost->BoostVector());
	return pBoosted;
}
inline TLorentzVector* boost( TLorentzVector pa, TLorentzVector pb, TLorentzVector p )
{
	TLorentzVector pTmp = pa+pb;
	TLorentzVector* pBoosted = (TLorentzVector*)p.Clone("");
	pBoosted->Boost(-1.*pTmp.BoostVector());
	return pBoosted;
}


void setFermions()
{	
	// FROM PYTHIA6.4 MANUAL
	////////////////////////////////////////////////////// id  mass               Q    I3    Y     B-L
	ui2f.insert( make_pair(1,  new fermion("dwn",          1,  0.33,            -f13, -f12, +f16, +f13 ) ) );
	ui2f.insert( make_pair(2,  new fermion("up",           2,  0.33,            +f23, +f12, +f16, +f13 ) ) );
	ui2f.insert( make_pair(3,  new fermion("strange",      3,  0.5,             -f13, -f12, +f16, +f13 ) ) );
	ui2f.insert( make_pair(4,  new fermion("charm",        4,  1.5,             +f23, +f12, +f16, +f13 ) ) );
	ui2f.insert( make_pair(5,  new fermion("bottom",       5,  4.8,             -f13, -f12, +f16, +f13 ) ) );
	ui2f.insert( make_pair(6,  new fermion("top",          6,  175,             +f23, +f12, +f16, +f13 ) ) );
	ui2f.insert( make_pair(11, new fermion("electron",     11, 0.511*MeV2GeV,   -1,   -f12, -f12, -1   ) ) );
	ui2f.insert( make_pair(12, new fermion("neutrino e",   12, 0.*MeV2GeV,       0,   +f12, -f12, -1   ) ) );
	ui2f.insert( make_pair(13, new fermion("muon",         13, 105.6*MeV2GeV,   -1,   -f12, -f12, -1   ) ) );
	ui2f.insert( make_pair(14, new fermion("neutrino mu",  14, 0.*MeV2GeV,       0,   +f12, -f12, -1   ) ) );
	ui2f.insert( make_pair(15, new fermion("tau",          15, 1776.82*MeV2GeV, -1,   -f12, -f12, -1   ) ) );
	ui2f.insert( make_pair(16, new fermion("neutrino tau", 16, 0.*MeV2GeV,       0,   +f12, -f12, -1   ) ) );
	
	s2f.insert( make_pair("dwn",   ui2f[1] ) );
	s2f.insert( make_pair("up",    ui2f[2] ) );
	s2f.insert( make_pair("str",   ui2f[3] ) );
	s2f.insert( make_pair("chm",   ui2f[4] ) );
	s2f.insert( make_pair("bot",   ui2f[5] ) );
	s2f.insert( make_pair("top",   ui2f[6] ) );
	s2f.insert( make_pair("elec",  ui2f[11] ) );
	s2f.insert( make_pair("nuel",  ui2f[12] ) );
	s2f.insert( make_pair("muon",  ui2f[13] ) );
	s2f.insert( make_pair("numu",  ui2f[14] ) );
	s2f.insert( make_pair("tau",   ui2f[15] ) );
	s2f.insert( make_pair("nutau", ui2f[16] ) );
}
inline string namef(unsigned int id) { return ui2f[id]->name; }
inline int    Ncf(unsigned int id)   { return ui2f[id]->Nc; }
inline double mf(unsigned int id)    { return ui2f[id]->M; }
inline double I3f(unsigned int id)   { return ui2f[id]->I3; }
inline double qf(unsigned int id)    { return ui2f[id]->Q; }
inline double Yf(unsigned int id)    { return ui2f[id]->Y; }
inline double BLf(unsigned int id)   { return ui2f[id]->BL; }



//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// COUPLINGS /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//// SM gamma, Z or Z'SSM couplings
inline double gG(unsigned int id)  { double g = qf(id);                    if(minZprimeFormat) g*=Qe; return g; }
inline double gZL(unsigned int id) { double g = (I3f(id)-qf(id)*sw2)/cwsw; if(minZprimeFormat) g*=Qe; return g; }
inline double gZR(unsigned int id) { double g = -qf(id)*sw2/cwsw;          if(minZprimeFormat) g*=Qe; return g; }
inline double gZH(unsigned int id, double h)
{
	if     (h==-f12) return gZL(id);
	else if(h==+f12) return gZR(id);
	else _FAT("unknown helicity: "<<h);
	return 0.;
}

//// Minimal Z' couplings
inline double gBL() { return gammaP*cos(thetaP)*Qe/cwsw; }
inline double gY()  { return gammaP*sin(thetaP)*Qe/cwsw; }
inline double gMinZPL(unsigned int idf)
{	
	// if      (idf==uQ  || idf==cQ  || idf==tQ)  return (f16*gY()+f13*gBL());
	// else if (idf==dQ  || idf==sQ  || idf==bQ)  return (f16*gY()+f13*gBL());
	// else if (idf==eL  || idf==mL  || idf==tL)  return (-f12*gY()-1.*gBL());
	// else if (idf==veL || idf==vmL || idf==vtL) return (-f12*gY()-1.*gBL());
	// else _FAT("unsupported fermion "<<idf);
	// return 0;

	return (Yf(idf)*gY()+BLf(idf)*gBL());
}
inline double gMinZPR(unsigned int idf)
{
	if      (idf==uQ  || idf==cQ  || idf==tQ)  return (f23*gY()+f13*gBL());
	else if (idf==dQ  || idf==sQ  || idf==bQ)  return (-f13*gY()+f13*gBL());
	else if (idf==eL  || idf==mL  || idf==tL)  return (-1.*gY()-1.*gBL());
	else if (idf==veL || idf==vmL || idf==vtL) return 0;
	else _FAT("unsupported fermion "<<idf);
	return 0;

	// if(idf==veL || idf==vmL || idf==vtL) return 0; 
	// return (Yf(idf)*gY()+BLf(idf)*gBL());
}
inline double gMinZPH(unsigned int idf, double h)
{
	if     (h==-f12) return gMinZPL(idf);
	else if(h==+f12) return gMinZPR(idf);
	else _FAT("unknown helicity: "<<h);
	return 0.;
}
inline double gMinZPV(unsigned int idf) { return (gMinZPL(idf)+gMinZPR(idf))*2*cwsw/Qe; }
inline double gMinZPA(unsigned int idf) { return (gMinZPL(idf)-gMinZPR(idf))*2*cwsw/Qe; }




//////////////////////////////////////////////////////////////////////////////
///////////////////////////// COUPLINGS SCALE ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//// ZP SSM couplings scale
static const dcomplex fgZPinit(1.,0.);
static       dcomplex fgZP(1.,0.);
inline void   setFgZP(double re, double im) { fgZP  = dcomplex(re,im); }
inline double fgZP2()                       { return real(fgZP*conj(fgZP)); }
inline void resetfgZP()                     { fgZP  = fgZPinit; }
inline dcomplex fgZPL(unsigned int idf)           { return fgZP*gZL(idf);   }
inline dcomplex fgZPR(unsigned int idf)           { return fgZP*gZR(idf);   }
inline dcomplex fgZPH(unsigned int idf, double h) { return fgZP*gZH(idf,h); }








///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// WIDTHS /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//// Z'SSM-->ffbar
inline double wZP2ffbar(unsigned int idf)
{
	double w = 0.;
	double mZP2 = mZP*mZP;
	double mf2 = mf(idf)*mf(idf);
	double gL = gZL(idf);
	double gL2 = gL*gL;
	double gR = gZR(idf);
	double gR2 = gR*gR;
	if(mZP<2.*mf(idf))      return 0.;
	if((1.-4.*mf2/mZP2)<0.) return 0.;
	w = Ncf(idf)*(alphaEM/6.)*mZP*sqrt(1.-4.*mf2/mZP2)*((gL2+gR2)+(mf2/mZP2)*(6.*gL*gR-gL2-gR2));
	if(minZprimeFormat) w *= (1./(24.*pi))/(alphaEM/6.);
	if(idf<=tQ)         w *= (1.+alphaST/pi); // QCD corrections as in Pythia
	if(doScale && doScaleWidth) w *= fgZP2();
	return w;
}
//// Minimal Z'-->ffbar
inline double wMinZP2ffbar(unsigned int idf)
{
	double w = 0.;
	double mZP2 = mZP*mZP;
	double mf2 = mf(idf)*mf(idf);
	double gL = gMinZPL(idf);
	double gL2 = gL*gL;
	double gR = gMinZPR(idf);
	double gR2 = gR*gR;
	if(mZP<2.*mf(idf))      return 0.;
	if((1.-4.*mf2/mZP2)<0.) return 0.;
	w = (Ncf(idf)/(24.*pi))*mZP*sqrt(1.-4.*mf2/mZP2)*((gL2+gR2)+(mf2/mZP2)*(6.*gL*gR-gL2-gR2));
	if(idf<=tQ) w *= (1.+alphaST/pi); // QCD corrections as in Pythia
	return w;
}
inline double wTotZP()
{
	double w = 0.;
	for(ui2fermion::iterator it=ui2f.begin() ; it!=ui2f.end() ; ++it) w += wZP2ffbar(it->first);
	// w = (doScale && doScaleWidth) ? fgZP2()*wZ0*mZP/mZ0 + wZP2ffbar(tQ) : wZ0*mZP/mZ0 + wZP2ffbar(tQ);
	return w;
}
inline double wTotMinZP()
{
	double w = 0.;
	for(ui2fermion::iterator it=ui2f.begin() ; it!=ui2f.end() ; ++it) w += wMinZP2ffbar(it->first);
	return w;
	///////////////////////////////////////
	// return 36.83; // !!!!!!!!!!!!!!!!!!!!! 3 TeV Z'_chi from Pythia
	///////////////////////////////////////
}





/////////////////////////////////////////////////////////////////
/////////////////////// AMPLITUDES //////////////////////////////
/////////////////////////////////////////////////////////////////

// Drell-Yan amplitudes
inline dcomplex hAG0(double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A = gG(idIn)*gG(idOut)/s;
	return A;
}
inline dcomplex hAZ0(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	double widthterm = (doFixedWidth) ? wZ0*mZ0 : s*(wZ0/mZ0);
	dcomplex A = gZH(idIn,hIn)*gZH(idOut,hOut)/(s-mZ2 + Im*widthterm);
	return A;
}
inline dcomplex hADY(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A = hAG0(s,idIn,idOut) + hAZ0(s,idIn,idOut,hIn,hOut);
	return A;
}
inline double hA2DY(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hADY(s,idIn,idOut,hIn,hOut);
			double angular = (1.+4.*hIn*hOut*cosTheta);
			double angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2DY(double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hADY(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}



// Z'SSM amplitudes
inline dcomplex hAZP0(double s, double w, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	double m2 = mZP*mZP;
	dcomplex gIn  = (doScale) ? fgZPH(idIn,hIn)   : gZH(idIn,hIn);
	dcomplex gOut = (doScale) ? fgZPH(idOut,hOut) : gZH(idOut,hOut);
	double widthterm = (doFixedWidth) ? w*mZP : s*(w/mZP);
	dcomplex A = gIn*gOut/(s-m2 + Im*widthterm);
	return A;
}
inline dcomplex hAZP(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	double w = wTotZP();
	dcomplex A = hADY(s,idIn,idOut,hIn,hOut) + hAZP0(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline dcomplex hAZPnoDY(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	double w = wTotZP();
	dcomplex A = hAZP0(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline double hA2ZP(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAZP(s,idIn,idOut,hIn,hOut);
			double angular = (1.+4.*hIn*hOut*cosTheta);
			double angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2ZP(double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAZP(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}
inline double hA2ZPnoDY(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAZPnoDY(s,idIn,idOut,hIn,hOut);
			double angular = (1.+4.*hIn*hOut*cosTheta);
			double angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2ZPnoDY(double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAZPnoDY(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}





// Minimal Z' amplitudes
inline dcomplex hAMinZP0(double s, double w, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	double m2 = mZP*mZP;
	dcomplex gIn  = gMinZPH(idIn,hIn);
	dcomplex gOut = gMinZPH(idOut,hOut);
	double widthterm = (doFixedWidth) ? w*mZP : s*(w/mZP);
	dcomplex A = gIn*gOut/(s-m2 + Im*widthterm);
	return A;
}
inline dcomplex hAMinZP(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	double w = wTotMinZP();
	dcomplex A = hADY(s,idIn,idOut,hIn,hOut) + hAMinZP0(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline dcomplex hAMinZPnoDY(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	double w = wTotMinZP();
	dcomplex A = hAMinZP0(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline double hA2MinZP(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAMinZP(s,idIn,idOut,hIn,hOut);
			double angular = (1.+4.*hIn*hOut*cosTheta);
			double angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2MinZP(double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAMinZP(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}
inline double hA2MinZPnoDY(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAMinZPnoDY(s,idIn,idOut,hIn,hOut);
			double angular = (1.+4.*hIn*hOut*cosTheta);
			double angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2MinZPnoDY(double s, unsigned int idIn, unsigned int idOut)
{
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			dcomplex A = hAMinZPnoDY(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}



inline void writeparameters(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	_ERR(1,"Event parameters: cosTheta="+str(cosTheta)+", s="+str(s)+", idIn="+str(idIn)+", idOut="+str(idOut));
}
inline void writeparameters(double s, unsigned int idIn, unsigned int idOut)
{
	_ERR(1,"Event parameters: s="+str(s)+", idIn="+str(idIn)+", idOut="+str(idOut));
}
inline void validateinput(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	if(s<0.)              { _ERR(1,"s<0., exitting now.");              exit(-1); }
	if(fabs(cosTheta)>1.) { _ERR(1,"fabs(cosTheta)>1., exitting now."); exit(-1); }
	if(idIn<=0)           { _ERR(1,"idIn<0, exitting now.");            exit(-1); }
	if(idOut<=0)          { _ERR(1,"idOut<0, exitting now.");           exit(-1); }
}
inline void validateinput(double s, unsigned int idIn, unsigned int idOut)
{
	if(s<0.)              { _ERR(1,"s<0., exitting now.");              exit(-1); }
	if(idIn<=0)           { _ERR(1,"idIn<0, exitting now.");            exit(-1); }
	if(idOut<=0)          { _ERR(1,"idOut<0, exitting now.");           exit(-1); }
}
inline void validateoutput(double N, double D)
{
	double R = N/D;
	if(std::isinf(R)) _FAT("value is infinity -> "+str(N)+"/"+str(D)+", exitting now.");
	if(std::isnan(R)) _FAT("value is not a number -> "+str(N)+"/"+str(D)+", exitting now.");
	if(R<=min_weight) _WRN(1,"value too small or negative -> "+str(N)+"/"+str(D));
	if(R>max_weight)  _WRN(1,"value is too large -> "+str(N)+"/"+str(D));
}
inline double weightZP(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(cosTheta,s,idIn,idOut);
	double N = hA2ZP(cosTheta,s,idIn,idOut);
	double D = hA2DY(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZP is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZP is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZP is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightZP(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2ZP(s,idIn,idOut);
	double D = hA2DY(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZP is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZP is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZP is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightZPnoDY(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(cosTheta,s,idIn,idOut);
	double N = hA2ZPnoDY(cosTheta,s,idIn,idOut);
	double D = hA2DY(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZPnoDY is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZPnoDY is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZPnoDY is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightZPnoDY(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2ZPnoDY(s,idIn,idOut);
	double D = hA2DY(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZPnoDY is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZPnoDY is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZPnoDY is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightMinZP(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(cosTheta,s,idIn,idOut);
	double N = hA2MinZP(cosTheta,s,idIn,idOut);
	double D = hA2DY(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2MinZP is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2MinZP is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2MinZP is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightMinZP(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2MinZP(s,idIn,idOut);
	double D = hA2DY(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2MinZP is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2MinZP is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2MinZP is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightMinZPnoDY(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(cosTheta,s,idIn,idOut);
	double N = hA2MinZPnoDY(cosTheta,s,idIn,idOut);
	double D = hA2DY(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2MinZPnoDY is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2MinZPnoDY is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2MinZPnoDY is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightMinZPnoDY(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2MinZPnoDY(s,idIn,idOut);
	double D = hA2DY(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2MinZPnoDY is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2MinZPnoDY is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2DY is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2DY is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2DY is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2MinZPnoDY is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}