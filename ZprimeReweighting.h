#include "std.h"

//// Some enums
enum pdg
{
	DWN=1, UP=2, STR=3, CHM=4, BOT=5, TOP=6,
	E=11, NUE=12, MU=13, NUMU=14, TAU=15, NUTAU=16,
	TAUPRIME=17, NUTAUPRIME=18,
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
static const double evv = 246; // GeV is the SM Higgs vacuum expectation value

static const double sw2     = 0.2312015; // ->PDG
static const double cw2     = 0.7687985; // ->PDG
static       double alphaEM = 0.008138; //<--at 3 TeV // 0.00729735; // ->PDG // IT MAY FLOAT !!!
static       double alphaST = 0.086; //<--at 3 TeV // 0.11856; // IT MAY FLOAT !!!
static const double Gmu     = 0.0000116633980690699; // GeV^-2

static const double GeV2mb  = 0.38937930419;      // 1/GeV^2 to mb (mili-barn) conversion constant
static const double GeV2nb  = 1e6*0.38937930419;  // 1/GeV^2 to nb (nano-barn) conversion constant
static const double GeV2pb  = 1e9*0.38937930419;  // 1/GeV^2 to pb (pico-barn) conversion constant
static const double GeV2fb  = 1e12*0.38937930419;  // 1/GeV^2 to fb (femto-barn) conversion constant

static bool doScale      = false;
static bool doScaleWidth = true; // turn on/off the g^2 scale in the BW width (default = true)
static double thetaP = 0.;
static double gammaP = 1.;



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

void setCouplingsScale(bool doscale) { doScale = doscale; }
void setScaleWidth(bool doscale)     { doScaleWidth = doscale; }

void setThetaPrime(double theta)
{
	if(theta<0 || theta>pi) _FAT("thetaP<0 || thetaP>pi: "<<theta);
	thetaP = theta;
}
void setGammaPrime(double gamma)
{
	if(gamma<0) _FAT("gammaP<0: "<<gamma);
	gammaP = gamma;
}
void setModelPrime(TString model)
{
	double theta = 0.;
	double gamma = 0.;
	if     (model=="psi") { theta = 0;                         gamma = 1.; }
	else if(model=="chi") { theta = -pi/2.;                    gamma = 1.; }
	else if(model=="eta") { theta = atan(-sqrt(5./3.))+pi/2.;  gamma = 1.; }
	else if(model=="I")   { theta = -asin(sqrt(5./8.));        gamma = 1.; }
	else _FAT("Model :"<<model<<", is not supported.");
	setThetaPrime(theta);
	setGammaPrime(gamma);
}
void setModelPrime(double theta, double gamma)
{
	setThetaPrime(theta);
	setGammaPrime(gamma);
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
// Z'/KK scale factor couplings:
// This should affect the type of the couplings (complex) and the widths which
// are proportional to (|fgX*gL|^2+|fgX*gR|^2) ~ |fgX|^2*(|gL|^2+|gR|^2)
static const dcomplex fgZPinit(1.,0.);
static       dcomplex fgZP(1.,0.);
inline void   setFgZP(double re, double im) { fgZP  = dcomplex(re,im); }
inline double fgZP2()                       { return real(fgZP*conj(fgZP)); }
inline void resetfgZP()                     { fgZP  = fgZPinit; }

//// SM
inline double gG(unsigned int id)  { return qf(id); }
inline double gZL(unsigned int id) { return (I3f(id)-qf(id)*sw2)/sqrt(sw2*cw2); }
inline double gZR(unsigned int id) { return -qf(id)*sw2/sqrt(sw2*cw2); }
inline double gZH(unsigned int id, double h)
{
	if     (h==-f12) return gZL(id);
	else if(h==+f12) return gZR(id);
	else _FAT("unknown helicity: "<<h);
	return 0.;
}

//// E6
inline double gVE6(unsigned int id)
{
	double gV = 0.;
	if     (id==s2f["nuel"]->id || id==s2f["numu"]->id || id==s2f["nutau"]->id) gV =1./6.*(sqrt(10.)*cos(thetaP)-3.*sqrt(6.)*sin(thetaP))*sqrt(sw2);
	else if(id==s2f["elec"]->id || id==s2f["muon"]->id || id==s2f["tau"]->id)   gV = -4./sqrt(6.)*sin(thetaP)*sqrt(sw2);
	else if(id==s2f["up"]->id   || id==s2f["chm"]->id  || id==s2f["top"]->id)   gV = 0.;
	else if(id==s2f["dwn"]->id  || id==s2f["str"]->id  || id==s2f["bot"]->id)   gV = +4./sqrt(6.)*sin(thetaP)*sqrt(sw2);
	else _FAT("id="+str(id)+" is not supported");
	return gV;
}
inline double gAE6(unsigned int id)
{
	double gA = 0.;
	if     (id==s2f["nuel"]->id || id==s2f["numu"]->id || id==s2f["nutau"]->id) gA = 1./6.*(sqrt(10.)*cos(thetaP)-3.*sqrt(6.)*sin(thetaP))*sqrt(sw2);
	else if(id==s2f["elec"]->id || id==s2f["muon"]->id || id==s2f["tau"]->id)   gA = 1./3.*(sqrt(10.)*cos(thetaP)-sqrt(6.)*sin(thetaP))*sqrt(sw2);
	else if(id==s2f["up"]->id   || id==s2f["chm"]->id  || id==s2f["top"]->id)   gA = 1./3.*(sqrt(10.)*cos(thetaP)+sqrt(6.)*sin(thetaP))*sqrt(sw2);
	else if(id==s2f["dwn"]->id  || id==s2f["str"]->id  || id==s2f["bot"]->id)   gA = 1./3.*(sqrt(10.)*cos(thetaP)-sqrt(6.)*sin(thetaP))*sqrt(sw2);
	else _FAT("id="+str(id)+" is not supported");
	return gA;
}
inline double gLE6(unsigned int id) { return (1./(2.*sqrt(cw2*sw2)))*(gVE6(id)+gAE6(id))/2.; }
inline double gRE6(unsigned int id) { return (1./(2.*sqrt(cw2*sw2)))*(gVE6(id)-gAE6(id))/2.; }
inline double gHE6(unsigned int id, double h)
{
	if     (h==-f12) return gLE6(id);
	else if(h==+f12) return gRE6(id);
	else _FAT("unknown helicity: "<<h);
	return 0.;
}

//// ZP (real methods for ZP are like for Z0)
inline dcomplex fgZPL(unsigned int id)           { return fgZP*gZL(id); }
inline dcomplex fgZPR(unsigned int id)           { return fgZP*gZR(id); }
inline dcomplex fgZPH(unsigned int id, double h) { return fgZP*gZH(id,h); }

//// ZP E6 
inline dcomplex fgE6L(unsigned int id)           { return fgZP*gLE6(id); }
inline dcomplex fgE6R(unsigned int id)           { return fgZP*gRE6(id); }
inline dcomplex fgE6H(unsigned int id, double h) { return fgZP*gHE6(id,h); }



///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// WIDTHS /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
static const double mZPinit = 2000.; // GeV
static double mZP = 2000.; // GeV


void setZPmass(double m) { mZP = m; }
void resetZPmass() { mZP = mZPinit; }

//// Z'-->ffbar
inline double wZP2ffbar(unsigned int id)
{
	double w = 0.;
	double mZP2 = mZP*mZP;
	double mf2 = mf(id)*mf(id);
	double gL = gZL(id);
	double gL2 = gL*gL;
	double gR = gZR(id);
	double gR2 = gR*gR;
	if(mZP<2.*mf(id))       return 0.;
	if((1.-4.*mf2/mZP2)<0.) return 0.;
	w = Ncf(id)*(alphaEM/6.)*mZP*sqrt(1.-4.*mf2/mZP2)*((gL2+gR2)+(mf2/mZP2)*(6.*gL*gR-gL2-gR2));
	if(id<=TOP) w *= (1.+alphaST/pi); // QCD corrections
	if(doScale && doScaleWidth) w *= fgZP2();
	// cout << "Gamma("<<namef(id)<<")=" << w << endl;
	return w;
}
//// Z'E6->ffbar
inline double wE62ffbar(unsigned int id)
{
	double w = 0.;
	double mZP2 = mZP*mZP;
	double mf2 = mf(id)*mf(id);
	double gL = gLE6(id);
	double gL2 = gL*gL;
	double gR = gRE6(id);
	double gR2 = gR*gR;
	if(mZP<2.*mf(id))       return 0.;
	if((1.-4.*mf2/mZP2)<0.) return 0.;
	w = Ncf(id)*(alphaEM/6.)*mZP*sqrt(1.-4.*mf2/mZP2)*((gL2+gR2)+(mf2/mZP2)*(6.*gL*gR-gL2-gR2));
	if(id<=TOP) w *= (1.+alphaST/pi); // QCD corrections
	if(doScale && doScaleWidth) w *= fgZP2();
	return w;
}
inline double wTotZP()
{
	double w = 0.;
	for(ui2fermion::iterator it=ui2f.begin() ; it!=ui2f.end() ; ++it) w += wZP2ffbar(it->first);
	// w = wZ0*mZP/mZ0 + wZP2ffbar(TOP);
	if(doScale && doScaleWidth) w *= fgZP2();
	return w;
}
inline double wTotE6()
{
	double w = 0.;
	for(ui2fermion::iterator it=ui2f.begin() ; it!=ui2f.end() ; ++it) w += wE62ffbar(it->first);
	if(doScale && doScaleWidth) w *= fgZP2();
	return w;
}



/////////////////////////////////////////////////////////////////
/////////////////////// AMPLITUDES //////////////////////////////
/////////////////////////////////////////////////////////////////
static const double       min_weight   = 1.e-30;
static const double       max_weight   = 1.e+10;
static bool               dokFactors   = false;
static bool               doFixedWidth = false;


void setkFactors(bool dokF) { dokFactors = dokF; }
void setFixedWidth(bool doFixed) { doFixedWidth = doFixed; }

inline dcomplex hAG0(double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	if(s<0.) return A;
	A = gG(idIn)*gG(idOut)/s;
	return A;
}
inline dcomplex hAZ0(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	if(s<0.) return A;
	double widthterm = (doFixedWidth) ? wZ0*mZ0 : s*(wZ0/mZ0);
	A = gZH(idIn,hIn)*gZH(idOut,hOut)/(s-mZ2 + Im*widthterm);
	return A;
}
inline dcomplex hASM(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	A = hAG0(s,idIn,idOut) + hAZ0(s,idIn,idOut,hIn,hOut);
	if(isnaninf(real(A*conj(A)))) _FAT("|hASM|^2 is nan/inf");
	return A;
}

inline dcomplex hAZP0(double s, double w, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	if(s<0.) return A;
	double mass = mZP;
	double m2 = mass*mass;
	dcomplex gIn  = (doScale) ? fgZPH(idIn,hIn)   : gZH(idIn,hIn);
	dcomplex gOut = (doScale) ? fgZPH(idOut,hOut) : gZH(idOut,hOut);
	double widthterm = (doFixedWidth) ? w*mass : s*(w/mass);
	A = gIn*gOut/(s-m2 + Im*widthterm);
	return A;
}
inline dcomplex hAZP(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	A = hASM(s,idIn,idOut,hIn,hOut); // the SM term
	double w = wTotZP();
	A += hAZP0(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline dcomplex hAZPnoSM(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	double w = wTotZP();
	A += hAZP0(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline dcomplex hAE60(double s, double w, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	if(s<0.) return A;
	double mass = mZP;
	double m2 = mass*mass;
	dcomplex gIn  = (doScale) ? fgE6H(idIn,hIn)   : gHE6(idIn,hIn);
	dcomplex gOut = (doScale) ? fgE6H(idOut,hOut) : gHE6(idOut,hOut);
	double widthterm = (doFixedWidth) ? w*mass : s*(w/mass);
	A = gIn*gOut/(s-m2 + Im*widthterm);
	return A;
}
inline dcomplex hAE6(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	A = hASM(s,idIn,idOut,hIn,hOut); // the SM term
	double w = wTotE6();
	A += hAE60(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline dcomplex hAE6noSM(double s, unsigned int idIn, unsigned int idOut, double hIn, double hOut)
{
	dcomplex A(0,0);
	double w = wTotE6();
	A += hAE60(s,w,idIn,idOut,hIn,hOut);
	return A;
}
inline double hA2SM(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	double angular = 0.;
	double angular2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hASM(s,idIn,idOut,hIn,hOut);
			angular = (1.+4.*hIn*hOut*cosTheta);
			angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2SM(double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hASM(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}
inline double hA2ZP(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	double angular = 0.;
	double angular2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAZP(s,idIn,idOut,hIn,hOut);
			angular = (1.+4.*hIn*hOut*cosTheta);
			angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2ZP(double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAZP(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}
inline double hA2ZPnoSM(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	double angular = 0.;
	double angular2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAZPnoSM(s,idIn,idOut,hIn,hOut);
			angular = (1.+4.*hIn*hOut*cosTheta);
			angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2ZPnoSM(double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAZPnoSM(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}
inline double hA2E6(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	double angular = 0.;
	double angular2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAE6(s,idIn,idOut,hIn,hOut);
			angular = (1.+4.*hIn*hOut*cosTheta);
			angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2E6(double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAE6(s,idIn,idOut,hIn,hOut);
			A2 += real(A*conj(A));
		}
	}
	return A2;
}
inline double hA2E6noSM(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	double angular = 0.;
	double angular2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAE6noSM(s,idIn,idOut,hIn,hOut);
			angular = (1.+4.*hIn*hOut*cosTheta);
			angular2 = angular*angular;
			A2 += real(A*conj(A))*angular2;
		}
	}
	return A2;
}
inline double hA2E6noSM(double s, unsigned int idIn, unsigned int idOut)
{
	dcomplex A(0,0);
	double A2 = 0.;
	for(double hIn=-f12 ; hIn<=+f12 ; hIn++)
	{
		for(double hOut=-f12 ; hOut<=+f12 ; hOut++)
		{
			A = hAE6noSM(s,idIn,idOut,hIn,hOut);
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
	double D = hA2SM(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZP is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZP is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZP is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightZP(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2ZP(s,idIn,idOut);
	double D = hA2SM(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZP is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZP is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZP is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightZPnoSM(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(cosTheta,s,idIn,idOut);
	double N = hA2ZPnoSM(cosTheta,s,idIn,idOut);
	double D = hA2SM(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZPnoSM is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZPnoSM is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZPnoSM is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightZPnoSM(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2ZPnoSM(s,idIn,idOut);
	double D = hA2SM(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2ZPnoSM is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2ZPnoSM is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2ZPnoSM is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightE6(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(cosTheta,s,idIn,idOut);
	double N = hA2E6(cosTheta,s,idIn,idOut);
	double D = hA2SM(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2E6 is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2E6 is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2E6 is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightE6(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2E6(s,idIn,idOut);
	double D = hA2SM(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2E6 is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2E6 is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2E6 is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightE6noSM(double cosTheta, double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(cosTheta,s,idIn,idOut);
	double N = hA2E6noSM(cosTheta,s,idIn,idOut);
	double D = hA2SM(cosTheta,s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2E6noSM is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2E6noSM is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2E6noSM is "+str(N)+", returning weight=0 for this event"); writeparameters(cosTheta,s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}
inline double weightE6noSM(double s, unsigned int idIn, unsigned int idOut)
{
	validateinput(s,idIn,idOut);
	double N = hA2E6noSM(s,idIn,idOut);
	double D = hA2SM(s,idIn,idOut);
	if(std::isinf(N)) {_ERR(1,"hA2E6noSM is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(N)) {_ERR(1,"hA2E6noSM is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isinf(D)) {_ERR(1,"hA2SM is inf, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(std::isnan(D)) {_ERR(1,"hA2SM is nan, returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(D<=0.)         {_ERR(1,"hA2SM is "+str(D)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	if(N<0.)          {_ERR(1,"hA2E6noSM is "+str(N)+", returning weight=0 for this event"); writeparameters(s,idIn,idOut); return 0.;}
	// validateoutput(N,D);
	return N/D;
}