// c++ classes
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <ctime>
// ROOT classes
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include <TMath.h>
#include "TVector3.h"
#include <TLegend.h>
#include "TCanvas.h"
#include <TPad.h>
#include "TStyle.h"
#include "TFile.h"
#include "TRandom.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"
#include "TStopwatch.h"
#include "TTimeStamp.h"

// constants
Double_t pi = TMath::Pi();
Double_t pi2 = pi*pi;
Double_t pi3 = pi2*pi;
Double_t twopi3 = 8.0*pi*pi*pi;
Double_t hbarc = 0.197327;
Double_t hbarc2 = hbarc*hbarc;
Double_t hbarc3 = hbarc2*hbarc;

// Not using 
const double rJpsi = 0.45;//fm
const double rJpsi2 = rJpsi*rJpsi;
Double_t mg = 0.0;


//==================== Parameters JPsi =======================//
Int_t isUps = 0.0;
Double_t mJpsi = 3.1;
Double_t mD = 1.869;
Double_t mQ = 1.87;
Double_t mQ2 = mQ*mQ;
Double_t ep0=0.64;
Double_t ep02=ep0*ep0;
Double_t epCon = 4.0*pi*pow(8./3.0,3)*pow(ep0, 3)/pow(mQ, 1.5);
Double_t redMass = mQ/2.0;
Double_t a0 = 1.0/TMath::Sqrt(mQ*ep0);
Double_t epConMod = (8.0*pi*16.0*16.0/27.0)*(a0/mQ);

//const Double_t nCC0_21 = (10.0*1658)/362.0; ;  
//const Double_t nJpsi0_21 = (0.18*1658)/362.0; ; 
//const Double_t nCC0_21_All = 22.0;  // used at 2.76 

Double_t Fac1 = 1.0;
const Double_t nCC0_21_All = Fac1*26.0; 
//const Double_t Fugicity = 0.45;   
const Double_t Fugicity = 0.9;   

const Double_t NCollMB = 392;
const Double_t NColl_0_5 = 1818;

const Double_t nCC0_21 = (nCC0_21_All*Fugicity*NColl_0_5)/NCollMB;  // Initial with shadowing 
const Double_t nJpsi0_21 = Fac1*(0.1223*NColl_0_5)/NCollMB;  // Initial pp to be used in RAA to be divided by 0.65 // now it is 0.932 (at three places) 

const Double_t JPsiIntSh_21 = 1.0/0.72;
const Double_t JPsiIntSh_1 = 1.0/0.72;

//const Double_t FormTau=0.89;
const Double_t FormTau=0.0;
const Double_t CNMSigmaMid=0.15;
const Double_t CNMSigmaFor=0.15; //fm2
//sigma in fm^2 (1 fm^2 = 10 mb)
Double_t SigmaPionDiss = 0.10;
int QQbarVar=1;

//AddMag = 0 --> Electric
//AddMag = 1 --> Electric + Magnetic
//Add Mag = -1 --> only Magnetic
 
Int_t AddMag = 1;


/*
//Parameters Y(1S)
Int_t isUps=1;
Double_t mJpsi = 9.46;
Double_t mD = 5.280;
Double_t mQ = 4.2;
Double_t mQ2 = mQ*mQ;
Double_t ep0=1.1;
Double_t ep02=ep0*ep0;
Double_t epCon = 4.0*pi*pow(8./3.0, 3)*pow(ep0, 3)/pow(mQ, 1.5);
Double_t redMass = mQ/2.0;

Double_t a0 = 1.0/TMath::Sqrt(mQ*ep0);
Double_t epConMod = (8.0*pi*16.0*16.0/27.0)*(a0/mQ);

const Double_t nCC0_21_All = 1.007;
const Double_t Fugicity = 1.0;   
const Double_t NCollMB = 392;
const Double_t NColl_0_5 = 1818;
const Double_t nCC0_21 = (nCC0_21_All*Fugicity*NColl_0_5)/NCollMB;
const Double_t nJpsi0_21 = (0.00201*NColl_0_5)/NCollMB;
const Double_t JPsiIntSh_21 = 1.0/0.8526;
const Double_t JPsiIntSh_1 = 1.0/0.8526; 
const Double_t FormTau=0.0;
const Double_t CNMSigmaMid=0.15;
const Double_t CNMSigmaFor=0.15;
Double_t SigmaPionDiss = 0.10;
int QQbarVar=4; 
Int_t AddMag = 1;
*/



/*
//Parameters Y(2S)
Int_t isUps=2;
Double_t mJpsi = 10.02;
Double_t mD = 5.280;
Double_t mQ = 4.2;
Double_t mQ2 = mQ*mQ;
Double_t ep0=1.10/4.0;
Double_t ep02=ep0*ep0;
Double_t epCon = 4.0*pi*pow(8./3.0, 3)*pow(ep0, 3)/pow(mQ, 1.5);
Double_t redMass = mQ/2.0;
Double_t a0 = 1.0/TMath::Sqrt(mQ*1.10);
Double_t epConMod = 16.0*(8.0*pi*16.0*16.0/27.0)*(a0/mQ);
const Double_t nCC0_21_All = 1.007;
const Double_t Fugicity = 1.0;   
const Double_t NCollMB = 392;
const Double_t NColl_0_5 = 1818;
const Double_t nCC0_21 = (nCC0_21_All*Fugicity*NColl_0_5)/NCollMB;
const Double_t nJpsi0_21 = (0.00201*NColl_0_5)/NCollMB;
const Double_t JPsiIntSh_21 = 1.0/0.8526;
const Double_t JPsiIntSh_1 = 1.0/0.8526; 
const Double_t FormTau=0.0;
const Double_t CNMSigmaMid=0.20;
const Double_t CNMSigmaFor=0.20;
Double_t SigmaPionDiss = 0.10;
int QQbarVar=5; 
Int_t AddMag = 1;
*/


/*
//Parameters Y(3S)
Int_t isUps=3;
Double_t mJpsi = 10.02;
Double_t mD = 5.280;
Double_t mQ = 4.2;
Double_t mQ2 = mQ*mQ;
Double_t ep0=1.10/4.0;
Double_t ep02=ep0*ep0;
Double_t epCon = 4.0*pi*pow(8./3.0, 3)*pow(ep0, 3)/pow(mQ, 1.5);
Double_t redMass = mQ/2.0;
const Double_t nCC0_21_All = 0.607;
const Double_t Fugicity = 0.45;   
const Double_t nCC0_21 = (nCC0_21_All*Fugicity*1658)/362.0;
const Double_t nJpsi0_21 = (0.00123*1658)/362.0; 
const Double_t JPsiIntSh_21 = 1.0/0.8880; // on 20.04.2015
const Double_t JPsiIntSh_1 = 1.0/0.856; 
Double_t a0 = 1.0/TMath::Sqrt(mQ*1.10);
Double_t epConMod = 16.0*(8.0*pi*16.0*16.0/27.0)*(a0/mQ);
const Double_t FormTau=0.0;
int QQbarVar=6; 
Int_t AddMag = 1;
*/



/*
//Parameters Chib(1P)
Int_t isUps=1;
Double_t mJpsi = 9.99;
Double_t mD = 5.280;
Double_t mQ = 4.2;
Double_t mQ2 = mQ*mQ;
Double_t ep0= 0.67;
Double_t ep02=ep0*ep0;
Double_t redMass = mQ/2.0;
const Double_t nCC0_21_All = 0.607;
const Double_t Fugicity = 0.45;   
const Double_t nCC0_21 = (nCC0_21_All*Fugicity*1658)/362.0;
const Double_t nJpsi0_21 = (0.00123*1658)/362.0; 
const Double_t JPsiIntSh_21 = 1.0/0.8880; // on 20.04.2015
const Double_t JPsiIntSh_1 = 1.0/0.856;
//Double_t gs2 = 16.0*pi*TMath::Sqrt(1.10/mQ)/3.0;
Double_t gs2 = 16.0*pi*TMath::Sqrt(4.0*0.67/mQ)/3.0;
Double_t epConMod = TMath::Power(2.0,7)*gs2/(9.0*mQ);
const Double_t FormTau=0.0;
int QQbarVar=7; 
Int_t AddMag =1;
*/




//Matter at extream condition
//T0 is 550 - 580 MeV and tau0 is 0.3
const Double_t mPi = 0.140; 
const Double_t RPb = 7.11;  //1.2*TMath::Power(208,1.0/3.0)

//const Double_t R05 = 0.92*RPb;
const Double_t R05 = 0.96*RPb;
const Double_t NPart05 = 384.0;
const Double_t NColl05 = 1819.0;
//===========================================//
//============= Variation ==================//
//==========================================//

//Double_t Fac = 0.5;
Double_t Fac = 1.0;
//Double_t Fac = 1.5;

//const Double_t tau0 = 0.1;  //         T0 = 0.744535 (GeV)
const Double_t tau0 = 0.3; //nominal T0 = 0.51603  (GeV)
//const Double_t tau0 = 0.6;   //nominal T0 = 0.409034 (GeV)

const Double_t dNChdEta = 1943; //0-5% at 5 TeV (ALICE)
//const Double_t SS=3.6*1.5*1600;
const Double_t am=5.0;
const Double_t SS=am*1.5*dNChdEta; //5 TeV

const Double_t Nf=2.5;
const Double_t aq = (7.0*Nf/60.0 + 16.0/90.0)*pi2;
const Double_t ah = 4.5*pi2/90.0;

Double_t aT = 0.1;
//Double_t aT = 0.1;
//Double_t z0=1.8*tau0; //0
//Double_t vZ=1.4;     //1.0
//for longitudnal
//Double_t aT = 0.0;   // 0
Double_t z0=0.0; //0
Double_t vZ=1.0;     //1.0

const Double_t VTau0 = (R05+0.5*aT*tau0*tau0)*(R05+0.5*aT*tau0*tau0)*(z0+vZ*tau0)*pi;

const Double_t ss05 = SS/VTau0;

const Double_t T0=TMath::Power(SS/(4.0*aq*VTau0),1.0/3.0)*hbarc;

//const Double_t TC = 0.170;
const Double_t TF = 0.140;

//const Double_t tauf = pow(T0/TC, 3.)*tau0;
const Double_t nPart0 = 384;  // NPart 0-5%
//const Double_t nColl0 = 1747; //NColl 0-5% at 2.76 TeV
const Double_t nColl0 = 1818; //NColl 0-5%

int NTau; double stepTau;
double Tau[10000], TempTau[10000], fQGP[10000];

///////////////////////////////
TH1D *HistRegenJpsiPt[1000];

// Dynamics 
Double_t Npart(int BinLow, int BinHigh);
Double_t NColl(int BinLow, int BinHigh);
Double_t NPartVsNColl(Double_t NPart);

Double_t Npart_276(int BinLow, int BinHigh);
Double_t NColl_276(int BinLow, int BinHigh);


//Double_t TAA(int BinLow, int BinHigh);

// Dissociation functions
Double_t IntDiss_All(Double_t PtMin);
Double_t IntDiss_All(Double_t PtMin, Double_t NPart);

Double_t IntDiss_NPartInt(Double_t CentMin, Double_t CentMax, Double_t PtRange, Double_t Pt);
Double_t IntDiss_PtInt(Double_t PtMin, Double_t NPart);
Double_t SigmaQuasiGluon(Int_t Flag, Double_t PJPsi, Double_t T);
//Double_t SigmaGluonDissS(Double_t s);
//Double_t SigmaGluonDissS_Med(Double_t s, Double_t T);
//Double_t EDiss(Double_t Flag, Double_t T );
//Double_t SigmaThermalAct(Double_t Flag, Double_t T);

//=== formation functions
Double_t SigmaFSMod(Double_t s);
Double_t FormFun(Double_t p1, Double_t p2, Double_t theta1, Double_t theta2, Double_t phi1, Double_t phi2);
Double_t FormRateMC_P_New(Double_t Temp, Int_t iTau);

Double_t FormRateMC_P(Double_t Temp, Double_t Pt); // formation function for plots

Double_t IntFormVsPt(Double_t PtMin, Double_t R0Cent, Double_t NPart);
Double_t IntFormVsPt(Double_t PtMin, Double_t R0Cent);
Double_t RegenratedQuarkoniaNPart(Double_t NCC, Double_t NJPsi, Double_t R0, Double_t IntSh, Int_t PtRange);
Double_t RegenratedQuarkoniaNPart(Double_t NCC, Double_t NJPsi, Double_t R0, Double_t IntSh, Int_t PtRange, Double_t NPart);

Double_t RegenratedQuarkoniaPt(Double_t CentMin, Double_t CentMax, Double_t Pt);


//============ jpsi pion =================//
Double_t fpion(Double_t pPi, Double_t T);
Double_t PionJPsiDiss(Double_t P1, Double_t T);
Double_t RhoPi(Double_t T);
Double_t PionDiss_All(Double_t PtMin);
Double_t PionDiss_PtInt(Double_t PtMin);
Double_t SigmaPionDissFunc(Double_t S);

Double_t CalculateTandf_LatticeEOS( Double_t ssCent, Double_t R0Cent);
Double_t SigmaQuasiGluon(Int_t Flag , Double_t P1, Double_t T);
Double_t SigmaGluonDissSMod(Double_t s);
Double_t fGluon(Double_t pg, Double_t T);
Double_t fcharm_thermal(Double_t p, Double_t T);
Double_t fcharm(Double_t p, Double_t T);
Double_t CNMVsNPart(Double_t NPart, TGraph *Shgrf);

//============ Quarkonia ELoss =================//

Double_t ELossVsPt(Double_t Pt, TF1 *PtFunc);

//======================= Data Functions =====================//
void Draw_ATLAS_JPsi_RaaVsPt();

void Draw_Clone_ATLAS_JPsi_RaaVsPt();

void Draw_ALICE_JPsi_RaaVsPt();
void Draw_ALICE_JPsi_RaaVsNpart();
void Draw_ATLAS_JPsi_RaaVsNpart();

void Draw_CMS_JPsi_RaaVsPt_0To100();


void Draw_CMS_JPsi_RaaVsPt_0To10();
void Draw_CMS_JPsi_RaaVsPt_10To30();
void Draw_CMS_JPsi_RaaVsPt_30To100();
void Draw_CMS_JPsi_RaaVsNpart();


void TotalCharmProductionCross_ALICE();
void TotalBeautyProductionCross_ALICE();


void Draw_CMS_Y1S_5TeV_RaaVsNpart();
void Draw_CMS_Y2S_5TeV_RaaVsNpart();

void Draw_CMS_Y1S_5TeV_RaaVsPt();
void Draw_CMS_Y2S_5TeV_RaaVsPt();

void Draw_ALICE_Y1S_5TeV_RaaVsNpart();


Double_t tsallis_fitting_function(Double_t* x, Double_t* par);
Double_t JPsi_fitting_function(Double_t* x, Double_t* par);


TGraphErrors *grf_GhostGraph(Int_t MarkerStyle, Int_t MarkerColor, Int_t LineColor);


//==================================== Lattice EOS ================================================//
TFile *fileEOS=new TFile("InRootFile/LatticeEOS_s95p-v1.2.root","R");

TGraph *grfSSVsTemp = (TGraph*)fileEOS->Get("grfSSVsTemp");
TGraph *TempVsFQGP = (TGraph*)fileEOS->Get("TempVsFQGP");
TGraph *TempVsFQGP2 = (TGraph*)fileEOS->Get("TempVsFQGP2");
//================================================================================================//
// pT upsilon from Pythia
//TFile *filejpsi=new TFile("InRootFile/JPsiPt.root","R");
TFile *filejpsi=new TFile("InRootFile/Acc_PP_JPsi.root","R");
TH1D *Jpsi_Pt = (TH1D*)filejpsi->Get("diMuonsPt_Gen");

TFile *fileUpsilon=new TFile("InRootFile/dimuonGenPt1s2sPbPb_Pt.root","R");
TH1D *Y1S_Pt = (TH1D*)fileUpsilon->Get("diMuonsPt_Gen1S");
TH1D *Y2S_Pt = (TH1D*)fileUpsilon->Get("diMuonsPt_Gen2S");


//================================== Shadowing from root file ==========================//
//TFile *fileShadowing = new TFile("InRootFile/shadowing_at_276tev_bMass946.root","R");
//TGraph *grf_shadowing_276tev_24_y_24 = (TGraph*)fileShadowing->Get("shadowing_276tev_24_y_24");

TFile *fileShadowing = new TFile("InRootFile/Shadowing_5_TeV/Bottom_502_MinBias.root","R");
TGraph *grf_shadowing_502tev_24_y_24_Y1S = (TGraph*)fileShadowing->Get("gr_shadowing");
TGraph *grf_shadowing_502tev_24_y_24_Y2S = (TGraph*)fileShadowing->Get("gr_shadowing");

TGraph *grf_shadowing_502tev_25_y_40_Y1S = (TGraph*)fileShadowing->Get("gr_shadowing_2540");
TGraph *grf_shadowing_502tev_25_y_40_Y2S = (TGraph*)fileShadowing->Get("gr_shadowing_2540");



//TFile *fileShadowing_pt = new TFile("InRootFile/17042015_Shadowing.root","R");            
TFile *fileShadowing_pt = new TFile("InRootFile/Shadowing_5_TeV/Charm_502_MinBias.root","R");            
//TH1D *HistJpsiRaaShVsPt_Y1 = (TH1D*)fileShadowing_pt->Get("HistJpsiRaaShVsPt_Y1");

TGraph *gr_shadowing_2540= (TGraph*)fileShadowing_pt->Get("gr_shadowing_2540");
//TH1D *HistJpsiRaaShVsPt_Y2440 = (TH1D*)fileShadowing_pt->Get("HistJpsiRaaShVsPt_Y2440");

TGraph *gr_shadowing_2424= (TGraph*)fileShadowing_pt->Get("gr_shadowing");
//TH1D *HistJpsiRaaShVsPt_Y21 = (TH1D*)fileShadowing_pt->Get("HistJpsiRaaShVsPt_Y21");


//TFile *fileShadowing_Y_pt = new TFile("InRootFile/Shadowing_5_TeV/Bottom_502_MinBias.root","R");      
//TH1D *HistUpsilonRaaShVsPt_Y21 = (TH1D*)fileShadowing_pt->Get("HistUpsilonRaaShVsPt_Y21");
//TGraph *gr_shadowing_Y_2424= (TGraph*)fileShadowing_Y_pt->Get("gr_shadowing");

/*
TFile *fileShadowingNPart = new TFile("InRootFile/2804_ShawQuarkonia_NPart.root","R");
TGraph *grShJPsiVsNPart_Y1= (TGraph*)fileShadowingNPart->Get("grShJPsiVsNPart_Y1");
TGraph *grShJPsiVsNPart_Y2440= (TGraph*)fileShadowingNPart->Get("grShJPsiVsNPart_Y2440");
TGraph *grShJPsiPtCutVsNPart_Y21= (TGraph*)fileShadowingNPart->Get("grShJPsiPtCutVsNPart_Y21");
TGraph *grShUpsilonVsNPartMid = (TGraph*)fileShadowingNPart->Get("grShUpsilonVsNPartMid");
TGraph *grShUpsilonVsNPartFor = (TGraph*)fileShadowingNPart->Get("grShUpsilonVsNPartFor");
*/

TFile *fileShadowingNPart_JPsi_2540 = new TFile("InRootFile/Shadowing_5_TeV/Charm_502_NPart_YPlus25_to_Plus40_Pt_00_50.root","R");
TGraph *grShJPsiVsNPart_Y2440= (TGraph*)fileShadowingNPart_JPsi_2540->Get("grShVsNpart_Full");

TFile *fileShadowingNPart_JPsi_2424 = new TFile("InRootFile/Shadowing_5_TeV/Charm_502_NPart_YMinus24_to_Plus24_Pt_65_50.root","R");
TGraph *grShJPsiPtCutVsNPart_Y24= (TGraph*)fileShadowingNPart_JPsi_2424->Get("grShVsNpart_Full");

TFile *fileShadowingNPart_JPsi_2020 = new TFile("InRootFile/Shadowing_5_TeV/Charm_502_NPart_YMinus20_to_Plus20_Pt_95_50.root","R");
TGraph *grShJPsiPtCutVsNPart_Y20= (TGraph*)fileShadowingNPart_JPsi_2020->Get("grShVsNpart_Full");
                                                                         
TFile *fileShadowingNPart_Y_2020 = new TFile("InRootFile/Shadowing_5_TeV/Bottom_502_NPart_YMinus24_to_Plus24_Pt_00_40.root","R");
TGraph *grShUpsilonVsNPartMid = (TGraph*)fileShadowingNPart_Y_2020->Get("grShVsNpart_Full");

TFile *fileShadowingNPart_Y_2540 = new TFile("InRootFile/Shadowing_5_TeV/Beauty_502_NPart_YTwoPointFour_to_Four.root","R");
TGraph *grShUpsilonVsNPartFor = (TGraph*)fileShadowingNPart_Y_2540->Get("grShVsNpart_Full");

//================================== Comover from root file ==========================//
TFile *fileJPsiComover = new TFile("InRootFile/OutJPsiPionCross.root","R");
TGraph *grfJPsiPionCross = (TGraph*)fileJPsiComover->Get("grCalcSigmaJPsiPion");


TFile *fileUpsilonComover = new TFile("InRootFile/OutUpsilonPionCross.root","R");
TGraph *grfUpsilonPionCross = (TGraph*)fileUpsilonComover->Get("grCalcSigmaJPsiPion");

TFile *fileUpsilon2SComover = new TFile("InRootFile/OutUpsilon2SPionCross.root","R");
TGraph *grfUpsilon2SPionCross = (TGraph*)fileUpsilon2SComover->Get("grCalcSigmaJPsiPion");


Double_t ErrorReCalculate(Double_t OldNumber, Double_t OldError, Double_t NewNumber);
Double_t ErrorReCalculateExp(Double_t NewNumber, Double_t ErrorPercent );

void QuarkSupp()
{

  cout<<endl<<endl;
  TTimeStamp ts1;
  cout << ts1.AsString() << endl; 
  cout<<endl<<endl;


  TStopwatch timer;


  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat("nmr");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.048,"xyz");
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadRightMargin(0.065);
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetTitleXOffset(1.15);
  //gStyle->SetTitleYOffset(1.2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetErrorX(0);   
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gROOT->ForceStyle();

  Int_t SaveTime =0;
  SaveTime =0;

  Int_t FORM =0;
  FORM =0;



  TLatex *tb = new TLatex();
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);



  //Double_t val = NPartVsNColl(397);
  //cout<<val<<endl;
  //return;

  

  // j/psi pT histogram  
  cout << "J/psi hist Integral = " << Jpsi_Pt->Integral() << endl;
  cout<<Jpsi_Pt->GetMean()<<"   "<<Jpsi_Pt->GetRMS()<<"  "<<Jpsi_Pt->GetSize()<<"  "<<Jpsi_Pt->FindBin(6.5)<<"  "<<Jpsi_Pt->GetBinWidth(1)<<"  "<<Jpsi_Pt->GetNbinsX()<<endl;
  
  Jpsi_Pt->Rebin(5); 
  
  cout << "J/psi hist Integral = " << Jpsi_Pt->Integral() << endl;
  cout<<Jpsi_Pt->GetMean()<<"   "<<Jpsi_Pt->GetRMS()<<"  "<<Jpsi_Pt->GetSize()<<"  "<<Jpsi_Pt->FindBin(6.5)<<"  "<<Jpsi_Pt->GetBinWidth(1)<<endl;
  
  for(int i=0;i<Jpsi_Pt->GetSize();i++)
    {
      cout<<" Jpsi Bin Centre "<<Jpsi_Pt->GetBinCenter(i)<<endl;
    }



  
  Y1S_Pt->Rebin(5); 
  cout << "Y1S hist Integral = " << Y1S_Pt->Integral() << endl;
  cout<<Y1S_Pt->GetMean()<<"   "<<Y1S_Pt->GetRMS()<<"  "<<Y1S_Pt->GetSize()<<"  "<<Y1S_Pt->FindBin(6.5)<<"  "<<Y1S_Pt->GetBinWidth(1)<<endl;
  

  Y2S_Pt->Rebin(5); 
  cout << "Y2S hist Integral = " << Y2S_Pt->Integral() << endl;
  cout<<Y2S_Pt->GetMean()<<"   "<<Y2S_Pt->GetRMS()<<"  "<<Y2S_Pt->GetSize()<<"  "<<Y2S_Pt->FindBin(6.5)<<"  "<<Y2S_Pt->GetBinWidth(1)<<endl;
 

  cout<<" T0 "<<T0<<"  aq "<<aq<<" ah "<<ah<<endl;
  cout<<" NPart_MB  "<<Npart(0,200)<<" NColl_MB  "<<NColl(0,200)<<endl;

 
  new TCanvas;
  TotalCharmProductionCross_ALICE();

  new TCanvas;
  TotalBeautyProductionCross_ALICE();



  /*
  Double_t ErrorReCalculate(Double_t OldNumber, Double_t OldError, Double_t NewNumber)
  cout<<" sigma pp ccbar  "<<ErrorReCalculate(4.11, 2.69, 6.754)<<"   "<<ErrorReCalculate(4.11, 2.50, 6.754)<<endl;
  cout<<"sigma PbPb ccbar "<<ErrorReCalculate(3.21, 2.1, 4.669)<<"   "<<ErrorReCalculate(3.21, 1.95, 4.669)<<endl;
  cout<<"N PbPb ccbar "<<ErrorReCalculate(18.21, 12.0, 26.23)<<"   "<<ErrorReCalculate(18.21, 11.0, 26.23)<<endl;
  cout<<endl<<endl;


  cout<<" sigma pp JPsi  "<<ErrorReCalculate(21.6, 10.6, 35.32)<<"   "<<ErrorReCalculate(21.6, 10.4, 35.32)<<endl;
  cout<<"sigma PbPb JPsi "<<ErrorReCalculate(16.83,8.26, 24.56)<<"   "<<ErrorReCalculate(16.83, 8.10, 24.56)<<endl;
  cout<<"N PbPb JPsi "<<ErrorReCalculate(0.0952, 0.047, 0.1381)<<"   "<<ErrorReCalculate(0.0952, 0.046, 0.1381)<<endl;
  cout<<endl<<endl;


  cout<<" sigma pp bb_bar  "<<ErrorReCalculate(110.5, 15.1, 210.3)<<"   "<<ErrorReCalculate(110.5, 14.2, 210.3)<<endl;
  cout<<"sigma PbPb bb_bar "<<ErrorReCalculate(100.5, 13.7, 179.3)<<"   "<<ErrorReCalculate(100.5, 12.9, 179.3)<<endl;
  cout<<"N PbPb bb_bar "<<ErrorReCalculate(0.57, 0.08, 1.007)<<"   "<<ErrorReCalculate(0.57, 0.07, 1.007)<<endl;
  cout<<endl<<endl;




  cout<<" sigma pp Y  "<<ErrorReCalculate(0.22, 0.07, 0.4206)<<"   "<<ErrorReCalculate(0.22, 0.06, 0.4206)<<endl;
  cout<<"sigma PbPb Y "<<ErrorReCalculate(0.199,0.063,0.3586)<<"   "<<ErrorReCalculate(0.199,0.054,0.3586)<<endl;
  cout<<"N PbPb Y "<<ErrorReCalculate(0.001123, 0.0004, 0.0020)<<"   "<<ErrorReCalculate(0.001123, 0.0003, 0.0020)<<endl;
  cout<<endl<<endl;
  */

  //Double_t ErrorReCalculateExp(Double_t NewNumber, Double_t ErrorPercent )

  Double_t CErrorPercent = (0.641509/6.75472);
  Double_t BErrorPercent = (46.217/210.302);

  cout<<"sigma PbPb ccbar "<<ErrorReCalculateExp(4.669,CErrorPercent)<<endl;
  cout<<"N PbPb ccbar     "<<ErrorReCalculateExp(26.23,CErrorPercent)<<endl;

  cout<<" sigma pp JPsi  "<<ErrorReCalculateExp(35.32,CErrorPercent)<<endl;
  cout<<"sigma PbPb JPsi "<<ErrorReCalculateExp(24.56,CErrorPercent)<<endl;
  cout<<"N PbPb JPsi     "<<ErrorReCalculateExp(0.1381,CErrorPercent)<<endl;



  cout<<"sigma PbPb bbbar "<<ErrorReCalculateExp(179.30,BErrorPercent)<<endl;
  cout<<"N PbPb bbbar     "<<ErrorReCalculateExp(1.007,BErrorPercent)<<endl;


  cout<<" sigma pp Y  "<<ErrorReCalculateExp(0.4206,BErrorPercent)<<endl;
  cout<<"sigma PbPb Y "<<ErrorReCalculateExp(0.3586,BErrorPercent)<<endl;
  cout<<"N PbPb Y     "<<ErrorReCalculateExp(0.0020,BErrorPercent)<<endl;


  return;


  new TCanvas;
  Draw_CMS_JPsi_RaaVsPt_0To100();
  Draw_Clone_ATLAS_JPsi_RaaVsPt();


  new TCanvas;
  Draw_CMS_JPsi_RaaVsPt_0To10();

  new TCanvas;
  Draw_CMS_JPsi_RaaVsPt_10To30();
  
  new TCanvas;
  Draw_CMS_JPsi_RaaVsPt_30To100();



  new TCanvas;
  Draw_CMS_JPsi_RaaVsNpart();

  
  //============= define the output rootfile ================//
  Char_t OutFileName[100];
  if(QQbarVar==1){sprintf(OutFileName,"JPsiCalculations.root");} 
  if(QQbarVar==4){sprintf(OutFileName,"Y1SCalculations.root");}
  if(QQbarVar==5){sprintf(OutFileName,"Y2SCalculations.root");}
  if(QQbarVar==6){sprintf(OutFileName,"Y3SCalculations.root");}
  if(QQbarVar==7){sprintf(OutFileName,"ChiBCalculations.root");}
  
  TFile *OutFile =new TFile(OutFileName,"RECREATE");




  //==================== Fit J/psi with Tsallis =================//

  TF1 *tsallis_fun = new TF1("tsallis_fun", tsallis_fitting_function, 0.0, 50.0, 5);
  tsallis_fun->SetParNames("dNdy", "n", "pzero","mass","Beta");
    
  tsallis_fun->SetParameters(11.0, 7.0, 2.0,3.1,0.45);
  
  //tsallis_fun->SetParLimits(0, 1000.0, 40000.0) ;
  tsallis_fun->SetParLimits(1, 1.0, 8.0) ;
  //tsallis_fun->SetParLimits(2, 0.5, 10.0); 
  
  //tsallis_fun->FixParameter(3, 3.1); 
  tsallis_fun->SetLineColor(2);
  tsallis_fun->SetLineStyle(1);
  tsallis_fun->SetLineWidth(2);
    
 
  new TCanvas;
  gPad->SetLogy(1);
  Jpsi_Pt->Fit(tsallis_fun, "ME", "", 1.0, 25.0);  
  Jpsi_Pt->SetMarkerStyle(20);
  Jpsi_Pt->Draw("P");
  tsallis_fun->Draw("same");
  
  
  TF1 *ptfit_fun = new TF1("ptfit_fun", JPsi_fitting_function, 0.0, 50.0, 4);
  ptfit_fun->SetParNames("#alpha_{1}", "#alpha_{2}", "#alpha_{3}","#alpha_{4}");
  ptfit_fun->SetLineColor(4);
  //ptfit_fun->SetParameters(3415,2.88,10.4,4.5);


  new TCanvas;
  gPad->SetLogy(1);
  Jpsi_Pt->Fit(ptfit_fun, "ME", "", 1.0, 25.0); 
  Jpsi_Pt->SetMarkerStyle(20);
  Jpsi_Pt->Draw("P");
  ptfit_fun->Draw("same");

  //return;

  //making all histograms for regenration 
  char NameHist[1000];
  for(int i=0;i<200;i++)
    {
      sprintf(NameHist,"HistRegenJpsiPt_%d",i);
      //cout<<NameHist<<endl;
      //HistRegenJpsiPt[i]=new TH1D(NameHist,NameHist,100,0.0,50.0);
      HistRegenJpsiPt[i]=new TH1D(NameHist,NameHist,50,0.0,25.0);
    }


  Double_t q0min=0.0;
  Double_t q0max = 3.0;
  if(QQbarVar ==4 )q0max = 6.0;

  Double_t q0step = 0.01;
  int NQ0 = (int)((q0max-q0min)/q0step); 
  Double_t Q0[2000], SigDQ0[2000], SigDQ0Mod[2000];
  
  for(int i=0; i<NQ0; i++) {

    Q0[i] = q0min + q0step*i;
    Double_t ss = 2.0*Q0[i]*mJpsi+mJpsi*mJpsi;

    SigDQ0[i] = SigmaGluonDissSMod(ss)*hbarc2*10.0;

    SigDQ0Mod[i] = SigmaGluonDissSMod(ss)*hbarc2*10.0;

    //Double_t SigmaGluonDissSMod(Double_t s)
    //Double_t q0 = (s-(mJpsi*mJpsi+mg*mg))/(2.0*mJpsi); 
    //SigDQ0[i] = SigmaD(1,Q0[i])*hbarc2*10.0;
    //SigDQ0Mod[i] = SigmaDMod(1,Q0[i])*hbarc2*10.0;
  }
  
  TGraph *grSigDQ0 = new TGraph(NQ0, Q0, SigDQ0);
  grSigDQ0->SetName("grSigDQ0");
  grSigDQ0->SetTitle("grSigDQ0");
  grSigDQ0->SetLineWidth(2);
  grSigDQ0->GetXaxis()->SetTitle("q^{0} (GeV)");
  grSigDQ0->GetYaxis()->SetTitle("#sigma (mb)");


  TGraph *grSigDQ0Mod = new TGraph(NQ0, Q0, SigDQ0Mod);
  grSigDQ0Mod->SetName("grSigDQ0Mod");
  grSigDQ0Mod->SetTitle("grSigDQ0Mod");
  grSigDQ0Mod->SetLineWidth(2);
  grSigDQ0Mod->SetLineColor(2);
  grSigDQ0Mod->GetXaxis()->SetTitle("q^{0} (GeV)");
  grSigDQ0Mod->GetYaxis()->SetTitle("#sigma (mb)");
  
  new TCanvas;
  gPad->SetTicks();
  grSigDQ0->Draw("AL");
  grSigDQ0Mod->Draw("Lsame");
  
  tb->DrawLatex(0.46,0.82,"g+J/#psi #rightarrow c+#bar{c}");
  
  gPad->SaveAs("Figures/Fig_SigmaDq0.pdf");
  gPad->SaveAs("Figures/Fig_SigmaDq0.png");
  

  grSigDQ0->Write();
  grSigDQ0Mod->Write();
  

  //return;

  // ================ dn/deta graph for making Temp as a function of nPart ==========================================//
  //===================== 0,2.5,5,7.5,10,20,30,40,50,60,70,80
  Double_t NPartdNdEta[20] = {398.0, 372.0, 346.0, 320.0, 263.0, 188.0, 131.0, 86.3, 53.6, 30.4, 15.6} ; 
  Double_t Err_NPartdNdEta[20] = {2.0, 3.0, 4.0, 4.0, 4.0, 3.0, 2.0, 1.7, 1.2, 0.8, 0.5} ; 
  
  Double_t dNdEtabyNpartby2[20] = {10.2, 9.9, 9.6, 9.4, 9.0, 8.4, 7.8, 7.4, 6.8, 6.3, 5.8} ; 
  Double_t Err_dNdEtabyNpartby2[20] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.5} ; 
  
  
  TGraphErrors *grdNDetaNpart = new TGraphErrors(11, NPartdNdEta, dNdEtabyNpartby2, Err_NPartdNdEta, Err_dNdEtabyNpartby2);
  grdNDetaNpart->SetLineWidth(2);
  grdNDetaNpart->SetMarkerStyle(20);

  grdNDetaNpart->GetYaxis()->SetRangeUser(3.5,12.0);
  grdNDetaNpart->GetYaxis()->SetTitle("#frac{dN}{d#eta}/(#frac{N_{Part}}{2})");
  grdNDetaNpart->GetXaxis()->SetTitle("N_{Part}");

  new TCanvas;
  gPad->SetLeftMargin(0.14);
  grdNDetaNpart->Draw("AP");
  
  gPad->SaveAs("Figures/Fig_DnDEtaVsNpartBy2.pdf");
  gPad->SaveAs("Figures/Fig_DnDEtaVsNpartBy2.png");

  
  Double_t TempNpart[400];
  Double_t TempNpart2[400];

  Double_t NpartVal[400];
  
  for(int i=0;i<=40;i++)
    {
      NpartVal[i]=i*10;
      
      TempNpart[i]=T0*TMath::Power(grdNDetaNpart->Eval(NpartVal[i])/grdNDetaNpart->Eval(Npart(0,10)),1.0/3); 
      //Double_t FacNPart = NpartVal[i]/Npart(0,10);
      //TempNpart[i]=T0*TMath::Power((grdNDetaNpart->Eval(NpartVal[i])*FacNPart)/grdNDetaNpart->Eval(Npart(0,10)),1.0/3); 

      TempNpart2[i]=T0*TMath::Power(NpartVal[i]/Npart(0,10),1.0/3); 
      
    }

  TGraph *grTempVsNpart = new TGraph(40, NpartVal,TempNpart);
  grTempVsNpart->SetLineWidth(2);
  grTempVsNpart->SetMarkerStyle(20);
  grTempVsNpart->SetMarkerColor(4);
  grTempVsNpart->GetXaxis()->SetTitle("N_{Part}");
  grTempVsNpart->GetYaxis()->SetTitle("T (GeV)");
  grTempVsNpart->GetYaxis()->SetRangeUser(0.1,1.0);
  
  TGraph *grTempVsNpart2 = new TGraph(40, NpartVal,TempNpart2);

  new TCanvas;
  gPad->SetTicks();
  grTempVsNpart->Draw("AP");
  grTempVsNpart2->Draw("Psame");
  gPad->SaveAs("Figures/Fig_TempVsNpart.pdf");
  gPad->SaveAs("Figures/Fig_TempVsNPart.png");
 
  //Double_t T0Cen = T0*TMath::Power((grdNDetaNpart->Eval(Npart(0,10)))/dNdEtabyNpartby2[0],1.0/3.0);
  //Double_t ssCen = ss05*(grdNDetaNpart->Eval(Npart(0,10))/dNdEtabyNpartby2[0]);  

  //return;
  
  double TauLatt[10000], TempTauLatt[10000];
  double fQGPLatt[10000];

  cout<<" Calculating Lattice equation of state "<<endl;
  NTau=0;
  CalculateTandf_LatticeEOS(ss05, R05);

  cout<<" Lattice EOS : "<<endl;
  
  cout<<" NTau "<< NTau<<" ss05 "<< ss05 <<endl;
    
  cout<<" Tau "<< "   "<< "Temp "<<endl;
  for(int i=0;i<NTau;i++){
    TauLatt[i]=Tau[i];
    TempTauLatt[i]=TempTau[i];
    fQGPLatt[i]=fQGP[i];
    cout<<i<<"  "<<Tau[i]<<"     "<<TempTau[i]<<endl;
  }

  

  TGraph *grTempVsTauLatt = new TGraph(NTau,TauLatt,TempTauLatt);
  grTempVsTauLatt->SetName("grTempVsTauLatt");
  grTempVsTauLatt->SetTitle("grTempVsTauLatt");
  grTempVsTauLatt->SetLineColor(1);
  grTempVsTauLatt->SetLineWidth(2);
  grTempVsTauLatt->SetLineStyle(9);
  grTempVsTauLatt->GetYaxis()->SetTitleOffset(1.4);
  grTempVsTauLatt->GetYaxis()->SetTitle("Temperature (GeV)");
  grTempVsTauLatt->GetXaxis()->SetTitle("#tau (fm)");

  new TCanvas;
  gPad->SetLeftMargin(0.14);
  gPad->SetTicks();
  grTempVsTauLatt->Draw("AL");


  TGraph *grFQGPVsTauLatt = new TGraph(NTau,TauLatt,fQGPLatt);
  grFQGPVsTauLatt->SetName("grFQGPVsTauLatt");
  grFQGPVsTauLatt->SetTitle("grFQGPVsTauLatt");
  grFQGPVsTauLatt->SetLineColor(1);
  grFQGPVsTauLatt->SetLineWidth(2);
  grFQGPVsTauLatt->SetLineStyle(9);
  grFQGPVsTauLatt->GetYaxis()->SetTitle("FQGP");
  grFQGPVsTauLatt->GetXaxis()->SetTitle("#tau (fm)");

  new TCanvas; 
  gPad->SetLeftMargin(0.14);
  gPad->SetTicks();
  grFQGPVsTauLatt->Draw("AL");



  Int_t NTau_020=0;
  double TauLatt_020[10000], TempTauLatt_020[10000];
  double fQGPLatt_020[10000];

  Double_t ss_020 = ss05*(grdNDetaNpart->Eval(Npart(40,100))/grdNDetaNpart->Eval(Npart(0,10)) );  
  Double_t R0_020 = R05*TMath::Power(Npart(40,100)/Npart(0,10),0.5);


  cout<<" Calculating Lattice equation of state for 0 -20 %"<<endl;
  NTau=0;
  CalculateTandf_LatticeEOS(ss_020, R0_020);
  NTau_020=NTau;
  
  cout<<" NTau "<<NTau_020<<" ss_020 "<< ss_020 <<endl;
    
  cout<<" Tau "<< "   "<< "Temp "<<endl;
  for(int i=0;i<NTau;i++){
    TauLatt_020[i]=Tau[i];
    TempTauLatt_020[i]=TempTau[i];
    fQGPLatt_020[i]=fQGP[i];
    cout<<i<<"  "<<Tau[i]<<"     "<<TempTau[i]<<endl;
  }

  

  TGraph *grTempVsTauLatt_020 = new TGraph(NTau_020,TauLatt_020,TempTauLatt_020);
  grTempVsTauLatt_020->SetName("grTempVsTauLatt_020");
  grTempVsTauLatt_020->SetTitle("grTempVsTauLatt_020");
  grTempVsTauLatt_020->SetLineColor(2);
  grTempVsTauLatt_020->SetLineWidth(2);
  grTempVsTauLatt_020->SetLineStyle(1);
  grTempVsTauLatt_020->GetYaxis()->SetTitleOffset(1.4);
  grTempVsTauLatt_020->GetYaxis()->SetTitle("Temperature (GeV)");
  grTempVsTauLatt_020->GetXaxis()->SetTitle("#tau (fm)");

  TLegend *legd_TempTau = new TLegend( 0.21,0.81,0.81,0.93);
  legd_TempTau->SetBorderSize(0);
  legd_TempTau->SetFillStyle(0);
  legd_TempTau->SetFillColor(0);
  legd_TempTau->SetTextSize(0.040);
  
  legd_TempTau->SetHeader("Lattice EOS, Cylindrical expansion");
  legd_TempTau->AddEntry(grTempVsTauLatt, "0-5%", "l");
  legd_TempTau->AddEntry(grTempVsTauLatt_020, "20-50%", "l");


  new TCanvas;
  gPad->SetLeftMargin(0.14);
  gPad->SetTicks();
  grTempVsTauLatt->Draw("AL");
  grTempVsTauLatt_020->Draw("Lsame");
  legd_TempTau->Draw("same");

  gPad->SaveAs("Figures/Fig_TauVsTemp.pdf");
  gPad->SaveAs("Figures/Fig_TauVsTemp.png");


  TGraph *grFQGPVsTauLatt_020 = new TGraph(NTau_020,TauLatt_020,fQGPLatt_020);
  grFQGPVsTauLatt_020->SetName("grFQGPVsTauLatt_020");
  grFQGPVsTauLatt_020->SetTitle("grFQGPVsTauLatt_020");
  grFQGPVsTauLatt_020->SetLineColor(2);
  grFQGPVsTauLatt_020->SetLineWidth(2);
  grFQGPVsTauLatt_020->SetLineStyle(1);
  grFQGPVsTauLatt_020->GetYaxis()->SetTitle("FQGP");
  grFQGPVsTauLatt_020->GetXaxis()->SetTitle("#tau (fm)");

  new TCanvas; 
  gPad->SetLeftMargin(0.14);
  gPad->SetTicks();
  grFQGPVsTauLatt->GetYaxis()->SetRangeUser(0.0,1.3);
  grFQGPVsTauLatt->Draw("AL"); 
  grFQGPVsTauLatt_020->Draw("Lsame");
 
  legd_TempTau->Draw("same");

  gPad->SaveAs("Figures/Fig_TauVsFQGP.pdf");
  gPad->SaveAs("Figures/Fig_TauVsFQGP.png");
  



  TGraph2D *g = new TGraph2D(NTau,TauLatt, TempTauLatt, fQGPLatt);
  g->SetTitle(";#tau(fm);T(GeV);f_{QGP}");
  g->SetMarkerColor(10);
  //g->GetXaxis()->SetTitle("#tau") ;
  g->GetXaxis()->CenterTitle() ;
  g->GetXaxis()->SetTitleOffset(1.4);
  //g->GetYaxis()->SetTitle("T(GeV)") ;
  g->GetYaxis()->CenterTitle() ;
  g->GetYaxis()->SetTitleOffset(1.8);
  g->GetYaxis()->SetNdivisions(505);

  //g->GetZaxis()->SetTitle("FQGP") ;
  g->GetZaxis()->CenterTitle() ;
  g->GetZaxis()->SetTitleOffset(1.2);
  g->GetZaxis()->SetNdivisions(505);

  TLegend *legd_TempTau2D = new TLegend( 0.25,0.76,0.67,0.85);
  legd_TempTau2D->SetBorderSize(0);
  legd_TempTau2D->SetFillStyle(0);
  legd_TempTau2D->SetFillColor(0);
  legd_TempTau2D->SetTextSize(0.040);
  
  //legd_TempTau2D->SetHeader("Lattice EOS, Cylindrical expansion");
  legd_TempTau2D->SetHeader("#splitline{PbPb #sqrt{s_{NN}} = 5.02 TeV}{0-5% Centrality}");
  //legd_TempTau2D->AddEntry(g, "0-5% Centrality", "");
  




  TCanvas *c1 = new TCanvas;
 
  //gPad->SetPalette(1);
  gPad->SetTicks(1);
  g->Draw("p0") ; 
  legd_TempTau2D->Draw("same");
  g->Draw("surf2 same") ; 
  gPad->SaveAs("Figures/Fig_TauVsFQGP2D_1.png");
  
  new TCanvas;
  g->Draw("p0") ; 
  legd_TempTau2D->Draw("same");
  

  
  TGraph2D *g2 = new TGraph2D(NTau_020, TauLatt_020, TempTauLatt_020, fQGPLatt_020);
  g2->SetTitle(";#tau(fm);T(GeV);f_{QGP}");
  g2->SetMarkerColor(10);
  g2->GetXaxis()->CenterTitle() ;
  g2->GetXaxis()->SetTitleOffset(1.4);
  g2->GetYaxis()->CenterTitle() ;
  g2->GetYaxis()->SetTitleOffset(1.8);
  g2->GetYaxis()->SetNdivisions(505);
  
  g2->GetZaxis()->CenterTitle() ;
  g2->GetZaxis()->SetTitleOffset(1.2);
  g2->GetZaxis()->SetNdivisions(505);

  TLegend *legd_TempTau2D_2 = new TLegend( 0.25,0.76,0.67,0.85);
  legd_TempTau2D_2->SetBorderSize(0);
  legd_TempTau2D_2->SetFillStyle(0);
  legd_TempTau2D_2->SetFillColor(0);
  legd_TempTau2D_2->SetTextSize(0.040);
  legd_TempTau2D_2->SetHeader("#splitline{PbPb #sqrt{s_{NN}} = 5.02 TeV}{20-50% Centrality}");
  

  new TCanvas;
  gPad->SetTicks(1);
  g2->Draw("p0") ; 
  legd_TempTau2D_2->Draw("same");
  g2->Draw("surf2 same") ; 
  gPad->SaveAs("Figures/Fig_TauVsFQGP2D_2.png");
  
  new TCanvas;
  g2->Draw("p0") ; 
  legd_TempTau2D_2->Draw("same");
  



  return;


  //============================================================================================================================//
  //==================================== DISSOCIATION RATES ====================================================================//
  //===========================================================================================================================//
  
  Double_t LambdaD_Khar[1000]={0};
  Double_t LambdaD_Khar36[1000]={0};  
  Double_t LambdaD_Khar86[1000]={0};
  
  Double_t FCharm_Temp[1000]={0.0};
  Double_t FCharm_Thermal_Temp[1000]={0.0};


  Double_t TempD[1000];  
  Double_t TempDMin=0.17;
  Double_t TempDMax=1.0;
  Double_t TempDStep=0.05;
  int NTempD= (int)((TempDMax-TempDMin)/TempDStep);


  for(int i = 0; i<=NTempD; i++) {
    TempD[i]=TempDMin+TempDStep*i;
    
    LambdaD_Khar[i]= SigmaQuasiGluon(1,0.0,TempD[i])/hbarc;
    LambdaD_Khar36[i]= SigmaQuasiGluon(1,3.6,TempD[i])/hbarc;
    LambdaD_Khar86[i]= SigmaQuasiGluon(1,8.6,TempD[i])/hbarc;

    FCharm_Temp[i]=fcharm(3.6,TempD[i]);
    FCharm_Thermal_Temp[i]=fcharm_thermal(3.6,TempD[i]);


  }  


  TGraph *DissRate_Khar = new TGraph(NTempD,TempD,LambdaD_Khar);
  DissRate_Khar->SetName("DissRateVsT_KharPt0");
  DissRate_Khar->SetTitle("DissRateVsT_KharPt0");
  DissRate_Khar->SetLineWidth(2);
  DissRate_Khar->SetLineColor(1);
  DissRate_Khar->SetLineStyle(1);
  DissRate_Khar->GetXaxis()->SetTitle("Temperature (GeV)");
  DissRate_Khar->GetYaxis()->SetTitle("Dissociation Rate #lambda_{D}#rho_{g} (fm^{-1})");



  TGraph *DissRate_Khar36 = new TGraph(NTempD,TempD,LambdaD_Khar36);
  DissRate_Khar36->SetName("DissRateVsT_KharPt36");
  DissRate_Khar36->SetTitle("DissRateVsT_KharPt36");
  DissRate_Khar36->SetLineWidth(2);
  DissRate_Khar36->SetLineColor(4);
  DissRate_Khar36->SetLineStyle(4);
  
  TGraph *DissRate_Khar86 = new TGraph(NTempD,TempD,LambdaD_Khar86);
  DissRate_Khar86->SetName("DissRateVsT_KharPt86");
  DissRate_Khar86->SetTitle("DissRateVsT_KharPt86");
  DissRate_Khar86->SetLineWidth(2);
  DissRate_Khar86->SetLineColor(6);
  DissRate_Khar86->SetLineStyle(6);

  // Diss rates
  TLegend *legd_dissT = new TLegend( 0.24,0.67,0.69,0.86);
  legd_dissT->SetBorderSize(0);
  legd_dissT->SetFillStyle(0);
  legd_dissT->SetFillColor(0);
  legd_dissT->SetTextSize(0.040);
  
  legd_dissT->AddEntry(DissRate_Khar, "p_{T}= 0", "l");
  legd_dissT->AddEntry(DissRate_Khar36, "p_{T}= 3.6 GeV/c", "l");
  legd_dissT->AddEntry(DissRate_Khar86, "p_{T}= 8.6 GeV/c", "l");
  
  new TCanvas;
  gPad->SetTicks();
  DissRate_Khar->Draw("AL");
  DissRate_Khar36->Draw("Lsame");
  DissRate_Khar86->Draw("Lsame");
  legd_dissT->Draw("same");
  gPad->SaveAs("Figures/Fig_DRateVsT.pdf");  
  gPad->SaveAs("Figures/Fig_DRateVsT.png");  


  TGraph *grf_FCharm_Temp = new TGraph(NTempD,TempD,FCharm_Temp);
  grf_FCharm_Temp->SetName("grf_FCharm_Temp");
  grf_FCharm_Temp->SetTitle("grf_FCharm_Temp");
  grf_FCharm_Temp->SetLineWidth(2);
  grf_FCharm_Temp->SetLineColor(1);
  grf_FCharm_Temp->SetLineStyle(1);
  grf_FCharm_Temp->GetXaxis()->SetTitle("Temperature (GeV)");
  grf_FCharm_Temp->GetYaxis()->SetTitle("fcharm");

  TGraph *grf_FCharm_Thermal_Temp = new TGraph(NTempD,TempD,FCharm_Thermal_Temp);
  grf_FCharm_Thermal_Temp->SetName("grf_FCharm_Thermal_Temp");
  grf_FCharm_Thermal_Temp->SetTitle("grf_FCharm_Thermal_Temp");
  grf_FCharm_Thermal_Temp->SetLineWidth(2);
  grf_FCharm_Thermal_Temp->SetLineColor(2);
  grf_FCharm_Thermal_Temp->SetLineStyle(1);

  TLegend *legd_fcharm = new TLegend( 0.24,0.67,0.69,0.86);
  legd_fcharm->SetBorderSize(0);
  legd_fcharm->SetFillStyle(0);
  legd_fcharm->SetFillColor(0);
  legd_fcharm->SetTextSize(0.040);
  
  legd_fcharm->AddEntry(grf_FCharm_Thermal_Temp, "thermal", "l");
  legd_fcharm->AddEntry(grf_FCharm_Temp, "tsallis", "l");


  new TCanvas;
  gPad->SetTicks();
  grf_FCharm_Temp->Draw("AL");
  grf_FCharm_Thermal_Temp->Draw("Lsame");
  legd_fcharm->Draw("Lsame");

  
 //========================= Diss Rate Vs pT ==========================//
  Double_t DissRateKharVsPt1[10000]={0};
  Double_t DissRateKharVsPt2[10000]={0};
  Double_t DissRateKharVsPt3[10000]={0};
  

  Double_t FCharm_Pt[1000]={0.0};
  Double_t FCharm_Thermal_Pt[1000]={0.0};

  Double_t PtD[10000];
  Double_t PtDMax=50.0;
  Double_t PtDMin=0.25;
  Double_t PtDStep=0.1;
    
  int NPtD=int((PtDMax-PtDMin)/PtDStep);
  
  for(int i = 0; i<=NPtD; i++) {
    
    PtD[i]=PtDMin+PtDStep*i;
    DissRateKharVsPt1[i]= SigmaQuasiGluon(1,PtD[i],0.200)/hbarc;
    DissRateKharVsPt2[i]= SigmaQuasiGluon(1,PtD[i],0.400)/hbarc;
    DissRateKharVsPt3[i]= SigmaQuasiGluon(1,PtD[i],0.600)/hbarc;


    FCharm_Pt[i]=fcharm(PtD[i],0.400);
    FCharm_Thermal_Pt[i]=fcharm_thermal(PtD[i],0.400);


  }
   

  //========================= Diss Rate Vs pT graphs ==========================//

  //DRate vs Pt graph
  TGraph *grDissRateVsPt1 = new TGraph(NPtD,PtD,DissRateKharVsPt1);
  grDissRateVsPt1->SetName("grDissRateVsPt_T200");
  grDissRateVsPt1->SetTitle("grDissRateVsPt_T200");
  grDissRateVsPt1->SetLineWidth(2);
  grDissRateVsPt1->SetLineColor(2);
  grDissRateVsPt1->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grDissRateVsPt1->GetYaxis()->SetTitle("Dissociation Rate #lambda_{D}#rho_{g} (fm^{-1})");
  
  
  TGraph *grDissRateVsPt2 = new TGraph(NPtD,PtD,DissRateKharVsPt2);
  grDissRateVsPt2->SetName("grDissRateVsPt_T400");
  grDissRateVsPt2->SetTitle("grDissRateVsPt_T400");
  grDissRateVsPt2->SetLineWidth(2);
  grDissRateVsPt2->SetLineStyle(6);  
  grDissRateVsPt2->SetLineColor(1);
  
  
  TGraph *grDissRateVsPt3 = new TGraph(NPtD,PtD,DissRateKharVsPt3);
  grDissRateVsPt3->SetName("grDissRateVsPt_T600");
  grDissRateVsPt3->SetTitle("grDissRateVsPt_T600");
  grDissRateVsPt3->SetLineWidth(2);
  grDissRateVsPt3->SetLineStyle(4);
  grDissRateVsPt3->SetLineColor(4);
  
  

  TLegend *legd_dissPt = new TLegend( 0.31,0.67,0.76,0.86);
  legd_dissPt->SetBorderSize(0);
  legd_dissPt->SetFillStyle(0);
  legd_dissPt->SetFillColor(0);
  legd_dissPt->SetTextSize(0.040);
  legd_dissPt->AddEntry(grDissRateVsPt1, "T=0.2 GeV", "l");
  legd_dissPt->AddEntry(grDissRateVsPt2, "T=0.4 GeV", "l");
  legd_dissPt->AddEntry(grDissRateVsPt3, "T=0.6 GeV", "l");


  new TCanvas; 
  gPad->SetTicks();
  grDissRateVsPt1->GetYaxis()->SetRangeUser(0.0,6.0);
  grDissRateVsPt1->Draw("AL");
  grDissRateVsPt2->Draw("Lsame");
  grDissRateVsPt3->Draw("Lsame");
  legd_dissPt->Draw("same");
  gPad->SaveAs("Figures/Fig_DRateVsPt.pdf");  
  gPad->SaveAs("Figures/Fig_DRateVsPt.png");  

 
  TGraph *grf_FCharm_Pt = new TGraph(NPtD,PtD,FCharm_Pt);
  grf_FCharm_Pt->SetName("grf_FCharm_Pt");
  grf_FCharm_Pt->SetTitle("grf_FCharm_Pt");
  grf_FCharm_Pt->SetLineWidth(2);
  grf_FCharm_Pt->SetLineColor(1);
  grf_FCharm_Pt->SetLineStyle(1);
  grf_FCharm_Pt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grf_FCharm_Pt->GetYaxis()->SetTitle("fcharm");

  TGraph *grf_FCharm_Thermal_Pt = new TGraph(NPtD,PtD,FCharm_Thermal_Pt);
  grf_FCharm_Thermal_Pt->SetName("grf_FCharm_Thermal_Pt");
  grf_FCharm_Thermal_Pt->SetTitle("grf_FCharm_Thermal_Pt");
  grf_FCharm_Thermal_Pt->SetLineWidth(2);
  grf_FCharm_Thermal_Pt->SetLineColor(2);
  grf_FCharm_Thermal_Pt->SetLineStyle(1);

  TLegend *legd_fcharm_pt = new TLegend( 0.24,0.67,0.69,0.86);
  legd_fcharm_pt->SetBorderSize(0);
  legd_fcharm_pt->SetFillStyle(0);
  legd_fcharm_pt->SetFillColor(0);
  legd_fcharm_pt->SetTextSize(0.040);
  
  legd_fcharm_pt->AddEntry(grf_FCharm_Thermal_Pt, "thermal", "l");
  legd_fcharm_pt->AddEntry(grf_FCharm_Pt, "tsallis", "l");


  new TCanvas;
  gPad->SetTicks();
  grf_FCharm_Pt->Draw("AL");
  grf_FCharm_Thermal_Pt->Draw("Lsame");
  legd_fcharm_pt->Draw("Lsame");


  //========================================================================================================//
  //==================================== Formation rate as a function of Temperature =======================//
  //========================================================================================================//

  cout<<" ==================== Formation Rate vs Temp ===================== " <<endl;
  cout<<"Temp: "<<" MC Int  "<< " Triple Int " << " (MC Int - Triple Int / MC Int)%"<<endl;
  
  Double_t LambdaFMC_P1[1000]={0};
  Double_t LambdaFMC_P2[1000]={0};
  Double_t LambdaFMC_P3[1000]={0};
      
  Double_t TempF[1000];  
  Double_t TempFMin=0.17;
  Double_t TempFMax=0.9;
  Double_t TempFStep=0.05;
  
  int NTempF= (int)((TempFMax-TempFMin)/TempFStep);
  
  
  for(int i = 0; i<=NTempF; i++) {
    TempF[i]=TempFMin+TempFStep*i;
    
    if(FORM ==1){
      
      LambdaFMC_P1[i]=FormRateMC_P(TempF[i],3.5);
      LambdaFMC_P2[i]=FormRateMC_P(TempF[i],4.5);
      LambdaFMC_P3[i]=FormRateMC_P(TempF[i],6.5);

      //LambdaFMC_P1[i]=FormRateMC_P(TempF[i],0.0);
      //LambdaFMC_P2[i]=FormRateMC_P(TempF[i],3.6);
      //LambdaFMC_P3[i]=FormRateMC_P(TempF[i],8.6);
      
    }
  }
      



  TGraph *grFormRateVsT_MCInt_P1 = new TGraph(NTempF,TempF,LambdaFMC_P1);
  grFormRateVsT_MCInt_P1->SetName("grFormRateVsT_MCInt_P1");
  grFormRateVsT_MCInt_P1->SetTitle("grFormRateVsT_MCInt_P1");
  grFormRateVsT_MCInt_P1->SetLineWidth(2);
  grFormRateVsT_MCInt_P1->SetLineColor(2);
  grFormRateVsT_MCInt_P1->SetLineStyle(1);
  
  TGraph *grFormRateVsT_MCInt_P2 = new TGraph(NTempF,TempF,LambdaFMC_P2);
  grFormRateVsT_MCInt_P2->SetName("grFormRateVsT_MCInt_P2");
  grFormRateVsT_MCInt_P2->SetTitle("grFormRateVsT_MCInt_P2");
  grFormRateVsT_MCInt_P2->SetLineWidth(2);
  grFormRateVsT_MCInt_P2->SetLineColor(1);
  grFormRateVsT_MCInt_P2->SetLineStyle(2);
  
  
  TGraph *grFormRateVsT_MCInt_P3 = new TGraph(NTempF,TempF,LambdaFMC_P3);
  grFormRateVsT_MCInt_P3->SetName("grFormRateVsT_MCInt_P3");
  grFormRateVsT_MCInt_P3->SetTitle("grFormRateVsT_MCInt_P3");
  grFormRateVsT_MCInt_P3->SetLineWidth(2);
  grFormRateVsT_MCInt_P3->SetLineColor(4);
  grFormRateVsT_MCInt_P3->SetLineStyle(4);
  
  TLegend *lgd_FRVsT = new TLegend( 0.51,0.78,0.84,0.92);
  lgd_FRVsT->SetBorderSize(0);
  lgd_FRVsT->SetFillStyle(0);
  lgd_FRVsT->SetFillColor(0);
  lgd_FRVsT->SetTextSize(0.040);
  lgd_FRVsT->AddEntry(grFormRateVsT_MCInt_P1, "p_{T} = 3.5 GeV/c", "L");
  lgd_FRVsT->AddEntry(grFormRateVsT_MCInt_P2, "p_{T} = 4.5 GeV/c", "L");
  lgd_FRVsT->AddEntry(grFormRateVsT_MCInt_P3, "p_{T} = 6.5 GeV/c", "L");
  
  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  grFormRateVsT_MCInt_P1->GetXaxis()->SetTitle("Temperature (GeV)");
  grFormRateVsT_MCInt_P1->GetYaxis()->SetTitleOffset(1.1);
  grFormRateVsT_MCInt_P1->GetYaxis()->SetTitle("d#lambda_{F}/dp_{T}(fm^{2}GeV^{-1})");
  grFormRateVsT_MCInt_P1->GetYaxis()->SetRangeUser(0.00001,1.0);
  grFormRateVsT_MCInt_P1->Draw("AL");
  grFormRateVsT_MCInt_P2->Draw("Lsame");
  grFormRateVsT_MCInt_P3->Draw("Lsame");
  
  
  lgd_FRVsT->Draw("same");
  //tb->DrawLatex(0.16,0.90,"(a)");
  gPad->SaveAs("Figures/Fig4a_FRateVsT.png");
  gPad->SaveAs("Figures/Fig4a_FRateVsT.pdf");
  gPad->SaveAs("Figures/Fig4a_FRateVsT.eps");
  
  
  




  //_______________________________________________________________________________________________________//
  //_______________________________ Formation rate as a function of pT ___________________________________//
  //_____________________________________________________________________________________________________//

  cout<<" _____________________ Formation Rate as a function of pT  ____________________________________" <<endl;
  

    
  cout<<" Formation Rate vs P_{J/psi}(GeV/c) :"<<endl;
  

  const int NPtF = 26;
  

  //Double_t PtF[NPtF]={0.2, 0.6, 1, 1.4, 1.8, 2.2, 2.6, 3, 3.4, 3.8, 4.2, 4.6, 5, 5.4, 5.8, 6.2, 6.6,  7, 7.4, 7.8,
  //		      8.2, 8.6, 9, 9.4, 9.8, 10.2, 10.6, 11, 11.4, 11.8};
  
  Double_t PtF[NPtF]={0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75};
    

  Double_t LambdaF_P1[NPtF]={0};
  Double_t LambdaF_P2[NPtF]={0};
  Double_t LambdaF_P3[NPtF]={0};
  


  //Double_t PtFMax=8.75;
  //Double_t PtFMin=0.25;
  //Double_t PtFStep=0.5;
  //int NPtF=int((PtFMax-PtFMin)/PtFStep);
  

  for(int i = 0; i<=NPtF; i++) {
    
    //PtF[i]=PtFMin+PtFStep*i;
    
    if(FORM ==1){
      LambdaF_P1[i]=FormRateMC_P(0.200,PtF[i]);
      LambdaF_P2[i]=FormRateMC_P(0.400,PtF[i]);
      LambdaF_P3[i]=FormRateMC_P(0.600,PtF[i]);
    }
    
    cout<<PtF[i]<<"  "<< LambdaF_P1[i] <<"  "<<LambdaF_P2[i]<<"   "<<LambdaF_P3[i]<<endl;
  }
  
  
  
  TGraph *grFormRate_TripleInt_P1 = new TGraph(NPtF,PtF,LambdaF_P1);
  grFormRate_TripleInt_P1->SetName("grFormRateVsP_MCInt_T1");
  grFormRate_TripleInt_P1->SetTitle("grFormRateVsP_MCInt_T1");
  grFormRate_TripleInt_P1->SetMarkerStyle(20);
  grFormRate_TripleInt_P1->SetLineWidth(2);
  grFormRate_TripleInt_P1->SetLineColor(2);
  grFormRate_TripleInt_P1->SetLineStyle(1);
  grFormRate_TripleInt_P1->GetXaxis()->SetTitle("P(GeV)");
  grFormRate_TripleInt_P1->GetYaxis()->SetTitle("FormRate(Triple Int)");
  
  TGraph *grFormRate_TripleInt_P2 = new TGraph(NPtF,PtF,LambdaF_P2);
  grFormRate_TripleInt_P2->SetName("grFormRateVsP_MCInt_T2");
  grFormRate_TripleInt_P2->SetTitle("grFormRateVsP_MCInt_T2");
  grFormRate_TripleInt_P2->SetLineWidth(2);
  grFormRate_TripleInt_P2->SetLineColor(1);
  grFormRate_TripleInt_P2->SetLineStyle(2);
  
  TGraph *grFormRate_TripleInt_P3 = new TGraph(NPtF,PtF,LambdaF_P3);
  grFormRate_TripleInt_P3->SetName("grFormRateVsP_MCInt_T3");
  grFormRate_TripleInt_P3->SetTitle("grFormRateVsP_MCInt_T3");
  grFormRate_TripleInt_P3->SetLineWidth(2);
  grFormRate_TripleInt_P3->SetLineColor(4);
  grFormRate_TripleInt_P3->SetLineStyle(4);
  
  TLegend *legd_formPt = new TLegend( 0.44,0.73,0.90,0.92);
  legd_formPt->SetBorderSize(0);
  legd_formPt->SetFillStyle(0);
  legd_formPt->SetFillColor(0);
  legd_formPt->SetTextSize(0.040);  
  
  legd_formPt->AddEntry(grFormRate_TripleInt_P1,"T=0.2 GeV", "L");
  legd_formPt->AddEntry(grFormRate_TripleInt_P2,"T=0.4 GeV", "L");
  legd_formPt->AddEntry(grFormRate_TripleInt_P3,"T=0.6 GeV", "L");
  
  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.14);
  grFormRate_TripleInt_P1->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grFormRate_TripleInt_P1->GetYaxis()->SetTitleOffset(1.2);
  grFormRate_TripleInt_P1->GetYaxis()->SetTitle("d#lambda_{F}/dp_{T} (fm^{2}GeV^{-1})");
  grFormRate_TripleInt_P1->GetYaxis()->SetRangeUser(0.00001,1.0);
  
  grFormRate_TripleInt_P1->Draw("AC");
  grFormRate_TripleInt_P2->Draw("Csame");
  grFormRate_TripleInt_P3->Draw("Csame");
  legd_formPt->Draw("same");
  //tb->DrawLatex(0.16,0.90,"(b)");
  gPad->SaveAs("Figures/Fig4b_FRateVsPt.png");
  gPad->SaveAs("Figures/Fig4b_FRateVsPt.eps");
  gPad->SaveAs("Figures/Fig4b_FRateVsPt.pdf");
  


  if(FORM==1)
    {
    
      grSigDQ0->Write();
      DissRate_Khar->Write();
      DissRate_Khar36->Write();
      DissRate_Khar86->Write();
      
      grDissRateVsPt1->Write();
      grDissRateVsPt2->Write();
      grDissRateVsPt3->Write();
      
      grFormRateVsT_MCInt_P1->Write();
      grFormRateVsT_MCInt_P2->Write();
      grFormRateVsT_MCInt_P3->Write();

      grFormRate_TripleInt_P1->Write();
      grFormRate_TripleInt_P2->Write();
      grFormRate_TripleInt_P3->Write();
      
      //OutFile->Write();
      //OutFile->Close();

    }



  //return;

 //==============================================================================================//
  //=================================  ALICE Calculations RAA Vs Pt  ===============================//
  //=============================================================================================//
  
  const int NNPtALICE =12;

  Double_t RAAALICEPt[NNPtALICE];
  Double_t PionDissALICEPt[NNPtALICE];
  Double_t ShadowDissALICEPt[NNPtALICE];
  Double_t GluonDissALICEPt[NNPtALICE];
  Double_t TotalDissALICEPt[NNPtALICE];
  Double_t ELossALICEPt[NNPtALICE];
  Double_t RAAALICEPt_ELoss[NNPtALICE];

  Double_t NJPsiRegenALICEPt[NNPtALICE];
  Double_t NJPsi0ALICEPt[NNPtALICE];
  Double_t FormationALICEPt[NNPtALICE];
  
  Double_t PtALICE[12]={0.25, 0.75, 1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 8.25, 10.25, 12.25, 16.25};

  Int_t BinMin = 0;
  Int_t BinMax = 40;

  if(QQbarVar==4 || QQbarVar==5 || QQbarVar==6 || QQbarVar==7){BinMax=180;}


  Double_t NPart0To20 = Npart(BinMin,BinMax);
  Double_t FacNaprtTemp0To20 = Npart(BinMin,BinMax)/Npart(0,10);

  Double_t ssALICE = ss05*(grdNDetaNpart->Eval(Npart(BinMin,BinMax))/grdNDetaNpart->Eval(Npart(0,10)));  
  Double_t R0ALICE = R05*TMath::Power(Npart(BinMin,BinMax)/Npart(0,10),0.5);
  
  Double_t NCCALICE=nCC0_21*NColl(BinMin,BinMax)/NColl(0,10);
  Double_t NJPsiALICE=nJpsi0_21*NColl(BinMin,BinMax)/NColl(0,10);
  
  cout<<endl<<endl;
  cout<<" ====== Calculating lattice EOS Bin Min "<< BinMin <<" BinMax "<<BinMax<<endl;
  NTau=0;
  CalculateTandf_LatticeEOS(ssALICE, R0ALICE);
  cout<<" NTau "<<NTau<<" ss "<< ssALICE <<"  "<<R0ALICE<<endl<<endl;

  //for(int i=0;i<NTau;i++){
  //cout<<Tau[i]<<"     "<<TempTau[i]<<"   "<<fQGP[i]<<endl;
  //}
 
 cout<<" Resetting the histograms in pT loop : "<<endl;
 for(int i=0;i<200;i++)
   {
     HistRegenJpsiPt[i]->Reset();
   }
 
 timer.Start();
 for(int i=0;i<NTau;i++){
   FormRateMC_P_New(TempTau[i],i);  
 }
 timer.Stop();
 
 cout<<" calculated formation in ALICE pT cpu sec: "<< timer.CpuTime()<<" real sec: "<<timer.RealTime()<<endl;
 cout<<" Calculating ALICE RAA Vs Pt "<<endl;
 
 //return;
 // NNPtALICE=1;
 
 
 new TCanvas;
 Jpsi_Pt->Draw();
 
 cout<<" JPsi pT integral "<<Jpsi_Pt->Integral()<<endl;
 
 for(int i=0; i<NNPtALICE; i++) {
 //for(int i=0; i<0; i++) {
   
   PionDissALICEPt[i] = PionDiss_All(PtALICE[i]);
   //TotalDissALICEPt[i]=IntDiss_All(PtALICE[i]);
   //TotalDissALICEPt[i]=IntDiss_All(PtALICE[i],NPart0To20);
 

   if(QQbarVar==1){TotalDissALICEPt[i]= IntDiss_NPartInt(BinMin,BinMax,1,PtALICE[i]);}
   if(QQbarVar==4 || QQbarVar==5 || QQbarVar==6 || QQbarVar==7){TotalDissALICEPt[i]= IntDiss_NPartInt(BinMin,BinMax,1,PtALICE[i]);}

   cout<<" ALICE Pt : "<<PtALICE[i]<<"  Total Dissociation "<<TotalDissALICEPt[i]<<endl;
   ELossALICEPt[i]=ELossVsPt(PtALICE[i], tsallis_fun);

   if(QQbarVar==1)ShadowDissALICEPt[i]=gr_shadowing_2540->Eval(PtALICE[i]);
   if(QQbarVar==4)ShadowDissALICEPt[i]=grf_shadowing_502tev_25_y_40_Y1S->Eval(PtALICE[i]);
   if(QQbarVar==5 || QQbarVar==6 || QQbarVar==7)ShadowDissALICEPt[i]=grf_shadowing_502tev_25_y_40_Y2S->Eval(PtALICE[i]);

   GluonDissALICEPt[i] =  TotalDissALICEPt[i]/(PionDissALICEPt[i]*ShadowDissALICEPt[i]);
   
   //NJPsiRegenALICEPt[i] = GluonDissALICEPt[i] * NCCALICE * NCCALICE * IntFormVsPt(PtALICE[i], R0ALICE, NPart0To20);
   //NJPsi0ALICEPt[i]= Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(PtALICE[i]))/Jpsi_Pt->GetBinWidth(0);
   //FormationALICEPt[i]= NJPsiRegenALICEPt[i]/(NJPsiALICE*JPsiIntSh_21*NJPsi0ALICEPt[i]);
   
   Double_t NewFormationALICE =0.0;
   
   if(QQbarVar==1){NewFormationALICE = RegenratedQuarkoniaPt(BinMin, BinMax, PtALICE[i]);} 

   FormationALICEPt[i]= NewFormationALICE;

   cout<<" Formation  : Current (no nPart Int) "<<  FormationALICEPt[i] <<" new "<<NewFormationALICE<<endl;

   // Total RAA 
   RAAALICEPt[i]= TotalDissALICEPt[i] + FormationALICEPt[i]; 
   //cout<<" 9 "<<endl;
   
   RAAALICEPt_ELoss[i]= ShadowDissALICEPt[i]*ELossALICEPt[i] ; 
   
   //cout<<"PT ALICE "<<PtALICE[i]<<"  "<<NJPsiRegenALICEPt[i]<<"   "<<(NJPsiALICE*JPsiIntSh_1*NJPsi0ALICEPt[i])<<endl;
 }
 
 
 cout<<" shadow "<<"   "<< " gluon "<<"   "<< " Pion "<< "   "<< " Diss  "<< "   "<<" For "<< "  "<< " RAA " <<endl;  
 
 for(int i=0; i<NNPtALICE; i++){ 
   cout<<ShadowDissALICEPt[i]<<"   "<<GluonDissALICEPt[i]<<"   "<<PionDissALICEPt[i]<<"   "<<
     ShadowDissALICEPt[i]*PionDissALICEPt[i]*GluonDissALICEPt[i]<<"  "<<
     FormationALICEPt[i]<<"  "<<RAAALICEPt[i]<<endl;
 }


 //RegenratedQuarkoniaNPart(e_t NCC, Double_t NJPsi, Double_t R0, Double_t IntSh) 
 //cout<<" value of regenration : "<<RegenratedQuarkoniaNPart() 

 // RAA vs pT graphs
 TGraph *grRAAALICEPt_M = new TGraph(NNPtALICE,PtALICE,RAAALICEPt);
 grRAAALICEPt_M->SetName("grRAAALICEPt_M");
 grRAAALICEPt_M->SetTitle("grRAAALICEPt_M");
 grRAAALICEPt_M->SetLineWidth(2);
 grRAAALICEPt_M->SetLineColor(1);
 grRAAALICEPt_M->SetLineStyle(1);
 grRAAALICEPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grRAAALICEPt_M->GetYaxis()->SetTitle("R_{AA}");
 grRAAALICEPt_M->GetYaxis()->SetRangeUser(0.0,1.0);
 
 
 TGraph *grRAAALICEPt_ELoss_M = new TGraph(NNPtALICE,PtALICE,RAAALICEPt_ELoss);
 grRAAALICEPt_ELoss_M->SetName("grRAAALICEPt_ELoss_M");
 grRAAALICEPt_ELoss_M->SetTitle("grRAAALICEPt_ELoss_M");
 grRAAALICEPt_ELoss_M->SetLineWidth(2);
 grRAAALICEPt_ELoss_M->SetLineColor(1);
 grRAAALICEPt_ELoss_M->SetLineStyle(1);
 grRAAALICEPt_ELoss_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grRAAALICEPt_ELoss_M->GetYaxis()->SetTitle("R_{AA}");
 grRAAALICEPt_ELoss_M->GetYaxis()->SetRangeUser(0.0,1.0);


 //Pion diss vs pT graphs
 TGraph *grPionDissALICEPt_M = new TGraph(NNPtALICE,PtALICE,PionDissALICEPt);
 grPionDissALICEPt_M->SetName("grPionDissALICEPt_M");
 grPionDissALICEPt_M->SetTitle("grPionDissALICEPt_M");
 grPionDissALICEPt_M->SetLineWidth(2);
 grPionDissALICEPt_M->SetLineColor(8);
 grPionDissALICEPt_M->SetLineStyle(2);
 grPionDissALICEPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grPionDissALICEPt_M->GetYaxis()->SetTitle("#pi Dissociation");
 grPionDissALICEPt_M->GetYaxis()->SetRangeUser(0.0,1.0);
 
 //Shadow diss vs pT graphs
 TGraph *grShadowDissALICEPt_M = new TGraph(NNPtALICE,PtALICE,ShadowDissALICEPt);
 grShadowDissALICEPt_M->SetName("grShadowDissALICEPt_M");
 grShadowDissALICEPt_M->SetTitle("grShadowDissALICEPt_M");
 grShadowDissALICEPt_M->SetLineWidth(2);
 grShadowDissALICEPt_M->SetLineColor(6);
 grShadowDissALICEPt_M->SetLineStyle(2);
 grShadowDissALICEPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grShadowDissALICEPt_M->GetYaxis()->SetTitle("R_{AA}^{Shadowing}");
 grShadowDissALICEPt_M->GetYaxis()->SetRangeUser(0.0,1.0);
 
 
 //Gluon diss vs pT graphs
 TGraph *grGluonDissALICEPt_M = new TGraph(NNPtALICE,PtALICE,GluonDissALICEPt);
 grGluonDissALICEPt_M->SetName("grGluonDissALICEPt_M");
 grGluonDissALICEPt_M->SetTitle("grGluonDissALICEPt_M");
 grGluonDissALICEPt_M->SetLineWidth(2);
 grGluonDissALICEPt_M->SetLineColor(4);
 grGluonDissALICEPt_M->SetLineStyle(6);
 grGluonDissALICEPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grGluonDissALICEPt_M->GetYaxis()->SetTitle("gluon Dissociation");
 grGluonDissALICEPt_M->GetYaxis()->SetRangeUser(0.0,1.0);
 

 //ELoss diss vs pT graphs
 TGraph *grELossALICEPt_M = new TGraph(NNPtALICE,PtALICE,ELossALICEPt);
 grELossALICEPt_M->SetName("grELossALICEPt_M");
 grELossALICEPt_M->SetTitle("grELossALICEPt_M");
 grELossALICEPt_M->SetLineWidth(2);
 grELossALICEPt_M->SetLineColor(kRed+3);
 grELossALICEPt_M->SetLineStyle(2);
 grELossALICEPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grELossALICEPt_M->GetYaxis()->SetTitle("R_{AA}^{ELoss}");
 grELossALICEPt_M->GetYaxis()->SetRangeUser(0.0,3.0);



 //Regenration vs pT graphs
 TGraph *grRegenALICEPt_M = new TGraph(NNPtALICE,PtALICE,FormationALICEPt);
 grRegenALICEPt_M->SetName("grRegenALICEPt_M");
 grRegenALICEPt_M->SetTitle("grRegenALICEPt_M");
 grRegenALICEPt_M->SetLineWidth(2);
 grRegenALICEPt_M->SetLineColor(2);
 grRegenALICEPt_M->SetLineStyle(4);
 grRegenALICEPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grRegenALICEPt_M->GetYaxis()->SetTitle("Regenration");
 grRegenALICEPt_M->GetYaxis()->SetRangeUser(0.0,1.0);
 
 TLegend *leg_ALICEJPsiRaaPt = new TLegend(0.16,0.61,0.82,0.93);
 leg_ALICEJPsiRaaPt->SetBorderSize(0);
 leg_ALICEJPsiRaaPt->SetFillStyle(0);
 leg_ALICEJPsiRaaPt->SetFillColor(0);
 leg_ALICEJPsiRaaPt->SetTextSize(0.030);


 TGraphErrors *ALICE_GhostGraph = grf_GhostGraph(20, 6, 6);
 if(QQbarVar==1){leg_ALICEJPsiRaaPt->AddEntry(ALICE_GhostGraph,"ALICE Data, 2.5<y^{J/#psi}<4.0, 0-20%","P");}
 else{leg_ALICEJPsiRaaPt->AddEntry(ALICE_GhostGraph,"ALICE Prediction, 2.5<y^{#varUpsilon}<4.0, 0-100%","P");}

 leg_ALICEJPsiRaaPt->AddEntry(grGluonDissALICEPt_M,"Gluon Dissociation","L");
 leg_ALICEJPsiRaaPt->AddEntry(grRegenALICEPt_M,"Formation","L");
 leg_ALICEJPsiRaaPt->AddEntry(grPionDissALICEPt_M,"Comover","L");
 leg_ALICEJPsiRaaPt->AddEntry(grShadowDissALICEPt_M,"CNM Effects","L");
 leg_ALICEJPsiRaaPt->AddEntry(grRAAALICEPt_M,"Total (R_{AA})","L");
 



 TLegend *leg_ALICEJPsiRaaPt_ELoss = new TLegend(0.16,0.61,0.82,0.93);
 leg_ALICEJPsiRaaPt_ELoss->SetBorderSize(0);
 leg_ALICEJPsiRaaPt_ELoss->SetFillStyle(0);
 leg_ALICEJPsiRaaPt_ELoss->SetFillColor(0);
 leg_ALICEJPsiRaaPt_ELoss->SetTextSize(0.030);

 if(QQbarVar==1){leg_ALICEJPsiRaaPt_ELoss->AddEntry(ALICE_GhostGraph,"ALICE Data, 2.5<y<4.0, 0-20%","P");}
 else{leg_ALICEJPsiRaaPt_ELoss->AddEntry(ALICE_GhostGraph,"ALICE Prediction, 2.5<y^{#varUpsilon}<4.0, 0-100%","P");}
 leg_ALICEJPsiRaaPt_ELoss->AddEntry(grShadowDissALICEPt_M,"CNM Effects","L");
 leg_ALICEJPsiRaaPt_ELoss->AddEntry(grELossALICEPt_M,"E Loss","L");
 leg_ALICEJPsiRaaPt_ELoss->AddEntry(grRAAALICEPt_M,"Total (R_{AA})","L");
 

 new TCanvas;
 gPad->SetTicks();
 if(QQbarVar==1){Draw_ALICE_JPsi_RaaVsPt();}
 if(QQbarVar==4 || QQbarVar==5 || QQbarVar==6 || QQbarVar==7 ){
   grRAAALICEPt_M->GetYaxis()->SetRangeUser(0,2.0);
   grRAAALICEPt_M->Draw("AC"); }
 
 grRAAALICEPt_M->Draw("Csame");
 grPionDissALICEPt_M->Draw("Csame");
 grShadowDissALICEPt_M->Draw("Csame");
 grGluonDissALICEPt_M->Draw("Csame");
 grRegenALICEPt_M->Draw("Csame");
 leg_ALICEJPsiRaaPt->Draw("same"); 
 gPad->SaveAs("Figures/Fig_ALICE_RAAPt.pdf");
 gPad->SaveAs("Figures/Fig_ALICE_RAAPt.png");
 
 

 new TCanvas;
 gPad->SetTicks();
 if(QQbarVar==1){Draw_ALICE_JPsi_RaaVsPt();}
 if(QQbarVar==4 || QQbarVar==5 || QQbarVar==6 || QQbarVar==7 ){
   grRAAALICEPt_M->GetYaxis()->SetRangeUser(0,2.0);
   grRAAALICEPt_M->Draw("AC"); }
 
 grRAAALICEPt_ELoss_M->Draw("Csame");
 grShadowDissALICEPt_M->Draw("Csame");
 grELossALICEPt_M->Draw("Csame");
 leg_ALICEJPsiRaaPt_ELoss->Draw("same"); 
 gPad->SaveAs("Figures/Fig_ALICE_RAAPt_ELoss.pdf");
 gPad->SaveAs("Figures/Fig_ALICE_RAAPt_ELoss.png");

 cout<<endl<<endl<<endl;
 
 

 //return;



  //==============================================================================================//
  //================================= CMS RAA Vs Pt Calculations ===================================//
  //=============================================================================================//
  
  //int NNPtCMS = 37;

  int NNPtCMS = 15;

  Double_t PionDissCMSPt[37]={0.0};
  Double_t ShadowDissCMSPt[37]={0.0};
  Double_t GluonDissCMSPt[37]={0.0};
  Double_t ELossCMSPt[37]={0.0};
  Double_t TotalDissCMSPt[37]={0.0};

  Double_t RAACMSPt[37]={0.0};
  Double_t RAACMSPt_ELoss[37]={0.0};

  Double_t NJPsiRegenCMSPt[37]={0.0};
  Double_t NJPsi0CMSPt[37]={0.0};
  Double_t FormationCMSPt[37]={0.0};
 

  //Double_t PtCMS[37]={5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,10.25,11.25,12.25,13.25,14.25,15.25,16.25,17.25,18.25,19.25,20.25,21.25,22.25, 23.25,
  //24.25,25.25,26.25,27.25,28.25,29.25,30.25,31.25,33.25,37.25,40.25,43.25,46.25,49.25};



  Double_t PtCMS[15]={5.25,6.25,7.25,8.25,10.25,12.25,14.25,16.25,20.25,24.25,28.25,32.25,36.25,44.25,49.25};
  

  int NNPtCMS_Upsilon = 26;
  Double_t PtCMS_Upsilon[26]={0.25,1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25,9.25,10.25,11.25,12.25,13.25,14.25,
			       15.25,16.25,17.25,18.25,19.25,20.25,21.25,22.25,23.25,26.25,30.25};

  //Double_t FacNaprtTempMB = Npart(0,200)/Npart(0,10);
  
  Double_t NPartMB = Npart(0,200);
  Double_t ssMB = ss05*(grdNDetaNpart->Eval(Npart(0,200))/grdNDetaNpart->Eval(Npart(0,10)));  
  //Double_t ssMB = ss05*Npart(0,200)/Npart(0,10);  
  Double_t R0MB = R05*TMath::Power(Npart(0,200)/Npart(0,10),0.5); 
  
  cout<<endl<<endl;
  cout<<" ====== Calculating lattice EOS MB "<<endl;
  NTau=0;
  CalculateTandf_LatticeEOS(ssMB, R0MB);
  
  cout<<" NTau "<<NTau<<" ssMB "<< ssMB <<endl;
  
  //for(int i=0;i<NTau;i++){
  //cout<<Tau[i]<<"     "<<TempTau[i]<<"   "<<fQGP[i]<<endl;
  //}

  

  //================== regenration at CMS =============================//
  Double_t NCCCMS=nCC0_21*NColl(0,200)/NColl(0,10);
  Double_t NJPsiCMS=nJpsi0_21*NColl(0,200)/NColl(0,10);
 
  cout<<" Resetting the histograms in pT loop : "<<endl;
  for(int i=0;i<200;i++)
    {
      HistRegenJpsiPt[i]->Reset();
    }
 
  timer.Start();
  for(int i=0;i<NTau;i++){
    FormRateMC_P_New(TempTau[i],i);  
  }
  timer.Stop();
  
  cout<<" CMS:calculated formation in 0-100% in cpu sec: "<< timer.CpuTime()<<" real sec: "<<timer.RealTime()<<endl;
  


  if(QQbarVar==1){cout<<" J/psi CMS RAA pT Calculations "<<endl;}

  if(QQbarVar==4){cout<<" Y(1S) CMS RAA pT Calculations "<<endl;}
  if(QQbarVar==5){cout<<" Y(2S) CMS RAA pT Calculations "<<endl;}
  
  cout<<"  Pt  "<<"  Pion     "<<"    CNM     "<<"    Gluon     Formation "<<" ELoss      RAA     "<<endl;

  if(QQbarVar == 4 || QQbarVar == 5 || QQbarVar== 6 || QQbarVar == 7 ){NNPtCMS=NNPtCMS_Upsilon;}

  for(int i=0; i<NNPtCMS; i++) {
  //for(int i=0; i<0; i++) {
    
    if(QQbarVar==4 || QQbarVar ==5 || QQbarVar== 6 || QQbarVar == 7 ){PtCMS[i] = PtCMS_Upsilon[i];}
    
    PionDissCMSPt[i]= PionDiss_All(PtCMS[i]);
    
    //TotalDissCMSPt[i]=IntDiss_All(PtCMS[i]);
    //TotalDissCMSPt[i]=IntDiss_All(PtCMS[i],NPartMB);
    
    TotalDissCMSPt[i]= IntDiss_NPartInt(0,200,2,PtCMS[i]);
    
    ELossCMSPt[i]=ELossVsPt(PtCMS[i], tsallis_fun);
    
    ShadowDissCMSPt[i]=1.0;
    
    if(QQbarVar ==1){ShadowDissCMSPt[i]=gr_shadowing_2424->Eval(PtCMS[i]);}
    
    if(QQbarVar ==4){ShadowDissCMSPt[i]=grf_shadowing_502tev_24_y_24_Y1S->Eval(PtCMS[i]);}
    if(QQbarVar ==5){ShadowDissCMSPt[i]=grf_shadowing_502tev_24_y_24_Y2S->Eval(PtCMS[i]);}
    if(QQbarVar ==6){ShadowDissCMSPt[i]=grf_shadowing_502tev_24_y_24_Y2S->Eval(PtCMS[i]);}
    if(QQbarVar ==7){ShadowDissCMSPt[i]=grf_shadowing_502tev_24_y_24_Y2S->Eval(PtCMS[i]);}
    
    GluonDissCMSPt[i]=TotalDissCMSPt[i]/(PionDissCMSPt[i]*ShadowDissCMSPt[i]);


    if(PtCMS[i]<20.0){
      
      //NJPsiRegenCMSPt[i]= GluonDissCMSPt[i]*NCCCMS*NCCCMS*IntFormVsPt(PtCMS[i],R0MB,NPartMB);
      //NJPsi0CMSPt[i]= Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(PtCMS[i]))/Jpsi_Pt->GetBinWidth(0);
      //FormationCMSPt[i]= NJPsiRegenCMSPt[i]/(NJPsiCMS*JPsiIntSh_21*NJPsi0CMSPt[i]);
    
      Double_t NewFormationCMS = RegenratedQuarkoniaPt(0, 200, PtCMS[i]); 
      FormationCMSPt[i]= NewFormationCMS;

      cout<<" pT "<<PtCMS[i]<<"  "<<NewFormationCMS<<endl;
    }
   // Total RAA 
   RAACMSPt[i]= ShadowDissCMSPt[i]*PionDissCMSPt[i]*GluonDissCMSPt[i] + FormationCMSPt[i]; 
   //cout<<" 9 "<<endl;

   //RAACMSPt[i]= ShadowDissCMSPt[i]*PionDissCMSPt[i]*GluonDissCMSPt[i];

    RAACMSPt_ELoss[i]=ShadowDissCMSPt[i]*ELossCMSPt[i];
    cout<<PtCMS[i]<<"   "<<PionDissCMSPt[i]<<"   "<<ShadowDissCMSPt[i]<<"   "<<GluonDissCMSPt[i]<<"   "<<FormationCMSPt[i]<<"  "<<ELossCMSPt[i]<<"   "<<RAACMSPt[i]<<endl;

}



  // RAA vs pT graphs
  TGraph *grRAACMSPt_M = new TGraph(NNPtCMS,PtCMS,RAACMSPt);
  grRAACMSPt_M->SetName("grRAACMSPt_M");
  grRAACMSPt_M->SetTitle("grRAACMSPt_M");
  grRAACMSPt_M->SetLineWidth(2);
  grRAACMSPt_M->SetLineColor(1);
  grRAACMSPt_M->SetLineStyle(1);
  grRAACMSPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRAACMSPt_M->GetYaxis()->SetTitle("R_{AA}");
  grRAACMSPt_M->GetYaxis()->SetRangeUser(0.0,3.0);


  // RAA vs pT graphs (ELoss)
  TGraph *grRAACMSPt_ELoss_M = new TGraph(NNPtCMS,PtCMS,RAACMSPt_ELoss);
  grRAACMSPt_ELoss_M->SetName("grRAACMSPt_ELoss_M");
  grRAACMSPt_ELoss_M->SetTitle("grRAACMSPt_ELoss_M");
  grRAACMSPt_ELoss_M->SetLineWidth(2);
  grRAACMSPt_ELoss_M->SetLineColor(1);
  grRAACMSPt_ELoss_M->SetLineStyle(1);
  grRAACMSPt_ELoss_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRAACMSPt_ELoss_M->GetYaxis()->SetTitle("R_{AA}");
  grRAACMSPt_ELoss_M->GetYaxis()->SetRangeUser(0.0,3.0);
  



  //Gluon diss vs pT graphs
  TGraph *grGluonDissCMSPt_M = new TGraph(NNPtCMS,PtCMS,GluonDissCMSPt);
  grGluonDissCMSPt_M->SetName("grGluonDissCMSPt_M");
  grGluonDissCMSPt_M->SetTitle("grGluonDissCMSPt_M");
  grGluonDissCMSPt_M->SetLineWidth(2);
  grGluonDissCMSPt_M->SetLineColor(4);
  grGluonDissCMSPt_M->SetLineStyle(6);
  grGluonDissCMSPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grGluonDissCMSPt_M->GetYaxis()->SetTitle("gluon Dissociation");
  grGluonDissCMSPt_M->GetYaxis()->SetRangeUser(0.0,3.0);
  
  
  //Pion diss vs pT graphs
  TGraph *grPionDissCMSPt_M = new TGraph(NNPtCMS,PtCMS,PionDissCMSPt);
  grPionDissCMSPt_M->SetName("grPionDissCMSPt_M");
  grPionDissCMSPt_M->SetTitle("grPionDissCMSPt_M");
  grPionDissCMSPt_M->SetLineWidth(2);
  grPionDissCMSPt_M->SetLineColor(8);
  grPionDissCMSPt_M->SetLineStyle(2);
  grPionDissCMSPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grPionDissCMSPt_M->GetYaxis()->SetTitle("#pi Dissociation");
  grPionDissCMSPt_M->GetYaxis()->SetRangeUser(0.0,3.0);

  //Shadow diss vs pT graphs
  TGraph *grShadowDissCMSPt_M = new TGraph(NNPtCMS,PtCMS,ShadowDissCMSPt);
  grShadowDissCMSPt_M->SetName("grShadowDissCMSPt_M");
  grShadowDissCMSPt_M->SetTitle("grShadowDissCMSPt_M");
  grShadowDissCMSPt_M->SetLineWidth(2);
  grShadowDissCMSPt_M->SetLineColor(6);
  grShadowDissCMSPt_M->SetLineStyle(2);
  grShadowDissCMSPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grShadowDissCMSPt_M->GetYaxis()->SetTitle("R_{AA}^{Shadowing}");
  grShadowDissCMSPt_M->GetYaxis()->SetRangeUser(0.0,3.0);

  
  //ELoss diss vs pT graphs
  TGraph *grELossCMSPt_M = new TGraph(NNPtCMS,PtCMS,ELossCMSPt);
  grELossCMSPt_M->SetName("grELossCMSPt_M");
  grELossCMSPt_M->SetTitle("grELossCMSPt_M");
  grELossCMSPt_M->SetLineWidth(2);
  grELossCMSPt_M->SetLineColor(kRed+3);
  grELossCMSPt_M->SetLineStyle(2);
  grELossCMSPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grELossCMSPt_M->GetYaxis()->SetTitle("R_{AA}^{ELoss}");
  grELossCMSPt_M->GetYaxis()->SetRangeUser(0.0,3.0);

  
  //Regenration vs pT graphs
  TGraph *grRegenCMSPt_M = new TGraph(NNPtCMS,PtCMS,FormationCMSPt);
  grRegenCMSPt_M->SetName("grRegenCMSPt_M");
  grRegenCMSPt_M->SetTitle("grRegenCMSPt_M");
  grRegenCMSPt_M->SetLineWidth(2);
  grRegenCMSPt_M->SetLineColor(2);
  grRegenCMSPt_M->SetLineStyle(4);
  grRegenCMSPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRegenCMSPt_M->GetYaxis()->SetTitle("Formation");
  grRegenCMSPt_M->GetYaxis()->SetRangeUser(0.0,3.0);

  TLegend *leg_CMSJPsiRaaPt = new TLegend(0.15,0.60,0.90,0.93);
  leg_CMSJPsiRaaPt->SetBorderSize(0);
  leg_CMSJPsiRaaPt->SetFillStyle(0);
  leg_CMSJPsiRaaPt->SetFillColor(0);
  leg_CMSJPsiRaaPt->SetTextSize(0.040);
     

  //TGraphErrors *grf_GhostGraph(Int_t MarkerStyle, Int_t MarkerColor, Int_t LineColor)
  TGraphErrors *CMS_GhostGraph = grf_GhostGraph(21, 2, 2);

  leg_CMSJPsiRaaPt->AddEntry(CMS_GhostGraph,"CMS Data |y|<2.4, 0-100%","P");
  leg_CMSJPsiRaaPt->AddEntry(grGluonDissCMSPt_M,"Gluon Dissociation","L");
  leg_CMSJPsiRaaPt->AddEntry(grPionDissCMSPt_M,"Comover","L");
  leg_CMSJPsiRaaPt->AddEntry(grRegenCMSPt_M,"Formation","L");
  leg_CMSJPsiRaaPt->AddEntry(grShadowDissCMSPt_M,"CNM Effects","L");
  leg_CMSJPsiRaaPt->AddEntry(grRAACMSPt_M,"Total (R_{AA})","L");


  new TCanvas;
  gPad->SetTicks();
  if(QQbarVar==1) Draw_CMS_JPsi_RaaVsPt_0To100();
  
  if(QQbarVar==4)Draw_CMS_Y1S_5TeV_RaaVsPt();
  if(QQbarVar==5)Draw_CMS_Y2S_5TeV_RaaVsPt();
  if(QQbarVar==6)Draw_CMS_Y2S_5TeV_RaaVsPt();
  if(QQbarVar==7)Draw_CMS_Y2S_5TeV_RaaVsPt();

  grRAACMSPt_M->Draw("Lsame");
  grPionDissCMSPt_M->Draw("Lsame");
  grShadowDissCMSPt_M->Draw("Lsame");
  grGluonDissCMSPt_M->Draw("Lsame");
  grRegenCMSPt_M->Draw("Lsame");
  leg_CMSJPsiRaaPt->Draw("same"); 
  
  gPad->SaveAs("Figures/Fig_CMS_JPsiRAAPt.pdf");
  gPad->SaveAs("Figures/Fig_CMS_JPsiRAAPt.png");


  TLegend *leg_CMSJPsiRaaPt_ELoss = new TLegend(0.15,0.60,0.90,0.93);
  leg_CMSJPsiRaaPt_ELoss->SetBorderSize(0);
  leg_CMSJPsiRaaPt_ELoss->SetFillStyle(0);
  leg_CMSJPsiRaaPt_ELoss->SetFillColor(0);
  leg_CMSJPsiRaaPt_ELoss->SetTextSize(0.040);
  
  leg_CMSJPsiRaaPt_ELoss->AddEntry(CMS_GhostGraph,"CMS Data |y|<2.4, 0-100%","P");
  leg_CMSJPsiRaaPt_ELoss->AddEntry(grShadowDissCMSPt_M,"CNM Effects","L");
  leg_CMSJPsiRaaPt_ELoss->AddEntry(grELossCMSPt_M,"E-Loss Effect","L");
  leg_CMSJPsiRaaPt_ELoss->AddEntry(grRAACMSPt_M,"Total (R_{AA})","L");
  
  new TCanvas;
  gPad->SetTicks();
  if(QQbarVar==1) Draw_CMS_JPsi_RaaVsPt_0To100();
  if(QQbarVar==4)Draw_CMS_Y1S_5TeV_RaaVsPt();
  if(QQbarVar==5)Draw_CMS_Y2S_5TeV_RaaVsPt();
  if(QQbarVar==6)Draw_CMS_Y2S_5TeV_RaaVsPt();
  if(QQbarVar==7)Draw_CMS_Y2S_5TeV_RaaVsPt();

  grRAACMSPt_ELoss_M->Draw("Lsame");
  grShadowDissCMSPt_M->Draw("Lsame");
  grELossCMSPt_M->Draw("Lsame");
  leg_CMSJPsiRaaPt_ELoss->Draw("same"); 

  gPad->SaveAs("Figures/Fig_CMS_JPsiRAAPt_ELoss.pdf");
  gPad->SaveAs("Figures/Fig_CMS_JPsiRAAPt_ELoss.png");


  

  cout<<endl<<endl<<endl;

  TGraph *grRAAATLASPt_M;
  TGraph *grRAAATLASPt_ELoss_M;
  TGraph *grGluonDissATLASPt_M;
  TGraph *grPionDissATLASPt_M;
  TGraph *grShadowDissATLASPt_M;
  TGraph *grELossATLASPt_M;
  TGraph *grRegenATLASPt_M;
  
  if(QQbarVar==1)
    {
      //==============================================================================================//
      //================================= ATLAS RAA Vs Pt Calculations ===================================//
      //=============================================================================================//
    

      // const int NNPtATLAS = 26;

      const int NNPtATLAS = 12;

  
      Double_t PionDissATLASPt[NNPtATLAS];
      Double_t ShadowDissATLASPt[NNPtATLAS];
      Double_t GluonDissATLASPt[NNPtATLAS];
      Double_t ELossATLASPt[NNPtATLAS];
      Double_t TotalDissATLASPt[NNPtATLAS];

      Double_t RAAATLASPt[NNPtATLAS];
      Double_t RAAATLASPt_ELoss[NNPtATLAS];
      
      Double_t NJPsiRegenATLASPt[NNPtATLAS];
      Double_t NJPsi0ATLASPt[NNPtATLAS];
      Double_t FormationATLASPt[NNPtATLAS];
      
      //Double_t PtATLAS[26]={9.25,10.25,11.25,12.25,13.25,14.25,15.25,16.25,17.25,18.25,19.25,20.25,21.25,22.25,
      //23.25,24.25,25.25,26.25,27.25,28.25,29.25,30.25,31.25,33.25,37.25,40.25};


      Double_t PtATLAS[12]={9.25,11.25,13.25,15.25,17.25,19.25,22.25,25.25,28.25,31.25,37.25,40.25};


      Double_t NPart0To80 = Npart(0,160);
      Double_t FacNaprtTemp0To80 = Npart(0,160)/Npart(0,10);
      Double_t ss0To80 = ss05*(grdNDetaNpart->Eval(Npart(0,160))/grdNDetaNpart->Eval(Npart(0,10)));  
      //Double_t ss0To80 = ss05*Npart(0,160)/Npart(0,10);  
      Double_t R00To80 = R05*TMath::Power(Npart(0,160)/Npart(0,10),0.5); 
      
      
      cout<<endl<<endl;
      cout<<" ====== Calculating lattice EOS 0To80 "<<endl;
      NTau=0;
      CalculateTandf_LatticeEOS(ss0To80, R00To80);
      
      cout<<" NTau "<<NTau<<" ss0To80 "<< ss0To80 <<" R00To80 "<<R00To80<<endl;
      
      //for(int i=0;i<NTau;i++){
      //cout<<Tau[i]<<"     "<<TempTau[i]<<"   "<<fQGP[i]<<endl;
      //}
      

      //================== regenration at ATLAS =============================//
      Double_t NCCATLAS=nCC0_21*NColl(0,160)/NColl(0,10);
      Double_t NJPsiATLAS=nJpsi0_21*NColl(0,160)/NColl(0,10);
      
      cout<<" ATLAS: Resetting the histograms in pT loop : "<<endl;
      for(int i=0;i<200;i++)
	{
	  HistRegenJpsiPt[i]->Reset();
	}
      
      timer.Start();
      for(int i=0;i<NTau;i++){
	FormRateMC_P_New(TempTau[i],i);  
      }
      timer.Stop();
      
      cout<<" ATLAS:calculated formation in 0-80% in cpu sec: "<< timer.CpuTime()<<" real sec: "<<timer.RealTime()<<endl;
      cout<<" J/psi ATLAS RAA pT Calculations for Npart :"<<NPart0To80<<endl;
      
         
      for(int i=0; i<NNPtATLAS; i++) {
	//for(int i=0; i<0; i++) {
	
	PionDissATLASPt[i]=PionDiss_All(PtATLAS[i]);
	

	//TotalDissATLASPt[i]=IntDiss_All(PtATLAS[i]);
	//TotalDissATLASPt[i]= IntDiss_All(PtATLAS[i],NPart0To80);
	TotalDissATLASPt[i]= IntDiss_NPartInt(0,160,3,PtATLAS[i]);
	
	ELossATLASPt[i]=ELossVsPt(PtATLAS[i], tsallis_fun);
	
	ShadowDissATLASPt[i]=1.0;
	
	if(QQbarVar ==1){ShadowDissATLASPt[i]=gr_shadowing_2424->Eval(PtATLAS[i]);}
	GluonDissATLASPt[i]= TotalDissATLASPt[i]/(PionDissATLASPt[i]*ShadowDissATLASPt[i]);

	if(PtATLAS[i]<20.0){
	  
	  //Double_t ForInt = IntFormVsPt(PtATLAS[i], R00To80, NPart0To80);
	  //NJPsiRegenATLASPt[i]=GluonDissATLASPt[i]*NCCATLAS*NCCATLAS*ForInt;
	  //NJPsi0ATLASPt[i]= Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(PtATLAS[i]))/Jpsi_Pt->GetBinWidth(0);
	  //FormationATLASPt[i]= NJPsiRegenATLASPt[i]/(NJPsiATLAS*JPsiIntSh_21*NJPsi0ATLASPt[i]);

	  Double_t NewFormationATLAS = RegenratedQuarkoniaPt(0, 160, PtATLAS[i]); 
	  FormationATLASPt[i]= NewFormationATLAS;

	  //cout<<GluonDissATLASPt[i]<<"   "<<NCCATLAS<<"   "<<ForInt<<endl;
	  //cout<<" NJPsiRegenATLASPt "<<NJPsiRegenATLASPt[i]<<"   "<<NJPsi0ATLASPt[i]<<"   "<<FormationATLASPt[i]<<endl;
	}

	// Total RAA 
	RAAATLASPt[i]= TotalDissATLASPt[i] + FormationATLASPt[i]; 
	RAAATLASPt_ELoss[i]=ShadowDissATLASPt[i]*ELossATLASPt[i];
	
      }
      
      
      cout<<"  Pt  "<<"  Pion     "<<"    CNM     "<<"    Gluon   "<<"    Formation   "<<"   ELoss    "<<"   RAA     "<<endl;  
      for(int i=0; i<NNPtATLAS; i++) {
	cout<<PtATLAS[i]<<"   "<<PionDissATLASPt[i]<<"   "<<ShadowDissATLASPt[i]<<"   "<<GluonDissATLASPt[i]<<"   "<<
	  FormationATLASPt[i]<<"   "<<ELossATLASPt[i]<<"   "<<RAAATLASPt[i]<<endl;
      }
      
      // RAA vs pT graphs
      grRAAATLASPt_M = new TGraph(NNPtATLAS,PtATLAS,RAAATLASPt);
      grRAAATLASPt_M->SetName("grRAAATLASPt_M");
      grRAAATLASPt_M->SetTitle("grRAAATLASPt_M");
      grRAAATLASPt_M->SetLineWidth(2);
      grRAAATLASPt_M->SetLineColor(1);
      grRAAATLASPt_M->SetLineStyle(1);
      grRAAATLASPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      grRAAATLASPt_M->GetYaxis()->SetTitle("R_{AA}");
      grRAAATLASPt_M->GetYaxis()->SetRangeUser(0.0,3.0);
      
      

      // RAA vs pT graphs (ELoss)
      grRAAATLASPt_ELoss_M = new TGraph(NNPtATLAS,PtATLAS,RAAATLASPt_ELoss);
      grRAAATLASPt_ELoss_M->SetName("grRAAATLASPt_ELoss_M");
      grRAAATLASPt_ELoss_M->SetTitle("grRAAATLASPt_ELoss_M");
      grRAAATLASPt_ELoss_M->SetLineWidth(2);
      grRAAATLASPt_ELoss_M->SetLineColor(1);
      grRAAATLASPt_ELoss_M->SetLineStyle(1);
      grRAAATLASPt_ELoss_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      grRAAATLASPt_ELoss_M->GetYaxis()->SetTitle("R_{AA}");
      grRAAATLASPt_ELoss_M->GetYaxis()->SetRangeUser(0.0,3.0);
      
     

      //Gluon diss vs pT graphs
      grGluonDissATLASPt_M = new TGraph(NNPtATLAS,PtATLAS,GluonDissATLASPt);
      grGluonDissATLASPt_M->SetName("grGluonDissATLASPt_M");
      grGluonDissATLASPt_M->SetTitle("grGluonDissATLASPt_M");
      grGluonDissATLASPt_M->SetLineWidth(2);
      grGluonDissATLASPt_M->SetLineColor(4);
      grGluonDissATLASPt_M->SetLineStyle(6);
      grGluonDissATLASPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      grGluonDissATLASPt_M->GetYaxis()->SetTitle("gluon Dissociation");
      grGluonDissATLASPt_M->GetYaxis()->SetRangeUser(0.0,3.0);
      
  
      //Pion diss vs pT graphs
      grPionDissATLASPt_M = new TGraph(NNPtATLAS,PtATLAS,PionDissATLASPt);
      grPionDissATLASPt_M->SetName("grPionDissATLASPt_M");
      grPionDissATLASPt_M->SetTitle("grPionDissATLASPt_M");
      grPionDissATLASPt_M->SetLineWidth(2);
      grPionDissATLASPt_M->SetLineColor(8);
      grPionDissATLASPt_M->SetLineStyle(2);
      grPionDissATLASPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      grPionDissATLASPt_M->GetYaxis()->SetTitle("#pi Dissociation");
      grPionDissATLASPt_M->GetYaxis()->SetRangeUser(0.0,3.0);
      
      //Shadow diss vs pT graphs
      grShadowDissATLASPt_M = new TGraph(NNPtATLAS,PtATLAS,ShadowDissATLASPt);
      grShadowDissATLASPt_M->SetName("grShadowDissATLASPt_M");
      grShadowDissATLASPt_M->SetTitle("grShadowDissATLASPt_M");
      grShadowDissATLASPt_M->SetLineWidth(2);
      grShadowDissATLASPt_M->SetLineColor(6);
      grShadowDissATLASPt_M->SetLineStyle(2);
      grShadowDissATLASPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      grShadowDissATLASPt_M->GetYaxis()->SetTitle("R_{AA}^{Shadowing}");
      grShadowDissATLASPt_M->GetYaxis()->SetRangeUser(0.0,3.0);
      
      
      

      //ELoss diss vs pT graphs
      grELossATLASPt_M = new TGraph(NNPtATLAS,PtATLAS,ELossATLASPt);
      grELossATLASPt_M->SetName("grELossATLASPt_M");
      grELossATLASPt_M->SetTitle("grELossATLASPt_M");
      grELossATLASPt_M->SetLineWidth(2);
      grELossATLASPt_M->SetLineColor(kRed+3);
      grELossATLASPt_M->SetLineStyle(2);
      grELossATLASPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      grELossATLASPt_M->GetYaxis()->SetTitle("R_{AA}^{ELoss}");
      grELossATLASPt_M->GetYaxis()->SetRangeUser(0.0,3.0);
      
      //Regenration vs pT graphs
      grRegenATLASPt_M = new TGraph(NNPtATLAS,PtATLAS,FormationATLASPt);
      grRegenATLASPt_M->SetName("grRegenATLASPt_M");
      grRegenATLASPt_M->SetTitle("grRegenATLASPt_M");
      grRegenATLASPt_M->SetLineWidth(2);
      grRegenATLASPt_M->SetLineColor(2);
      grRegenATLASPt_M->SetLineStyle(4);
      grRegenATLASPt_M->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      grRegenATLASPt_M->GetYaxis()->SetTitle("Formation");
      grRegenATLASPt_M->GetYaxis()->SetRangeUser(0.0,3.0);
      
      TLegend *leg_ATLASJPsiRaaPt = new TLegend(0.15,0.69,0.90,0.92);
      leg_ATLASJPsiRaaPt->SetBorderSize(0);
      leg_ATLASJPsiRaaPt->SetFillStyle(0);
      leg_ATLASJPsiRaaPt->SetFillColor(0);
      leg_ATLASJPsiRaaPt->SetTextSize(0.040);
      

      TGraphErrors *ATLAS_GhostGraph = grf_GhostGraph(24, 4, 4);
      
      leg_ATLASJPsiRaaPt->AddEntry(ATLAS_GhostGraph,"ATLAS Data |y|<2.0, 0-80%","P");
      leg_ATLASJPsiRaaPt->AddEntry(grGluonDissATLASPt_M,"Gluon Dissociation","L");
      leg_ATLASJPsiRaaPt->AddEntry(grPionDissATLASPt_M,"Comover","L");
      leg_ATLASJPsiRaaPt->AddEntry(grRegenATLASPt_M,"Formation","L");
      leg_ATLASJPsiRaaPt->AddEntry(grShadowDissATLASPt_M,"CNM Effects","L");
      leg_ATLASJPsiRaaPt->AddEntry(grRAAATLASPt_M,"Total (R_{AA})","L");
      
      new TCanvas;
      gPad->SetTicks();
      Draw_ATLAS_JPsi_RaaVsPt();
      grRAAATLASPt_M->Draw("Lsame");
      grPionDissATLASPt_M->Draw("Lsame");
      grShadowDissATLASPt_M->Draw("Lsame");
      grGluonDissATLASPt_M->Draw("Lsame");
      grRegenATLASPt_M->Draw("Lsame");
      leg_ATLASJPsiRaaPt->Draw("same"); 
      
      gPad->SaveAs("Figures/Fig_ATLAS_JPsiRAAPt.pdf");
      gPad->SaveAs("Figures/Fig_ATLAS_JPsiRAAPt.png");
      
    
      TLegend *leg_ATLASJPsiRaaPt_ELoss = new TLegend(0.15,0.69,0.90,0.92);
      leg_ATLASJPsiRaaPt_ELoss->SetBorderSize(0);
      leg_ATLASJPsiRaaPt_ELoss->SetFillStyle(0);
      leg_ATLASJPsiRaaPt_ELoss->SetFillColor(0);
      leg_ATLASJPsiRaaPt_ELoss->SetTextSize(0.040);
      
      leg_ATLASJPsiRaaPt_ELoss->AddEntry(ATLAS_GhostGraph,"ATLAS Data |y|<2.0, 0-80%","P");
      leg_ATLASJPsiRaaPt_ELoss->AddEntry(grShadowDissATLASPt_M,"CNM Effects","L");
      leg_ATLASJPsiRaaPt_ELoss->AddEntry(grELossATLASPt_M,"E-Loss Effect","L");
      leg_ATLASJPsiRaaPt_ELoss->AddEntry(grRAAATLASPt_M,"Total (R_{AA})","L");


  
      new TCanvas;
      gPad->SetTicks();
      Draw_ATLAS_JPsi_RaaVsPt();
      grRAAATLASPt_ELoss_M->Draw("Lsame");
      grShadowDissATLASPt_M->Draw("Lsame");
      grELossATLASPt_M->Draw("Lsame");
      leg_ATLASJPsiRaaPt_ELoss->Draw("same"); 
      gPad->SaveAs("Figures/Fig_ATLAS_JPsiRAAPt_ELoss.pdf");
      gPad->SaveAs("Figures/Fig_ATLAS_JPsiRAAPt_ELoss.png");
      
      cout<<endl<<endl<<endl;
      
    }


 
  //return;

 //========================================================================================//
 //=======================================================================================//
 // ======================== Calculations for centrality ===========================//
 //=====================================================================================//
 //====================================================================================//
  
  //Double_t nPartMin = 0.0;
  //Double_t nPartMax = 500;
  //Double_t nstep = 50;
  //int NN = (int)((nPartMax-nPartMin)/nstep);
  
  //Double_t nPart[100]={Npart(0,10),Npart(10,20),Npart(20,30),Npart(30,40),Npart(40,50),Npart(50,60),
  //		       Npart(60,70),Npart(70,80),Npart(80,90),Npart(90,100),Npart(100,120),
  //		       Npart(120,200)};
  
  //Double_t N_Coll[100]={NColl(0,10),NColl(10,20),NColl(20,30),NColl(30,40),NColl(40,50),NColl(50,60),
  //		       NColl(60,70),NColl(70,80),NColl(80,90),NColl(90,100),NColl(100,120),
  //		       NColl(120,200)};

  
 
 const int  NN=5;  
 
 Double_t nPart[5]={ Npart(0,10),  Npart(10,40), Npart(40,80), Npart(80,120), Npart(120,200)};
 Double_t N_Coll[5]={ NColl(0,10),  NColl(10,40), NColl(40,80), NColl(80,120), NColl(120,200)};

 Double_t jpsiFormCMS[NN]={0.0}, jpsiDissCMS[NN]={0.0}, nJpsiFCMS[NN]={0.0}, RAACMS[NN]={0.0}, CNMEffectsCMS[NN]={0.0},PionDissCMS[NN]={0.0},GluonDissCMS[NN]={0.0};

 
 // Loop over centrality
 
 cout<<endl;
 cout<<" CMS calculations Npart      ========================     "<<endl;
 cout<<endl;
 
 for(int i=0; i<NN; i++) {
 //for(int i=0; i<0; i++) {  
   Double_t nColl = N_Coll[i];
   
   Double_t nJpsi = nJpsi0_21*nColl/nColl0;
   
   Double_t nCC = nCC0_21*nColl/nColl0;

   Double_t FacNaprtTempCMS = nPart[i]/Npart(0,10);

   Double_t ssCent = ss05*(grdNDetaNpart->Eval(nPart[i])/grdNDetaNpart->Eval(Npart(0,10)));
   //Double_t ssCent = ss05*nPart[i]/Npart(0,10);
   Double_t R0Cent = R05*TMath::Power(nPart[i]/Npart(0,10),0.5);   
 
   cout<<" calculating the Temp Vs Tau in Cent Loop CMS "<<endl;
   
   NTau=0;
   CalculateTandf_LatticeEOS(ssCent, R0Cent);

   cout<<" calculated the Temp Vs Tau in Cent Loop"<<endl;
   cout<<" NTau Value = "<<NTau<<endl;
   
   //for(int j=0;j<NTau;j++){
   //cout<<"j  "<<j<<"   "<<Tau[j]<<" Temp : "<<TempTau[j]<<" FQGP  "<< fQGP[j]<<endl;
   //}
   
   cout<<" Value of isUps: "<<isUps<<" Remaining steps: "<<(NN-i)<<endl;
   if(isUps == 1|| isUps == 2 || isUps == 3 ){cout<<" save Time =============================== "<<endl;}
   
   

   if(SaveTime ==0){
     cout<<" CMS:Resetting the histograms in Cent loop : "<<endl;
     for(int j=0;j<200;j++)
       {
	 HistRegenJpsiPt[j]->Reset();
       }
     
     if(isUps ==0){
       timer.Start();  
       for(int j=0;j<NTau;j++){
	 FormRateMC_P_New(TempTau[j],j);  
       }
       timer.Stop();
       cout<<" Calculated the Formation histograms in centrality loop : cpu sec ==============: "<<
	 timer.CpuTime()<< "  real sec: "<<timer.RealTime()<<endl;
     }
   }

   
   //for (int l= 0; l<NTau; l++) {
   //cout<<i<<"  "<<NTau<<" trying to acsess Histo Array :  "<<HistRegenJpsiPt[l]->GetName()<<endl;
   //cout<<" bin content  "<<HistRegenJpsiPt[l]->GetBinContent(HistRegenJpsiPt[1]->FindBin(6.25))<<endl;
   //}
   
   //====================================================================================//
   //====================   CMS ========================================================//
   //==================================================================================//



   Double_t intPiondissCMS = 0.0;
   
   if(isUps==0){intPiondissCMS = PionDiss_PtInt(6.25);}
   if(isUps==1){intPiondissCMS = PionDiss_PtInt(0.25);}
   if(isUps==2){intPiondissCMS = PionDiss_PtInt(0.25);}
   if(isUps==3){intPiondissCMS = PionDiss_PtInt(0.25);}
   PionDissCMS[i]=intPiondissCMS;
     
   cout<<" calculated Pion Diss: "<<PionDissCMS[i]<<endl;

   Double_t intCNMEffectCMS = 0.0;
   if(isUps==0){intCNMEffectCMS = CNMVsNPart(nPart[i],grShJPsiPtCutVsNPart_Y24);}
   if(isUps==1){intCNMEffectCMS = CNMVsNPart(nPart[i],grShUpsilonVsNPartMid);}
   if(isUps==2){intCNMEffectCMS = CNMVsNPart(nPart[i],grShUpsilonVsNPartMid);}
   if(isUps==3){intCNMEffectCMS = CNMVsNPart(nPart[i],grShUpsilonVsNPartMid);}
   CNMEffectsCMS[i]=intCNMEffectCMS;

   cout<<" calculated Shadowing: "<< CNMEffectsCMS[i] <<endl;

   Double_t intdissCMS = 0.0;
   if(isUps==0){intdissCMS = IntDiss_PtInt(6.25,nPart[i]);}
   if(isUps==1){intdissCMS = IntDiss_PtInt(0.25,nPart[i]);}
   if(isUps==2){intdissCMS = IntDiss_PtInt(0.25,nPart[i]);}
   if(isUps==3){intdissCMS = IntDiss_PtInt(0.25,nPart[i]);}
   jpsiDissCMS[i] = intdissCMS;
   
   GluonDissCMS[i]= jpsiDissCMS[i]/(PionDissCMS[i]*CNMEffectsCMS[i]);
   
   cout<<" calculated Gluon Dissociation: "<<  GluonDissCMS[i] <<endl;

   cout<<" Total Dissociation : "<<jpsiDissCMS[i]<<endl;

   Double_t intformCMS=0.0;
   
   //intformCMS = RegenratedQuarkoniaNPart(nCC,nJpsi,R0Cent,JPsiIntSh_21,1);
   //RegenratedQuarkoniaNPart(Double_t NCC, Double_t NJPsi, Double_t R0, Double_t IntSh, Int_t PtRange, Double_t NPart);
   
   intformCMS = RegenratedQuarkoniaNPart(nCC,nJpsi,R0Cent,JPsiIntSh_21,1,nPart[i]);
   
   //jpsiFormCMS[i] = intdissCMS*intformCMS*nCC*nCC/(nJpsi);

   jpsiFormCMS[i] = intformCMS;
  
   nJpsiFCMS[i] = (jpsiDissCMS[i] + jpsiFormCMS[i] );
   RAACMS[i] = nJpsiFCMS[i];

   cout<<" calculated Formation: "<<  jpsiFormCMS[i] <<endl;

  }
 
 cout<<endl;
 for(int i=0; i<NN; i++) {
   cout <<"RAA CMS  nPart "<<nPart[i]<<" pion   "<<PionDissCMS[i]<<" CNM  "<<CNMEffectsCMS[i]<<" Gluon  "<< GluonDissCMS[i] <<" Regen "<<jpsiFormCMS[i]<<" Raa  "<<RAACMS[i]<<endl;
   cout<<jpsiDissCMS[i]<<"   "<<jpsiFormCMS[i]<<endl;


 }
 cout<<endl;




  //====================== Calculation graphs for CMS ==========================//
  
  //PionDiss graph
  TGraph *grPionDissCMS_M = new TGraph(NN,nPart,PionDissCMS);
  grPionDissCMS_M->SetName("grPionDissCMS_M");
  grPionDissCMS_M->SetTitle("grPionDissCMS_M");
  grPionDissCMS_M->SetLineWidth(2);
  grPionDissCMS_M->SetLineColor(8);
  grPionDissCMS_M->SetLineStyle(2);
  grPionDissCMS_M->GetXaxis()->SetTitle("N_{Part}");
  grPionDissCMS_M->GetYaxis()->SetTitle("PionDiss");
  
  //CNM graph
  TGraph *grCNMCMS_M = new TGraph(NN,nPart,CNMEffectsCMS);
  grCNMCMS_M->SetName("grCNMCMS_M");
  grCNMCMS_M->SetTitle("grCNMCMS_M");
  grCNMCMS_M->SetLineWidth(2);
  grCNMCMS_M->SetLineColor(2);
  grCNMCMS_M->SetLineStyle(2);
  grCNMCMS_M->GetXaxis()->SetTitle("N_{Part}");
  grCNMCMS_M->GetYaxis()->SetTitle("CNM");
  
  
  TGraph *grDissCMS_M = new TGraph(NN,nPart,GluonDissCMS);
  grDissCMS_M->SetName("grDissCMS_M");
  grDissCMS_M->SetTitle("grDissCMS_M");
  grDissCMS_M->SetLineWidth(2);
  grDissCMS_M->SetLineColor(4);
  grDissCMS_M->SetLineStyle(6);
  
  
  TGraph *grFormCMS_M = new TGraph(NN,nPart,jpsiFormCMS);
  grFormCMS_M->SetName("grFormCMS_M");
  grFormCMS_M->SetTitle("grFormCMS_M");
  grFormCMS_M->SetLineWidth(2);
  grFormCMS_M->SetLineColor(6);
  grFormCMS_M->SetLineStyle(4);
  
  
  // RAA graphs
  TGraph *grRAACMS_M = new TGraph(NN,nPart,RAACMS);
  grRAACMS_M->SetName("grRAACMS_M");
  grRAACMS_M->SetTitle("grRAACMS_M");
  grRAACMS_M->SetLineWidth(2);
  grRAACMS_M->GetXaxis()->SetTitle("N_{Part}");
  grRAACMS_M->GetYaxis()->SetTitle("R_{AA}");
  
  TLegend *legd7 = new TLegend(0.17,0.66,0.63,0.93);
  legd7->SetBorderSize(0);
  legd7->SetFillStyle(0);
  legd7->SetFillColor(0);
  legd7->SetTextSize(0.040);
  
  TGraphErrors *CMS_NPart_GhostGraph = grf_GhostGraph(21, 2, 2);
  
  if(QQbarVar==1){legd7->AddEntry(CMS_NPart_GhostGraph,"CMS Data, |y|<2.4, 6.5<p_{T}<30","P");}
  
  if(QQbarVar==4){legd7->AddEntry(CMS_NPart_GhostGraph,"CMS Data, |y|^{#varUpsilon(1S)}<2.4, p_{T}<30","P");}
  if(QQbarVar==5 || QQbarVar==6 || QQbarVar==7){legd7->AddEntry(CMS_NPart_GhostGraph,"CMS Data, |y|^{#varUpsilon(2S)}<2.4, p_{T}<30","P");}

  legd7->AddEntry(grDissCMS_M,"Gluon Dissociation","L");
  legd7->AddEntry(grPionDissCMS_M,"Comover","L");
  legd7->AddEntry(grCNMCMS_M,"CNM Effects","L");
  legd7->AddEntry(grFormCMS_M,"Regeneration","L");
  legd7->AddEntry(grRAACMS_M,"Total (R_{AA})","L");
  
    

  TLegend *legd7_NoData = new TLegend(0.17,0.66,0.63,0.93);
  legd7_NoData->SetBorderSize(0);
  legd7_NoData->SetFillStyle(0);
  legd7_NoData->SetFillColor(0);
  legd7_NoData->SetTextSize(0.040);
  
  legd7_NoData->AddEntry(grDissCMS_M,"Gluon Dissociation","L");
  legd7_NoData->AddEntry(grPionDissCMS_M,"Comover","L");
  legd7_NoData->AddEntry(grCNMCMS_M,"CNM Effects","L");
  legd7_NoData->AddEntry(grFormCMS_M,"Regeneration","L");
  legd7_NoData->AddEntry(grRAACMS_M,"Total (R_{AA})","L");


  new TCanvas;
  gPad->SetTicks();
  if(QQbarVar==1)Draw_CMS_JPsi_RaaVsNpart();
  
  if(QQbarVar==4)Draw_CMS_Y1S_5TeV_RaaVsNpart();
  if(QQbarVar==5)Draw_CMS_Y2S_5TeV_RaaVsNpart();
  if(QQbarVar==6)Draw_CMS_Y2S_5TeV_RaaVsNpart();
  if(QQbarVar==7)Draw_CMS_Y2S_5TeV_RaaVsNpart();
  
  grCNMCMS_M->Draw("Lsame");
  grPionDissCMS_M->Draw("sameL");
  grFormCMS_M->Draw("sameL");
  grDissCMS_M->Draw("sameL");
  grRAACMS_M->Draw("sameL");
  legd7->Draw("same");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_CMS.pdf");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_CMS.png");
    
  
  new TCanvas;
  gPad->SetTicks();
  grCNMCMS_M->GetYaxis()->SetRangeUser(0.0,2.0);
  grCNMCMS_M->Draw("AL");
  grPionDissCMS_M->Draw("sameL");
  grFormCMS_M->Draw("sameL");
  grDissCMS_M->Draw("sameL");
  grRAACMS_M->Draw("sameL");
  legd7_NoData->Draw("same");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_CMS_NoData.pdf");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_CMS_NoData.png");
  
  cout<<endl<<endl<<endl;


  Double_t jpsiFormATLAS[NN]={0.0}, jpsiDissATLAS[NN]={0.0}, nJpsiFATLAS[NN]={0.0}; 
  Double_t RAAATLAS[NN]={0.0}, CNMEffectsATLAS[NN]={0.0},PionDissATLAS[NN]={0.0}, GluonDissATLAS[NN]={0.0};
        

  // Loop over centrality
  

  for(int i=0; i<NN; i++) {
  //for(int i=0; i<0; i++) {
    
    if(QQbarVar==1){

      
      Double_t nColl = N_Coll[i];
      Double_t nJpsi = nJpsi0_21*nColl/nColl0;
      Double_t nCC = nCC0_21*nColl/nColl0;
      
      Double_t FacNaprtTempATLAS = nPart[i]/Npart(0,10);
      
      Double_t ssCent = ss05*(grdNDetaNpart->Eval(nPart[i])/grdNDetaNpart->Eval(Npart(0,10)));
      //Double_t ssCent = ss05*nPart[i]/Npart(0,10);
      Double_t R0Cent = R05*TMath::Power(nPart[i]/Npart(0,10),0.5);   
      
      cout<<" calculating the Temp Vs Tau in Cent Loop ATLAS "<<endl;
      NTau=0;
      CalculateTandf_LatticeEOS(ssCent, R0Cent);
      cout<<" calculated the Temp Vs Tau in Cent Loop"<<endl;
      cout<<" NTau Value = "<<NTau<<endl;
      
      //for(int j=0;j<=NTau;j++){
      //cout<<" Temp : "<<TempTau[j]<<" FQGP  "<< fQGP[j]<<endl;
      //}
      
      cout<<" Value of isUps: "<<isUps<<" Remaining steps: "<<(NN-i)<<endl;
      if(isUps == 1|| isUps == 2 || isUps == 3 ){cout<<" save Time =============================== "<<endl;}
      
      
      //====================================================================================//
      //====================   ATLAS ========================================================//
      //==================================================================================//
      Double_t intPiondissATLAS = 0.0;
      
      if(isUps==0){intPiondissATLAS = PionDiss_PtInt(9.25);}
      if(isUps==1){intPiondissATLAS = PionDiss_PtInt(0.25);}
      if(isUps==2){intPiondissATLAS = PionDiss_PtInt(0.25);}
      PionDissATLAS[i]=intPiondissATLAS;
      
      cout<<" calculated Pion Diss: "<<PionDissATLAS[i]<<endl;
      
      Double_t intCNMEffectATLAS = 0.0;
      if(isUps==0){intCNMEffectATLAS = CNMVsNPart(nPart[i],grShJPsiPtCutVsNPart_Y20);}
      if(isUps==1){intCNMEffectATLAS = CNMVsNPart(nPart[i],grShUpsilonVsNPartMid);}
      if(isUps==2){intCNMEffectATLAS = CNMVsNPart(nPart[i],grShUpsilonVsNPartMid);}
      CNMEffectsATLAS[i]=intCNMEffectATLAS;
      
      cout<<" calculated Shadowing: "<< CNMEffectsATLAS[i] <<endl;
      
      Double_t intdissATLAS = 0.0;
      if(isUps==0){intdissATLAS = IntDiss_PtInt(9.25,nPart[i]);}
      
      if(isUps==1){intdissATLAS =IntDiss_PtInt(0.25,nPart[i]);}
      if(isUps==2){intdissATLAS =IntDiss_PtInt(0.25,nPart[i]);}
      jpsiDissATLAS[i] = intdissATLAS;
      
      GluonDissATLAS[i]= jpsiDissATLAS[i]/(PionDissATLAS[i]*CNMEffectsATLAS[i]);
      
      cout<<" calculated Gluon Dissociation: "<<  jpsiDissATLAS[i] <<endl;
      
      
      Double_t intformATLAS=0.0;
      intformATLAS = RegenratedQuarkoniaNPart(nCC,nJpsi,R0Cent,JPsiIntSh_21,2,nPart[i]);
      jpsiFormATLAS[i] = intformATLAS;
      nJpsiFATLAS[i] = (jpsiDissATLAS[i] + jpsiFormATLAS[i] );
      RAAATLAS[i] = nJpsiFATLAS[i];
      
      cout<<" calculated Formation: "<<  jpsiFormATLAS[i] <<endl;
      
      cout<<"RAA ATLAS "<<endl;
      cout <<"nPart "<<nPart[i]<<" pion   "<<PionDissATLAS[i]<<" cnm  "<<CNMEffectsATLAS[i]<<" raa  "<<RAAATLAS[i]<<" gluon "<< RAAATLAS[i]/(PionDissATLAS[i]*CNMEffectsATLAS[i])<<endl;
    }
    
  }


  //====================== NPart Calculation graphs for ATLAS ==========================//
  
  //PionDiss graph
  TGraph *grPionDissATLAS_M = new TGraph(NN,nPart,PionDissATLAS);
  grPionDissATLAS_M->SetName("grPionDissATLAS_M");
  grPionDissATLAS_M->SetTitle("grPionDissATLAS_M");
  grPionDissATLAS_M->SetLineWidth(2);
  grPionDissATLAS_M->SetLineColor(8);
  grPionDissATLAS_M->SetLineStyle(2);
  grPionDissATLAS_M->GetXaxis()->SetTitle("N_{Part}");
  grPionDissATLAS_M->GetYaxis()->SetTitle("PionDiss");
  
  //CNM graph
  TGraph *grCNMATLAS_M = new TGraph(NN,nPart,CNMEffectsATLAS);
  grCNMATLAS_M->SetName("grCNMATLAS_M");
  grCNMATLAS_M->SetTitle("grCNMATLAS_M");
  grCNMATLAS_M->SetLineWidth(2);
  grCNMATLAS_M->SetLineColor(2);
  grCNMATLAS_M->SetLineStyle(2);
  grCNMATLAS_M->GetXaxis()->SetTitle("N_{Part}");
  grCNMATLAS_M->GetYaxis()->SetTitle("CNM");
  
  
  TGraph *grDissATLAS_M = new TGraph(NN,nPart,GluonDissATLAS);
  grDissATLAS_M->SetName("grDissATLAS_M");
  grDissATLAS_M->SetTitle("grDissATLAS_M");
  grDissATLAS_M->SetLineWidth(2);
  grDissATLAS_M->SetLineColor(4);
  grDissATLAS_M->SetLineStyle(6);
  
  
  TGraph *grFormATLAS_M = new TGraph(NN,nPart,jpsiFormATLAS);
  grFormATLAS_M->SetName("grFormATLAS_M");
  grFormATLAS_M->SetTitle("grFormATLAS_M");
  grFormATLAS_M->SetLineWidth(2);
  grFormATLAS_M->SetLineColor(6);
  grFormATLAS_M->SetLineStyle(4);
  
  
  // RAA graphs
  TGraph *grRAAATLAS_M = new TGraph(NN,nPart,RAAATLAS);
  grRAAATLAS_M->SetName("grRAAATLAS_M");
  grRAAATLAS_M->SetTitle("grRAAATLAS_M");
  grRAAATLAS_M->SetLineWidth(2);
  grRAAATLAS_M->GetXaxis()->SetTitle("N_{Part}");
  grRAAATLAS_M->GetYaxis()->SetTitle("R_{AA}");
  

  TLegend *legd8 = new TLegend(0.17,0.66,0.63,0.93);
  legd8->SetBorderSize(0);
  legd8->SetFillStyle(0);
  legd8->SetFillColor(0);
  legd8->SetTextSize(0.040);
  

  TGraphErrors *ATLAS_NPart_GhostGraph = grf_GhostGraph(24, 4, 4);

  legd8->AddEntry(ATLAS_NPart_GhostGraph,"ATLAS Data, |y|<2.0, 9<p_{T}<40","P");
  legd8->AddEntry(grDissATLAS_M,"Gluon Dissociation","L");
  legd8->AddEntry(grPionDissATLAS_M,"Comover","L");
  legd8->AddEntry(grCNMATLAS_M,"CNM Effects","L");
  legd8->AddEntry(grRAAATLAS_M,"Total (R_{AA})","L");
  
    
  new TCanvas;
  gPad->SetTicks();
  Draw_ATLAS_JPsi_RaaVsNpart();
  grCNMATLAS_M->Draw("Lsame");
  grPionDissATLAS_M->Draw("sameL");
  grFormATLAS_M->Draw("sameL");
  grDissATLAS_M->Draw("sameL");
  grRAAATLAS_M->Draw("sameL");
  legd8->Draw("same");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_ATLAS.pdf");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_ATLAS.png");
    
    

  TLegend *legd8_NoData = new TLegend(0.16,0.66,0.86,0.93);
  legd8_NoData->SetBorderSize(0);
  legd8_NoData->SetFillStyle(0);
  legd8_NoData->SetFillColor(0);
  legd8_NoData->SetTextSize(0.040);
  

  legd8_NoData->AddEntry(ATLAS_NPart_GhostGraph,"ATLAS Data, |y| < 2.0, 9 < p_{T} < 40","P");
  legd8_NoData->AddEntry(grDissATLAS_M,"Gluon Dissociation","L");
  legd8_NoData->AddEntry(grPionDissATLAS_M,"Comover","L");
  legd8_NoData->AddEntry(grCNMATLAS_M,"CNM Effects","L");
  legd8_NoData->AddEntry(grRAAATLAS_M,"Total (R_{AA})","L");


  new TCanvas;
  gPad->SetTicks();
  grCNMATLAS_M->GetYaxis()->SetRangeUser(0.0,2.0);
  grCNMATLAS_M->Draw("AL");
  grPionDissATLAS_M->Draw("sameL");
  grFormATLAS_M->Draw("sameL");
  grDissATLAS_M->Draw("sameL");
  grRAAATLAS_M->Draw("sameL");
  legd8_NoData->Draw("same");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_ATLAS_NoData.pdf");
  gPad->SaveAs("Figures/Fig_RAA_Npart_HighPt_ATLAS_NoData.png");
    

  cout<<endl<<endl<<endl;







  Double_t jpsiFormALICE[NN]={0.0}, jpsiDissALICE[NN]={0.0}, nJpsiFALICE[NN]={0.0}, RAAALICE[NN]={0.0};
  Double_t CNMEffectsALICE[NN]={0.0},PionDissALICE[NN]={0.0},GluonDissALICE[NN]={0.0};

  //============== centrality calculations ALICE

  
  for(int i=0; i<NN; i++) {
  //for(int i=0; i<0; i++) {

   Double_t nColl = N_Coll[i];
   Double_t nJpsi = nJpsi0_21*nColl/NColl(0,10);
   Double_t nCC = nCC0_21*nColl/NColl(0,10);

   Double_t FacNaprtTempALICE = nPart[i]/Npart(0,10);   
   Double_t ssCent = ss05*(grdNDetaNpart->Eval(nPart[i])/grdNDetaNpart->Eval(Npart(0,10)));
   //Double_t ssCent = ss05*nPart[i]/Npart(0,10);
   Double_t R0Cent = R05*TMath::Power(nPart[i]/Npart(0,10),0.5);   
   
   cout<<"NCOLL: "<<nColl<<"  "<<nJpsi<<"  "<<nCC<<"  "<<ssCent<<"  "<<R0Cent<<endl;
   
   cout<<" calculating the Temp Vs Tau in Cent Loop ALICE "<<endl;
   NTau=0;
   CalculateTandf_LatticeEOS(ssCent, R0Cent);
   cout<<" calculated the Temp Vs Tau in Cent Loop"<<endl;
   cout<<" NTau Value = "<<NTau<<endl;
   //for(int j=0;j<=NTau;j++){
   //cout<<" Temp : "<<TempTau[j]<<" FQGP  "<< fQGP[j]<<endl;
   //}
   
   cout<<" Value of isUps: "<<isUps<<" Remaining steps: "<<(NN-i)<<endl;
   if(isUps == 1|| isUps == 2 || isUps == 3 ){cout<<" save Time =============================== "<<endl;}
    
    
    if(SaveTime ==0){
      cout<<" Resetting the histograms in Cent loop : ALICE"<<endl;
      for(int j=0;j<200;j++)
	{
	  HistRegenJpsiPt[j]->Reset();
	}
      
      if(isUps ==0){
	timer.Start();  
	for(int j=0;j<NTau;j++){
	  FormRateMC_P_New(TempTau[j],j);  
	}
	timer.Stop();
	cout<<" Calculated the Formation part in centrality loop : cpu sec ==============: "<<timer.CpuTime()<< "  real sec: "<<timer.RealTime()<<endl;
      }
    }
    
    ////   ALICE  calculations with centrality///////////
    Double_t intdissALICE=0; 
    Double_t intformALICE=0;
    
    //if(isUps==0){
      
    if(QQbarVar==1) CNMEffectsALICE[i] = CNMVsNPart(nPart[i],grShJPsiVsNPart_Y2440);
    
    if(QQbarVar==4) CNMEffectsALICE[i] = CNMVsNPart(nPart[i],grShUpsilonVsNPartFor);
    if(QQbarVar==5) CNMEffectsALICE[i] = CNMVsNPart(nPart[i],grShUpsilonVsNPartFor);
    if(QQbarVar==6) CNMEffectsALICE[i] = CNMVsNPart(nPart[i],grShUpsilonVsNPartFor);
    if(QQbarVar==7) CNMEffectsALICE[i] = CNMVsNPart(nPart[i],grShUpsilonVsNPartFor);
    cout<<" calculated CNM low pT: "<<CNMEffectsALICE[i]<<endl;
    

    PionDissALICE[i] = PionDiss_PtInt(0.25);
    cout<<" calculated Pion Diss low pT: "<<PionDissALICE[i]<<endl;
      
    
    if(QQbarVar==1) intdissALICE = IntDiss_PtInt(0.25,nPart[i]);
    if(QQbarVar==4 || QQbarVar==5 || QQbarVar==6 || QQbarVar==7)intdissALICE = IntDiss_PtInt(0.25,nPart[i]);

    GluonDissALICE[i]=intdissALICE/(CNMEffectsALICE[i]*PionDissALICE[i]);
    
    cout<<" calculated total Diss low pT: "<<intdissALICE<<endl;
    
    if(SaveTime==0){
      cout<<" Running Big loop "<<endl;
      timer.Start();
      //intformALICE = RegenratedQuarkoniaNPart(nCC,nJpsi,R0Cent,JPsiIntSh_21,0);
      intformALICE = RegenratedQuarkoniaNPart(nCC, nJpsi, R0Cent, JPsiIntSh_21, 0, nPart[i]);
      timer.Stop();
      cout<<" Regenrating Quarkonia ALICE: "<<intformALICE<<" in cpu sec: "<<timer.CpuTime()<<endl;
    }
    //}
    
    jpsiDissALICE[i] = intdissALICE;
    jpsiFormALICE[i] = intformALICE;
    
    //RAA ALICE
    nJpsiFALICE[i] = (jpsiDissALICE[i] + jpsiFormALICE[i] );
    RAAALICE[i] = nJpsiFALICE[i];
    
    
    
  }
  
  
  cout<<endl<<endl;
  cout<<" NPart "<<"    "<<" CNM "<<"    "<<" Pion "<<"    "<<" Total Diss "<<"    "<<"  Formation  "<<"  RAA "<<endl;
  for(int i=0; i<NN; i++) {
    cout << nPart[i]<<"  "<<CNMEffectsALICE[i]<<"  "<<PionDissALICE[i]<<"   "<<jpsiDissALICE[i]<<"   "<<jpsiFormALICE[i]<<"  "<<RAAALICE[i]<<endl;
  }
  cout<<endl<<endl;







  //====================== Calculation graphs for ALICE ==========================//
  //PionDiss graph
  TGraph *grPionDissALICE_M = new TGraph(NN,nPart,PionDissALICE);
  grPionDissALICE_M->SetName("grPionDissALICE_M");    
  grPionDissALICE_M->SetTitle("grPionDissALICE_M");    
  grPionDissALICE_M->SetLineWidth(2);
  grPionDissALICE_M->SetLineColor(8);
  grPionDissALICE_M->SetLineStyle(2);
  grPionDissALICE_M->GetXaxis()->SetTitle("N_{Part}");
  grPionDissALICE_M->GetYaxis()->SetTitle("PionDiss");
  
  //CNM graph
  TGraph *grCNMALICE_M = new TGraph(NN,nPart,CNMEffectsALICE);
  grCNMALICE_M->SetName("grCNMALICE_M");
  grCNMALICE_M->SetTitle("grCNMALICE_M");
  grCNMALICE_M->SetLineWidth(2);
  grCNMALICE_M->SetLineColor(2);
  grCNMALICE_M->SetLineStyle(2);
  grCNMALICE_M->GetXaxis()->SetTitle("N_{Part}");
  grCNMALICE_M->GetYaxis()->SetTitle("CNM");
  
  
  
  TGraph *grDissALICE_M = new TGraph(NN,nPart,GluonDissALICE);
  grDissALICE_M->SetName("grDissALICE_M");
  grDissALICE_M->SetTitle("grDissALICE_M");
  grDissALICE_M->SetLineWidth(2);
  grDissALICE_M->SetLineColor(4);
  grDissALICE_M->SetLineStyle(6);
  
  
  
  TGraph *grFormALICE_M = new TGraph(NN,nPart,jpsiFormALICE);
  grFormALICE_M->SetName("grFormALICE_M");
  grFormALICE_M->SetTitle("grFormALICE_M");
  grFormALICE_M->SetLineWidth(2);
  grFormALICE_M->SetLineColor(2);
  grFormALICE_M->SetLineStyle(4);
  
  
  // RAA graphs
  TGraph *grRAAALICE_M = new TGraph(NN,nPart,RAAALICE);
  grRAAALICE_M->SetName("grRAAALICE_M");
  grRAAALICE_M->SetTitle("grRAAALICE_M");
  grRAAALICE_M->SetLineWidth(2);
  grRAAALICE_M->GetXaxis()->SetTitle("N_{Part}");
  grRAAALICE_M->GetYaxis()->SetTitle("R_{AA} ");
  
    
  TLegend *legd9 = new TLegend(0.17,0.56,0.54,0.91);
  legd9->SetBorderSize(0);
  legd9->SetFillStyle(0);
  legd9->SetFillColor(0);
  legd9->SetTextSize(0.035);
  

  TGraphErrors *ALICE_NPart_GhostGraph = grf_GhostGraph(20, 6, 6);
  if(QQbarVar==4 || QQbarVar== 5 || QQbarVar== 6 || QQbarVar== 7)
    {
      ALICE_NPart_GhostGraph->SetMarkerColor(2);
      ALICE_NPart_GhostGraph->SetMarkerStyle(21);
      ALICE_NPart_GhostGraph->SetLineColor(2);
    }


  if(QQbarVar==1)legd9->AddEntry(ALICE_NPart_GhostGraph,"ALICE Data, 2.5 < y^{J/#psi} < 4.0, 0.3 < p_{T} < 8","P");
  
  if(QQbarVar==4 || QQbarVar== 5 || QQbarVar== 6 || QQbarVar== 7){
    legd9->AddEntry(ALICE_NPart_GhostGraph,"ALICE Data, 2.5 < y^{#varUpsilon} < 4.0, p_{T} < 12","P");
  }
  
  legd9->AddEntry(grDissALICE_M,"Gluon Dissociation","L");
  legd9->AddEntry(grCNMALICE_M,"CNM Effects","L");
  legd9->AddEntry(grPionDissALICE_M,"Comover","L");
  legd9->AddEntry(grFormALICE_M,"Formation","L");
  legd9->AddEntry(grRAAALICE_M,"Total","L");
  
  
  new TCanvas;
  gPad->SetTicks();
  
  if(QQbarVar==1) Draw_ALICE_JPsi_RaaVsNpart();
  if(QQbarVar==4 || QQbarVar== 5 || QQbarVar== 6 || QQbarVar== 7 )Draw_ALICE_Y1S_5TeV_RaaVsNpart();
  
  grFormALICE_M->Draw("sameL");
  grDissALICE_M->Draw("sameL");
  grRAAALICE_M->Draw("sameL");
  grCNMALICE_M->Draw("sameL"); 
  grPionDissALICE_M->Draw("sameL"); 
  legd9->Draw("same");
  gPad->SaveAs("Figures/Fig_RAA_Npart_ALICE.pdf");
  gPad->SaveAs("Figures/Fig_RAA_Npart_ALICE.png");
  

  TLegend *legd9_NoData = new TLegend(0.17,0.56,0.54,0.91);
  legd9_NoData->SetBorderSize(0);
  legd9_NoData->SetFillStyle(0);
  legd9_NoData->SetFillColor(0);
  legd9_NoData->SetTextSize(0.035);
  
  legd9_NoData->AddEntry(grDissALICE_M,"Gluon Dissociation","L");
  legd9_NoData->AddEntry(grCNMALICE_M,"CNM Effects","L");
  legd9_NoData->AddEntry(grPionDissALICE_M,"Comover","L");
  legd9_NoData->AddEntry(grFormALICE_M,"Formation","L");
  legd9_NoData->AddEntry(grRAAALICE_M,"Total","L");

  
  new TCanvas;
  grFormALICE_M->GetYaxis()->SetRangeUser(0.0,2.0);
  grFormALICE_M->Draw("AL");
  grDissALICE_M->Draw("sameL");
  grRAAALICE_M->Draw("sameL");
  grCNMALICE_M->Draw("sameL"); 
  grPionDissALICE_M->Draw("sameL");
  legd9_NoData->Draw("same"); 
  gPad->SaveAs("Figures/Fig_RAA_Npart_ALICE_NoData.pdf");
  gPad->SaveAs("Figures/Fig_RAA_Npart_ALICE_NoData.png");
  
  /*
  Char_t OutFileName[100];
  if(QQbarVar==1){sprintf(OutFileName,"JPsiCalculations.root");} 
  if(QQbarVar==4){sprintf(OutFileName,"Y1SCalculations.root");}
  if(QQbarVar==5){sprintf(OutFileName,"Y2SCalculations.root");}
  if(QQbarVar==6){sprintf(OutFileName,"Y3SCalculations.root");}
  if(QQbarVar==7){sprintf(OutFileName,"ChiBCalculations.root");}
  
  TFile *OutFile =new TFile(OutFileName,"RECREATE");
  */

  grSigDQ0->Write();
  grSigDQ0Mod->Write();
 
  grTempVsNpart->Write();
  
  grTempVsTauLatt->Write();
  grFQGPVsTauLatt->Write();
  grTempVsTauLatt_020->Write();
  grFQGPVsTauLatt_020->Write();
  
  DissRate_Khar->Write();
  DissRate_Khar36->Write();
  DissRate_Khar86->Write();
  
  grf_FCharm_Temp->Write();
  grf_FCharm_Thermal_Temp->Write();
  
  grDissRateVsPt1->Write();
  grDissRateVsPt2->Write();
  grDissRateVsPt3->Write();
  
  grf_FCharm_Pt->Write();
  grf_FCharm_Thermal_Pt->Write();
  
  grRAACMSPt_M->Write();
  grRAACMSPt_ELoss_M->Write();
  grGluonDissCMSPt_M->Write();
  grPionDissCMSPt_M->Write();
  grShadowDissCMSPt_M->Write();
  grELossCMSPt_M->Write();
  grRegenCMSPt_M->Write();
  
  if(QQbarVar==1){
    grRAAATLASPt_M->Write();
    grRAAATLASPt_ELoss_M->Write();
    grGluonDissATLASPt_M->Write();
    grPionDissATLASPt_M->Write();
    grShadowDissATLASPt_M->Write();
    grELossATLASPt_M->Write();
    grRegenATLASPt_M->Write();
  }
  
  grRAAALICEPt_M->Write();
  grRAAALICEPt_ELoss_M->Write();
  grPionDissALICEPt_M->Write();
  grShadowDissALICEPt_M->Write();
  grGluonDissALICEPt_M->Write();
  grELossALICEPt_M->Write();
  grRegenALICEPt_M->Write();
  
  
  grPionDissCMS_M->Write();
  grCNMCMS_M->Write();
  grDissCMS_M->Write();
  grFormCMS_M->Write();
  grRAACMS_M->Write();
  
  grPionDissATLAS_M->Write();
  grCNMATLAS_M->Write();
  grDissATLAS_M->Write();
  grFormATLAS_M->Write();
  grRAAATLAS_M->Write();
  
  
  grPionDissALICE_M->Write();
  grCNMALICE_M->Write();
  grDissALICE_M->Write();
  grFormALICE_M->Write();
  grRAAALICE_M->Write();
  
  
  
  OutFile->Write();
  OutFile->Close();
  
  
  TTimeStamp ts2;
  cout<<endl<<endl;
  cout << ts1.AsString() << "    "<<ts2.AsString() <<endl; 
  cout<<endl<<endl;
  

}





//===============================================================================================//
//==================== JPsi pion dissociation ===================================================//
//===============================================================================================//

Double_t PionDiss_PtInt(Double_t PtMin)
{
  //for 500 bins and step = 0.1
  //PtMin=0.05;
  //PtMin=6.45;
  
  //for 100 bins and step = 0.5
  //PtMin=0.25;
  //PtMin=6.25;

  Double_t Ptmax=0;
  if(isUps==0)Ptmax=20.75;
  if(isUps==1 || isUps ==2 || isUps ==3)Ptmax=40.25;

  Double_t Ptstep=0.5;
  int NN_Pt= (int)((Ptmax-PtMin)/Ptstep);
  Double_t sum=0;
  Double_t sumPt=0;
  Double_t yield =0;
  for(int i=0;i<=NN_Pt;i++) {
    Double_t Pt=PtMin+i*Ptstep;
    if(isUps==0){yield=Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(Pt));}
    if(isUps==1){yield=Y1S_Pt->GetBinContent(Y1S_Pt->FindBin(Pt));}
    if(isUps==2){yield=Y2S_Pt->GetBinContent(Y2S_Pt->FindBin(Pt));}
    if(isUps==3){yield=Y2S_Pt->GetBinContent(Y2S_Pt->FindBin(Pt));}
    sum= sum + (yield*PionDiss_All(Pt));
    sumPt= sumPt + yield;  
  }
  return sum/sumPt;
}




Double_t PionDiss_All(Double_t PtMin)
{
  Double_t sumd = 0.0;
  for (int i= 0; i<NTau; i++) {
    sumd = sumd + ( (1-fQGP[i])*PionJPsiDiss(PtMin,TempTau[i]));
  }
  return exp(-sumd*stepTau/hbarc);
}


Double_t PionJPsiDiss(Double_t P1, Double_t T)
{
  double m1Square =  mJpsi*mJpsi;
  double m2Square =  mPi*mPi;
  Double_t E1 = TMath::Sqrt( P1*P1 + m1Square);

  Double_t P2Step = 0.1;
  Double_t P2Start = 0.00001;
  Double_t P2End = 10.0;
  int P2N =(int)((P2End - P2Start)/P2Step);
  
  Double_t CosThetaStep = 0.1;
  Double_t CosThetaStart = 1.0;
  Double_t CosThetaEnd = -1.0;
  int CosThetaN =(int)((CosThetaStart - CosThetaEnd)/CosThetaStep);
  
  Double_t SumP2=0;  
  
  Double_t ThetaVal=4.0*mD*mD;
      
  for(int i =0; i <= P2N ; i++) {
    Double_t P2 = P2Start + i*P2Step;
    Double_t E2 = TMath::Sqrt(P2*P2 + m2Square);
    double fPion = fpion(P2,T);

    Double_t SumCosTheta=0;
    for(int k =0; k <= CosThetaN ; k++)  {
      Double_t CosTheta = CosThetaStart - k*CosThetaStep;
      Double_t fourP1DotP2 = E1*E2-P1*P2*CosTheta;
      Double_t S = m1Square + m2Square + 2*fourP1DotP2;
      if(S <= ThetaVal ) continue;
      Double_t VRel =TMath::Sqrt(fourP1DotP2*fourP1DotP2 - m1Square*m2Square)/(E1*E2);
      
      //Pion - Quarkonia Cross Section

      
      //Double_t SGD = 0.1/hbarc2;
      //Double_t SGD = SigmaPionDiss/hbarc2;
      //Quarkonia-Pion Cross Section
      
      //sigma in (GeV)^{-2}
      Double_t SGD = SigmaPionDissFunc(S)/(10.0*hbarc2);
    
      SumCosTheta = SumCosTheta + SGD*VRel;


    }//theta
    
    double IntCosTheta = SumCosTheta*CosThetaStep;
    SumP2 = SumP2 + P2*P2*fPion*IntCosTheta;  
  }//P2

  Double_t IntP2 = (SumP2*P2Step)/(4.0*pi2);
  return IntP2;
}

Double_t RhoPi(Double_t T)
{
 
  Double_t P2Step = 0.1;
  Double_t P2Start = 0.00001;
  Double_t P2End = 10.0;
  int P2N =(int)((P2End - P2Start)/P2Step);
  Double_t SumP2=0;  
  
  for(int i =0; i <= P2N ; i++) {
    Double_t P2 = P2Start + i*P2Step;
    SumP2 = SumP2 + P2*P2*fpion(P2,T);
  }
  Double_t IntP2 = SumP2*P2Step/(2.0*pi2);
  return IntP2;
}




Double_t fpion(Double_t pPi, Double_t T)
{
  Double_t gi=3.0;
  Double_t fg = gi / (TMath::Exp( sqrt(pPi*pPi + mPi*mPi) / T ) - 1 );
  //Double_t fg = gi/TMath::Exp(sqrt(pPi*pPi + mPi*mPi)/T);
  
  return fg;
}


Double_t SigmaPionDissFunc(Double_t S)
{
  Double_t RootS = TMath::Sqrt(S);
  Double_t Sigma =0;  

  if(QQbarVar==1){
    Sigma = grfJPsiPionCross->Eval(RootS);
  }
  
  if(QQbarVar==4){
    Sigma = grfUpsilonPionCross->Eval(RootS);
  }

 if(QQbarVar==5){
    Sigma = grfUpsilon2SPionCross->Eval(RootS);
  }

 if(QQbarVar==6){
   Sigma = grfUpsilon2SPionCross->Eval(RootS);
 }

 if(QQbarVar==7){
   Sigma = grfUpsilon2SPionCross->Eval(RootS);
 }
  
  return Sigma; // it return as mb

}






// Dissociation rates ///////////////////////////////////////////////
//Double_t NPartVsNColl(Double_t ValNPart)
Double_t IntDiss_NPartInt(Double_t CentMin, Double_t CentMax, Double_t PtRange, Double_t Pt)
{
  
  Double_t CentStep=10;
  
  //Double_t NPartMax = Npart(CentMin*2,CentMin*2+1); 
  //Double_t NPartMin = Npart(CentMax*2-1,CentMax*2); 
  
  Int_t NN_Cent= (Int_t)((CentMax-CentMin)/CentStep);
  
  Double_t sum=0;
  Double_t sumPt=0;
  Double_t yy = 0.0;
  
  Double_t yield =0; // This is NColl*N_{Q}^{pp}(pT)
                     //N_{Q}^{pp}(pT) --> comes from JPsi pT histogram but cancel in this case

  Double_t ShadowingFactor=0;
  Double_t PionDissFactor=0;
  Double_t GluonDissFactor=0;

  //if(QQbarVar ==1){ShadowDissATLASPt[i]=gr_shadowing_2424->Eval(PtATLAS[i]);}

  if(PtRange==1 && isUps ==0){ yy = gr_shadowing_2540->Eval(Pt);} //ALICE
  if(PtRange==2 && isUps ==0){ yy = gr_shadowing_2424->Eval(Pt);} //CMS
  if(PtRange==3 && isUps ==0){ yy = gr_shadowing_2424->Eval(Pt);} //ATLAS 
  
  if(PtRange==1 && isUps ==1){ yy = grf_shadowing_502tev_25_y_40_Y1S->Eval(Pt);}
  if(PtRange==1 && isUps ==2){ yy = grf_shadowing_502tev_25_y_40_Y2S->Eval(Pt);}
  if(PtRange==1 && isUps ==3){ yy = grf_shadowing_502tev_25_y_40_Y2S->Eval(Pt);}


  if(PtRange==2 && isUps ==1){ yy = grf_shadowing_502tev_24_y_24_Y1S->Eval(Pt);}
  if(PtRange==2 && isUps ==2){ yy = grf_shadowing_502tev_24_y_24_Y2S->Eval(Pt);}
  if(PtRange==2 && isUps ==3){ yy = grf_shadowing_502tev_24_y_24_Y2S->Eval(Pt);}


  if(PtRange==3 && isUps ==1){ yy = grf_shadowing_502tev_24_y_24_Y1S->Eval(Pt);}
  if(PtRange==3 && isUps ==2){ yy = grf_shadowing_502tev_24_y_24_Y2S->Eval(Pt);}
  if(PtRange==3 && isUps ==3){ yy = grf_shadowing_502tev_24_y_24_Y2S->Eval(Pt);}


  for(Int_t i=0;i<=NN_Cent;i++) {
    
    Int_t Cent1 = CentMin + i*CentStep;
    Int_t Cent2 = CentMin + (i+1)*CentStep;

    Double_t NPart=Npart(Cent1,Cent2);    

    //yield=NPartVsNColl(NPart);

    yield=NColl(Cent1,Cent2);
    
    /*
    if(PtRange==1 && isUps ==0){ yy = CNMVsNPart(NPart,grShJPsiVsNPart_Y2440);}
    if(PtRange==2 && isUps ==0){ yy = CNMVsNPart(NPart,grShJPsiPtCutVsNPart_Y24);}
    if(PtRange==3 && isUps ==0){ yy = CNMVsNPart(NPart,grShJPsiPtCutVsNPart_Y20);}
    
    if(PtRange==1 && isUps ==1){ yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}
    if(PtRange==1 && isUps ==2){ yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}
    if(PtRange==1 && isUps ==3){ yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}
    */

    ShadowingFactor = yy;
    PionDissFactor = PionDiss_All(Pt);
    GluonDissFactor = IntDiss_All(Pt,NPart);

    sum = sum + (ShadowingFactor*yield*GluonDissFactor*PionDissFactor);

    //sum = sum + (ShadowingFactor*GluonDissFactor*PionDissFactor);

    sumPt = sumPt+yield;  
    
    //cout<<" NPart "<<NPart<<" Pt "<<Pt<<"   "<<"  "<< ShadowingFactor <<"  "<<PionDissFactor<<" "<< GluonDissFactor <<"  "<<sum<<"   "<<sumPt<<endl;
  
  }
  
  
  //cout<<" RAA "<<sum/sumPt<<endl;

  return sum/sumPt;
  //return sum;
}


Double_t IntDiss_PtInt(Double_t PtMin, Double_t NPart)
{
  //for 500 bins and step = 0.1
  //PtMin=0.05;
  //PtMin=6.45;
  
  //for 100 bins and step = 0.5
  //PtMin=0.25;
  //PtMin=6.25;
  
  Double_t Ptmax=0;
  if(isUps==0)Ptmax=40.75;
  if(isUps==0 && PtMin ==0.25)Ptmax=9.25;
  if(isUps==1|| isUps ==2 || isUps ==3)Ptmax=40.25;
  
  Double_t Ptstep=1.0;
  int NN_Pt= (int)((Ptmax-PtMin)/Ptstep);
  Double_t sum=0;
  Double_t sumPt=0;
  Double_t yield =0;
  Double_t yy = 0.0;
  

  /*
  if(PtMin==0.25 &&  isUps ==0){yy = 1.0 - ((0.07/111)*(NPart-2));}
  if(PtMin==6.25 && isUps ==0){ yy = 1.0 + ((0.034/111)*(NPart-2));}
  if(PtMin==9.25 && isUps ==0){ yy = 1.0 + ((0.034/111)*(NPart-2));}

  if(PtMin==0.25 &&  isUps ==1){ yy = 1.0 - ((0.045/111)*(NPart-2));}
  if(PtMin==0.25 &&  isUps ==2){ yy = 1.0 - ((0.045/111)*(NPart-2));}
  if(PtMin==0.25 &&  isUps ==3){ yy = 1.0 - ((0.045/111)*(NPart-2));}
  */

  if(PtMin==0.25 &&  isUps ==0){ yy = CNMVsNPart(NPart,grShJPsiVsNPart_Y2440);}
  if(PtMin==6.25 && isUps ==0){ yy = CNMVsNPart(NPart,grShJPsiPtCutVsNPart_Y24);}
  if(PtMin==9.25 && isUps ==0){ yy = CNMVsNPart(NPart,grShJPsiPtCutVsNPart_Y20);}

  if(PtMin==0.25 &&  isUps ==1){yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}
  if(PtMin==0.25 &&  isUps ==2){yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}
  if(PtMin==0.25 &&  isUps ==3){yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}

 
  if(PtMin==0.25 &&  isUps ==1){yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}
  if(PtMin==0.25 &&  isUps ==2){yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}
  if(PtMin==0.25 &&  isUps ==3){yy = CNMVsNPart(NPart,grShUpsilonVsNPartMid);}

  //grShUpsilonVsNPartFor



  Double_t ShadowingFactor=0;
  Double_t PionDissFactor=0;
  Double_t GluonDissFactor=0;

  for(int i=0;i<=NN_Pt;i++) {
    
    Double_t Pt=PtMin+i*Ptstep;

    if(isUps==0){yield=Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(Pt));} // This should be multiplied by Ncoll (which cancel in numrator and denomenator in this case).
    
    if(isUps==1){yield=Y1S_Pt->GetBinContent(Y1S_Pt->FindBin(Pt));}
    if(isUps==2){yield=Y2S_Pt->GetBinContent(Y2S_Pt->FindBin(Pt));}
    if(isUps==3){yield=Y2S_Pt->GetBinContent(Y2S_Pt->FindBin(Pt));}
    
    ShadowingFactor=yy;
    PionDissFactor = PionDiss_All(Pt);
    GluonDissFactor = IntDiss_All(Pt,NPart);
    //sum = sum + (yy*yield*IntDiss_All(Pt)*PionDissFactor);
    //sum = sum + (ShadowingFactor*yield*GluonDissFactor*PionDissFactor);
    sum = sum + (ShadowingFactor*yield*GluonDissFactor*PionDissFactor*Pt);
    //cout<<" pT "<<Pt<<"   "<<sum <<"  "<< ShadowingFactor <<"  "<<yield<<"  "<<PionDissFactor<<" "<< GluonDissFactor <<endl;
    //sumPt = sumPt+yield;  
    sumPt = sumPt+(yield*Pt);  
  

  }
  return sum/sumPt;
}




Double_t IntDiss_All(Double_t PtMin, Double_t NPart)
{
  Double_t sumd = 0.0;
  for (int i= 0; i<NTau; i++) {
    
    Double_t theta =1.0;
    if(Tau[i] < FormTau * (PtMin/mJpsi) ) theta =0;
    Double_t NPartRatio = TMath::Sqrt(NPart/nPart0);
    sumd = sumd + NPartRatio*fQGP[i]*theta*SigmaQuasiGluon(1,PtMin,TempTau[i]);
  }
  return exp(-sumd*stepTau/hbarc);
}



Double_t IntDiss_All(Double_t PtMin)
{
  Double_t sumd = 0.0;
  for (int i= 0; i<NTau; i++) {
    Double_t theta =1.0;
    if(Tau[i] < FormTau * (PtMin/mJpsi) ) theta =0;
    sumd = sumd +fQGP[i]*theta*SigmaQuasiGluon(1,PtMin,TempTau[i]);
  }
  return exp(-sumd*stepTau/hbarc);
}





Double_t SigmaQuasiGluon(Int_t Flag , Double_t P1, Double_t T)
{
  double m1Square =  mJpsi*mJpsi;
  double m2Square =  mg*mg;
 
  Double_t E1 = TMath::Sqrt( P1*P1 + m1Square);


  Double_t P2Step = 0.01;
  Double_t P2Start = 0.00001;
  Double_t P2End = 8.0;
  int P2N =(int)((P2End - P2Start)/P2Step);
  
  Double_t CosThetaStep = 0.1;
  Double_t CosThetaStart = 1.0;
  Double_t CosThetaEnd = -1.0;
  int CosThetaN =(int)((CosThetaStart - CosThetaEnd)/CosThetaStep);
  
  Double_t SumP2=0;  
  
  for(int i =0; i <= P2N ; i++) {
    Double_t P2 = P2Start + i*P2Step;
    Double_t E2 = TMath::Sqrt(P2*P2 + m2Square);
    double fgluon = fGluon(P2,T);

    Double_t SumCosTheta=0;
    for(int k =0; k <= CosThetaN ; k++)  {
      Double_t CosTheta = CosThetaStart - k*CosThetaStep;
      Double_t fourP1DotP2 = E1*E2-P1*P2*CosTheta;
      Double_t S = m1Square + m2Square + 2*fourP1DotP2;
      
      //if(S<ThetaVal) continue;
      
      Double_t VRel =TMath::Sqrt(fourP1DotP2*fourP1DotP2 - m1Square*m2Square)/(E1*E2);
      Double_t SGD = 0;
      
      
      //Modified cross section
      if(Flag ==1) {SGD = SigmaGluonDissSMod(S);}
      //if(Flag ==1) {SGD = SigmaGluonDissS(S);}
      
      SumCosTheta = SumCosTheta + SGD*VRel;
    }//theta
    
    
    double IntCosTheta = SumCosTheta*CosThetaStep;
    
    ////cout<<"P2 :"<<P2<<"   "<<P2*P2*fgluon*IntCosTheta<<endl;

    SumP2 = SumP2 + P2*P2*fgluon*IntCosTheta;  
  

  }//P2

  Double_t IntP2 = (SumP2*P2Step)/(4.0*pi2);
  return IntP2;
}






Double_t SigmaGluonDissSMod(Double_t s)
{
  Double_t q0 = (s-(mJpsi*mJpsi+mg*mg))/(2.0*mJpsi); 

  Double_t SigmaDq0 = 0.0;
  Double_t SigmaDq0_Mag = 0.0;


  if( (QQbarVar == 1 || QQbarVar == 4)  && q0 > ep0) {SigmaDq0 = epConMod*pow( (q0/ep0-1.0), 1.5)/pow(q0/ep0,5.0);}  

  if( (QQbarVar == 1 || QQbarVar == 4)  && q0 > ep0) 
    {
      Double_t Const = 8.0/3.0;
      Double_t gsSq = 16.0*pi*sqrt(ep0)/(3.0*sqrt(mQ));
      SigmaDq0_Mag = Const*gsSq*pow(ep0,2.5)*sqrt(q0-ep0)/(pow(q0,3.0)*mQ*mQ);
    }  

  if( (QQbarVar == 2 || QQbarVar == 5 ||  QQbarVar == 6)  && q0 > ep0 )  {SigmaDq0 = epConMod*pow( (q0/ep0 - 1.0), 1.5) * pow( (q0/ep0 - 3.0), 2.0)  /pow((q0/ep0),7.0);}  

  if( (QQbarVar == 2 || QQbarVar == 5 || QQbarVar == 6)  && q0 > ep0) 
    {
      Double_t Const = 32.0/3.0;
      Double_t gsSq = 16.0*pi*sqrt(4.0*ep0)/(3.0*sqrt(mQ));
      SigmaDq0_Mag = Const*gsSq*pow(ep0,2.5)*sqrt(q0-ep0)*pow((q0-2.0*ep0),2.0)/(pow(q0,5.0)*mQ*mQ);
    }  

  if( (QQbarVar == 7) && q0 > ep0 )  {SigmaDq0 = (epConMod * TMath::Power(ep0,7.0/2.0)*TMath::Sqrt(q0-ep0)*(9.0*q0*q0 - 20.0*q0*ep0 + 12.0*ep0*ep0))/TMath::Power(q0,7.0);}  

  //if((QQbarVar == 7) && q0 > ep0 )  {
  //Double_t const1 = 2.0*pi/3.0;
  //Double_t const2 = 32.0*32.0/(3.0*3.0);
  //Double_t const3 = TMath::Power(4.0*ep0*mQ,-0.5)/mQ;
  //Double_t xx = q0/ep0;
  //SigmaDq0 =  const1*const2*const3* 4.0*TMath::Sqrt(xx-1)*(9.0*xx*xx - 20.0*xx + 12) /TMath::Power(xx,7.0);
  //}  

  if( (QQbarVar == 7)  && q0 > ep0) 
    {
      Double_t Const = pow(2.0,7.0)/9.0;
      Double_t gsSq = 16.0*pi*sqrt(1.10)/(3.0*sqrt(mQ));
      SigmaDq0_Mag = Const*gsSq*pow(ep0,3.5)*pow((q0-ep0),1.5)/(pow(q0,5.0)*mQ*mQ);
    }  

  
  Double_t TotalSigma = Fac*SigmaDq0;

  if(AddMag==1) {TotalSigma = Fac*(SigmaDq0 + SigmaDq0_Mag);}
  if(AddMag==-1) {TotalSigma = Fac*(SigmaDq0_Mag);}

  return TotalSigma;

}








Double_t fGluon(Double_t pg, Double_t T)
{
  Double_t gi=16.0;
  Double_t fg = gi / (TMath::Exp( sqrt(pg*pg + mg*mg) / T ) - 1 );
  return fg;
}


 /*
Double_t fGluon(Double_t p, Double_t T)
{

  double nn = 6.0;

  double Tem = T;
  double Mass = 0.0 ;

  double MT = sqrt(p*p + Mass*Mass) ;

  //  double fac = (1/(2.0*pi))*(nn-1.0)*(nn-2.0)/((nn*Tem + Mass*(nn-1.0))*(nn*Tem+Mass));
  //  double Norm = fac*pow((nn*Tem/(nn*Tem+Mass)),-nn);

  double dNdy =16.0;
  double Norm = 1.0;
 
  double tsallis = dNdy * Norm * pow( (1 + MT/(nn*Tem)), -nn) ;

  return tsallis ;
  
}
 */



Double_t fcharm_thermal(Double_t p, Double_t T)
{
  Double_t gi=1.0;
  Double_t fq = gi /(TMath::Exp(sqrt(p*p + mQ2)/T) + 1.0);
  //Double_t fq = gi/(TMath::Exp(p/T)+1.0);
  return fq;
}



Double_t fcharm(Double_t p, Double_t T)
{

  double nn = 12.0;

  double Tem = T;
  double Mass = mQ ;

  double MT = sqrt(p*p + Mass*Mass) ;

  //  double fac = (1/(2.0*pi))*(nn-1.0)*(nn-2.0)/((nn*Tem + Mass*(nn-1.0))*(nn*Tem+Mass));
  //  double Norm = fac*pow((nn*Tem/(nn*Tem+Mass)),-nn);

  double dNdy =1.0;
  double Norm = 1.0;
 
  double tsallis = dNdy * Norm * pow( (1 + MT/(nn*Tem)), -nn) ;

  return tsallis ;
  
}


/////////////////////////////////////////////////

Double_t RegenratedQuarkoniaPt(Double_t CentMin, Double_t CentMax, Double_t Pt) 
{
  
  Double_t CentStep=10;
  Int_t NN_Cent= (Int_t)((CentMax-CentMin)/CentStep);
  
  
  Double_t RegenQuarkNum=0.0;
  Double_t RegenQuarkDeno=0.0;
  Double_t RegenQuark=0.0;
  
  for(Int_t i=0;i<NN_Cent;i++) {
    
    Int_t Cent1 = CentMin + i*CentStep;
    Int_t Cent2 = CentMin + (i+1)*CentStep;

    //cout<<"Bin1 "<<Cent1<<" Bin2  "<<Cent2<<endl;

    
    Double_t NPart  =Npart(Cent1,Cent2);    
    Double_t CollNo =NColl(Cent1,Cent2);
       
    Double_t R0ALICE = R05*TMath::Power(NPart/NPart05,0.5); 
    Double_t NCCALICE=nCC0_21*CollNo/NColl05;
   
    Double_t NJPsiALICE=nJpsi0_21*CollNo/NColl05;
   
    Double_t GluonDissALICEPt = IntDiss_All(Pt,NPart);
    Double_t IntForm = IntFormVsPt(Pt, R0ALICE, NPart);
    

    Double_t NJPsiRegenALICEPt = GluonDissALICEPt*NCCALICE*NCCALICE*IntForm;
    
    //cout<<"  gluon diss "<<GluonDissALICEPt<<"   "<<NCCALICE*NCCALICE<<"   "<<IntForm<< "   " <<NPart<<endl;

    RegenQuarkNum = RegenQuarkNum + NJPsiRegenALICEPt;   
    
    Double_t NJPsi0ALICEPt= Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(Pt))/Jpsi_Pt->GetBinWidth(0);
         
    RegenQuarkDeno = RegenQuarkDeno + ( NJPsiALICE*JPsiIntSh_21*NJPsi0ALICEPt);

    //cout<<"PT: "<< Pt <<"  "<< NJPsiRegenALICEPt <<"  "<<( NJPsiALICE* JPsiIntSh_21 * NJPsi0ALICEPt)<<endl;

}

  //cout<<" Num: "<<RegenQuarkNum<<" Deno: "<<RegenQuarkDeno<<endl;

  RegenQuark=(RegenQuarkNum/RegenQuarkDeno);

  return  RegenQuark;
}






/////////////////////////////////////////////////
Double_t RegenratedQuarkoniaNPart(Double_t NCC, Double_t NJPsi, Double_t R0, Double_t IntSh, Int_t PtRange, Double_t NPart) 
{
  
  int NN=17; 
  Double_t Pt[17]={0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25};
  Double_t PtTemp[8]={6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};
  Double_t PtTemp2[3]={9.25,9.75,10.25};

  if(isUps ==0 && PtRange==1)
    {
      NN=8;
      for(int i=0; i<NN; i++) {
	Pt[i]=PtTemp[i];
      }
    }

if(isUps ==0 && PtRange==2)
    {
      NN=3;
      for(int i=0; i<NN; i++) {
	Pt[i]=PtTemp2[i];
      }
    }


  Double_t RegenQuarkNum=0.0;
  Double_t RegenQuarkDeno=0.0;
  Double_t RegenQuark=0.0;

  for(int i=0; i<NN; i++) {
    
    Double_t GluonDissALICEPt=IntDiss_All(Pt[i],NPart);
    Double_t IntForm = IntFormVsPt(Pt[i],R0,NPart);
    
    Double_t NJPsiRegenALICEPt=GluonDissALICEPt*NCC*NCC*IntForm;
    
    //cout<<GluonDissALICEPt<<"   "<<NCC*NCC<<"   "<<IntForm<<endl;

    //RegenQuarkNum = RegenQuarkNum + (NJPsiRegenALICEPt);   
    RegenQuarkNum = RegenQuarkNum + (NJPsiRegenALICEPt*Pt[i]);   
    
    Double_t NJPsi0ALICEPt= Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(Pt[i]))/Jpsi_Pt->GetBinWidth(0);
         
    //RegenQuarkDeno = RegenQuarkDeno + ( NJPsi*IntSh * NJPsi0ALICEPt);
    RegenQuarkDeno = RegenQuarkDeno + ( NJPsi*IntSh * NJPsi0ALICEPt*Pt[i]);

    //cout<<"PT: "<< Pt[i] <<"  "<< NJPsiRegenALICEPt<<"  "<<( NJPsi* IntSh * NJPsi0ALICEPt)<<endl;

}


  //cout<<" Num: "<<RegenQuarkNum<<" Deno: "<<RegenQuarkDeno<<endl;

  RegenQuark=(RegenQuarkNum/RegenQuarkDeno);

  return  RegenQuark;
}




Double_t RegenratedQuarkoniaNPart(Double_t NCC, Double_t NJPsi, Double_t R0, Double_t IntSh, Int_t PtRange) 
{
  
  int NN=17; 
  Double_t Pt[17]={0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25};
  Double_t PtTemp[8]={6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};
  Double_t PtTemp2[3]={9.25,9.75,10.25};
  
  if(PtRange==1)
    {
      NN=8;
      for(int i=0; i<NN; i++) {
	Pt[i]=PtTemp[i];
      }
    }


  if(PtRange==2)
    {
      NN=3;
      for(int i=0; i<NN; i++) {
	Pt[i]=PtTemp2[i];
      }
    }

  Double_t RegenQuarkNum=0.0;
  Double_t RegenQuarkDeno=0.0;
  Double_t RegenQuark=0.0;

  for(int i=0; i<NN; i++) {
    
    Double_t GluonDissALICEPt=IntDiss_All(Pt[i]);
    Double_t IntForm = IntFormVsPt(Pt[i],R0);
    Double_t NJPsiRegenALICEPt=GluonDissALICEPt*NCC*NCC*IntForm;
    
    //cout<<GluonDissALICEPt<<"   "<<NCC*NCC<<"   "<<IntForm<<endl;

    RegenQuarkNum = RegenQuarkNum + NJPsiRegenALICEPt;   
    
    Double_t NJPsi0ALICEPt= Jpsi_Pt->GetBinContent(Jpsi_Pt->FindBin(Pt[i]))/Jpsi_Pt->GetBinWidth(0);
         
    RegenQuarkDeno = RegenQuarkDeno + ( NJPsi* IntSh * NJPsi0ALICEPt);

    //cout<<"PT: "<< Pt[i] <<"  "<< NJPsiRegenALICEPt<<"  "<<( NJPsi* IntSh * NJPsi0ALICEPt)<<endl;

}


  //cout<<" Num: "<<RegenQuarkNum<<" Deno: "<<RegenQuarkDeno<<endl;

  RegenQuark=(RegenQuarkNum/RegenQuarkDeno);

  return  RegenQuark;
}









Double_t IntFormVsPt(Double_t PtMin, Double_t R0Cent, Double_t NPart)
{

  Double_t IntDrate[10000];
  Double_t sumd = 0.0;
  
  for (int i= 0; i<NTau; i++) {
    Double_t theta =1.0;
    if(Tau[i] < FormTau * (PtMin/mJpsi) ) theta =0;
    Double_t FacNPart = TMath::Sqrt(NPart/nPart0);
    sumd = sumd +  FacNPart*fQGP[i]*theta*SigmaQuasiGluon(1,PtMin,TempTau[i]);
    IntDrate[i] =  exp(-sumd*stepTau/hbarc); 
  }
  
  Double_t sum=0;

  for (int i= 0; i<NTau; i++) {
    Double_t DeltaY=1.0;
    Double_t VTau=Tau[i]*pi*(R0Cent+0.5*aT*Tau[i]*Tau[i])*(R0Cent+0.5*aT*Tau[i]*Tau[i])*DeltaY;
    
    ///////////// this is how I add
    //cout<<i<<"  "<<NTau<<" trying to acsess Histo Array :  "<<HistRegenJpsiPt[i]->GetName()<<endl;
    //cout<<" bin content  "<<HistRegenJpsiPt[i]->GetBinContent(HistRegenJpsiPt[i]->FindBin(PtMin))<<endl;
    
    Double_t ForRateVal = HistRegenJpsiPt[i]->GetBinContent(HistRegenJpsiPt[1]->FindBin(PtMin));
    
   
    //Double_t Frate = (fQGP[i]*ForRateVal)/(VTau*IntDrate[i]); 

    Double_t Frate = (ForRateVal)/(VTau*IntDrate[i]); 
   
    //cout<<" PtMin "<<PtMin<<"  "<<HistRegenJpsiPt[1]->FindBin(PtMin)<<endl;
    //cout<<" ForRateVal "<<ForRateVal<<"   "<<IntDrate[i]<<endl;

 
    sum = sum + Frate;
  }
  //cout<<" sum "<<sum<<endl;
  //cout<<" going back to main program "<<endl;
  return sum*stepTau*hbarc2;
}












Double_t IntFormVsPt(Double_t PtMin, Double_t R0Cent)
{

  Double_t IntDrate[10000];
  Double_t sumd = 0.0;
  
  for (int i= 0; i<NTau; i++) {
    Double_t theta =1.0;
    if(Tau[i] < FormTau * (PtMin/mJpsi) ) theta =0;
    sumd = sumd +  fQGP[i]*theta*SigmaQuasiGluon(1,PtMin,TempTau[i]);
    IntDrate[i] =  exp(-sumd*stepTau/hbarc); 
  }
  
  Double_t sum=0;

  for (int i= 0; i<NTau; i++) {
    Double_t DeltaY=1.0;
    Double_t VTau=Tau[i]*pi*(R0Cent+0.5*aT*Tau[i]*Tau[i])*(R0Cent+0.5*aT*Tau[i]*Tau[i])*DeltaY;
    
    ///////////// this is how I add
    //cout<<i<<"  "<<NTau<<" trying to acsess Histo Array :  "<<HistRegenJpsiPt[i]->GetName()<<endl;
    //cout<<" bin content  "<<HistRegenJpsiPt[i]->GetBinContent(HistRegenJpsiPt[i]->FindBin(PtMin))<<endl;
    
    Double_t ForRateVal = HistRegenJpsiPt[i]->GetBinContent(HistRegenJpsiPt[1]->FindBin(PtMin));
    
   
    //Double_t Frate = (fQGP[i]*ForRateVal)/(VTau*IntDrate[i]); 

    Double_t Frate = (ForRateVal)/(VTau*IntDrate[i]); 
   
    //cout<<" PtMin "<<PtMin<<"  "<<HistRegenJpsiPt[1]->FindBin(PtMin)<<endl;
    //cout<<" ForRateVal "<<ForRateVal<<"   "<<IntDrate[i]<<endl;

 
    sum = sum + Frate;
  }
  //cout<<" sum "<<sum<<endl;
  //cout<<" going back to main program "<<endl;
  return sum*stepTau*hbarc2;
}




///This is uesd now for formation rate Aug 2014
Double_t FormRateMC_P_New(Double_t Temp, Int_t iTau) 
{
  
  Double_t PMin= 0.0001;
  Double_t PMax= 6.0;

  Double_t ThetaMin=0.0;
  Double_t ThetaMax=TMath::Pi();
  Double_t PhiMin=0.0;
  Double_t PhiMax=2.0*TMath::Pi();
  
  //int nTrials = 500000;
  int nTrials = 100000;
  
  Double_t Sum=0;
  Double_t Sum1=0;

  for(int i =1; i<= nTrials; i++) {
    
    Double_t RanNo1, RanNo2;
    Double_t RanNo3, RanNo4, RanNo5, RanNo6;	
    
    RanNo1=gRandom->Rndm();
    RanNo2=gRandom->Rndm();
    RanNo3=gRandom->Rndm();
    RanNo4=gRandom->Rndm();
    RanNo5=gRandom->Rndm();
    RanNo6=gRandom->Rndm();
    
    Double_t p1=PMin+(PMax-PMin)*RanNo1;
    Double_t p2=PMin+(PMax-PMin)*RanNo2;
    Double_t theta1= ThetaMin + (ThetaMax -ThetaMin )*RanNo3;
    Double_t theta2= ThetaMin + (ThetaMax -ThetaMin )*RanNo4;
    Double_t phi1= PhiMin + (PhiMax -PhiMin )*RanNo5;
    Double_t phi2= PhiMin + (PhiMax -PhiMin )*RanNo6;
      
    Double_t Px1,Py1,Pz1,Px2,Py2,Pz2;
	
    Px1=p1*sin(theta1)*cos(phi1);
    Py1=p1*sin(theta1)*sin(phi1);
    Pz1=p1*cos(theta1);
    TVector3 PC;
    PC.SetXYZ(Px1,Py1,Pz1);
	
    Px2=p2*sin(theta2)*cos(phi2);
    Py2=p2*sin(theta2)*sin(phi2);
    Pz2=p2*cos(theta2);
    TVector3 PCbar;
    PCbar.SetXYZ(Px2,Py2,Pz2);
    
    TVector3 PCCbar = PC + PCbar;
    Double_t PJpsi=PCCbar.Mag();

    Double_t funNum = FormFun(p1,p2,theta1,theta2,phi1,phi2)*p1*p1*p2*p2*sin(theta1)*sin(theta2)*fcharm(p1,Temp)*fcharm(p2,Temp);
    Double_t funDeno = p1*p1*p2*p2*sin(theta1)*sin(theta2)*fcharm(p1,Temp)*fcharm(p2,Temp);
    
    Sum = Sum + funNum;
    Sum1 = Sum1 + funDeno;
    HistRegenJpsiPt[iTau]->Fill(PJpsi,funNum);
  }
  
  Double_t step = ((PMax- PMin)*(PMax- PMin)*(ThetaMax-ThetaMin)*(ThetaMax-ThetaMin)*(PhiMax-PhiMin)*(PhiMax-PhiMin))/nTrials;
  
  //Double_t ForRate  = step*Sum;
  Double_t ForRateDeno  = step*Sum1; 
  HistRegenJpsiPt[iTau]->Scale(hbarc2*step/(ForRateDeno*HistRegenJpsiPt[iTau]->GetBinWidth(0))); 
  return 0;

}




Double_t FormFun(Double_t p1, Double_t p2, Double_t theta1, Double_t theta2, Double_t phi1, Double_t phi2)
{
  Double_t p1dotp2 = p1*p2*TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi1-phi2) + p1*p2*TMath::Cos(theta1)*TMath::Cos(theta2);
  Double_t E1 = sqrt(p1*p1 + mQ2);
  Double_t E2 = sqrt(p2*p2 + mQ2);
  Double_t s = 2*mQ2 + 2*E1*E2 - 2.0*p1dotp2; 
  //  ccbar to Jpsi cross section
  //  Double_t FJPsi = 0;
  //  if(qSquare > 0 && (4*mD*mD - 4*mQ2 - qSquare) > 0) FJPsi = 0.01*sigmapp ; // sigmapp ???
  
  ///Modified cross section
  //Double_t FJPsi = SigmaFS(s);
  Double_t FJPsi = SigmaFSMod(s);
  
  //Double_t FJPsi = 1.0;
  // Relative Velocity 
  Double_t num = s-4.0*mQ2;
  Double_t RelVel = sqrt(s)*sqrt(num)/(2.0*E1*E2);
  Double_t allfactors = RelVel *FJPsi;

  return allfactors;
}




Double_t SigmaFSMod(Double_t s)
{
  
  Double_t q0 = (s-mJpsi*mJpsi)/(2.0*mJpsi); 
  Double_t SigmaF = 0.0;

  //cout<<" q0 "<<q0<<" 4.0 * MQ2 "<<4*mQ2<<endl;
  
  if(q0<ep0) return SigmaF;
  if(s<=4.0*mQ2) return SigmaF;
  
  Double_t kJpsi2 = pow(s-mJpsi*mJpsi, 2.0);
  Double_t kcc2 = s*(s-4.0*mQ2);
  
  Double_t SigmaD =0.0; 
  


  //SigmaD = epConMod*pow( (q0/ep0-1.0), 1.5)/pow(q0/ep0,5.0);
  //if( (QQbarVar == 1 || QQbarVar == 4)  && q0 > ep0 ) {SigmaD = epConMod*pow( (q0/ep0-1.0), 1.5)/pow(q0/ep0,5.0);}  
  //if( (QQbarVar == 2 || QQbarVar == 5)  && q0 > ep0 ) {SigmaD = epConMod*pow( (q0/ep0 - 1.0), 1.5) * pow( (q0/ep0 - 3.0), 2.0)  /pow((q0/ep0),7.0);}  

  SigmaD = SigmaGluonDissSMod(s);
  //Double_t SigmaGluonDissSMod(Double_t s)


  //Double_t Fac = 1.0;

  //SigmaF = Fac*SigmaD*48.0*kJpsi2/(36.0*kcc2);
  SigmaF = SigmaD*48.0*kJpsi2/(36.0*kcc2);
  
  //cout<<SigmaF<<" sigmaF "<<endl;
  return SigmaF;
}



Double_t CNMVsNPart(Double_t NPart, TGraph *Shgrf)
{
  Double_t yy = 0.0; 
  yy=Shgrf->Eval(NPart);
  return yy;
}



Double_t ELossVsPt(Double_t Pt, TF1 *InFunc)
{
  

  //  InPtHist->Smooth(4);

  //Double_t aa = 0.2;
  //Double_t nn = 0.6;
  //Double_t DeltaPt = aa*TMath::Power((Pt),nn);


  Double_t aa = 0.724; //0-100% --> 0.7
  Double_t CC = 2.5;
  Double_t nn = 0.55;
  
  Double_t DeltaPt = 0.0;

  if(Pt > CC) DeltaPt = aa*TMath::Power((Pt-CC),nn);

  
  Double_t Raa=1.0;

  Double_t Yld_PbPb = InFunc->Eval(Pt+DeltaPt);
  Double_t Yld_PP = InFunc->Eval(Pt);

  Raa = Yld_PbPb/Yld_PP;

  return Raa;
}











Double_t CalculateTandf_LatticeEOS(Double_t ssCent, Double_t R0Cent)
{
  stepTau = 0.1;
  NTau=0;
  Double_t CutTemp = 0.0;
  do{
    Tau[NTau] = tau0 + stepTau*NTau;


    Double_t VTau0Cent =  (R0Cent+0.5*aT*tau0*tau0)*(R0Cent+0.5*aT*tau0*tau0)*(z0+vZ*tau0)*pi;
    Double_t VTauCent =  (R0Cent+0.5*aT*Tau[NTau]*Tau[NTau])*(R0Cent+0.5*aT*Tau[NTau]*Tau[NTau])*(z0+vZ*Tau[NTau])*pi;


    Double_t ssTau = (hbarc3*ssCent)*(VTau0Cent/VTauCent);
    TempTau[NTau]=grfSSVsTemp->Eval(ssTau);
    Double_t tt= TempTau[NTau];
    //fQGP[NTau]=TempVsFQGP->Eval(tt);
    fQGP[NTau]=TempVsFQGP2->Eval(tt);
    CutTemp = TempTau[NTau];
    //cout<<"Lattice: NTau "<<NTau<<" ss "<<ssTau<<"  "<<Tau[NTau]<<"   "<<TempTau[NTau]<<"  "<<fQGP[NTau]<<endl;
    NTau++;
  }while(CutTemp > TF);
  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
Double_t Npart(int BinLow, int BinHigh)
{

  const int nbins = 200;
  Double_t NpartArray[nbins]={401.99, 398.783, 396.936, 392.71, 387.901, 383.593, 377.914, 374.546, 367.507, 361.252, 356.05, 352.43, 345.701, 341.584, 
			      335.148, 330.581, 325.135, 320.777, 315.074, 310.679, 306.687, 301.189, 296.769, 291.795, 287.516, 283.163, 277.818, 274.293, 
			      269.29, 265.911, 260.574, 256.586, 252.732, 249.194, 245.011, 241.292, 236.715, 232.55, 229.322, 225.328, 221.263, 218.604, 214.728, 
			      210.554, 206.878, 203.924, 200.84, 196.572, 193.288, 189.969, 186.894, 183.232, 180.24, 177.36, 174.008, 171.222, 168.296, 165.319, 
			      162.013, 158.495, 156.05, 154.218, 150.559, 148.455, 145.471, 142.496, 139.715, 137.395, 134.469, 131.926, 129.817, 127.045, 124.467, 
			      122.427, 119.698, 117.607, 114.543, 112.662, 110.696, 108.294, 105.777, 103.544, 101.736, 99.943, 97.4951, 95.4291, 93.2148, 91.2133, 
			      89.5108, 87.2103, 85.7498, 83.5134, 81.9687, 79.7456, 78.1684, 76.4873, 74.7635, 72.761, 71.0948, 69.6102, 67.7806, 66.2215, 64.5813, 
			      63.0269, 61.4325, 59.8065, 58.2423, 57.2432, 55.8296, 54.2171, 52.8809, 51.3254, 49.9902, 48.6927, 47.5565, 46.136, 44.8382, 43.6345, 
			      42.3964, 41.4211, 39.9681, 39.178, 37.9341, 36.9268, 35.5626, 34.5382, 33.6912, 32.8156, 31.6695, 30.6552, 29.7015, 28.8655, 27.9609, 
			      27.0857, 26.105, 25.3163, 24.4872, 23.6394, 23.0484, 22.2774, 21.4877, 20.5556, 19.9736, 19.3296, 18.5628, 17.916, 17.2928, 16.6546, 16.1131, 
			      15.4013, 14.8264, 14.3973, 13.7262, 13.2853, 12.8253, 12.2874, 11.7558, 11.2723, 10.8829, 10.4652, 9.96477, 9.6368, 9.09316, 8.84175, 
			      8.48084, 8.05694, 7.64559, 7.29709, 7.07981, 6.70294, 6.45736, 6.10284, 5.91788, 5.5441, 5.33311, 5.06641, 4.96415, 4.6286, 4.38214, 
			      4.2076, 4.01099, 3.81054, 3.63854, 3.43403, 3.23244, 3.08666, 2.86953, 2.74334, 2.62787, 2.48354, 2.38115, 2.26822, 2.23137, 2.1665, 
			      2.14264, 2.10636, 2.07358, 2.05422, 2.04126, 2.00954};
  
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NpartArray[i];
  }
  Double_t NPart = sum/(BinHigh-BinLow);
  return NPart;
}


Double_t NColl(int BinLow, int BinHigh)
{ 
  const int nbins = 200;
  Double_t NCollArray[nbins]={1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 
			   1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 
			   1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 
			   751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 
			   505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 
			   328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 
			   205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 
			   122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 
			   68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 
			   36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 
			   18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 
			   8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 
			   3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 
			   1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NCollArray[i];
  }
  Double_t NColl = sum/(BinHigh-BinLow);
  return NColl;
}




Double_t NPartVsNColl(Double_t ValNPart)
{
  const int NN =100;
  Int_t CentMin =0;
  Int_t CentMax = 100;
  Int_t CentStep = (CentMax-CentMin)/NN;


  //Int_t NN = (CentMax-CentMin)/CentStep;

  Double_t ANPart[NN]={0.0};
  Double_t ANColl[NN]={0.0};



  for(int i =0;i<NN;i++)
    {
      Int_t Cent1 = CentMin + i*CentStep;
      Int_t Cent2 = CentMin + (i+1)*CentStep;
      //cout<<" Cent1 "<<Cent1<<"  Cent2 "<<Cent2<<endl;

      ANPart[i]=Npart(2.0*Cent1,2.0*Cent2);
      ANColl[i]=NColl(2.0*Cent1,2.0*Cent2);

      //cout<<ANPart[i]<<"   "<<ANColl[i]<<endl;

    }


  TGraph *grf_NPart_NColl = new TGraph(NN,ANPart,ANColl);
  grf_NPart_NColl->SetMarkerStyle(20);

  //new TCanvas;
  //grf_NPart_NColl->Draw("APL");

  Double_t ValNColl = grf_NPart_NColl->Eval(ValNPart);
  //cout<<" NColl "<<ValNColl<<endl;

  return ValNColl;


}














void Draw_ATLAS_JPsi_RaaVsPt()
{

  //================= ATLAS JPsi Raa Vs Pt ===================================//
 
  const int nbinsPtATLAS = 8;
  Double_t PtATLASLow[nbinsPtATLAS]={9.0,10.0,11.0,12.0,14.0,15.0,20.0,30.0};
  Double_t PtATLASHigh[nbinsPtATLAS]={10.0,11.0,12.0,14.0,15.0,20.0,30.0,40.0};
  


  Double_t PtATLAS[nbinsPtATLAS]={0.0};
  Double_t ErrPtATLAS[nbinsPtATLAS]={0.0};
  

  Double_t RaaPtATLAS[nbinsPtATLAS] = {0.271875, 0.303125, 0.303125, 0.284375, 0.309375, 0.33125, 0.403125, 0.521875}; 
  Double_t RaaPtStatErrATLAS[nbinsPtATLAS] = {0.0,0.0,0.0,0.0,0.0,0.0,0.03,0.08};
  Double_t RaaPtSystErrATLASLow[nbinsPtATLAS] = {0.2, 0.2375, 0.253125, 0.234375, 0.259375, 0.284375, 0.340625, 0.446875};
  Double_t RaaPtSystErrATLASHigh[nbinsPtATLAS] = {0.346875, 0.36875, 0.359375, 0.334375, 0.359375, 0.3875, 0.4625, 0.6};
  

  for(int j=0;j<nbinsPtATLAS;j++){
    PtATLAS[j]=0.5*(PtATLASLow[j]+PtATLASHigh[j]);
  }

  TGraphErrors *grRaaPtATLAS = new TGraphErrors(nbinsPtATLAS, PtATLAS, RaaPtATLAS, ErrPtATLAS, RaaPtStatErrATLAS);  
  grRaaPtATLAS->SetMarkerStyle(24);
  grRaaPtATLAS->SetMarkerSize(1.4);
  grRaaPtATLAS->SetMarkerColor(kBlue);
  grRaaPtATLAS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtATLAS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtATLAS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtATLAS->GetXaxis();
  Xaxis2->SetLimits(6.0,42.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtATLAS->Draw("AP");

  TLine *lh4 = new TLine(6.0,1.0,42.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"ATLAS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"Prompt J/#psi, |y| < 2");  
  //tb->DrawLatex(0.16,0.74,"0-80% centrality");  
  
  TBox *RaaPtJPsiATLASSys[nbinsPtATLAS];
  for(int j=0;j<nbinsPtATLAS;j++){
    //PtATLAS[j]=PtATLASLow[j]+PtATLASHigh[j];
    //  RaaPtJPsiATLASSys[j] = new TBox(PtATLASLow[j],  RaaPtATLAS[j]-RaaPtSystErrATLAS[j], PtATLASHigh[j],  RaaPtATLAS[j]+RaaPtSystErrATLAS[j]);
    RaaPtJPsiATLASSys[j] = new TBox(PtATLASLow[j],  RaaPtSystErrATLASLow[j], PtATLASHigh[j],  RaaPtSystErrATLASHigh[j]);
  }
  
  for(int j=0;j<nbinsPtATLAS;j++){
    RaaPtJPsiATLASSys[j]->SetFillStyle(0000);
    RaaPtJPsiATLASSys[j]->SetLineColor(kBlue);
    RaaPtJPsiATLASSys[j]->Draw("same"); 
  }
  
  /*
  TBox *ATLASGlobalSysJPsiPt;
  ATLASGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  ATLASGlobalSysJPsiPt->SetFillStyle(3001);
  ATLASGlobalSysJPsiPt->SetLineColor(4);
  ATLASGlobalSysJPsiPt->SetFillColor(4);
  ATLASGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtATLAS,"ATLAS Data", "P");  

  }





void Draw_Clone_ATLAS_JPsi_RaaVsPt()
{

  //================= ATLAS JPsi Raa Vs Pt ===================================//
 
  const int nbinsPtATLAS = 8;
  Double_t PtATLASLow[nbinsPtATLAS]={9.0,10.0,11.0,12.0,14.0,15.0,20.0,30.0};
  Double_t PtATLASHigh[nbinsPtATLAS]={10.0,11.0,12.0,14.0,15.0,20.0,30.0,40.0};
  


  Double_t PtATLAS[nbinsPtATLAS]={0.0};
  Double_t ErrPtATLAS[nbinsPtATLAS]={0.0};
  

  Double_t RaaPtATLAS[nbinsPtATLAS] = {0.271875, 0.303125, 0.303125, 0.284375, 0.309375, 0.33125, 0.403125, 0.521875}; 
  Double_t RaaPtStatErrATLAS[nbinsPtATLAS] = {0.0,0.0,0.0,0.0,0.0,0.0,0.03,0.08};
  Double_t RaaPtSystErrATLASLow[nbinsPtATLAS] = {0.2, 0.2375, 0.253125, 0.234375, 0.259375, 0.284375, 0.340625, 0.446875};
  Double_t RaaPtSystErrATLASHigh[nbinsPtATLAS] = {0.346875, 0.36875, 0.359375, 0.334375, 0.359375, 0.3875, 0.4625, 0.6};
  

  for(int j=0;j<nbinsPtATLAS;j++){
    PtATLAS[j]=0.5*(PtATLASLow[j]+PtATLASHigh[j]);
  }

  TGraphErrors *grRaaPtATLAS = new TGraphErrors(nbinsPtATLAS, PtATLAS, RaaPtATLAS, ErrPtATLAS, RaaPtStatErrATLAS);  
  grRaaPtATLAS->SetMarkerStyle(24);
  grRaaPtATLAS->SetMarkerSize(1.4);
  grRaaPtATLAS->SetMarkerColor(kBlue);
  grRaaPtATLAS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtATLAS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtATLAS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtATLAS->GetXaxis();
  //Xaxis2->SetLimits(6.0,42.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtATLAS->Draw("Psame");

  grRaaPtATLAS->GetXaxis()->SetRangeUser(6.0,42.0);

  TLine *lh4 = new TLine(6.0,1.0,42.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"ATLAS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"Prompt J/#psi, |y| < 2");  
  //tb->DrawLatex(0.16,0.74,"0-80% centrality");  
  
  TBox *RaaPtJPsiATLASSys[nbinsPtATLAS];
  for(int j=0;j<nbinsPtATLAS;j++){
    //PtATLAS[j]=PtATLASLow[j]+PtATLASHigh[j];
    //  RaaPtJPsiATLASSys[j] = new TBox(PtATLASLow[j],  RaaPtATLAS[j]-RaaPtSystErrATLAS[j], PtATLASHigh[j],  RaaPtATLAS[j]+RaaPtSystErrATLAS[j]);
    RaaPtJPsiATLASSys[j] = new TBox(PtATLASLow[j],  RaaPtSystErrATLASLow[j], PtATLASHigh[j],  RaaPtSystErrATLASHigh[j]);
  }
  
  for(int j=0;j<nbinsPtATLAS;j++){
    RaaPtJPsiATLASSys[j]->SetFillStyle(0000);
    RaaPtJPsiATLASSys[j]->SetLineColor(kBlue);
    RaaPtJPsiATLASSys[j]->Draw("same"); 
  }
  
  /*
  TBox *ATLASGlobalSysJPsiPt;
  ATLASGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  ATLASGlobalSysJPsiPt->SetFillStyle(3001);
  ATLASGlobalSysJPsiPt->SetLineColor(4);
  ATLASGlobalSysJPsiPt->SetFillColor(4);
  ATLASGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtATLAS,"ATLAS Data", "P");  

  }
















void Draw_ALICE_JPsi_RaaVsPt()
{

  //================= ALICE JPsi Raa Vs Pt ===================================//
 
  const int nbinsPtALICE = 11;

  Double_t PtALICE[nbinsPtALICE]={0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11.0};
  Double_t ErrPtALICEPlus[nbinsPtALICE]={0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,0.5, 1.0};
  Double_t ErrPtALICEMinus[nbinsPtALICE]={0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,0.5, 1.0};
  
  Double_t RaaPtALICE[nbinsPtALICE] = {0.75, 0.73, 0.64, 0.51, 0.47, 0.38, 0.32, 0.3, 0.35,  0.23, 0.35 }; 
  
  Double_t RaaPtStatErrALICEMinus[nbinsPtALICE] = {0.03, 0.02, 0.02, 0.01, 0.02, 0.01, 0.01, 0.02, 0.03, 0.03, 0.03};
  Double_t RaaPtStatErrALICEPlus[nbinsPtALICE] = {0.03, 0.02, 0.02, 0.01, 0.02, 0.01, 0.01, 0.02, 0.03, 0.03, 0.03};
  
  Double_t RaaPtSyst1ErrALICEMinus[nbinsPtALICE] = {0.07, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.04, 0.06, 0.05, 0.08};
  Double_t RaaPtSyst1ErrALICEPlus[nbinsPtALICE] = {0.07, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.04, 0.06, 0.05, 0.08};
  
  Double_t RaaPtSyst2ErrALICEMinus[nbinsPtALICE] = {0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01};
  Double_t RaaPtSyst2ErrALICEPlus[nbinsPtALICE] = {0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01};
  
  Double_t RaaPtSystErrALICEMinus[nbinsPtALICE]={0.0};
  Double_t RaaPtSystErrALICEPlus[nbinsPtALICE]={0.0};
  
  for(int j=0;j<nbinsPtALICE;j++){
    
    RaaPtSystErrALICEMinus[j] = TMath::Sqrt((RaaPtSyst1ErrALICEMinus[j] * RaaPtSyst1ErrALICEMinus[j])+(RaaPtSyst2ErrALICEMinus[j] * RaaPtSyst2ErrALICEMinus[j]));
    RaaPtSystErrALICEPlus[j] = TMath::Sqrt((RaaPtSyst1ErrALICEPlus[j] * RaaPtSyst1ErrALICEPlus[j])+(RaaPtSyst2ErrALICEPlus[j] * RaaPtSyst2ErrALICEPlus[j]));

  }


  TGraphAsymmErrors *grRaaPtALICE = new TGraphAsymmErrors(nbinsPtALICE, PtALICE, RaaPtALICE, ErrPtALICEMinus, ErrPtALICEPlus, RaaPtStatErrALICEMinus, RaaPtStatErrALICEPlus);  
  grRaaPtALICE->SetMarkerStyle(20);
  grRaaPtALICE->SetMarkerSize(1.4);
  grRaaPtALICE->SetMarkerColor(6);
  grRaaPtALICE->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtALICE->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtALICE->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtALICE->GetXaxis();
  Xaxis2->SetLimits(0.0,12.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtALICE->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,12.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");

  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"ALICE Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"J/#psi, 2.5 < y < 4.0, 0-20% centrality");  
  //  tb->DrawLatex(0.16,0.74,"");  

  //tb->DrawLatex(0.55,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  //tb->DrawLatex(0.55,0.16,"#varUpsilon #rightarrow #mu^{+} #mu^{-}, p_{T}^{#varUpsilon} > 0.0 GeV/c");  
  
  TBox *RaaPtJPsiALICESys[nbinsPtALICE];
  for(int j=0;j<nbinsPtALICE;j++){
    //PtALICE[j]=PtALICELow[j]+PtALICEHigh[j];
    //  RaaPtJPsiALICESys[j] = new TBox(PtALICELow[j],  RaaPtALICE[j]-RaaPtSystErrALICE[j], PtALICEHigh[j],  RaaPtALICE[j]+RaaPtSystErrALICE[j]);
    RaaPtJPsiALICESys[j] = new TBox(PtALICE[j]-0.25,  RaaPtALICE[j]-RaaPtSystErrALICEMinus[j], PtALICE[j]+0.25,  RaaPtALICE[j]+RaaPtSystErrALICEPlus[j]);
  }
  
  for(int j=0;j<nbinsPtALICE;j++){
    RaaPtJPsiALICESys[j]->SetFillStyle(0000);
    RaaPtJPsiALICESys[j]->SetLineColor(6);
    RaaPtJPsiALICESys[j]->Draw("same"); 
  }
  
  /*
  TBox *ALICEGlobalSysJPsiPt;
  ALICEGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  ALICEGlobalSysJPsiPt->SetFillStyle(3001);
  ALICEGlobalSysJPsiPt->SetLineColor(4);
  ALICEGlobalSysJPsiPt->SetFillColor(4);
  ALICEGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtALICE,"ALICE Data", "P");  

  }





void Draw_ATLAS_JPsi_RaaVsNpart()
{

  //================= ATLAS JPsi Raa Vs Npart ===================================//
 
  const int nbinsNPartATLAS = 7;

  Double_t NPartATLAS[nbinsNPartATLAS]={21.9341176471, 53.3756862745, 86.2290196078, 131.141960784, 188.86745098, 264.136470588, 358.605490196};
  Double_t ErrNPartATLAS[nbinsNPartATLAS]={0.0};
  

  Double_t RaaNPartATLAS[nbinsNPartATLAS] = {0.865098039216, 0.623529411765, 0.570196078431, 0.457254901961, 0.331764705882, 0.262745098039, 0.193725490196}; 

  Double_t RaaNPartStatErrATLASLow[nbinsNPartATLAS] = {0.802352941176, 0.579607843137, 0.538823529412, 0.432156862745, 0.316078431373, 0.212549019608, 0.181176470588};
  Double_t RaaNPartStatErrATLASHigh[nbinsNPartATLAS] = {0.927843137255, 0.673725490196, 0.604705882353, 0.479215686275, 0.353725490196, 0.309803921569, 0.212549019608};

  Double_t RaaNPartSystErrATLASLow[nbinsNPartATLAS] = {0.742745098039, 0.538823529412, 0.491764705882, 0.39137254902, 0.287843137255, 0.225098039216, 0.16862745098};
  Double_t RaaNPartSystErrATLASHigh[nbinsNPartATLAS] = {0.990588235294, 0.714509803922, 0.651764705882, 0.52, 0.385098039216, 0.297254901961, 0.228235294118};
  


  Double_t RaaNPartStatErrATLAS[nbinsNPartATLAS] = {0.0};
  Double_t RaaNPartSystErrATLAS[nbinsNPartATLAS]={0.0};

  for(int j=0;j<nbinsNPartATLAS;j++){
    RaaNPartStatErrATLAS[j]= 0.5*(RaaNPartStatErrATLASHigh[j] - RaaNPartStatErrATLASLow[j]);
    RaaNPartSystErrATLAS[j]= 0.5*(RaaNPartSystErrATLASHigh[j] - RaaNPartSystErrATLASLow[j]);
    
  }



  TGraphErrors *grRaaNPartATLAS = new TGraphErrors(nbinsNPartATLAS, NPartATLAS, RaaNPartATLAS, ErrNPartATLAS, RaaNPartStatErrATLAS);  
  grRaaNPartATLAS->SetMarkerStyle(24);
  grRaaNPartATLAS->SetMarkerSize(1.4);
  grRaaNPartATLAS->SetMarkerColor(kBlue);
  grRaaNPartATLAS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaNPartATLAS->GetXaxis()->SetTitle("N_{Part}");
  grRaaNPartATLAS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaNPartATLAS->GetXaxis();
  Xaxis2->SetLimits(0.0,400.0);

  Xaxis2->SetMoreLogLabels();

  grRaaNPartATLAS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,400.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");

  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  // tb->DrawLatex(0.16,0.86,"ATLAS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"Prompt J/#psi, |y| < 2, 9 < p_{T} < 40 GeV/c");  
  //  tb->DrawLatex(0.16,0.74,"");  

  //tb->DrawLatex(0.55,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  //tb->DrawLatex(0.55,0.16,"#varUpsilon #rightarrow #mu^{+} #mu^{-}, p_{T}^{#varUpsilon} > 0.0 GeV/c");  
  
  TBox *RaaNPartJPsiATLASSys[nbinsNPartATLAS];
  for(int j=0;j<nbinsNPartATLAS;j++){
    //NPartATLAS[j]=NPartATLASLow[j]+NPartATLASHigh[j];
    //  RaaNPartJPsiATLASSys[j] = new TBox(NPartATLASLow[j],  RaaNPartATLAS[j]-RaaNPartSystErrATLAS[j], NPartATLASHigh[j],  RaaNPartATLAS[j]+RaaNPartSystErrATLAS[j]);
    RaaNPartJPsiATLASSys[j] = new TBox(NPartATLAS[j]-5,  RaaNPartATLAS[j]-RaaNPartSystErrATLAS[j], NPartATLAS[j]+5,  RaaNPartATLAS[j]+RaaNPartSystErrATLAS[j]);
  }
  
  for(int j=0;j<nbinsNPartATLAS;j++){
    RaaNPartJPsiATLASSys[j]->SetFillStyle(0000);
    RaaNPartJPsiATLASSys[j]->SetLineColor(kBlue);
    RaaNPartJPsiATLASSys[j]->Draw("same"); 
  }
  
  /*
  TBox *ATLASGlobalSysJPsiNPart;
  ATLASGlobalSysJPsiNPart = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  ATLASGlobalSysJPsiNPart->SetFillStyle(3001);
  ATLASGlobalSysJPsiNPart->SetLineColor(4);
  ATLASGlobalSysJPsiNPart->SetFillColor(4);
  ATLASGlobalSysJPsiNPart->Draw("same"); 
  */
  //lgd->AddEntry(grRaaNPartATLAS,"ATLAS Data", "P");  

  }




void Draw_ALICE_JPsi_RaaVsNpart()
{

  //================= ALICE JPsi Raa Vs Npart ===================================//
 
  const int nbinsNPartALICE = 9;

  Double_t NPartALICE[nbinsNPartALICE]={359.0, 263.0, 188.0, 131.0, 86.3, 53.6, 30.4, 15.6, 7.6 };
  Double_t ErrNPartALICEPlus[nbinsNPartALICE]={0.0};
  Double_t ErrNPartALICEMinus[nbinsNPartALICE]={0.0};
  

  Double_t RaaNPartALICE[nbinsNPartALICE] = {0.63, 0.64, 0.67, 0.62, 0.69, 0.72, 0.81, 0.92, 0.9}; 

  Double_t RaaNPartStatErrALICEMinus[nbinsNPartALICE] = {0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.04};
  Double_t RaaNPartStatErrALICEPlus[nbinsNPartALICE] = {0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.04};

  Double_t RaaNPartSyst1ErrALICEMinus[nbinsNPartALICE] = {0.03, 0.02, 0.02, 0.03, 0.03, 0.04, 0.05,0.07, 0.09};
  Double_t RaaNPartSyst1ErrALICEPlus[nbinsNPartALICE] = {0.03, 0.02, 0.02, 0.03, 0.03, 0.04, 0.05,0.07, 0.09};

  Double_t RaaNPartSyst2ErrALICEMinus[nbinsNPartALICE] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06,0.07, 0.07};
  Double_t RaaNPartSyst2ErrALICEPlus[nbinsNPartALICE] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06,0.07, 0.07};
  


  Double_t RaaNPartSystErrALICEMinus[nbinsNPartALICE]={0.0};
  Double_t RaaNPartSystErrALICEPlus[nbinsNPartALICE]={0.0};

  for(int j=0;j<nbinsNPartALICE;j++){
   
    RaaNPartSystErrALICEMinus[j] = TMath::Sqrt((RaaNPartSyst1ErrALICEMinus[j] * RaaNPartSyst1ErrALICEMinus[j])+(RaaNPartSyst2ErrALICEMinus[j] * RaaNPartSyst2ErrALICEMinus[j]));
    RaaNPartSystErrALICEPlus[j] = TMath::Sqrt((RaaNPartSyst1ErrALICEPlus[j] * RaaNPartSyst1ErrALICEPlus[j])+(RaaNPartSyst2ErrALICEPlus[j] * RaaNPartSyst2ErrALICEPlus[j]));


  }


  TGraphAsymmErrors *grRaaNPartALICE = new TGraphAsymmErrors(nbinsNPartALICE, NPartALICE, RaaNPartALICE, ErrNPartALICEMinus, ErrNPartALICEPlus, RaaNPartStatErrALICEMinus, RaaNPartStatErrALICEPlus);  
  grRaaNPartALICE->SetMarkerStyle(20);
  grRaaNPartALICE->SetMarkerSize(1.4);
  grRaaNPartALICE->SetMarkerColor(6);
  grRaaNPartALICE->GetYaxis()->SetRangeUser(0,2.0);
  grRaaNPartALICE->GetXaxis()->SetTitle("N_{Part}");
  grRaaNPartALICE->GetYaxis()->SetTitle("R_{AA}");

  TAxis *Xaxis2 = grRaaNPartALICE->GetXaxis();
  Xaxis2->SetLimits(0.0,400.0);
  Xaxis2->SetMoreLogLabels();

  grRaaNPartALICE->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,400.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");

  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  
  //tb->DrawLatex(0.16,0.86,"ALICE Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"J/#psi, 2.5 < y < 4.0, 0.3 < p_{T} < 8 GeV/c");  
  //tb->DrawLatex(0.16,0.74,"");  
  //tb->DrawLatex(0.55,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  //tb->DrawLatex(0.55,0.16,"#varUpsilon #rightarrow #mu^{+} #mu^{-}, p_{T}^{#varUpsilon} > 0.0 GeV/c");  
  
  TBox *RaaNPartJPsiALICESys[nbinsNPartALICE];
  for(int j=0;j<nbinsNPartALICE;j++){
    //NPartALICE[j]=NPartALICELow[j]+NPartALICEHigh[j];
    //  RaaNPartJPsiALICESys[j] = new TBox(NPartALICELow[j],  RaaNPartALICE[j]-RaaNPartSystErrALICE[j], NPartALICEHigh[j],  RaaNPartALICE[j]+RaaNPartSystErrALICE[j]);
    RaaNPartJPsiALICESys[j] = new TBox(NPartALICE[j]-5,  RaaNPartALICE[j]-RaaNPartSystErrALICEMinus[j], NPartALICE[j]+5,  RaaNPartALICE[j]+RaaNPartSystErrALICEPlus[j]);
  }
  
  for(int j=0;j<nbinsNPartALICE;j++){
    RaaNPartJPsiALICESys[j]->SetFillStyle(0000);
    RaaNPartJPsiALICESys[j]->SetLineColor(6);
    RaaNPartJPsiALICESys[j]->Draw("same"); 
  }
  
  /*
  TBox *ALICEGlobalSysJPsiNPart;
  ALICEGlobalSysJPsiNPart = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  ALICEGlobalSysJPsiNPart->SetFillStyle(3001);
  ALICEGlobalSysJPsiNPart->SetLineColor(4);
  ALICEGlobalSysJPsiNPart->SetFillColor(4);
  ALICEGlobalSysJPsiNPart->Draw("same"); 
  */
  //lgd->AddEntry(grRaaNPartALICE,"ALICE Data", "P");  

  }











void Draw_CMS_Y1S_5TeV_RaaVsPt()
{

  //================= CMS JPsi Raa Vs Pt ===================================//
  const int nbinsPtCMS = 6;
  Double_t PtCMSLow[nbinsPtCMS]={0.0,2.0,4.0,6.0,9.0,12.0};
  Double_t PtCMSHigh[nbinsPtCMS]={2.0,4.0,6.0,9.0,12.0,30.0};
  
  Double_t PtCMS[nbinsPtCMS]={0.0};
  Double_t ErrPtCMS[nbinsPtCMS]={0.0};
  Double_t RaaPtStatErrCMS[nbinsPtCMS]={0.0};

  Double_t RaaPtCMS[nbinsPtCMS] = {0.310,0.336,0.359,0.395,0.421,0.415}; 
  
  Double_t RaaPtStatErrCMSLow[nbinsPtCMS] = {0.282,0.313,0.349,0.372,0.403,0.392};
  Double_t RaaPtStatErrCMSHigh[nbinsPtCMS] = {0.341,0.362,0.377,0.413,0.444,0.441};
  
  Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.246,0.308,0.331,0.364,0.377,0.382};
  Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.372,0.367,0.390,0.423,0.469,0.451};
  

  for(int j=0;j<nbinsPtCMS;j++){
    PtCMS[j]=0.5*(PtCMSLow[j]+PtCMSHigh[j]);
    cout<<PtCMS[j]<<endl;
  
    RaaPtStatErrCMS[j]=0.5*( (RaaPtCMS[j]-RaaPtStatErrCMSLow[j]) + (RaaPtStatErrCMSHigh[j]-RaaPtCMS[j])  );

  }

  TGraphErrors *grRaaPtCMS = new TGraphErrors(nbinsPtCMS, PtCMS, RaaPtCMS, ErrPtCMS, RaaPtStatErrCMS);  
  grRaaPtCMS->SetMarkerStyle(20);
  grRaaPtCMS->SetMarkerSize(1.4);
  grRaaPtCMS->SetMarkerColor(kRed+2);
  grRaaPtCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtCMS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,32.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,32.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"#varUpsilon(1S), |y| < 2.4");  
  //tb->DrawLatex(0.16,0.74,"0-100% centrality");  
  
  TBox *RaaPtJPsiCMSSys[nbinsPtCMS];
  for(int j=0;j<nbinsPtCMS;j++){
    //PtCMS[j]=PtCMSLow[j]+PtCMSHigh[j];
    //  RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtCMS[j]-RaaPtSystErrCMS[j], PtCMSHigh[j],  RaaPtCMS[j]+RaaPtSystErrCMS[j]);
    RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtSystErrCMSHigh[j]);
  }
  
  for(int j=0;j<nbinsPtCMS;j++){
    RaaPtJPsiCMSSys[j]->SetFillStyle(0000);
    RaaPtJPsiCMSSys[j]->SetLineColor(kRed+2);
    RaaPtJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiPt;
  CMSGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiPt->SetFillStyle(3001);
  CMSGlobalSysJPsiPt->SetLineColor(4);
  CMSGlobalSysJPsiPt->SetFillColor(4);
  CMSGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtCMS,"CMS Data", "P");  

  }








void Draw_CMS_Y2S_5TeV_RaaVsPt()
{

  //================= CMS JPsi Raa Vs Pt ===================================//
  const int nbinsPtCMS = 3;
  Double_t PtCMSLow[nbinsPtCMS]={0.0,4.0,9.0};
  Double_t PtCMSHigh[nbinsPtCMS]={4.0,9.0,30.0};
  
  Double_t PtCMS[nbinsPtCMS]={0.0};
  Double_t ErrPtCMS[nbinsPtCMS]={0.0};
  Double_t RaaPtStatErrCMS[nbinsPtCMS]={0.0};

  Double_t RaaPtCMS[nbinsPtCMS] = {0.100,0.123,0.131}; 
  
  Double_t RaaPtStatErrCMSLow[nbinsPtCMS] = {0.064,0.064,0.087};
  Double_t RaaPtStatErrCMSHigh[nbinsPtCMS] = {0.141,0.174,0.172};
  
  Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.028,0.079,0.097};
  Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.179,0.159,0.162};
  

  for(int j=0;j<nbinsPtCMS;j++){
    PtCMS[j]=0.5*(PtCMSLow[j]+PtCMSHigh[j]);
    cout<<PtCMS[j]<<endl;
  
    RaaPtStatErrCMS[j]=0.5*( (RaaPtCMS[j]-RaaPtStatErrCMSLow[j]) + (RaaPtStatErrCMSHigh[j]-RaaPtCMS[j])  );

  }

  TGraphErrors *grRaaPtCMS = new TGraphErrors(nbinsPtCMS, PtCMS, RaaPtCMS, ErrPtCMS, RaaPtStatErrCMS);  
  grRaaPtCMS->SetMarkerStyle(21);
  grRaaPtCMS->SetMarkerSize(1.4);
  grRaaPtCMS->SetMarkerColor(kBlue+2);
  grRaaPtCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtCMS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,32.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,32.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"#varUpsilon(2S), |y| < 2.4");  
  //tb->DrawLatex(0.16,0.74,"0-100% centrality");  
  
  TBox *RaaPtJPsiCMSSys[nbinsPtCMS];
  for(int j=0;j<nbinsPtCMS;j++){
    //PtCMS[j]=PtCMSLow[j]+PtCMSHigh[j];
    //  RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtCMS[j]-RaaPtSystErrCMS[j], PtCMSHigh[j],  RaaPtCMS[j]+RaaPtSystErrCMS[j]);
    RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtSystErrCMSHigh[j]);
  }
  
  for(int j=0;j<nbinsPtCMS;j++){
    RaaPtJPsiCMSSys[j]->SetFillStyle(0000);
    RaaPtJPsiCMSSys[j]->SetLineColor(kBlue+2);
    RaaPtJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiPt;
  CMSGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiPt->SetFillStyle(3001);
  CMSGlobalSysJPsiPt->SetLineColor(4);
  CMSGlobalSysJPsiPt->SetFillColor(4);
  CMSGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtCMS,"CMS Data", "P");  

  }



void Draw_CMS_Y1S_5TeV_RaaVsNpart()
{
  //cat tt | awk '{printf("%.3f,",$2)}'
  //cat tt | awk 'NR%2==0' | awk '{printf("%.3f,",$2)}'
  //================= CMS Y1S Raa Vs Npart ===================================//
 
  const int nbinsNPartCMS = 9;

  Double_t NPartCMS[nbinsNPartCMS]={Npart(140,200),Npart(120,140),Npart(100,120),Npart(80,100),Npart(60,80),Npart(40,60),Npart(20,40),Npart(10,20),Npart(0,10)};
  Double_t ErrNPartCMSPlus[nbinsNPartCMS]={0.0};
  Double_t ErrNPartCMSMinus[nbinsNPartCMS]={0.0};
  

  Double_t RaaNPartCMS[nbinsNPartCMS] = {0.797,0.931,0.616,0.525,0.484,0.405,0.329,0.322,0.319}; 

  Double_t RaaNPartStatErrCMSMinus_temp[nbinsNPartCMS] = {0.660,0.844,0.564,0.487,0.455,0.382,0.311,0.296,0.298};
  Double_t RaaNPartStatErrCMSPlus_temp[nbinsNPartCMS] = {0.929,1.019,0.671,0.560,0.516,0.431,0.343,0.343,0.342};

  Double_t RaaNPartSystErrCMSMinus[nbinsNPartCMS] = {0.640,0.809,0.543,0.470,0.440,0.373,0.299,0.293,0.296};
  Double_t RaaNPartSystErrCMSPlus[nbinsNPartCMS] = {0.961,1.057,0.686,0.575,0.531,0.437,0.355,0.348,0.348};

   
  Double_t RaaNPartStatErrCMSMinus[nbinsNPartCMS] = {0.0};
  Double_t RaaNPartStatErrCMSPlus[nbinsNPartCMS] = {0.0};


  
  for(int j=0;j<nbinsNPartCMS;j++){
   
    RaaNPartStatErrCMSMinus[j] = RaaNPartCMS[j] - RaaNPartStatErrCMSMinus_temp[j];
    RaaNPartStatErrCMSPlus[j] = RaaNPartStatErrCMSPlus_temp[j] - RaaNPartCMS[j];

  



  }


  TGraphAsymmErrors *grRaaNPartCMS = new TGraphAsymmErrors(nbinsNPartCMS, NPartCMS, RaaNPartCMS, ErrNPartCMSMinus, ErrNPartCMSPlus, RaaNPartStatErrCMSMinus, RaaNPartStatErrCMSPlus);  
  grRaaNPartCMS->SetMarkerStyle(20);
  grRaaNPartCMS->SetMarkerSize(1.4);
  grRaaNPartCMS->SetMarkerColor(kRed+2);
  grRaaNPartCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaNPartCMS->GetXaxis()->SetTitle("N_{Part}");
  grRaaNPartCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaNPartCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,400.0);

  Xaxis2->SetMoreLogLabels();

  grRaaNPartCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,400.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");

  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"#varUpsilon(1S), |y| < 2.4, p_{T} < 30 GeV/c");  
  //  tb->DrawLatex(0.16,0.74,"");  

  //tb->DrawLatex(0.55,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  //tb->DrawLatex(0.55,0.16,"#varUpsilon #rightarrow #mu^{+} #mu^{-}, p_{T}^{#varUpsilon} > 0.0 GeV/c");  
  
  TBox *RaaNPartJPsiCMSSys[nbinsNPartCMS];
  for(int j=0;j<nbinsNPartCMS;j++){
    //NPartCMS[j]=NPartCMSLow[j]+NPartCMSHigh[j];
    //  RaaNPartJPsiCMSSys[j] = new TBox(NPartCMSLow[j],  RaaNPartCMS[j]-RaaNPartSystErrCMS[j], NPartCMSHigh[j],  RaaNPartCMS[j]+RaaNPartSystErrCMS[j]);
    RaaNPartJPsiCMSSys[j] = new TBox(NPartCMS[j]-5,  RaaNPartSystErrCMSMinus[j], NPartCMS[j]+5,  RaaNPartSystErrCMSPlus[j]);
  }
  
  for(int j=0;j<nbinsNPartCMS;j++){
    RaaNPartJPsiCMSSys[j]->SetFillStyle(0000);
    RaaNPartJPsiCMSSys[j]->SetLineColor(kRed+2);
    RaaNPartJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiNPart;
  CMSGlobalSysJPsiNPart = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiNPart->SetFillStyle(3001);
  CMSGlobalSysJPsiNPart->SetLineColor(4);
  CMSGlobalSysJPsiNPart->SetFillColor(4);
  CMSGlobalSysJPsiNPart->Draw("same"); 
  */
  //lgd->AddEntry(grRaaNPartCMS,"CMS Data", "P");  

  }




void Draw_CMS_Y2S_5TeV_RaaVsNpart()
{
  //cat tt | awk '{printf("%.3f,",$2)}'
  //cat tt | awk 'NR%2==0' | awk '{printf("%.3f,",$2)}'
  //================= CMS Y1S Raa Vs Npart ===================================//
 
  const int nbinsNPartCMS = 9;

  Double_t NPartCMS[nbinsNPartCMS]={Npart(140,200),Npart(120,140),Npart(100,120),Npart(80,100),Npart(60,80),Npart(40,60),Npart(20,40),Npart(10,20),Npart(0,10)};
  Double_t ErrNPartCMSPlus[nbinsNPartCMS]={0.0};
  Double_t ErrNPartCMSMinus[nbinsNPartCMS]={0.0};
  

  Double_t RaaNPartCMS[nbinsNPartCMS] = {0.517,0.526,0.199,0.213,0.178,0.145,0.092,0.115,0.027}; 

  Double_t RaaNPartStatErrCMSMinus_temp[nbinsNPartCMS] = {0.205,0.354,0.088,0.140,0.119,0.064,0.048,0.062,0.0};
  Double_t RaaNPartStatErrCMSPlus_temp[nbinsNPartCMS] = {0.841,0.698,0.307,0.295,0.236,0.230,0.136,0.168,0.071};

  Double_t RaaNPartSystErrCMSMinus[nbinsNPartCMS] = {0.365,0.400,0.164,0.190,0.160,0.131,0.083,0.103,0.010};
  Double_t RaaNPartSystErrCMSPlus[nbinsNPartCMS] = {0.681,0.651,0.237,0.242,0.198,0.160,0.101,0.124,0.036};

   
  Double_t RaaNPartStatErrCMSMinus[nbinsNPartCMS] = {0.0};
  Double_t RaaNPartStatErrCMSPlus[nbinsNPartCMS] = {0.0};


  
  for(int j=0;j<nbinsNPartCMS;j++){
   
    RaaNPartStatErrCMSMinus[j] = RaaNPartCMS[j] - RaaNPartStatErrCMSMinus_temp[j];
    RaaNPartStatErrCMSPlus[j] = RaaNPartStatErrCMSPlus_temp[j] - RaaNPartCMS[j];

  



  }


  TGraphAsymmErrors *grRaaNPartCMS = new TGraphAsymmErrors(nbinsNPartCMS, NPartCMS, RaaNPartCMS, ErrNPartCMSMinus, ErrNPartCMSPlus, RaaNPartStatErrCMSMinus, RaaNPartStatErrCMSPlus);  
  grRaaNPartCMS->SetMarkerStyle(21);
  grRaaNPartCMS->SetMarkerSize(1.4);
  grRaaNPartCMS->SetMarkerColor(kBlue+2);
  grRaaNPartCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaNPartCMS->GetXaxis()->SetTitle("N_{Part}");
  grRaaNPartCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaNPartCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,400.0);

  Xaxis2->SetMoreLogLabels();

  grRaaNPartCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,400.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");

  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"#varUpsilon(1S), |y| < 2.4, p_{T} < 30 GeV/c");  
  //  tb->DrawLatex(0.16,0.74,"");  

  //tb->DrawLatex(0.55,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  //tb->DrawLatex(0.55,0.16,"#varUpsilon #rightarrow #mu^{+} #mu^{-}, p_{T}^{#varUpsilon} > 0.0 GeV/c");  
  
  TBox *RaaNPartJPsiCMSSys[nbinsNPartCMS];
  for(int j=0;j<nbinsNPartCMS;j++){
    //NPartCMS[j]=NPartCMSLow[j]+NPartCMSHigh[j];
    //  RaaNPartJPsiCMSSys[j] = new TBox(NPartCMSLow[j],  RaaNPartCMS[j]-RaaNPartSystErrCMS[j], NPartCMSHigh[j],  RaaNPartCMS[j]+RaaNPartSystErrCMS[j]);
    RaaNPartJPsiCMSSys[j] = new TBox(NPartCMS[j]-5,  RaaNPartSystErrCMSMinus[j], NPartCMS[j]+5,  RaaNPartSystErrCMSPlus[j]);
  }
  
  for(int j=0;j<nbinsNPartCMS;j++){
    RaaNPartJPsiCMSSys[j]->SetFillStyle(0000);
    RaaNPartJPsiCMSSys[j]->SetLineColor(kBlue+2);
    RaaNPartJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiNPart;
  CMSGlobalSysJPsiNPart = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiNPart->SetFillStyle(3001);
  CMSGlobalSysJPsiNPart->SetLineColor(4);
  CMSGlobalSysJPsiNPart->SetFillColor(4);
  CMSGlobalSysJPsiNPart->Draw("same"); 
  */
  //lgd->AddEntry(grRaaNPartCMS,"CMS Data", "P");  

  }




void Draw_CMS_JPsi_RaaVsPt_0To100()
{

  //================= CMS JPsi Raa Vs Pt ===================================//
 
  const int nbinsPtCMS = 9;
  Double_t PtCMSLow[nbinsPtCMS]={6.5,7.5,8.5,9.5,11.0,13.0,15.0,20.0,30.0};
  Double_t PtCMSHigh[nbinsPtCMS]={7.5,8.5,9.5,11.0,13.0,15.0,20.0,30.0,50.0};
  

  Double_t PtCMS[nbinsPtCMS]={0.0};
  Double_t ErrPtCMS[nbinsPtCMS]={0.0};
  

  Double_t RaaPtCMS[nbinsPtCMS] = {0.348, 0.346, 0.333, 0.331, 0.341, 0.361, 0.361, 0.452, 0.518}; 
  Double_t RaaPtStatErrCMS[nbinsPtCMS] = {0.012, 0.010, 0.009, 0.008, 0.008, 0.011, 0.010, 0.019, 0.043};
  

  Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};
  Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};
  

  for(int j=0;j<nbinsPtCMS;j++){
    PtCMS[j]=0.5*(PtCMSLow[j]+PtCMSHigh[j]);
  }

  TGraphErrors *grRaaPtCMS = new TGraphErrors(nbinsPtCMS, PtCMS, RaaPtCMS, ErrPtCMS, RaaPtStatErrCMS);  
  grRaaPtCMS->SetMarkerStyle(21);
  grRaaPtCMS->SetMarkerSize(1.2);
  grRaaPtCMS->SetMarkerColor(kRed);
  grRaaPtCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtCMS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,52.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,52.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"Prompt J/#psi, |y| < 2");  
  //tb->DrawLatex(0.16,0.74,"0-80% centrality");  
  

  //Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};
  //Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};

  TBox *RaaPtJPsiCMSSys[nbinsPtCMS];
  for(int j=0;j<nbinsPtCMS;j++){
    //PtCMS[j]=PtCMSLow[j]+PtCMSHigh[j];
    RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtCMS[j]-RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtCMS[j]+RaaPtSystErrCMSHigh[j]);
    //RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtSystErrCMSHigh[j]);
  }
  
  for(int j=0;j<nbinsPtCMS;j++){
    RaaPtJPsiCMSSys[j]->SetFillStyle(0000);
    RaaPtJPsiCMSSys[j]->SetLineColor(kRed);
    RaaPtJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiPt;
  CMSGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiPt->SetFillStyle(3001);
  CMSGlobalSysJPsiPt->SetLineColor(4);
  CMSGlobalSysJPsiPt->SetFillColor(4);
  CMSGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtCMS,"CMS Data", "P");  

}






void Draw_CMS_JPsi_RaaVsPt_0To10()
{

  //================= CMS JPsi Raa Vs Pt ===================================//
 
  const int nbinsPtCMS = 8;
  Double_t PtCMSLow[nbinsPtCMS]={6.5,7.5,8.5,9.5,11.0,13.0,15.0,20.0};
  Double_t PtCMSHigh[nbinsPtCMS]={7.5,8.5,9.5,11.0,13.0,15.0,20.0,50.0};
  

  Double_t PtCMS[nbinsPtCMS]={0.0};
  Double_t ErrPtCMS[nbinsPtCMS]={0.0};
  
  Double_t RaaPtCMS[nbinsPtCMS] = {0.254, 0.235, 0.242, 0.234, 0.237, 0.259, 0.277, 0.324}; 
  Double_t RaaPtStatErrCMS[nbinsPtCMS] = {0.017, 0.012, 0.012, 0.010, 0.010, 0.014, 0.013, 0.020};
  
  Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.020, 0.015, 0.011, 0.011, 0.012, 0.021, 0.012, 0.016};
  Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.020, 0.015, 0.011, 0.011, 0.012, 0.021, 0.012, 0.016};
  

  for(int j=0;j<nbinsPtCMS;j++){
    PtCMS[j]=0.5*(PtCMSLow[j]+PtCMSHigh[j]);
  }

  TGraphErrors *grRaaPtCMS = new TGraphErrors(nbinsPtCMS, PtCMS, RaaPtCMS, ErrPtCMS, RaaPtStatErrCMS);  
  grRaaPtCMS->SetMarkerStyle(21);
  grRaaPtCMS->SetMarkerSize(1.2);
  grRaaPtCMS->SetMarkerColor(kRed+2);
  grRaaPtCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtCMS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,52.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,52.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"Prompt J/#psi, |y| < 2");  
  tb->DrawLatex(0.75,0.90,"0-10%");  
  

  //Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};
  //Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};

  TBox *RaaPtJPsiCMSSys[nbinsPtCMS];
  for(int j=0;j<nbinsPtCMS;j++){
    //PtCMS[j]=PtCMSLow[j]+PtCMSHigh[j];
    RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtCMS[j]-RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtCMS[j]+RaaPtSystErrCMSHigh[j]);
    //RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtSystErrCMSHigh[j]);
  }
  
  for(int j=0;j<nbinsPtCMS;j++){
    RaaPtJPsiCMSSys[j]->SetFillStyle(0000);
    RaaPtJPsiCMSSys[j]->SetLineColor(kRed+2);
    RaaPtJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiPt;
  CMSGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiPt->SetFillStyle(3001);
  CMSGlobalSysJPsiPt->SetLineColor(4);
  CMSGlobalSysJPsiPt->SetFillColor(4);
  CMSGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtCMS,"CMS Data", "P");  

}






void Draw_CMS_JPsi_RaaVsPt_10To30()
{

  //================= CMS JPsi Raa Vs Pt ===================================//
 
  const int nbinsPtCMS = 8;
  Double_t PtCMSLow[nbinsPtCMS]={6.5,7.5,8.5,9.5,11.0,13.0,15.0,20.0};
  Double_t PtCMSHigh[nbinsPtCMS]={7.5,8.5,9.5,11.0,13.0,15.0,20.0,50.0};
  

  Double_t PtCMS[nbinsPtCMS]={0.0};
  Double_t ErrPtCMS[nbinsPtCMS]={0.0};
  
  Double_t RaaPtCMS[nbinsPtCMS] = {0.353, 0.340, 0.329, 0.330, 0.341, 0.356, 0.341, 0.478}; 
  Double_t RaaPtStatErrCMS[nbinsPtCMS] = {0.018, 0.014, 0.013, 0.012, 0.012, 0.015, 0.014, 0.025};
  
  Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.024, 0.018, 0.015, 0.014, 0.012, 0.014, 0.019, 0.021};
  Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.024, 0.018, 0.015, 0.014, 0.012, 0.014, 0.019, 0.021};
  

  for(int j=0;j<nbinsPtCMS;j++){
    PtCMS[j]=0.5*(PtCMSLow[j]+PtCMSHigh[j]);
  }

  TGraphErrors *grRaaPtCMS = new TGraphErrors(nbinsPtCMS, PtCMS, RaaPtCMS, ErrPtCMS, RaaPtStatErrCMS);  
  grRaaPtCMS->SetMarkerStyle(20);
  grRaaPtCMS->SetMarkerSize(1.4);
  grRaaPtCMS->SetMarkerColor(kBlue+2);
  grRaaPtCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtCMS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,52.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,52.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"Prompt J/#psi, |y| < 2");  
  //tb->DrawLatex(0.16,0.74,"0-80% centrality");  
   tb->DrawLatex(0.75,0.90,"10-30%"); 

  //Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};
  //Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};

  TBox *RaaPtJPsiCMSSys[nbinsPtCMS];
  for(int j=0;j<nbinsPtCMS;j++){
    //PtCMS[j]=PtCMSLow[j]+PtCMSHigh[j];
    RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtCMS[j]-RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtCMS[j]+RaaPtSystErrCMSHigh[j]);
    //RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtSystErrCMSHigh[j]);
  }
  
  for(int j=0;j<nbinsPtCMS;j++){
    RaaPtJPsiCMSSys[j]->SetFillStyle(0000);
    RaaPtJPsiCMSSys[j]->SetLineColor(kBlue+4);
    RaaPtJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiPt;
  CMSGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiPt->SetFillStyle(3001);
  CMSGlobalSysJPsiPt->SetLineColor(4);
  CMSGlobalSysJPsiPt->SetFillColor(4);
  CMSGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtCMS,"CMS Data", "P");  

}





void Draw_CMS_JPsi_RaaVsPt_30To100()
{

  //================= CMS JPsi Raa Vs Pt ===================================//
 
  const int nbinsPtCMS = 8;
  Double_t PtCMSLow[nbinsPtCMS]={6.5,7.5,8.5,9.5,11.0,13.0,15.0,20.0};
  Double_t PtCMSHigh[nbinsPtCMS]={7.5,8.5,9.5,11.0,13.0,15.0,20.0,50.0};
  

  Double_t PtCMS[nbinsPtCMS]={0.0};
  Double_t ErrPtCMS[nbinsPtCMS]={0.0};
  
  Double_t RaaPtCMS[nbinsPtCMS] = {0.516, 0.565, 0.512, 0.536, 0.564, 0.581, 0.586, 0.739}; 
  Double_t RaaPtStatErrCMS[nbinsPtCMS] = {0.024, 0.025, 0.022, 0.021, 0.023, 0.028, 0.028, 0.046};
  
  Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.034, 0.028, 0.024, 0.022, 0.022, 0.023, 0.022, 0.029};
  Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.034, 0.028, 0.024, 0.022, 0.022, 0.023, 0.022, 0.029};
  

  for(int j=0;j<nbinsPtCMS;j++){
    PtCMS[j]=0.5*(PtCMSLow[j]+PtCMSHigh[j]);
  }

  TGraphErrors *grRaaPtCMS = new TGraphErrors(nbinsPtCMS, PtCMS, RaaPtCMS, ErrPtCMS, RaaPtStatErrCMS);  
  grRaaPtCMS->SetMarkerStyle(34);
  grRaaPtCMS->SetMarkerSize(1.6);
  grRaaPtCMS->SetMarkerColor(kGreen+2);
  grRaaPtCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaPtCMS->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  grRaaPtCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaPtCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,52.0);

  Xaxis2->SetMoreLogLabels();

  grRaaPtCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,52.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"Prompt J/#psi, |y| < 2");  
  //tb->DrawLatex(0.16,0.74,"0-80% centrality");  
  tb->DrawLatex(0.75,0.90,"30-100%"); 

  //Double_t RaaPtSystErrCMSLow[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};
  //Double_t RaaPtSystErrCMSHigh[nbinsPtCMS] = {0.033, 0.026, 0.021, 0.019, 0.018, 0.019, 0.019, 0.024, 0.038};

  TBox *RaaPtJPsiCMSSys[nbinsPtCMS];
  for(int j=0;j<nbinsPtCMS;j++){
    //PtCMS[j]=PtCMSLow[j]+PtCMSHigh[j];
    RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtCMS[j]-RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtCMS[j]+RaaPtSystErrCMSHigh[j]);
    //RaaPtJPsiCMSSys[j] = new TBox(PtCMSLow[j],  RaaPtSystErrCMSLow[j], PtCMSHigh[j],  RaaPtSystErrCMSHigh[j]);
  }
  
  for(int j=0;j<nbinsPtCMS;j++){
    RaaPtJPsiCMSSys[j]->SetFillStyle(0000);
    RaaPtJPsiCMSSys[j]->SetLineColor(kGreen+4);
    RaaPtJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiPt;
  CMSGlobalSysJPsiPt = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiPt->SetFillStyle(3001);
  CMSGlobalSysJPsiPt->SetLineColor(4);
  CMSGlobalSysJPsiPt->SetFillColor(4);
  CMSGlobalSysJPsiPt->Draw("same"); 
  */
  //lgd->AddEntry(grRaaPtCMS,"CMS Data", "P");  

}






void Draw_CMS_JPsi_RaaVsNpart()
{

  //================= CMS Y1S Raa Vs Npart ===================================//
 
  const int nbinsNPartCMS = 13;

  Double_t NPartCMS[nbinsNPartCMS]={Npart(0,10),Npart(10,20),Npart(20,30),Npart(30, 40),Npart(40, 50),Npart(50, 60),Npart(60, 70),Npart(70, 80),Npart(80, 90),
				    Npart(90, 100),Npart(100, 120),Npart(120, 140),Npart(140, 200)};
  Double_t ErrNPartCMSPlus[nbinsNPartCMS]={0.0};
  Double_t ErrNPartCMSMinus[nbinsNPartCMS]={0.0};
  

  Double_t RaaNPartCMS[nbinsNPartCMS] = {0.223, 0.270, 0.304, 0.330, 0.376, 0.397, 0.457, 0.482, 0.569, 0.605, 0.654, 0.739, 0.800}; 

  Double_t RaaNPartStatErrCMSMinus_temp[nbinsNPartCMS] = {0.005, 0.007, 0.008, 0.009, 0.011, 0.013, 0.014, 0.017, 0.023, 0.026, 0.025, 0.041, 0.066};
  Double_t RaaNPartStatErrCMSPlus_temp[nbinsNPartCMS] = {0.005, 0.007, 0.008, 0.009, 0.011, 0.013, 0.014, 0.017, 0.023, 0.026, 0.025, 0.041, 0.066};

  Double_t RaaNPartSystErrCMSMinus[nbinsNPartCMS] = {0.014, 0.017, 0.020, 0.022, 0.025, 0.028, 0.033, 0.038, 0.047, 0.056, 0.067, 0.091, 0.091};
  Double_t RaaNPartSystErrCMSPlus[nbinsNPartCMS] = {0.014, 0.017, 0.018, 0.021, 0.025, 0.028, 0.033, 0.038, 0.047, 0.056, 0.069, 0.099, 0.138};
   
  Double_t RaaNPartStatErrCMSMinus[nbinsNPartCMS] = {0.0};
  Double_t RaaNPartStatErrCMSPlus[nbinsNPartCMS] = {0.0};

  
  for(int j=0;j<nbinsNPartCMS;j++){
  
    RaaNPartStatErrCMSMinus[j] = RaaNPartStatErrCMSMinus_temp[j];
    RaaNPartStatErrCMSPlus[j] = RaaNPartStatErrCMSPlus_temp[j];
  
  }

  TGraphAsymmErrors *grRaaNPartCMS = new TGraphAsymmErrors(nbinsNPartCMS, NPartCMS, RaaNPartCMS, ErrNPartCMSMinus, ErrNPartCMSPlus, RaaNPartStatErrCMSMinus, RaaNPartStatErrCMSPlus);  
  grRaaNPartCMS->SetMarkerStyle(21);
  grRaaNPartCMS->SetMarkerSize(1.4);
  grRaaNPartCMS->SetMarkerColor(2);
  grRaaNPartCMS->GetYaxis()->SetRangeUser(0,2.0);
  grRaaNPartCMS->GetXaxis()->SetTitle("N_{Part}");
  grRaaNPartCMS->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaNPartCMS->GetXaxis();
  Xaxis2->SetLimits(0.0,400.0);
  Xaxis2->SetMoreLogLabels();

  grRaaNPartCMS->Draw("AP");

  TLine *lh4 = new TLine(0.0,1.0,400.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");

  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"CMS Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"#varUpsilon(1S), |y| < 2.4, p_{T} < 30 GeV/c");  
  //  tb->DrawLatex(0.16,0.74,"");  

  //tb->DrawLatex(0.55,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  //tb->DrawLatex(0.55,0.16,"#varUpsilon #rightarrow #mu^{+} #mu^{-}, p_{T}^{#varUpsilon} > 0.0 GeV/c");  
  
  TBox *RaaNPartJPsiCMSSys[nbinsNPartCMS];
  for(int j=0;j<nbinsNPartCMS;j++){
    //NPartCMS[j]=NPartCMSLow[j]+NPartCMSHigh[j];
    //  RaaNPartJPsiCMSSys[j] = new TBox(NPartCMSLow[j],  RaaNPartCMS[j]-RaaNPartSystErrCMS[j], NPartCMSHigh[j],  RaaNPartCMS[j]+RaaNPartSystErrCMS[j]);
    RaaNPartJPsiCMSSys[j] = new TBox(NPartCMS[j]-5,  RaaNPartCMS[j]-RaaNPartSystErrCMSMinus[j], NPartCMS[j]+5,  RaaNPartCMS[j]+ RaaNPartSystErrCMSPlus[j]);
  }
  
  for(int j=0;j<nbinsNPartCMS;j++){
    RaaNPartJPsiCMSSys[j]->SetFillStyle(0000);
    RaaNPartJPsiCMSSys[j]->SetLineColor(2);
    RaaNPartJPsiCMSSys[j]->Draw("same"); 
  }
  
  /*
  TBox *CMSGlobalSysJPsiNPart;
  CMSGlobalSysJPsiNPart = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
  CMSGlobalSysJPsiNPart->SetFillStyle(3001);
  CMSGlobalSysJPsiNPart->SetLineColor(4);
  CMSGlobalSysJPsiNPart->SetFillColor(4);
  CMSGlobalSysJPsiNPart->Draw("same"); 
  */
  //lgd->AddEntry(grRaaNPartCMS,"CMS Data", "P");  

  }




void Draw_ALICE_Y1S_5TeV_RaaVsNpart()
{
  //cat tt | awk '{printf("%.3f,",$2)}'
  //cat tt | awk 'NR%2==0' | awk '{printf("%.3f,",$2)}'
  //================= ALICE Y1S Raa Vs Npart ===================================//
  
  const int nbinsNPartALICE = 4;
  
  Double_t NPartALICE[nbinsNPartALICE]={27, 108, 226, 360};
  Double_t ErrNPartALICEPlus[nbinsNPartALICE]={0.0};
  Double_t ErrNPartALICEMinus[nbinsNPartALICE]={0.0};
  
  
  Double_t RaaNPartALICE[nbinsNPartALICE] = {0.652, 0.556, 0.353, 0.346}; 
  
  Double_t RaaNPartStatErrALICEMinus_temp[nbinsNPartALICE] ={0.523,0.470,0.315,0.311};
  Double_t RaaNPartStatErrALICEPlus_temp[nbinsNPartALICE] =  {0.787,0.645,0.390,0.384};
  
  Double_t RaaNPartSystErrALICEMinus[nbinsNPartALICE] = {0.700,0.590,0.373,0.368};
  Double_t RaaNPartSystErrALICEPlus[nbinsNPartALICE] = {0.605,0.523,0.331,0.322};
  
  
  Double_t RaaNPartStatErrALICEMinus[nbinsNPartALICE] = {0.0};
  Double_t RaaNPartStatErrALICEPlus[nbinsNPartALICE] = {0.0};

  
  for(int j=0;j<nbinsNPartALICE;j++){
    
    RaaNPartStatErrALICEMinus[j] = RaaNPartALICE[j] - RaaNPartStatErrALICEMinus_temp[j];
    RaaNPartStatErrALICEPlus[j] = RaaNPartStatErrALICEPlus_temp[j] - RaaNPartALICE[j];
    
    
  }
  
  
  //cout<<endl;
  //for(int j=0;j<nbinsNPartALICE;j++){ cout<<" Y(1S) ALICE NPart Stat Errors "<< RaaNPartStatErrALICEMinus[j] << "   " << RaaNPartStatErrALICEPlus[j] <<endl; }
  //cout<<endl;
  
  TGraphAsymmErrors *grRaaNPartALICE = new TGraphAsymmErrors(nbinsNPartALICE, NPartALICE, RaaNPartALICE, ErrNPartALICEMinus, ErrNPartALICEPlus, RaaNPartStatErrALICEMinus, RaaNPartStatErrALICEPlus);  
  grRaaNPartALICE->SetMarkerStyle(21);
  grRaaNPartALICE->SetMarkerSize(1.4);
  grRaaNPartALICE->SetMarkerColor(kRed);
  grRaaNPartALICE->SetLineColor(kRed);
  
  grRaaNPartALICE->GetYaxis()->SetRangeUser(0,2.0);
  grRaaNPartALICE->GetXaxis()->SetTitle("N_{Part}");
  grRaaNPartALICE->GetYaxis()->SetTitle("R_{AA}");
  
  TAxis *Xaxis2 = grRaaNPartALICE->GetXaxis();
  Xaxis2->SetLimits(0.0,400.0);
  
  Xaxis2->SetMoreLogLabels();
  
  grRaaNPartALICE->Draw("AP");
  
  TLine *lh4 = new TLine(0.0,1.0,400.0,1.0);
  lh4->SetLineColor(1);
  lh4->SetLineStyle(1);
  lh4->SetLineWidth(2);
  lh4->Draw("same");
  
  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  //tb->DrawLatex(0.16,0.86,"ALICE Data, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //tb->DrawLatex(0.16,0.80,"#varUpsilon(1S), 2.5 < y < 4, p_{T} < 12 GeV/c");  
  //  tb->DrawLatex(0.16,0.74,"");  
  
  //tb->DrawLatex(0.55,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  //tb->DrawLatex(0.55,0.16,"#varUpsilon #rightarrow #mu^{+} #mu^{-}, p_{T}^{#varUpsilon} > 0.0 GeV/c");  
  
  TBox *RaaNPartJPsiALICESys[nbinsNPartALICE];
  for(int j=0;j<nbinsNPartALICE;j++){
    //NPartALICE[j]=NPartALICELow[j]+NPartALICEHigh[j];
    //  RaaNPartJPsiALICESys[j] = new TBox(NPartALICELow[j],  RaaNPartALICE[j]-RaaNPartSystErrALICE[j], NPartALICEHigh[j],  RaaNPartALICE[j]+RaaNPartSystErrALICE[j]);
    RaaNPartJPsiALICESys[j] = new TBox(NPartALICE[j]-5,  RaaNPartSystErrALICEMinus[j], NPartALICE[j]+5,  RaaNPartSystErrALICEPlus[j]);
  }
  
  for(int j=0;j<nbinsNPartALICE;j++){
    RaaNPartJPsiALICESys[j]->SetFillStyle(0000);
    RaaNPartJPsiALICESys[j]->SetLineColor(kRed);
    RaaNPartJPsiALICESys[j]->Draw("same"); 
  }
  
  /*
    TBox *ALICEGlobalSysJPsiNPart;
    ALICEGlobalSysJPsiNPart = new TBox(18-0.2, 1 - 0.083, 18+0.2, 1 + 0.083);
    ALICEGlobalSysJPsiNPart->SetFillStyle(3001);
    ALICEGlobalSysJPsiNPart->SetLineColor(4);
    ALICEGlobalSysJPsiNPart->SetFillColor(4);
    ALICEGlobalSysJPsiNPart->Draw("same"); 
  */
  //lgd->AddEntry(grRaaNPartALICE,"ALICE Data", "P");  
  
}










//////////////////////////////////////////////////////////////////////////////////////////////////
Double_t Npart_276(int BinLow, int BinHigh)
{
  Double_t NpartArray[40]={393.622,368.96,342.32,316.49,293.49,271.98,249.65,230.53,212.28,194.50,178.54,
			 163.25,149.05,135.92,123.28,111.67,100.79,90.71,80.93,72.60,64.15,56.61,49.95,
			 43.39,37.83,32.70,27.86,23.79,20.20,16.85,14.04,11.60,9.55,7.72,6.44,4.96,4.22,
			 3.50,3.17,2.79};
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NpartArray[i];
  }
  Double_t NPart = sum/(BinHigh-BinLow);
  return NPart;
}

Double_t NColl_276(int BinLow, int BinHigh)
{
  Double_t NCollArray[40]={1747.8600, 1567.5300, 1388.3900, 1231.7700, 1098.2000, 980.4390, 861.6090, 766.0420, 676.5150, 593.4730,
			   521.9120, 456.5420, 398.5460, 346.6470, 299.3050, 258.3440, 221.2160, 188.6770, 158.9860, 134.7000,
			   112.5470, 93.4537, 77.9314, 63.5031, 52.0469, 42.3542, 33.9204, 27.3163, 21.8028, 17.2037,
			   13.5881, 10.6538, 8.3555, 6.4089, 5.1334, 3.7322, 3.0663, 2.4193, 2.1190, 1.7695};
  
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NCollArray[i];
  }
  Double_t NColl = sum/(BinHigh-BinLow);
  return NColl;
}



void TotalCharmProductionCross_ALICE()
{

   
  const int nbins = 2;

  Double_t Energy[nbins]={2.76,7.0};
  Double_t ErrEnrgy[nbins]={0.0};
  
  Double_t Sigma[nbins] = {4.8,8.5}; 
  Double_t ErrSigma[nbins] = {0.8,0.5}; 

  Double_t SigmaU[nbins] = {5.6,9.0}; 
  Double_t SigmaD[nbins] = {4.0,8.0}; 

  TGraphErrors *grSigmaVsEnergy = new TGraphErrors(nbins, Energy, Sigma, ErrEnrgy, ErrSigma);  
  grSigmaVsEnergy->SetMarkerStyle(24);
  grSigmaVsEnergy->SetMarkerSize(1.4);
  grSigmaVsEnergy->SetMarkerColor(kBlue);
  grSigmaVsEnergy->GetYaxis()->SetRangeUser(0,10.0);
  grSigmaVsEnergy->GetXaxis()->SetTitle("s^{1/2}");
  grSigmaVsEnergy->GetYaxis()->SetTitle("#sigma_{cc}");
  
  TAxis *Xaxis2 = grSigmaVsEnergy->GetXaxis();
  Xaxis2->SetLimits(0.0,10.0);
  Xaxis2->SetMoreLogLabels();

  TGraph *grSigmaVsEnergyU = new TGraphErrors(nbins, Energy, SigmaU);  
  TGraph *grSigmaVsEnergyD = new TGraphErrors(nbins, Energy, SigmaD);  

  grSigmaVsEnergy->Draw("AP");
  grSigmaVsEnergyU->Draw("Lsame");
  grSigmaVsEnergyD->Draw("Lsame");
  cout<<endl<<endl;

  Double_t ErrU = grSigmaVsEnergyU->Eval(5.0) - grSigmaVsEnergy->Eval(5.0);
  Double_t ErrD = -grSigmaVsEnergyD->Eval(5.0) + grSigmaVsEnergy->Eval(5.0);

  cout<<" Total ccbar cross-section at five TeV is : "<<grSigmaVsEnergy->Eval(5.0)<<" mb "<<" pm "<<ErrU<<"  "<<ErrD<<endl;
  cout<<endl<<endl;

  }



void TotalBeautyProductionCross_ALICE()
{

  
  const int nbins = 2;

  Double_t Energy[nbins]={2.76,7.0};
  Double_t ErrEnrgy[nbins]={0.0};
  
  Double_t Sigma[nbins] = {130.0,282.0}; 
  Double_t ErrSigma[nbins] = {15.1,74.0}; 

  Double_t SigmaU[nbins] = {145.1,356.0}; 
  Double_t SigmaD[nbins] = {115.0,208.0}; 

  TGraphErrors *grSigmaVsEnergy = new TGraphErrors(nbins, Energy, Sigma, ErrEnrgy, ErrSigma);  
  grSigmaVsEnergy->SetMarkerStyle(24);
  grSigmaVsEnergy->SetMarkerSize(1.4);
  grSigmaVsEnergy->SetMarkerColor(kBlue);
  grSigmaVsEnergy->GetYaxis()->SetRangeUser(50,400.0);
  grSigmaVsEnergy->GetXaxis()->SetTitle("s^{1/2}");
  grSigmaVsEnergy->GetYaxis()->SetTitle("#sigma_{bb}");
  
  TAxis *Xaxis2 = grSigmaVsEnergy->GetXaxis();
  Xaxis2->SetLimits(0.0,10.0);
  Xaxis2->SetMoreLogLabels();

  TGraph *grSigmaVsEnergyU = new TGraphErrors(nbins, Energy, SigmaU);  
  TGraph *grSigmaVsEnergyD = new TGraphErrors(nbins, Energy, SigmaD);  

  grSigmaVsEnergy->Draw("AP");
  grSigmaVsEnergyU->Draw("Lsame");
  grSigmaVsEnergyD->Draw("Lsame");
  cout<<endl<<endl;

  Double_t ErrU = grSigmaVsEnergyU->Eval(5.0) - grSigmaVsEnergy->Eval(5.0);
  Double_t ErrD = -grSigmaVsEnergyD->Eval(5.0) + grSigmaVsEnergy->Eval(5.0);

  cout<<" Total bbbar cross-section at five TeV is : "<<grSigmaVsEnergy->Eval(5.0)<<" mub "<<" pm "<<ErrU<<"  "<<ErrD<<endl;
  cout<<endl<<endl;

  /*
  grSigmaVsEnergy->Draw("AP");
  cout<<endl<<endl;
  cout<<" Total bbbar cross-section at five TeV is : "<<grSigmaVsEnergy->Eval(5.0)<<" mub "<<endl;
  cout<<endl<<endl;
  */
  }






Double_t tsallis_fitting_function(Double_t* x, Double_t* par)
{
  double pT = x[0] ; 

  double dNdy = par[0];
  double nn = par[1];
  double pzero = par[2] ;
  double Mass = par[3] ; 
  double Beta = par[4] ; 
  
  //double Yrap = par[5] ;
  //  double gamma = 1.0/sqrt(1-Beta*Beta);
  
    double MT = sqrt(pT*pT + Mass*Mass) ; 
    
    //  double fac = (1/(2.0*pi))*(nn-1.0)*(nn-2.0)/((nn*Tem + Mass*(nn-1.0))*(nn*Tem+Mass));
    // double Norm = fac*pow((nn*Tem/(nn*Tem+Mass)),-nn);
    // double Norm = pT/MT;
    double Norm = 1.0;
    
    double tsallis = dNdy * Norm * pT *pow((TMath::Exp(-Beta*pT) + MT/(pzero)), -nn) ;

  //   double tsallis = dNdy * Norm * pow((1 + pT/pzero), -nn) ;

  return tsallis ;
   
}


Double_t JPsi_fitting_function(Double_t* x, Double_t* par)
{
  double pT = x[0] ; 

  double alpha1 = par[0];
  double alpha2 = par[1];
  double alpha3 = par[2] ;
  double alpha4 = par[3] ; 
  
  double dNdPt = alpha1*TMath::Power(pT,alpha2)/TMath::Power((alpha3 + pT*pT),alpha4);
  
  return dNdPt;
   
}


TGraphErrors *grf_GhostGraph(Int_t MarkerStyle, Int_t MarkerColor, Int_t LineColor)
{
  const int NN =1;
  Double_t XX[NN]={0.0};
  Double_t Err_XX[NN]={0.0};

  Double_t YY[NN]={0.0};
  Double_t Err_YY[NN]={0.0};


  TGraphErrors *GhostGraph = new TGraphErrors(NN,XX,YY,Err_XX,Err_YY);
  
  GhostGraph->SetMarkerStyle(MarkerStyle);
  GhostGraph->SetMarkerColor(MarkerColor);
  GhostGraph->SetLineColor(LineColor);

  return GhostGraph;


}


//===============================================================================//
//============================ Form Rate Function for plots =====================//
//===============================================================================//

Double_t FormRateMC_P(Double_t Temp, Double_t Pt) 
{

  TH1D *Hist = new TH1D("HistPt","HistPt",50,0.0,25.0);
  Double_t PMin=0.0001;
  Double_t PMax=6.0;
  
  Double_t ThetaMin=0.0;
  Double_t ThetaMax=TMath::Pi();
  Double_t PhiMin=0.0;
  Double_t PhiMax=2.0*TMath::Pi();
  
  int nTrials = 500000;
  
  Double_t Sum=0;
  Double_t Sum1=0;

  for(int i =1; i<= nTrials; i++) {
    
    Double_t RanNo1, RanNo2;
    Double_t RanNo3, RanNo4, RanNo5, RanNo6;	
    
    RanNo1=gRandom->Rndm();
    RanNo2=gRandom->Rndm();
    RanNo3=gRandom->Rndm();
    RanNo4=gRandom->Rndm();
    RanNo5=gRandom->Rndm();
    RanNo6=gRandom->Rndm();
    
    Double_t p1=PMin+(PMax-PMin)*RanNo1;
    Double_t p2=PMin+(PMax-PMin)*RanNo2;
    Double_t theta1= ThetaMin + (ThetaMax -ThetaMin )*RanNo3;
    Double_t theta2= ThetaMin + (ThetaMax -ThetaMin )*RanNo4;
    Double_t phi1= PhiMin + (PhiMax -PhiMin )*RanNo5;
    Double_t phi2= PhiMin + (PhiMax -PhiMin )*RanNo6;
      
    Double_t Px1,Py1,Pz1,Px2,Py2,Pz2;
	
    Px1=p1*sin(theta1)*cos(phi1);
    Py1=p1*sin(theta1)*sin(phi1);
    Pz1=p1*cos(theta1);
    TVector3 PC;
    PC.SetXYZ(Px1,Py1,Pz1);
	
    Px2=p2*sin(theta2)*cos(phi2);
    Py2=p2*sin(theta2)*sin(phi2);
    Pz2=p2*cos(theta2);
    TVector3 PCbar;
    PCbar.SetXYZ(Px2,Py2,Pz2);
    
    TVector3 PCCbar = PC + PCbar;
    Double_t PJpsi=PCCbar.Mag();

    Double_t funNum = FormFun(p1,p2,theta1,theta2,phi1,phi2)*p1*p1*p2*p2*sin(theta1)*sin(theta2)*fcharm(p1,Temp)*fcharm(p2,Temp);
    Double_t funDeno = p1*p1*p2*p2*sin(theta1)*sin(theta2)*fcharm(p1,Temp)*fcharm(p2,Temp);
    
    Sum = Sum + funNum;
    Sum1 = Sum1 + funDeno;
    
    Hist->Fill(PJpsi,funNum);
  }
  
  Double_t step = ((PMax- PMin)*(PMax- PMin)*(ThetaMax-ThetaMin)*(ThetaMax-ThetaMin)*(PhiMax-PhiMin)*(PhiMax-PhiMin))/nTrials;
  
  //Double_t ForRate  = step*Sum;
  Double_t ForRateDeno  = step*Sum1; 
  
  Hist->Scale(hbarc2*step/(ForRateDeno*Hist->GetBinWidth(0))); 
  Double_t ForRateAtPt = Hist->GetBinContent(Hist->FindBin(Pt));
  Hist->Delete();
  return ForRateAtPt;
  
}




Double_t ErrorReCalculate(Double_t OldNumber, Double_t OldError, Double_t NewNumber)
{

  Double_t OldErrorInPercent = OldError/OldNumber;

  Double_t NewError = NewNumber*OldErrorInPercent;

  return NewError;




} 




Double_t ErrorReCalculateExp(Double_t NewNumber, Double_t ErrorPercent )
{

  Double_t NewError = NewNumber*ErrorPercent;

  return NewError;




} 
