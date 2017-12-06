#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
using namespace std;



/// The large arrays that were here are now GONE.
/// Instead, we have this helper that holds the
/// information of all our histograms.

class HistoFactory{
 public:
  std::vector<std::string> vars;
  std::vector<int> nBins;
  std::vector<double> minBin;
  std::vector<double> maxBin;
  void setHisto(std::string s, int n, double min, double max) {
    vars.push_back(s);
    nBins.push_back(n);
    minBin.push_back(min);
    maxBin.push_back(max);
  }
};

/// EDBRHistoMaker is the class that analyzes the flat
/// TTree that comes out from the NTuple dumper module.
/// It doesn't make analysis; it just makes plots.
/// There are a few switches to do different plots:
/// SIGNAL - BACKGROUND,
/// MUON - ELECTRON, etc...

class EDBRHistoMaker {
 public:
  EDBRHistoMaker(TTree *tree=0,
         TFile *fileTMP=0,
                 TH1F* hR1=0,
		 bool wantElectrons=true,
		 bool wantMuons=true,
		 bool wantSideband=true, 
		 bool wantSignal=true,
		 bool wantFullRange=false,
		 int  wantNXJets=1,
		 bool isZZchannel=1);
  virtual ~EDBRHistoMaker();

  /// This is the tree structure. This comes directly from MakeClass
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   // Declaration of leaf types
  
   Double_t        scalef;
   Int_t              nVtx;
   Double_t        theWeight;
   Double_t        lumiWeight;
   Double_t        pileupWeight;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu2;
   Int_t           HLT_Ele2;
   Double_t        nump;
   Double_t        numm;
   Double_t        npT;
   Int_t              lep;
   Double_t        ptVlepJEC;
   Double_t        mtVlepJECnew;
   Double_t        ptlep1;
   Double_t        etalep1;
   Int_t           nlooseeles;
   Int_t           nloosemus;
   Double_t        MET_et;
   Double_t        photonet;
   Double_t        photoneta;
   Double_t        photonsieie;
   Int_t              isprompt;
   Double_t        jet1pt;
   Double_t        jet1eta;
   Double_t        jet2pt;
   Double_t        jet2eta;
   Double_t        Mjj;
   Double_t        Mla;
   Double_t        zepp;
   Double_t        deltaeta;
   Double_t        drla;
   Double_t        drj1a;
   Double_t        drj2a;
   Double_t        drj1l;
   Double_t        drj2l;
   Double_t        j1metPhi;
   Double_t        j2metPhi;
   Double_t        Dphiwajj;
   Double_t        jet1icsv;
   Double_t        jet2icsv;

 


    TFile*  fileTMP_;
    TH1F* hR1_;
    
   // List of branches
   TBranch        *b_scalef;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_theWeight;   //!
   TBranch        *b_lumiWeight;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu2;   //!
   TBranch        *b_HLT_Ele2;
   TBranch        *b_nump;   //!
   TBranch        *b_numm;   //!
   TBranch        *b_npT;   //!
   TBranch        *b_lep;   //!
   TBranch        *b_ptVlepJEC;   //!
   TBranch        *b_mtVlepJECnew;   //!
   TBranch        *b_ptlep1;   //!
   TBranch        *b_etalep1;   //!
   TBranch        *b_nlooseeles;   //!
   TBranch        *b_nloosemus;   //!
   TBranch        *b_MET_et;   //!
   TBranch        *b_photonet;   //!
   TBranch        *b_photoneta;   //!
   TBranch        *b_photonsieie;   //!!
   TBranch        *b_isprompt;   //!
   TBranch        *b_jet1pt;   //!
   TBranch        *b_jet1eta;   //!
   TBranch        *b_jet2pt;   //!
   TBranch        *b_jet2eta;   //!
   TBranch        *b_Mjj;   //!
   TBranch        *b_Mla;   //!
   TBranch        *b_zepp;   //!
   TBranch        *b_deltaeta;   //!
   TBranch        *b_drla;
   TBranch        *b_drj1a;
   TBranch        *b_drj2a;
   TBranch        *b_drj1l;
   TBranch        *b_drj2l;
   TBranch        *b_j1metPhi;
   TBranch        *b_j2metPhi;
   TBranch        *b_Dphiwajj;
   TBranch        *b_jet1icsv;
   TBranch        *b_jet2icsv;


    
  // Basic functions directly from MakeClass
  Int_t    GetEntry(Long64_t entry);
  Long64_t LoadTree(Long64_t entry);
  void     Init(TTree *tree);
  void     Loop(std::string outFileName);

  // Our added functions
  void createAllHistos();
  void printAllHistos();
  void saveAllHistos(std::string outFileName);

  void setWantElectrons(bool doele=false){wantElectrons_=doele;}
  void setWantMuons(bool domu=false){wantMuons_=domu;}
  void setWantSideband(bool dosb=false){wantSideband_=dosb;}
  void setWantSignal(bool dosig=false){wantSignal_=dosig;}
  void setWantNXJets(int nxj=1){wantNXJets_=nxj;}
  void setUnitaryWeights(bool setuniw=false){setUnitaryWeights_=setuniw;}


  int check ( double pt, vector<double> * ptZ  )
  {
    int goodw=1;
    for(unsigned int i =0; i< ptZ->size(); i++)
      {   
	//printf("Comparing %g and %g\n",pt,ptZ->at(i));
	if(pt==ptZ->at(i)) { goodw=0; break;}
	//else {printf("I think they're different\n");}
      }   

    return goodw;
  }

  // Our added variables
  int nVars;
  bool wantElectrons_;
  bool wantMuons_;
  bool wantSideband_;
  bool wantSignal_;
  bool wantFullRange_;
  bool setUnitaryWeights_;
  bool debug_;
  int wantNXJets_;

  double sidebandVHMassLow_;
  double sidebandVHMassHigh_;
  double signalVHMassLow_;
  double signalVHMassHigh_;
  bool isZZchannel_;

  // The histograms
  HistoFactory hs;
  std::map<std::string,TH1D*> theHistograms;
  TH2D *hmjmzz; 
  TH1D *hmzzNEW;
};

void EDBRHistoMaker::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
   fChain->SetBranchAddress("scalef", &scalef, &b_scalef);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("lumiWeight", &lumiWeight, &b_lumiWeight);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu2", &HLT_Mu2, &b_HLT_Mu2);
   fChain->SetBranchAddress("HLT_Ele2",&HLT_Ele2,&b_HLT_Ele2);
   fChain->SetBranchAddress("nump", &nump, &b_nump);
   fChain->SetBranchAddress("numm", &numm, &b_numm);
   fChain->SetBranchAddress("npT", &npT, &b_npT);
   fChain->SetBranchAddress("lep", &lep, &b_lep);
   fChain->SetBranchAddress("ptVlepJEC", &ptVlepJEC, &b_ptVlepJEC);
   fChain->SetBranchAddress("mtVlepJECnew", &mtVlepJECnew, &b_mtVlepJECnew);
   fChain->SetBranchAddress("ptlep1", &ptlep1, &b_ptlep1);
   fChain->SetBranchAddress("etalep1", &etalep1, &b_etalep1);
   fChain->SetBranchAddress("nlooseeles", &nlooseeles, &b_nlooseeles);
   fChain->SetBranchAddress("nloosemus", &nloosemus, &b_nloosemus);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("photonet", &photonet, &b_photonet);
   fChain->SetBranchAddress("photoneta", &photoneta, &b_photoneta);
   fChain->SetBranchAddress("isprompt", &isprompt, &b_isprompt);
   fChain->SetBranchAddress("jet1pt", &jet1pt, &b_jet1pt);
   fChain->SetBranchAddress("jet1eta", &jet1eta, &b_jet1eta);
   fChain->SetBranchAddress("jet2pt", &jet2pt, &b_jet2pt);
   fChain->SetBranchAddress("jet2eta", &jet2eta, &b_jet2eta);
   fChain->SetBranchAddress("Mjj", &Mjj, &b_Mjj);
   fChain->SetBranchAddress("Mla", &Mla, &b_Mla);
   fChain->SetBranchAddress("zepp", &zepp, &b_zepp);
   fChain->SetBranchAddress("deltaeta", &deltaeta, &b_deltaeta);
   fChain->SetBranchAddress("drla",&drla,&b_drla);
   fChain->SetBranchAddress("drj1a",&drj1a,&b_drj1a);
   fChain->SetBranchAddress("drj2a",&drj2a,&b_drj2a);
   fChain->SetBranchAddress("drj1l",&drj1l,&b_drj1l);
   fChain->SetBranchAddress("drj2l",&drj2l,&b_drj2l);
   fChain->SetBranchAddress("j1metPhi",&j1metPhi,&b_j1metPhi);
   fChain->SetBranchAddress("j2metPhi",&j2metPhi,&b_j2metPhi);
   fChain->SetBranchAddress("Dphiwajj",&Dphiwajj,&b_Dphiwajj);
   fChain->SetBranchAddress("jet1icsv",&jet1icsv,&b_jet1icsv);
   fChain->SetBranchAddress("jet2icsv",&jet2icsv,&b_jet2icsv);

//   fChain->SetBranchAddress("photonsieie", &photonsieie, &b_photonsieie);
//   fChain->SetBranchAddress("photonchiso", &photonchiso, &b_photonchiso);
//   fChain->SetBranchAddress("HLT_Ele1", &HLT_Ele1, &b_HLT_Ele1);
//   fChain->SetBranchAddress("HLT_Mu1", &HLT_Mu1, &b_HLT_Mu1);
//   fChain->SetBranchAddress("photonphoiso", &photonphoiso, &b_photonphoiso);
//   fChain->SetBranchAddress("photonnhiso", &photonnhiso, &b_photonnhiso);
 

}

EDBRHistoMaker::EDBRHistoMaker(TTree* tree,
                               TFile* fileTMP,
                               TH1F* hR1,
			       bool wantElectrons,
			       bool wantMuons,
			       bool wantSideband,
			       bool wantSignal,
			       bool wantFullRange,
			       int  wantNXJets,
			       bool isZZchannel){
  fChain = 0;

  // Which category do we want to analyze?
  wantElectrons_ = wantElectrons;
  wantMuons_ = wantMuons;
  wantSideband_ = wantSideband;
  wantSignal_ = wantSignal;
  wantFullRange_ = wantFullRange;
  wantNXJets_ = wantNXJets;
  isZZchannel_ = isZZchannel;
  fileTMP_ = fileTMP;
  hR1_ = hR1;
    

  debug_ = true;
  Init(tree);
  createAllHistos();
  printAllHistos();
}

EDBRHistoMaker::~EDBRHistoMaker() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}                  

Int_t EDBRHistoMaker::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t EDBRHistoMaker::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  return centry;
}

//-------------------------
// Infrastructure functions
//-------------------------

void EDBRHistoMaker::createAllHistos() {

  /// This part substitutes the big arrays that used to be 
  /// in the beginning of this file.
  /// Much simpler to create histos now: just add them to
  /// hs with hs.setHisto(name,nbins,min,max);
    hs.setHisto("nVtx", 20, 0., 50.);
    hs.setHisto("ptlep1", 20, 25,200);
    hs.setHisto("etalep1", 20,-2.1,2.1);
    hs.setHisto("mtVlepJECnew", 20,30,200);
    hs.setHisto("ptVlepJEC", 20,-10,1000);
    hs.setHisto("photonet", 20,20,150);
    hs.setHisto("photoneta", 20,-1.442,1.442);
    hs.setHisto("jet1pt", 20,40,400.);  
    hs.setHisto("jet1eta", 20,-5.,5.);
    hs.setHisto("jet2pt", 20,30,400.);  
    hs.setHisto("jet2eta", 20,-5.,5.);
    hs.setHisto("Mjj", 20,200,400);
    hs.setHisto("zepp", 20,0.,10.);
    hs.setHisto("deltaeta", 20,-10.,10.);
    hs.setHisto("MET_et",20,35,200);
    hs.setHisto("Dphiwajj",20,0,3.2);
    hs.setHisto("Mla",20,-10,200);
    hs.setHisto("drla",20,0.5,10);
    hs.setHisto("drj1a",20,0.5,10);
    hs.setHisto("drj2a",20,0.5,10);
    hs.setHisto("drj1l",20,0.5,10);
    hs.setHisto("drj2l",20,0.5,10);
    hs.setHisto("j1metPhi",20,0.4,3.2);
    hs.setHisto("j2metPhi",20,0.4,3.2);
//    hs.setHisto("photonsieie", 20,0,0.02);
//    hs.setHisto("photonphoiso", 40,0,10.);
//    hs.setHisto("photonchiso", 20,0,5.);
//    hs.setHisto("photonnhiso", 40,0,20.);


  char buffer[256];
  char buffer2[256];

  nVars = hs.vars.size();

  for(int i = 0; i!= nVars; ++i) {
    sprintf(buffer,"%s_mu",hs.vars[i].c_str());
//    sprintf(buffer,"%s_el",hs.vars[i].c_str());
    sprintf(buffer2,"%s;%s;Number of events;",hs.vars[i].c_str(),hs.vars[i].c_str());
    TH1D* histogram = new TH1D(buffer,
			       buffer2,
			       hs.nBins[i],
			       hs.minBin[i],
			       hs.maxBin[i]);
    histogram->SetStats(kFALSE);
    histogram->SetDirectory(0);
    histogram->Sumw2();
    theHistograms[hs.vars[i]] = histogram;
  }

}

void EDBRHistoMaker::printAllHistos() {
  printf("We have %i histograms \n",int(theHistograms.size()));
  typedef std::map<std::string, TH1D*>::iterator it_type;
  for(it_type iterator = theHistograms.begin(); iterator != theHistograms.end(); iterator++) {

    iterator->second->Print();
    // Repeat if you also want to iterate through the second map.
  }
}

void EDBRHistoMaker::saveAllHistos(std::string outFileName) {

  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");

  for(int i = 0; i!=nVars; ++i) {
    std::string name = hs.vars[i];
    const TH1D* thisHisto = this->theHistograms[name];
    thisHisto->Write();
  }
  outFile->Close();
}

///----------------------------------------------------------------
/// This is the important function, the loop over all events.
/// Here we fill the histograms according to cuts, weights,
/// and can also filter out events on an individual basis.
///----------------------------------------------------------------

void EDBRHistoMaker::Loop(std::string outFileName){

 // TFile * input1 = new TFile ("test.root");
//  TH1F* hR1= (TH1F*)input1->Get("hRatio");
  if (fChain == 0) {std::cout<<outFileName<<" =filename "<<std::endl; return;}
else{
 // double lumivalue = 0.57831;
    int numbe_out = 0;
    double actualWeight;
  Long64_t nentries = fChain->GetEntriesFast();
    std::cout<<"nentries"<<nentries<<std::endl;
  Long64_t npp = fChain->GetEntries("theWeight>0.");
  Long64_t nmm = fChain->GetEntries("theWeight<0.");
  Double_t nn;
  std::cout<< "numberofnp:" << npp << "  numberofnm:" <<nmm << std::endl;
  Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //    for (Long64_t jentry=0; jentry<10000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%1000000==0) std::cout << "Entry num " << jentry << std::endl;

 
   //  int bin = hR1_->FindBin(nVtx);
   // pileupWeight = hR1_->GetBinContent(bin);
 

    if(theWeight>0) nn=1;
    else nn= -1;

   // if(npp>0) actualWeight = lumiWeight*pileupWeight/(npp-nmm)*nn*lumivalue;
   // else actualWeight = lumiWeight*pileupWeight/nentries*lumivalue;//0.00810255;//0.0521;//0.04024;//0.00559;
      actualWeight = lumiWeight*pileupWeight*scalef;
      //      actualWeight =scalef;
    if(setUnitaryWeights_) {
      if(jentry==0)printf("Unitary weights set!\n");
      actualWeight=1.0;
    }


    // We get the histogram from the map by string and fill it.
    // We could wrap all the fills in the this->eventPassesCut()
    // to fill histograms only for the events which pass a given
    // cut (Sideband / SignalRegion, Muon / Electron, 
    // Single / Double jet ...) 

    // Remember: bool eventPassesCut(double ptZ_threshold, double ptlep1_threshold );
    //if(eventPassesCut(20, 20)){


   //   int bveto=1;
     // for(int i=0; i<8; i++)  {

	//if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.890  && fabs(ak4jet_eta[i])<2.4  && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8) {bveto=0; break;}

 //     }

	//Int_t nLooseLep=nLooseEle+nLooseMu;
///        if(lep==13 && mtVlepJECnew>50. && ptlep1>30. && fabs(etalep1)<2.1 && MET_et>30. && nlooseeles==0 && nloosemus<2);
//        if(lep==11  && photonet>15.  && photonet<45. && HLT_Ele1>0 && fabs(photoneta)<1.4442 && ptlep1>30. && fabs(etalep1)<2.1 && ptlep2>30. && fabs(etalep2)<2.1 && nlooseeles<3 && nloosemus==0);
           TString filename = fileTMP_->GetName();

	if(lep==13 && nlooseeles==0 && nloosemus<2 && HLT_Mu3 ==1 && mtVlepJECnew>30 && ptlep1>25 && fabs(etalep1)<2.1 && photonet>20 && photonet<150 && fabs(photoneta)<1.4442 && drla>0.5 && jet1pt>40 && jet2pt>30 && fabs(jet1eta)<4.7 && fabs(jet2eta)<4.7 && Mjj>200 && Mjj<400 && drj1a>0.5 && drj2a>0.5 && drj1l>0.5 && drj2l>0.5 && j1metPhi>0.4 && j2metPhi>0.4  && MET_et>35);  //  && zepp <1.3 && fabs(deltaeta)>2.4   && jet1icsv<0.8484 && jet2icsv<0.8484 && ptVlepJEC>0
                  else continue;


            int iswjets = 0;
            int isnotwets = 0;
            int iszjets = 0;
            int isttjets = 0;

            if (filename.Contains("WJets") && isprompt!= 1){
                iswjets = 1;
             //   std::cout<<"WJETs"<<iswjets<<std::endl;
                }
            if (filename.Contains("DY") && isprompt!= 1){
                iszjets = 1;
                }
            if (filename.Contains("TTJets") && isprompt!= 1){
                isttjets = 1;
                }
            if (!(filename.Contains("WJets")) && !(filename.Contains("DY")) && !(filename.Contains("TTJets")) ) {
                isnotwets = 1;
                }
 

 
            if (isnotwets>0 || iswjets>0 || iszjets>0 || isttjets>0) {


if(Mjj>=4000)Mjj=3999;
if(photonet>=400)photonet=399;
if(nVtx>=50)nVtx=49;
if(ptlep1>=200)ptlep1=199;
if(mtVlepJECnew>=200)mtVlepJECnew=199;
if(ptVlepJEC>=300)ptVlepJEC=299;
if(jet1pt>=400)jet1pt=399;
if(jet2pt>=300)jet2pt=299;
if(MET_et>=200)MET_et=199;
if(Mla>=200)Mla=199;
      (theHistograms["nVtx"])->Fill(nVtx,actualWeight);
      (theHistograms["ptlep1"])->Fill(ptlep1,actualWeight);
      (theHistograms["etalep1"])->Fill(etalep1,actualWeight);
      (theHistograms["mtVlepJECnew"])->Fill(mtVlepJECnew,actualWeight);
      (theHistograms["ptVlepJEC"])->Fill(ptVlepJEC,actualWeight);
      (theHistograms["photonet"])->Fill(photonet,actualWeight);
      (theHistograms["photoneta"])->Fill(photoneta,actualWeight);
      (theHistograms["jet1pt"])->Fill(jet1pt,actualWeight);
      (theHistograms["jet1eta"])->Fill(jet1eta,actualWeight);
      (theHistograms["jet2pt"])->Fill(jet2pt,actualWeight);
      (theHistograms["jet2eta"])->Fill(jet2eta,actualWeight);
      (theHistograms["Mjj"])->Fill(Mjj,actualWeight);
      (theHistograms["zepp"])->Fill(zepp,actualWeight);
      (theHistograms["deltaeta"])->Fill(deltaeta,actualWeight);
      (theHistograms["MET_et"])->Fill(MET_et,actualWeight);
      (theHistograms["Dphiwajj"])->Fill(Dphiwajj,actualWeight);
      (theHistograms["Mla"])->Fill(Mla,actualWeight);
      (theHistograms["drla"])->Fill(drla,actualWeight);
      (theHistograms["drj1a"])->Fill(drj1a,actualWeight);
      (theHistograms["drj2a"])->Fill(drj2a,actualWeight);
      (theHistograms["drj1l"])->Fill(drj1l,actualWeight);
      (theHistograms["drj2l"])->Fill(drj2l,actualWeight);
      (theHistograms["j1metPhi"])->Fill(j1metPhi,actualWeight);
      (theHistograms["j2metPhi"])->Fill(j2metPhi,actualWeight);
//      (theHistograms["photonsieie"])->Fill(photonsieie,actualWeight);
//      (theHistograms["photonphoiso"])->Fill(photonphoiso,actualWeight);
//      (theHistograms["photonchiso"])->Fill(photonchiso,actualWeight);
//      (theHistograms["photonnhiso"])->Fill(photonnhiso,actualWeight);


 
    //  (theHistograms["philep1"])->Fill(philep1,actualWeight);
    //  (theHistograms["MET_et"])->Fill(MET_et,actualWeight);
     // }//end if eventPassesCut
                }
  }//end loop over entries
   cout << "after cut: " << numbe_out << "*actualweight" << actualWeight << " result " << numbe_out*actualWeight << endl; 
  //std::cout<<"From makeHisto: the histo with #vtx has "<<(theHistograms["nVtx"])->GetEntries()<<" entries"<<std::endl;
  this->saveAllHistos(outFileName);
}
}
