#include <string>
#include <vector>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
//#include "CMSLabels.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLatex.h" 
#include "CMS_lumi.C"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
 

class EDBRHistoPlotter {
public:

  EDBRHistoPlotter(std::string nameInDir,
                   std::vector<std::string> nameFileDATA,
                   std::vector<std::string> nameFileMC,
                   std::vector<std::string> nameFileMCSig,
                   double targetLumi,
                   int wantNXJets,
                   int flavour,
                   bool isZZchannel,
                   bool scaleToData,
                   bool scaleOnlyWJets,
                   bool makeRatio,
                   bool isSignalStackOnBkg,
                   std::vector<double> kFactorsMC,
                   std::vector<double> kFactorsMCSig)
  {
    std::cout << "Plotter constructor" << std::endl;
    nameInDir_      = nameInDir;
    fileNamesMC     = nameFileMC;
    fileNamesMCSig  = nameFileMCSig;
    fileNamesDATA   = nameFileDATA;
    kFactorsMC_     = kFactorsMC;
    kFactorsSig_    = kFactorsMCSig;
    targetLumi_     = targetLumi;
    wantNXJets_     = wantNXJets;
    flavour_        = flavour;
    isZZchannel_    = isZZchannel;
    scaleToData_    = scaleToData;
    scaleOnlyWJets_ = scaleOnlyWJets;
    makeRatio_      = makeRatio;
    isSignalStackOnBkg_ = isSignalStackOnBkg;
    debug_          = true;
    if (fileNamesDATA.size() != 0)
      isDataPresent_ = true;
    else
      isDataPresent_ = false;
    std::cout << "Check" << std::endl;
    EDBRColors.resize(20, kBlack);
    EDBRLineColors.resize(20, kBlack);
    std::cout << "Check" << std::endl;
    labels.resize(0);
    labelsSig.resize(0);
    std::cout << "Check" << std::endl;
    makeLabels();
    std::cout << "Check" << std::endl;
   if (fileNamesMCSig.size() != kFactorsSig_.size()) {
      cout << "======> Mismatch in size of input MC Sig arrays !!! " << fileNamesMCSig.size() << "  " << kFactorsSig_.size() << endl;
    }

    printf("Target lumi is %g fb-1\n", targetLumi);
    std::cout << "k factors for MC backgrounds are: " << std::endl;
    int myKindex = 0;
    for (std::vector<double>::iterator it = kFactorsMC_.begin(); it != kFactorsMC_.end(); ++it) {
      std::cout << *it << " for " << fileNamesMC.at(myKindex) << std::endl ;
      myKindex++;
    }
    myKindex = 0;
    for (std::vector<double>::iterator it = kFactorsSig_.begin(); it != kFactorsSig_.end(); ++it) {
      std::cout << *it << " for " << fileNamesMCSig.at(myKindex) << std::endl ;
      myKindex++;
    }

    std::cout << std::endl;
  }//end constructor

  virtual ~EDBRHistoPlotter()
  {
  }//end destructor

  /// Members
  std::vector<std::string> fileNamesMC;
  std::vector<std::string> fileNamesMCSig;
  std::vector<std::string> fileNamesDATA;
  std::vector<std::string> labels;
  std::vector<std::string> labelsSig;
  std::vector<TFile*>      filesMC;
  std::vector<TFile*>      filesMCSig;
  std::vector<TFile*>      filesDATA;
  std::vector<TH1D*>       histosMC;
  std::vector<TH1D*>       histosMCSig;
  std::vector<TH1D*>       histosMCSigOrig;
  std::vector<TH1D*>       histosDATA;
  std::vector<int>         EDBRColors;
  std::vector<int>         EDBRLineColors;
  std::vector<double>      kFactorsMC_;
  std::vector<double>      kFactorsSig_;

  std::string nameInDir_;
  std::string nameOutDir_;
  double dataIntegral_;
  double targetLumi_;
  int    wantNXJets_;
  int    flavour_;
  bool   isZZchannel_;
  bool   scaleToData_;
  bool   scaleOnlyWJets_;
  bool   makeRatio_;
  bool   isSignalStackOnBkg_;
  bool   isDataPresent_;
  bool   debug_;

  /// Functions
  void cleanupMC();
  void cleanupMCSig();
  void cleanupDATA();
  void makeLabels();
  void makeStackPlots(std::string histoName);
  void setOutDir(std::string outDirNew);

  /// set debug mode
  void setDebug(bool debug)
  {
    debug_ = debug;
  }

  /// get reasonable colors for stacks.
  int getFillColor(int index)
  {
    if (index < 20)return EDBRColors[index];
    return kWhite;
  }

  /// set reasonable colors for stacks.
  void setFillColor(std::vector<int> colorList)
  {
    unsigned int ind = 0;
    while (ind < 20 && ind < colorList.size()) {
      EDBRColors.at(ind) = colorList.at(ind);
      ind++;
    }
  }

  /// get reasonable colors for stacks.
  int getLineColor(int index)
  {
    if (index < 20)return EDBRLineColors[index];
    return kBlack;
  }

  /// set reasonable colors for stacks.
  void setLineColor(std::vector<int> colorList)
  {
    unsigned int ind = 0;
    while (ind < 20 && ind < colorList.size()) {
      EDBRLineColors.at(ind) = colorList.at(ind);
      ind++;
    }
  }

};

void EDBRHistoPlotter::cleanupMC()
{
  for (size_t i = 0; i != filesMC.size(); ++i) {
    filesMC.at(i)->Close();
  }
  filesMC.clear();

  for (size_t i = 0; i != histosMC.size(); ++i) {
    histosMC.at(i)->Delete();
  }
  histosMC.clear();
}

void EDBRHistoPlotter::cleanupMCSig()
{
  for (size_t i = 0; i != filesMCSig.size(); ++i) {
    filesMCSig.at(i)->Close();
  }
  filesMCSig.clear();

  for (size_t i = 0; i != histosMCSig.size(); ++i) {
    histosMCSig.at(i)->Delete();
    histosMCSigOrig.at(i)->Delete();
  }
  histosMCSig.clear();
  histosMCSigOrig.clear();
}

void EDBRHistoPlotter::cleanupDATA()
{
  for (size_t i = 0; i != filesDATA.size(); ++i) {
    filesDATA.at(i)->Close();
  }
  filesDATA.clear();

  for (size_t i = 0; i != histosDATA.size(); ++i) {
    histosDATA.at(i)->Delete();
  }
  histosDATA.clear();
}

void EDBRHistoPlotter::makeLabels()
{
  for (size_t i = 0; i != fileNamesMC.size(); i++) {
    TString s1 = fileNamesMC.at(i);
    TString s2 = "_.";
    TObjArray* tokens = s1.Tokenize(s2);
    std::string aLabel  = ((TObjString*)(tokens->At(1)))->String().Data();
    std::string aLabel2 = ((TObjString*)(tokens->At(2)))->String().Data();
    labels.push_back(aLabel);
//    labels.push_back(aLabel + aLabel2);
  }
  std::cout << "Labels MC done" << std::endl;

  for (size_t i = 0; i != fileNamesMCSig.size(); i++) {
    TString s1 = fileNamesMCSig.at(i);
    TString s2 = "_.";
    TObjArray* tokens = s1.Tokenize(s2);

    std::string aLabelType = ((TObjString*)(tokens->At(1)))->String().Data();

    std::string aLabelCoupling = ((TObjString*)(tokens->At(2)))->String().Data();

  //  std::string aLabelMass = ((TObjString*)(tokens->At(3)))->String().Data();
    //  std::cout << "Right3?" << std::endl;

    std::cout << s1.Data();

    std::string aLabel = aLabelType + aLabelCoupling;

    labelsSig.push_back(aLabel);

  }
  std::cout << "Labels MC signal done" << std::endl;

}

///set output directories for plots.
void EDBRHistoPlotter::setOutDir(std::string outDirNew)
{
  char buffer[256];
  nameOutDir_ = outDirNew;

  sprintf(buffer, "%s/pdf", nameOutDir_.c_str());
  printf("%s\n", buffer);
  gSystem->mkdir(buffer, true);

  sprintf(buffer, "%s/root", nameOutDir_.c_str());
  printf("%s\n", buffer);
  gSystem->mkdir(buffer, true);

  sprintf(buffer, "%s/png", nameOutDir_.c_str());
  printf("%s\n", buffer);
  gSystem->mkdir(buffer, true);
}

void EDBRHistoPlotter::makeStackPlots(std::string histoName)
{

  cleanupMC();
  cleanupMCSig();
  cleanupDATA();

  //printf("Making histo %s\n",histoName.c_str());
  std::cout << "\rMaking histo " << histoName.c_str() << std::endl;

  TCanvas* cv = new TCanvas(("cv_" + histoName).c_str(), ("cv_" + histoName).c_str(), 700, 600);

  //create 3 pads in the canvas
  TPad* fPads1 = NULL;
  TPad* fPads2 = NULL;
//  TPad* fPads3 = NULL;

  if (makeRatio_ && isDataPresent_) {
    fPads1 = new TPad("pad1", "", 0.00, 0.25, 0.99, 0.99);
    fPads2 = new TPad("pad2", "", 0.00, 0.00, 0.99, 0.25);
    //fPads2 = new TPad("pad2", "", 0.00, 0.20, 0.99, 0.40);
    //fPads3 = new TPad("pad3", "", 0.00, 0.00, 0.99, 0.20);
    fPads1->SetFillColor(0);
    fPads1->SetLineColor(0);
    fPads2->SetFillColor(0);
    fPads2->SetLineColor(0);
      fPads1->SetBottomMargin(0);
      fPads2->SetTopMargin(0);
      fPads2->SetBottomMargin(0.3);

    //fPads3->SetFillColor(0);
    //fPads3->SetLineColor(0);
    fPads1->Draw();
    fPads2->Draw();
    //fPads3->Draw();
  }


  //============ Data vs MC plots ==============

  if (makeRatio_ && isDataPresent_) {
    fPads1->cd();
  }
 // fPads1->cd();
  ///--------------------
  /// Make the DATA stack
  ///--------------------

  TH1D* sumDATA = NULL;

  for (size_t i = 0; i != fileNamesDATA.size(); ++i) {
    filesDATA.push_back(TFile::Open((nameInDir_ +
                                     fileNamesDATA.at(i)).c_str()));
  }

  for (size_t i = 0; i != filesDATA.size(); ++i) {
    TH1D* histo = (TH1D*)(filesDATA.at(i)->Get(histoName.c_str())->Clone("clone"));
    histo->SetDirectory(0);
    histosDATA.push_back(histo);
  }

  if (histosDATA.size() != 0) {
    sumDATA = (TH1D*)(histosDATA.at(0)->Clone("masterDATA"));
    sumDATA->Reset();
    sumDATA->SetDirectory(0);
  }

  for (size_t i = 0; i != histosDATA.size(); ++i) {
    sumDATA->Add(histosDATA.at(i));
  }

  double sumDataIntegral = 0;
  double sumDataIntegralerr;
  int sumDataIntegraln = sumDATA->GetNbinsX();
  if (isDataPresent_)
    sumDataIntegral = sumDATA->IntegralAndError(1,sumDataIntegraln,sumDataIntegralerr);


  ///------------------
  /// Make the MC stack
  ///------------------

  TH1D* sumMC = NULL;
  TH1D* hdiff = NULL;
  double sumBkgAtTargetLumi = 0;

  for (size_t i = 0; i != fileNamesMC.size(); ++i) {
    filesMC.push_back(TFile::Open((nameInDir_ +
                                   fileNamesMC.at(i)).c_str()));
  }

  double sumBkgOther = 0.;
  double sumBkgOthererr = 0;
  double sumBkgOthererrpre;
  int sumBkgOthern;

  double sumWJets = 0.;
  double sumWJetserr;
  int sumWJetsn;

  double sumWAJJ = 0;


  for (size_t i = 0; i != filesMC.size(); ++i) {
    TH1D* histo = (TH1D*)(filesMC.at(i)->Get(histoName.c_str())->Clone(labels.at(i).c_str()));
    histo->Scale(kFactorsMC_.at(i));

    sumBkgAtTargetLumi += (histo->Integral() * targetLumi_);

    histo->Scale(targetLumi_);
    TString filename = filesMC.at(i)->GetName();
    if(filename.Contains("WAJJ"))sumWAJJ=histo->Integral();
    if (filename.Contains("WA") && !(filename.Contains("WAJJ"))){sumWJetsn=histo->GetNbinsX(); sumWJets = histo->IntegralAndError(1,sumWJetsn,sumWJetserr);}
    else {sumBkgOthern=histo->GetNbinsX();  sumBkgOther += (histo->IntegralAndError(1,sumBkgOthern,sumBkgOthererrpre)); sumBkgOthererr=sqrt(sumBkgOthererr*sumBkgOthererr+sumBkgOthererrpre*sumBkgOthererrpre);}
  }


   double WJetsScaleFactor = 1.;
   double WJetsScaleFactorerr;
   double WAJJratio;

  if (scaleOnlyWJets_) {
    WJetsScaleFactor = (sumDataIntegral - sumBkgOther) / sumWJets;
    WJetsScaleFactorerr = 1/sumWJets*sqrt(sumDataIntegralerr*sumDataIntegralerr+sumBkgOthererr*sumBkgOthererr+sumWJetserr*sumWJetserr*(sumDataIntegral-sumBkgOther)*(sumDataIntegral-sumBkgOther)/sumWJets/sumWJets);
    WAJJratio= sumWAJJ/(sumBkgOther+sumWJets);

    cout << "WAScaleFactor= " << WJetsScaleFactor  <<" +- "<< WJetsScaleFactorerr << endl;
ofstream myfile("WAscalefactor.txt");
myfile << "WA scalefactor in control region: " << WJetsScaleFactor <<" +- "<< WJetsScaleFactorerr << endl;
myfile << "signal ratio in control region: "<< WAJJratio <<endl;
myfile.close();

  }

  for (size_t i = 0; i != filesMC.size(); ++i) {
    TH1D* histo = (TH1D*)(filesMC.at(i)->Get(histoName.c_str())->Clone(labels.at(i).c_str()));
    histo->SetDirectory(0);
    histo->SetFillColor(getFillColor(i));
    histo->SetLineColor(0);

    TString filename = filesMC.at(i)->GetName();
    histo->Scale(kFactorsMC_.at(i));
//    if (filename.Contains("WA") && !(filename.Contains("WAJJ")))histo->Scale(WJetsScaleFactor);    //normalize the WA to data
    histosMC.push_back(histo);
  }



 
  if (histosMC.size() != 0) {
    sumMC = (TH1D*)(histosMC.at(0)->Clone("masterMC"));
    sumMC->Reset();
    sumMC->SetDirectory(0);
    hdiff = (TH1D*)(histosMC.at(0)->Clone("masterMC"));
    hdiff->Reset();
    hdiff->SetDirectory(0);
  }

  /// Do we normalize to data or to lumi?
  /// NOTICE THAT THIS DEPENDS ON THE HISTOGRAMS HAVING BEING
  /// CORRECTLY FILLED WITH PUweight*LumiWeight*GenWeight
  for (size_t is = 0; is != histosMC.size(); is++) {
    if (scaleToData_ && isDataPresent_) {
      histosMC.at(is)->Scale(targetLumi_ *sumDataIntegral / sumBkgAtTargetLumi);
    } else {
      histosMC.at(is)->Scale(targetLumi_);
    }
  }

  THStack* hs = new THStack("hs", "");

  // Make a histogram just for the sum
  for (size_t i = 0; i != histosMC.size(); ++i) {

    histosMC.at(i)->SetFillColor(getFillColor(i));
    histosMC.at(i)->SetLineColor(kBlack);
    histosMC.at(i)->SetLineWidth(1);	
    sumMC->Add(histosMC.at(i));
    hs->Add(histosMC.at(i));
  }


  sumMC->SetFillStyle(0);
//  sumMC->SetLineColor(kBlack);
//  sumMC->SetLineColor(0);


  if (scaleToData_ && isDataPresent_) {
    std::cout << "===> Residual DATA/MC Scale Factor is: " << sumDataIntegral / sumBkgAtTargetLumi << std::endl;
  }

  ///-------------------------------
  /// Add the MC signal to the stack
  ///-------------------------------

  for (size_t i = 0; i != fileNamesMCSig.size(); ++i) {
    filesMCSig.push_back(TFile::Open((nameInDir_ +
                                      fileNamesMCSig.at(i)).c_str()));
  }

  for (size_t i = 0; i != filesMCSig.size(); ++i) {


    TH1D* histo = (TH1D*)(filesMCSig.at(i)->Get(histoName.c_str())->Clone(labelsSig.at(i).c_str()));
    TH1D* histoOrig = (TH1D*)(filesMCSig.at(i)->Get(histoName.c_str())->Clone(labelsSig.at(i).c_str()));
    histo->SetDirectory(0);
    histo->SetLineColor(kBlack);
    histo->SetFillColor(2);
    histo->SetLineWidth(1);


    histoOrig->SetDirectory(0);
    histoOrig->SetLineColor(getLineColor(i));
    histoOrig->SetFillColor(getLineColor(i));
    if (i % 2 == 0)histoOrig->SetFillStyle(3004);
    else histoOrig->SetFillStyle(3005);
    //histo->Scale(kFactor_); //============= SCALE FACTORS FOR SIGNAL? ==== FIXME

     histosMCSig.push_back(histo);
     histosMCSigOrig.push_back(histoOrig);
  }

  //scale the MC signal histogram
  if (histosMCSig.size() != kFactorsSig_.size()) cout << "+++++++++++++++++ Mismatch in size of input MC Sig arrays !!!" << endl;



  for (size_t is = 0; is != histosMCSig.size(); is++) {
 

    histosMCSig.at(is)->Scale(targetLumi_ * kFactorsSig_.at(is));
    histosMCSigOrig.at(is)->Scale(targetLumi_ * kFactorsSig_.at(is));

    hs->Add(histosMCSig.at(is));
    sumMC->Add(histosMCSig.at(is)); 
    //add the signal to the total background
    // histosMCSig.at(is)->Add(sumMC);
  }


  ///-----------------------------------
  /// Draw both MC and DATA in the stack
  ///-----------------------------------

  hs->Draw("HIST");  //
 // hs->GetXaxis()->SetTitle(histoName.c_str());
  hs->GetYaxis()->SetTitle("Events/bin");//40.24pb-1");
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleOffset(0.9);
//  hs->GetYaxis()->CenterTitle();
  
  double maximumMC = 1.5 * sumMC->GetMaximum();
  double maximumDATA = -100;
  if (isDataPresent_)
    maximumDATA = 1.5 * sumDATA->GetMaximum();
  double maximumForStack = -100;
  if (isDataPresent_)
    maximumForStack = (maximumMC > maximumDATA ? maximumMC : maximumDATA);
  else
    maximumForStack = maximumMC;
  hs->SetMaximum(maximumForStack);
  // Some hacks for better aestetics
  // Extra vertical space in eta plots
  hs->SetMinimum(1.0);
    if (isDataPresent_) {
          sumDATA->SetMarkerColor(1);
          sumDATA->SetMarkerStyle(20);
          sumDATA->Draw("SAME EP"); 
 }

   // histosMCSig.at(0)->Draw("SAME HIST");
      
  // For the legend, we have to tokenize the name "histos_XXX.root"
  TLegend* leg = new TLegend(0.75, 0.55, 0.88, 0.85);
  TLegend* leg2 = new TLegend(0.5, 0.6, 0.72, 0.85);  //every 0.055 a draw in y axis
  leg->SetTextSize(0.04);
  leg2->SetTextSize(0.04);

  leg->SetMargin(0.4);
  if (isDataPresent_)
    leg2->AddEntry(sumDATA, "Data", "ep");
  for (size_t i = 0; i != histosMC.size(); ++i){
	if(i<6){leg->AddEntry(histosMC.at(i), labels.at(i).c_str(), "f"); }
	else{leg2->AddEntry(histosMC.at(i), labels.at(i).c_str(), "f");}
						}
 
  if (histosMCSig.size() > 0) {
    char rescalingLabel[64];
    for (size_t i = 0; i != histosMCSig.size(); ++i) {
      sprintf(rescalingLabel, " (x%g)", kFactorsSig_.at(i));
      std::string rescalingStr(rescalingLabel);
      if (kFactorsSig_.at(i) != 2.0){
        if(i==0) leg2->AddEntry(histosMCSig.at(i), "EWK W#gamma+2Jets", "f");
	}
      else leg2->AddEntry(histosMCSig.at(i), (labelsSig.at(i)).c_str(), "f");
    }
  }

  leg->SetFillStyle(0);
  leg2->SetFillStyle(0);
//  leg->SetFillColor(kWhite);

  // Nice labels
 // TMathText* l = makeCMSPreliminaryTop(13, 0.50, 0.935);
  //  TMathText* l = makeCMSLumi(13,551.7,0.6,0.935);
  //l->Draw();
  //============ Data-MC/Error ==============THstack error band=============

  TLine* lineAtZero = NULL;
  TLine* lineAtPlusTwo = NULL;
  TLine* lineAtMinusTwo = NULL;
  if (makeRatio_ && isDataPresent_) {


    double thisYmin = 0.5;  //-5
    double thisYmax = 1.5;  //5 

//    TVectorD nsigma_x(sumDATA->GetNbinsX());
//    TVectorD nsigma_y(sumDATA->GetNbinsX());

int binnumber=sumMC->GetNbinsX();
double nsigma_x[binnumber],nsigma_y[binnumber];
double Data[binnumber],Bkg[binnumber],eData[binnumber],eBkg[binnumber],x[binnumber],halfwidth[binnumber],ex[binnumber],ey[binnumber],sigma[binnumber];


    for (int ibin = 0; ibin != binnumber; ++ibin) 
	{

      Data[ibin] = sumDATA->GetBinContent(ibin + 1);
      Bkg[ibin] = sumMC->GetBinContent(ibin + 1);
      eData[ibin] = sumDATA->GetBinError(ibin + 1);
      eBkg[ibin] = sumMC->GetBinError(ibin + 1);
      x[ibin] = sumDATA->GetBinCenter(ibin + 1);
	halfwidth[ibin]=sumMC->GetBinWidth(ibin + 1)/2;

      double diff = Data[ibin] - Bkg[ibin];
	double difference = Data[ibin]/Bkg[ibin];
      sigma[ibin] = sqrt((eData[ibin] * eData[ibin]) + (eBkg[ibin] * eBkg[ibin]));
	ex[ibin]=0;ey[ibin]=sqrt(eData[ibin]*eData[ibin]/Bkg[ibin]/Bkg[ibin]+Data[ibin]*Data[ibin]*eBkg[ibin]*eBkg[ibin]/Bkg[ibin]/Bkg[ibin]/Bkg[ibin]/Bkg[ibin]);

      if (sigma[ibin] != 0.0 && Data[ibin] != 0.0) {
        nsigma_x[ibin] = x[ibin];
        nsigma_y[ibin] = difference;   //diff / sigma;
      } else {
        nsigma_x[ibin] = +999999;
        nsigma_y[ibin] = 0;
	      }
	 }

    fPads1->cd();
//double hsmax=hs->GetMaximum();  //maximumForStack   min=0.1
//double hsmin=hs->GetMinimum();
double xlow=x[0]-halfwidth[0];  double xhigh=x[binnumber-1]+halfwidth[binnumber-1];
TH2F *hempty=new TH2F("hempty","hempty",binnumber,xlow,xhigh,10000,0.1,maximumForStack);

	TGraphAsymmErrors *gr=new TGraphAsymmErrors(binnumber,x,Bkg,halfwidth,halfwidth,eBkg,eBkg);
	gr->SetFillStyle(3004);
	hempty->Draw("same");
	gr->Draw("2same");
	leg2->AddEntry(gr, "MC unc", "f");
	leg->Draw();
	  leg2->Draw("same");
    fPads1->Update();


    fPads2->cd();

//    fPads2->SetGridx();
    fPads2->SetGridy();

//    if (nsigma_x.GetNoElements() != 0) {
      //TGraph *nsigmaGraph = new TGraph(nsigma_x, nsigma_y);
TGraphErrors *nsigmaGraph = new TGraphErrors(binnumber, nsigma_x, nsigma_y, ex, ey);
      nsigmaGraph->SetTitle("");
      nsigmaGraph->GetYaxis()->SetRangeUser(thisYmin, thisYmax);
      nsigmaGraph->GetYaxis()->SetTitle("Data/MC");     //("(Data-Bkg)/#sigma");
      nsigmaGraph->GetYaxis()->CenterTitle();
      nsigmaGraph->GetYaxis()->SetTitleOffset(0.27);
      nsigmaGraph->GetYaxis()->SetTitleSize(0.15);
      nsigmaGraph->GetYaxis()->SetLabelSize(0.06);
      nsigmaGraph->GetXaxis()->SetTitle(histoName.c_str());
        nsigmaGraph->GetXaxis()->SetTitleSize(0.1);
      nsigmaGraph->GetXaxis()->SetLimits(sumMC->GetXaxis()->GetXmin() , sumMC->GetXaxis()->GetXmax());
      nsigmaGraph->GetXaxis()->SetRangeUser(sumMC->GetXaxis()->GetXmin() , sumMC->GetXaxis()->GetXmax());
      nsigmaGraph->GetXaxis()->SetTitleOffset(0.9);
      nsigmaGraph->GetXaxis()->SetLabelSize(0.08);
      nsigmaGraph->SetMarkerStyle(8);
      nsigmaGraph->SetMarkerSize(1.);
      nsigmaGraph->Draw("ape");
//    }

    fPads2->Update();

//    lineAtZero = new TLine(sumMC->GetXaxis()->GetXmin(), 0, sumMC->GetXaxis()->GetXmax(), 0);
//    lineAtZero->SetLineColor(2);
//    lineAtZero->Draw();
    lineAtPlusTwo = new TLine(sumMC->GetXaxis()->GetXmin(), 1, sumMC->GetXaxis()->GetXmax(), 1);
    lineAtPlusTwo->SetLineColor(2);
    lineAtPlusTwo->SetLineStyle(2);
    lineAtPlusTwo->SetLineWidth(2);
    lineAtPlusTwo->Draw();
//    lineAtMinusTwo = new TLine(sumMC->GetXaxis()->GetXmin(), -2, sumMC->GetXaxis()->GetXmax(), -2);
//    lineAtMinusTwo->SetLineColor(2);
//    lineAtMinusTwo->SetLineStyle(2);
//    lineAtMinusTwo->Draw();
  }

//  char bufferstr[20];
//  sprintf(bufferstr, "%.3f",targetLumi_);
//  CMS_lumi(fPads1,4,0, bufferstr);
  CMS_lumi(fPads1,4,0, "35.861");

  // Save the picture
  char buffer[256];
  cv->SetLogy(false);
  sprintf(buffer, "%s/png/can_%s.png", nameOutDir_.c_str(), histoName.c_str());
  cv->SaveAs(buffer);
  sprintf(buffer, "%s/root/can_%s.root", nameOutDir_.c_str(), histoName.c_str());
  cv->SaveAs(buffer);
  sprintf(buffer, "%s/pdf/can_%s.pdf", nameOutDir_.c_str(), histoName.c_str());
  cv->SaveAs(buffer);
  if (makeRatio_ && isDataPresent_) {
    fPads1->cd();
    fPads1->SetLogy(true);
  } else {
    cv->SetLogy(true);
    fPads1->SetLogy(true);
  }
  //-- resize y axis --
  hs->SetMaximum(100 * maximumForStack);
  //
  cv->SetLogy(true);
  cv->Update();
  cv->Modified(true);
  cv->Draw();
  sprintf(buffer, "%s/png/LOG_can_%s.png", nameOutDir_.c_str(), histoName.c_str());
  cv->SaveAs(buffer);
 // sprintf(buffer, "%s/root/LOG_can_%s.root", nameOutDir_.c_str(), histoName.c_str());
 // cv->SaveAs(buffer);
  sprintf(buffer, "%s/pdf/LOG_can_%s.pdf", nameOutDir_.c_str(), histoName.c_str());
  cv->SaveAs(buffer);

 }



