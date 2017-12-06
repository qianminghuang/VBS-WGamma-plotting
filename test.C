TH1F* test(TTree * treeD,TTree * treeMC){
 

 float inix= 0;
 float finx= 60;
 float nbinx= 60.0; 
 
int nVtx;
treeD->SetBranchAddress("nVtx", &nVtx);

Long64_t nentries = treeD->GetEntriesFast();

TH1F *hD= new TH1F("hD","h500",nbinx,inix,finx);
TH1F *hRatio= new TH1F("hRatio","hRatio500",nbinx,inix,finx);
 hD->Sumw2();
 hRatio->Sumw2();
    for(int i2=0; i2<= nentries; i2++)
  //      for(int i2=0; i2<= 100; i2++)
{
    treeD->GetEntry(i2);
       hD->Fill(nVtx);
       hRatio->Fill(nVtx);
 }

cout<< "nDCount: " << endl;
double nDError;
cout<< "nDCount1: " << endl;
double nDCount = hD->IntegralAndError(1,nbinx-1,nDError);
cout<< "nDCount: " << nDCount<<"nDError: "<<nDError<<endl;

nDCount = hD->IntegralAndError(1,nbinx,nDError);
cout<< "nDCount: " << nDCount<<"nDError: "<<nDError<<endl;

treeMC->SetBranchAddress("nVtx", &nVtx);

Long64_t nentriesMC = treeMC->GetEntriesFast();
TH1F *h1= new TH1F("h1","h500",nbinx,inix,finx);
h1->Sumw2();
    for(int i=0; i<= nentriesMC; i++)
   //     for(int i=0; i<= 100; i++)
{
    treeMC->GetEntry(i);
    h1->Fill(nVtx);
}

double n1Error;
double n1Count = h1->IntegralAndError(1,nbinx-1,n1Error);
cout<<"n1Count: "<<n1Count<<"n1Error: "<<n1Error<<endl;


hRatio->Scale(n1Count/nDCount);
hRatio->Divide(h1);  
 
    return hRatio;



}
