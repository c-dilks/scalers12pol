//void ListRuns(TString fn="counts_bd10.root") {
void ListRuns(TString fn="counts_bd11.root") {
//void ListRuns(TString fn="rtree.root") {
  TFile *f=new TFile(fn.Data(),"READ");
  TTree *t=(TTree*) f->Get("sca");
  //TTree *t=(TTree*) f->Get("rellum");
  Int_t r;
  t->SetBranchAddress("runnum",&r);
  for(int i=0;i<t->GetEntries(); i++) {
    t->GetEntry(i);
    printf("%d\n",r);
  };
};
