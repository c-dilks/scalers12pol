// plots total bXings / bXing rate vs. run time
// -- should be equal

void nbx_check(const char * filename="sums.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sum");
  Int_t i_max = tr->GetMaximum("i");
  Int_t fi_max = tr->GetMaximum("fi");

  TH2D * tt = new TH2D("tt","total bXings / bXing rate vs. run time == #tau vs. t",
    100,0,4000,100,0,4000);
  TH1D * diff_vs_i = new TH1D("diff_vs_i","t-#tau vs. run index",i_max,1,i_max+1);
  TH1D * diff_vs_fi = new TH1D("diff_vs_fi","average t-#tau in a fill vs. fill index",fi_max,1,fi_max+1);
  TH1D * quot_vs_i = new TH1D("quot_vs_i","t/#tau vs. run index",i_max,1,i_max+1);
  TH1D * quot_vs_fi = new TH1D("quot_vs_fi","average t/#tau in a fill vs. fill index",fi_max,1,fi_max+1);

  tr->Project("tt","tau:t");
  tr->Project("diff_vs_i","i","t-tau");
  tr->Project("diff_vs_fi","fi","(t-tau)/num_runs");
  tr->Project("quot_vs_i","i","t/tau");
  tr->Project("quot_vs_fi","fi","(t/tau)/num_runs");

  tt->Fit("pol1","","",0,4000);
  c1->Close();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TCanvas * ctt = new TCanvas("ctt","ctt",700,500);
  ctt->SetGrid(1,1);
  tt->Draw("colz");
  ctt->Print("nbx_check/ctt.png","png");

  TCanvas * diff_canv = new TCanvas("diff_canv","diff_canv",1300,700);
  diff_canv->Divide(1,2);
  for(Int_t x=1; x<=2; x++) diff_canv->GetPad(x)->SetGrid(1,1);
  diff_canv->cd(1); diff_vs_i->Draw();
  diff_canv->cd(2); diff_vs_fi->Draw();
  diff_canv->Print("nbx_check/diff_canv.png","png");

  TCanvas * quot_canv = new TCanvas("quot_canv","quot_canv",1300,700);
  quot_canv->Divide(1,2);
  for(Int_t x=1; x<=2; x++) quot_canv->GetPad(x)->SetGrid(1,1);
  quot_canv->cd(1); quot_vs_i->Draw();
  quot_canv->cd(2); quot_vs_fi->Draw();
  quot_canv->Print("nbx_check/quot_canv.png","png");
};

