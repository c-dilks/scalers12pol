// draws fill index vs. run index, to get a sense
// of how many runs there are per fill

void draw_fill_vs_run(const char * filename="counts.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sca");
  Double_t i_min = tr->GetMinimum("i");
  Double_t i_max = tr->GetMaximum("i");
  Int_t i_bins = (Int_t)(i_max-i_min+1);
  Double_t fi_min = tr->GetMinimum("fi");
  Double_t fi_max = tr->GetMaximum("fi");
  Int_t fi_bins = (Int_t)(fi_max-fi_min+1);
  TH2F * ff = new TH2F("ff","fi :: i",i_bins,i_min,i_max,fi_bins,fi_min,fi_max);
  ff->GetXaxis()->SetTitle("i");
  ff->GetYaxis()->SetTitle("fi");
  TCanvas * cc = new TCanvas("cc","cc",700,500);
  cc->SetGrid(1,1);
  gStyle->SetOptStat(0);
  tr->Draw("fi:i>>ff","","colz");
};
