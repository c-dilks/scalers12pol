// checks correlation between run time and nbx shift problem

void nbx_time_correlation(const char * filename="counts.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sca");

  Int_t i,bx;
  Double_t t,tot_bx;
  tr->SetBranchAddress("i",&i);
  tr->SetBranchAddress("bx",&bx);
  tr->SetBranchAddress("t",&t);
  tr->SetBranchAddress("tot_bx",&tot_bx);

  Int_t NRUNS_tmp = tr->GetMaximum("i");
  const Int_t NRUNS = NRUNS_tmp;
  Double_t ratio[NRUNS];
  Double_t numer[NRUNS];
  Double_t denom[NRUNS];
  Double_t time[NRUNS];

  for(Int_t e=0; e<tr->GetEntries(); e++)
  {
    tr->GetEntry(e);
    if(bx==0) denom[i]=tot_bx;
    else if(bx==119) numer[i]=tot_bx;
    time[i]=t;
  };


  TH2D * td = new TH2D("td","N_{bx}^{(119)} / N_{bx}^{(0)}  vs. run time",100,0,4000,50,0.9996,1);
  for(Int_t e=0; e<NRUNS; e++) 
  {
    ratio[e] = numer[e] / denom[e];
    td->Fill(time[e],ratio[e]);
  };

  TCanvas * cc = new TCanvas("cc","cc",1200,500);
  td->Draw("colz");

  TProfile * td_pfx = (TProfile*) td->ProfileX();
  td_pfx->SetLineColor(kRed);
  td_pfx->SetLineWidth(2);
  td_pfx->Draw("same");
};
  
