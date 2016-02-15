// THIS ONE IS FOR THE SPIN PATTERN TESTS DONE ON AUG. 5th 2015
//
// draws distribution of asymmetry determined from bunch fit
//
// requires the following files, all of which are produced from RunPatterns:
//  - fit_result.[numer].[denom].root
//  - pats/fit_result.[numer].[denom].pat*.root

void DrawProjections2(const char * testchar="",
                      const char * numer="zdcx",
                      const char * denom="vpdx")
{
  const Int_t NBINS = 90;
  gStyle->SetOptFit(1);
  char filename[256];
  if(!strcmp(testchar,""))
    sprintf(filename,"fit_result.%s.%s.root",numer,denom);
  else 
    sprintf(filename,"tests/%s/fit_result.%s.%s.root",
      testchar,numer,denom);
  TFile * infile = new TFile(filename,"READ");
  TH1D * hist_all = (TH1D*) infile->Get("/asymmetry/asym_a3_v_run");

  TFile * patfile[8];
  char patfile_n[8][256];
  Int_t pat_arr[8] = {13,14,23,24,31,32,41,42};
  //Color_t colours[8] = {kOrange,kRed,kMagenta,kBlue,kCyan+1,kGreen+1,kYellow+2,kViolet-6}; 
  Color_t colours[8] = {kBlue,kRed,kGreen+1,kMagenta,kBlue,kGreen+1,kRed,kMagenta}; // pairs
  TH1D * hist_pat[8];
  for(Int_t pp=0; pp<8; pp++)
  {
    if(!strcmp(testchar,""))
      sprintf(patfile_n[pp],"pats/fit_result.%s.%s.pat%d.root",numer,denom,pat_arr[pp]);
    else
      sprintf(patfile_n[pp],"tests/%s/pats/fit_result.%s.%s.pat%d.root",
        testchar,numer,denom,pat_arr[pp]);
    patfile[pp] = new TFile(patfile_n[pp],"READ");
    hist_pat[pp] = (TH1D*) patfile[pp]->Get("/asymmetry/asym_a3_v_run");
  };

  // bin boundaries
  Double_t LBOUND=10;
  Double_t UBOUND=-10;
  for(Int_t pp=0; pp<8; pp++)
  {
    LBOUND = (hist_pat[pp]->GetMinimum() < LBOUND) ? hist_pat[pp]->GetMinimum():LBOUND;
    UBOUND = (hist_pat[pp]->GetMaximum() > UBOUND) ? hist_pat[pp]->GetMaximum():UBOUND;
  };
  //LBOUND -= (UBOUND-LBOUND)*0.05;
  //UBOUND += (UBOUND-LBOUND)*0.05;
  LBOUND = -0.6e-2; // override
  UBOUND = 0.6e-2;  // override
  printf("LBOUND=%f UBOUND=%f\n",LBOUND,UBOUND);


  char omit_str[64];
  if(!strcmp(testchar,"")) strcpy(omit_str,"");
  else sprintf(omit_str," :: test %s",testchar);

  char asym_dist_all_n[256];
  char asym_dist_all_t[256];
  sprintf(asym_dist_all_n,"asym_dist_all");
  //sprintf(asym_dist_all_t,"Run 12 S_{LL} distribution for %s/%s%s",numer,denom,omit_str);
  sprintf(asym_dist_all_t,"Run 12 S_{LL} distribution for rate-safe %s/%s%s",numer,denom,omit_str);
  TH1D * asym_dist_all = new TH1D(asym_dist_all_n,asym_dist_all_t,NBINS,LBOUND,UBOUND);
  asym_dist_all->SetLineWidth(2);
  asym_dist_all->SetLineColor(kBlack);

  char asym_dist_pat_n[8][256];
  char asym_dist_pat_t[8][256];
  TH1D * asym_dist_pat[8];
  for(Int_t pp=0; pp<8; pp++)
  {
    sprintf(asym_dist_pat_n[pp],"asym_dist_pat%d",pp);
    //sprintf(asym_dist_pat_t[pp],"Run 12 S_{LL} distribution for %s/%s%s",numer,denom,omit_str);
    sprintf(asym_dist_pat_t[pp],"Run 12 S_{LL} distribution for rate-safe %s/%s%s",numer,denom,omit_str);
    asym_dist_pat[pp] = new TH1D(asym_dist_pat_n[pp],asym_dist_pat_t[pp],NBINS,LBOUND,UBOUND);
    asym_dist_pat[pp]->SetLineColor(colours[pp]);
  };

  Double_t bc;

  for(Int_t b=1; b<=hist_all->GetNbinsX(); b++)
  {
    bc = hist_all->GetBinContent(b);
    asym_dist_all->Fill(bc);
  };

  
  for(Int_t pp=0; pp<8; pp++)
  {
    for(Int_t b=1; b<=hist_pat[pp]->GetNbinsX(); b++)
    {
      bc = hist_pat[pp]->GetBinContent(b);
      if(bc!=0) asym_dist_pat[pp]->Fill(bc);
    };
  };

  // run 13 fits (2 gaussians)
  TF1 * gaus1 = new TF1("gaus1","gaus",LBOUND,0);
  TF1 * gaus2 = new TF1("gaus2","gaus",0,UBOUND);
  //asym_dist_all->Fit(gaus1,"","",LBOUND,0);
  //asym_dist_all->Fit(gaus2,"","",0,UBOUND);
  //c1->Close();


  TCanvas * cc = new TCanvas("cc","cc",1600,400);
  cc->SetGrid(1,0);
  gStyle->SetOptStat(1100);
  gStyle->SetStatFontSize(0.1);
  Float_t size = 0.05;
  asym_dist_all->GetXaxis()->SetLabelSize(size);
  asym_dist_all->GetYaxis()->SetLabelSize(size);
  asym_dist_all->Draw();
  gaus1->Draw("same");
  gaus2->Draw("same");
  for(Int_t pp=0; pp<8; pp++)
  {
    asym_dist_pat[pp]->GetXaxis()->SetLabelSize(size);
    asym_dist_pat[pp]->GetYaxis()->SetLabelSize(size);
    asym_dist_pat[pp]->Draw("same");
  };
  asym_dist_all->Draw("same");

  Double_t sigma1 = gaus1->GetParameter(2);
  Double_t sigma2 = gaus2->GetParameter(2);
  Double_t sigma = (sigma1>sigma2)?sigma1:sigma2;
  Double_t mean = asym_dist_all->GetMean();
  Double_t sys = sigma + fabs(mean);
  printf("sigma1 = %f\nsigma2 = %f\n",sigma1,sigma2);
  printf("sigma = %f\nmean = %f\nsystematic = %f\n",sigma,mean,sys);


  char printname[64];
  if(!strcmp(testchar,"")) strcpy(printname,"projection.png");
  else sprintf(printname,"projection.test%s.png",testchar);
  cc->Print(printname,"png");
    
};
