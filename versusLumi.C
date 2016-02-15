// plots correction factor vs. instantaneous luminosity

void versusLumi(const char * filename = "lum_perrun_FMSJP2.txt")
{
  // read rellum tree and build lumi tree
  TFile * rellum_file = new TFile("rdat_i.root","READ");
  TTree * rtr = (TTree*) rellum_file->Get("rtr");
  TTree * ltr = new TTree();
  ltr->ReadFile(filename,"runnum/I:time_start/D:time_stop/D:fill/I:lumi/F:prescale/F:livetime_trig/F:base/C:livetime_base/F:fom/F:fomp4/F");
  Int_t runnum_l, runnum_r;
  Float_t lumi;
  Double_t time_start, time_stop;
  rtr->SetBranchAddress("runnum",&runnum_r);
  ltr->SetBranchAddress("runnum",&runnum_l);
  ltr->SetBranchAddress("lumi",&lumi);
  ltr->SetBranchAddress("time_start",&time_start);
  ltr->SetBranchAddress("time_stop",&time_stop);


  // strings
  char tbit[3][4];
  sprintf(tbit[0],"bbc");
  sprintf(tbit[1],"zdc");
  sprintf(tbit[2],"vpd");
  char cbit[3][4];
  sprintf(cbit[0],"e");
  sprintf(cbit[1],"w");
  sprintf(cbit[2],"x");
  char sbit[5][4];
  sprintf(sbit[0],"nn");
  sprintf(sbit[1],"np");
  sprintf(sbit[2],"pn");
  sprintf(sbit[3],"pp");
  sprintf(sbit[4],"all");


  // set run index minimum & maximum
  Int_t var_l = rtr->GetMinimum("i");
  Int_t var_h = rtr->GetMaximum("i");
  var_h++;
  Int_t var_bins = var_h-var_l;



  // lumi bins
  Int_t lumi_bins = 100;

  Float_t lumit_l = 0;
  Float_t lumit_h = 0;
  for(Int_t l=0; l<ltr->GetEntries(); l++)
  {
    ltr->GetEntry(l);
    lumit_h = (lumi/(time_stop-time_start) > lumit_h) ? lumi/(time_stop-time_start):lumit_h;
  };

  Float_t lumi_l = 0;
  Float_t lumi_h = ltr->GetMaximum("lumi");;

  //lumi_h = 0.06; // override
  //lumit_h = 0.03e-3; // override

  // obtain 1d hists
  TH1D * fac_hist[3][3][4]; // [tbit] [cbit] [spinbit]
  char fac_hist_name[3][3][4][64]; // [tbit] [cbit] [spinbit]
  TCanvas * fac_canvas[3]; // [tbit]
  char fac_canvas_name[3][64]; // [tbit];
  Double_t fac_l[3][3][5];
  Double_t fac_h[3][3][5];
  TH1D * raw_hist[3][3][4]; // [tbit] [cbit] [spinbit]
  char raw_hist_name[3][3][4][64]; // [tbit] [cbit] [spinbit]
  TCanvas * raw_canvas[3]; // [tbit]
  char raw_canvas_name[3][64]; // [tbit];
  Double_t raw_l[3][3][5];
  Double_t raw_h[3][3][5];
  TH1D * mul_hist[3][3][4]; // [tbit] [cbit] [spinbit]
  char mul_hist_name[3][3][4][64]; // [tbit] [cbit] [spinbit]
  TCanvas * mul_canvas[3]; // [tbit]
  char mul_canvas_name[3][64]; // [tbit];
  Double_t mul_l[3][3][5];
  Double_t mul_h[3][3][5];


  for(Int_t t=0; t<3; t++)
  {
    sprintf(fac_canvas_name[t],"c_fac_%s",tbit[t]);
    sprintf(raw_canvas_name[t],"c_raw_%s",tbit[t]);
    sprintf(mul_canvas_name[t],"c_mul_%s",tbit[t]);
    fac_canvas[t] = (TCanvas*) rellum_file->Get(fac_canvas_name[t]);
    raw_canvas[t] = (TCanvas*) rellum_file->Get(raw_canvas_name[t]);
    mul_canvas[t] = (TCanvas*) rellum_file->Get(mul_canvas_name[t]);
    for(Int_t c=0; c<3; c++)
    {
      for(Int_t s=0; s<4; s++)
      {
        sprintf(fac_hist_name[t][c][s],"fac_%s%s_%s",tbit[t],cbit[c],sbit[s]);
        sprintf(raw_hist_name[t][c][s],"raw_%s%s_%s",tbit[t],cbit[c],sbit[s]);
        sprintf(mul_hist_name[t][c][s],"mul_%s%s_%s",tbit[t],cbit[c],sbit[s]);
        fac_hist[t][c][s] = (TH1D*) fac_canvas[t]->GetPad(c+1)->GetPrimitive(fac_hist_name[t][c][s]);
        raw_hist[t][c][s] = (TH1D*) raw_canvas[t]->GetPad(c+1)->GetPrimitive(raw_hist_name[t][c][s]);
        mul_hist[t][c][s] = (TH1D*) mul_canvas[t]->GetPad(c+1)->GetPrimitive(mul_hist_name[t][c][s]);
        fac_l[t][c][s] = fac_hist[t][c][s]->GetMinimum();
        fac_h[t][c][s] = fac_hist[t][c][s]->GetMaximum();
        raw_l[t][c][s] = raw_hist[t][c][s]->GetMinimum();
        raw_h[t][c][s] = raw_hist[t][c][s]->GetMaximum();
        mul_l[t][c][s] = mul_hist[t][c][s]->GetMinimum();
        mul_h[t][c][s] = mul_hist[t][c][s]->GetMaximum();
      };
    };
  };

  // set spinbit=4 minima & maxima to be minimum minimum and maximum maximum
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      fac_l[t][c][4] = 1e12;
      fac_h[t][c][4] = 0;
      raw_l[t][c][4] = 1e12;
      raw_h[t][c][4] = 0;
      mul_l[t][c][4] = 1e12;
      mul_h[t][c][4] = 0;
      for(Int_t s=0; s<4; s++)
      {
        fac_l[t][c][4] = (fac_l[t][c][s]<fac_l[t][c][4]) ? fac_l[t][c][s]:fac_l[t][c][4];
        fac_h[t][c][4] = (fac_h[t][c][s]>fac_h[t][c][4]) ? fac_h[t][c][s]:fac_h[t][c][4];
        raw_l[t][c][4] = (raw_l[t][c][s]<raw_l[t][c][4]) ? raw_l[t][c][s]:raw_l[t][c][4];
        raw_h[t][c][4] = (raw_h[t][c][s]>raw_h[t][c][4]) ? raw_h[t][c][s]:raw_h[t][c][4];
        mul_l[t][c][4] = (mul_l[t][c][s]<mul_l[t][c][4]) ? mul_l[t][c][s]:mul_l[t][c][4];
        mul_h[t][c][4] = (mul_h[t][c][s]>mul_h[t][c][4]) ? mul_h[t][c][s]:mul_h[t][c][4];
      };
    };
  };

  
  // vs. lumi plots
  TH2D * fac_vs_lumi[3][3][5];
  char fac_vs_lumi_title[3][3][5][128];
  char fac_vs_lumi_name[3][3][5][128];
  Double_t fac;
  TH2D * raw_vs_lumi[3][3][5];
  char raw_vs_lumi_title[3][3][5][128];
  char raw_vs_lumi_name[3][3][5][128];
  Double_t raw;
  TH2D * mul_vs_lumi[3][3][5];
  char mul_vs_lumi_title[3][3][5][128];
  char mul_vs_lumi_name[3][3][5][128];
  Double_t mul;

  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      printf("t=%d c=%d\n",t,c);
      for(Int_t s=0; s<5; s++)
      {
        sprintf(fac_vs_lumi_name[t][c][s],"fac_%s%s_%s_vs_lumit",tbit[t],cbit[c],sbit[s]);
        sprintf(fac_vs_lumi_title[t][c][s],"fac_%s%s_%s vs. lumi/t",tbit[t],cbit[c],sbit[s]);
        fac_vs_lumi[t][c][s] = new TH2D(fac_vs_lumi_name[t][c][s],fac_vs_lumi_title[t][c][s],
          lumi_bins,lumit_l,lumit_h,lumi_bins,fac_l[t][c][s],fac_h[t][c][s]);
        sprintf(raw_vs_lumi_name[t][c][s],"raw_%s%s_%s_vs_lumi",tbit[t],cbit[c],sbit[s]);
        sprintf(raw_vs_lumi_title[t][c][s],"raw_%s%s_%s vs. lumi",tbit[t],cbit[c],sbit[s]);
        raw_vs_lumi[t][c][s] = new TH2D(raw_vs_lumi_name[t][c][s],raw_vs_lumi_title[t][c][s],
          lumi_bins,lumi_l,lumi_h,lumi_bins,raw_l[t][c][s],raw_h[t][c][s]);
        sprintf(mul_vs_lumi_name[t][c][s],"mul_%s%s_%s_vs_lumi",tbit[t],cbit[c],sbit[s]);
        sprintf(mul_vs_lumi_title[t][c][s],"mul_%s%s_%s vs. lumi",tbit[t],cbit[c],sbit[s]);
        mul_vs_lumi[t][c][s] = new TH2D(mul_vs_lumi_name[t][c][s],mul_vs_lumi_title[t][c][s],
          lumi_bins,lumi_l,lumi_h,lumi_bins,mul_l[t][c][s],mul_h[t][c][s]);
      }
      for(Int_t s=0; s<4; s++)
      {
        for(Int_t b=var_l; b<var_h; b++)
        {
          fac = fac_hist[t][c][s]->GetBinContent(b);
          raw = raw_hist[t][c][s]->GetBinContent(b);
          mul = mul_hist[t][c][s]->GetBinContent(b);
          rtr->GetEntry(b-1); // tree index = run index - 1, where bin number = run index
          for(Int_t l=0; l<ltr->GetEntries(); l++)
          {
            ltr->GetEntry(l);
            if(runnum_l == runnum_r) 
            {
              fac_vs_lumi[t][c][s]->Fill(lumi/(time_stop-time_start),fac);
              fac_vs_lumi[t][c][4]->Fill(lumi/(time_stop-time_start),fac);
              raw_vs_lumi[t][c][s]->Fill(lumi,raw);
              raw_vs_lumi[t][c][4]->Fill(lumi,raw);
              mul_vs_lumi[t][c][s]->Fill(lumi,mul);
              mul_vs_lumi[t][c][4]->Fill(lumi,mul);
            };
          };
        };
      };
    };
  };




  TFile * outfile = new TFile("lumi.root","RECREATE");
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      for(Int_t s=0; s<5; s++) fac_vs_lumi[t][c][s]->Write();
    };
  };
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      raw_vs_lumi[t][c][4]->Write();
    };
  };
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      mul_vs_lumi[t][c][4]->Write();
    };
  };
};







