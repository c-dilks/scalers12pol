// makes some plots to explore properties of rate-safe counting

void RateSafeTests(const char * filename="rootfiles/all.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * m = (TTree*) infile->Get("matx");




  // open matx tree
  Int_t i,fi,runnum,fill,tbit,cbit,bx;
  Double_t t,raw,acc,mul,fac,rsc,fac2;
  Bool_t kicked;
  m->SetBranchAddress("i",&i);
  m->SetBranchAddress("fi",&fi);
  m->SetBranchAddress("runnum",&runnum);
  m->SetBranchAddress("fill",&fill);
  m->SetBranchAddress("tbit",&tbit);
  m->SetBranchAddress("cbit",&cbit);
  m->SetBranchAddress("bx",&bx);
  m->SetBranchAddress("t",&t);
  m->SetBranchAddress("raw",&raw);
  m->SetBranchAddress("acc",&acc);
  m->SetBranchAddress("mul",&mul);
  m->SetBranchAddress("fac",&fac);
  m->SetBranchAddress("rsc",&rsc);
  m->SetBranchAddress("fac2",&fac2);
  m->SetBranchAddress("kicked",&kicked);


  // get maxima for plots
  Double_t max_vpd_mul_rate=0; // maximum of multiples corrected vpdx / t
  Double_t max_vpd_rsc_rate=0; // maximum of rate-safe vpd / t
  Double_t max_rat_mul=0; // maximum of multiples corrected zdcx/vpdx
  Double_t max_rat_rsc=0; // maximum of rate-safe zdc/vpd
  Double_t bbcx,vpdx,zdcx,bbc_rsc,vpd_rsc,zdc_rsc,time;
  Int_t i_tmp=0;
  for(Int_t e=0; e<m->GetEntries(); e++)
  {
    m->GetEntry(e);
    if(i_tmp==0) i_tmp=i;

    if(cbit==2)
    {
      if(i!=i_tmp || e+1==m->GetEntries())
      {
        i_tmp = i;
        max_vpd_mul_rate = (vpdx/time > max_vpd_mul_rate) ? vpdx/time : max_vpd_mul_rate;
        max_vpd_rsc_rate = (vpd_rsc/time > max_vpd_rsc_rate) ? vpd_rsc/time : max_vpd_rsc_rate;
        max_rat_mul = (zdcx/vpdx > max_rat_mul) ? zdcx/vpdx : max_rat_mul;
        max_rat_rsc = (zdc_rsc/vpd_rsc > max_rat_rsc) ? zdc_rsc/vpd_rsc : max_rat_rsc;
        bbcx=zdcx=vpdx=bbc_rsc=zdc_rsc=vpd_rsc=0; 
      };

      if(i==i_tmp && kicked==0)
      {
        time=t;
        if(tbit==0)
        {
          bbcx += mul;
          bbc_rsc += rsc;
        }
        else if(tbit==1)
        {
          zdcx += mul;
          zdc_rsc += rsc;
        }
        else if(tbit==2)
        {
          vpdx += mul;
          vpd_rsc += rsc;
        };
      };
    };
  };
  printf("max_vpd_mul_rate=%f\n",max_vpd_mul_rate);
  printf("max_vpd_rsc_rate=%f\n",max_vpd_rsc_rate);
  printf("max_rat_mul=%f\n",max_rat_mul);
  printf("max_rat_rsc=%f\n",max_rat_rsc);




  // define 2-d distributions
  const Int_t NBINS = 250;
  TH2D * rat_mul_dist = new TH2D("rat_mul_dist",
    "N_{ZDCX}^{mul} / N_{VPDX}^{mul} vs. N_{VPDX}^{mul} / t; N_{VPDX}^{mul} / t; N_{ZDCX}^{mul} / N_{VPDX}^{mul}",
    NBINS,0,max_vpd_mul_rate,NBINS,0,max_rat_mul);
  TH2D * rat_rsc_dist = new TH2D("rat_rsc_dist",
    "#epsilon_{ZDCE}#epsilon_{ZDCW}N_{ZDC}^{rsc} / #epsilon_{VPDE}#epsilon_{VPDW}N_{VPD}^{rsc} vs. #epsilon_{VPDE}#epsilon_{VPDW}N_{VPD}^{rsc} / t; #epsilon_{VPDE}#epsilon_{VPDW}N_{VPD}^{rsc} / t; #epsilon_{ZDCE}#epsilon_{ZDCW}N_{ZDC}^{rsc} / #epsilon_{VPDE}#epsilon_{VPDW}N_{VPD}^{rsc} ",
    NBINS,0,max_vpd_rsc_rate,NBINS,0,max_rat_rsc);


  // fill 2-d distributions
  Double_t vpd_mul_rate,vpd_rsc_rate,rat_mul,rat_rsc;
  i_tmp = 0;
  for(Int_t e=0; e<m->GetEntries(); e++)
  {
    m->GetEntry(e);
    if(i_tmp==0) i_tmp=i;

    if(cbit==2)
    {
      if(i!=i_tmp || e+1==m->GetEntries())
      {
        i_tmp=i;
        vpd_mul_rate = vpdx/time;
        vpd_rsc_rate = vpd_rsc/time;
        rat_mul = zdcx/vpdx;
        rat_rsc = zdc_rsc/vpd_rsc;
        rat_mul_dist->Fill(vpd_mul_rate,rat_mul);
        rat_rsc_dist->Fill(vpd_rsc_rate,rat_rsc);
        bbcx=zdcx=vpdx=bbc_rsc=zdc_rsc=vpd_rsc=0; 
      };

      if(i==i_tmp && kicked==0)
      {
        time=t;
        if(tbit==0)
        {
          bbcx += mul;
          bbc_rsc += rsc;
        }
        else if(tbit==1)
        {
          zdcx += mul;
          zdc_rsc += rsc;
        }
        else if(tbit==2)
        {
          vpdx += mul;
          vpd_rsc += rsc;
        };
      };
    };
  };

  TCanvas * c_rat_mul_dist = new TCanvas("c_rat_mul_dist","c_rat_mul_dist",1200,450);
  c_rat_mul_dist->SetGrid(1,1);
  rat_mul_dist->Draw("colz");
  TCanvas * c_rat_rsc_dist = new TCanvas("c_rat_rsc_dist","c_rat_rsc_dist",1200,450);
  c_rat_rsc_dist->SetGrid(1,1);
  rat_rsc_dist->Draw("colz");
};
