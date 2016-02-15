// computes correlations between triggers 

void correlations(const char * filename="counts.root")
{
  // read data tree
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sca");
  Int_t index,t,runnum;
  Int_t bbce,bbcw,bbcx;
  Int_t zdce,zdcw,zdcx;
  Int_t vpde,vpdw,vpdx;
  tr->SetBranchAddress("t",&t);
  tr->SetBranchAddress("i",&index);
  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("bbce",&bbce);
  tr->SetBranchAddress("bbcw",&bbcw);
  tr->SetBranchAddress("bbcx",&bbcx);
  tr->SetBranchAddress("zdce",&zdce);
  tr->SetBranchAddress("zdcw",&zdcw);
  tr->SetBranchAddress("zdcx",&zdcx);
  tr->SetBranchAddress("vpde",&vpde);
  tr->SetBranchAddress("vpdw",&vpdw);
  tr->SetBranchAddress("vpdx",&vpdx);


  // determine maximum trigger rates
  /*
  Float_t bbce_max_rate, bbcw_max_rate, bbcx_max_rate;
  Float_t zdce_max_rate, zdcw_max_rate, zdcx_max_rate;
  Float_t vpde_max_rate, vpdw_max_rate, vpdx_max_rate;
  Float_t bbce_rate, bbcw_rate, bbcx_rate;
  Float_t zdce_rate, zdcw_rate, zdcx_rate;
  Float_t vpde_rate, vpdw_rate, vpdx_rate;
  bbce_max_rate = bbcw_max_rate = bbcx_max_rate = 0;
  zdce_max_rate = zdcw_max_rate = zdcx_max_rate = 0;
  vpde_max_rate = vpdw_max_rate = vpdx_max_rate = 0;
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    bbce_rate = ((Float_t)bbce)/((Float_t)t);
    bbcw_rate = ((Float_t)bbcw)/((Float_t)t);
    bbcx_rate = ((Float_t)bbcx)/((Float_t)t);
    zdce_rate = ((Float_t)zdce)/((Float_t)t);
    zdcw_rate = ((Float_t)zdcw)/((Float_t)t);
    zdcx_rate = ((Float_t)zdcx)/((Float_t)t);
    vpde_rate = ((Float_t)vpde)/((Float_t)t);
    vpdw_rate = ((Float_t)vpdw)/((Float_t)t);
    vpdx_rate = ((Float_t)vpdx)/((Float_t)t);

    bbce_max_rate = (bbce_rate > bbce_max_rate) ? bbce_rate : bbce_max_rate;
    bbcw_max_rate = (bbcw_rate > bbcw_max_rate) ? bbcw_rate : bbcw_max_rate;
    bbcx_max_rate = (bbcx_rate > bbcx_max_rate) ? bbcx_rate : bbcx_max_rate;
    zdce_max_rate = (zdce_rate > zdce_max_rate) ? zdce_rate : zdce_max_rate;
    zdcw_max_rate = (zdcw_rate > zdcw_max_rate) ? zdcw_rate : zdcw_max_rate;
    zdcx_max_rate = (zdcx_rate > zdcx_max_rate) ? zdcx_rate : zdcx_max_rate;
    vpde_max_rate = (vpde_rate > vpde_max_rate) ? vpde_rate : vpde_max_rate;
    vpdw_max_rate = (vpdw_rate > vpdw_max_rate) ? vpdw_rate : vpdw_max_rate;
    vpdx_max_rate = (vpdx_rate > vpdx_max_rate) ? vpdx_rate : vpdx_max_rate;
  };
  Float_t factor=0.1;
  bbce_max_rate += factor*bbce_max_rate;  
  bbcw_max_rate += factor*bbcw_max_rate;
  bbcx_max_rate += factor*bbcx_max_rate;
  zdce_max_rate += factor*zdce_max_rate;
  zdcw_max_rate += factor*zdcw_max_rate;
  zdcx_max_rate += factor*zdcx_max_rate;
  vpde_max_rate += factor*vpde_max_rate;
  vpdw_max_rate += factor*vpdw_max_rate;
  vpdx_max_rate += factor*vpdx_max_rate;


  // initialise correlation histograms
  const Int_t NBINS = 100;
  TH2F * bbc_ew = new TH2F("bbc_ew","BBCE vs BBCW Rate Correlation",NBINS,0,bbcw_max_rate,NBINS,0,bbce_max_rate);
  TH2F * bbc_ex = new TH2F("bbc_ex","BBCE vs BBCX Rate Correlation",NBINS,0,bbcx_max_rate,NBINS,0,bbce_max_rate);
  TH2F * bbc_wx = new TH2F("bbc_wx","BBCW vs BBCX Rate Correlation",NBINS,0,bbcx_max_rate,NBINS,0,bbcw_max_rate);
  TH2F * zdc_ew = new TH2F("zdc_ew","ZDCE vs ZDCW Rate Correlation",NBINS,0,zdcw_max_rate,NBINS,0,zdce_max_rate);
  TH2F * zdc_ex = new TH2F("zdc_ex","ZDCE vs ZDCX Rate Correlation",NBINS,0,zdcx_max_rate,NBINS,0,zdce_max_rate);
  TH2F * zdc_wx = new TH2F("zdc_wx","ZDCW vs ZDCX Rate Correlation",NBINS,0,zdcx_max_rate,NBINS,0,zdcw_max_rate);
  TH2F * vpd_ew = new TH2F("vpd_ew","VPDE vs VPDW Rate Correlation",NBINS,0,vpdw_max_rate,NBINS,0,vpde_max_rate);
  TH2F * vpd_ex = new TH2F("vpd_ex","VPDE vs VPDX Rate Correlation",NBINS,0,vpdx_max_rate,NBINS,0,vpde_max_rate);
  TH2F * vpd_wx = new TH2F("vpd_wx","VPDW vs VPDX Rate Correlation",NBINS,0,vpdx_max_rate,NBINS,0,vpdw_max_rate);

  TH2F * bz_e = new TH2F("bz_e","BBCE vs ZDCE Rate Correlation",NBINS,0,zdce_max_rate,NBINS,0,bbce_max_rate);
  TH2F * bz_w = new TH2F("bz_w","BBCW vs ZDCW Rate Correlation",NBINS,0,zdcw_max_rate,NBINS,0,bbcw_max_rate);
  TH2F * bz_x = new TH2F("bz_x","BBCX vs ZDCX Rate Correlation",NBINS,0,zdcx_max_rate,NBINS,0,bbcx_max_rate);
  TH2F * bv_e = new TH2F("bv_e","BBCE vs VPDE Rate Correlation",NBINS,0,vpde_max_rate,NBINS,0,bbce_max_rate);
  TH2F * bv_w = new TH2F("bv_w","BBCW vs VPDW Rate Correlation",NBINS,0,vpdw_max_rate,NBINS,0,bbcw_max_rate);
  TH2F * bv_x = new TH2F("bv_x","BBCX vs VPDX Rate Correlation",NBINS,0,vpdx_max_rate,NBINS,0,bbcx_max_rate);
  TH2F * vz_e = new TH2F("vz_e","VPDE vs ZDCE Rate Correlation",NBINS,0,zdce_max_rate,NBINS,0,vpde_max_rate);
  TH2F * vz_w = new TH2F("vz_w","VPDW vs ZDCW Rate Correlation",NBINS,0,zdcw_max_rate,NBINS,0,vpdw_max_rate);
  TH2F * vz_x = new TH2F("vz_x","VPDX vs ZDCX Rate Correlation",NBINS,0,zdcx_max_rate,NBINS,0,vpdx_max_rate);

  // projections
  tr->Project("bbc_ew","bbce/t:bbcw/t","bbce>0 && bbcw>0 && t>0");
  tr->Project("bbc_ex","bbce/t:bbcx/t","bbce>0 && bbcx>0 && t>0");
  tr->Project("bbc_wx","bbcw/t:bbcx/t","bbcw>0 && bbcx>0 && t>0");
  tr->Project("zdc_ew","zdce/t:zdcw/t","zdce>0 && zdcw>0 && t>0");
  tr->Project("zdc_ex","zdce/t:zdcx/t","zdce>0 && zdcx>0 && t>0");
  tr->Project("zdc_wx","zdcw/t:zdcx/t","zdcw>0 && zdcx>0 && t>0");
  tr->Project("vpd_ew","vpde/t:vpdw/t","vpde>0 && vpdw>0 && t>0");
  tr->Project("vpd_ex","vpde/t:vpdx/t","vpde>0 && vpdx>0 && t>0");
  tr->Project("vpd_wx","vpdw/t:vpdx/t","vpdw>0 && vpdx>0 && t>0");
  tr->Project("bz_e","bbce/t:zdce/t","bbce>0 && zdce>0 && t>0");
  tr->Project("bz_w","bbcw/t:zdcw/t","bbcw>0 && zdcw>0 && t>0");
  tr->Project("bz_x","bbcx/t:zdcx/t","bbcx>0 && zdcx>0 && t>0");
  tr->Project("bv_e","bbce/t:vpde/t","bbce>0 && vpde>0 && t>0");
  tr->Project("bv_w","bbcw/t:vpdw/t","bbcw>0 && vpdw>0 && t>0");
  tr->Project("bv_x","bbcx/t:vpdx/t","bbcx>0 && vpdx>0 && t>0");
  tr->Project("vz_e","vpde/t:zdce/t","vpde>0 && zdce>0 && t>0");
  tr->Project("vz_w","vpdw/t:zdcw/t","vpdw>0 && zdcw>0 && t>0");
  tr->Project("vz_x","vpdx/t:zdcx/t","vpdx>0 && zdcx>0 && t>0");
  */


  // define run total variables (total scaler counts in a single run)
  Int_t IMAX_tmp = tr->GetMaximum("i");
  const Int_t IMAX = IMAX_tmp;
  Long_t bbce_total[IMAX];
  Long_t bbcw_total[IMAX];
  Long_t bbcx_total[IMAX];
  Long_t zdce_total[IMAX];
  Long_t zdcw_total[IMAX];
  Long_t zdcx_total[IMAX];
  Long_t vpde_total[IMAX];
  Long_t vpdw_total[IMAX];
  Long_t vpdx_total[IMAX];
  Int_t time[IMAX];
  for(Int_t i=0; i<IMAX; i++)
  {
    bbce_total[i] = 0;
    bbcw_total[i] = 0;
    bbcx_total[i] = 0;
    zdce_total[i] = 0;
    zdcw_total[i] = 0;
    zdcx_total[i] = 0;
    vpde_total[i] = 0;
    vpdw_total[i] = 0;
    vpdx_total[i] = 0;
  };


  // tree loop --> fills run totals arrays
  Int_t ii;
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    ii = index-1; // run index starts at 1; arrays start at 0
    time[ii] = t;
    bbce_total[ii] += bbce;
    bbcw_total[ii] += bbcw;
    bbcx_total[ii] += bbcx;
    zdce_total[ii] += zdce;
    zdcw_total[ii] += zdcw;
    zdcx_total[ii] += zdcx;
    vpde_total[ii] += vpde;
    vpdw_total[ii] += vpdw;
    vpdx_total[ii] += vpdx;
  };


  // compute trigger rates & and max trigger rates
  Float_t bbce_rate[IMAX], bbcw_rate[IMAX], bbcx_rate[IMAX];
  Float_t zdce_rate[IMAX], zdcw_rate[IMAX], zdcx_rate[IMAX];
  Float_t vpde_rate[IMAX], vpdw_rate[IMAX], vpdx_rate[IMAX];
  Float_t bbce_rate_max, bbcw_rate_max, bbcx_rate_max;
  Float_t zdce_rate_max, zdcw_rate_max, zdcx_rate_max;
  Float_t vpde_rate_max, vpdw_rate_max, vpdx_rate_max;
  bbce_rate_max = bbcw_rate_max = bbcx_rate_max = 0;
  zdce_rate_max = zdcw_rate_max = zdcx_rate_max = 0;
  vpde_rate_max = vpdw_rate_max = vpdx_rate_max = 0;
  for(Int_t i=0; i<IMAX; i++)
  {
    bbce_rate[i] = ((Float_t) bbce_total[i])/((Float_t) time[i]);
    bbcw_rate[i] = ((Float_t) bbcw_total[i])/((Float_t) time[i]);
    bbcx_rate[i] = ((Float_t) bbcx_total[i])/((Float_t) time[i]);
    zdce_rate[i] = ((Float_t) zdce_total[i])/((Float_t) time[i]);
    zdcw_rate[i] = ((Float_t) zdcw_total[i])/((Float_t) time[i]);
    zdcx_rate[i] = ((Float_t) zdcx_total[i])/((Float_t) time[i]);
    vpde_rate[i] = ((Float_t) vpde_total[i])/((Float_t) time[i]);
    vpdw_rate[i] = ((Float_t) vpdw_total[i])/((Float_t) time[i]);
    vpdx_rate[i] = ((Float_t) vpdx_total[i])/((Float_t) time[i]);

    bbce_rate_max = (bbce_rate[i] > bbce_rate_max) ? bbce_rate[i] : bbce_rate_max;
    bbcw_rate_max = (bbcw_rate[i] > bbcw_rate_max) ? bbcw_rate[i] : bbcw_rate_max;
    bbcx_rate_max = (bbcx_rate[i] > bbcx_rate_max) ? bbcx_rate[i] : bbcx_rate_max;
    zdce_rate_max = (zdce_rate[i] > zdce_rate_max) ? zdce_rate[i] : zdce_rate_max;
    zdcw_rate_max = (zdcw_rate[i] > zdcw_rate_max) ? zdcw_rate[i] : zdcw_rate_max;
    zdcx_rate_max = (zdcx_rate[i] > zdcx_rate_max) ? zdcx_rate[i] : zdcx_rate_max;
    zdce_rate_max = (zdce_rate[i] > zdce_rate_max) ? zdce_rate[i] : zdce_rate_max;
    zdcw_rate_max = (zdcw_rate[i] > zdcw_rate_max) ? zdcw_rate[i] : zdcw_rate_max;
    zdcx_rate_max = (zdcx_rate[i] > zdcx_rate_max) ? zdcx_rate[i] : zdcx_rate_max;
  };
  Float_t factor=0.1;
  bbce_rate_max += factor*bbce_rate_max;
  bbcw_rate_max += factor*bbcw_rate_max;
  bbcx_rate_max += factor*bbcx_rate_max;
  zdce_rate_max += factor*zdce_rate_max;
  zdcw_rate_max += factor*zdcw_rate_max;
  zdcx_rate_max += factor*zdcx_rate_max;
  vpde_rate_max += factor*vpde_rate_max;
  vpdw_rate_max += factor*vpdw_rate_max;
  vpdx_rate_max += factor*vpdx_rate_max;
    

  // initialise correlation histograms
  const Int_t NBINS = 100;
  TH2F * bbc_ew = new TH2F("bbc_ew","BBCE vs BBCW Rate Correlation",NBINS,0,bbcw_rate_max,NBINS,0,bbce_rate_max);
  TH2F * bbc_ex = new TH2F("bbc_ex","BBCE vs BBCX Rate Correlation",NBINS,0,bbcx_rate_max,NBINS,0,bbce_rate_max);
  TH2F * bbc_wx = new TH2F("bbc_wx","BBCW vs BBCX Rate Correlation",NBINS,0,bbcx_rate_max,NBINS,0,bbcw_rate_max);
  TH2F * zdc_ew = new TH2F("zdc_ew","ZDCE vs ZDCW Rate Correlation",NBINS,0,zdcw_rate_max,NBINS,0,zdce_rate_max);
  TH2F * zdc_ex = new TH2F("zdc_ex","ZDCE vs ZDCX Rate Correlation",NBINS,0,zdcx_rate_max,NBINS,0,zdce_rate_max);
  TH2F * zdc_wx = new TH2F("zdc_wx","ZDCW vs ZDCX Rate Correlation",NBINS,0,zdcx_rate_max,NBINS,0,zdcw_rate_max);
  TH2F * vpd_ew = new TH2F("vpd_ew","VPDE vs VPDW Rate Correlation",NBINS,0,vpdw_rate_max,NBINS,0,vpde_rate_max);
  TH2F * vpd_ex = new TH2F("vpd_ex","VPDE vs VPDX Rate Correlation",NBINS,0,vpdx_rate_max,NBINS,0,vpde_rate_max);
  TH2F * vpd_wx = new TH2F("vpd_wx","VPDW vs VPDX Rate Correlation",NBINS,0,vpdx_rate_max,NBINS,0,vpdw_rate_max);

  TH2F * bz_e = new TH2F("bz_e","BBCE vs ZDCE Rate Correlation",NBINS,0,zdce_rate_max,NBINS,0,bbce_rate_max);
  TH2F * bz_w = new TH2F("bz_w","BBCW vs ZDCW Rate Correlation",NBINS,0,zdcw_rate_max,NBINS,0,bbcw_rate_max);
  TH2F * bz_x = new TH2F("bz_x","BBCX vs ZDCX Rate Correlation",NBINS,0,zdcx_rate_max,NBINS,0,bbcx_rate_max);
  TH2F * bv_e = new TH2F("bv_e","BBCE vs VPDE Rate Correlation",NBINS,0,vpde_rate_max,NBINS,0,bbce_rate_max);
  TH2F * bv_w = new TH2F("bv_w","BBCW vs VPDW Rate Correlation",NBINS,0,vpdw_rate_max,NBINS,0,bbcw_rate_max);
  TH2F * bv_x = new TH2F("bv_x","BBCX vs VPDX Rate Correlation",NBINS,0,vpdx_rate_max,NBINS,0,bbcx_rate_max);
  TH2F * vz_e = new TH2F("vz_e","VPDE vs ZDCE Rate Correlation",NBINS,0,zdce_rate_max,NBINS,0,vpde_rate_max);
  TH2F * vz_w = new TH2F("vz_w","VPDW vs ZDCW Rate Correlation",NBINS,0,zdcw_rate_max,NBINS,0,vpdw_rate_max);
  TH2F * vz_x = new TH2F("vz_x","VPDX vs ZDCX Rate Correlation",NBINS,0,zdcx_rate_max,NBINS,0,vpdx_rate_max);


  // fill correlation histograms
  for(Int_t i=0; i<IMAX; i++)
  {
    bbc_ew->Fill(bbcw_rate[i],bbce_rate[i]);
    bbc_ex->Fill(bbcx_rate[i],bbce_rate[i]);
    bbc_wx->Fill(bbcx_rate[i],bbcw_rate[i]);

    zdc_ew->Fill(zdcw_rate[i],zdce_rate[i]);
    zdc_ex->Fill(zdcx_rate[i],zdce_rate[i]);
    zdc_wx->Fill(zdcx_rate[i],zdcw_rate[i]);

    vpd_ew->Fill(vpdw_rate[i],vpde_rate[i]);
    vpd_ex->Fill(vpdx_rate[i],vpde_rate[i]);
    vpd_wx->Fill(vpdx_rate[i],vpdw_rate[i]);

    bz_e->Fill(zdce_rate[i],bbce_rate[i]);
    bz_w->Fill(zdcw_rate[i],bbcw_rate[i]);
    bz_x->Fill(zdcx_rate[i],bbcx_rate[i]);
    
    bv_e->Fill(vpde_rate[i],bbce_rate[i]);
    bv_w->Fill(vpdw_rate[i],bbcw_rate[i]);
    bv_x->Fill(vpdx_rate[i],bbcx_rate[i]);
    
    vz_e->Fill(zdce_rate[i],vpde_rate[i]);
    vz_w->Fill(zdcw_rate[i],vpdw_rate[i]);
    vz_x->Fill(zdcx_rate[i],vpdx_rate[i]);
  };



  // compute correlation coefficients
  char bbc_ew_text[64]; 
  char bbc_ex_text[64]; 
  char bbc_wx_text[64]; 
  char zdc_ew_text[64]; 
  char zdc_ex_text[64]; 
  char zdc_wx_text[64]; 
  char vpd_ew_text[64]; 
  char vpd_ex_text[64]; 
  char vpd_wx_text[64]; 
  char bz_e_text[64];
  char bz_w_text[64];
  char bz_x_text[64];
  char bv_e_text[64];
  char bv_w_text[64];
  char bv_x_text[64];
  char vz_e_text[64];
  char vz_w_text[64];
  char vz_x_text[64];
  sprintf(bbc_ew_text,"BBC E:W -- C=%0.5f",bbc_ew->GetCorrelationFactor());
  sprintf(bbc_ex_text,"BBC E:X -- C=%0.5f",bbc_ex->GetCorrelationFactor());
  sprintf(bbc_wx_text,"BBC W:X -- C=%0.5f",bbc_wx->GetCorrelationFactor());
  sprintf(zdc_ew_text,"ZDC E:W -- C=%0.5f",zdc_ew->GetCorrelationFactor());
  sprintf(zdc_ex_text,"ZDC E:X -- C=%0.5f",zdc_ex->GetCorrelationFactor());
  sprintf(zdc_wx_text,"ZDC W:X -- C=%0.5f",zdc_wx->GetCorrelationFactor());
  sprintf(vpd_ew_text,"VPD E:W -- C=%0.5f",vpd_ew->GetCorrelationFactor());
  sprintf(vpd_ex_text,"VPD E:X -- C=%0.5f",vpd_ex->GetCorrelationFactor());
  sprintf(vpd_wx_text,"VPD W:X -- C=%0.5f",vpd_wx->GetCorrelationFactor());
  sprintf(bz_e_text,"BBCE:ZDCE -- C=%0.5f",bz_e->GetCorrelationFactor());
  sprintf(bz_w_text,"BBCW:ZDCW -- C=%0.5f",bz_w->GetCorrelationFactor());
  sprintf(bz_x_text,"BBCX:ZDCX -- C=%0.5f",bz_x->GetCorrelationFactor());
  sprintf(bv_e_text,"BBCE:VPDE -- C=%0.5f",bv_e->GetCorrelationFactor());
  sprintf(bv_w_text,"BBCW:VPDW -- C=%0.5f",bv_w->GetCorrelationFactor());
  sprintf(bv_x_text,"BBCX:VPDX -- C=%0.5f",bv_x->GetCorrelationFactor());
  sprintf(vz_e_text,"VPDE:ZDCE -- C=%0.5f",vz_e->GetCorrelationFactor());
  sprintf(vz_w_text,"VPDW:ZDCW -- C=%0.5f",vz_w->GetCorrelationFactor());
  sprintf(vz_x_text,"VPDX:ZDCX -- C=%0.5f",vz_x->GetCorrelationFactor());

  Float_t starty=0.8;
  Float_t startx=0.3;
  Float_t interval=0.1;
  TLatex * corr_latex = new TLatex(startx-0.05,starty,"Correlation Coefficients");
  TLatex * bbc_ew_latex = new TLatex(startx,starty-1*interval,bbc_ew_text);
  TLatex * bbc_ex_latex = new TLatex(startx,starty-2*interval,bbc_ex_text);
  TLatex * bbc_wx_latex = new TLatex(startx,starty-3*interval,bbc_wx_text);
  TLatex * zdc_ew_latex = new TLatex(startx,starty-1*interval,zdc_ew_text);
  TLatex * zdc_ex_latex = new TLatex(startx,starty-2*interval,zdc_ex_text);
  TLatex * zdc_wx_latex = new TLatex(startx,starty-3*interval,zdc_wx_text);
  TLatex * vpd_ew_latex = new TLatex(startx,starty-1*interval,vpd_ew_text);
  TLatex * vpd_ex_latex = new TLatex(startx,starty-2*interval,vpd_ex_text);
  TLatex * vpd_wx_latex = new TLatex(startx,starty-3*interval,vpd_wx_text);
  TLatex * bz_e_latex = new TLatex(startx,starty-1*interval,bz_e_text);
  TLatex * bz_w_latex = new TLatex(startx,starty-2*interval,bz_w_text);
  TLatex * bz_x_latex = new TLatex(startx,starty-3*interval,bz_x_text);
  TLatex * bv_e_latex = new TLatex(startx,starty-1*interval,bv_e_text);
  TLatex * bv_w_latex = new TLatex(startx,starty-2*interval,bv_w_text);
  TLatex * bv_x_latex = new TLatex(startx,starty-3*interval,bv_x_text);
  TLatex * vz_e_latex = new TLatex(startx,starty-1*interval,vz_e_text);
  TLatex * vz_w_latex = new TLatex(startx,starty-2*interval,vz_w_text);
  TLatex * vz_x_latex = new TLatex(startx,starty-3*interval,vz_x_text);


  // draw correlation histograms
  gStyle->SetOptStat(0);
  TCanvas * canv_corr_bbc = new TCanvas("canv_corr_bbc","canv_corr_bbc",1100,700); canv_corr_bbc->Divide(2,2);
  TCanvas * canv_corr_zdc = new TCanvas("canv_corr_zdc","canv_corr_zdc",1100,700); canv_corr_zdc->Divide(2,2);
  TCanvas * canv_corr_vpd = new TCanvas("canv_corr_vpd","canv_corr_vpd",1100,700); canv_corr_vpd->Divide(2,2);
  TCanvas * canv_corr_bz = new TCanvas("canv_corr_bz","canv_corr_bz",1100,700); canv_corr_bz->Divide(2,2);
  TCanvas * canv_corr_bv = new TCanvas("canv_corr_bv","canv_corr_bv",1100,700); canv_corr_bv->Divide(2,2);
  TCanvas * canv_corr_vz = new TCanvas("canv_corr_vz","canv_corr_vz",1100,700); canv_corr_vz->Divide(2,2);

  for(Int_t i=1; i<=4; i++)
  {
    canv_corr_bbc->GetPad(i)->SetGrid(1,1);
    canv_corr_zdc->GetPad(i)->SetGrid(1,1);
    canv_corr_vpd->GetPad(i)->SetGrid(1,1);
    canv_corr_bz->GetPad(i)->SetGrid(1,1);
    canv_corr_bv->GetPad(i)->SetGrid(1,1);
    canv_corr_vz->GetPad(i)->SetGrid(1,1);
  };

  canv_corr_bbc->cd(1); bbc_ew->Draw("colz");
  canv_corr_bbc->cd(3); bbc_ex->Draw("colz");
  canv_corr_bbc->cd(4); bbc_wx->Draw("colz");
  canv_corr_zdc->cd(1); zdc_ew->Draw("colz");
  canv_corr_zdc->cd(3); zdc_ex->Draw("colz");
  canv_corr_zdc->cd(4); zdc_wx->Draw("colz");
  canv_corr_vpd->cd(1); vpd_ew->Draw("colz");
  canv_corr_vpd->cd(3); vpd_ex->Draw("colz");
  canv_corr_vpd->cd(4); vpd_wx->Draw("colz");
  canv_corr_bz->cd(1); bz_e->Draw("colz");
  canv_corr_bz->cd(3); bz_w->Draw("colz");
  canv_corr_bz->cd(4); bz_x->Draw("colz");
  canv_corr_bv->cd(1); bv_e->Draw("colz");
  canv_corr_bv->cd(3); bv_w->Draw("colz");
  canv_corr_bv->cd(4); bv_x->Draw("colz");
  canv_corr_vz->cd(1); vz_e->Draw("colz");
  canv_corr_vz->cd(3); vz_w->Draw("colz");
  canv_corr_vz->cd(4); vz_x->Draw("colz");

  canv_corr_bbc->cd(2); corr_latex->Draw(); bbc_ew_latex->Draw(); bbc_ex_latex->Draw(); bbc_wx_latex->Draw();
  canv_corr_zdc->cd(2); corr_latex->Draw(); zdc_ew_latex->Draw(); zdc_ex_latex->Draw(); zdc_wx_latex->Draw();
  canv_corr_vpd->cd(2); corr_latex->Draw(); vpd_ew_latex->Draw(); vpd_ex_latex->Draw(); vpd_wx_latex->Draw();
  canv_corr_bz->cd(2); corr_latex->Draw(); bz_e_latex->Draw(); bz_w_latex->Draw(); bz_x_latex->Draw();
  canv_corr_bv->cd(2); corr_latex->Draw(); bv_e_latex->Draw(); bv_w_latex->Draw(); bv_x_latex->Draw();
  canv_corr_vz->cd(2); corr_latex->Draw(); vz_e_latex->Draw(); vz_w_latex->Draw(); vz_x_latex->Draw();


};
