/* - computes relative luminosity for BBC,ZDC,VPD w.r.t. parameter "var"
 *   - useful "var"'s
 *     - run index "i"
 *     - fill index "fi"
 *     - bXing number "bx"
 *  - the file rdat.root is populated with TCanvases for scaler counts and rellums
 *    (see documentation for further details of what's contained in output rdat.root files)
 *    - "raw" = raw scaler counts
 *    - "acc" = accidentals corrected scaler counts (stage 1)
 *    - "mul" = multiples corrected scaler counts (stage 2)
 *    - "fac" = correction factor (mul / raw)
 *    - "rsc" = manion's rate-safe counts (written Omega * rsc, where
 *              Omega is the product of efficiency and acceptance for E & W scalers)
 *    - "rsr" = rate-safe counts ratio (Omega * rsc / raw coincidences) [cf. "fac"]
 *    - R1-9 = 9 possible relative luminosities
 *    - R*_zdc_minus_vpd = difference between rellum for zdc and vpd
 * - can write pngs to png_rellum subdirectory (run "rellum_all" to do this
 *   for the interesting independent vars)
 *
 * - printPNGs = true will produce png files (run as background!)
 * - drawLog = true will draw all scaler counts plots on log y scale
 * - zoomIn = true will "zoom in" on abort gaps scaler counts by setting
 *            a maximum in the distributions -- this is more useful for
 *            looking at structure in linear bXing distributions
 *            (only enabled if var==bx)
 * - specificFill = if 0 --> plot all fills
 *                  if >0 --> plot only fill no. "specificFill"
 * - specificRun = if 0 --> plot all runs
 *                 if >0 --> plot only run no. "specificFill"
 *
 *                 note: if both specificFill and specificRun>0, code will exit
 *                  
 */

void rellum4(const char * var="i",Bool_t printPNGs=0, 
             Bool_t drawLog=0, Int_t zoomIn=0, 
             Int_t specificFill=0, Int_t specificRun=0)
{
  // read counts.root file
  TFile * infile = new TFile("counts.root","READ");
  TTree * tr = (TTree*) infile->Get("sca");
  char outname[128];
  sprintf(outname,"rdat_%s.root",var);
  if(specificFill==0 && specificRun==0) TFile * outfile = new TFile(outname,"RECREATE");
  //tr->Print();

  if(specificFill>0 && specificRun>0)
  {
    fprintf(stderr,"ERROR: both specificFill and specificRun specified\n");
    return;
  };

  

  // zeroAborts boolean -- if true, will draw plots that lack "spin_cut" (i.e. plots
  // for all four possible spin combination; this is useful for looking at abort
  // gap scaler counts in bXing distributions; this boolean does not affect rellum calculation
  Bool_t zeroAborts;
  if(!strcmp(var,"bx")) zeroAborts=0;
  else zeroAborts=1;


  // define independent variable bounds
  Int_t var_l = tr->GetMinimum(var);
  Int_t var_h = tr->GetMaximum(var);
  var_h++; // fencepost
  Int_t var_bins = var_h-var_l;
  if(!strcmp(var,"t")) var_bins=800;
  const Int_t var_bins_const = var_bins;
  printf("var_bins=%d var_l=%d var_h=%d\n",var_bins,var_l,var_h);


  // some style variables
  const Double_t FSIZE = 0.08;
  const Double_t LWIDTH = 2;



  // ARRAY DEFINITIONS: 
  // from here on, all distributions are identified by 3-d arrays:
  // name[trigger bit][combination bit][spin bit]
  //
  // int   triggerbit     combbit     spinbit
  // ---   ----------     -------     -------
  // 0     BBC            E           B- Y-
  // 1     ZDC            W           B- Y+
  // 2     VPD            X           B+ Y-
  // 3                                B+ Y+ 
  //
  //                                  spinbit=4 = all of them



  // trigger bit character strings (tbit)
  char tbit[3][4];
  sprintf(tbit[0],"bbc");
  sprintf(tbit[1],"zdc");
  sprintf(tbit[2],"vpd");


  // combination bit character strings (cbit)
  char cbit[3][4];
  sprintf(cbit[0],"e");
  sprintf(cbit[1],"w");
  sprintf(cbit[2],"x");


  // spin bit character strings (sbit) ( p & n ... for th1 names)
  char sbit[5][4];
  sprintf(sbit[0],"nn");
  sprintf(sbit[1],"np");
  sprintf(sbit[2],"pn");
  sprintf(sbit[3],"pp");
  sprintf(sbit[4],"all");

  // spin bit character strings (nbit) ( + & - ... for th1 titles)
  char nbit[5][4];
  sprintf(nbit[0],"--");
  sprintf(nbit[1],"-+");
  sprintf(nbit[2],"+-");
  sprintf(nbit[3],"++");
  sprintf(nbit[4],"all");


  // set branch addresses
  Int_t index,runnum,fill,fi,bx;
  Double_t N[3][3]; // [tbit] [cbit]
  Double_t tot_bx,time;
  Int_t blue,yell;
  Bool_t kicked;
  tr->SetBranchAddress("i",&index);
  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("fill",&fill);
  tr->SetBranchAddress("fi",&fi);
  tr->SetBranchAddress("t",&time);
  tr->SetBranchAddress("bx",&bx);
  tr->SetBranchAddress("blue",&blue);
  tr->SetBranchAddress("yell",&yell);
  tr->SetBranchAddress("tot_bx",&tot_bx);
  // bbc
    tr->SetBranchAddress("bbce",&N[0][0]); // written explicitly for
    tr->SetBranchAddress("bbcw",&N[0][1]); // obfuscation clarifaction
    tr->SetBranchAddress("bbcx",&N[0][2]);
  // zdc 
    tr->SetBranchAddress("zdce",&N[1][0]);
    tr->SetBranchAddress("zdcw",&N[1][1]);
    tr->SetBranchAddress("zdcx",&N[1][2]); 
  // vpd
    tr->SetBranchAddress("vpde",&N[2][0]);
    tr->SetBranchAddress("vpdw",&N[2][1]);
    tr->SetBranchAddress("vpdx",&N[2][2]); 
  tr->SetBranchAddress("kicked",&kicked);



  // get run index or fill index if specific run/fill specified
  Int_t specificI=-1;
  Int_t specificFI=-1;
  Int_t specificRunMatx,specificFillMatx; // (for matx tree)
  Double_t specificT;
  fill=0;
  runnum=0;
  if(specificFill>0 || specificRun>0)
  {
    for(Int_t jj=0; jj<tr->GetEntries(); jj++)
    {
      tr->GetEntry(jj);
      if(fill==specificFill || runnum==specificRun) 
      {
        specificFI=fi;
        specificI=index;
        specificT=time;
        specificRunMatx=runnum;
        specificFillMatx=fill;
      };
    };
    if(specificFI==-1 && specificI==-1)
    {
      fprintf(stderr,"ERROR: specificRun or specificFill not in counts tr\n");
      return;
    };
  };
  if(specificFill>0) printf("fill=%d index=%d\n",specificFill,specificFI);
  else if(specificRun>0) printf("run=%d index=%d\n",specificRun,specificI);



  // build time array, fill array, runnum array  (for rate dependence; calculated iff var=="i")
  Int_t index_tmp=0;
  Int_t time_array[var_bins_const]; // fenceposting: array index = run index - 1
  Int_t fill_array[var_bins_const];
  Int_t runnum_array[var_bins_const];
  if(!strcmp(var,"i"))
  {
    for(Int_t i=0; i<tr->GetEntries(); i++)
    {
      tr->GetEntry(i);
      if(index_tmp != index)
      {
        time_array[index-1] = time;
        fill_array[index-1] = fill;
        runnum_array[index-1] = runnum;
        index_tmp = index;
      };
    };
  };


  // fill kicked bXing array
  Bool_t kicked_arr[120];
  for(Int_t i=0; i<120; i++) kicked_arr[i]=0;
  if(specificRun>0 || specificFI>0)
  {
    for(Int_t i=0; i<tr->GetEntries(); i++)
    {
      tr->GetEntry(i);
      if((index==specificI && specificRun>0) || (fi==specificFI && specificFill>0))
        kicked_arr[bx]=kicked;
    };
  };
  //for(Int_t i=0; i<120; i++) printf("bx=%d\tkicked=%d\n",i,(Int_t)(kicked_arr[bx]));


  // fill polarization data array (only if var=="i")
  // -- used in systematic error computation -- DEPRECATED
  /*
  Float_t polar_b_array[var_bins_const];
  Float_t polar_y_array[var_bins_const];
  Int_t pol_fill;
  Float_t b_pol,y_pol;
  pol_tr->SetBranchAddress("fill",&pol_fill);
  pol_tr->SetBranchAddress("b_pol",&b_pol);
  pol_tr->SetBranchAddress("y_pol",&y_pol);
  if(!strcmp(var,"i"))
  {
    for(Int_t i=0; i<pol_tr->GetEntries(); i++)
    {
      pol_tr->GetEntry(i);
      for(Int_t j=0; j<var_bins_const; j++)
      {
        if(fill_array[j] == pol_fill)
        {
          polar_b_array[j] = b_pol;
          polar_y_array[j] = y_pol;
        };
      };
    };
    //for(j=0; j<var_bins_const; j++) printf("%d %d %f %f\n",runnum_array[j],fill_array[j],polar_b_array[j],polar_y_array[j]);
  };
  */



  // define distributions
  /*   raw = raw (uncorrected) distributions
   *   acc = accidentals corrected distributions
   *   mul = multiples corrected distributions
   *   fac = correction factor (multiples / raw)
   *   rsc = Omega * rate-safe counts (Omega := acceptance * efficinecy of E & W)
   *   rsr = rate-safe counts ratio (Omega * rate-safe counts / raw coincidences)
   *   tot = tot_bx (total scaler counts) distributions
   */
  char raw_n[3][3][5][256];  // [tbit] [cbit] [sbit] [char buffer]
  char raw_t[3][3][5][256];
  char acc_n[3][3][5][256];
  char acc_t[3][3][5][256];
  char mul_n[3][3][5][256];
  char mul_t[3][3][5][256];
  char fac_n[3][3][5][256];
  char fac_t[3][3][5][256];
  char rsc_n[3][5][256]; // [tbit] [sbit] [char buffer]
  char rsc_t[3][5][256]; 
  char rsr_n[3][5][256];
  char rsr_t[3][5][256];
  char tot_n[5][256]; // [sbit] [char buffer]
  char tot_t[5][256];
  TH1D * raw_d[3][3][5]; // raw scaler counts
  TH1D * acc_d[3][3][5]; // accidentals corrected scaler counts (stage 1)
  TH1D * mul_d[3][3][5]; // multiples corrected scaler counts (stage 2)
  TH1D * fac_d[3][3][5]; // correction factor (multiples corrected / raw)
  TH1D * rsc_d[3][5]; // Omega * rate-safe counts
  TH1D * rsr_d[3][5]; // rate-safe counts ratio  ( Omega * rate-safe counts / raw coincidences)
  TH1D * tot_d[5]; // [sbit]
  char leg[50]; 
  if (zeroAborts) sprintf(leg,"(Grn:-- Orn:-+ Red:+- Blue:++)");
  else sprintf(leg,"(Grn:-- Orn:-+ Red:+- Blue:++ Blk:all)");
  char extra_t[16];
  if (specificFill==0 && specificRun==0) sprintf(extra_t,"");
  else if(specificFill>0) sprintf(extra_t," -- F%d",specificFill);
  else if(specificRun>0) sprintf(extra_t," -- R%d",specificRun);
  for(Int_t s=0; s<5; s++)
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        sprintf(raw_n[t][c][s],"raw_%s%s_%s",tbit[t],cbit[c],sbit[s]);
        sprintf(raw_t[t][c][s],"raw %s%s vs. %s %s%s",tbit[t],cbit[c],var,leg,extra_t);
        sprintf(acc_n[t][c][s],"acc_%s%s_%s",tbit[t],cbit[c],sbit[s]);
        sprintf(acc_t[t][c][s],"accidentals corrected %s%s vs %s %s%s",tbit[t],cbit[c],var,leg,extra_t);
        sprintf(mul_n[t][c][s],"mul_%s%s_%s",tbit[t],cbit[c],sbit[s]);
        sprintf(mul_t[t][c][s],"multiples corrected %s%s vs. %s %s%s",tbit[t],cbit[c],var,leg,extra_t);
        sprintf(fac_n[t][c][s],"fac_%s%s_%s",tbit[t],cbit[c],sbit[s]);
        sprintf(fac_t[t][c][s],"correction factor (mult/raw) for %s%s vs. %s %s%s",tbit[t],cbit[c],var,leg,extra_t);
        raw_d[t][c][s] = new TH1D(raw_n[t][c][s],raw_t[t][c][s],var_bins,var_l,var_h);
        acc_d[t][c][s] = new TH1D(acc_n[t][c][s],acc_t[t][c][s],var_bins,var_l,var_h);
        mul_d[t][c][s] = new TH1D(mul_n[t][c][s],mul_t[t][c][s],var_bins,var_l,var_h);
        fac_d[t][c][s] = new TH1D(fac_n[t][c][s],fac_t[t][c][s],var_bins,var_l,var_h);
      };
      sprintf(rsc_n[t][s],"rsc_%s_%s",tbit[t],sbit[s]);
      sprintf(rsc_t[t][s],"#Omega * N_{%s}^{rsc} vs. %s %s%s",tbit[t],var,leg,extra_t);
      sprintf(rsr_n[t][s],"rsr_%s_%s",tbit[t],sbit[s]);
      sprintf(rsr_t[t][s],"#Omega * N_{%s}^{rsc} / N_{%s}^{raw} vs. %s %s%s",tbit[t],tbit[t],var,leg,extra_t);
      rsc_d[t][s] = new TH1D(rsc_n[t][s],rsc_t[t][s],var_bins,var_l,var_h);
      rsr_d[t][s] = new TH1D(rsr_n[t][s],rsr_t[t][s],var_bins,var_l,var_h);
    };
    sprintf(tot_n[s],"tot_%s",sbit[s]);
    sprintf(tot_t[s],"total scaler counts N_{bx}^{%s}",nbit[s]);
    tot_d[s] = new TH1D(tot_n[s],tot_t[s],var_bins,var_l,var_h);
  };
  char spin_pat_t[64]; sprintf(spin_pat_t,"R_{spin} vs. %s",var);
  TH1D * spin_pat_same = new TH1D("spin_pat_same",spin_pat_t,var_bins,var_l,var_h);
  TH1D * spin_pat_diff = new TH1D("spin_pat_diff",spin_pat_t,var_bins,var_l,var_h);
  TH1D * spin_pat_rel = new TH1D("spin_pat_rel",spin_pat_t,var_bins,var_l,var_h);



  // define projection cuts
  char raw_cut[3][3][5][512]; // [tbit] [cbit] [sbit]
  char tot_cut[5][512];

  char spin_cut[5][128]; // (includes cut out of kicked bunches)
   sprintf(spin_cut[0],"blue==-1 && yell==-1 && !kicked");
   sprintf(spin_cut[1],"blue==-1 && yell==1 && !kicked");
   sprintf(spin_cut[2],"blue==1 && yell==-1 && !kicked");
   sprintf(spin_cut[3],"blue==1 && yell==1 && !kicked");
   sprintf(spin_cut[4],"!kicked"); // no cut (so abort gaps aren't zeroed)

  char spec_cut[128]; 
  if(specificFill>0) sprintf(spec_cut,"fill==%d",specificFill);
  else if(specificRun>0) sprintf(spec_cut,"runnum==%d",specificRun);
  else sprintf(spec_cut,"1");


  // project raw scaler data into dists
  for(Int_t s=0; s<5; s++)
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        sprintf(raw_cut[t][c][s],"%s%s*(%s&&%s)",tbit[t],cbit[c],spin_cut[s],spec_cut);
        tr->Project(raw_n[t][c][s],var,raw_cut[t][c][s]);
      };
    };
    sprintf(tot_cut[s],"tot_bx*(%s&&%s)",spin_cut[s],spec_cut);
    tr->Project(tot_n[s],var,tot_cut[s]);
  };
  tr->Project("spin_pat_same",var,"blue==yell && blue!=0 && yell!=0");
  tr->Project("spin_pat_diff",var,"blue!=yell && blue!=0 && yell!=0");
  spin_pat_rel->Divide(spin_pat_same,spin_pat_diff,1.0,1.0);



  // accidentals and multiples corrections and rate-safe corrections
  Double_t nn_raw[3][3][5]; // scaled counts     [tbit] [cbit] [sbit]
  Double_t p_scal[3][3][5]; // scale probabilities  (scaled counts / total bXings)
  Double_t p_phys[3][3][5]; // physical process probabilities  (accidentals corrections)
  Double_t nn_acc[3][3][5]; // counts corrected for accidentals (stage 1)
  Double_t nn_mul[3][3][5]; // counts corrected for multiples (stage 2)
  Double_t nn_fac[3][3][5]; // correction factor (mult/raw)
  Double_t nn_rsc[3][5]; // manion's rate-safe counts (Omega * lambda * Nbx) [tbit][sbit] (Omega=epsilon*epsilon)
  Double_t nn_rsr[3][5]; // correction factor (Omega * rcs / raw coin)
  Double_t nn_tot[5]; // total scaler counts (N_bx) 
  for(Int_t b=1; b<=var_bins; b++)
  {
    for(Int_t s=0; s<5; s++)
    {
      nn_tot[s] = tot_d[s]->GetBinContent(b);
      for(Int_t t=0; t<3; t++)
      {
        // get raw scaler counts for each trigger bit combination (E,W,X)
        for(Int_t c=0; c<3; c++)
        {
          nn_raw[t][c][s] = raw_d[t][c][s]->GetBinContent(b);
        };
        // only compute corrections for nonzero counts
        if(nn_tot[s]>0 && nn_raw[t][0][s]>0 && nn_raw[t][1][s]>0 && nn_raw[t][2][s]>0)
        {
          // compute scale probabilities
          for(Int_t c=0; c<3; c++) p_scal[t][c][s] = nn_raw[t][c][s] / nn_tot[s];

          // compute physical process probabilities
          p_phys[t][0][s] = (nn_raw[t][0][s] - nn_raw[t][2][s]) / (nn_tot[s] - nn_raw[t][1][s]);
          p_phys[t][1][s] = (nn_raw[t][1][s] - nn_raw[t][2][s]) / (nn_tot[s] - nn_raw[t][0][s]);
          p_phys[t][2][s] = (nn_raw[t][2][s] - (nn_raw[t][0][s]*nn_raw[t][1][s])/nn_tot[s]) /
                        (nn_tot[s] + nn_raw[t][2][s] - nn_raw[t][0][s] - nn_raw[t][1][s]);

          // rate-safe counts (manion's method) 
          // -- using scale probabilities for computation of rate-safe counts
          // -- if I use physical process probabilities (to correct for accidentals, sometimes the 
          //    number of rate-safe counts becomes negative, and VPD differs a lot from ZDC & BBC)
          nn_rsc[t][s] = nn_tot[s] * log( (1-p_scal[t][2][s]) / ((1-p_scal[t][0][s])*(1-p_scal[t][1][s])) );

          // fill rate-safe counts plots
          if(nn_raw[t][2][s]>0) nn_rsr[t][s] = nn_rsc[t][s] / nn_raw[t][2][s];
          else nn_rsr[t][s] = 0;
          rsc_d[t][s]->SetBinContent(b,nn_rsc[t][s]);
          rsr_d[t][s]->SetBinContent(b,nn_rsr[t][s]);

          // compute accidentals and multiples corrected counts and fill dists
          for(Int_t c=0; c<3; c++)
          {
            nn_acc[t][c][s] = p_phys[t][c][s] * nn_tot[s];
            nn_mul[t][c][s] = -1 * nn_tot[s] * log(1 - p_phys[t][c][s]);
            if(nn_raw[t][c][s]>0) nn_fac[t][c][s] = nn_mul[t][c][s] / nn_raw[t][c][s];
            else nn_fac[t][c][s] = 0;
            acc_d[t][c][s]->SetBinContent(b,nn_acc[t][c][s]);
            mul_d[t][c][s]->SetBinContent(b,nn_mul[t][c][s]);
            fac_d[t][c][s]->SetBinContent(b,nn_fac[t][c][s]);
          };
        };
      };
    };
  };

  
  // compute spinbit consistency (for var==fi only)
  // -- DEPRECATED
  TH1D * spinbit_dev[3][3];
  char spinbit_dev_t[3][3][64];
  char spinbit_dev_n[3][3][64];
  Double_t spinbit_mean;
  Double_t dev_calc;
  if(!strcmp(var,"fi") || !strcmp(var,"i"))
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        sprintf(spinbit_dev_t[t][c],"%s%s spinbit deviation vs. %s",tbit[t],cbit[c],var);
        sprintf(spinbit_dev_n[t][c],"spinbit_dev_%d_%d",t,c);
        spinbit_dev[t][c] = new TH1D(spinbit_dev_n[t][c],spinbit_dev_t[t][c],var_bins,var_l,var_h);
        for(Int_t b=1; b<=raw_d[t][c][0]->GetNbinsX(); b++)
        {
          spinbit_mean=0;
          for(Int_t s=0; s<4; s++)
            spinbit_mean += raw_d[t][c][s]->GetBinContent(b);
          spinbit_mean /= 4.0;
          dev_calc = 0;
          for(Int_t s=0; s<4; s++)
          {
            dev_calc += pow(raw_d[t][c][s]->GetBinContent(b) - spinbit_mean, 2); 
          };
          dev_calc = sqrt(dev_calc/3.0);
          spinbit_dev[t][c]->SetBinContent(b,dev_calc);
        };
      };
    };
  };



  // statistical uncertainties for acc+mul corrected counts -- ( just sqrt(counts) statistics)
  Double_t BC;
  Double_t unc;
  TH1D * mul_unc_d[3][3][4]; // [tbit] [cbit] [spinbit] -- mul count uncertainty
  char mul_unc_t[3][3][4][256];
  char mul_unc_n[3][3][4][32];
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      for(Int_t s=0; s<4; s++)
      {
        sprintf(mul_unc_t[t][c][s],"%s%s mul uncertainty (spin %s) vs. %s",tbit[t],cbit[c],nbit[s],var);
        sprintf(mul_unc_n[t][c][s],"mul_unc_%s%s_s%d",tbit[t],cbit[c],s);
        mul_unc_d[t][c][s] = new TH1D(mul_unc_n[t][c][s],mul_unc_t[t][c][s],var_bins,var_l,var_h);

        for(Int_t b=1; b<=var_bins; b++)
        {
          BC = mul_d[t][c][s]->GetBinContent(b);
          unc = sqrt(BC);
          mul_unc_d[t][c][s]->SetBinContent(b,unc);
        };
      };
    };
  };


  // staatistical uncertaintites for rate-safe corrected counts -- FORMULA
  Double_t RC[3]; // [cbit] // total number of raw counts
  Double_t TC; // total number of bXings
  Double_t zeta[3]; // [cbit] // zeta[c] := 1-RC[c]/TC
  Double_t corr_xe,corr_xw,corr_we; // pearson correlation coefficient
  TH1D * rsc_unc_d[3][4]; // [tbit] [spinbit] -- rate-safe count uncertainty for each spinbit
  char rsc_unc_t[3][4][256];
  char rsc_unc_n[3][4][32];
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t s=0; s<4; s++)
    {
      sprintf(rsc_unc_t[t][s],"%s rsc uncertainty (spin %s) vs. %s",tbit[t],nbit[s],var);
      sprintf(rsc_unc_n[t][s],"rsc_unc_%s_s%d",tbit[t],s);
      rsc_unc_d[t][s] = new TH1D(rsc_unc_n[t][s],rsc_unc_t[t][s],var_bins,var_l,var_h);

      for(Int_t b=1; b<=var_bins; b++)
      {
        TC = tot_d[s]->GetBinContent(b);
        for(Int_t c=0; c<3; c++)
        {
          RC[c] = raw_d[t][c][s]->GetBinContent(b);
          zeta[c] = 1-RC[c]/TC;
        };

        // pearson correlation coefficients --- IS THIS ASSUMPTION VALID !?!?!?!?!?!?! (see CheckCorrelations.C --> corr.root ....)
        corr_xe = 1; 
        corr_xw = 1;
        corr_we = 1; 

        unc = RC[0]/zeta[0] + RC[1]/zeta[1] + RC[2]/zeta[2]
              - 2 * corr_xe * sqrt( RC[2]*RC[0] / (zeta[2]*zeta[0]) )
              - 2 * corr_xw * sqrt( RC[2]*RC[1] / (zeta[2]*zeta[1]) )
              + 2 * corr_we * sqrt( RC[1]*RC[0] / (zeta[1]*zeta[0]) );
        if(unc<0) printf("==================== WARNING: RSC Im(unc)!=0\n");
        unc = sqrt(unc);
        rsc_unc_d[t][s]->SetBinContent(b,unc);
      };
    };
  };


  /*
  TH1D * Rerr_mul_d[3][3][10]; // [tbit] [cbit] [rellum] -- rellum error for that using mul counts
  char Rerr_mul_t[3][3][10][256];
  char Rerr_mul_n[3][3][10][256];
  Double_t LL[4]; // number of mul counts for each spinbit
  for(Int_t r=1; r<10; r++)
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        sprintf(Rerr_mul_t[t][c][r],"R%d error %s%s using mul vs. %s",r,tbit[t],cbit[c],var);
        sprintf(Rerr_mul_n[t][c][r],"Rerr_mul_%s%s_R%d",tbit[t],cbit[c],r);
        Rerr_mul_d[t][c][r] = new TH1D(Rerr_mul_n[t][c][r],Rerr_mul_t[t][c][r],var_bins,var_l,var_h);

        for(Int_t b=1; b<=var_bins; b++)
        {
          for(Int_t s=0; s<4; s++) LL[s] = mul_d[t][c][s]->GetBinContent(b);

          // rellum uncertainty propagation for acc+mul corrections -- FORMULA; assumes error on counts is sqrt(counts) 
          if(r==1)
            unc = ( (LL[1] + LL[3]) * (LL[0] + LL[1] + LL[2] + LL[3]) ) / pow((LL[0] + LL[2]), 3);
          else if(r==2)
            unc = ( (LL[2] + LL[3]) * (LL[0] + LL[1] + LL[2] + LL[3]) ) / pow((LL[0] + LL[1]), 3);
          else if(r==3)
            unc = ( (LL[0] + LL[3]) * (LL[0] + LL[1] + LL[2] + LL[3]) ) / pow((LL[1] + LL[2]), 3);
          else if(r==4) unc = ( LL[3] * (LL[0]+LL[3])) / pow(LL[0], 3);
          else if(r==5) unc = ( LL[1] * (LL[0]+LL[1])) / pow(LL[0], 3);
          else if(r==6) unc = ( LL[2] * (LL[0]+LL[2])) / pow(LL[0], 3);
          else if(r==7) unc = ( LL[3] * (LL[2]+LL[3])) / pow(LL[2], 3);
          else if(r==8) unc = ( LL[1] * (LL[1]+LL[2])) / pow(LL[2], 3);
          else if(r==9) unc = ( LL[3] * (LL[1]+LL[3])) / pow(LL[1], 3);

          unc = sqrt(unc);
          //unc = sqrt(fabs(unc)); // TESTING; only needed if counts go negative from "bad" corrections

          Rerr_mul_d[t][c][r]->SetBinContent(b,unc);
        };
      };
    };
  };
  */

  // propagate counts uncertainties to rellum errors
  TH1D * Rerr_mul_d[3][3][10]; // [tbit] [cbit] [rellum] -- rellum error for that using mul counts
  char Rerr_mul_t[3][3][10][256];
  char Rerr_mul_n[3][3][10][64];
  TH1D * Rerr_rsc_d[3][10]; // [tbit] [rellum] -- rellum error for that using rsc counts
  char Rerr_rsc_t[3][10][256];
  char Rerr_rsc_n[3][10][64];
  Double_t LL[4]; 
  Double_t SS[4];
  for(Int_t r=1; r<10; r++)
  {
    for(Int_t t=0; t<3; t++)
    {
      // loop over cbit one extra time; if c==3 then this loop does the uncertainty propagation
      // using rate-safe corrected counts uncertainties; if c<3, then we use sqrt(mul) uncertainties
      for(Int_t c=0; c<4; c++)
      {
        printf("r%d t%d c%d\n",r,t,c);
        if(c<3)
        {
          sprintf(Rerr_mul_t[t][c][r],"R%d error %s%s using mul vs. %s",r,tbit[t],cbit[c],var);
          sprintf(Rerr_mul_n[t][c][r],"Rerr_mul_%s%s_R%d",tbit[t],cbit[c],r);
          Rerr_mul_d[t][c][r] = new TH1D(Rerr_mul_n[t][c][r],Rerr_mul_t[t][c][r],var_bins,var_l,var_h);
        }
        else
        {
          sprintf(Rerr_rsc_t[t][r],"R%d error %s using rsc vs. %s",r,tbit[t],var);
          sprintf(Rerr_rsc_n[t][r],"Rerr_rsc_%s_R%d",tbit[t],r);
          Rerr_rsc_d[t][r] = new TH1D(Rerr_rsc_n[t][r],Rerr_rsc_t[t][r],var_bins,var_l,var_h);
        };
        

        for(Int_t b=1; b<=var_bins; b++)
        {
          for(Int_t s=0; s<4; s++)
          {
            if(c<3)
            {
              LL[s] = mul_d[t][c][s]->GetBinContent(b);
              SS[s] = mul_unc_d[t][c][s]->GetBinContent(b);
            }
            else
            {
              LL[s] = rsc_d[t][s]->GetBinContent(b);
              SS[s] = rsc_unc_d[t][s]->GetBinContent(b);
            };
            //printf("t%d c%d s%d LL=%f SS=%f\n",t,c,s,LL[s],SS[s]);
          };

          // rellum uncertainty propagation -- FORMULA; sqrt taken after if chain
          if(r==1)
            unc = ( ( pow(SS[1],2) + pow(SS[3],2) ) * pow(LL[0] + LL[2],2) + 
                    ( pow(SS[0],2) + pow(SS[2],2) ) * pow(LL[1] + LL[3],2) ) / 
                    ( pow(LL[0] + LL[2],4) );
          else if(r==2)
            unc = ( ( pow(SS[2],2) + pow(SS[3],2) ) * pow(LL[0] + LL[1],2) + 
                    ( pow(SS[0],2) + pow(SS[1],2) ) * pow(LL[2] + LL[3],2) ) / 
                    ( pow(LL[0] + LL[1],4) );
          else if(r==3)
            unc = ( ( pow(SS[0],2) + pow(SS[3],2) ) * pow(LL[1] + LL[2],2) + 
                    ( pow(SS[1],2) + pow(SS[2],2) ) * pow(LL[0] + LL[3],2) ) / 
                    ( pow(LL[1] + LL[2],4) );
          else if(r==4) 
            unc = ( pow(SS[3],2) * pow(LL[0],2) +
                    pow(SS[0],2) * pow(LL[3],2) ) /
                  ( pow(LL[0],4) );
          else if(r==5)
            unc = ( pow(SS[1],2) * pow(LL[0],2) +
                    pow(SS[0],2) * pow(LL[1],2) ) /
                  ( pow(LL[0],4) );
          else if(r==6)
            unc = ( pow(SS[2],2) * pow(LL[0],2) +
                    pow(SS[0],2) * pow(LL[2],2) ) /
                  ( pow(LL[0],4) );
          else if(r==7)
            unc = ( pow(SS[3],2) * pow(LL[2],2) +
                    pow(SS[2],2) * pow(LL[3],2) ) /
                  ( pow(LL[2],4) );
          else if(r==8)
            unc = ( pow(SS[2],2) * pow(LL[1],2) +
                    pow(SS[1],2) * pow(LL[2],2) ) /
                  ( pow(LL[2],4) );
          else if(r==9)
            unc = ( pow(SS[3],2) * pow(LL[1],2) +
                    pow(SS[1],2) * pow(LL[3],2) ) /
                  ( pow(LL[1],4) );

          unc = sqrt(unc);
          //unc = sqrt(fabs(unc)); // TESTING; only needed if counts go negative from "bad" corrections

          if(c<3) Rerr_mul_d[t][c][r]->SetBinContent(b,unc);
          else Rerr_rsc_d[t][r]->SetBinContent(b,unc);
        };
      };
    };
  };



  // set colours and font sizes
  //   -- cool colours (green & blue) for ++ & --
  //   -- warm colours (red & orange) for +- & -+
  //   IF COLORS CHANGED HERE, NEED TO CHANGE TITLES TOO
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      if(!strcmp(var,"bx"))
      {
        for(Int_t s=0; s<4; s++)
        {
          raw_d[t][c][s]->SetFillColor(kBlue);
          acc_d[t][c][s]->SetFillColor(kBlue);
          mul_d[t][c][s]->SetFillColor(kBlue);
          fac_d[t][c][s]->SetFillColor(kBlue);
          raw_d[t][c][s]->SetLineColor(kBlue);
          acc_d[t][c][s]->SetLineColor(kBlue);
          mul_d[t][c][s]->SetLineColor(kBlue);
          fac_d[t][c][s]->SetLineColor(kBlue);
          if(c==0)
          {
            rsc_d[t][s]->SetFillColor(kBlue);
            rsc_d[t][s]->SetLineColor(kBlue);
            rsr_d[t][s]->SetFillColor(kBlue);
            rsr_d[t][s]->SetLineColor(kBlue);
          };
        };
        raw_d[t][c][4]->SetFillColor(kBlack);
        acc_d[t][c][4]->SetFillColor(kBlack);
        mul_d[t][c][4]->SetFillColor(kBlack);
        fac_d[t][c][4]->SetFillColor(kBlack);
        raw_d[t][c][4]->SetLineColor(kBlack);
        acc_d[t][c][4]->SetLineColor(kBlack);
        mul_d[t][c][4]->SetLineColor(kBlack);
        fac_d[t][c][4]->SetLineColor(kBlack);
        if(c==0)
        {
          rsc_d[t][s]->SetFillColor(kBlack);
          rsc_d[t][s]->SetLineColor(kBlack);
          rsr_d[t][s]->SetFillColor(kBlack);
          rsr_d[t][s]->SetLineColor(kBlack);
        };
      }
      else
      {
        raw_d[t][c][0]->SetLineColor(kGreen+2);
        raw_d[t][c][1]->SetLineColor(kOrange+7);
        raw_d[t][c][2]->SetLineColor(kRed);
        raw_d[t][c][3]->SetLineColor(kBlue);
        raw_d[t][c][4]->SetLineColor(kBlack);
        acc_d[t][c][0]->SetLineColor(kGreen+2);
        acc_d[t][c][1]->SetLineColor(kOrange+7);
        acc_d[t][c][2]->SetLineColor(kRed);
        acc_d[t][c][3]->SetLineColor(kBlue);
        acc_d[t][c][4]->SetLineColor(kBlack);
        mul_d[t][c][0]->SetLineColor(kGreen+2);
        mul_d[t][c][1]->SetLineColor(kOrange+7);
        mul_d[t][c][2]->SetLineColor(kRed);
        mul_d[t][c][3]->SetLineColor(kBlue);
        mul_d[t][c][4]->SetLineColor(kBlack);
        fac_d[t][c][0]->SetLineColor(kGreen+2);
        fac_d[t][c][1]->SetLineColor(kOrange+7);
        fac_d[t][c][2]->SetLineColor(kRed);
        fac_d[t][c][3]->SetLineColor(kBlue);
        fac_d[t][c][4]->SetLineColor(kBlack);
        if(c==0)
        {
          rsc_d[t][0]->SetLineColor(kGreen+2);
          rsc_d[t][1]->SetLineColor(kOrange+7);
          rsc_d[t][2]->SetLineColor(kRed);
          rsc_d[t][3]->SetLineColor(kBlue);
          rsc_d[t][4]->SetLineColor(kBlack);
          rsr_d[t][0]->SetLineColor(kGreen+2);
          rsr_d[t][1]->SetLineColor(kOrange+7);
          rsr_d[t][2]->SetLineColor(kRed);
          rsr_d[t][3]->SetLineColor(kBlue);
          rsr_d[t][4]->SetLineColor(kBlack);
        };
      };
      for(Int_t s=0; s<5; s++)
      {
        raw_d[t][c][s]->GetXaxis()->SetLabelSize(FSIZE);
        raw_d[t][c][s]->GetYaxis()->SetLabelSize(FSIZE);
        acc_d[t][c][s]->GetXaxis()->SetLabelSize(FSIZE);
        acc_d[t][c][s]->GetYaxis()->SetLabelSize(FSIZE);
        mul_d[t][c][s]->GetXaxis()->SetLabelSize(FSIZE);
        mul_d[t][c][s]->GetYaxis()->SetLabelSize(FSIZE);
        fac_d[t][c][s]->GetXaxis()->SetLabelSize(FSIZE);
        fac_d[t][c][s]->GetYaxis()->SetLabelSize(FSIZE);
        if(c==0)
        {
          rsc_d[t][s]->GetXaxis()->SetLabelSize(FSIZE);
          rsc_d[t][s]->GetYaxis()->SetLabelSize(FSIZE);
          rsr_d[t][s]->GetXaxis()->SetLabelSize(FSIZE);
          rsr_d[t][s]->GetYaxis()->SetLabelSize(FSIZE);
        };
      };
    };
  };


  // zoom in on abort gaps (zoomIn)
  // implemented to see if there's fill-dependent afterpulsing effect
  // in abort gaps when viewing bXing distributions on linear scale
  Double_t raw_zoom[3][3]; // [tbit] [cbit]
  Double_t acc_zoom[3][3];
  Double_t mul_zoom[3][3];
  Double_t fac_zoom[3][3];
  Double_t rsc_zoom[3]; // [tbit]
  Double_t rsr_zoom[3];
  Double_t raw_bcc;
  Double_t acc_bcc;
  Double_t mul_bcc;
  Double_t fac_bcc;
  Double_t rsc_bcc;
  Double_t rsr_bcc;
  Int_t bset_l[2] = {32,112};
  Int_t bset_h[2] = {40,120}; 
  if(zoomIn && !strcmp(var,"bx"))
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        raw_zoom[t][c]=0;
        acc_zoom[t][c]=0;
        mul_zoom[t][c]=0;
        fac_zoom[t][c]=0;
        if(c==2)
        {
          rsc_zoom[t]=0;
          rsr_zoom[t]=0;
        };
        // loop over abort gap bXings 32-40 & 112-120
        for(Int_t bbb=0; bbb<2; bbb++)
        {
          for(Int_t b=bset_l[bbb]; b<=bset_h[bbb]; b++)
          {
            raw_bcc = raw_d[t][c][4]->GetBinContent(b);
            acc_bcc = acc_d[t][c][4]->GetBinContent(b);
            mul_bcc = mul_d[t][c][4]->GetBinContent(b);
            fac_bcc = fac_d[t][c][4]->GetBinContent(b);
            raw_zoom[t][c] = (raw_bcc > raw_zoom[t][c]) ? raw_bcc : raw_zoom[t][c];
            acc_zoom[t][c] = (acc_bcc > acc_zoom[t][c]) ? acc_bcc : acc_zoom[t][c];
            mul_zoom[t][c] = (mul_bcc > mul_zoom[t][c]) ? mul_bcc : mul_zoom[t][c];
            fac_zoom[t][c] = (fac_bcc > fac_zoom[t][c]) ? fac_bcc : fac_zoom[t][c];
            if(c==2)
            {
              rsc_bcc = rsc_d[t][4]->GetBinContent(b);
              rsr_bcc = rsr_d[t][4]->GetBinContent(b);
              rsc_zoom[t] = (rsc_bcc > rsc_zoom[t]) ? rsc_bcc : rsc_zoom[t];
              rsr_zoom[t] = (rsr_bcc > rsr_zoom[t]) ? rsr_bcc : rsr_zoom[t];
            };
          };
        };
        raw_zoom[t][c]*=1.1;
        acc_zoom[t][c]*=1.1;
        mul_zoom[t][c]*=1.1;
        fac_zoom[t][c]*=1.1;
        if(c==2)
        {
          rsc_zoom[t]*=1.1;
          rsr_zoom[t]*=1.1;
        };
        for(Int_t s=0; s<5; s++)
        {
          raw_d[t][c][s]->SetMaximum(raw_zoom[t][c]);
          acc_d[t][c][s]->SetMaximum(acc_zoom[t][c]);
          mul_d[t][c][s]->SetMaximum(mul_zoom[t][c]);
          fac_d[t][c][s]->SetMaximum(fac_zoom[t][c]);
          if(c==2)
          {
            rsc_d[t][s]->SetMaximum(rsc_zoom[t]);
            rsr_d[t][s]->SetMaximum(rsr_zoom[t]);
          };
        };
      };
    };
  };



  // compute relative luminosities using multiples corrected counts
  // -- for all arrays, we only fill entries 1-9 (skip 0th entry)
  //    so that the index matches rellum definition
  // R1 = (N++ + N-+) / (N+- + N--)
  // R2 = (N++ + N+-) / (N-+ + N--)
  // R3 = (N++ + N--) / (N+- + N-+)
  // R4 = N++ / N--
  // R5 = N-+ / N--
  // R6 = N+- / N--
  // R7 = N++ / N+-
  // R8 = N-+ / N+-
  // R9 = N++ / N-+
  TH1D * R_mul_d[3][3][10]; // [tbit] [cbit] [rellum] // relative luminosity using multiples corrections
  char R_mul_n[3][3][10][256];
  char R_mul_t[3][3][10][256];
  TH1D * R_rsc_d[3][10]; // [tbit] [rellum] // relative luminosity using rate-safe corrections
  char R_rsc_n[3][10][256];
  char R_rsc_t[3][10][256];
  Double_t mmm[4]; // [sbit] // mul counts
  Double_t rrr[10]; // [rellum]
  char rellum_equ[10][128];
    strcpy(rellum_equ[1],"R1 = (N++ + N-+) / (N+- + N--)");
    strcpy(rellum_equ[2],"R2 = (N++ + N+-) / (N-+ + N--)");
    strcpy(rellum_equ[3],"R3 = (N++ + N--) / (N+- + N-+)");
    strcpy(rellum_equ[4],"R4 = N++ / N--");
    strcpy(rellum_equ[5],"R5 = N-+ / N--");
    strcpy(rellum_equ[6],"R6 = N+- / N--");
    strcpy(rellum_equ[7],"R7 = N++ / N+-");
    strcpy(rellum_equ[8],"R8 = N-+ / N+-");
    strcpy(rellum_equ[9],"R9 = N++ / N-+");
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      // histogram initialisation
      for(Int_t r=1; r<=9; r++) 
      {
        sprintf(R_mul_t[t][c][r],"%s for %s%s (multiples corrected) vs. %s",rellum_equ[r],tbit[t],cbit[c],var);
        sprintf(R_mul_n[t][c][r],"R%d_mul_%s%s",r,tbit[t],cbit[c]);
        R_mul_d[t][c][r] = new TH1D(R_mul_n[t][c][r],R_mul_t[t][c][r],var_bins,var_l,var_h);
        if(c==2)
        {
          sprintf(R_rsc_t[t][r],"%s for %s (rate-safe corrected) vs. %s",rellum_equ[r],tbit[t],var);
          sprintf(R_rsc_n[t][r],"R%d_rsc_%s",r,tbit[t]);
          R_rsc_d[t][r] = new TH1D(R_rsc_n[t][r],R_rsc_t[t][r],var_bins,var_l,var_h);
        };
      };
      // compute relative luminosities, first using multiples (uu=0), then using rsc (uu=1)
      for(Int_t uu=0; uu<2; uu++)
      {
        for(Int_t b=1; b<=var_bins; b++)
        {
          for(Int_t s=0; s<5; s++) 
          {
            if(uu==0) mmm[s] = mul_d[t][c][s]->GetBinContent(b);
            else if(uu==1) mmm[s] = rsc_d[t][s]->GetBinContent(b);
            nn_tot[s] = tot_d[s]->GetBinContent(b);
          };
          if(nn_tot[0]*nn_tot[1]*nn_tot[2]*nn_tot[3]>0)
          {
            rrr[1] = (mmm[3] + mmm[1]) / (mmm[2] + mmm[0]);
            rrr[2] = (mmm[3] + mmm[2]) / (mmm[1] + mmm[0]);
            rrr[3] = (mmm[3] + mmm[0]) / (mmm[2] + mmm[1]);
            rrr[4] = mmm[3] / mmm[0];
            rrr[5] = mmm[1] / mmm[0];
            rrr[6] = mmm[2] / mmm[0];
            rrr[7] = mmm[3] / mmm[2];
            rrr[8] = mmm[1] / mmm[2];
            rrr[9] = mmm[3] / mmm[1];
            for(Int_t r=1; r<=9; r++) 
            {
              if(uu==0)
              {
                R_mul_d[t][c][r]->SetBinContent(b,rrr[r]);
                R_mul_d[t][c][r]->SetBinError(b,Rerr_mul_d[t][c][r]->GetBinContent(b));
              }
              else if(uu==1 && c==2)
              {
                R_rsc_d[t][r]->SetBinContent(b,rrr[r]);
                R_rsc_d[t][r]->SetBinError(b,Rerr_rsc_d[t][r]->GetBinContent(b));
              };
            };
          };
        };
      };
    };
  };




  // COMPARISONS AND AVERAGES OF RELLUMS --------------------------------------
  // means of rellum over {E,W,X}  (not done for rsc, since rsc is only one value for cbit==2)
  char mean_R_n[3][10][128]; // [tbit] [rellum]
  char mean_R_t[3][10][128];
  TH1D * mean_R[3][10];
  Double_t ave;
  Float_t unc_b[3];
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      sprintf(mean_R_n[t][r],"mean_R%d_%s",r,tbit[t]);
      sprintf(mean_R_t[t][r],"mean R%d via %s (multiples corrected) over e,w,x",r,tbit[t]);
      mean_R[t][r] = new TH1D(mean_R_n[t][r],mean_R_t[t][r],var_bins,var_l,var_h);
      for(Int_t b=1; b<=var_bins; b++)
      {
        ave=0;
        for(Int_t c=0; c<3; c++) ave += R_mul_d[t][c][r]->GetBinContent(b);
        ave /= 3.0;
        mean_R[t][r]->SetBinContent(b,ave);

        // propagate error bars into mean -- FORMULA
        unc_b[0] = Rerr_mul_d[t][0][r]->GetBinContent(b);
        unc_b[1] = Rerr_mul_d[t][1][r]->GetBinContent(b);
        unc_b[2] = Rerr_mul_d[t][2][r]->GetBinContent(b);
        unc = 1/3.0 * sqrt( pow(unc_b[0],2) + pow(unc_b[1],2) + pow(unc_b[2],2) );
        mean_R[t][r]->SetBinError(b,unc);
      };
    };
  };

  // deviations from mean of R*
  char dev_R_n[3][3][10][128]; // [tbit] [cbit] [rellum]
  char dev_R_t[3][3][10][256]; 
  TH1D * dev_R[3][3][10];
  Double_t dev;
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      for(Int_t r=1; r<10; r++)
      {
        sprintf(dev_R_n[t][c][r],"dev_R%d_%s%s",r,tbit[t],cbit[c]);
        sprintf(dev_R_t[t][c][r],"deviation = R%d - mean R%d for %s%s (multiples corrected)",r,r,tbit[t],cbit[c]);
        dev_R[t][c][r] = new TH1D(dev_R_n[t][c][r],dev_R_t[t][c][r],var_bins,var_l,var_h);
        dev_R[t][c][r]->Add(R_mul_d[t][c][r],mean_R[t][r],1.0,-1.0);
      };
    };
  };
  

  // compare zdc and vpd
  char D_mul_n[3][10][128]; // [cbit] [rellum]
  char D_mul_t[3][10][256];
  TH1D * D_mul_d[3][10];
  char D_rsc_n[10][128]; // [rellum]
  char D_rsc_t[10][256];
  TH1D * D_rsc_d[10];
  for(Int_t c=0; c<3; c++)
  {
    for(Int_t r=1; r<10; r++)
    {
      sprintf(D_mul_n[c][r],"delta_mul_zdc%s_vpd%s_%d",cbit[c],cbit[c],r);
      sprintf(D_mul_t[c][r],"R%d(zdc%s) minus R%d(vpd%s) via multiples corrections vs. %s",r,cbit[c],r,cbit[c],var);
      D_mul_d[c][r] = new TH1D(D_mul_n[c][r],D_mul_t[c][r],var_bins,var_l,var_h);
      D_mul_d[c][r]->Add(R_mul_d[1][c][r],R_mul_d[2][c][r],1.0,-1.0);
    };
  };
  for(Int_t r=1; r<10; r++)
  {
    sprintf(D_rsc_n[r],"delta_rsc_zdc_vpd_%d",r);
    sprintf(D_rsc_t[r],"R%d(zdc) minus R%d(vpd) via rate-safe corrections vs. %s",r,r,var);
    D_rsc_d[r] = new TH1D(D_rsc_n[r],D_rsc_t[r],var_bins,var_l,var_h);
    D_rsc_d[r]->Add(R_rsc_d[1][r],R_rsc_d[2][r],1.0,-1.0);
  };



  // compare singles bit combinations (east minus west, east minus coin, west minus coin)
  // xbit definition: 0=E-W, 1=E-X, 2=W-X
  char SD_n[3][3][10][128]; // [xbit] [tbit] [rellum]
  char SD_t[3][3][10][256];
  TH1D * SD_d[3][3][10]; 
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      sprintf(SD_n[0][t][r],"e_minus_w_diff_%s_%d",tbit[t],r);
      sprintf(SD_n[1][t][r],"e_minus_x_diff_%s_%d",tbit[t],r);
      sprintf(SD_n[2][t][r],"w_minus_x_diff_%s_%d",tbit[t],r);

      sprintf(SD_t[0][t][r],"R%d(%se)-R%d(%sw) vs. %s",r,tbit[t],r,tbit[t],var);
      sprintf(SD_t[1][t][r],"R%d(%se)-R%d(%sx) vs. %s",r,tbit[t],r,tbit[t],var);
      sprintf(SD_t[2][t][r],"R%d(%sw)-R%d(%sx) vs. %s",r,tbit[t],r,tbit[t],var);

      for(Int_t x=0; x<3; x++) SD_d[x][t][r] = new TH1D(SD_n[x][t][r],SD_t[x][t][r],var_bins,var_l,var_h);
      
      SD_d[0][t][r]->Add(R_mul_d[t][0][r],R_mul_d[t][1][r],1.0,-1.0); // E-W
      SD_d[1][t][r]->Add(R_mul_d[t][0][r],R_mul_d[t][2][r],1.0,-1.0); // E-X
      SD_d[2][t][r]->Add(R_mul_d[t][1][r],R_mul_d[t][2][r],1.0,-1.0); // W-X
    };
  };


  // compare rate-safe corrections to multiples corrections on coincidences
  char RD_n[3][10][128]; // [tbit] [rellum]
  char RD_t[3][10][256];
  TH1D * RD_d[3][10];
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      sprintf(RD_n[t][r],"rsc_minus_mul_%s_%d",tbit[t],r);
      sprintf(RD_t[t][r],"R%d(%s,rsc) - R%d(%sx,mul) vs. %s",r,tbit[t],r,tbit[t],var);
      RD_d[t][r] = new TH1D(RD_n[t][r],RD_t[t][r],var_bins,var_l,var_h);
      RD_d[t][r]->Add(R_rsc_d[t][r],R_mul_d[t][2][r],1.0,-1.0);
    };
  };



  // set more font sizes
  for(Int_t r=1; r<10; r++)
  {
    for(Int_t c=0; c<3; c++)
    {
      D_mul_d[c][r]->GetXaxis()->SetLabelSize(FSIZE);
      D_mul_d[c][r]->GetYaxis()->SetLabelSize(FSIZE);
      D_mul_d[c][r]->SetLineWidth(LWIDTH);
      if(c==2)
      {
        D_rsc_d[r]->GetXaxis()->SetLabelSize(FSIZE);
        D_rsc_d[r]->GetYaxis()->SetLabelSize(FSIZE);
        D_rsc_d[r]->SetLineWidth(LWIDTH);
      };
      for(Int_t t=0; t<3; t++)
      {
        R_mul_d[t][c][r]->GetXaxis()->SetLabelSize(FSIZE);
        R_mul_d[t][c][r]->GetYaxis()->SetLabelSize(FSIZE);
        R_mul_d[t][c][r]->SetLineWidth(LWIDTH);
        if(c==2)
        {
          R_rsc_d[t][r]->GetXaxis()->SetLabelSize(FSIZE);
          R_rsc_d[t][r]->GetYaxis()->SetLabelSize(FSIZE);
          R_rsc_d[t][r]->SetLineWidth(LWIDTH);
        };
        dev_R[t][c][r]->GetXaxis()->SetLabelSize(FSIZE);
        dev_R[t][c][r]->GetYaxis()->SetLabelSize(FSIZE);
        dev_R[t][c][r]->SetLineWidth(LWIDTH);
      };
    };
    for(Int_t t=0; t<3; t++)
    {
      mean_R[t][r]->GetXaxis()->SetLabelSize(FSIZE);
      mean_R[t][r]->GetYaxis()->SetLabelSize(FSIZE);
      mean_R[t][r]->SetLineWidth(LWIDTH);
      RD_d[t][r]->GetXaxis()->SetLabelSize(FSIZE);
      RD_d[t][r]->GetYaxis()->SetLabelSize(FSIZE);
      RD_d[t][r]->SetLineWidth(LWIDTH);
      for(Int_t x=0; x<3; x++)
      {
        SD_d[x][t][r]->GetXaxis()->SetLabelSize(FSIZE);
        SD_d[x][t][r]->GetYaxis()->SetLabelSize(FSIZE);
        SD_d[x][t][r]->SetLineWidth(LWIDTH);
      };
    };
  };


  // constant fits
  for(Int_t r=1; r<10; r++)
  {
    for(Int_t c=0; c<3; c++)
    {
      D_mul_d[c][r]->Fit("pol0","Q","",var_l,var_h);
      if(c==2) D_rsc_d[r]->Fit("pol0","Q","",var_l,var_h);
      for(Int_t t=0; t<3; t++)
      {
        R_mul_d[t][c][r]->Fit("pol0","Q","",var_l,var_h);
        if(c==2) R_rsc_d[t][r]->Fit("pol0","Q","",var_l,var_h);
        dev_R[t][c][r]->Fit("pol0","Q","",var_l,var_h);
      };
    };
    for(Int_t t=0; t<3; t++)
    {
      mean_R[t][r]->Fit("pol0","Q","",var_l,var_h);
      RD_d[t][r]->Fit("pol0","Q","",var_l,var_h);
      for(Int_t x=0; x<3; x++) SD_d[x][t][r]->Fit("pol0","Q","",var_l,var_h);
    };
  };


  
  // compute rate dependence of correction factors and "compare" plots (comparing mul, raw, and rsc)
  // !!!! computation runs iff var=="i"
  TH2D * rate_dep[3][3][10]; // [tbit] [cbit] [rellum] -- R3 vs. rate
  char rate_dep_n[3][3][10][128];
  char rate_dep_t[3][3][10][256];
  TH2D * rate_fac[3][3][5]; // [tbit] [cbit] [sbit] -- correction factor (mul/raw) vs. rate
  char rate_fac_n[3][3][5][128];
  char rate_fac_t[3][3][5][256];
  TProfile * rate_fac_pfx[3][3][5];
  TH2D * rate_rsr[3][5]; // [tbit] [sbit] -- (omega * rsc / rawx) vs. rate  (rawx := raw coincidences)
  char rate_rsr_n[3][5][256];
  char rate_rsr_t[3][5][512];
  TProfile * rate_rsr_pfx[3][5];
  TH2D * mul_compare_raw[3][5]; // [tbit] [sbit] -- (mul - raw) / mul vs. (mul/Nbx)
  char mul_compare_raw_n[3][5][256];
  char mul_compare_raw_t[3][5][512];
  TProfile * mul_compare_raw_pfx[3][5];
  // the following plots may not make much sense, since they compare Omega*rsc to raw and to mul
  TH2D * rsc_compare_raw[3][5]; // [tbit] [sbit] -- (Omega*rsc - raw) / (Omega*rsc) vs. (Omega*rsc/Nbx) 
  char rsc_compare_raw_n[3][5][256];
  char rsc_compare_raw_t[3][5][512];
  TProfile * rsc_compare_raw_pfx[3][5];
  TH2D * rsc_compare_mul[3][5]; // [tbit] [sbit] -- (Omega*rsc - mul) / (Omega*rsc) vs. (Omega*rsc/Nbx)
  char rsc_compare_mul_n[3][5][256];
  char rsc_compare_mul_t[3][5][512];
  TProfile * rsc_compare_mul_pfx[3][5];

  Double_t R_min[3][3][10];
  Double_t R_max[3][3][10];
  Double_t F_min[3][3][5];
  Double_t F_max[3][3][5];
  Double_t counts[3][3]; // [tbit] [cbit] 
  Double_t rate_array[3][3][var_bins_const]; // [tbit] [cbit] [run index]
  Double_t rate_array_max[3][3];
  Double_t zzz;
  Double_t zz1,zz2,zz3;
  Double_t maxx,maxx_tmp;
  Int_t ncbins=250;
  if(!strcmp(var,"i"))
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        // fill rate_array and compute max rates for each detector and cbit
        rate_array_max[t][c] = 0;
        for(Int_t b=1; b<var_bins; b++)
        {
          counts[t][c] = 0;
          for(Int_t s=0; s<4; s++) 
            counts[t][c] += mul_d[t][c][s]->GetBinContent(b);
          rate_array[t][c][b-1] = counts[t][c] / time_array[b-1];
          if(rate_array[t][c][b-1] > rate_array_max[t][c]) 
            rate_array_max[t][c] = rate_array[t][c][b-1];
        };

        // fill R3 rate dependence 2d hists
        for(Int_t r=1; r<10; r++)
        {
          sprintf(rate_dep_n[t][c][r],"rate_dep_R%d_%s%s",r,tbit[t],cbit[c]);
          sprintf(rate_dep_t[t][c][r],"R%d vs. %s%s corrected rate",r,tbit[t],cbit[c]);
          R_min[t][c][r] = R_mul_d[t][c][r]->GetMinimum();
          R_max[t][c][r] = R_mul_d[t][c][r]->GetMaximum();
          R_min[t][c][r] -= R_min[t][c][r] * 0.1;
          R_max[t][c][r] += R_max[t][c][r] * 0.1;
          rate_dep[t][c][r] = new TH2D(rate_dep_n[t][c][r],rate_dep_t[t][c][r],
              100,0,rate_array_max[t][c],
              100,R_min[t][c][r],R_max[t][c][r]);
          for(Int_t b=1; b<var_bins; b++)
          {
            zzz = R_mul_d[t][c][r]->GetBinContent(b);
            rate_dep[t][c][r]->Fill(rate_array[t][c][b-1],zzz);
          };
        };

        // fill correction factor (mul/raw) rate dependence plots
        for(Int_t s=0; s<4; s++)
        {
          sprintf(rate_fac_n[t][c][s],"rate_fac_%s%s_s%d",tbit[t],cbit[c],s);
          sprintf(rate_fac_t[t][c][s],"%s%s mul/raw vs. %s%s corrected rate -- %s",tbit[t],cbit[c],tbit[t],cbit[c],leg);
          F_min[t][c][s] = fac_d[t][c][s]->GetMinimum();
          F_max[t][c][s] = fac_d[t][c][s]->GetMaximum();
          F_min[t][c][s] -= F_min[t][c][s] * 0.1;
          F_max[t][c][s] += F_max[t][c][s] * 0.1;
          rate_fac[t][c][s] = new TH2D(rate_fac_n[t][c][s],rate_fac_t[t][c][s],
            10,0,rate_array_max[t][c],
            100,F_min[t][c][s],F_max[t][c][s]);
          for(Int_t b=1; b<var_bins; b++)
          {
            zzz = fac_d[t][c][s]->GetBinContent(b);
            rate_fac[t][c][s]->Fill(rate_array[t][c][b-1],zzz);
          };
          rate_fac_pfx[t][c][s] = rate_fac[t][c][s]->ProfileX();
        };

        // fill correction factor (Omega*rsc/rawx) rate dependence plots
        if(c==2)
        {
          for(Int_t s=0; s<4; s++)
          {
            sprintf(rate_rsr_n[t][s],"rate_rsr_%s_s%d",tbit[t],s);
            sprintf(rate_rsr_t[t][s],
              "%s #Omega*rsc/rawx vs. %sx acc+mul corrected rate -- %s",tbit[t],tbit[t],leg);
            rate_rsr[t][s] = new TH2D(rate_rsr_n[t][s],rate_rsr_t[t][s],
              10,0,rate_array_max[t][2],
              100,0,10);
            for(Int_t b=1; b<var_bins; b++)
            {
              zzz = rsr_d[t][s]->GetBinContent(b);
              rate_rsr[t][s]->Fill(rate_array[t][2][b-1],zzz);
            };
            rate_rsr_pfx[t][s] = rate_rsr[t][s]->ProfileX();
          };
        };


        // fill mul_compare_raw -- (mul - raw) / mul vs. (mul/Nbx)
        if(c==2)
        {
          for(Int_t s=0; s<4; s++)
          {
            sprintf(mul_compare_raw_n[t][s],"mul_compare_raw_%s_s%d",tbit[t],s);
            sprintf(mul_compare_raw_t[t][s],
              "(N_{%sx}^{mul} - N_{%sx}^{raw}) / N_{%sx}^{mul} vs. N_{%sx}^{mul}/N_{bx} -- %s",
              tbit[t],tbit[t],tbit[t],tbit[t],leg);

            maxx=0;
            for(Int_t b=1; b<var_bins; b++)
            {
              maxx_tmp = mul_d[t][c][s]->GetBinContent(b) / tot_d[s]->GetBinContent(b);
              maxx = (maxx_tmp > maxx) ? maxx_tmp : maxx;
            };

            mul_compare_raw[t][s] = new TH2D(mul_compare_raw_n[t][s],mul_compare_raw_t[t][s],
              ncbins,0,maxx,
              2*ncbins,-1,1);

            for(Int_t b=1; b<var_bins; b++)
            {
              zz1 = mul_d[t][c][s]->GetBinContent(b);
              zz2 = raw_d[t][c][s]->GetBinContent(b);
              zz3 = zz1 / tot_d[s]->GetBinContent(b);
              zzz = (zz1 - zz2) / zz1;
              mul_compare_raw[t][s]->Fill(zz3,zzz);
            };
            mul_compare_raw_pfx[t][s] = mul_compare_raw[t][s]->ProfileX();
          };
        };

        // fill rsc_compare_raw -- (Omega*rsc - raw) / (Omega*rsc) vs. (Omega*rsc/Nbx)
        if(c==2)
        {
          for(Int_t s=0; s<4; s++)
          {
            sprintf(rsc_compare_raw_n[t][s],"rsc_compare_raw_%s_s%d",tbit[t],s);
            sprintf(rsc_compare_raw_t[t][s],
              "(#Omega*N_{%s}^{rsc} - N_{%sx}^{raw})/(#Omega*N_{%s}^{rsc}) vs. #Omega*N_{%s}^{rsc}/N_{bx} -- %s",
              tbit[t],tbit[t],tbit[t],tbit[t],leg);

            maxx=0;
            for(Int_t b=1; b<var_bins; b++)
            {
              maxx_tmp = rsc_d[t][s]->GetBinContent(b) / tot_d[s]->GetBinContent(b);
              maxx = (maxx_tmp > maxx) ? maxx_tmp : maxx;
            };

            rsc_compare_raw[t][s] = new TH2D(rsc_compare_raw_n[t][s],rsc_compare_raw_t[t][s],
              ncbins,0,maxx,
              2*ncbins,-1,1);

            for(Int_t b=1; b<var_bins; b++)
            {
              zz1 = rsc_d[t][s]->GetBinContent(b);
              zz2 = raw_d[t][c][s]->GetBinContent(b);
              zz3 = zz1 / tot_d[s]->GetBinContent(b);
              zzz = (zz1 - zz2) / zz1;
              rsc_compare_raw[t][s]->Fill(zz3,zzz);
            };
            rsc_compare_raw_pfx[t][s] = rsc_compare_raw[t][s]->ProfileX();
          };
        };

        // fill rsc_compare_mul -- (Omega*rsc - mul) / (Omega*rsc) vs. (Omega*rsc/Nbx)
        if(c==2)
        {
          for(Int_t s=0; s<4; s++)
          {
            sprintf(rsc_compare_mul_n[t][s],"rsc_compare_mul_%s_s%d",tbit[t],s);
            sprintf(rsc_compare_mul_t[t][s],
              "(#Omega*N_{%s}^{rsc} - N_{%sx}^{mul})/(#Omega*N_{%s}^{rsc}) vs. #Omega*N_{%s}^{rsc}/N_{bx} -- %s",
              tbit[t],tbit[t],tbit[t],tbit[t],leg);

            maxx=0;
            for(Int_t b=1; b<var_bins; b++)
            {
              maxx_tmp = rsc_d[t][s]->GetBinContent(b) / tot_d[s]->GetBinContent(b);
              maxx = (maxx_tmp > maxx) ? maxx_tmp : maxx;
            };

            rsc_compare_mul[t][s] = new TH2D(rsc_compare_mul_n[t][s],rsc_compare_mul_t[t][s],
              ncbins,0,maxx,
              2*ncbins,-1,1);

            for(Int_t b=1; b<var_bins; b++)
            {
              zz1 = rsc_d[t][s]->GetBinContent(b);
              zz2 = mul_d[t][c][s]->GetBinContent(b);
              zz3 = zz1 / tot_d[s]->GetBinContent(b);
              zzz = (zz1 - zz2) / zz1;
              rsc_compare_mul[t][s]->Fill(zz3,zzz);
            };
            rsc_compare_mul_pfx[t][s] = rsc_compare_mul[t][s]->ProfileX();
          };
        };

      };
    }; // eo bit loops
  }; // eo if var=="i"


  // change tick marks (divisions) for bXing distributions
  if(!strcmp(var,"bx"))
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        for(Int_t s=0; s<5; s++)
        {
          raw_d[t][c][s]->SetNdivisions(524);
          acc_d[t][c][s]->SetNdivisions(524);
          mul_d[t][c][s]->SetNdivisions(524);
          fac_d[t][c][s]->SetNdivisions(524);
        };
      };
    };
  };

  Int_t sf; // drawing scale factor
  if(printPNGs) sf=3;
  else sf=1;



  // DRAW OUTPUT CANVASES
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetTitleFontSize(FSIZE);

  // controls drawing order whether we plot
  // fifth spinbit (sum over ++ -- +- -+) or not
  Int_t first_draw,first_same_draw;
  if(zeroAborts)
  {
    first_draw=0;
    first_same_draw=1;
  }
  else
  {
    first_draw=4;
    first_same_draw=0;
  };

  TCanvas * c_raw[3];
  char c_raw_n[3][32]; // [tbit] [char buffer]
  for(Int_t t=0; t<3; t++) 
  {
    sprintf(c_raw_n[t],"c_raw_%s",tbit[t]);
    c_raw[t] = new TCanvas(c_raw_n[t],c_raw_n[t],1100*sf,700*sf);
    c_raw[t]->Divide(1,3);
    for(Int_t ccc=1; ccc<=3; ccc++) 
    {
      c_raw[t]->GetPad(ccc)->SetGrid(1,1);
      if(drawLog) c_raw[t]->GetPad(ccc)->SetLogy();
      c_raw[t]->cd(ccc);
      raw_d[t][ccc-1][first_draw]->Draw();
      for(Int_t s=first_same_draw; s<4; s++) raw_d[t][ccc-1][s]->Draw("same"); 
    };
  };

  TCanvas * c_acc[3];
  char c_acc_n[3][32]; // [tbit] [char buffer]
  for(Int_t t=0; t<3; t++) 
  {
    sprintf(c_acc_n[t],"c_acc_%s",tbit[t]);
    c_acc[t] = new TCanvas(c_acc_n[t],c_acc_n[t],1100*sf,700*sf);
    c_acc[t]->Divide(1,3);
    for(Int_t ccc=1; ccc<=3; ccc++) 
    {
      if(drawLog) c_acc[t]->GetPad(ccc)->SetLogy();
      c_acc[t]->GetPad(ccc)->SetGrid(1,1);
      c_acc[t]->cd(ccc);
      acc_d[t][ccc-1][first_draw]->Draw();
      for(Int_t s=first_same_draw; s<4; s++) acc_d[t][ccc-1][s]->Draw("same");
    };
  };

  // scale all multiple corrected plots to be same
  if(specificRun>0)
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        for(Int_t s=0; s<5; s++)
        {
          if(mul_d[t][c][s]->GetMaximum() > 30e6)
            mul_d[t][c][s]->GetYaxis()->SetRangeUser(0,60e6);
          else
            mul_d[t][c][s]->GetYaxis()->SetRangeUser(0,30e6);
        };
      };
    };
  };


  TCanvas * c_mul[3];
  char c_mul_n[3][32]; // [tbit] [char buffer]
  for(Int_t t=0; t<3; t++) 
  {
    sprintf(c_mul_n[t],"c_mul_%s",tbit[t]);
    c_mul[t] = new TCanvas(c_mul_n[t],c_mul_n[t],1100*sf,700*sf);
    c_mul[t]->Divide(1,3);
    for(Int_t ccc=1; ccc<=3; ccc++) 
    {
      c_mul[t]->GetPad(ccc)->SetGrid(1,1);
      if(drawLog) c_mul[t]->GetPad(ccc)->SetLogy();
      c_mul[t]->cd(ccc);
      mul_d[t][ccc-1][first_draw]->Draw();
      for(Int_t s=first_same_draw; s<4; s++) mul_d[t][ccc-1][s]->Draw("same");
    };
  };
  


  TCanvas * c_fac[3];
  char c_fac_n[3][32]; // [tbit] [char buffer]
  for(Int_t t=0; t<3; t++) 
  {
    sprintf(c_fac_n[t],"c_fac_%s",tbit[t]);
    c_fac[t] = new TCanvas(c_fac_n[t],c_fac_n[t],1100*sf,700*sf);
    c_fac[t]->Divide(1,3);
    for(Int_t ccc=1; ccc<=3; ccc++) 
    {
      c_fac[t]->GetPad(ccc)->SetGrid(1,1);
      if(drawLog) c_fac[t]->GetPad(ccc)->SetLogy();
      c_fac[t]->cd(ccc);
      fac_d[t][ccc-1][first_draw]->Draw();
      for(Int_t s=first_same_draw; s<4; s++) fac_d[t][ccc-1][s]->Draw("same");
    };
  };


  TCanvas * c_rsc; // bbc,zdc,vpd all on same canvas
  char c_rsc_n[32];
  strcpy(c_rsc_n,"c_rsc");
  c_rsc = new TCanvas(c_rsc_n,c_rsc_n,1100*sf,700*sf);
  c_rsc->Divide(1,3);
  for(Int_t t=0; t<3; t++)
  {
    c_rsc->GetPad(t+1)->SetGrid(1,1);
    if(drawLog) c_rsc->GetPad(t+1)->SetLogy();
    c_rsc->cd(t+1);
    rsc_d[t][first_draw]->Draw();
    for(Int_t s=first_same_draw; s<4; s++) rsc_d[t][s]->Draw("same");
  };


  TCanvas * c_rsr; // bbc,zdc,vpd all on same canvas
  char c_rsr_n[32];
  strcpy(c_rsr_n,"c_rsr");
  c_rsr = new TCanvas(c_rsr_n,c_rsr_n,1100*sf,700*sf);
  c_rsr->Divide(1,3);
  for(Int_t t=0; t<3; t++)
  {
    c_rsr->GetPad(t+1)->SetGrid(1,1);
    if(drawLog) c_rsr->GetPad(t+1)->SetLogy();
    c_rsr->cd(t+1);
    rsr_d[t][first_draw]->Draw();
    for(Int_t s=first_same_draw; s<4; s++) rsr_d[t][s]->Draw("same");
  };

  
  TCanvas * c_R[3][10]; // [tbit] [rellum]
  char c_R_n[3][10][32]; // [tbit] [rellum] [char buffer]
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      sprintf(c_R_n[t][r],"c_R%d_%s",r,tbit[t]);
      c_R[t][r] = new TCanvas(c_R_n[t][r],c_R_n[t][r],1100*sf,940*sf);
      c_R[t][r]->Divide(1,4);
      for(Int_t ccc=1; ccc<=3; ccc++) 
      {
        c_R[t][r]->GetPad(ccc)->SetGrid(1,1);
        c_R[t][r]->cd(ccc);
        R_mul_d[t][ccc-1][r]->Draw();
      };
      c_R[t][r]->GetPad(4)->SetGrid(1,1);
      c_R[t][r]->cd(4);
      R_rsc_d[t][r]->Draw();
    };
  };


  TCanvas * c_D[10]; // [rellum]
  char c_D_n[10][32]; // [rellum] [char buffer]
  for(Int_t r=1; r<10; r++)
  {
    sprintf(c_D_n[r],"c_R%d_zdc_minus_vpd",r);
    c_D[r] = new TCanvas(c_D_n[r],c_D_n[r],1100*sf,940*sf);
    c_D[r]->Divide(1,4);
    for(Int_t ccc=1; ccc<=3; ccc++)
    {
      c_D[r]->GetPad(ccc)->SetGrid(1,1);
      c_D[r]->cd(ccc);
      D_mul_d[ccc-1][r]->Draw();
    };
    c_D[r]->GetPad(4)->SetGrid(1,1);
    c_D[r]->cd(4);
    D_rsc_d[r]->Draw();
  };

  TCanvas * c_RD[10]; // [rellum]
  char c_RD_n[10][32];
  for(Int_t r=1; r<10; r++)
  {
    sprintf(c_RD_n[r],"c_R%d_rsc_minus_mul",r);
    c_RD[r] = new TCanvas(c_RD_n[r],c_RD_n[r],1100*sf,700*sf);
    c_RD[r]->Divide(1,3);
    for(Int_t ttt=1; ttt<=3; ttt++)
    {
      c_RD[r]->GetPad(ttt)->SetGrid(1,1);
      c_RD[r]->cd(ttt);
      RD_d[ttt-1][r]->Draw();
    };
  };

  TCanvas * c_SD[3][10]; // [xbit] [rellum]
  char c_SD_n[3][10][32]; 
  for(Int_t r=1; r<10; r++)
  {
    sprintf(c_SD_n[0][r],"c_R%d_east_minus_west",r);
    sprintf(c_SD_n[1][r],"c_R%d_east_minus_coin",r);
    sprintf(c_SD_n[2][r],"c_R%d_west_minus_coin",r);
    for(Int_t x=0; x<3; x++)
    {
      c_SD[x][r] = new TCanvas(c_SD_n[x][r],c_SD_n[x][r],1100*sf,700*sf);
      c_SD[x][r]->Divide(1,3);
      for(Int_t ccc=1; ccc<=3; ccc++)
      {
        c_SD[x][r]->GetPad(ccc)->SetGrid(1,1);
        c_SD[x][r]->cd(ccc);
        SD_d[x][ccc-1][r]->Draw();
      };
    };
  };

  TCanvas * c_mean[10]; // [rellum]
  char c_mean_n[10][32];
  for(Int_t r=1; r<10; r++)
  {
    sprintf(c_mean_n[r],"c_mean_R%d",r);
    c_mean[r] = new TCanvas(c_mean_n[r],c_mean_n[r],1100*sf,700*sf);
    c_mean[r]->Divide(1,3);
    for(Int_t ccc=1; ccc<=3; ccc++)
    {
      c_mean[r]->GetPad(ccc)->SetGrid(1,1);
      c_mean[r]->cd(ccc);
      mean_R[ccc-1][r]->Draw();
    };
  };

  TCanvas * c_dev[3][10]; // [tbit] [rellum]
  char c_dev_n[3][10][32];
  for(Int_t r=1; r<10; r++)
  {
    for(Int_t t=0; t<3; t++)
    {
      sprintf(c_dev_n[t][r],"c_deviation_R%d_%s",r,tbit[t]);
      c_dev[t][r] = new TCanvas(c_dev_n[t][r],c_dev_n[t][r],1100*sf,700*sf);
      c_dev[t][r]->Divide(1,3);
      for(Int_t ccc=1; ccc<=3; ccc++)
      {
        c_dev[t][r]->GetPad(ccc)->SetGrid(1,1);
        c_dev[t][r]->cd(ccc);
        dev_R[t][ccc-1][r]->Draw();
      };
    };
  };


  TCanvas * c_spin_pat = new TCanvas("c_spin_pat","c_spin_pat",1100*sf,700*sf);
  spin_pat_rel->GetXaxis()->SetLabelSize(FSIZE);
  spin_pat_rel->GetYaxis()->SetLabelSize(FSIZE);
  spin_pat_rel->Draw();



  TCanvas * c_rate_fac[3]; // [tbit]
  char c_rate_fac_n[3][32];
  if(!strcmp(var,"i"))
  {
    for(Int_t t=0; t<3; t++)
    {
      sprintf(c_rate_fac_n[t],"c_rate_fac_%s",tbit[t]);
      c_rate_fac[t] = new TCanvas(c_rate_fac_n[t],c_rate_fac_n[t],800*sf,800*sf);
      c_rate_fac[t]->Divide(2,2);
      for(Int_t c=0; c<3; c++)
      {
        c_rate_fac[t]->cd(c+1);
        c_rate_fac[t]->GetPad(c+1)->SetGrid(1,1);
        rate_fac_pfx[t][c][0]->SetLineColor(kGreen+2);
        rate_fac_pfx[t][c][1]->SetLineColor(kOrange+7);
        rate_fac_pfx[t][c][2]->SetLineColor(kRed);
        rate_fac_pfx[t][c][3]->SetLineColor(kBlue);
        rate_fac_pfx[t][c][0]->Draw();
        for(Int_t s=1; s<4; s++) rate_fac_pfx[t][c][s]->Draw("same");
      };
    };
  };

  TCanvas * c_rate_rsr = new TCanvas("c_rate_rsr","c_rate_rsr",800*sf,800*sf);
  if(!strcmp(var,"i"))
  {
    c_rate_rsr->Divide(2,2);
    for(Int_t t=0; t<3; t++) 
    {
      c_rate_rsr->GetPad(t+1)->SetGrid(1,1);
      rate_rsr_pfx[t][0]->SetLineColor(kGreen+2);
      rate_rsr_pfx[t][1]->SetLineColor(kOrange+7);
      rate_rsr_pfx[t][2]->SetLineColor(kRed);
      rate_rsr_pfx[t][3]->SetLineColor(kBlue);
      c_rate_rsr->cd(t+1);
      rate_rsr_pfx[t][0]->Draw();
      for(Int_t s=1; s<4; s++) rate_rsr_pfx[t][s]->Draw("same");
    };
  };

  TCanvas * c_mul_compare_raw = new TCanvas("c_mul_compare_raw","c_mul_compare_raw",800*sf,800*sf);
  if(!strcmp(var,"i"))
  {
    c_mul_compare_raw->Divide(2,2);
    for(Int_t t=0; t<3; t++) 
    {
      c_mul_compare_raw->GetPad(t+1)->SetGrid(1,1);
      mul_compare_raw_pfx[t][0]->SetLineColor(kGreen+2);
      mul_compare_raw_pfx[t][1]->SetLineColor(kOrange+7);
      mul_compare_raw_pfx[t][2]->SetLineColor(kRed);
      mul_compare_raw_pfx[t][3]->SetLineColor(kBlue);
      c_mul_compare_raw->cd(t+1);
      mul_compare_raw_pfx[t][0]->Draw();
      for(Int_t s=1; s<4; s++) mul_compare_raw_pfx[t][s]->Draw("same");
    };
  };
  
  TCanvas * c_rsc_compare_raw = new TCanvas("c_rsc_compare_raw","c_rsc_compare_raw",800*sf,800*sf);
  if(!strcmp(var,"i"))
  {
    c_rsc_compare_raw->Divide(2,2);
    for(Int_t t=0; t<3; t++) 
    {
      c_rsc_compare_raw->GetPad(t+1)->SetGrid(1,1);
      rsc_compare_raw_pfx[t][0]->SetLineColor(kGreen+2);
      rsc_compare_raw_pfx[t][1]->SetLineColor(kOrange+7);
      rsc_compare_raw_pfx[t][2]->SetLineColor(kRed);
      rsc_compare_raw_pfx[t][3]->SetLineColor(kBlue);
      c_rsc_compare_raw->cd(t+1);
      rsc_compare_raw_pfx[t][0]->Draw();
      for(Int_t s=1; s<4; s++) rsc_compare_raw_pfx[t][s]->Draw("same");
    };
  };

  TCanvas * c_rsc_compare_mul = new TCanvas("c_rsc_compare_mul","c_rsc_compare_mul",800*sf,800*sf);
  if(!strcmp(var,"i"))
  {
    c_rsc_compare_mul->Divide(2,2);
    for(Int_t t=0; t<3; t++) 
    {
      c_rsc_compare_mul->GetPad(t+1)->SetGrid(1,1);
      rsc_compare_mul_pfx[t][0]->SetLineColor(kGreen+2);
      rsc_compare_mul_pfx[t][1]->SetLineColor(kOrange+7);
      rsc_compare_mul_pfx[t][2]->SetLineColor(kRed);
      rsc_compare_mul_pfx[t][3]->SetLineColor(kBlue);
      c_rsc_compare_mul->cd(t+1);
      rsc_compare_mul_pfx[t][0]->Draw();
      for(Int_t s=1; s<4; s++) rsc_compare_mul_pfx[t][s]->Draw("same");
    };
  };



  /*
  // draw only ZDC and VPD for these plots (for FMS meeting 11.10.14) 
  //
  TCanvas * c_mul_compare_raw = new TCanvas("c_mul_compare_raw","c_mul_compare_raw",800*sf,400*sf);
  if(!strcmp(var,"i"))
  {
    c_mul_compare_raw->Divide(2,1);
    for(Int_t t=1; t<3; t++) 
    {
      c_mul_compare_raw->GetPad(t)->SetGrid(1,1);
      mul_compare_raw_pfx[t][0]->SetLineColor(kGreen+2);
      mul_compare_raw_pfx[t][1]->SetLineColor(kOrange+7);
      mul_compare_raw_pfx[t][2]->SetLineColor(kRed);
      mul_compare_raw_pfx[t][3]->SetLineColor(kBlue);
      c_mul_compare_raw->cd(t);
      mul_compare_raw_pfx[t][0]->Draw();
      for(Int_t s=1; s<4; s++) mul_compare_raw_pfx[t][s]->Draw("same");
    };
  };
  
  TCanvas * c_rsc_compare_raw = new TCanvas("c_rsc_compare_raw","c_rsc_compare_raw",800*sf,400*sf);
  if(!strcmp(var,"i"))
  {
    c_rsc_compare_raw->Divide(2,1);
    for(Int_t t=1; t<3; t++) 
    {
      c_rsc_compare_raw->GetPad(t)->SetGrid(1,1);
      rsc_compare_raw_pfx[t][0]->SetLineColor(kGreen+2);
      rsc_compare_raw_pfx[t][1]->SetLineColor(kOrange+7);
      rsc_compare_raw_pfx[t][2]->SetLineColor(kRed);
      rsc_compare_raw_pfx[t][3]->SetLineColor(kBlue);
      c_rsc_compare_raw->cd(t);
      rsc_compare_raw_pfx[t][0]->Draw();
      for(Int_t s=1; s<4; s++) rsc_compare_raw_pfx[t][s]->Draw("same");
    };
  };
  */





  // write 
  if(specificFill==0 && specificRun==0)
  {
    c_spin_pat->Write();
    for(Int_t t=0; t<3; t++) c_raw[t]->Write();
    for(Int_t t=0; t<3; t++) c_acc[t]->Write();
    for(Int_t t=0; t<3; t++) c_mul[t]->Write();
    for(Int_t t=0; t<3; t++) c_fac[t]->Write();
    c_rsc->Write();
    c_rsr->Write();
    for(Int_t r=1; r<10; r++)
    {
      for(Int_t t=0; t<3; t++)
      {
        c_R[t][r]->Write();
      };
    };
    for(Int_t r=1; r<10; r++) c_mean[r]->Write();
    for(Int_t r=1; r<10; r++) c_D[r]->Write();
    for(Int_t r=1; r<10; r++) 
    {
      for(Int_t x=0; x<3; x++)
      {
        c_SD[x][r]->Write();
      };
    };
    for(Int_t r=1; r<10; r++) c_RD[r]->Write();
    for(Int_t r=1; r<10; r++)
    {
      for(Int_t t=0; t<3; t++)
      {
        c_dev[t][r]->Write();
      };
    };
    if(!strcmp(var,"i"))
    {
      for(Int_t r=1; r<10; r++)
      {
        for(Int_t t=0; t<3; t++)
        {
          for(Int_t c=0; c<3; c++)
          {
            rate_dep[t][c][r]->Write();
          };
        };
      };
      /*
      for(Int_t t=0; t<3; t++)
      {
        for(Int_t c=0; c<3; c++)
        {
          for(Int_t s=0; s<4; s++)
          {
            rate_fac[t][c][s]->Write();
            rate_fac_pfx[t][c][s]->Write();
          };
        };
      };
      */
      if(!strcmp(var,"i"))
      {
        for(Int_t t=0; t<3; t++)
        {
          c_rate_fac[t]->Write();
        };
      };
      c_rate_rsr->Write();
      c_mul_compare_raw->Write();
      c_rsc_compare_raw->Write();
      c_rsc_compare_mul->Write();
    };
  };


  
  // produce output tree (iff var=="i")
  TTree * rtr = new TTree("rtr","rtr");
  Float_t RR[3][10]; // [tbit] [rellum]
  Float_t RR_err[3][10]; // [tbit] [rellum]
  Float_t RR_rscerr[3][10]; // [tbit] [rellum]
  Float_t RRe[3][10]; // [tbit] [rellum]
  Float_t RRw[3][10]; // [tbit] [rellum]
  Float_t RRx[3][10]; // [tbit] [rellum]
  Float_t RRrsc[3][10]; // [tbit] [rellum]
  Float_t d_vz[3]; // [cbit] --- diagnostic (only for R3!)
  Float_t d_xx[3][3]; // [xbit] [tbit]
  Float_t R_LL[3]; // [cbit]
  if(!strcmp(var,"i") && specificFill==0 && specificRun==0)
  {
    rtr->Branch("i",&index,"i/I");
    rtr->Branch("runnum",&runnum,"runnum/I");
    rtr->Branch("fill",&fill,"fill/I");
    rtr->Branch("t",&time,"t/D");

    // n.b. it's probably better to write arrays to branches, but that's not
    //      how the original downstream code was designed to interpret this tree
    rtr->Branch("R1_bbce",&(RRe[0][1]),"R1_bbce/F"); // bbce relative luminosity (using mul)
    rtr->Branch("R2_bbce",&(RRe[0][2]),"R2_bbce/F");
    rtr->Branch("R3_bbce",&(RRe[0][3]),"R3_bbce/F");
    rtr->Branch("R4_bbce",&(RRe[0][4]),"R4_bbce/F");
    rtr->Branch("R5_bbce",&(RRe[0][5]),"R5_bbce/F");
    rtr->Branch("R6_bbce",&(RRe[0][6]),"R6_bbce/F");
    rtr->Branch("R7_bbce",&(RRe[0][7]),"R7_bbce/F");
    rtr->Branch("R8_bbce",&(RRe[0][8]),"R8_bbce/F");
    rtr->Branch("R9_bbce",&(RRe[0][9]),"R9_bbce/F");
    rtr->Branch("R1_zdce",&(RRe[1][1]),"R1_zdce/F"); // zdce relative luminosity (using mul)
    rtr->Branch("R2_zdce",&(RRe[1][2]),"R2_zdce/F");
    rtr->Branch("R3_zdce",&(RRe[1][3]),"R3_zdce/F");
    rtr->Branch("R4_zdce",&(RRe[1][4]),"R4_zdce/F");
    rtr->Branch("R5_zdce",&(RRe[1][5]),"R5_zdce/F");
    rtr->Branch("R6_zdce",&(RRe[1][6]),"R6_zdce/F");
    rtr->Branch("R7_zdce",&(RRe[1][7]),"R7_zdce/F");
    rtr->Branch("R8_zdce",&(RRe[1][8]),"R8_zdce/F");
    rtr->Branch("R9_zdce",&(RRe[1][9]),"R9_zdce/F");
    rtr->Branch("R1_vpde",&(RRe[2][1]),"R1_vpde/F"); // vpde relative luminosity (using mul)
    rtr->Branch("R2_vpde",&(RRe[2][2]),"R2_vpde/F");
    rtr->Branch("R3_vpde",&(RRe[2][3]),"R3_vpde/F");
    rtr->Branch("R4_vpde",&(RRe[2][4]),"R4_vpde/F");
    rtr->Branch("R5_vpde",&(RRe[2][5]),"R5_vpde/F");
    rtr->Branch("R6_vpde",&(RRe[2][6]),"R6_vpde/F");
    rtr->Branch("R7_vpde",&(RRe[2][7]),"R7_vpde/F");
    rtr->Branch("R8_vpde",&(RRe[2][8]),"R8_vpde/F");
    rtr->Branch("R9_vpde",&(RRe[2][9]),"R9_vpde/F");

    rtr->Branch("R1_bbcw",&(RRw[0][1]),"R1_bbcw/F"); // bbcw relative luminosity (using mul)
    rtr->Branch("R2_bbcw",&(RRw[0][2]),"R2_bbcw/F");
    rtr->Branch("R3_bbcw",&(RRw[0][3]),"R3_bbcw/F");
    rtr->Branch("R4_bbcw",&(RRw[0][4]),"R4_bbcw/F");
    rtr->Branch("R5_bbcw",&(RRw[0][5]),"R5_bbcw/F");
    rtr->Branch("R6_bbcw",&(RRw[0][6]),"R6_bbcw/F");
    rtr->Branch("R7_bbcw",&(RRw[0][7]),"R7_bbcw/F");
    rtr->Branch("R8_bbcw",&(RRw[0][8]),"R8_bbcw/F");
    rtr->Branch("R9_bbcw",&(RRw[0][9]),"R9_bbcw/F");
    rtr->Branch("R1_zdcw",&(RRw[1][1]),"R1_zdcw/F"); // zdcw relative luminosity (using mul)
    rtr->Branch("R2_zdcw",&(RRw[1][2]),"R2_zdcw/F");
    rtr->Branch("R3_zdcw",&(RRw[1][3]),"R3_zdcw/F");
    rtr->Branch("R4_zdcw",&(RRw[1][4]),"R4_zdcw/F");
    rtr->Branch("R5_zdcw",&(RRw[1][5]),"R5_zdcw/F");
    rtr->Branch("R6_zdcw",&(RRw[1][6]),"R6_zdcw/F");
    rtr->Branch("R7_zdcw",&(RRw[1][7]),"R7_zdcw/F");
    rtr->Branch("R8_zdcw",&(RRw[1][8]),"R8_zdcw/F");
    rtr->Branch("R9_zdcw",&(RRw[1][9]),"R9_zdcw/F");
    rtr->Branch("R1_vpdw",&(RRw[2][1]),"R1_vpdw/F"); // vpdw relative luminosity (using mul)
    rtr->Branch("R2_vpdw",&(RRw[2][2]),"R2_vpdw/F");
    rtr->Branch("R3_vpdw",&(RRw[2][3]),"R3_vpdw/F");
    rtr->Branch("R4_vpdw",&(RRw[2][4]),"R4_vpdw/F");
    rtr->Branch("R5_vpdw",&(RRw[2][5]),"R5_vpdw/F");
    rtr->Branch("R6_vpdw",&(RRw[2][6]),"R6_vpdw/F");
    rtr->Branch("R7_vpdw",&(RRw[2][7]),"R7_vpdw/F");
    rtr->Branch("R8_vpdw",&(RRw[2][8]),"R8_vpdw/F");
    rtr->Branch("R9_vpdw",&(RRw[2][9]),"R9_vpdw/F");

    rtr->Branch("R1_bbcx",&(RRx[0][1]),"R1_bbcx/F"); // bbcx relative luminosity (using mul)
    rtr->Branch("R2_bbcx",&(RRx[0][2]),"R2_bbcx/F");
    rtr->Branch("R3_bbcx",&(RRx[0][3]),"R3_bbcx/F");
    rtr->Branch("R4_bbcx",&(RRx[0][4]),"R4_bbcx/F");
    rtr->Branch("R5_bbcx",&(RRx[0][5]),"R5_bbcx/F");
    rtr->Branch("R6_bbcx",&(RRx[0][6]),"R6_bbcx/F");
    rtr->Branch("R7_bbcx",&(RRx[0][7]),"R7_bbcx/F");
    rtr->Branch("R8_bbcx",&(RRx[0][8]),"R8_bbcx/F");
    rtr->Branch("R9_bbcx",&(RRx[0][9]),"R9_bbcx/F");
    rtr->Branch("R1_zdcx",&(RRx[1][1]),"R1_zdcx/F"); // zdcx relative luminosity (using mul)
    rtr->Branch("R2_zdcx",&(RRx[1][2]),"R2_zdcx/F");
    rtr->Branch("R3_zdcx",&(RRx[1][3]),"R3_zdcx/F");
    rtr->Branch("R4_zdcx",&(RRx[1][4]),"R4_zdcx/F");
    rtr->Branch("R5_zdcx",&(RRx[1][5]),"R5_zdcx/F");
    rtr->Branch("R6_zdcx",&(RRx[1][6]),"R6_zdcx/F");
    rtr->Branch("R7_zdcx",&(RRx[1][7]),"R7_zdcx/F");
    rtr->Branch("R8_zdcx",&(RRx[1][8]),"R8_zdcx/F");
    rtr->Branch("R9_zdcx",&(RRx[1][9]),"R9_zdcx/F");
    rtr->Branch("R1_vpdx",&(RRx[2][1]),"R1_vpdx/F"); // vpdx relative luminosity (using mul)
    rtr->Branch("R2_vpdx",&(RRx[2][2]),"R2_vpdx/F");
    rtr->Branch("R3_vpdx",&(RRx[2][3]),"R3_vpdx/F");
    rtr->Branch("R4_vpdx",&(RRx[2][4]),"R4_vpdx/F");
    rtr->Branch("R5_vpdx",&(RRx[2][5]),"R5_vpdx/F");
    rtr->Branch("R6_vpdx",&(RRx[2][6]),"R6_vpdx/F");
    rtr->Branch("R7_vpdx",&(RRx[2][7]),"R7_vpdx/F");
    rtr->Branch("R8_vpdx",&(RRx[2][8]),"R8_vpdx/F");
    rtr->Branch("R9_vpdx",&(RRx[2][9]),"R9_vpdx/F");

    rtr->Branch("R1_bbcrsc",&(RRrsc[0][1]),"R1_bbcrsc/F"); // bbc_rsc relative luminosity (using rsc)
    rtr->Branch("R2_bbcrsc",&(RRrsc[0][2]),"R2_bbcrsc/F");
    rtr->Branch("R3_bbcrsc",&(RRrsc[0][3]),"R3_bbcrsc/F");
    rtr->Branch("R4_bbcrsc",&(RRrsc[0][4]),"R4_bbcrsc/F");
    rtr->Branch("R5_bbcrsc",&(RRrsc[0][5]),"R5_bbcrsc/F");
    rtr->Branch("R6_bbcrsc",&(RRrsc[0][6]),"R6_bbcrsc/F");
    rtr->Branch("R7_bbcrsc",&(RRrsc[0][7]),"R7_bbcrsc/F");
    rtr->Branch("R8_bbcrsc",&(RRrsc[0][8]),"R8_bbcrsc/F");
    rtr->Branch("R9_bbcrsc",&(RRrsc[0][9]),"R9_bbcrsc/F");
    rtr->Branch("R1_zdcrsc",&(RRrsc[1][1]),"R1_zdcrsc/F"); // zdc_rsc relative luminosity (using rsc)
    rtr->Branch("R2_zdcrsc",&(RRrsc[1][2]),"R2_zdcrsc/F");
    rtr->Branch("R3_zdcrsc",&(RRrsc[1][3]),"R3_zdcrsc/F");
    rtr->Branch("R4_zdcrsc",&(RRrsc[1][4]),"R4_zdcrsc/F");
    rtr->Branch("R5_zdcrsc",&(RRrsc[1][5]),"R5_zdcrsc/F");
    rtr->Branch("R6_zdcrsc",&(RRrsc[1][6]),"R6_zdcrsc/F");
    rtr->Branch("R7_zdcrsc",&(RRrsc[1][7]),"R7_zdcrsc/F");
    rtr->Branch("R8_zdcrsc",&(RRrsc[1][8]),"R8_zdcrsc/F");
    rtr->Branch("R9_zdcrsc",&(RRrsc[1][9]),"R9_zdcrsc/F");
    rtr->Branch("R1_vpdrsc",&(RRrsc[2][1]),"R1_vpdrsc/F"); // vpd_rsc relative luminosity (using rsc)
    rtr->Branch("R2_vpdrsc",&(RRrsc[2][2]),"R2_vpdrsc/F");
    rtr->Branch("R3_vpdrsc",&(RRrsc[2][3]),"R3_vpdrsc/F");
    rtr->Branch("R4_vpdrsc",&(RRrsc[2][4]),"R4_vpdrsc/F");
    rtr->Branch("R5_vpdrsc",&(RRrsc[2][5]),"R5_vpdrsc/F");
    rtr->Branch("R6_vpdrsc",&(RRrsc[2][6]),"R6_vpdrsc/F");
    rtr->Branch("R7_vpdrsc",&(RRrsc[2][7]),"R7_vpdrsc/F");
    rtr->Branch("R8_vpdrsc",&(RRrsc[2][8]),"R8_vpdrsc/F");
    rtr->Branch("R9_vpdrsc",&(RRrsc[2][9]),"R9_vpdrsc/F");

    // -- mean rellum branches (means over e,w,x)
    rtr->Branch("R1_bbc_mean",&(RR[0][1]),"R1_bbc_mean/F"); // mean bbc relative luminosity
    rtr->Branch("R2_bbc_mean",&(RR[0][2]),"R2_bbc_mean/F");
    rtr->Branch("R3_bbc_mean",&(RR[0][3]),"R3_bbc_mean/F");
    rtr->Branch("R4_bbc_mean",&(RR[0][4]),"R4_bbc_mean/F");
    rtr->Branch("R5_bbc_mean",&(RR[0][5]),"R5_bbc_mean/F");
    rtr->Branch("R6_bbc_mean",&(RR[0][6]),"R6_bbc_mean/F");
    rtr->Branch("R7_bbc_mean",&(RR[0][7]),"R7_bbc_mean/F");
    rtr->Branch("R8_bbc_mean",&(RR[0][8]),"R8_bbc_mean/F");
    rtr->Branch("R9_bbc_mean",&(RR[0][9]),"R9_bbc_mean/F");
    rtr->Branch("R1_zdc_mean",&(RR[1][1]),"R1_zdc_mean/F"); // mean zdc relative luminosity
    rtr->Branch("R2_zdc_mean",&(RR[1][2]),"R2_zdc_mean/F");
    rtr->Branch("R3_zdc_mean",&(RR[1][3]),"R3_zdc_mean/F");
    rtr->Branch("R4_zdc_mean",&(RR[1][4]),"R4_zdc_mean/F");
    rtr->Branch("R5_zdc_mean",&(RR[1][5]),"R5_zdc_mean/F");
    rtr->Branch("R6_zdc_mean",&(RR[1][6]),"R6_zdc_mean/F");
    rtr->Branch("R7_zdc_mean",&(RR[1][7]),"R7_zdc_mean/F");
    rtr->Branch("R8_zdc_mean",&(RR[1][8]),"R8_zdc_mean/F");
    rtr->Branch("R9_zdc_mean",&(RR[1][9]),"R9_zdc_mean/F");
    rtr->Branch("R1_vpd_mean",&(RR[2][1]),"R1_vpd_mean/F"); // mean vpd relative luminosity
    rtr->Branch("R2_vpd_mean",&(RR[2][2]),"R2_vpd_mean/F");
    rtr->Branch("R3_vpd_mean",&(RR[2][3]),"R3_vpd_mean/F");
    rtr->Branch("R4_vpd_mean",&(RR[2][4]),"R4_vpd_mean/F");
    rtr->Branch("R5_vpd_mean",&(RR[2][5]),"R5_vpd_mean/F");
    rtr->Branch("R6_vpd_mean",&(RR[2][6]),"R6_vpd_mean/F");
    rtr->Branch("R7_vpd_mean",&(RR[2][7]),"R7_vpd_mean/F");
    rtr->Branch("R8_vpd_mean",&(RR[2][8]),"R8_vpd_mean/F");
    rtr->Branch("R9_vpd_mean",&(RR[2][9]),"R9_vpd_mean/F");


    // mean rellum errors
    rtr->Branch("R1_bbc_mean_err",&(RR_err[0][1]),"R1_bbc_mean_err/F"); // bbc_{e,w,x mean} rellum error
    rtr->Branch("R2_bbc_mean_err",&(RR_err[0][2]),"R2_bbc_mean_err/F");
    rtr->Branch("R3_bbc_mean_err",&(RR_err[0][3]),"R3_bbc_mean_err/F");
    rtr->Branch("R4_bbc_mean_err",&(RR_err[0][4]),"R4_bbc_mean_err/F");
    rtr->Branch("R5_bbc_mean_err",&(RR_err[0][5]),"R5_bbc_mean_err/F");
    rtr->Branch("R6_bbc_mean_err",&(RR_err[0][6]),"R6_bbc_mean_err/F");
    rtr->Branch("R7_bbc_mean_err",&(RR_err[0][7]),"R7_bbc_mean_err/F");
    rtr->Branch("R8_bbc_mean_err",&(RR_err[0][8]),"R8_bbc_mean_err/F");
    rtr->Branch("R9_bbc_mean_err",&(RR_err[0][9]),"R9_bbc_mean_err/F");
    rtr->Branch("R1_zdc_mean_err",&(RR_err[1][1]),"R1_zdc_mean_err/F"); // zdc_{e,w,x mean} rellum error
    rtr->Branch("R2_zdc_mean_err",&(RR_err[1][2]),"R2_zdc_mean_err/F");
    rtr->Branch("R3_zdc_mean_err",&(RR_err[1][3]),"R3_zdc_mean_err/F");
    rtr->Branch("R4_zdc_mean_err",&(RR_err[1][4]),"R4_zdc_mean_err/F");
    rtr->Branch("R5_zdc_mean_err",&(RR_err[1][5]),"R5_zdc_mean_err/F");
    rtr->Branch("R6_zdc_mean_err",&(RR_err[1][6]),"R6_zdc_mean_err/F");
    rtr->Branch("R7_zdc_mean_err",&(RR_err[1][7]),"R7_zdc_mean_err/F");
    rtr->Branch("R8_zdc_mean_err",&(RR_err[1][8]),"R8_zdc_mean_err/F");
    rtr->Branch("R9_zdc_mean_err",&(RR_err[1][9]),"R9_zdc_mean_err/F");
    rtr->Branch("R1_vpd_mean_err",&(RR_err[2][1]),"R1_vpd_mean_err/F"); // vpd_{e,w,x mean} rellum error
    rtr->Branch("R2_vpd_mean_err",&(RR_err[2][2]),"R2_vpd_mean_err/F");
    rtr->Branch("R3_vpd_mean_err",&(RR_err[2][3]),"R3_vpd_mean_err/F");
    rtr->Branch("R4_vpd_mean_err",&(RR_err[2][4]),"R4_vpd_mean_err/F");
    rtr->Branch("R5_vpd_mean_err",&(RR_err[2][5]),"R5_vpd_mean_err/F");
    rtr->Branch("R6_vpd_mean_err",&(RR_err[2][6]),"R6_vpd_mean_err/F");
    rtr->Branch("R7_vpd_mean_err",&(RR_err[2][7]),"R7_vpd_mean_err/F");
    rtr->Branch("R8_vpd_mean_err",&(RR_err[2][8]),"R8_vpd_mean_err/F");
    rtr->Branch("R9_vpd_mean_err",&(RR_err[2][9]),"R9_vpd_mean_err/F");

    // rsc rellum error bars
    rtr->Branch("R1_bbc_rsc_err",&(RR_rscerr[0][1]),"R1_bbc_rsc_err/F"); // bbc_rsc rellum error
    rtr->Branch("R2_bbc_rsc_err",&(RR_rscerr[0][2]),"R2_bbc_rsc_err/F");
    rtr->Branch("R3_bbc_rsc_err",&(RR_rscerr[0][3]),"R3_bbc_rsc_err/F");
    rtr->Branch("R4_bbc_rsc_err",&(RR_rscerr[0][4]),"R4_bbc_rsc_err/F");
    rtr->Branch("R5_bbc_rsc_err",&(RR_rscerr[0][5]),"R5_bbc_rsc_err/F");
    rtr->Branch("R6_bbc_rsc_err",&(RR_rscerr[0][6]),"R6_bbc_rsc_err/F");
    rtr->Branch("R7_bbc_rsc_err",&(RR_rscerr[0][7]),"R7_bbc_rsc_err/F");
    rtr->Branch("R8_bbc_rsc_err",&(RR_rscerr[0][8]),"R8_bbc_rsc_err/F");
    rtr->Branch("R9_bbc_rsc_err",&(RR_rscerr[0][9]),"R9_bbc_rsc_err/F");
    rtr->Branch("R1_zdc_rsc_err",&(RR_rscerr[1][1]),"R1_zdc_rsc_err/F"); // zdc_rsc rellum error
    rtr->Branch("R2_zdc_rsc_err",&(RR_rscerr[1][2]),"R2_zdc_rsc_err/F");
    rtr->Branch("R3_zdc_rsc_err",&(RR_rscerr[1][3]),"R3_zdc_rsc_err/F");
    rtr->Branch("R4_zdc_rsc_err",&(RR_rscerr[1][4]),"R4_zdc_rsc_err/F");
    rtr->Branch("R5_zdc_rsc_err",&(RR_rscerr[1][5]),"R5_zdc_rsc_err/F");
    rtr->Branch("R6_zdc_rsc_err",&(RR_rscerr[1][6]),"R6_zdc_rsc_err/F");
    rtr->Branch("R7_zdc_rsc_err",&(RR_rscerr[1][7]),"R7_zdc_rsc_err/F");
    rtr->Branch("R8_zdc_rsc_err",&(RR_rscerr[1][8]),"R8_zdc_rsc_err/F");
    rtr->Branch("R9_zdc_rsc_err",&(RR_rscerr[1][9]),"R9_zdc_rsc_err/F");
    rtr->Branch("R1_vpd_rsc_err",&(RR_rscerr[2][1]),"R1_vpd_rsc_err/F"); // vpd_rsc rellum error
    rtr->Branch("R2_vpd_rsc_err",&(RR_rscerr[2][2]),"R2_vpd_rsc_err/F");
    rtr->Branch("R3_vpd_rsc_err",&(RR_rscerr[2][3]),"R3_vpd_rsc_err/F");
    rtr->Branch("R4_vpd_rsc_err",&(RR_rscerr[2][4]),"R4_vpd_rsc_err/F");
    rtr->Branch("R5_vpd_rsc_err",&(RR_rscerr[2][5]),"R5_vpd_rsc_err/F");
    rtr->Branch("R6_vpd_rsc_err",&(RR_rscerr[2][6]),"R6_vpd_rsc_err/F");
    rtr->Branch("R7_vpd_rsc_err",&(RR_rscerr[2][7]),"R7_vpd_rsc_err/F");
    rtr->Branch("R8_vpd_rsc_err",&(RR_rscerr[2][8]),"R8_vpd_rsc_err/F");
    rtr->Branch("R9_vpd_rsc_err",&(RR_rscerr[2][9]),"R9_vpd_rsc_err/F");


    // -- diagnostic branches (only for R3)
    rtr->Branch("d_vz_e",&(d_vz[0]),"d_vz_e/F"); // ZDCE - VPDE  ( [cbit] ) (using multiples corrections)
    rtr->Branch("d_vz_w",&(d_vz[1]),"d_vz_w/F"); // ZDCW - VPDW
    rtr->Branch("d_vz_x",&(d_vz[2]),"d_vz_x/F"); // ZDCX - VPDX

    rtr->Branch("d_ew_bbc",&(d_xx[0][0]),"d_ew_bbc/F"); // BBCE - BBCW  ( [xbit] [tbit] )
    rtr->Branch("d_ew_zdc",&(d_xx[0][1]),"d_ew_zdc/F"); // ZDCE - ZDCW
    rtr->Branch("d_ew_vpd",&(d_xx[0][2]),"d_ew_vpd/F"); // VPDE - VPDW

    rtr->Branch("d_ex_bbc",&(d_xx[1][0]),"d_ex_bbc/F"); // BBCE - BBCX 
    rtr->Branch("d_ex_zdc",&(d_xx[1][1]),"d_ex_zdc/F"); // ZDCE - ZDCX
    rtr->Branch("d_ex_vpd",&(d_xx[1][2]),"d_ex_vpd/F"); // VPDE - VPDX

    rtr->Branch("d_wx_bbc",&(d_xx[2][0]),"d_wx_bbc/F"); // BBCW - BBCX 
    rtr->Branch("d_wx_zdc",&(d_xx[2][1]),"d_wx_zdc/F"); // ZDCW - ZDCX
    rtr->Branch("d_wx_vpd",&(d_xx[2][2]),"d_wx_vpd/F"); // VPDW - VPDX


    
    for(Int_t b=1; b<=var_bins; b++)
    {
      index = b;
      runnum = runnum_array[b-1];
      fill = fill_array[b-1];
      time = time_array[b-1];
      //printf("%d\n",time);
      for(Int_t t=0; t<3; t++)
      {
        for(Int_t r=1; r<10; r++)
        {
          RR[t][r] = mean_R[t][r]->GetBinContent(b);
          RRe[t][r] = R_mul_d[t][0][r]->GetBinContent(b);
          RRw[t][r] = R_mul_d[t][1][r]->GetBinContent(b);
          RRx[t][r] = R_mul_d[t][2][r]->GetBinContent(b);
          RRrsc[t][r] = R_rsc_d[t][r]->GetBinContent(b);
          RR_err[t][r] = mean_R[t][r]->GetBinError(b);
          RR_rscerr[t][r] = R_rsc_d[t][r]->GetBinError(b);
        };
      };
      for(Int_t c=0; c<3; c++)
      {
        d_vz[c] = D_mul_d[c][3]->GetBinContent(b);
      };
      for(Int_t t=0; t<3; t++)
      {
        for(Int_t x=0; x<3; x++)
        {
          d_xx[x][t] = SD_d[x][t][3]->GetBinContent(b);
        };
      };


      rtr->Fill();
    };
    rtr->Write();
  };

  // print pngs
  if(printPNGs)
  {
    char pngdir[32]; 
    if(specificFill==0 && specificRun==0) sprintf(pngdir,"png_rellum");
    else if(specificFill>0) sprintf(pngdir,"pdf_bXings_fills/%d",specificFill);
    else if(specificRun>0) sprintf(pngdir,"pdf_bXings_runs/%d",specificRun);
    char mkdir[64];
    sprintf(mkdir,".! mkdir -p %s",pngdir);
    gROOT->ProcessLine(mkdir);
    char c_raw_png[3][256];
    char c_acc_png[3][256];
    char c_mul_png[3][256];
    char c_fac_png[3][256];
    char c_rsc_png[256];
    char c_rsr_png[256];
    char c_R_png[3][10][256];
    char c_dev_png[3][10][256];
    char c_mean_png[10][256];
    char c_RD_png[10][256];
    char c_D_png[10][256];
    char c_SD_png[3][10][256];
    char c_rate_fac_png[3][256];
    char c_rate_rsr_png[256];
    char c_mul_compare_raw_png[256];
    char c_rsc_compare_raw_png[256];
    char c_rsc_compare_mul_png[256];
    char file_type[4];
    if(specificFill>0 || specificRun>0) sprintf(file_type,"pdf"); // make pdfs instead of pngs for specific fills
    else sprintf(file_type,"png");
    for(Int_t t=0; t<3; t++) 
    {
      sprintf(c_raw_png[t],"%s/raw_%s_%s.%s",pngdir,tbit[t],var,file_type);
      sprintf(c_acc_png[t],"%s/acc_%s_%s.%s",pngdir,tbit[t],var,file_type);
      sprintf(c_mul_png[t],"%s/mul_%s_%s.%s",pngdir,tbit[t],var,file_type);
      sprintf(c_fac_png[t],"%s/fac_%s_%s.%s",pngdir,tbit[t],var,file_type);
      sprintf(c_rate_fac_png[t],"%s/rate_fac_%s.%s",pngdir,tbit[t],file_type);
      c_raw[t]->Print(c_raw_png[t],file_type);
      c_acc[t]->Print(c_acc_png[t],file_type);
      c_mul[t]->Print(c_mul_png[t],file_type);
      c_fac[t]->Print(c_fac_png[t],file_type);
      if(specificFill==0 && specificRun==0)
      {
        for(Int_t r=1; r<10; r++) 
        {
          sprintf(c_R_png[t][r],"%s/R%d_%s_%s.png",pngdir,r,tbit[t],var);
          sprintf(c_dev_png[t][r],"%s/deviation_R%d_%s_%s.png",pngdir,r,tbit[t],var);
          c_R[t][r]->Print(c_R_png[t][r],"png");
          c_dev[t][r]->Print(c_dev_png[t][r],"png");
        };
      };
    };
    sprintf(c_rsc_png,"%s/rsc_%s.%s",pngdir,var,file_type);
    sprintf(c_rsr_png,"%s/rsr_%s.%s",pngdir,var,file_type);
    c_rsc->Print(c_rsc_png,file_type);
    c_rsr->Print(c_rsr_png,file_type);
    if(specificFill==0 && specificRun==0)
    {
      for(Int_t r=1; r<10; r++)
      {
        sprintf(c_D_png[r],"%s/zdc_minus_vpd_R%d_%s.png",pngdir,r,var);
        sprintf(c_SD_png[0][r],"%s/east_minus_west_R%d_%s.png",pngdir,r,var);
        sprintf(c_SD_png[1][r],"%s/east_minus_coin_R%d_%s.png",pngdir,r,var);
        sprintf(c_SD_png[2][r],"%s/west_minus_coin_R%d_%s.png",pngdir,r,var);
        sprintf(c_mean_png[r],"%s/mean_R%d_%s.png",pngdir,r,var);
        sprintf(c_RD_png[r],"%s/rsc_minus_mul_R%d_%s.png",pngdir,r,var);
        c_D[r]->Print(c_D_png[r],"png");
        for(Int_t x=0; x<3; x++)
          c_SD[x][r]->Print(c_SD_png[x][r],"png");
        c_mean[r]->Print(c_mean_png[r],"png");
        c_RD[r]->Print(c_RD_png[r],"png");
      };
      char spin_pat_png[64];
      sprintf(spin_pat_png,"%s/spin_pat_%s.png",pngdir,var);
      c_spin_pat->Print(spin_pat_png,"png");
    };
    if(specificFill==0 && specificRun==0 && !strcmp(var,"i"))
    {
      for(Int_t t=0; t<3; t++)
      {
        c_rate_fac[t]->Print(c_rate_fac_png[t],"png");
      };
      sprintf(c_rate_rsr_png,"%s/rate_rsr.png",pngdir);
      sprintf(c_mul_compare_raw_png,"%s/mul_compare_raw.png",pngdir);
      sprintf(c_rsc_compare_raw_png,"%s/rsc_compare_raw.png",pngdir);
      sprintf(c_rsc_compare_mul_png,"%s/rsc_compare_mul.png",pngdir);
      c_rate_rsr->Print(c_rate_rsr_png,"png");
      c_mul_compare_raw->Print(c_mul_compare_raw_png,"png");
      c_rsc_compare_raw->Print(c_rsc_compare_raw_png,"png");
      c_rsc_compare_mul->Print(c_rsc_compare_mul_png,"png");
    };
  };


  // print cuts
  for(Int_t s=0; s<5; s++)
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        printf("raw_cut[%d][%d][%d]=%s\n",t,c,s,raw_cut[t][c][s]);
      };
    };
  };
  for(Int_t s=0; s<4; s++) printf("tot_cut[%d]=%s\n",s,tot_cut[s]);
  printf("%s created\n",outname);
  printf("it's best to use the TBrowser to look at objects in this file\n");



  // builds matrix tree "matx"
  // (named after bXing vs. run index weighted by, e.g., rellum plot,
  //  which is studied as a matrix using SVD)
  TFile * matrix_file;
  char matrix_file_n[128];
  if(specificRun>0 || specificFill>0)
  {
    // set filename
    if(specificRun>0)
      sprintf(matrix_file_n,"matrix/rootfiles/matxR%d.root",specificRun);
    else
      sprintf(matrix_file_n,"matrix/rootfiles/matxF%d.root",specificFill);
    matrix_file = new TFile(matrix_file_n,"RECREATE");
  }
  Double_t raw_cont;
  Double_t acc_cont;
  Double_t mul_cont;
  Double_t fac_cont;
  Double_t rsc_cont;
  Double_t rsr_cont;
  Double_t R_cont[10];
  Int_t tbit_set,cbit_set,bx_set;
  Bool_t kicked_set;
  TTree * matx = new TTree("matx","matx");
  matx->Branch("i",&specificI,"i/I");
  matx->Branch("fi",&specificFI,"fi/I");
  matx->Branch("runnum",&specificRunMatx,"runnum/I");
  matx->Branch("fill",&specificFillMatx,"fill/I");
  matx->Branch("t",&specificT,"t/D");
  matx->Branch("tbit",&tbit_set,"tbit/I"); // (see definition above)
  matx->Branch("cbit",&cbit_set,"cbit/I");
  matx->Branch("bx",&bx_set,"bx/I");
  matx->Branch("kicked",&kicked_set,"kicked/O");
  matx->Branch("raw",&raw_cont,"raw/D");
  matx->Branch("acc",&acc_cont,"acc/D");
  matx->Branch("mul",&mul_cont,"mul/D");
  matx->Branch("fac",&fac_cont,"fac/D");
  matx->Branch("rsc",&rsc_cont,"rsc/D"); // same for all cbits
  matx->Branch("rsr",&rsr_cont,"rsr/D"); // same for all cbits
  /*
  matx->Branch("R1",&R_cont[1],"R1/D");
  matx->Branch("R2",&R_cont[2],"R2/D");
  matx->Branch("R3",&R_cont[3],"R3/D");
  matx->Branch("R4",&R_cont[4],"R4/D");
  matx->Branch("R5",&R_cont[5],"R5/D");
  matx->Branch("R6",&R_cont[6],"R6/D");
  matx->Branch("R7",&R_cont[7],"R7/D");
  matx->Branch("R8",&R_cont[8],"R8/D");
  matx->Branch("R9",&R_cont[9],"R9/D");
  */
  if(specificRun>0 || specificFill>0)
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        tbit_set = t;
        cbit_set = c;
        // loop through histograms (all should have same number of bins)
        for(Int_t cc=1; cc<=mul_d[t][c][4]->GetNbinsX(); cc++)
        {
          bx_set = cc;
          kicked_set = kicked_arr[bx_set-1];
          raw_cont = raw_d[t][c][4]->GetBinContent(cc);
          acc_cont = acc_d[t][c][4]->GetBinContent(cc);
          mul_cont = mul_d[t][c][4]->GetBinContent(cc);
          fac_cont = fac_d[t][c][4]->GetBinContent(cc);
          rsc_cont = rsc_d[t][4]->GetBinContent(cc);
          rsr_cont = rsr_d[t][4]->GetBinContent(cc);
          /*
          for(Int_t rr=1; rr<=9; rr++)
            R_cont[rr] = R_mul_d[t][c][rr]->GetBinContent(cc);
            */
          
          matx->Fill();


          // print columns: run index - fill index - tbit - cbit - bx -
          //                raw - acc - mul - fac - R1 - ... - R9
          /*
          gSystem->RedirectOutput("matrix/matrix_out","a");
          printf("%d %d %f %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
            specificI,specificFI,specificT,t,c,cc,
            raw_cont,acc_cont,mul_cont,fac_cont,
            R_cont[1],R_cont[2],R_cont[3],R_cont[4],R_cont[5],R_cont[6],R_cont[7],R_cont[8],R_cont[9]);
          gSystem->RedirectOutput(0);
          */
        };
      };
    };
    matx->Write("matx");
    sprintf(matrix_file_n,"%s written.\n",matrix_file_n);
    printf(matrix_file_n);
  };


  /*
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t s=0; s<4; s++)
    {
      mul_compare_raw[t][s]->Write();
    };
  };
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t s=0; s<4; s++)
    {
      rsc_compare_raw[t][s]->Write();
    };
  };
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t s=0; s<4; s++)
    {
      rsc_compare_mul[t][s]->Write();
    };
  };
  */


  if(!strcmp(var,"fi") || !strcmp(var,"i"))
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        spinbit_dev[t][c]->Write();
      };
    };
  };
};
