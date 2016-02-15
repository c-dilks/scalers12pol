// bunch fitting algorithm
// -- see equ/bunch_fitter for details
//
// reads hadded matx root file "rootfiles/all.root"
//
// tbit
//  0 - bbc
//  1 - zdc
//  2 - vpd
// cbit
//  0 - east
//  1 - west
//  2 - coin

void BF(const char * var="mul",
        Int_t numer_tbit=1, Int_t numer_cbit=2, Int_t denom_tbit=2, Int_t denom_cbit=2,
        Bool_t evaluateChiSquare=false, Int_t specificPattern=0)
{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  // open hadded matx tree and check its ordering by running DrawMatrix.C
  TFile * infile = new TFile("rootfiles/all.root","READ");
  TTree * matx = (TTree*) infile->Get("matx");
  TFile * countsfile = new TFile("../counts.root","READ");
  TTree * counts = (TTree*) countsfile->Get("sca");
  TFile * sumsfile = new TFile("../sums.root","READ");
  TTree * sums_tr = (TTree*) sumsfile->Get("sum");
  TFile * polfile = new TFile("../pol.root","READ");
  TTree * pol_tr = (TTree*) polfile->Get("pol");
  //printf("executing DrawMatrix.C to check ordering of matx tree\n");
  //gROOT->ProcessLine(".x DrawMatrix.C");


  // trigger bit character strings (tbit)
  char tbit_s[3][4];
  sprintf(tbit_s[0],"bbc");
  sprintf(tbit_s[1],"zdc");
  sprintf(tbit_s[2],"vpd");


  // combination bit character strings (cbit)
  char cbit_s[3][4];
  sprintf(cbit_s[0],"e");
  sprintf(cbit_s[1],"w");
  sprintf(cbit_s[2],"x");

  char detstr[16];
  printf("\n--------------------------\n");
  printf("ratio r^i = %s%s / %s%s\n",
    tbit_s[numer_tbit],cbit_s[numer_cbit],
    tbit_s[denom_tbit],cbit_s[denom_cbit]);
  sprintf(detstr,"%s%s/%s%s",
    tbit_s[numer_tbit],cbit_s[numer_cbit],
    tbit_s[denom_tbit],cbit_s[denom_cbit]);
  printf("--------------------------\n\n");

  // set matx branch addresses
  Int_t i,fi,tbit,cbit,bx;
  Int_t runnum,fill;
  Double_t t,val;
  matx->SetBranchAddress("i",&i);
  matx->SetBranchAddress("fi",&fi);
  matx->SetBranchAddress("runnum",&runnum);
  matx->SetBranchAddress("fill",&fill);
  matx->SetBranchAddress("t",&t);
  matx->SetBranchAddress("tbit",&tbit);
  matx->SetBranchAddress("cbit",&cbit);
  matx->SetBranchAddress("bx",&bx);
  if(!strcmp(var,"raw")) matx->SetBranchAddress("raw",&val);
  else if(!strcmp(var,"acc")) matx->SetBranchAddress("acc",&val);
  else if(!strcmp(var,"mul")) matx->SetBranchAddress("mul",&val);
  else if(!strcmp(var,"fac")) matx->SetBranchAddress("fac",&val);
  else if(!strcmp(var,"rsc")) matx->SetBranchAddress("rsc",&val);
  else 
  {
    fprintf(stderr,"ERROR: val not valid\n");
    return;
  };

  
  // set pol tree branch addresses
  Int_t fill_p;
  Float_t b_pol,y_pol,b_pol_e,y_pol_e;
  pol_tr->SetBranchAddress("fill",&fill_p);
  pol_tr->SetBranchAddress("b_pol",&b_pol);
  pol_tr->SetBranchAddress("y_pol",&y_pol);
  pol_tr->SetBranchAddress("b_pol_e",&b_pol_e);
  pol_tr->SetBranchAddress("y_pol_e",&y_pol_e);


  // check to see if matx tree is sorted by run index from hadd and 
  // build run number & fill number arrays
  // -- if it's not, you'll need to write a sorting routine
  printf("checking matx...\n");
  Int_t imax_tmp = matx->GetMaximum("i");
  const Int_t IMAX = imax_tmp; // max number of runs
  const Int_t NB = 120; // max number of bunches
  Int_t runnum_arr[IMAX];
  Int_t fill_arr[IMAX];
  Int_t fi_arr[IMAX];
  Int_t index_tmp=0;
  for(Int_t ii=0; ii<matx->GetEntries(); ii++)
  {
    matx->GetEntry(ii);
    if(i != index_tmp)
    {
      if(i == index_tmp+1) index_tmp = i;
      else 
      {
        fprintf(stderr,"ERROR: hadded matx tree not correctly sorted\n");
        return;
      };
      runnum_arr[i-1] = runnum;
      fill_arr[i-1] = fill;
      fi_arr[i-1] = fi;
    };
  };


  // define ratio array "rat" and corresponding uncertainties "unc"
  // -- from documentation:  r^i = rat;  sigma_{r^i} = unc
  Double_t rat[IMAX][NB]; // rat = num / den
  Double_t den[IMAX][NB]; // denominator
  Double_t num[IMAX][NB]; // numerator
  Double_t unc[IMAX][NB]; // variances of rat
  printf("computing r^i for each run and bXing...\n");
  for(Int_t ii=0; ii<IMAX; ii++)
  {
    for(Int_t bb=0; bb<NB; bb++)
    {
      num[ii][bb]=0;
      den[ii][bb]=0;
    };
  };
  for(Int_t k=0; k<matx->GetEntries(); k++)
  {
    matx->GetEntry(k);
    if(tbit==numer_tbit && cbit==numer_cbit) num[i-1][bx-1] = val;
    if(tbit==denom_tbit && cbit==denom_cbit) den[i-1][bx-1] = val;
  };
  for(Int_t ii=0; ii<IMAX; ii++)
  {
    for(Int_t bb=0; bb<NB; bb++)
    {
      if(num[ii][bb]>0 && den[ii][bb]>0) 
      {
        rat[ii][bb]=num[ii][bb]/den[ii][bb];
        unc[ii][bb]=rat[ii][bb]*sqrt(1/num[ii][bb]+1/den[ii][bb]);
      }
      else rat[ii][bb]=0;
    };
  };

  
  // rat vs. bXing
  TH1D * rat_v_bx[IMAX];
  TObjArray * rat_v_bx_arr = new TObjArray();
  char rat_v_bx_n[IMAX][32];
  char rat_v_bx_t[IMAX][128];
  for(Int_t ii=0; ii<IMAX; ii++)
  {
    sprintf(rat_v_bx_n[ii],"rat_v_bx_run%d",ii);
    sprintf(rat_v_bx_t[ii],"%s vs. bXing for i=%d",detstr,ii);
    rat_v_bx[ii] = new TH1D(rat_v_bx_n[ii],rat_v_bx_t[ii],NB,0,NB);
    for(Int_t bb=0; bb<NB; bb++)
    {
      rat_v_bx[ii]->SetBinContent(bb+1,rat[ii][bb]);
      rat_v_bx[ii]->SetBinError(bb+1,unc[ii][bb]);
    };
    rat_v_bx_arr->AddLast(rat_v_bx[ii]);
  };


  // set counts branch addresses and fill blue/yell/kicked arrays
  Int_t i_c,bx_c,blue,yell;
  Bool_t kicked;
  Int_t blue_arr[IMAX][NB]; // spin arrays
  Int_t yell_arr[IMAX][NB];
  Bool_t kicked_arr[IMAX][NB];
  counts->SetBranchAddress("i",&i_c);
  counts->SetBranchAddress("bx",&bx_c);
  counts->SetBranchAddress("blue",&blue);
  counts->SetBranchAddress("yell",&yell);
  counts->SetBranchAddress("kicked",&kicked);
  for(Int_t k=0; k<counts->GetEntries(); k++)
  {
    counts->GetEntry(k);
    blue_arr[i_c-1][bx_c]=blue;
    yell_arr[i_c-1][bx_c]=yell;
    kicked_arr[i_c-1][bx_c]=kicked;
  };


  // set sums tree branch address and fill pattern number arrays
  Int_t i_s,pattern;
  Int_t pattern_arr[IMAX];
  sums_tr->SetBranchAddress("i",&i_s);
  sums_tr->SetBranchAddress("pattern",&pattern);
  for(Int_t ii=0; ii<sums_tr->GetEntries(); ii++)
  {
    sums_tr->GetEntry(ii);
    pattern_arr[i_s-1]=pattern;
  };

  
  // fill spin array H 
  // NOTE: H is automatically set to zero if kicked==true;
  //       and if H==0 for a bXing, it's not included in sigma functions
  //
  // -- see equ/bunch_fitting/H_table.pdf for H_a definitions
  //
  Int_t H[IMAX][10][NB]; // [run] [asym] [bx]
  Int_t hb,hy;
  for(Int_t ii=0; ii<IMAX; ii++)
  {
    for(Int_t bb=0; bb<NB; bb++)
    {
      hb = blue_arr[ii][bb];
      hy = yell_arr[ii][bb];
      if(kicked_arr[ii][bb]==0 && abs(hb*hy)==1)
      {
        H[ii][1][bb] = hy;
        H[ii][2][bb] = hb;
        H[ii][3][bb] = hb * hy;
        H[ii][4][bb] = (hb + hy) / 2;
        H[ii][5][bb] = (1 - hb) * hy / 2;
        H[ii][6][bb] = (1 - hy) * hb / 2;
        H[ii][7][bb] = (1 + hb) * hy / 2;
        H[ii][8][bb] = (hy - hb) / 2;
        H[ii][9][bb] = (1 + hy) * hb / 2;
      }
      else for(Int_t aa=1; aa<10; aa++) H[ii][aa][bb]=0;
    };
  };


  // H vs. bXing
  TH1D * H_v_bx[10][IMAX]; // [asym] [run index]
  TObjArray * H_v_bx_arr[10];
  char H_v_bx_n[10][IMAX][32];
  char H_v_bx_t[10][IMAX][128];
  for(Int_t aa=1; aa<10; aa++)
  {
    H_v_bx_arr[aa] = new TObjArray();
    for(Int_t ii=0; ii<IMAX; ii++)
    {
      sprintf(H_v_bx_n[aa][ii],"H%d_v_bx_%d",aa,ii);
      sprintf(H_v_bx_t[aa][ii],"H%d vs. bXing for i=%d",aa,ii);
      H_v_bx[aa][ii] = new TH1D(H_v_bx_n[aa][ii],H_v_bx_t[aa][ii],NB,0,NB);
      for(Int_t bb=0; bb<NB; bb++)
      {
        H_v_bx[aa][ii]->SetBinContent(bb+1,H[ii][aa][bb]);
      };
      H_v_bx_arr[aa]->AddLast(H_v_bx[aa][ii]);
    };
  };


  // fill Hsum arrays (sum of H over bXings for each run & asymmetry number)
  Int_t Hsum[10][IMAX]; // [asym] [run]
  for(Int_t aa=1; aa<10; aa++)
  {
    for(Int_t ii=0; ii<IMAX; ii++)
    {
      Hsum[aa][ii] = 0;
      for(Int_t bb=0; bb<NB; bb++)
      {
        if(kicked_arr[ii][bb]==0) Hsum[aa][ii] += H[ii][aa][bb];
      };
    };
  };


  // Hsum vs. run index
  TH1D * Hsum_dist[10]; //asym
  char Hsum_dist_n[10][32];
  char Hsum_dist_t[10][128];
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(Hsum_dist_n[aa],"Hsum_a%d",aa);
    sprintf(Hsum_dist_t[aa],"Sum of H_{%d} over filled bXings vs. run index",aa);
    Hsum_dist[aa] = new TH1D(Hsum_dist_n[aa],Hsum_dist_t[aa],IMAX,1,IMAX+1);
    for(Int_t ii=0; ii<IMAX; ii++) Hsum_dist[aa]->SetBinContent(ii+1,Hsum[aa][ii]);
  };


  // compute sigma functions for each run, summing over bXings
  // -- if specificPattern!=0, only let sigma function be nonzero for fills with pattern
  //    no. equal to specificPattern
  Double_t sigma_id[10][IMAX]; // [asym] [run]
  Double_t sigma_r[10][IMAX];
  Double_t sigma_h[10][IMAX];
  Double_t sigma_hr[10][IMAX];
  printf("computing sigma functions...\n");
  if(specificPattern!=0) printf(" [+] specificPattern=%d\n",specificPattern);
  for(Int_t aa=1; aa<10; aa++)
  {
    for(Int_t ii=0; ii<IMAX; ii++)
    {
      sigma_id[aa][ii]=0;
      sigma_r[aa][ii]=0;
      sigma_h[aa][ii]=0;
      sigma_hr[aa][ii]=0;
      if(specificPattern==0 || pattern_arr[ii]==specificPattern)
      {
        for(Int_t bb=0; bb<NB; bb++)
        {
          if(kicked_arr[ii][bb]==0 && H[ii][aa][bb]!=0 && unc[ii][bb]>0)
          {
            sigma_id[aa][ii] += 1 / pow(unc[ii][bb],2);
            sigma_r[aa][ii] += rat[ii][bb] / pow(unc[ii][bb],2);
            sigma_h[aa][ii] += H[ii][aa][bb] / pow(unc[ii][bb],2);
            sigma_hr[aa][ii] += H[ii][aa][bb] * rat[ii][bb] / pow(unc[ii][bb],2);
          };
        };
      };
    };
  };


  // compute sigma_sum functions, which are sums of sigma functions over runs
  Double_t sigma_sum_id[10];
  Double_t sigma_sum_r[10];
  Double_t sigma_sum_h[10];
  Double_t sigma_sum_hr[10];
  for(Int_t aa=1; aa<10; aa++)
  {
    sigma_sum_id[aa] = 0;
    sigma_sum_r[aa] = 0;
    sigma_sum_h[aa] = 0;
    sigma_sum_hr[aa] = 0;
    for(Int_t ii=0; ii<IMAX; ii++)
    {
      sigma_sum_id[aa] += sigma_id[aa][ii];
      sigma_sum_r[aa] += sigma_r[aa][ii];
      sigma_sum_h[aa] += sigma_h[aa][ii];
      sigma_sum_hr[aa] += sigma_hr[aa][ii];
    };
  };


  // sigma function vs. run
  TH1D * sigma_id_dist[10]; // [asym]
  TH1D * sigma_r_dist[10];
  TH1D * sigma_h_dist[10];
  TH1D * sigma_hr_dist[10];
  char sigma_id_dist_n[10][16];
  char sigma_r_dist_n[10][16];
  char sigma_h_dist_n[10][16];
  char sigma_hr_dist_n[10][16];
  char sigma_id_dist_t[10][128];
  char sigma_r_dist_t[10][128];
  char sigma_h_dist_t[10][128];
  char sigma_hr_dist_t[10][128];
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(sigma_id_dist_n[aa],"sigma_id_a%d",aa);
    sprintf(sigma_r_dist_n[aa],"sigma_r_a%d",aa);
    sprintf(sigma_h_dist_n[aa],"sigma_h_a%d",aa);
    sprintf(sigma_hr_dist_n[aa],"sigma_hr_a%d",aa);
    sprintf(sigma_id_dist_t[aa],"#Sigma(1) vs. run // asym=%d // r=%s",aa,detstr);
    sprintf(sigma_r_dist_t[aa],"#Sigma(r) vs. run // asym=%d // r=%s",aa,detstr);
    sprintf(sigma_h_dist_t[aa],"#Sigma(H) vs. run // asym==%d // r=%s",aa,detstr);
    sprintf(sigma_hr_dist_t[aa],"#Sigma(Hr) vs. run // asym==%d // r=%s",aa,detstr);
    sigma_id_dist[aa] = new TH1D(sigma_id_dist_n[aa],sigma_id_dist_t[aa],IMAX,1,IMAX+1);
    sigma_r_dist[aa] = new TH1D(sigma_r_dist_n[aa],sigma_r_dist_t[aa],IMAX,1,IMAX+1);
    sigma_h_dist[aa] = new TH1D(sigma_h_dist_n[aa],sigma_h_dist_t[aa],IMAX,1,IMAX+1);
    sigma_hr_dist[aa] = new TH1D(sigma_hr_dist_n[aa],sigma_hr_dist_t[aa],IMAX,1,IMAX+1);
    for(Int_t ii=0; ii<IMAX; ii++)
    {
      sigma_id_dist[aa]->SetBinContent(ii+1,sigma_id[aa][ii]);
      sigma_r_dist[aa]->SetBinContent(ii+1,sigma_r[aa][ii]);
      sigma_h_dist[aa]->SetBinContent(ii+1,sigma_h[aa][ii]);
      sigma_hr_dist[aa]->SetBinContent(ii+1,sigma_hr[aa][ii]);
    };
  };


  // constant vs. run index and epsilon vs. run index
  // -- if specificPattern!=0, only compute constant & epsilon for fills with pattern
  //    no. equal to specificPattern
  TH1D * cons_dist[10]; // [asym]
  char cons_dist_n[10][16];
  char cons_dist_t[10][64];
  TH1D * epsi_dist[10]; // [asym]
  char epsi_dist_n[10][16];
  char epsi_dist_t[10][64];
  Double_t Sid,Sr,Sh,Shr,cons,epsi;
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(cons_dist_n[aa],"cons_a%d_v_run",aa);
    sprintf(epsi_dist_n[aa],"epsi_a%d_v_run",aa);
    sprintf(cons_dist_t[aa],"c_{%d} vs. run index // r=%s",aa,detstr);
    sprintf(epsi_dist_t[aa],"#varepsilon_{%d} vs. run index // r=%s",aa,detstr);
    cons_dist[aa] = new TH1D(cons_dist_n[aa],cons_dist_t[aa],IMAX,1,IMAX+1);
    epsi_dist[aa] = new TH1D(epsi_dist_n[aa],epsi_dist_t[aa],IMAX,1,IMAX+1);
    for(Int_t ii=0; ii<IMAX; ii++)
    {
      if(specificPattern==0 || pattern_arr[ii]==specificPattern)
      {
        Sid = sigma_id[aa][ii];
        Sr = sigma_r[aa][ii];
        Sh = sigma_h[aa][ii];
        Shr = sigma_hr[aa][ii];
        cons = ( Sh * Shr - Sid * Sr) / ( pow(Sh,2) - pow(Sid,2) );
        epsi = ( Sh * Sr - Sid * Shr ) / ( Sh * Shr - Sid * Sr );
        cons_dist[aa]->SetBinContent(ii+1,cons);
        epsi_dist[aa]->SetBinContent(ii+1,epsi);
      };
    };
  };


  // compute statistical errors on constant and epsilon, vs. run index
  Double_t cons_err,epsi_err;
  for(Int_t aa=1; aa<10; aa++)
  {
    for(Int_t ii=0; ii<IMAX; ii++)
    {
      cons_err = 0; 
      epsi_err = 0;
      if(specificPattern==0 || pattern_arr[ii]==specificPattern)
      {
        // first sum the squares...
        for(Int_t bb=0; bb<NB; bb++)
        {
          if(kicked_arr[ii][bb]==0 && H[ii][aa][bb]!=0 && unc[ii][bb]>0)
          {
            epsi_err += pow(1/unc[ii][bb] * 
              ( ( pow(sigma_h[aa][ii],2) - pow(sigma_id[aa][ii],2) ) * 
                ( sigma_hr[aa][ii] - H[ii][aa][bb]*sigma_r[aa][ii] ) ) / 
              ( pow(sigma_hr[aa][ii]*sigma_h[aa][ii] - sigma_id[aa][ii]*sigma_r[aa][ii],2) ),2);
            cons_err += pow(1/unc[ii][bb] *
              ( H[ii][aa][bb]*sigma_h[aa][ii] - sigma_id[aa][ii] ) /
              ( pow(sigma_h[aa][ii],2) - pow(sigma_id[aa][ii],2) ),2);
          };
        };
        // ... then take the square root of the sum of squares
        cons_err = sqrt(cons_err);
        epsi_err = sqrt(epsi_err);
        cons_dist[aa]->SetBinError(ii+1,cons_err);
        epsi_dist[aa]->SetBinError(ii+1,epsi_err);
      };
    };
  };


  // asymmetry vs. run index for a==1,2,3
  // -- asymmetry = epsilon / polarization
  TH1D * asym_dist[4]; // [asym]
  char asym_dist_n[4][16];
  char asym_dist_t[4][64];
  sprintf(asym_dist_t[1],"%s%s A_{L}^{Y} w.r.t. %s%s vs. run index",
    tbit_s[numer_tbit],cbit_s[numer_cbit],tbit_s[denom_tbit],cbit_s[denom_cbit]);
  sprintf(asym_dist_t[2],"%s%s A_{L}^{B} w.r.t. %s%s vs. run index",
    tbit_s[numer_tbit],cbit_s[numer_cbit],tbit_s[denom_tbit],cbit_s[denom_cbit]);
  sprintf(asym_dist_t[3],"%s%s A_{LL} w.r.t. %s%s vs. run index",
    tbit_s[numer_tbit],cbit_s[numer_cbit],tbit_s[denom_tbit],cbit_s[denom_cbit]);
  Float_t b_pol_use,y_pol_use,b_pol_e_use,y_pol_e_use;
  Double_t bc,be,asymmetry,asymmetry_e;
  for(Int_t aa=1; aa<4; aa++)
  {
    sprintf(asym_dist_n[aa],"asym_a%d_v_run",aa);
    asym_dist[aa] = new TH1D(asym_dist_n[aa],asym_dist_t[aa],IMAX,1,IMAX+1);
    for(Int_t bb=1; bb<=epsi_dist[aa]->GetNbinsX(); bb++)
    {
      bc = epsi_dist[aa]->GetBinContent(bb);
      be = epsi_dist[aa]->GetBinError(bb);
      for(Int_t xx=0; xx<pol_tr->GetEntries(); xx++)
      {
        pol_tr->GetEntry(xx);
        if(fill_p == fill_arr[bb-1])
        {
          b_pol_use = b_pol;
          y_pol_use = y_pol;
          b_pol_e_use = b_pol_e;
          y_pol_e_use = y_pol_e;
        };
      };
      if(be>0)
      {
        if(aa==1)
        {
          asymmetry = bc / y_pol_use;
          asymmetry_e = asymmetry * sqrt( pow(be/bc,2) + pow(y_pol_e_use/y_pol_use,2) );
        }
        else if(aa==2) 
        {
          asymmetry = bc / b_pol_use;
          asymmetry_e = asymmetry * sqrt( pow(be/bc,2) + pow(b_pol_e_use/b_pol_use,2) );
        }
        else if(aa==3) 
        {
          asymmetry = bc / (y_pol_use*b_pol_use);
          asymmetry_e = asymmetry * sqrt( pow(be/bc,2) + pow(b_pol_e_use/b_pol_use,2) + pow(y_pol_e_use/y_pol_use,2) ); 
        };
        asym_dist[aa]->SetBinContent(bb,asymmetry);
        asym_dist[aa]->SetBinError(bb,asymmetry_e);

        //if(aa==1) printf("%d %d %f %f\n",bb,fill_arr[bb-1],b_pol_use,y_pol_use);
      };
    };
  };



  // binary search algorithm for parameter range such that
  Double_t optimal_cons[10]; // [asym]  // constant fit of cons vs. i 
  Double_t optimal_epsi[10]; // [asym]  // constant fit of epsi vs. i
  Double_t error_cons[10]; // [asym] // constant fit error of cons vs. i
  Double_t error_epsi[10]; // [asym] // constant fit error of epsi vs. i
  Double_t lb,ub; // binary search lower & upper bounds
  Double_t cx,ex,chi2_eval;
  for(Int_t aa=1; aa<10; aa++)
  {
    // don't do the constant fits if specificPattern!=0
    if(specificPattern==0)
    {
      cons_dist[aa]->Fit("pol0","Q","",1,IMAX+1);
      epsi_dist[aa]->Fit("pol0","Q","",1,IMAX+1);
      if(aa<4) asym_dist[aa]->Fit("pol0","Q","",1,IMAX+1);
      //optimal_cons[aa] = cons_dist[aa]->GetFunction("pol0")->GetParameter(0);
      //optimal_epsi[aa] = epsi_dist[aa]->GetFunction("pol0")->GetParameter(0);
      error_cons[aa] = cons_dist[aa]->GetFunction("pol0")->GetParError(0);
      error_epsi[aa] = epsi_dist[aa]->GetFunction("pol0")->GetParError(0);
    }
    else
    {
      error_cons[aa] = 1;  // set "fake error" if specificPattern!=0
      error_epsi[aa] = 1;
    };

    // determine overall chi2 minimum for all bXings and runs
    Sid = sigma_sum_id[aa];
    Sr = sigma_sum_r[aa];
    Sh = sigma_sum_h[aa];
    Shr = sigma_sum_hr[aa];
    optimal_cons[aa] = ( Sh * Shr - Sid * Sr) / ( pow(Sh,2) - pow(Sid,2) );
    optimal_epsi[aa] = ( Sh * Sr - Sid * Shr ) / ( Sh * Shr - Sid * Sr );
  };
  // determine minimum chi2 value
  Double_t chi2_min[10];
  for(Int_t aa=1; aa<10; aa++)
  {
    chi2_min[aa] = CalcChi2(aa,optimal_cons[aa],optimal_epsi[aa],IMAX,H,rat,unc);
  };
  Double_t ub_cons[10];  // upper / lower bounds for binary search 
  Double_t ub_epsi[10];
  Double_t lb_cons[10];
  Double_t lb_epsi[10];
  Double_t max_cons[10]; // max range of cons
  Double_t max_epsi[10]; // max range of epsi
  Double_t range_cons[10]; // max_cons - optimal_cons
  Double_t range_epsi[10]; // max_epsi - optimal_epsi
  Double_t max_delta = 5; // chi2 @ chi2 profile max range = chi2_min +/- max_delta
  Double_t max_delta_pad = 0.001; // how close chi2_check-chi2_min is to max_delta in
                                //   order to halt binary search
  Double_t chi2_check;
  Int_t NN;
  Double_t mp;
  if(evaluateChiSquare)
  {
    for(Int_t aa=1; aa<10; aa++)
    {
      // cons upper bound
      chi2_check=0;
      NN = 0;
      while(chi2_check - chi2_min[aa] < max_delta)
      {
        printf("[+] chi2_check-chi2_min = %f, max_delta = %f\n",
          chi2_check-chi2_min[aa],max_delta);
        NN++;
        chi2_check = CalcChi2(aa,
                              optimal_cons[aa] + NN * error_cons[aa],
                              optimal_epsi[aa],
                              IMAX,H,rat,unc);
      };
      printf("[+] chi2_check-chi2_min = %f, max_delta = %f\n",
        chi2_check-chi2_min[aa],max_delta);
      ub_cons[aa] = optimal_cons[aa] + NN * error_cons[aa];

      // cons binary search
      lb_cons[aa] = optimal_cons[aa];
      while(fabs(max_delta - (chi2_check - chi2_min[aa])) >= max_delta_pad)
      {
        mp = lb_cons[aa] + (ub_cons[aa] - lb_cons[aa])/2;
        chi2_check = CalcChi2(aa,mp,optimal_epsi[aa],IMAX,H,rat,unc);
        printf(" --> | max_delta - (chi2_check-chi2_min) | = %f\n",
          fabs(max_delta - (chi2_check-chi2_min[aa])));
        if(chi2_check - chi2_min[aa] > max_delta) ub_cons[aa] = mp;
        else lb_cons[aa] = mp;
      };
      max_cons[aa] = mp;
      printf(" ------> chi2_check-chi2_min = %f\n",chi2_check-chi2_min[aa]);
      printf("         optimal_cons=%f max_cons=%f\n",optimal_cons[aa],max_cons[aa]);
      

      // epsi upper bound
      chi2_check=0;
      NN = 0;
      while(chi2_check - chi2_min[aa] < max_delta)
      {
        NN++;
        chi2_check = CalcChi2(aa,
                              optimal_cons[aa],
                              optimal_epsi[aa] + NN * error_epsi[aa],
                              IMAX,H,rat,unc);
      };
      ub_epsi[aa] = optimal_epsi[aa] + NN * error_epsi[aa];

      // epsi binary search
      lb_epsi[aa] = optimal_epsi[aa];
      while(fabs(max_delta - (chi2_check - chi2_min[aa])) >= max_delta_pad)
      {
        mp = lb_epsi[aa] + (ub_epsi[aa] - lb_epsi[aa])/2;
        chi2_check = CalcChi2(aa,optimal_cons[aa],mp,IMAX,H,rat,unc);
        printf(" --> | max_delta - (chi2_check-chi2_min) | = %f\n",
          fabs(max_delta - (chi2_check-chi2_min[aa])));
        if(chi2_check - chi2_min[aa] > max_delta) ub_epsi[aa] = mp;
        else lb_epsi[aa] = mp;
      };
      max_epsi[aa] = mp;
      printf(" ------> chi2_check-chi2_min = %f\n",chi2_check-chi2_min[aa]);
      printf("         optimal_epsi=%f max_epsi=%f\n",optimal_epsi[aa],max_epsi[aa]);

      range_cons[aa] = fabs(max_cons[aa] - optimal_cons[aa]);
      range_epsi[aa] = fabs(max_epsi[aa] - optimal_epsi[aa]);

      printf("[***] max_cons - optimal_cons = %.15f\n",range_cons[aa]);
      printf("[***] max_epsi - optimal_epsi = %.15f\n",range_epsi[aa]);

    };
  };


  // evaluate chi2 by varying fit parameters
  const Int_t NC = 15; // number of bins around optimal to evaluate chi2
  TH1D * chi2_vs_cons[10];
  char chi2_vs_cons_n[10][32];
  char chi2_vs_cons_t[10][128];
  TH1D * chi2_vs_epsi[10];
  char chi2_vs_epsi_n[10][32];
  char chi2_vs_epsi_t[10][128];
  TH2D * chi2_vs_both[10];
  char chi2_vs_both_n[10][32];
  char chi2_vs_both_t[10][128];
  if(evaluateChiSquare)
  {
    for(Int_t aa=1; aa<10; aa++)
    {
      printf("evaluating chi2 profile for asym %d...\n",aa);
      sprintf(chi2_vs_cons_n[aa],"chi2_vs_cons_a%d",aa);
      sprintf(chi2_vs_epsi_n[aa],"chi2_vs_epsi_a%d",aa);
      sprintf(chi2_vs_both_n[aa],"chi2_vs_both_a%d",aa);
      sprintf(chi2_vs_cons_t[aa],
        "#chi^{2}_{%d} vs. c_{%d} // #varepsilon_{%d}=%f // r=%s;c_{%d}",
        aa,aa,aa,optimal_epsi[aa],detstr,aa);
      sprintf(chi2_vs_epsi_t[aa],
        "#chi^{2}_{%d} vs. #varepsilon_{%d} // c_{%d}=%f // r=%s;#varepsilon_{%d}",
        aa,aa,aa,optimal_cons[aa],detstr,aa);
      sprintf(chi2_vs_both_t[aa],
        "#chi^{2}_{%d} vs. c_{%d} vs. #varepsilon_{%d} // r=%s;#varepsilon_{%d};c_{%d}",
        aa,aa,aa,detstr,aa,aa);
      chi2_vs_cons[aa] = new TH1D(chi2_vs_cons_n[aa],chi2_vs_cons_t[aa],
        NC,optimal_cons[aa]-range_cons[aa],optimal_cons[aa]+range_cons[aa]);
      chi2_vs_epsi[aa] = new TH1D(chi2_vs_epsi_n[aa],chi2_vs_epsi_t[aa],
        NC,optimal_epsi[aa]-range_epsi[aa],optimal_epsi[aa]+range_epsi[aa]);
      chi2_vs_both[aa] = new TH2D(chi2_vs_both_n[aa],chi2_vs_both_t[aa],
        NC,optimal_epsi[aa]-range_epsi[aa],optimal_epsi[aa]+range_epsi[aa],
        NC,optimal_cons[aa]-range_cons[aa],optimal_cons[aa]+range_cons[aa]);


      // vary cons holding epsi fixed
      printf(" stage 1/3\n");
      ex = optimal_epsi[aa];
      for(Int_t bin=1; bin<=NC; bin++)
      {
        cx = chi2_vs_cons[aa]->GetBinCenter(bin);
        chi2_eval = CalcChi2(aa,cx,ex,IMAX,H,rat,unc);
        printf(" -- chi2_eval-chi2_min=%f\n",chi2_eval-chi2_min[aa]);
        chi2_vs_cons[aa]->SetBinContent(bin,chi2_eval);
      };

      // vary epsi holding cons fixed
      printf(" stage 2/3\n");
      cx = optimal_cons[aa];
      for(Int_t bin=1; bin<=NC; bin++)
      {
        ex = chi2_vs_epsi[aa]->GetBinCenter(bin);
        chi2_eval = CalcChi2(aa,cx,ex,IMAX,H,rat,unc);
        printf(" -- chi2_eval-chi2_min=%f\n",chi2_eval-chi2_min[aa]);
        chi2_vs_epsi[aa]->SetBinContent(bin,chi2_eval);
      };

      // vary both cons and epsi
      printf(" stage 3/3\n");
      for(Int_t binx=1; binx<=NC; binx++)
      {
        for(Int_t biny=1; biny<=NC; biny++)
        {
          ex = chi2_vs_both[aa]->GetXaxis()->GetBinCenter(binx);
          cx = chi2_vs_both[aa]->GetYaxis()->GetBinCenter(biny);
          chi2_eval = CalcChi2(aa,cx,ex,IMAX,H,rat,unc);
          printf(" -- chi2_eval-chi2_min=%f\n",chi2_eval-chi2_min[aa]);
          chi2_vs_both[aa]->SetBinContent(binx,biny,chi2_eval);
        };
      };
    };
  };
      


  // write output
  const char outfile_name[64];
  if(specificPattern==0)
    sprintf(outfile_name,"fit_result.%s%s.%s%s.root",
      tbit_s[numer_tbit],cbit_s[numer_cbit],
      tbit_s[denom_tbit],cbit_s[denom_cbit]);
  else
    sprintf(outfile_name,"pats/fit_result.%s%s.%s%s.pat%d.root",
      tbit_s[numer_tbit],cbit_s[numer_cbit],
      tbit_s[denom_tbit],cbit_s[denom_cbit],specificPattern);
  TFile * outfile = new TFile(outfile_name,"RECREATE");
  rat_v_bx_arr->Write("rat_v_bx_arr",TObject::kSingleKey);

  outfile->mkdir("H_v_bx_arrs");
  outfile->cd("/H_v_bx_arrs");
  char H_v_bx_arr_n[10][32];
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(H_v_bx_arr_n[aa],"H_a%d_v_bx_arr",aa);
    H_v_bx_arr[aa]->Write(H_v_bx_arr_n[aa],TObject::kSingleKey);
  };

  outfile->mkdir("Hsum");
  outfile->mkdir("sigma_id");
  outfile->mkdir("sigma_r");
  outfile->mkdir("sigma_h");
  outfile->mkdir("sigma_hr");
  outfile->mkdir("constant");
  outfile->mkdir("epsilon");
  outfile->mkdir("asymmetry");
  for(Int_t aa=1; aa<10; aa++) 
  {
    outfile->cd("/Hsum"); Hsum_dist[aa]->Write();
    outfile->cd("/sigma_id"); sigma_id_dist[aa]->Write();
    outfile->cd("/sigma_r"); sigma_r_dist[aa]->Write();
    outfile->cd("/sigma_h"); sigma_h_dist[aa]->Write();
    outfile->cd("/sigma_hr"); sigma_hr_dist[aa]->Write();
    outfile->cd("/constant"); cons_dist[aa]->Write();
    outfile->cd("/epsilon"); epsi_dist[aa]->Write();
    if(aa<4)
    {
      outfile->cd("/asymmetry");
      asym_dist[aa]->Write();
    };
  };
  if(evaluateChiSquare)
  {
    outfile->mkdir("chi2_profile_constant");
    outfile->mkdir("chi2_profile_epsilon");
    outfile->mkdir("chi2_profile_2d");
    for(Int_t aa=1; aa<10; aa++) 
    {
      outfile->cd("chi2_profile_constant"); chi2_vs_cons[aa]->Write();
      outfile->cd("chi2_profile_epsilon"); chi2_vs_epsi[aa]->Write();
      outfile->cd("chi2_profile_2d"); chi2_vs_both[aa]->Write();
    };
  };


  // bunch fit results tree
  TTree * restr  = new TTree();
  Int_t i_res,fi_res,pattern_res,aa_res;
  Double_t cons_res,cons_err_res;
  Double_t epsi_res,epsi_err_res;
  Double_t asym_res,asym_err_res;
  restr->Branch("i",&i_res,"i/I");
  restr->Branch("fi",&fi_res,"fi/I");
  restr->Branch("aa",&aa_res,"aa/I"); // asymmetry number
  restr->Branch("pattern",&pattern_res,"pattern/I");
  restr->Branch("cons",&cons_res,"cons/D");
  restr->Branch("cons_err",&cons_err_res,"cons_err/D");
  restr->Branch("epsi",&epsi_res,"epsi/D");
  restr->Branch("epsi_err",&epsi_err_res,"epsi_err/D");
  restr->Branch("asym",&asym_res,"asym/D");
  restr->Branch("asym_err",&asym_err_res,"asym_err/D");

  for(Int_t aa=1; aa<4; aa++) 
  {
    for(Int_t bb=1; bb<=IMAX; bb++)
    {
      i_res = bb;
      fi_res = fi_arr[bb-1];
      pattern_res = pattern_arr[bb-1];
      aa_res = aa;

      cons_res = cons_dist[aa]->GetBinContent(bb);
      cons_err_res = cons_dist[aa]->GetBinError(bb);
      epsi_res = epsi_dist[aa]->GetBinContent(bb);
      epsi_err_res = epsi_dist[aa]->GetBinError(bb);
      asym_res = asym_dist[aa]->GetBinContent(bb);
      asym_err_res = asym_dist[aa]->GetBinError(bb);

      restr->Fill();
    };
  };

  outfile->cd();
  restr->Write("restr");


  printf("%s created\n",outfile_name);
}


//Double_t CalcChi2(Double_t cx0, Double_t ex0, 
                  //Int_t imax0, Int_t * H0, Double_t * rat0, Double_t * unc0)
Double_t CalcChi2(Int_t aa0, Double_t cx0, Double_t ex0, 
                  Int_t imax0, 
                  Int_t H0[][10][120], 
                  Double_t rat0[][120], 
                  Double_t unc0[][120])
{
  Double_t value = 0;
  for(Int_t ii=0; ii<imax0; ii++)
  {
    for(Int_t bb=0; bb<120; bb++)
    {
      if(H0[aa0][ii][bb]!=0 && unc0[ii][bb]>0)
      {
        value += pow((cx0*(1+H0[aa0][ii][bb]*ex0)-rat0[ii][bb])/unc0[ii][bb], 2);
      };
    };
  };
  return value;
};
