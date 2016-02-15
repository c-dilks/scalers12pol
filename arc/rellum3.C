void rellum3(const char * var="i", const char * trigger="zdc",Bool_t printPNGs=0)
{
  // read counts.root file
  TFile * infile = new TFile("counts.root","READ");
  TTree * s = (TTree*) infile->Get("sca");
  //s->Print();

  // define independent variable bounds
  Int_t var_l = s->GetMinimum(var);
  Int_t var_h = s->GetMaximum(var);
  var_h++; // fencepost
  Int_t var_bins = var_h-var_l;
  if(!strcmp(var,"t")) var_bins=800;


  // set up trigger names for reading
  if(strcmp(trigger,"bbc") && 
     strcmp(trigger,"zdc") &&
     strcmp(trigger,"vpd"))
  { 
    fprintf(stderr,"error: \"%s\" is an invalid trigger name\n",trigger);
    return;
  };
  char trigger_e[16];
  char trigger_w[16];
  char trigger_x[16];
  sprintf(trigger_e,"%se",trigger);
  sprintf(trigger_w,"%sw",trigger);
  sprintf(trigger_x,"%sx",trigger);

  // set branch addresses
  Int_t index,runnum,fill,bx;
  Double_t TTe,TTw,TTx;
  Int_t blue,yell;
  s->SetBranchAddress("i",&index);
  s->SetBranchAddress("runnum",&runnum);
  s->SetBranchAddress("fill",&fill);
  s->SetBranchAddress("bx",&bx);
  s->SetBranchAddress(trigger_e,&TTe);
  s->SetBranchAddress(trigger_w,&TTw);
  s->SetBranchAddress(trigger_x,&TTx);
  s->SetBranchAddress("blue",&blue);
  s->SetBranchAddress("yell",&yell);
  

  // define distributions "raw" -- raw number of triggers binned by independent variable 
  //                                for same and diff spin cases
  char TTe_same_raw_n[128]; sprintf(TTe_same_raw_n,"%s_same",trigger_e);
  char TTe_diff_raw_n[128]; sprintf(TTe_diff_raw_n,"%s_diff",trigger_e);
  char TTe_both_raw_n[128]; sprintf(TTe_both_raw_n,"%s_both",trigger_e);
  char TTw_same_raw_n[128]; sprintf(TTw_same_raw_n,"%s_same",trigger_w);
  char TTw_diff_raw_n[128]; sprintf(TTw_diff_raw_n,"%s_diff",trigger_w);
  char TTw_both_raw_n[128]; sprintf(TTw_both_raw_n,"%s_both",trigger_w);
  char TTx_same_raw_n[128]; sprintf(TTx_same_raw_n,"%s_same",trigger_x);
  char TTx_diff_raw_n[128]; sprintf(TTx_diff_raw_n,"%s_diff",trigger_x);
  char TTx_both_raw_n[128]; sprintf(TTx_both_raw_n,"%s_both",trigger_x);
  char trig_same_raw_n[128]; sprintf(trig_same_raw_n,"%s_same","trig");
  char trig_diff_raw_n[128]; sprintf(trig_diff_raw_n,"%s_diff","trig");
  char trig_both_raw_n[128]; sprintf(trig_both_raw_n,"%s_both","trig");
  char spin_same_raw_n[128]; sprintf(spin_same_raw_n,"%s_same","spin");
  char spin_diff_raw_n[128]; sprintf(spin_diff_raw_n,"%s_diff","spin");
  char spin_both_raw_n[128]; sprintf(spin_both_raw_n,"%s_both","spin");
  char totBx_raw_n[128]; sprintf(totBx_raw_n,"%s_same","totBx");
  char TTe_raw_t[128]; sprintf(TTe_raw_t,"%s weighted by raw %s scalers",var,trigger_e); 
  char TTw_raw_t[128]; sprintf(TTw_raw_t,"%s weighted by raw %s scalers",var,trigger_w);
  char TTx_raw_t[128]; sprintf(TTx_raw_t,"%s weighted by raw %s scalers",var,trigger_x);
  char trig_raw_t[128]; sprintf(trig_raw_t,"%s weighted by raw trig scalers",var);
  char spin_raw_t[128]; sprintf(spin_raw_t,"%s weighted by no. of buckets",var);
  char totBx_raw_t[128]; sprintf(totBx_raw_t,"%s weighted by total bXings",var);

  TH1D * TTe_same_raw = new TH1D(TTe_same_raw_n,TTe_raw_t,var_bins,var_l,var_h);
  TH1D * TTe_diff_raw = new TH1D(TTe_diff_raw_n,TTe_raw_t,var_bins,var_l,var_h);
  TH1D * TTe_both_raw = new TH1D(TTe_both_raw_n,TTe_raw_t,var_bins,var_l,var_h);
  TH1D * TTw_same_raw = new TH1D(TTw_same_raw_n,TTw_raw_t,var_bins,var_l,var_h);
  TH1D * TTw_diff_raw = new TH1D(TTw_diff_raw_n,TTw_raw_t,var_bins,var_l,var_h);
  TH1D * TTw_both_raw = new TH1D(TTw_both_raw_n,TTw_raw_t,var_bins,var_l,var_h);
  TH1D * TTx_same_raw = new TH1D(TTx_same_raw_n,TTx_raw_t,var_bins,var_l,var_h);
  TH1D * TTx_diff_raw = new TH1D(TTx_diff_raw_n,TTx_raw_t,var_bins,var_l,var_h);
  TH1D * TTx_both_raw = new TH1D(TTx_both_raw_n,TTx_raw_t,var_bins,var_l,var_h);
  TH1D * trig_same_raw = new TH1D(trig_same_raw_n,trig_raw_t,var_bins,var_l,var_h);
  TH1D * trig_diff_raw = new TH1D(trig_diff_raw_n,trig_raw_t,var_bins,var_l,var_h);
  TH1D * trig_both_raw = new TH1D(trig_both_raw_n,trig_raw_t,var_bins,var_l,var_h);
  TH1D * spin_same_raw = new TH1D(spin_same_raw_n,spin_raw_t,var_bins,var_l,var_h);
  TH1D * spin_diff_raw = new TH1D(spin_diff_raw_n,spin_raw_t,var_bins,var_l,var_h);
  TH1D * spin_both_raw = new TH1D(spin_both_raw_n,spin_raw_t,var_bins,var_l,var_h);
  TH1D * totBx_raw = new TH1D(totBx_raw_n,totBx_raw_t,var_bins,var_l,var_h);
  

  // define corrected distributions "corr" -- to be filled later after projections
  // -- these are either accidentals only or accidentals + multiples counts
  // -- for correction factor, see "factor" histograms below
  char TTe_same_corr_n[128]; sprintf(TTe_same_corr_n,"%s_same_corr",trigger_e);
  char TTe_diff_corr_n[128]; sprintf(TTe_diff_corr_n,"%s_diff_corr",trigger_e);
  char TTe_both_corr_n[128]; sprintf(TTe_both_corr_n,"%s_both_corr",trigger_e);
  char TTw_same_corr_n[128]; sprintf(TTw_same_corr_n,"%s_same_corr",trigger_w);
  char TTw_diff_corr_n[128]; sprintf(TTw_diff_corr_n,"%s_diff_corr",trigger_w);
  char TTw_both_corr_n[128]; sprintf(TTw_both_corr_n,"%s_both_corr",trigger_w);
  char TTx_same_corr_n[128]; sprintf(TTx_same_corr_n,"%s_same_corr",trigger_x);
  char TTx_diff_corr_n[128]; sprintf(TTx_diff_corr_n,"%s_diff_corr",trigger_x);
  char TTx_both_corr_n[128]; sprintf(TTx_both_corr_n,"%s_both_corr",trigger_x);
  char TTe_corr_t[128]; sprintf(TTe_corr_t,"%s weighted by corrected %s scalers",var,trigger_e); 
  char TTw_corr_t[128]; sprintf(TTw_corr_t,"%s weighted by corrected %s scalers",var,trigger_w);
  char TTx_corr_t[128]; sprintf(TTx_corr_t,"%s weighted by corrected %s scalers",var,trigger_x);
  TH1D * TTe_same_corr = new TH1D(TTe_same_corr_n,TTe_corr_t,var_bins,var_l,var_h);
  TH1D * TTe_diff_corr = new TH1D(TTe_diff_corr_n,TTe_corr_t,var_bins,var_l,var_h);
  TH1D * TTe_both_corr = new TH1D(TTe_both_corr_n,TTe_corr_t,var_bins,var_l,var_h);
  TH1D * TTw_same_corr = new TH1D(TTw_same_corr_n,TTw_corr_t,var_bins,var_l,var_h);
  TH1D * TTw_diff_corr = new TH1D(TTw_diff_corr_n,TTw_corr_t,var_bins,var_l,var_h);
  TH1D * TTw_both_corr = new TH1D(TTw_both_corr_n,TTw_corr_t,var_bins,var_l,var_h);
  TH1D * TTx_same_corr = new TH1D(TTx_same_corr_n,TTx_corr_t,var_bins,var_l,var_h);
  TH1D * TTx_diff_corr = new TH1D(TTx_diff_corr_n,TTx_corr_t,var_bins,var_l,var_h);
  TH1D * TTx_both_corr = new TH1D(TTx_both_corr_n,TTx_corr_t,var_bins,var_l,var_h);


  // --- projections over indepenedent variable, weighted by raw scalers
  char TTe_same_raw_cut[256];
  char TTe_diff_raw_cut[256];
  char TTe_both_raw_cut[256];
  char TTw_same_raw_cut[256];
  char TTw_diff_raw_cut[256];
  char TTw_both_raw_cut[256];
  char TTx_same_raw_cut[256];
  char TTx_diff_raw_cut[256];
  char TTx_both_raw_cut[256];
  char trig_same_raw_cut[256];
  char trig_diff_raw_cut[256];
  char trig_both_raw_cut[256];
  char spin_same_raw_cut[256];
  char spin_diff_raw_cut[256];
  char spin_both_raw_cut[256];
  char totBx_raw_cut[256];

  sprintf(TTe_same_raw_cut,"%s*(blue*yell!=0 && blue==yell)",trigger_e);
  sprintf(TTe_diff_raw_cut,"%s*(blue*yell!=0 && blue!=yell)",trigger_e);
  sprintf(TTe_both_raw_cut,"%s*(blue*yell!=0)",trigger_e);
  sprintf(TTw_same_raw_cut,"%s*(blue*yell!=0 && blue==yell)",trigger_w);
  sprintf(TTw_diff_raw_cut,"%s*(blue*yell!=0 && blue!=yell)",trigger_w);
  sprintf(TTw_both_raw_cut,"%s*(blue*yell!=0)",trigger_w);
  sprintf(TTx_same_raw_cut,"%s*(blue*yell!=0 && blue==yell)",trigger_x);
  sprintf(TTx_diff_raw_cut,"%s*(blue*yell!=0 && blue!=yell)",trigger_x);
  sprintf(TTx_both_raw_cut,"%s*(blue*yell!=0)",trigger_x);
  sprintf(trig_same_raw_cut,"(%s+%s+%s)*(blue*yell!=0 && blue==yell)",
                              trigger_e,trigger_w,trigger_x);
  sprintf(trig_diff_raw_cut,"(%s+%s+%s)*(blue*yell!=0 && blue!=yell)",
                              trigger_e,trigger_w,trigger_x);
  sprintf(trig_both_raw_cut,"(%s+%s+%s)*(blue*yell!=0)",
                              trigger_e,trigger_w,trigger_x);
  sprintf(spin_same_raw_cut,"1*(blue*yell!=0 && blue==yell)");
  sprintf(spin_diff_raw_cut,"1*(blue*yell!=0 && blue!=yell)");
  sprintf(spin_both_raw_cut,"1*(blue*yell!=0)");
  sprintf(totBx_raw_cut,"tot_bx*(blue*yell!=0)"); // only bXings with possible collision

  s->Project(TTe_same_raw_n,var,TTe_same_raw_cut);
  s->Project(TTe_diff_raw_n,var,TTe_diff_raw_cut);
  s->Project(TTe_both_raw_n,var,TTe_both_raw_cut);
  s->Project(TTw_same_raw_n,var,TTw_same_raw_cut);
  s->Project(TTw_diff_raw_n,var,TTw_diff_raw_cut);
  s->Project(TTw_both_raw_n,var,TTw_both_raw_cut);
  s->Project(TTx_same_raw_n,var,TTx_same_raw_cut);
  s->Project(TTx_diff_raw_n,var,TTx_diff_raw_cut);
  s->Project(TTx_both_raw_n,var,TTx_both_raw_cut);
  s->Project(trig_same_raw_n,var,trig_same_raw_cut);
  s->Project(trig_diff_raw_n,var,trig_diff_raw_cut);
  s->Project(trig_both_raw_n,var,trig_both_raw_cut);
  s->Project(spin_same_raw_n,var,spin_same_raw_cut);
  s->Project(spin_diff_raw_n,var,spin_diff_raw_cut);
  s->Project(spin_both_raw_n,var,spin_both_raw_cut);
  s->Project(totBx_raw_n,var,totBx_raw_cut);


  // accidentals and multiples corrections
  Double_t ne, nw, nx, total_bx; // scaled counts
  Double_t pe, pw, px; // physical process probabilities
  Double_t ae, aw, ax; // counts corrected for accidentals (stage 1)
  Double_t me, mw, mx; // counts corrected for multiples (stage 2)
  for(Int_t n=1; n<=var_bins; n++)
  {
    // scaled numbers
    ne = TTe_both_raw->GetBinContent(n);
    nw = TTw_both_raw->GetBinContent(n);
    nx = TTx_both_raw->GetBinContent(n);
    total_bx = totBx_raw->GetBinContent(n);


    // physical process probabilities
    pe = (ne-nx)/(total_bx-nw);
    pw = (nw-nx)/(total_bx-ne);
    px = (nx-(ne*nw)/total_bx)/(total_bx+nx-ne-nw);

    /*
    printf("probabilities: pe=%f,pe);
    printf("probabilities: pw=%f,pw);
    printf("probabilities: px=%f,px);
    */

    // accidentals corrections (stage 1)
    ae = pe*total_bx;
    aw = pw*total_bx;
    ax = px*total_bx;
    

    // multiples corrections (stage 2)
    me = -1*total_bx*log(1-pe);
    mw = -1*total_bx*log(1-pw);
    mx = -1*total_bx*log(1-px);


    // fill corrected dists
    /*
    TTe_both_corr->SetBinContent(n,ae); // accidentals only
    TTw_both_corr->SetBinContent(n,aw);
    TTx_both_corr->SetBinContent(n,ax);
    */
    TTe_both_corr->SetBinContent(n,me); // accidentals + multiples
    TTw_both_corr->SetBinContent(n,mw);
    TTx_both_corr->SetBinContent(n,mx);
  };
  /*
  new TCanvas(); TTe_both_corr->Draw();
  new TCanvas(); TTe_both_raw->Draw();
  */
  printf("raw max: %.2f  corr max: %.2f\n",TTe_both_raw->GetMaximum(),TTe_both_corr->GetMaximum());


  // correction factor calculation -------------------------------------------- is this right////????
  char TTe_factor_n[128]; sprintf(TTe_factor_n,"%s_factor",trigger_e);
  char TTw_factor_n[128]; sprintf(TTw_factor_n,"%s_factor",trigger_w);
  char TTx_factor_n[128]; sprintf(TTx_factor_n,"%s_factor",trigger_x);
  char TTe_factor_t[128]; sprintf(TTe_factor_t,"%s correction factor vs. %s",trigger_e,var);
  char TTw_factor_t[128]; sprintf(TTw_factor_t,"%s correction factor vs. %s",trigger_w,var);
  char TTx_factor_t[128]; sprintf(TTx_factor_t,"%s correction factor vs. %s",trigger_x,var);
  TH1D * TTe_factor = new TH1D(TTe_factor_n,TTe_factor_t,var_bins,var_l,var_h);
  TH1D * TTw_factor = new TH1D(TTw_factor_n,TTw_factor_t,var_bins,var_l,var_h);
  TH1D * TTx_factor = new TH1D(TTx_factor_n,TTx_factor_t,var_bins,var_l,var_h);
  TTe_factor->Divide(TTe_both_corr,TTe_both_raw,1.0,1.0);
  TTw_factor->Divide(TTw_both_corr,TTw_both_raw,1.0,1.0);
  TTx_factor->Divide(TTx_both_corr,TTx_both_raw,1.0,1.0);
  
  Int_t sf; // drawing scale factor
  if(printPNGs) sf=3;
  else sf=1;
  TCanvas * canv_factor = new TCanvas("canv_factor","canv_factor",1100*sf,700*sf);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.08);
  TTe_factor->GetYaxis()->SetLabelSize(0.08);
  TTe_factor->GetXaxis()->SetLabelSize(0.08);
  TTw_factor->GetYaxis()->SetLabelSize(0.08);
  TTw_factor->GetXaxis()->SetLabelSize(0.08);
  TTx_factor->GetYaxis()->SetLabelSize(0.08);
  TTx_factor->GetXaxis()->SetLabelSize(0.08);
  canv_factor->Divide(1,3);
  for(Int_t i=1; i<=3; i++) canv_factor->cd(i)->SetGrid(1,1);
  canv_factor->cd(1); TTe_factor->Draw();
  canv_factor->cd(2); TTw_factor->Draw();
  canv_factor->cd(3); TTx_factor->Draw();
  char canv_factor_print[128];
  sprintf(canv_factor_print,"pngs/canv_factor_%s_vs_%s.png",trigger,var);
  if(printPNGs) canv_factor->Print(canv_factor_print,"png");
  return;



  // compute raw relative luminosity R3
  char RTTe_raw_t[128]; sprintf(RTTe_raw_t,"raw R_{%s} :: %s",trigger_e,var);
  char RTTw_raw_t[128]; sprintf(RTTw_raw_t,"raw R_{%s} :: %s",trigger_w,var);
  char RTTx_raw_t[128]; sprintf(RTTx_raw_t,"raw R_{%s} :: %s",trigger_x,var);
  char Rtrig_raw_t[128]; sprintf(Rtrig_raw_t,"raw R_{trig} :: %s",var);
  char Rspin_t[128]; sprintf(Rspin_t,"R_{bx} :: %s",var);
  TH1D * RTTe_raw = new TH1D("RTTe_raw",RTTe_raw_t,var_bins,var_l,var_h);
  TH1D * RTTw_raw = new TH1D("RTTw_raw",RTTw_raw_t,var_bins,var_l,var_h);
  TH1D * RTTx_raw = new TH1D("RTTx_raw",RTTx_raw_t,var_bins,var_l,var_h);
  TH1D * Rtrig_raw = new TH1D("Rtrig_raw",Rtrig_raw_t,var_bins,var_l,var_h);
  TH1D * Rspin = new TH1D("Rspin",Rspin_t,var_bins,var_l,var_h);
  RTTe_raw->Divide(TTe_same_raw,TTe_diff_raw,1.0,1.0);
  RTTw_raw->Divide(TTw_same_raw,TTw_diff_raw,1.0,1.0);
  RTTx_raw->Divide(TTx_same_raw,TTx_diff_raw,1.0,1.0);
  Rtrig_raw->Divide(trig_same_raw,trig_diff_raw,1.0,1.0);
  Rspin->Divide(spin_same_raw,spin_diff_raw,1.0,1.0);


  // compute corrected relative luminosity R3
  char RTTe_corr_t[128]; sprintf(RTTe_corr_t,"corrected R_{%s} :: %s",trigger_e,var);
  char RTTw_corr_t[128]; sprintf(RTTw_corr_t,"corrected R_{%s} :: %s",trigger_w,var);
  char RTTx_corr_t[128]; sprintf(RTTx_corr_t,"corrected R_{%s} :: %s",trigger_x,var);
  char Rtrig_corr_t[128]; sprintf(Rtrig_corr_t,"corrected R_{trig} :: %s",var);
  TH1D * RTTe_corr = new TH1D("RTTe_corr",RTTe_corr_t,var_bins,var_l,var_h);
  TH1D * RTTw_corr = new TH1D("RTTw_corr",RTTw_corr_t,var_bins,var_l,var_h);
  TH1D * RTTx_corr = new TH1D("RTTx_corr",RTTx_corr_t,var_bins,var_l,var_h);
  TH1D * Rtrig_corr = new TH1D("Rtrig_corr",Rtrig_corr_t,var_bins,var_l,var_h);
  RTTe_corr->Divide(TTe_same_corr,TTe_diff_corr,1.0,1.0);
  RTTw_corr->Divide(TTw_same_corr,TTw_diff_corr,1.0,1.0);
  RTTx_corr->Divide(TTx_same_corr,TTx_diff_corr,1.0,1.0);


  // Q variable Q:=|R-1|
  /*
  char QTTe_t[128]; sprintf(QTTe_t,"Q_{%s} :: %s",trigger_e,var);
  char QTTw_t[128]; sprintf(QTTw_t,"Q_{%s} :: %s",trigger_w,var);
  char QTTx_t[128]; sprintf(QTTx_t,"Q_{%s} :: %s",trigger_x,var);
  char Qtrig_t[128]; sprintf(Qtrig_t,"Q_{trig} :: %s",var);
  char Qspin_t[128]; sprintf(Qspin_t,"Q_{spin} :: %s",var);
  TH1D * QTTe = new TH1D("QTTe",QTTe_t,var_bins,var_l,var_h);
  TH1D * QTTw = new TH1D("QTTw",QTTw_t,var_bins,var_l,var_h);
  TH1D * QTTx = new TH1D("QTTx",QTTx_t,var_bins,var_l,var_h);
  TH1D * Qtrig = new TH1D("Qtrig",Qtrig_t,var_bins,var_l,var_h);
  TH1D * Qspin = new TH1D("Qspin",Qspin_t,var_bins,var_l,var_h);
  for(Int_t r=1; r<=RTTe->GetNbinsX(); r++)
  {
    QTTe->SetBinContent(r,fabs(RTTe->GetBinContent(r)-1));
    QTTw->SetBinContent(r,fabs(RTTw->GetBinContent(r)-1));
    QTTx->SetBinContent(r,fabs(RTTx->GetBinContent(r)-1));
    Qtrig->SetBinContent(r,fabs(Rtrig->GetBinContent(r)-1));
    Qspin->SetBinContent(r,fabs(Rspin->GetBinContent(r)-1));
  };

  // Average rates (same spin)
  char ASTTe_t[128]; sprintf(ASTTe_t,"AS_{%s} :: %s",trigger_e,var);
  char ASTTw_t[128]; sprintf(ASTTw_t,"AS_{%s} :: %s",trigger_w,var);
  char ASTTx_t[128]; sprintf(ASTTx_t,"AS_{%s} :: %s",trigger_x,var);
  char AStrig_t[128]; sprintf(AStrig_t,"AS_{trig} :: %s",var);
  TH1D * ASTTe = new TH1D("ASTTe",ASTTe_t,var_bins,var_l,var_h);
  TH1D * ASTTw = new TH1D("ASTTw",ASTTw_t,var_bins,var_l,var_h);
  TH1D * ASTTx = new TH1D("ASTTx",ASTTx_t,var_bins,var_l,var_h);
  TH1D * AStrig = new TH1D("AStrig",AStrig_t,var_bins,var_l,var_h);
  ASTTe->Divide(TTe_same_raw,spin_same_raw,1.0,1.0);
  ASTTw->Divide(TTw_same_raw,spin_same_raw,1.0,1.0);
  ASTTx->Divide(TTx_same_raw,spin_same_raw,1.0,1.0);
  AStrig->Divide(trig_same_raw,spin_same_raw,1.0,1.0);


  // Average rates (diff spin)
  char ADTTe_t[128]; sprintf(ADTTe_t,"AD_{%s} :: %s",trigger_e,var);
  char ADTTw_t[128]; sprintf(ADTTw_t,"AD_{%s} :: %s",trigger_w,var);
  char ADTTx_t[128]; sprintf(ADTTx_t,"AD_{%s} :: %s",trigger_x,var);
  char ADtrig_t[128]; sprintf(ADtrig_t,"AD_{trig} :: %s",var);
  TH1D * ADTTe = new TH1D("ADTTe",ADTTe_t,var_bins,var_l,var_h);
  TH1D * ADTTw = new TH1D("ADTTw",ADTTw_t,var_bins,var_l,var_h);
  TH1D * ADTTx = new TH1D("ADTTx",ADTTx_t,var_bins,var_l,var_h);
  TH1D * ADtrig = new TH1D("ADtrig",ADtrig_t,var_bins,var_l,var_h);
  ADTTe->Divide(TTe_diff_raw,spin_diff_raw,1.0,1.0);
  ADTTw->Divide(TTw_diff_raw,spin_diff_raw,1.0,1.0);
  ADTTx->Divide(TTx_diff_raw,spin_diff_raw,1.0,1.0);
  ADtrig->Divide(trig_diff_raw,spin_diff_raw,1.0,1.0);

  // sum of S's for different and same spins (SS)
  char SSTTe_t[128]; sprintf(SSTTe_t,"SS_{%s} :: %s",trigger_e,var);
  char SSTTw_t[128]; sprintf(SSTTw_t,"SS_{%s} :: %s",trigger_w,var);
  char SSTTx_t[128]; sprintf(SSTTx_t,"SS_{%s} :: %s",trigger_x,var);
  char SStrig_t[128]; sprintf(SStrig_t,"SS_{trig} :: %s",var);
  char SSspin_t[128]; sprintf(SSspin_t,"SS_{bx} :: %s",var);
  TH1D * SSTTe = new TH1D("SSTTe",SSTTe_t,var_bins,var_l,var_h);
  TH1D * SSTTw = new TH1D("SSTTw",SSTTw_t,var_bins,var_l,var_h);
  TH1D * SSTTx = new TH1D("SSTTx",SSTTx_t,var_bins,var_l,var_h);
  TH1D * SStrig = new TH1D("SStrig",SStrig_t,var_bins,var_l,var_h);
  TH1D * SSspin = new TH1D("SSspin",SSspin_t,var_bins,var_l,var_h);
  SSTTe->Add(TTe_same_raw,TTe_diff_raw);
  SSTTw->Add(TTw_same_raw,TTw_diff_raw);
  SSTTx->Add(TTx_same_raw,TTx_diff_raw);
  SStrig->Add(trig_same_raw,trig_diff_raw);
  SSspin->Add(spin_same_raw,spin_diff_raw);

  // Average rates (both spin)
  char ABTTe_t[128]; sprintf(ABTTe_t,"AB_{%s} :: %s",trigger_e,var);
  char ABTTw_t[128]; sprintf(ABTTw_t,"AB_{%s} :: %s",trigger_w,var);
  char ABTTx_t[128]; sprintf(ABTTx_t,"AB_{%s} :: %s",trigger_x,var);
  char ABtrig_t[128]; sprintf(ABtrig_t,"AB_{trig} :: %s",var);
  TH1D * ABTTe = new TH1D("ABTTe",ABTTe_t,var_bins,var_l,var_h);
  TH1D * ABTTw = new TH1D("ABTTw",ABTTw_t,var_bins,var_l,var_h);
  TH1D * ABTTx = new TH1D("ABTTx",ABTTx_t,var_bins,var_l,var_h);
  TH1D * ABtrig = new TH1D("ABtrig",ABtrig_t,var_bins,var_l,var_h);
  ABTTe->Divide(SSTTe,SSspin,1.0,1.0);
  ABTTw->Divide(SSTTw,SSspin,1.0,1.0);
  ABTTx->Divide(SSTTx,SSspin,1.0,1.0);
  ABtrig->Divide(SStrig,SSspin,1.0,1.0);
  */

  /*
  // draw var weighted by total scalers (S++ :: var)
  TCanvas * canv_same = new TCanvas("canv_same","canv_same",1200,700);
  canv_same->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_same->GetPad(j)->SetGrid(1,1);
  canv_same->cd(1); trig_same_raw->Draw();
  canv_same->cd(4); spin_same_raw->Draw();
  canv_same->cd(2); TTe_same_raw->Draw();
  canv_same->cd(3); TTw_same_raw->Draw();
  canv_same->cd(5); TTx_same_raw->Draw();
  
  // draw var weighted by total scalers (S+- :: var)
  TCanvas * canv_diff = new TCanvas("canv_diff","canv_diff",1200,700);
  canv_diff->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_diff->GetPad(j)->SetGrid(1,1);
  canv_diff->cd(1); trig_diff_raw->Draw();
  canv_diff->cd(4); spin_diff_raw->Draw();
  canv_diff->cd(2); TTe_diff_raw->Draw();
  canv_diff->cd(3); TTw_diff_raw->Draw();
  canv_diff->cd(5); TTx_diff_raw->Draw();

  // draw average rate (AS :: var)
  TCanvas * canv_AS = new TCanvas("canv_AS","canv_AS",1200,700);
  canv_AS->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_AS->GetPad(j)->SetGrid(1,1);
  canv_AS->cd(1); AStrig->Draw();
  canv_AS->cd(2); ASTTe->Draw();
  canv_AS->cd(3); ASTTw->Draw();
  canv_AS->cd(5); ASTTx->Draw();
  
  // draw average rate (AD :: var)
  TCanvas * canv_AD = new TCanvas("canv_AD","canv_AD",1200,700);
  canv_AD->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_AD->GetPad(j)->SetGrid(1,1);
  canv_AD->cd(1); ADtrig->Draw();
  canv_AD->cd(2); ADTTe->Draw();
  canv_AD->cd(3); ADTTw->Draw();
  canv_AD->cd(5); ADTTx->Draw();
  */
  
  // draw average rate (AB :: var)
  /*
  TCanvas * canv_AB = new TCanvas("canv_AB","canv_AB",1200,700);
  canv_AB->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_AB->GetPad(j)->SetGrid(1,1);
  canv_AB->cd(1); ABtrig->Draw();
  canv_AB->cd(2); ABTTe->Draw();
  canv_AB->cd(3); ABTTw->Draw();
  canv_AB->cd(5); ABTTx->Draw();
  */
  
  // draw normaliser (total bx's binned in independent variable "var")
  /*
  TCanvas * canv_N = new TCanvas("canv_N","canv_N",600,400);
  canv_N->SetGrid(1,1);
  spin_same_raw->Draw();
  */


  /*
  // draw Q graphs (Q=|R-1| :: var)
  TCanvas * canv_Q = new TCanvas("canv_Q","canv_Q",1200,700);
  canv_Q->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_Q->GetPad(j)->SetGrid(1,1);
  canv_Q->cd(1); Qtrig->Draw();
  canv_Q->cd(4); Qspin->Draw();
  canv_Q->cd(2); QTTe->Draw();
  canv_Q->cd(3); QTTw->Draw();
  canv_Q->cd(5); QTTx->Draw();
  */

  // draw raw rellum graphs
  TCanvas * canv_R_raw = new TCanvas("canv_R_raw","canv_R_raw",1200,700);
  canv_R_raw->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_R_raw->GetPad(j)->SetGrid(1,1);
  canv_R_raw->cd(1); Rtrig_raw->Draw();
  canv_R_raw->cd(4); Rspin->Draw();
  canv_R_raw->cd(2); RTTe_raw->Draw();
  canv_R_raw->cd(3); RTTw_raw->Draw();
  canv_R_raw->cd(5); RTTx_raw->Draw();

  // draw corrected rellum graphs
  TCanvas * canv_R_corr = new TCanvas("canv_R_corr","canv_R_corr",1200,700);
  canv_R_corr->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_R_corr->GetPad(j)->SetGrid(1,1);
  canv_R_corr->cd(4); Rspin->Draw();
  canv_R_corr->cd(2); RTTe_corr->Draw();
  canv_R_corr->cd(3); RTTw_corr->Draw();
  canv_R_corr->cd(5); RTTx_corr->Draw();

  return;
};
