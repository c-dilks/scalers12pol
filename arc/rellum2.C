void rellum2(const char * var="i", const char * trigger="zdc")
{
  // read counts.root file
  TFile * infile = new TFile("counts.root","READ");
  TTree * s = (TTree*) infile->Get("sca");
  s->Print();

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
  

  // define distributions
  char TTe_same_hist_n[32]; sprintf(TTe_same_hist_n,"%s_same",trigger_e);
  char TTe_diff_hist_n[32]; sprintf(TTe_diff_hist_n,"%s_diff",trigger_e);
  char TTw_same_hist_n[32]; sprintf(TTw_same_hist_n,"%s_same",trigger_w);
  char TTw_diff_hist_n[32]; sprintf(TTw_diff_hist_n,"%s_diff",trigger_w);
  char TTx_same_hist_n[32]; sprintf(TTx_same_hist_n,"%s_same",trigger_x);
  char TTx_diff_hist_n[32]; sprintf(TTx_diff_hist_n,"%s_diff",trigger_x);
  char trig_same_hist_n[32]; sprintf(trig_same_hist_n,"%s_same","trig");
  char trig_diff_hist_n[32]; sprintf(trig_diff_hist_n,"%s_diff","trig");
  char spin_same_hist_n[32]; sprintf(spin_same_hist_n,"%s_same","spin");
  char spin_diff_hist_n[32]; sprintf(spin_diff_hist_n,"%s_diff","spin");
  char TTe_hist_t[32]; sprintf(TTe_hist_t,"%s weighted by %s rate",var,trigger_e); 
  char TTw_hist_t[32]; sprintf(TTw_hist_t,"%s weighted by %s rate",var,trigger_w);
  char TTx_hist_t[32]; sprintf(TTx_hist_t,"%s weighted by %s rate",var,trigger_x);
  char trig_hist_t[32]; sprintf(trig_hist_t,"%s weighted by trig rate",var);
  char spin_hist_t[32]; sprintf(spin_hist_t,"%s weighted by N_bx",var);

  TH1F * TTe_same_hist = new TH1F(TTe_same_hist_n,TTe_hist_t,var_bins,var_l,var_h);
  TH1F * TTe_diff_hist = new TH1F(TTe_diff_hist_n,TTe_hist_t,var_bins,var_l,var_h);
  TH1F * TTw_same_hist = new TH1F(TTw_same_hist_n,TTw_hist_t,var_bins,var_l,var_h);
  TH1F * TTw_diff_hist = new TH1F(TTw_diff_hist_n,TTw_hist_t,var_bins,var_l,var_h);
  TH1F * TTx_same_hist = new TH1F(TTx_same_hist_n,TTx_hist_t,var_bins,var_l,var_h);
  TH1F * TTx_diff_hist = new TH1F(TTx_diff_hist_n,TTx_hist_t,var_bins,var_l,var_h);
  TH1F * trig_same_hist = new TH1F(trig_same_hist_n,trig_hist_t,var_bins,var_l,var_h);
  TH1F * trig_diff_hist = new TH1F(trig_diff_hist_n,trig_hist_t,var_bins,var_l,var_h);
  TH1F * spin_same_hist = new TH1F(spin_same_hist_n,spin_hist_t,var_bins,var_l,var_h);
  TH1F * spin_diff_hist = new TH1F(spin_diff_hist_n,spin_hist_t,var_bins,var_l,var_h);
  

  // --- projections --------------------------- ****

  char TTe_same_hist_cut[256];
  char TTe_diff_hist_cut[256];
  char TTw_same_hist_cut[256];
  char TTw_diff_hist_cut[256];
  char TTx_same_hist_cut[256];
  char TTx_diff_hist_cut[256];
  char trig_same_hist_cut[256];
  char trig_diff_hist_cut[256];
  char spin_same_hist_cut[256];
  char spin_diff_hist_cut[256];

  sprintf(TTe_same_hist_cut,"%s/t*(blue*yell!=0 && blue==yell)",trigger_e);
  sprintf(TTe_diff_hist_cut,"%s/t*(blue*yell!=0 && blue!=yell)",trigger_e);
  sprintf(TTw_same_hist_cut,"%s/t*(blue*yell!=0 && blue==yell)",trigger_w);
  sprintf(TTw_diff_hist_cut,"%s/t*(blue*yell!=0 && blue!=yell)",trigger_w);
  sprintf(TTx_same_hist_cut,"%s/t*(blue*yell!=0 && blue==yell)",trigger_x);
  sprintf(TTx_diff_hist_cut,"%s/t*(blue*yell!=0 && blue!=yell)",trigger_x);
  sprintf(trig_same_hist_cut,"(%s+%s+%s)/t*(blue*yell!=0 && blue==yell)",
                              trigger_e,trigger_w,trigger_x);
  sprintf(trig_diff_hist_cut,"(%s+%s+%s)/t*(blue*yell!=0 && blue!=yell)",
                              trigger_e,trigger_w,trigger_x);
  sprintf(spin_same_hist_cut,"1*(blue*yell!=0 && blue==yell)");
  sprintf(spin_diff_hist_cut,"1*(blue*yell!=0 && blue!=yell)");

  s->Project(TTe_same_hist_n,var,TTe_same_hist_cut);
  s->Project(TTe_diff_hist_n,var,TTe_diff_hist_cut);
  s->Project(TTw_same_hist_n,var,TTw_same_hist_cut);
  s->Project(TTw_diff_hist_n,var,TTw_diff_hist_cut);
  s->Project(TTx_same_hist_n,var,TTx_same_hist_cut);
  s->Project(TTx_diff_hist_n,var,TTx_diff_hist_cut);
  s->Project(trig_same_hist_n,var,trig_same_hist_cut);
  s->Project(trig_diff_hist_n,var,trig_diff_hist_cut);
  s->Project(spin_same_hist_n,var,spin_same_hist_cut);
  s->Project(spin_diff_hist_n,var,spin_diff_hist_cut);

  char RTTe_t[32]; sprintf(RTTe_t,"R_{%s} :: %s",trigger_e,var);
  char RTTw_t[32]; sprintf(RTTw_t,"R_{%s} :: %s",trigger_w,var);
  char RTTx_t[32]; sprintf(RTTx_t,"R_{%s} :: %s",trigger_x,var);
  char Rtrig_t[32]; sprintf(Rtrig_t,"R_{trig} :: %s",var);
  char Rspin_t[32]; sprintf(Rspin_t,"R_{bx} :: %s",var);
  TH1F * RTTe = new TH1F("RTTe",RTTe_t,var_bins,var_l,var_h);
  TH1F * RTTw = new TH1F("RTTw",RTTw_t,var_bins,var_l,var_h);
  TH1F * RTTx = new TH1F("RTTx",RTTx_t,var_bins,var_l,var_h);
  TH1F * Rtrig = new TH1F("Rtrig",Rtrig_t,var_bins,var_l,var_h);
  TH1F * Rspin = new TH1F("Rspin",Rspin_t,var_bins,var_l,var_h);

  RTTe->Divide(TTe_same_hist,TTe_diff_hist,1.0,1.0);
  RTTw->Divide(TTw_same_hist,TTw_diff_hist,1.0,1.0);
  RTTx->Divide(TTx_same_hist,TTx_diff_hist,1.0,1.0);
  Rtrig->Divide(trig_same_hist,trig_diff_hist,1.0,1.0);
  Rspin->Divide(spin_same_hist,spin_diff_hist,1.0,1.0);

  // Q variable Q:=|R-1|
  char QTTe_t[32]; sprintf(QTTe_t,"Q_{%s} :: %s",trigger_e,var);
  char QTTw_t[32]; sprintf(QTTw_t,"Q_{%s} :: %s",trigger_w,var);
  char QTTx_t[32]; sprintf(QTTx_t,"Q_{%s} :: %s",trigger_x,var);
  char Qtrig_t[32]; sprintf(Qtrig_t,"Q_{trig} :: %s",var);
  char Qspin_t[32]; sprintf(Qspin_t,"Q_{spin} :: %s",var);
  TH1F * QTTe = new TH1F("QTTe",QTTe_t,var_bins,var_l,var_h);
  TH1F * QTTw = new TH1F("QTTw",QTTw_t,var_bins,var_l,var_h);
  TH1F * QTTx = new TH1F("QTTx",QTTx_t,var_bins,var_l,var_h);
  TH1F * Qtrig = new TH1F("Qtrig",Qtrig_t,var_bins,var_l,var_h);
  TH1F * Qspin = new TH1F("Qspin",Qspin_t,var_bins,var_l,var_h);
  for(Int_t r=1; r<=RTTe->GetNbinsX(); r++)
  {
    QTTe->SetBinContent(r,fabs(RTTe->GetBinContent(r)-1));
    QTTw->SetBinContent(r,fabs(RTTw->GetBinContent(r)-1));
    QTTx->SetBinContent(r,fabs(RTTx->GetBinContent(r)-1));
    Qtrig->SetBinContent(r,fabs(Rtrig->GetBinContent(r)-1));
    Qspin->SetBinContent(r,fabs(Rspin->GetBinContent(r)-1));
  };

  // Average rates (same spin)
  char ASTTe_t[32]; sprintf(ASTTe_t,"AS_{%s} :: %s",trigger_e,var);
  char ASTTw_t[32]; sprintf(ASTTw_t,"AS_{%s} :: %s",trigger_w,var);
  char ASTTx_t[32]; sprintf(ASTTx_t,"AS_{%s} :: %s",trigger_x,var);
  char AStrig_t[32]; sprintf(AStrig_t,"AS_{trig} :: %s",var);
  TH1F * ASTTe = new TH1F("ASTTe",ASTTe_t,var_bins,var_l,var_h);
  TH1F * ASTTw = new TH1F("ASTTw",ASTTw_t,var_bins,var_l,var_h);
  TH1F * ASTTx = new TH1F("ASTTx",ASTTx_t,var_bins,var_l,var_h);
  TH1F * AStrig = new TH1F("AStrig",AStrig_t,var_bins,var_l,var_h);
  ASTTe->Divide(TTe_same_hist,spin_same_hist,1.0,1.0);
  ASTTw->Divide(TTw_same_hist,spin_same_hist,1.0,1.0);
  ASTTx->Divide(TTx_same_hist,spin_same_hist,1.0,1.0);
  AStrig->Divide(trig_same_hist,spin_same_hist,1.0,1.0);


  // Average rates (diff spin)
  char ADTTe_t[32]; sprintf(ADTTe_t,"AD_{%s} :: %s",trigger_e,var);
  char ADTTw_t[32]; sprintf(ADTTw_t,"AD_{%s} :: %s",trigger_w,var);
  char ADTTx_t[32]; sprintf(ADTTx_t,"AD_{%s} :: %s",trigger_x,var);
  char ADtrig_t[32]; sprintf(ADtrig_t,"AD_{trig} :: %s",var);
  TH1F * ADTTe = new TH1F("ADTTe",ADTTe_t,var_bins,var_l,var_h);
  TH1F * ADTTw = new TH1F("ADTTw",ADTTw_t,var_bins,var_l,var_h);
  TH1F * ADTTx = new TH1F("ADTTx",ADTTx_t,var_bins,var_l,var_h);
  TH1F * ADtrig = new TH1F("ADtrig",ADtrig_t,var_bins,var_l,var_h);
  ADTTe->Divide(TTe_diff_hist,spin_diff_hist,1.0,1.0);
  ADTTw->Divide(TTw_diff_hist,spin_diff_hist,1.0,1.0);
  ADTTx->Divide(TTx_diff_hist,spin_diff_hist,1.0,1.0);
  ADtrig->Divide(trig_diff_hist,spin_diff_hist,1.0,1.0);

  // sum of S's for different and same spins (SS)
  char SSTTe_t[32]; sprintf(SSTTe_t,"SS_{%s} :: %s",trigger_e,var);
  char SSTTw_t[32]; sprintf(SSTTw_t,"SS_{%s} :: %s",trigger_w,var);
  char SSTTx_t[32]; sprintf(SSTTx_t,"SS_{%s} :: %s",trigger_x,var);
  char SStrig_t[32]; sprintf(SStrig_t,"SS_{trig} :: %s",var);
  char SSspin_t[32]; sprintf(SSspin_t,"SS_{bx} :: %s",var);
  TH1F * SSTTe = new TH1F("SSTTe",SSTTe_t,var_bins,var_l,var_h);
  TH1F * SSTTw = new TH1F("SSTTw",SSTTw_t,var_bins,var_l,var_h);
  TH1F * SSTTx = new TH1F("SSTTx",SSTTx_t,var_bins,var_l,var_h);
  TH1F * SStrig = new TH1F("SStrig",SStrig_t,var_bins,var_l,var_h);
  TH1F * SSspin = new TH1F("SSspin",SSspin_t,var_bins,var_l,var_h);
  SSTTe->Add(TTe_same_hist,TTe_diff_hist);
  SSTTw->Add(TTw_same_hist,TTw_diff_hist);
  SSTTx->Add(TTx_same_hist,TTx_diff_hist);
  SStrig->Add(trig_same_hist,trig_diff_hist);
  SSspin->Add(spin_same_hist,spin_diff_hist);

  // Average rates (both spin)
  char ABTTe_t[32]; sprintf(ABTTe_t,"AB_{%s} :: %s",trigger_e,var);
  char ABTTw_t[32]; sprintf(ABTTw_t,"AB_{%s} :: %s",trigger_w,var);
  char ABTTx_t[32]; sprintf(ABTTx_t,"AB_{%s} :: %s",trigger_x,var);
  char ABtrig_t[32]; sprintf(ABtrig_t,"AB_{trig} :: %s",var);
  TH1F * ABTTe = new TH1F("ABTTe",ABTTe_t,var_bins,var_l,var_h);
  TH1F * ABTTw = new TH1F("ABTTw",ABTTw_t,var_bins,var_l,var_h);
  TH1F * ABTTx = new TH1F("ABTTx",ABTTx_t,var_bins,var_l,var_h);
  TH1F * ABtrig = new TH1F("ABtrig",ABtrig_t,var_bins,var_l,var_h);
  ABTTe->Divide(SSTTe,SSspin,1.0,1.0);
  ABTTw->Divide(SSTTw,SSspin,1.0,1.0);
  ABTTx->Divide(SSTTx,SSspin,1.0,1.0);
  ABtrig->Divide(SStrig,SSspin,1.0,1.0);

  /*
  // draw var weighted by total rate (S++ :: var)
  TCanvas * canv_same = new TCanvas("canv_same","canv_same",1200,700);
  canv_same->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_same->GetPad(j)->SetGrid(1,1);
  canv_same->cd(1); trig_same_hist->Draw();
  canv_same->cd(4); spin_same_hist->Draw();
  canv_same->cd(2); TTe_same_hist->Draw();
  canv_same->cd(3); TTw_same_hist->Draw();
  canv_same->cd(5); TTx_same_hist->Draw();
  
  // draw var weighted by total rate (S+- :: var)
  TCanvas * canv_diff = new TCanvas("canv_diff","canv_diff",1200,700);
  canv_diff->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_diff->GetPad(j)->SetGrid(1,1);
  canv_diff->cd(1); trig_diff_hist->Draw();
  canv_diff->cd(4); spin_diff_hist->Draw();
  canv_diff->cd(2); TTe_diff_hist->Draw();
  canv_diff->cd(3); TTw_diff_hist->Draw();
  canv_diff->cd(5); TTx_diff_hist->Draw();

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
  TCanvas * canv_AB = new TCanvas("canv_AB","canv_AB",1200,700);
  canv_AB->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_AB->GetPad(j)->SetGrid(1,1);
  canv_AB->cd(1); ABtrig->Draw();
  canv_AB->cd(2); ABTTe->Draw();
  canv_AB->cd(3); ABTTw->Draw();
  canv_AB->cd(5); ABTTx->Draw();
  
  // draw normaliser (total bx's binned in independent variable "var")
  TCanvas * canv_N = new TCanvas("canv_N","canv_N",600,400);
  canv_N->SetGrid(1,1);
  spin_same_hist->Draw();


  // draw ratio graphs (R=S++/S+- :: var)
  TCanvas * canv_R = new TCanvas("canv_R","canv_R",1200,700);
  canv_R->Divide(3,2);
  for(Int_t j=1; j<=6; j++) canv_R->GetPad(j)->SetGrid(1,1);
  canv_R->cd(1); Rtrig->Draw();
  canv_R->cd(4); Rspin->Draw();
  canv_R->cd(2); RTTe->Draw();
  canv_R->cd(3); RTTw->Draw();
  canv_R->cd(5); RTTx->Draw();

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
  return;
};
