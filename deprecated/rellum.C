// old version 1 of the relative luminosity calculation 
// -- no corrections for accidentals or multiples here

void rellum(const char * counts_file = "counts.root", const char * trigger="zdc")
{
  // read counts.root file
  TFile * infile = new TFile(counts_file,"READ");
  TTree * s = (TTree*) infile->Get("sca");
  s->Print();

  // set up trigger names for reading
  // (note: variable names refer to bbc, but can be used for zdc
  //  or vpd, since these trigger name strings determine what's 
  //  read from the sca tree)
  Int_t max_trig=0;
  if(!strcmp(trigger,"bbc")) max_trig=300e6;
  else if(!strcmp(trigger,"zdc")) max_trig=50e6;
  else if(!strcmp(trigger,"vpd")) max_trig=150e6;
  else
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
  Double_t bbce,bbcw,bbcx;
  Int_t blue,yell;
  s->SetBranchAddress("i",&index);
  s->SetBranchAddress("runnum",&runnum);
  s->SetBranchAddress("fill",&fill);
  s->SetBranchAddress("bx",&bx);
  s->SetBranchAddress(trigger_e,&bbce);
  s->SetBranchAddress(trigger_w,&bbcw);
  s->SetBranchAddress(trigger_x,&bbcx);
  s->SetBranchAddress("blue",&blue);
  s->SetBranchAddress("yell",&yell);
  
  // determine maximum number of runs
  Int_t i_max_tmp = s->GetMaximum("i");
  const Int_t i_max = i_max_tmp;

  // define distributions
  char bbce_same_hist_n[32]; sprintf(bbce_same_hist_n,"%s_same",trigger_e);
  char bbce_diff_hist_n[32]; sprintf(bbce_diff_hist_n,"%s_diff",trigger_e);
  char bbcw_same_hist_n[32]; sprintf(bbcw_same_hist_n,"%s_same",trigger_w);
  char bbcw_diff_hist_n[32]; sprintf(bbcw_diff_hist_n,"%s_diff",trigger_w);
  char bbcx_same_hist_n[32]; sprintf(bbcx_same_hist_n,"%s_same",trigger_x);
  char bbcx_diff_hist_n[32]; sprintf(bbcx_diff_hist_n,"%s_diff",trigger_x);
  char trig_same_hist_n[32]; sprintf(trig_same_hist_n,"%s_same","trig");
  char trig_diff_hist_n[32]; sprintf(trig_diff_hist_n,"%s_diff","trig");
  char bbce_hist_t[32]; sprintf(bbce_hist_t,"%s G:same P:diff",trigger_e);
  char bbcw_hist_t[32]; sprintf(bbcw_hist_t,"%s G:same P:diff",trigger_w);
  char bbcx_hist_t[32]; sprintf(bbcx_hist_t,"%s G:same P:diff",trigger_x);
  char trig_hist_t[32]; sprintf(trig_hist_t,"%s G:same P:diff","trig");

  TH1F * bbce_same_hist = new TH1F(bbce_same_hist_n,bbce_hist_t,300,0,max_trig);
  TH1F * bbce_diff_hist = new TH1F(bbce_diff_hist_n,bbce_hist_t,300,0,max_trig);
  TH1F * bbcw_same_hist = new TH1F(bbcw_same_hist_n,bbcw_hist_t,300,0,max_trig);
  TH1F * bbcw_diff_hist = new TH1F(bbcw_diff_hist_n,bbcw_hist_t,300,0,max_trig);
  TH1F * bbcx_same_hist = new TH1F(bbcx_same_hist_n,bbcx_hist_t,300,0,max_trig);
  TH1F * bbcx_diff_hist = new TH1F(bbcx_diff_hist_n,bbcx_hist_t,300,0,max_trig);
  TH1F * trig_same_hist = new TH1F(trig_same_hist_n,trig_hist_t,300,0,max_trig*3);
  TH1F * trig_diff_hist = new TH1F(trig_diff_hist_n,trig_hist_t,300,0,max_trig*3);

  // arrays for graphs
  Long_t bbce_same[i_max], bbce_diff[i_max];
  Long_t bbcw_same[i_max], bbcw_diff[i_max];
  Long_t bbcx_same[i_max], bbcx_diff[i_max];
  Long_t trig_same[i_max], trig_diff[i_max];
  Int_t spin_same[i_max], spin_diff[i_max];
  Float_t indexArr[i_max];
  Int_t runnumArr[i_max];
  for(Int_t j=0; j<i_max; j++)
  {
    bbce_same[j] = bbce_diff[j] = 0;
    bbcw_same[j] = bbcw_diff[j] = 0;
    bbcx_same[j] = bbcx_diff[j] = 0;
    trig_same[j] = trig_diff[j] = 0;
    spin_same[j] = spin_diff[j] = 0;
    indexArr[j] = (Float_t)(j+1);
  };

  // loop through tree and fill arrays
  Int_t ix,trig;
  for(Int_t q=0; q<s->GetEntries(); q++)
  {
    s->GetEntry(q);
    ix = index - 1; 
    runnumArr[ix] = runnum;
    if(blue!=0 && yell!=0)
    {
      // TRIG -- the main variable name for relative luminosity computation

      // naive trig set to sum of east & west triggers
      // trig = bbce + bbcw;
      // naive inclusion of coincidence triggers
      trig = bbce + bbcw + bbcx;

      if(blue==yell) 
      {
        bbce_same[ix] += bbce;
        bbcw_same[ix] += bbcw;
        bbcx_same[ix] += bbcx;
        trig_same[ix] += trig;
        spin_same[ix]++;
        bbce_same_hist->Fill(bbce);
        bbcw_same_hist->Fill(bbcw);
        bbcx_same_hist->Fill(bbcx);
        trig_same_hist->Fill(trig);
      }
      else if(blue!=yell)
      {
        bbce_diff[ix] += bbce;
        bbcw_diff[ix] += bbcw;
        bbcx_diff[ix] += bbcx;
        trig_diff[ix] += trig;
        spin_diff[ix]++;
        bbce_diff_hist->Fill(bbce);
        bbcw_diff_hist->Fill(bbcw);
        bbcx_diff_hist->Fill(bbcx);
        trig_diff_hist->Fill(trig);
      };

    };
  };

  // compute relative luminosity
  Float_t bbce_R[i_max];
  Float_t bbcw_R[i_max];
  Float_t bbcx_R[i_max];
  Float_t trig_R[i_max];
  Float_t spin_R[i_max];
  char rdat[32];
  strcpy(rdat,"rdat");
  char rdat_clear[128];
  sprintf(rdat_clear,"touch %s; rm %s; touch %s",rdat,rdat,rdat);
  for(Int_t j=0; j<i_max; j++)
  {

    // check for overfill on trig
    if(trig_same[j]<0 || trig_diff[j]<0) 
      printf("WARNING: run index %d has counts greater than max allowable long int\n",j);

    bbce_R[j] = ((Float_t)bbce_same[j])/((Float_t)bbce_diff[j]);
    bbcw_R[j] = ((Float_t)bbcw_same[j])/((Float_t)bbcw_diff[j]);
    bbcx_R[j] = ((Float_t)bbcx_same[j])/((Float_t)bbcx_diff[j]);
    trig_R[j] = ((Float_t)trig_same[j])/((Float_t)trig_diff[j]);
    spin_R[j] = ((Float_t)spin_same[j])/((Float_t)spin_diff[j]);
    gSystem->RedirectOutput(rdat);
    printf("%d %d %f\n",indexArr[j],runnumArr[j],trig_R[j]);
    gSystem->RedirectOutput(0);
  };

  // plot graphs from arrays
  //
  TGraph * gr_bbce = new TGraph(i_max,indexArr,bbce_R);
  char gr_bbce_t[32]; sprintf(gr_bbce_t,"R_%s :: index",trigger_e);
  gr_bbce->SetTitle(gr_bbce_t);
  gr_bbce->GetXaxis()->SetTitle("run index");
  gr_bbce->GetYaxis()->SetTitle("L_{++}/L_{+-}");
  gr_bbce->SetMarkerStyle(kFullDotMedium);
  gr_bbce->SetMarkerColor(kBlue+2);

  TGraph * gr_bbcw = new TGraph(i_max,indexArr,bbcw_R);
  char gr_bbcw_t[32]; sprintf(gr_bbcw_t,"R_%s :: index",trigger_w);
  gr_bbcw->SetTitle(gr_bbcw_t);
  gr_bbcw->GetXaxis()->SetTitle("run index");
  gr_bbcw->GetYaxis()->SetTitle("L_{++}/L_{+-}");
  gr_bbcw->SetMarkerStyle(kFullDotMedium);
  gr_bbcw->SetMarkerColor(kBlue+2);

  TGraph * gr_bbcx = new TGraph(i_max,indexArr,bbcx_R);
  char gr_bbcx_t[32]; sprintf(gr_bbcx_t,"R_%s :: index",trigger_x);
  gr_bbcx->SetTitle(gr_bbcx_t);
  gr_bbcx->GetXaxis()->SetTitle("run index");
  gr_bbcx->GetYaxis()->SetTitle("L_{++}/L_{+-}");
  gr_bbcx->SetMarkerStyle(kFullDotMedium);
  gr_bbcx->SetMarkerColor(kBlue+2);

  TGraph * gr_trig = new TGraph(i_max,indexArr,trig_R);
  gr_trig->SetTitle("R_trig :: index");
  gr_trig->GetXaxis()->SetTitle("run index");
  gr_trig->GetYaxis()->SetTitle("L_{++}/L_{+-}");
  gr_trig->SetMarkerStyle(kFullDotMedium);
  gr_trig->SetMarkerColor(kBlue+2);

  TGraph * gr_spin = new TGraph(i_max,indexArr,spin_R);
  gr_spin->SetTitle("R_spin :: index");
  gr_spin->GetXaxis()->SetTitle("run index");
  gr_spin->GetYaxis()->SetTitle("L_{++}/L_{+-}");
  gr_spin->SetMarkerStyle(kFullDotMedium);
  gr_spin->SetMarkerColor(kBlue+2);

  gr_bbce->Fit("pol0","","",0,i_max);
  gr_bbcw->Fit("pol0","","",0,i_max);
  gr_bbcx->Fit("pol0","","",0,i_max);
  gr_trig->Fit("pol0","","",0,i_max);
  gr_spin->Fit("pol0","","",0,i_max);

  // set distribution colours
  bbce_same_hist->SetLineColor(kGreen+3);
  bbcw_same_hist->SetLineColor(kGreen+3);
  bbcx_same_hist->SetLineColor(kGreen+3);
  trig_same_hist->SetLineColor(kGreen+3);
  bbce_diff_hist->SetLineColor(kMagenta+2);
  bbcw_diff_hist->SetLineColor(kMagenta+2);
  bbcx_diff_hist->SetLineColor(kMagenta+2);
  trig_diff_hist->SetLineColor(kMagenta+2);

  // draw distributions
  TCanvas * canv_dists = new TCanvas("canv_dists","canv_dists",1200,700);
  canv_dists->Divide(3,2);
  //for(Int_t j=1; j<=6; j++) canv_dists->GetPad(j)->SetLogy();
  canv_dists->cd(1); trig_same_hist->Draw(); trig_diff_hist->Draw("same");
  canv_dists->cd(2); bbce_same_hist->Draw(); bbce_diff_hist->Draw("same");
  canv_dists->cd(3); bbcw_same_hist->Draw(); bbcw_diff_hist->Draw("same");
  canv_dists->cd(4); bbcx_same_hist->Draw(); bbcx_diff_hist->Draw("same");

  // draw ratio graphs (R=L++/L+- :: index)
  TCanvas * canv_R = new TCanvas("canv_R","canv_R",1200,700);
  canv_R->Divide(3,2);
  for(Int_t j=1; j<=6; j++)
  {
    //canv_R->GetPad(j)->SetLogy();
    canv_R->GetPad(j)->SetGrid(1,1);
  };
  canv_R->cd(1); gr_trig->Draw("AP");
  canv_R->cd(4); gr_spin->Draw("AP");
  canv_R->cd(2); gr_bbce->Draw("AP");
  canv_R->cd(3); gr_bbcw->Draw("AP");
  canv_R->cd(5); gr_bbcx->Draw("AP");


  // single windows (pres)
  gStyle->SetOptFit(1);
  new TCanvas(); trig_same_hist->Draw(); trig_diff_hist->Draw("same");
  new TCanvas(); bbce_same_hist->Draw(); bbce_diff_hist->Draw("same");
  new TCanvas(); bbcw_same_hist->Draw(); bbcw_diff_hist->Draw("same");
  new TCanvas(); bbcx_same_hist->Draw(); bbcx_diff_hist->Draw("same");
  new TCanvas(); gr_trig->Draw("AP");
  new TCanvas(); gr_spin->Draw("AP");
  new TCanvas(); gr_bbce->Draw("AP");
  new TCanvas(); gr_bbcw->Draw("AP");
  new TCanvas(); gr_bbcx->Draw("AP");

  TFile * tf = new TFile("out.root","RECREATE");
  gr_bbce->Write("gr_bbce");
  gr_bbcw->Write("gr_bbcw");
  gr_bbcx->Write("gr_bbcx");
  gr_spin->Write("gr_spin");
  gr_trig->Write("gr_trig");
  bbce_same_hist->Write();
  bbce_diff_hist->Write();
  bbcw_same_hist->Write();
  bbcw_diff_hist->Write();
  bbcx_same_hist->Write();
  bbcx_diff_hist->Write();
  trig_same_hist->Write();
  trig_diff_hist->Write();
};
