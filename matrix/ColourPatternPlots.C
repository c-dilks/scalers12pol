// produces coloured "vs. run index" plots, where the colours denote the different spin patterns;
//
// this is used to see if there are any correlations with spin pattern and bunch fit results

void ColourPatternPlots(const char * numerDet="zdce", const char * denomDet="vpdx")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  // define pattern numbers
  //  -- there are 8 patterns, but NPAT = 8 + 1, since the zeroth index will
  //     represent the "full" plots (where "full" = all patterns)
  const Int_t NPAT=9;
  Int_t pat[NPAT] = {0,13,14,23,24,31,32,41,42};
  char pat_str[NPAT][256];
  char beampat[5][20];
  strcpy(beampat[1],"+ - + - - + - +");
  strcpy(beampat[2],"- + - + + - + -");
  strcpy(beampat[3],"+ + - - + + - -");
  strcpy(beampat[4],"- - + + - - + +");
  // old way
  /*
  sprintf(pat_str[1],"B[%s]  Y[%s]  (%d)",beampat[1],beampat[3],pat[1]); 
  sprintf(pat_str[2],"B[%s]  Y[%s]  (%d)",beampat[1],beampat[4],pat[2]); 
  sprintf(pat_str[3],"B[%s]  Y[%s]  (%d)",beampat[2],beampat[3],pat[3]); 
  sprintf(pat_str[4],"B[%s]  Y[%s]  (%d)",beampat[2],beampat[4],pat[4]); 
  sprintf(pat_str[5],"B[%s]  Y[%s]  (%d)",beampat[3],beampat[1],pat[5]); 
  sprintf(pat_str[6],"B[%s]  Y[%s]  (%d)",beampat[3],beampat[2],pat[6]); 
  sprintf(pat_str[7],"B[%s]  Y[%s]  (%d)",beampat[4],beampat[1],pat[7]); 
  sprintf(pat_str[8],"B[%s]  Y[%s]  (%d)",beampat[4],beampat[2],pat[8]); 
  */
  // new way (just has more details)
  sprintf(pat_str[1],"B_{%d}[%s]  Y_{%d}[%s]  ( S O O S O S S O )",pat[1]/10,beampat[1],pat[1]%10,beampat[3]); 
  sprintf(pat_str[2],"B_{%d}[%s]  Y_{%d}[%s]  ( O S S O S O O S )",pat[2]/10,beampat[1],pat[2]%10,beampat[4]); 
  sprintf(pat_str[3],"B_{%d}[%s]  Y_{%d}[%s]  ( O S S O S O O S )",pat[3]/10,beampat[2],pat[3]%10,beampat[3]); 
  sprintf(pat_str[4],"B_{%d}[%s]  Y_{%d}[%s]  ( S O O S O S S O )",pat[4]/10,beampat[2],pat[4]%10,beampat[4]); 
  sprintf(pat_str[5],"B_{%d}[%s]  Y_{%d}[%s]  ( S O O S O S S O )",pat[5]/10,beampat[3],pat[5]%10,beampat[1]); 
  sprintf(pat_str[6],"B_{%d}[%s]  Y_{%d}[%s]  ( O S S O S O O S )",pat[6]/10,beampat[3],pat[6]%10,beampat[2]); 
  sprintf(pat_str[7],"B_{%d}[%s]  Y_{%d}[%s]  ( O S S O S O O S )",pat[7]/10,beampat[4],pat[7]%10,beampat[1]); 
  sprintf(pat_str[8],"B_{%d}[%s]  Y_{%d}[%s]  ( S O O S O S S O )",pat[8]/10,beampat[4],pat[8]%10,beampat[2]); 

  const Int_t NASY=10;

  // open root files
  char infile_n[NPAT][128];
  TFile * infile[NPAT];
  for(Int_t pp=0; pp<NPAT; pp++)
  {
    if(pp==0) sprintf(infile_n[pp],"fit_result.%s.%s.root",numerDet,denomDet);
    else sprintf(infile_n[pp],"pats/fit_result.%s.%s.pat%d.root",numerDet,denomDet,pat[pp]);
    infile[pp] = new TFile(infile_n[pp],"READ");
  };


  // load dists
  TH1D * sigma_id_dist[NPAT][NASY];
  TH1D * sigma_r_dist[NPAT][NASY];
  TH1D * sigma_h_dist[NPAT][NASY];
  TH1D * sigma_hr_dist[NPAT][NASY];
  TH1D * cons_dist[NPAT][NASY];
  TH1D * epsi_dist[NPAT][NASY];
  TH1D * asym_dist[NPAT][4];
  char sigma_id_dist_n[NPAT][NASY][64];
  char sigma_r_dist_n[NPAT][NASY][64];
  char sigma_h_dist_n[NPAT][NASY][64];
  char sigma_hr_dist_n[NPAT][NASY][64];
  char cons_dist_n[NPAT][NASY][64];
  char epsi_dist_n[NPAT][NASY][64];
  char asym_dist_n[NPAT][4][64];
  for(Int_t pp=0; pp<NPAT; pp++)
  {
    for(Int_t aa=1; aa<NASY; aa++)
    {
      sprintf(sigma_id_dist_n[pp][aa],"/sigma_id/sigma_id_a%d",aa);
      sprintf(sigma_r_dist_n[pp][aa],"/sigma_r/sigma_r_a%d",aa);
      sprintf(sigma_h_dist_n[pp][aa],"/sigma_h/sigma_h_a%d",aa);
      sprintf(sigma_hr_dist_n[pp][aa],"/sigma_hr/sigma_hr_a%d",aa);
      sprintf(cons_dist_n[pp][aa],"/constant/cons_a%d_v_run",aa);
      sprintf(epsi_dist_n[pp][aa],"/epsilon/epsi_a%d_v_run",aa);
      if(aa<4) sprintf(asym_dist_n[pp][aa],"/asymmetry/asym_a%d_v_run",aa);
      sigma_id_dist[pp][aa] = (TH1D*) infile[pp]->Get(sigma_id_dist_n[pp][aa]);
      sigma_r_dist[pp][aa] = (TH1D*) infile[pp]->Get(sigma_r_dist_n[pp][aa]);
      sigma_h_dist[pp][aa] = (TH1D*) infile[pp]->Get(sigma_h_dist_n[pp][aa]);
      sigma_hr_dist[pp][aa] = (TH1D*) infile[pp]->Get(sigma_hr_dist_n[pp][aa]);
      cons_dist[pp][aa] = (TH1D*) infile[pp]->Get(cons_dist_n[pp][aa]);
      epsi_dist[pp][aa] = (TH1D*) infile[pp]->Get(epsi_dist_n[pp][aa]);
      if(aa<4) asym_dist[pp][aa] = (TH1D*) infile[pp]->Get(asym_dist_n[pp][aa]);
      if(aa<4) printf("asym_dist[%d][%d] @ %p\n",pp,aa,(void*)asym_dist[pp][aa]);

      sigma_id_dist[pp][aa]->GetXaxis()->SetLabelSize(0.08);
      sigma_r_dist[pp][aa]->GetXaxis()->SetLabelSize(0.08);
      sigma_h_dist[pp][aa]->GetXaxis()->SetLabelSize(0.08);
      sigma_hr_dist[pp][aa]->GetXaxis()->SetLabelSize(0.08);
      cons_dist[pp][aa]->GetXaxis()->SetLabelSize(0.08);
      epsi_dist[pp][aa]->GetXaxis()->SetLabelSize(0.08);
      if(aa<4) asym_dist[pp][aa]->GetXaxis()->SetLabelSize(0.08);

      sigma_id_dist[pp][aa]->GetYaxis()->SetLabelSize(0.08);
      sigma_r_dist[pp][aa]->GetYaxis()->SetLabelSize(0.08);
      sigma_h_dist[pp][aa]->GetYaxis()->SetLabelSize(0.08);
      sigma_hr_dist[pp][aa]->GetYaxis()->SetLabelSize(0.08);
      cons_dist[pp][aa]->GetYaxis()->SetLabelSize(0.08);
      epsi_dist[pp][aa]->GetYaxis()->SetLabelSize(0.08);
      if(aa<4) asym_dist[pp][aa]->GetYaxis()->SetLabelSize(0.08);

      printf("sigma_id_dist[%d][%d] @ %p\n",pp,aa,(void*)sigma_id_dist[pp][aa]);
      printf("sigma_r_dist[%d][%d] @ %p\n",pp,aa,(void*)sigma_r_dist[pp][aa]);
      printf("sigma_h_dist[%d][%d] @ %p\n",pp,aa,(void*)sigma_h_dist[pp][aa]);
      printf("sigma_hr_dist[%d][%d] @ %p\n",pp,aa,(void*)sigma_hr_dist[pp][aa]);
      printf("cons_dist[%d][%d] @ %p\n",pp,aa,(void*)cons_dist[pp][aa]);
      printf("epsi_dist[%d][%d] @ %p\n",pp,aa,(void*)epsi_dist[pp][aa]);
    };
  };


  // set plot colours
  Color_t colours[NPAT] = {kBlack,kOrange,kRed,kMagenta,kBlue,kCyan+1,kGreen+1,kYellow+2,kViolet-6}; // individual colours
  //Color_t colours[NPAT] = {kBlack,kRed,kMagenta,kBlue,kGreen+1,kRed,kBlue,kMagenta,kGreen+1}; // partner colours
  for(Int_t pp=0; pp<NPAT; pp++)
  {
    for(Int_t aa=1; aa<NASY; aa++)
    {
      sigma_id_dist[pp][aa]->SetLineColor(colours[pp]);
      sigma_r_dist[pp][aa]->SetLineColor(colours[pp]);
      sigma_h_dist[pp][aa]->SetLineColor(colours[pp]);
      sigma_hr_dist[pp][aa]->SetLineColor(colours[pp]);
      cons_dist[pp][aa]->SetLineColor(colours[pp]);
      epsi_dist[pp][aa]->SetLineColor(colours[pp]);
      if(aa<4) asym_dist[pp][aa]->SetLineColor(colours[pp]);
      sigma_id_dist[pp][aa]->SetLineWidth(2);
      sigma_r_dist[pp][aa]->SetLineWidth(2);
      sigma_h_dist[pp][aa]->SetLineWidth(2);
      sigma_hr_dist[pp][aa]->SetLineWidth(2);
      cons_dist[pp][aa]->SetLineWidth(2);
      epsi_dist[pp][aa]->SetLineWidth(2);
      if(aa<4) asym_dist[pp][aa]->SetLineWidth(2);
    };
  };

  
  // build TCanvases; three asymmetry numbers per canvas: 1-3, 4-6, 7-9
  Int_t sf = 1;
  TCanvas * sigma_id_canv[3];
  TCanvas * sigma_r_canv[3];
  TCanvas * sigma_h_canv[3];
  TCanvas * sigma_hr_canv[3];
  TCanvas * cons_canv[3];
  TCanvas * epsi_canv[3];
  TCanvas * asym_canv;
  char sigma_id_canv_n[3][32];
  char sigma_r_canv_n[3][32];
  char sigma_h_canv_n[3][32];
  char sigma_hr_canv_n[3][32];
  char cons_canv_n[3][32];
  char epsi_canv_n[3][32];
  char asym_canv_n[32];
  char asy_str[3][8];
  strcpy(asy_str[0],"1-3");
  strcpy(asy_str[1],"4-6");
  strcpy(asy_str[2],"7-9");
  for(Int_t cc=0; cc<3; cc++)
  {
    printf("cc=%d\n",cc);
    sprintf(sigma_id_canv_n[cc],"sigma_id_canv_a%s",asy_str[cc]);
    sprintf(sigma_r_canv_n[cc],"sigma_r_canv_a%s",asy_str[cc]);
    sprintf(sigma_h_canv_n[cc],"sigma_h_canv_a%s",asy_str[cc]);
    sprintf(sigma_hr_canv_n[cc],"sigma_hr_canv_a%s",asy_str[cc]);
    sprintf(cons_canv_n[cc],"cons_canv_a%s",asy_str[cc]);
    sprintf(epsi_canv_n[cc],"epsi_canv_a%s",asy_str[cc]);
    sigma_id_canv[cc] = new TCanvas(sigma_id_canv_n[cc],sigma_id_canv_n[cc],1100*sf,700*sf);
    sigma_h_canv[cc] = new TCanvas(sigma_h_canv_n[cc],sigma_h_canv_n[cc],1100*sf,700*sf);
    sigma_r_canv[cc] = new TCanvas(sigma_r_canv_n[cc],sigma_r_canv_n[cc],1100*sf,700*sf);
    sigma_hr_canv[cc] = new TCanvas(sigma_hr_canv_n[cc],sigma_hr_canv_n[cc],1100*sf,700*sf);
    cons_canv[cc] = new TCanvas(cons_canv_n[cc],cons_canv_n[cc],1100*sf,700*sf);
    epsi_canv[cc] = new TCanvas(epsi_canv_n[cc],epsi_canv_n[cc],1100*sf,700*sf);
    char drawtype[8];

    sigma_id_canv[cc]->Divide(1,3);
    sigma_r_canv[cc]->Divide(1,3);
    sigma_h_canv[cc]->Divide(1,3);
    sigma_hr_canv[cc]->Divide(1,3);
    cons_canv[cc]->Divide(1,3);
    epsi_canv[cc]->Divide(1,3);
    
    for(Int_t dd=1; dd<=3; dd++)
    {
      printf(" dd=%d\n",dd);
      sigma_id_canv[cc]->GetPad(dd)->SetGrid(1,1);
      sigma_id_canv[cc]->cd(dd);
      for(Int_t pp=0; pp<NPAT; pp++)
      {
        if(pp==0) strcpy(drawtype,"E");
        else strcpy(drawtype,"ESAME");
        sigma_id_dist[pp][cc*3+dd]->Draw(drawtype);
      };

      sigma_r_canv[cc]->GetPad(dd)->SetGrid(1,1);
      sigma_r_canv[cc]->cd(dd);
      for(Int_t pp=0; pp<NPAT; pp++)
      {
        if(pp==0) strcpy(drawtype,"E");
        else strcpy(drawtype,"ESAME");
        sigma_r_dist[pp][cc*3+dd]->Draw(drawtype);
      };

      sigma_h_canv[cc]->GetPad(dd)->SetGrid(1,1);
      sigma_h_canv[cc]->cd(dd);
      for(Int_t pp=0; pp<NPAT; pp++)
      {
        if(pp==0) strcpy(drawtype,"E");
        else strcpy(drawtype,"ESAME");
        sigma_h_dist[pp][cc*3+dd]->Draw(drawtype);
      };

      sigma_hr_canv[cc]->GetPad(dd)->SetGrid(1,1);
      sigma_hr_canv[cc]->cd(dd);
      for(Int_t pp=0; pp<NPAT; pp++)
      {
        if(pp==0) strcpy(drawtype,"E");
        else strcpy(drawtype,"ESAME");
        sigma_hr_dist[pp][cc*3+dd]->Draw(drawtype);
      };

      cons_canv[cc]->GetPad(dd)->SetGrid(1,1);
      cons_canv[cc]->cd(dd);
      for(Int_t pp=0; pp<NPAT; pp++)
      {
        if(pp==0) strcpy(drawtype,"E");
        else strcpy(drawtype,"ESAME");
        cons_dist[pp][cc*3+dd]->Draw(drawtype);
      };

      epsi_canv[cc]->GetPad(dd)->SetGrid(1,1);
      epsi_canv[cc]->cd(dd);
      for(Int_t pp=0; pp<NPAT; pp++)
      {
        if(pp==0) strcpy(drawtype,"E");
        else strcpy(drawtype,"ESAME");
        epsi_dist[pp][cc*3+dd]->Draw(drawtype);
      };
    };
  };
  // asymmetry canvas
  sprintf(asym_canv_n,"asym_canv_a1-3");
  asym_canv = new TCanvas(asym_canv_n,asym_canv_n,1100*sf,700*sf);
  asym_canv->Divide(1,3);
  for(Int_t dd=1; dd<=3; dd++) 
  {
    asym_canv->GetPad(dd)->SetGrid(1,1);
    asym_canv->cd(dd);
    for(Int_t pp=0; pp<NPAT; pp++)
    {
      if(pp==0) strcpy(drawtype,"E");
      else strcpy(drawtype,"ESAME");
      asym_dist[pp][dd]->Draw(drawtype);
    };
  };
      


  
  // make legend canvas
  TCanvas * leg_canv = new TCanvas("leg_canv","leg_canv",1000,500);
  TLegend * leg = new TLegend(0,0,1,1);
  for(Int_t pp=1; pp<NPAT; pp++) leg->AddEntry(epsi_dist[pp][1],pat_str[pp],"l");
  leg_canv->cd();
  leg->Draw();
    



  // write canvases to root file
  TFile * outfile;
  char outfile_n[64];
  sprintf(outfile_n,"colour.%s.%s.root",numerDet,denomDet);
  outfile = new TFile(outfile_n,"RECREATE");
  for(Int_t cc=0; cc<3; cc++) sigma_id_canv[cc]->Write(sigma_id_canv_n[cc]);
  for(Int_t cc=0; cc<3; cc++) sigma_r_canv[cc]->Write(sigma_r_canv_n[cc]);
  for(Int_t cc=0; cc<3; cc++) sigma_h_canv[cc]->Write(sigma_h_canv_n[cc]);
  for(Int_t cc=0; cc<3; cc++) sigma_hr_canv[cc]->Write(sigma_hr_canv_n[cc]);
  for(Int_t cc=0; cc<3; cc++) cons_canv[cc]->Write(cons_canv_n[cc]);
  for(Int_t cc=0; cc<3; cc++) epsi_canv[cc]->Write(epsi_canv_n[cc]);
  asym_canv->Write(asym_canv_n);
  leg_canv->Write("legend");
  printf("%s written\n",outfile_n);
};
