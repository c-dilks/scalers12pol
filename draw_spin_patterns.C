// draws spin patterns for all the fills
// -- run in background!
// -- execute after accumulate!
// --> produces pattern.pdf, which is graphical representation of spin patterns

void draw_spin_patterns(Bool_t onlyDrawEmptyBuckets=false)
{
  TFile * infile = new TFile("counts.root","READ");
  TTree * tr = (TTree*) infile->Get("sca");

  // count number of fills
  // add fill numbers to fill_arr
  Int_t fill,fill_tmp;
  Int_t fill_arr[1000];
  tr->SetBranchAddress("fill",&fill);
  fill_tmp=0;
  Int_t nfills_tmp=0;
  for(Int_t k=0; k<tr->GetEntries(); k++)
  {
    tr->GetEntry(k);
    if(fill!=fill_tmp)
    {
      fill_arr[nfills_tmp] = fill;
      fill_tmp=fill;
      nfills_tmp++;
    };
  };
  const Int_t NFILLS = nfills_tmp;
  printf("%d\n",NFILLS);

  // define spin pat histograms and project
  TH2D * spinpat_blue[NFILLS];
  TH2D * spinpat_yell[NFILLS];
  char spinpat_blue_n[NFILLS][32];
  char spinpat_yell_n[NFILLS][32];
  char cut[NFILLS][32];
  for(Int_t k=0; k<NFILLS; k++)
  {
    sprintf(spinpat_blue_n[k],"blue %d",fill_arr[k]);
    sprintf(spinpat_yell_n[k],"yell %d",fill_arr[k]);
    if(onlyDrawEmptyBuckets) sprintf(cut[k],"fill==%d && blue*yell==0",fill_arr[k]);
    else sprintf(cut[k],"fill==%d",fill_arr[k]);
    spinpat_blue[k] = new TH2D(spinpat_blue_n[k],spinpat_blue_n[k],3,-1.5,1.5,120,-0.5,119.5);
    spinpat_yell[k] = new TH2D(spinpat_yell_n[k],spinpat_yell_n[k],3,-1.5,1.5,120,-0.5,119.5);
    spinpat_blue[k]->GetYaxis()->SetLabelSize(0.07);
    spinpat_yell[k]->GetYaxis()->SetLabelSize(0.07);
    tr->Project(spinpat_blue_n[k],"bx:blue",cut[k]);
    tr->Project(spinpat_yell_n[k],"bx:yell",cut[k]);
  };

  // recolour spinpatterns
  Int_t bc_pos_blue,bc_pos_yell; // binx=3
  Int_t bc_neg_blue,bc_neg_yell; // binx=1
  Int_t bc_nul_blue,bc_nul_yell; // binx=2
  Int_t binx,biny;
  Int_t blue_state,yell_state;
  TH2D * spinpat_blue_col[NFILLS];
  TH2D * spinpat_yell_col[NFILLS];
  char spinpat_blue_col_n[NFILLS][32];
  char spinpat_yell_col_n[NFILLS][32];
  for(Int_t k=0; k<NFILLS; k++)
  {
    sprintf(spinpat_blue_col_n[k],"blue_%d",fill_arr[k]);
    sprintf(spinpat_yell_col_n[k],"yell_%d",fill_arr[k]);
    spinpat_blue_col[k] = (TH2D*) spinpat_blue[k]->Clone(spinpat_blue_col_n[k]);
    spinpat_yell_col[k] = (TH2D*) spinpat_yell[k]->Clone(spinpat_yell_col_n[k]);
    for(biny=1; biny<=spinpat_blue[k]->GetNbinsY(); biny++)
    {
      bc_neg_blue = spinpat_blue[k]->GetBinContent(1,biny);
      bc_nul_blue = spinpat_blue[k]->GetBinContent(2,biny);
      bc_pos_blue = spinpat_blue[k]->GetBinContent(3,biny);
      bc_neg_yell = spinpat_yell[k]->GetBinContent(1,biny);
      bc_nul_yell = spinpat_yell[k]->GetBinContent(2,biny);
      bc_pos_yell = spinpat_yell[k]->GetBinContent(3,biny);

      if(bc_neg_blue>0) blue_state=-1;
      else if(bc_nul_blue>0) blue_state=0;
      else if(bc_pos_blue>0) blue_state=1;
      if(bc_neg_yell>0) yell_state=-1;
      else if(bc_nul_yell>0) yell_state=0;
      else if(bc_pos_yell>0) yell_state=1;

      if(blue_state*yell_state==0) 
      {
        spinpat_blue_col[k]->SetBinContent(blue_state+2,biny,1);
        spinpat_yell_col[k]->SetBinContent(yell_state+2,biny,1);
      }
      else if(blue_state==yell_state)
      {
        spinpat_blue_col[k]->SetBinContent(blue_state+2,biny,2);
        spinpat_yell_col[k]->SetBinContent(yell_state+2,biny,2);
      }
      else if(blue_state!=yell_state)
      {
        spinpat_blue_col[k]->SetBinContent(blue_state+2,biny,3);
        spinpat_yell_col[k]->SetBinContent(yell_state+2,biny,3);
      }

      // resize for box drawing (bc_neg_blue et al are, if nonzero, 
      // the number of runs for a fill #; here we change that
      // to 10, so that all boxes are the same size.. this makes
      // the drawing easier to read)
      spinpat_blue[k]->SetBinContent(blue_state+2,biny,10);
      spinpat_yell[k]->SetBinContent(yell_state+2,biny,10);
    };
  };
  

  // draw loop
  TCanvas * cc = new TCanvas("cc","cc",1000,2000);
  cc->Divide(4,1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleY(1.05);
  gStyle->SetTitleFontSize(0.12);
  Int_t cd_num;
  Int_t pg_num=0;
  Int_t j;
  for(j=0; j<NFILLS; j++)
  {
    cd_num = 2*(j%2)+1;
    spinpat_blue[j]->SetMarkerSize(3);
    spinpat_yell[j]->SetMarkerSize(3);
    spinpat_blue[j]->SetNdivisions(524,"y");
    spinpat_yell[j]->SetNdivisions(524,"y");
    spinpat_blue_col[j]->SetNdivisions(524,"y");
    spinpat_yell_col[j]->SetNdivisions(524,"y");


    for(Int_t kk=1; kk<=4; kk++) cc->GetPad(kk)->SetGrid(0,1);
    


    cc->cd(cd_num); 
     spinpat_blue_col[j]->Draw("col");
     spinpat_blue[j]->Draw("boxsame");
    cc->cd(cd_num+1); 
     spinpat_yell_col[j]->Draw("col");
     spinpat_yell[j]->Draw("boxsame");
    if(j%2)
    {
      if(pg_num) cc->Print("pattern.pdf","pdf");
      else cc->Print("pattern.pdf(","pdf");
      cc->Clear();
      cc->Divide(4,1);
    };
  };
  printf("%d %d\n",j,fill);

  if(j%2) cc->Print("pattern.pdf)","pdf");
  else { cc->Clear(); cc->Print("pattern.pdf)","pdf"); };
};
