// draws ver* and hor* vs. bXing number for checking for bad bXings / abort gap alignment
// -- if byRun = true, it does it by run number instead of by fill number

void ZDCpolBxCheck(Int_t board=10, Bool_t byRun=false)
{
  TString fname = Form("counts_bd%d.root",board);
  TFile * ff = new TFile(fname.Data(),"READ");
  TTree * tr = (TTree*) ff->Get("sca");
  TString outfile_n = Form("bxs_bd%d.root",board);
  TFile * outfile = new TFile(outfile_n.Data(),"RECREATE");

  TH1D * vh[8];
  TH1D * hh[8];
  TString vhn[8];
  TString hhn[8];
  int x;
  for(x=0;x<8;x++) {
    vhn[x] = Form("vert%d",x);
    hhn[x] = Form("hori%d",x);
    vh[x] = new TH1D(vhn[x].Data(),vhn[x].Data(),120,0,120);
    hh[x] = new TH1D(hhn[x].Data(),hhn[x].Data(),120,0,120);
    vh[x]->SetLineColor(kBlack); vh[x]->SetFillColor(kBlack);
    hh[x]->SetLineColor(kBlack); hh[x]->SetFillColor(kBlack);
  };


  Int_t IMAX_tmp;
  TString invar = byRun? "i":"fi";
  IMAX_tmp = tr->GetMaximum(invar.Data()) + 1;
  const Int_t IMAX = IMAX_tmp;


  Int_t invar_arr[IMAX];
  Int_t index,number;
  tr->SetBranchAddress(invar,&index);
  if(byRun) tr->SetBranchAddress("runnum",&number);
  else tr->SetBranchAddress("fill",&number);
  for(int xx=0; xx<tr->GetEntries(); xx++) {
    tr->GetEntry(xx);
    invar_arr[index] = number;
  };
  for(int xx=0; xx<IMAX; xx++) {
    printf("%d %d\n",xx,invar_arr[xx]);
  };

  TString vcut,hcut;

  TCanvas * cc = new TCanvas("cc","cc",1400,1000);
  cc->Divide(2,8);
  for(int qq=1;qq<=16;qq++) cc->GetPad(qq)->SetGrid(1,0); 
  TString pdfname=Form("bxs_bd%d.pdf",board);
  TString pdfnamel=pdfname+"(";
  TString pdfnamer=pdfname+")";
  TString titl;
  TString cct;
  TString extra_cut = "";
  for(int i=1; i<IMAX; i++) {
    printf("%d/%d\n",i,IMAX-1);
    cct = Form("cc%d",i);
    for(x=0;x<8;x++) {
      //vcut = Form("ver%d*(%s==%d && goodTAC==1)",x,invar.Data(),i);
      //hcut = Form("hor%d*(%s==%d && goodTAC==1)",x,invar.Data(),i);
      vcut = Form("ver%d*(%s==%d%s)",x,invar.Data(),i,extra_cut.Data());
      hcut = Form("hor%d*(%s==%d%s)",x,invar.Data(),i,extra_cut.Data());
      titl = Form("%s %s=%d (%d)",vhn[x].Data(),invar.Data(),i,invar_arr[i]); vh[x]->SetTitle(titl.Data());
      titl = Form("%s %s=%d (%d)",hhn[x].Data(),invar.Data(),i,invar_arr[i]); hh[x]->SetTitle(titl.Data());
      tr->Project(vhn[x].Data(),"bx",vcut.Data());
      tr->Project(hhn[x].Data(),"bx",hcut.Data());
      cc->cd(2*x+1); vh[x]->Draw();
      cc->cd(2*x+2); hh[x]->Draw();
    };
    if(i==1) cc->Print(pdfnamel.Data(),"pdf");
    else if(i+1==IMAX) cc->Print(pdfnamer.Data(),"pdf");
    else cc->Print(pdfname.Data(),"pdf");
    cc->Write(cct.Data());
  };
};
