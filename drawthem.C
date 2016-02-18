// draws all tcanvases from bxs_bd*.root (opens 50 windows!!)

void drawthem(Int_t bd=10)
{
  TString fname = Form("bxs_bd%d.root",bd);
  TFile * f = new TFile(fname.Data(),"READ");
  TString n;
  const Int_t NFILLS=49;
  TCanvas * c[NFILLS];
  for(int x=0; x<NFILLS; x++) {
    n=Form("cc%d",x);
    c[x]=(TCanvas*) f->Get(n.Data());
    c[x]->SetName(n.Data());
    c[x]->ToggleEventStatus();
    c[x]->Draw();
  };
};


