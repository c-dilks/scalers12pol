// prints run index vs. run and fill index vs. fill

void print_run_table(const char * filename="sums.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sum");

  tr->SetScanField(10000);

  gROOT->ProcessLine(".! touch tab_run.dat; rm tab_run.dat");
  gROOT->ProcessLine(".! touch tab_fill.dat; rm tab_fill.dat");

  gSystem->RedirectOutput("tab_run.dat");
  tr->Scan("i:runnum:fi:fill");
  gSystem->RedirectOutput(0);
}
