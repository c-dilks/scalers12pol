void DrawConsistency()
{
  TFile * cons_file = new TFile("cons_file.root","RECREATE");
  TTree * tr = new TTree();
  // cons is run number or fill number, depending on what study you're doing
  tr->ReadFile("cons_study","cons/I:tbit/I:cbit/I:spinbit/I:p0/D:p0_err/D");
  tr->Write("tr");
}
