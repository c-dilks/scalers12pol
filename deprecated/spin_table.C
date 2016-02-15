// DEPRECATED --> use spin_cogging bash script in stead of this
// why this produces incorrect spin patterns needs to be investigated 
// prior to analyses of spin-dependent data!!
// -- see log.07.01.14 for further details 
// -- see cogging.ods for cogging info
//
// reads spin patterns from $FMSTXT/spinpat CDEV files, corrects
// for cogging at STAR, and outputs ./spinpat/[fill#].spin files
// with columns [bx][blue spin][yell spin]

void spin_table(Int_t first_fill=17333, Int_t last_fill=17601, 
                Int_t first_run=14096078, Int_t last_run=14161020)
{
  // set environment
  TString str=gSystem->ExpandPathName("$SETFMSENV");
  if(str!="SETFMSENV"){printf("source SetFMSEnv first");exit();};
  TString FMSROOT =gSystem->ExpandPathName("${FMSROOT}");
  TString START_MACRO=FMSROOT+"/start.C";
  gROOT->Macro(START_MACRO);
  p_files->Print();

  // remove old *.spin files
  gROOT->ProcessLine(".! rm spinpat/*.spin");

  // read list of fills
  Fill * fills = new Fill(p_files,first_fill,last_fill,first_run,last_run);
  //gROOT->ProcessLine(".! cat fill.txt | awk '{print $2}' | uniq > fill_list.txt");
  TTree * fill_tree = new TTree("fill_tree","fill_tree");
  fill_tree->ReadFile("fill.txt","runnum/I:fill/I");
  Int_t fill_num;
  fill_tree->SetBranchAddress("fill",&fill_num);
  //fill_tree->Scan("fill");

  // create *.spin files
  Int_t blue_spin;
  Int_t yell_spin;
  char outfile[64];
  Int_t fill_num_tmp=0;
  for(Int_t i=0; i<fill_tree->GetEntries(); i++)
  {
    fill_tree->GetEntry(i);
    if(fill_num!=fill_num_tmp)
    {
      fill_num_tmp = fill_num;
      fills->SetFillNumber(fill_num);
      sprintf(outfile,"spinpat/%d.spin",fill_num);
      // setting bx bounds 0-119, equivalent to OFile Bunchid7bit
      // -- additionally, fill->YellowSpin(120) always returns 0 in run13
      //    while bx 0 is nontrivial
      for(Int_t bx=0; bx<120; bx++)
      {
        blue_spin = fills->BlueSpin(bx);
        yell_spin = fills->YellowSpin(bx);
        gSystem->RedirectOutput(outfile);
        printf("%d %d %d\n",bx,blue_spin,yell_spin);
        gSystem->RedirectOutput(0);
      };
      printf("fill %d\n",fill_num);
    };
  };

  printf("spinpat/*.spin files created\n");
};
