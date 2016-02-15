// checks the counts.root file produced from the new bunch kicking algorithm (bunch_kicker)

void bunch_kicker_check(const char * filename="counts.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sca");

  Int_t fill,fill_tmp;
  Int_t index;
  tr->SetBranchAddress("fill",&fill);
  tr->SetBranchAddress("i",&index);
  Int_t sb[4];
  Int_t sb_kicked[4];
  char cut[4][64];
  char cut_kicked[4][64];

  /*
  TFile * outfile = new TFile("bktest.root","RECREATE");
  TTree * outr = new TTree();
  outr->Branch("fill",&fill,"fill/I");
  outr->Branch("sb0",&(sb[0]),"sb0/I");
  outr->Branch("sb1",&(sb[1]),"sb1/I");
  outr->Branch("sb2",&(sb[2]),"sb2/I");
  outr->Branch("sb3",&(sb[3]),"sb3/I");
  outr->Branch("sb0_kicked",&(sb_kicked[0]),"sb0_kicked/I");
  outr->Branch("sb1_kicked",&(sb_kicked[1]),"sb1_kicked/I");
  outr->Branch("sb2_kicked",&(sb_kicked[2]),"sb2_kicked/I");
  outr->Branch("sb3_kicked",&(sb_kicked[3]),"sb3_kicked/I");
  */

  fill_tmp=0;
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    if(fill!=fill_tmp)
    {
      sprintf(cut[0],"i==%d && blue==-1 && yell==-1",index);
      sprintf(cut[1],"i==%d && blue==-1 && yell==1",index);
      sprintf(cut[2],"i==%d && blue==1 && yell==-1",index);
      sprintf(cut[3],"i==%d && blue==1 && yell==1",index);
      for(Int_t s=0; s<4; s++) 
      {
        sprintf(cut_kicked[s],"%s && !kicked",cut[s]);
        sb[s] = tr->GetEntries(cut[s]);
        sb_kicked[s] = tr->GetEntries(cut_kicked[s]);
      };
      //outr->Fill();
      fill_tmp = fill;
      printf("%d\n",fill);
      printf("  %d  %d  %d  %d\n",sb[0],sb[1],sb[2],sb[3]);
      printf("  %d  %d  %d  %d\n",sb_kicked[0],sb_kicked[1],sb_kicked[2],sb_kicked[3]);
      if(sb_kicked[0]!=sb_kicked[1] || 
         sb_kicked[0]!=sb_kicked[2] || 
         sb_kicked[0]!=sb_kicked[3]) printf("                  ^^ not equalized\n");
      printf("\n");
    };
  };
  //outr->Write("outr");
}
