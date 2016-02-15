// NOTE: this is used by bunch_kicker
//
// generates 10 unique random numbers from 0 to max-1


void RNG(Int_t max,Int_t remove, Int_t spin)
{
  if(max<remove)
  {
    fprintf(stderr,"max is less than remove=number to generate\n");
    return;
  };
  TRandom3 * r3 = new TRandom3(0); // Mersenne twistor method, with seed set by UUID
  Int_t arr[120];
  Bool_t uniq=0;
  char rand_file_name[32];
  sprintf(rand_file_name,"rand_file_%d",spin);
  for(Int_t i=0; i<remove; i++)
  {
    while(!uniq)
    {
      arr[i] = (Int_t) ((r3->Rndm())*max);
      uniq=1;
      for(Int_t j=0; j<i; j++)
      {
        if(arr[j] == arr[i]) uniq=0;
      };
    };
    uniq=0;
    gSystem->RedirectOutput(rand_file_name,"a");
    printf("%d\n",arr[i]);
    gSystem->RedirectOutput(0);
  };
};
