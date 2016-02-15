/* Chris Perkins - <cwperkins@lbl.gov> - 2011115
 *
 * Example reader for ScalerII Histogrammed Data Files
 *
 */

#include <stdio.h>
#include <fcntl.h>
#include <string.h>

struct hist_hdr {
  unsigned int version;
  unsigned int bdnum;
  unsigned int runnum;
  unsigned int count;
  unsigned int num_words;
  unsigned int rs_count_lo;
  unsigned int rs_count_hi;
  unsigned int status1;
  unsigned int status2;
  unsigned int status3;
};
struct sca_val{
  int bx;                        // bunch crossing
  unsigned long long valBBC[8];
  unsigned long long valZDC[8];
  unsigned long long valVPD[4];
  unsigned long long valSum;
};

int main(int argc, char *argv[]) {
  struct hist_hdr hist_hdr_l;
  char filename_l[100];
  int fd;
  int bytes_read;
  int i;
  int j;
  unsigned int chn_cnt[3];
  unsigned long long channel;
  unsigned long long count;
  unsigned long long sum = 0;
  int debug_l = 0;

  char outfile[100];
  FILE *fout;
  struct sca_val scaler[128];
  int bunch;
  int chnl_bbc;
  int chnl_zdc;
  int chnl_vpd;
  //  unsigned long long val[120];

  //Initialize scaler
  for(i=0;i<120;i++)
  {
    scaler[i].bx=i;
    for(j=0;j<8;j++)
    {
      scaler[i].valBBC[j]=0;
      scaler[i].valZDC[j]=0;
      if(j<4)
        scaler[i].valVPD[j]=0;
    }
    scaler[i].valSum=0;
  }

  printf("Usage: scaler2_reader [filename] [debug]\n");

  if(argc > 1) strcpy(filename_l, argv[1]);
  else sprintf(filename_l, "out.dat");
  if(argc > 2) debug_l = atoi(argv[2]);
  printf("Opening datafile: %s\n", filename_l);

  // Open File
  fd = open(filename_l, O_RDONLY);
  if(fd == -1) {
    printf("Error opening datafile: %s\n", filename_l);
    return -1;
  }

  // Read Version
  bytes_read = read(fd, &hist_hdr_l, 4);
  if(hist_hdr_l.version == 1) {    // Data Format Version 1
    // Read Rest of Header
    bytes_read = read(fd, &(hist_hdr_l.bdnum), sizeof(struct hist_hdr) - 4);
    if(bytes_read != (sizeof(struct hist_hdr) - 4)) {
      printf("Error reading header\n");
      close(fd);
      return -1;
    }

    // Print Header
    printf("Runnumber = %8.8d, BoardNumber = %2.2d\n", hist_hdr_l.runnum, hist_hdr_l.bdnum);
    printf("  Data Format Version = %d, Count = %d\n", hist_hdr_l.version, hist_hdr_l.count);
    printf("  Num Data Words = %d (0x%8.8x)\n", hist_hdr_l.num_words, hist_hdr_l.num_words);
    printf("  Num Channels = %d (0x%8.8x),  RS Count = 0x%8.8x%8.8x\n", hist_hdr_l.num_words / 3, hist_hdr_l.num_words / 3, hist_hdr_l.rs_count_hi, hist_hdr_l.rs_count_lo);
    printf("  Status Words = 0x%8.8x  0x%8.8x  0x%8.8x\n", hist_hdr_l.status1, hist_hdr_l.status2, hist_hdr_l.status3);


    sprintf(outfile,"datfiles/run%d_%d.dat", hist_hdr_l.runnum,hist_hdr_l.bdnum);
    fout=fopen(outfile,"w");

    // Read Channels and Counts
    for(i=0; i < (hist_hdr_l.num_words / 3) ; i++) {
      bytes_read = read(fd, chn_cnt, 12);
      if(bytes_read != 12) {
        printf("Error reading Channel/Count\n");
        close(fd);
        printf("%d\n", i);
        return -1;
      }
      channel = ((unsigned long long)(chn_cnt[2]) << 16) | ((unsigned long long)(chn_cnt[1]) >> 16);
      count = (((unsigned long long)(chn_cnt[1]) & 0xffff) << 32) | (unsigned long long)(chn_cnt[0]);

      /*------Analyze Channel------*/
      //Read channel
      bunch = (int) (channel >> 25) & 0x7f;
      printf("bx=%d\n",bunch);
      //printf("%d\n",bunch);
      if(bunch < 120)
      {

        //BBC Coincidence
        chnl_bbc = (int) channel & 0x3;
        chnl_zdc = (int) (channel >> 10) & 0x7;
        chnl_vpd = (int) (channel >> 4) & 0x3;
        if(channel & 0x2000) chnl_bbc = chnl_bbc + 4;
        //      printf("%d\n", chnl_bbc);
        scaler[bunch].valBBC[chnl_bbc] += count;
        scaler[bunch].valZDC[chnl_zdc] += count;
        scaler[bunch].valVPD[chnl_vpd] += count;
        // val[bunch] += count;
        //channel sum
        scaler[bunch].valSum += count;
      }else
      {
        printf("Bunch Crossing %d Out of Range\n",bunch);
      }
      /*------Done------*/
      sum += count;
      if(debug_l) printf("Channel 0x%12.12llx = 0x%12.12llx Counts\n", channel, count);
    }
    printf("\n");
    printf("RS Count = 0x%8.8x%8.8x, Count Sum = 0x%16.16llx\n", hist_hdr_l.rs_count_hi, hist_hdr_l.rs_count_lo, sum);
  } else {
    printf("Unknown Data File Format Version: %d\n", hist_hdr_l.version);
  }

  close(fd);


  for(i=0;i<120;i++)
  {
    fprintf(fout, "%d ", scaler[i].bx);
    for(j=0;j<8;j++)
    {
      fprintf(fout, "%lld ", scaler[i].valBBC[j]);
    }
    //      fprintf(fout, "%lld\n", scaler[i].valSum);

    for(j=0;j<8;j++)
    {
      fprintf(fout, "%lld ", scaler[i].valZDC[j]);
    }
    //      fprintf(fout, "%lld\n", scaler[i].valSum);

    for(j=0;j<4;j++)
    {
      fprintf(fout, "%lld ", scaler[i].valVPD[j]);
    }
    fprintf(fout, "%lld\n", scaler[i].valSum);
  } 

  fclose(fout);

  return 0;
}
