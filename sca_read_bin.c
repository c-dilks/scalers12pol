/**********************************************************************
 *
 *  to read STAR scaler binary data files
 *
 *  expects open file_descriptor as input
 *  returns arguments:
 *     number of channels, channel value, and counts in that channel
 *  return value == 1 if data good, 0 otherwise
 *
 *  note that chn=0 means val is the RHIC strobe counter
 *
 *  hank crawford 510-486-6962 hjcrawford@lbl.gov
 *
 *  4 March 2004
 *
 *********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "sca_read.h"

int sca_read_bin(int debug, FILE *infile, int *num, int *chn, unsigned long long *val) {

  static int nwrd=0;
  static int nbuf=0;
  static int numloc=0;

  int nby;
  unsigned int lowd, hiwd;
  unsigned int vlo, vhi;
  int i, *p;

  static unsigned long long sum=0;

  if (nwrd==0) {
    if (0 == (nby = fread(sca_rawdat,1,4*SCA_MAXWRDS,infile)) ) return(0);
    /*
     * 2004-09-11 - EPSichtermann@lbl.gov
     *
     * the for-loop below converts between big and little endian so that
     * the scaler data can be read at pdsf
     */
  
    for ( i = SCA_MAXWRDS, p = sca_rawdat; i-- ; p++ ) {
  	*p=(((*p>>24)&0xff)|((*p&0xff)<<24)|((*p>>8)&0xff00)|((*p&0xff00)<<8));
    }
   
    nbuf++;
    sca_board   = sca_rawdat[0];
    sca_runnum  = sca_rawdat[1];
    sca_count   = sca_rawdat[2];
    sca_maxwrds = sca_rawdat[3];

    if (debug) 
      printf("SCA_READ_BIN: buffer=%3d bd=%2d run=%4d count=%2d wds=%4d\n",
	   nbuf,sca_board, sca_runnum, sca_count, sca_maxwrds);

    if (sca_maxwrds>SCA_MAXWRDS)
      printf("#SCA_READ_BIN: maxwrds %d too large in buffer %d\n",
	     nbuf, sca_maxwrds);
    nwrd=4;
  }
  
  numloc++;

  lowd = sca_rawdat[nwrd];
  hiwd = sca_rawdat[nwrd+1];
  
  *chn = (hiwd & 0xffffff00)>>8;
  vhi    = hiwd  & 0x000000ff;
  vlo    = lowd;
  *val   = ((long long)vhi<<32) + vlo;

  if ((*chn==0)&&(debug>1)) 
    printf("SCA_READ_BIN: *chn==0 buf=%d word=%d hi=0x%8x lo=0x%8x\n", 
		      nbuf, nwrd,hiwd, lowd);
  if ((numloc>1)&&(*chn!=0)) sum += *val;


  if (debug>1) {
    printf("\n SCA_READ_BIN: num= %d lowd=%8x hiwd=%8x chn=%6x"
	   " vlo=%8x vhi=%8x val=%lld sum=%lld\n",
	   *num,lowd, hiwd,*chn,vlo,vhi,*val,sum);
    if (debug>2) if (getc(stdin)=='q') exit(0);
  }

  if (nwrd<(SCA_MAXWRDS-2)) nwrd+=2;
  else {
    if (debug) {
      printf("SCA_READ_BIN: end of buffer %3d: num=%d nwrd=%d\n",
	     nbuf,numloc,nwrd);
      if (getc(stdin)=='q') exit(0);
    }
    nwrd=0;
  }
  *num = numloc;
 
  // if (*chn==0 && *val==0) return(0);

  return(1);
}
