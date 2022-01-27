/* evallf.c caluclates a factor proportional to twice the negative log likelihood */
/* for multivariate normal data with missing values.                              */
/*                                                                                */
/* evallf.c is designed to be dynamically loaded into R version 1.1.1 or later    */
/* and used to find the MLE of mu and sigma for the MVN distribution              */
/*                                                                                */
/* instead of sigmahat (the etimate of the variance-covariance matrix) itself,    */
/* the code works on the inverse of the cholesky factor of sigma, del.            */
/* (cf. Pinheiro and Bates [2000], _Mixed-effects models in S and S-PLUS_,        */  
/* Springer.)  A maximum of 50 variables are allowed.                             */
/*                                                                                */
/* written by K. Gross, Nov 2000, with guidance from D. Bates.                    */

#include <math.h>
#define MAX_VARS 50 /*  maximum number of variables */

void evallf(double data[], int *nv, int freq[], int *nbl, int pa[], double pars[], double* val)
{
  int iblock, datacounter=0, nvars = *nv, nblocks = *nbl;
  int i,j,k,pcount,acount,muplace,idata;
  double del[MAX_VARS][MAX_VARS]={0.0}, subdel[MAX_VARS][MAX_VARS]={0.0},c,d;
  double newcol1[MAX_VARS]={0.0},newcol2[MAX_VARS]={0.0};
  double diagsum;
  double pmu[MAX_VARS]={0.0}, prod[MAX_VARS]={0.0};

  *val=0.0;

  /* construct the full del */

  for(i=0;i<nvars;i++)
    {
      del[i][i]=exp(pars[nvars+i]);
    }
  
  k=0;
  for(j=1;j<nvars;j++)
    {
      for(i=0;i<j;i++)
	{
	  del[i][j]=pars[2*nvars+k];
	  k++;
	}
    }

  /* loop through the blocks in the data set */

  for(iblock=0;iblock<nblocks;iblock++)
    {
      /* find the new precision matrix, `subdel', corresponding to the */
      /* pattern of missingness in the present block   */

      /* shuffle the rows, putting the rows corresponging */
      /* to present observations on top.                  */

      pcount=0;  /* running tallies of the number of variables (p)resent and (a)bsent */
      acount=0;

      for(i=0;i<nvars;i++)  
	{
	  if (pa[iblock*nvars+i]>0.5)  /* pa=1 means the variable is present, pa=0 means absent */
	    {
	      for(j=0;j<nvars;j++)
		{
		  subdel[pcount][j]=del[i][j];
		}
	      pcount++;
	    }
	  else
	    {
	      for(j=0;j<nvars;j++)
		{
		  subdel[nvars-acount-1][j]=del[i][j];
		}
	      acount++;
	    }
	}

      /* zero out elements below the main diagonal */
      /* using Givens rotations on the right.      */

      for(i=nvars-1;i>-1;i--) /* start with the bottom row and move up */
	{
	  for(j=0;j<i;j++)    /* in each row, start with column 1 and move right to the main daigonal */
	    {
	      /* zero out subdel[i][j] */
	      if((subdel[i][j]<.000001) & (subdel[i][j]>-.000001))
		subdel[i][j]=0;
	      else
		{
		  c = subdel[i][j]/sqrt(subdel[i][j]*subdel[i][j] + subdel[i][j+1]*subdel[i][j+1]);
		  d = subdel[i][j+1]/sqrt(subdel[i][j]*subdel[i][j] + subdel[i][j+1]*subdel[i][j+1]);
	      
		  for(k=0;k<nvars;k++)  /* calculate the new columns of subdel */
		    {
		      newcol1[k]=d*subdel[k][j]-c*subdel[k][j+1];
		      newcol2[k]=c*subdel[k][j]+d*subdel[k][j+1];
		    }

		  for(k=0;k<nvars;k++)  /* update subdel */                             
		    {
		      subdel[k][j]=newcol1[k];
		      subdel[k][j+1]=newcol2[k];
		    }
		  subdel[i][j]=0.0;
		}
	    }
	}

      /* flip column signs so that diagonal elements are all positive */

      for(i=0;i<pcount;i++)
	{
	  if (subdel[i][i]<0)
	    {
	      for(j=0;j<i+1;j++)
		{
		  subdel[j][i] = -1 * subdel[j][i];
		}
	    }
	}
      
      /* extract the relevant mu parameters */
      
      muplace=0;
      for(i=0;i<nvars;i++)
	{
	  if(pa[iblock*nvars+i]>0.5)
	    {
	      pmu[muplace]=pars[i];
	      muplace++;
	    }
	}

      /* calculate the likelihood */

      diagsum=0.0;   /* calculate the part of the likelihood due to the log of the determinant */
      for(i=0;i<pcount;i++)       
	diagsum += log(subdel[i][i]);
      
      *val -= 2*freq[iblock]*diagsum;  

      for(idata=0;idata<freq[iblock];idata++) /* calculate the part of the likelihood due to (y-mu)' * sigma^-1 * y-mu */
	{
	  for(j=0;j<pcount;j++)
	    prod[j]=0.0;

	  for(j=0;j<pcount;j++)
	    {
	      for(k=0;k<j+1;k++)
	        prod[j] += (data[datacounter+idata*pcount+k]-pmu[k])*subdel[k][j];
	    }  

	  for(j=0;j<pcount;j++)
	    *val += prod[j]*prod[j];
	}
    datacounter += pcount*freq[iblock];  

    }
}
