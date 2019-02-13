/*
 This software was developed by employees of the US Naval Research Laboratory
(NRL), an agency of the Federal Government. Pursuant to title 17 section 105 of
the United States Code, works of NRL employees are not not subject to copyright
protection, and this software is in the public domain. FiPy is an experimental
system. NRL assumes no responsibility whatsoever for its use by other parties,
and makes no guarantees, expressed or implied, about its quality, reliability,
or any other characteristic. We would appreciate acknowledgment if the software
is used.

To the extent that NRL may hold copyright in countries other than the United
States, you are hereby granted the non-exclusive irrevocable and unconditional
right to print, publish, prepare derivative works and distribute this software,
in any medium, or authorize others to do so on your behalf, on a royalty-free
basis throughout the world.

You may improve, modify, and create derivative works of the software or any
portion of the software, and you may copy and distribute such modifications or
works. Modified works should carry a notice stating that you changed the
software and should note the date and nature of any such change. Please
explicitly acknowledge the US Naval Research Laboratory as the original source.

This software can be redistributed and/or modified freely provided that any
derivative works bear some notice that they are derived from it, and any
modified versions bear some notice that they have been modified.

 
Author: David Rowenhorst; The US Naval Research Laboratory
Date: 5 Dec 2018   
 

Description: This contains the primary code for performning the Non-Local 
Pattern Averaging Reindexing (NLPAR) algorithm.  
 */


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <map>
//#include <limits>
#include "idl_export.h" // this is only needed for converting IDL types.  
//#include "cNL_Means.h"
#include <math.h>
#ifdef _OPENMP
    #include <omp.h>
#endif


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void getpairid(long long int* index0, long long int* index1, long long unsigned int* pairid) {
	unsigned int* p;
	p = (unsigned int*)(pairid);

	if (*index0 < *index1){
    p[0] = *index0; 
    p[1] = *index1;
  } else{
    p[1] = *index0;
    p[0] = *index1;
  }
}

int SigEst ( float* patterns, long long* patDim, long long* szim, long long* xy_start_end, 
      long long* win_sz, float* sigma, int* nthread){
  //int nthread=1;
  long long int ncol, nrow, patStride;
  long long int istart, iend, jstart, jend;
  long long int npatch_el, nwin_el;
  double  patStrideFlt;

  
  
  

  ncol =  szim[0];
  nrow =  szim[1]; 
  patStride = patDim[0];
  patStrideFlt = (double) patDim[0];
  //printf("%i \t %i \n", (int) ncol,(int) nrow);
  istart = (long long) xy_start_end[0];
  iend = (long long) xy_start_end[1]+1;
  jstart = (long long) xy_start_end[2];
  jend = (long long) xy_start_end[3]+1;
  

  #ifdef _OPENMP
    //  omp_set_dynamic(1);
    if (nthread[0] <= 0){
      nthread[0] = MIN( (int)floor(((float)(iend-istart))/10.0), omp_get_num_procs());
      nthread[0] = MAX(1, nthread[0]);
    }
    omp_set_num_threads(nthread[0]);
       //omp_set_num_threads(MIN(10, (iend-istart)));
       //omp_set_num_threads(10);
  #endif
  

  
  nwin_el =((long long) (2*win_sz[0])+1)*(2*win_sz[1]+1);
  //printf("%i \n", nwin_el);
  //patch_r = (long) (patch_sz[0]-1)/2;
  
  

#pragma omp parallel default(shared)
  {   
  // allocate a bunch of variables that will only be used inside the parallel
  // contruction
  //long long unsigned int pairid; 
  //unsigned int* pairidfirst;
  double d,d2, sum, mind, sigma0, sigma1; 
  
  long long int  index0, index1, counter, index01, w0indx; 
  long long int m,n, q, i,j;
  long long int winstart_x, winend_x, winstart_y, winend_y;
  bool test;
  //double inf =std::numeric_limits<double>::infinity();
  

#pragma omp for schedule(static)
  for ( i=istart; i< iend;i++){
    for (j = jstart;j < jend;j++){
      
      index0 = (i+j*ncol)*patStride;
      //printf("%i \n", i+j*ncol);
      //sigma[i+j*ncol]=0.0;

      //for (q=0; q < patStride;q++){
      //  image_out[index0+q] = 0.0;
      //}
      //printf(" %i \t %i \n", (int) i,(int) j);
      //printf(" %i \n", index0);
      
      
      winstart_x = MAX((i-win_sz[0]),0) - MAX( (i+win_sz[0] - (ncol-1)), 0);
      winend_x = MIN((i+win_sz[0]),(ncol-1)) + MAX( (win_sz[0]-i), 0) +1;
    
      winstart_y = MAX((j-win_sz[1]),0) - MAX((j+win_sz[1] - (nrow-1)), 0);
      winend_y = MIN((j+win_sz[1]),(nrow-1)) + MAX( (win_sz[1]-j), 0) +1;
      //printf("%i \t %i \t%i \t %i \n", winstart_x, winend_x, winstart_y, winend_y);
      //printf("%i \n", nwin_el);
      //counter = 0;
      mind =  INFINITY;
      for (n = winstart_y; n < winend_y; n++){
        for (m = winstart_x; m < winend_x; m++){    
          //w[counter] = 0.0;
          //counter += 1;
          index1 = (m+n*ncol)*patStride;
          
          if (index0 == index1) { continue; } else {
            /*if (index0 < index1)
            {
              p[0] = index0; //assign values to 
              p[1] = index1;
            }
            else
            {
              p[1] = index0;
              p[0] = index1;
            }
            */
            //getpairid(&index0, &index1, &pairid);

            //#pragma omp critical 
            //{test = pointdist.find(pairid) == pointdist.end();}
            
            
          
          d2 = 0.0;
          for (q=0; q < patStride;q++){
            d =  patterns[index0+q] - patterns[index1+q];
            d2 += d * d;
          }
          if (d2 > 1.0){ //sometimes EDAX collects the same pattern twice
            d2 /= (patStrideFlt*2.0);
            d2 = sqrt(d2);
            if (d2 < mind){
              mind = d2;
            }
          }
          }


          //d2 -= patStrideFlt*(sigma1 + sigma0);
          //d2 /= (sigma1 + sigma0)*0.5;
          //d2 /= (sigma1 + sigma0)*sqrt(2.0*patStrideFlt);
          //printf("%f \n", d2);
          //#pragma omp critical
          //{pointdist[pairid] = d2 ;}
          
       }
      }//end window scan

      sigma[i+j*ncol]=mind;
              
    } 
  } //kjhj
} //end of parallel construction
  return 1;
}


int cNLPARV5 ( float* patterns, long long* patDim, long long* szim, long long* xy_start_end, 
		  long long* win_sz, float* sigma, float* lambda, float* dthresh, 
      float* pat_out, float* w0, int* nthread){
	//int nthread=1;
	long long int ncol, nrow, patStride;
  long long int istart, iend, jstart, jend;
	long long int npatch_el, nwin_el;
	double lambda2, patStrideFlt;

	
	// pre calc a lambda factor
  lambda2 = 1.0 / (double) (lambda[0]*lambda[0]);

  ncol =  szim[0];
	nrow =  szim[1]; 
	patStride = patDim[0];
  patStrideFlt = (double) patDim[0];
	//printf("%i \t %i \n", (int) ncol,(int) nrow);
  istart = (long long) xy_start_end[0];
  iend = (long long) xy_start_end[1]+1;
  jstart = (long long) xy_start_end[2];
  jend = (long long) xy_start_end[3]+1;
  
  // calculate how many threads should be used.  
  //This has been optimized only on a small dataset. 
  // YRMV
  #ifdef _OPENMP
		//	omp_set_dynamic(1);
    if (nthread[0] <= 0){
      nthread[0] = MIN( (int)floor(((float)(iend-istart))/10.0), omp_get_num_procs());
      nthread[0] = MAX(1, nthread[0]);
    }
    omp_set_num_threads(nthread[0]);
	     //omp_set_num_threads(MIN(10, (iend-istart)));
		   //omp_set_num_threads(10);
	#endif
  // number of patterns within the search window.   
	nwin_el =((long long) (2*win_sz[0])+1)*(2*win_sz[1]+1);

	
	

#pragma omp parallel default(shared)
  {		
	// allocate a bunch of variables that will only be used inside the parallel
	// contruction
	long long unsigned int pairid; 
	unsigned int* pairidfirst;
	double d,d2, sum, mind, sigma0, sigma1; 
	double  *w;
	long long int  index0, index1, counter, index01, w0indx; 
	long long int m,n, q, i,j;
	long long int winstart_x, winend_x, winstart_y, winend_y;
  bool test;
  
  // must allocate within the parallel construction - must also free within
  w = (double *) malloc(nwin_el * sizeof(double));
  // going to use a simple hash to store the distances between two patterns. 
  // This avoids the double calculation of the distances	
  pairidfirst = (unsigned int*)(&pairid);
  std::map<long long unsigned int, double> pointdist;

#pragma omp for schedule(static)
	for ( i=istart; i< iend;i++){
		for (j = jstart;j < jend;j++){
			
			index0 = (i+j*ncol)*patStride;
			//printf("%i \n", i+j*ncol);
      sigma0 = sigma[i+j*ncol]*sigma[i+j*ncol];
      // intialize pat_out to 0
			for (q=0; q < patStride;q++){
				pat_out[index0+q] = 0.0;
			}
			
			
			// calculate the start and end of the sliding window.  
			winstart_x = MAX((i-win_sz[0]),0) - MAX( (i+win_sz[0] - (ncol-1)), 0);
			winend_x = MIN((i+win_sz[0]),(ncol-1)) + MAX( (win_sz[0]-i), 0) +1;
		
			winstart_y = MAX((j-win_sz[1]),0) - MAX((j+win_sz[1] - (nrow-1)), 0);
			winend_y = MIN((j+win_sz[1]),(nrow-1)) + MAX( (win_sz[1]-j), 0) +1;
			
			counter = 0;
			for (n = winstart_y; n < winend_y; n++){
				for (m = winstart_x; m < winend_x; m++){		
					w[counter] = 0.0;
					counter += 1;
					index1 = (m+n*ncol)*patStride;
          
					if (index0 == index1) { continue; }
            
  					getpairid(&index0, &index1, &pairid);

            // test and see if the distance has been calculated.  
            {test = pointdist.find(pairid) == pointdist.end();}
            // true is that there is no previous value. 
            if(test){
              sigma1 = sigma[m+n*ncol]*sigma[m+n*ncol];
            	d2 = 0.0;
              for (q=0; q < patStride;q++){
  							d =  patterns[index0+q] - patterns[index1+q];
  							d2 += d * d;
  						}
              // normalize the distance. 
  						d2 -= patStrideFlt*(sigma1 + sigma0); 
              //d2 /= (sigma1 + sigma0)*0.5;
              d2 /= (sigma1 + sigma0)*sqrt(2.0*patStrideFlt);
              //printf("%f \n", d2);
              //#pragma omp critical
              {pointdist[pairid] = d2 ;}
            }

			 }
			}//end window scan
      // now calculate the weights
			sum = 0.0;
			counter = 0;
			for (n = winstart_y; n < winend_y; n++){
				for (m = winstart_x; m < winend_x; m++){
					index1 = (m + n*ncol)*patStride;
					if (index1 == index0){
						w[counter] = 1.0;
            w0indx = counter;
					} else {
						getpairid(&index0, &index1, &pairid);
						w[counter] = pointdist[pairid];
						w[counter] = MAX(w[counter]-(*dthresh), 0.0);
						w[counter] = exp(-1.0*w[counter]*lambda2);
					}
					sum += w[counter];
					counter +=1; 
				}
			}
			
      // normalize the weights to total = 1.0
			for (q=0; q < nwin_el; q++){
				w[q] = w[q]/(sum);
			}	
      //write out the weight factor for the pattern of interest. 
      w0[(i+j*ncol)] = (float) w[w0indx];
			
      // now apply the weights and store in pat_out
			counter = 0;
			for (n = winstart_y; n < winend_y; n++){
				for (m = winstart_x; m < winend_x; m++){
					index1 = (m + n*ncol)*patStride;

					for (q=0; q < patStride;q++){
						pat_out[index0+q] +=  patterns[index1+q] * (float)w[counter];
					}
					counter += 1;
				}
			}
							
		}	
	} //kjhj
	
	//printf("%i \t %i \t %i \n", index0, index1, index01);

	//free(d2);
	//d2 = NULL;
	free(w); 
	w = NULL;
} //end of parallel construction
	return 1;
}

#ifdef __cplusplus  /* Need to avoind name mangling for IDL */
extern "C"
{
#endif

int cNLPARV5_wrap (int argc, void* argv[]) {
	/* 
	 *
	 * argc must equal 11
	 * argv[0] (float*) pointer to orginal pattern array (floating)
	 * argv[1] (IDL_LONG*) pointer to the number of elements in 1 pattern.  
	 * argv[2] (IDL_LONG*) array w size scanned image (ncol, nrows)
	 * argv[3] (IDL_LONG*) array with the starting and stoping points of the nl_means (xstart, xend, ystart, yend)
	 * argv[4] (IDL_LONG*) the 1/2 width size of the search window
	 * argv[5] (float*) sigma array
	 * argv[6] (float*) h
	 * argv[7] (float*) threshold distance that will be treated with highest weighting factor
	 * argv[8] (float*) output array
   * argv[9] (float*) output weight 
   * argv[10] (IDL_LONG*) number of threads to feed openmp
	 */

	 
	if(argc != 11) return 0;
	//long long nMaskindex[1], ii;
	//long long *maskindex;
	long long patDim[1], imageDim[2], xystartend[4], windowSize[2], ii;
	int returnval, nthread;
	//IDL_LONG *maskindexptr;
  IDL_LONG *temp; 
  IDL_LONG *temp2;
  

  temp = (IDL_LONG*) argv[1];
  patDim[0] = (long long) temp[0];
  temp = (IDL_LONG*) argv[2];
  imageDim[0] = (long long) temp[0];
  imageDim[1] = (long long) temp[1];
  temp = (IDL_LONG*) argv[3];
  for (ii = 0; ii < 4; ii++) xystartend[ii] = (long long) temp[ii];
  //printf( " %i \t %i \t %i \t %i \n", xystartend[0], xystartend[1], xystartend[2], xystartend[3]);
  temp = (IDL_LONG*) argv[4];
	windowSize[0] = (long long) temp[0];
	windowSize[1] = (long long) temp[1];
	//printf( " %i \t %i \n", windowSize[0], windowSize[1]);

	//temp2 = (IDL_LONG*) argv[5];
	//nMaskindex[0] = (long long) *temp2;
	//printf(" %i \n", nMaskindex[0]);
	
	//maskindex = (long long *) malloc(nMaskindex[0] * sizeof(long long));
	//maskindexptr = (IDL_LONG*) argv[4];
	//for (ii = 0; ii < nMaskindex[0]; ii++) {maskindex[ii] = (long long) (maskindexptr[ii]);}

	//printf(" %i \n", maskindex[nMaskindex[0]-1]);
  temp = (IDL_LONG*) argv[10];
  nthread = (int) temp[0];

	returnval = cNLPARV5( (float*) argv[0], patDim, imageDim, xystartend,  windowSize , 
		    (float*) argv[5], (float*) argv[6], (float*) argv[7], (float*) argv[8], (float*) argv[9], &nthread);

  
	return returnval;

	
}

int cSigEst_wrap (int argc, void* argv[]) {
  /* 
   *
   * argc must equal 7
   * argv[0] (float*) pointer to orginal pattern array (floating)
   * argv[1] (IDL_LONG*) pointer to the number of elements in 1 pattern.  
   * argv[2] (IDL_LONG*) array w size scanned image (ncol, nrows)
   * argv[3] (IDL_LONG*) array with the starting and stoping points of the nl_means (xstart, xend, ystart, yend)
   * argv[4] (IDL_LONG*) the 1/2 width size of the search window
   * argv[5] (float*) sigma array output
   * argv[6] (IDL_LONG*) number of threads to feed openmp
   */

   
  if(argc != 7) return 0;
  
  long long patDim[1], imageDim[2], xystartend[4], windowSize[2], ii;
  int returnval, nthread;
  
  IDL_LONG *temp; 
  IDL_LONG *temp2;
  

  temp = (IDL_LONG*) argv[1];
  patDim[0] = (long long) temp[0];
  temp = (IDL_LONG*) argv[2];
  imageDim[0] = (long long) temp[0];
  imageDim[1] = (long long) temp[1];
  temp = (IDL_LONG*) argv[3];
  for (ii = 0; ii < 4; ii++) xystartend[ii] = (long long) temp[ii];
  
  temp = (IDL_LONG*) argv[4];
  windowSize[0] = (long long) temp[0];
  windowSize[1] = (long long) temp[1];
  
  temp = (IDL_LONG*) argv[6];
  nthread = (int) temp[0];

  returnval = SigEst( (float*) argv[0], patDim, imageDim, xystartend,  windowSize, 
        (float*) argv[5], &nthread);

  
  return returnval;

  
}


#ifdef __cplusplus /* Need to avoind name mangling for IDL */
}
#endif

