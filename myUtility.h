#ifndef _MYUTILITY_H_
#define _MYUTILITY_H_

#include "stdio.h"
#include "math.h"
//#include "conio.h"
#include "stdlib.h"

float*** alloc3DMatrix( int d1, int d2, int d3 )
{
	float ***arr = new float **[d1];
	for (int i = 0; i < d1; i++) {
		arr[i] = new float *[d2];
		for (int j = 0; j < d2; j++) {
		arr[i][j] = new float[d3];
		memset( arr[i][j], 0, sizeof(float)*d3 );
		}
	}

	return arr;
}

void free3DMatrix(float*** arr, int d1, int d2, int d3)
{
	for (int i = 0; i < d1; i++) {
		for (int j = 0; j < d2; j++) {
		free(arr[i][j]);
		}
		free(arr[i]);
	}
	free(arr);
}


void intepolateTransferFunctionFromControlPoint( float* ctfa, int bins, float* ctPointV, float* ctPointA )
{
	//calculate transfer function
	int cpNow = 0;
	for( int i=0; i<bins; i++ ){
		float vxl = i/(float)(bins-1.0);
		//current vxl > vxl in current control point, move to next control point
		if( vxl > ctPointV[cpNow+1] )cpNow++;
		//ctfr[i] = interpolation1D(vxl, g_ctrlPtV[cpNow], g_ctrlPtV[cpNow+1],
		//								g_ctrlPtR[cpNow], g_ctrlPtR[cpNow+1]);//r
		//ctfg[i] = interpolation1D(vxl, g_ctrlPtV[cpNow], g_ctrlPtV[cpNow+1],
		//								g_ctrlPtG[cpNow], g_ctrlPtG[cpNow+1]);//g
		//ctfb[i] = interpolation1D(vxl, g_ctrlPtV[cpNow], g_ctrlPtV[cpNow+1],
		//								g_ctrlPtB[cpNow], g_ctrlPtB[cpNow+1]);//b

		float p, ps, pt, vs, vt;
		p = vxl;
		ps = ctPointV[cpNow];
		pt = ctPointV[cpNow+1];
		vs = ctPointA[cpNow];
		vt = ctPointA[cpNow+1];

		float ratio = (p-ps)/(pt-ps);
		ctfa[i] = ((vt-vs)*ratio + vs);
	}
}

float computeEMD( int bins, float* h1, float* h2 )
{
	float emd = 0;
	float h1cdf = 0;
	float h2cdf = 0;
	for( int j=0; j<bins; j++ ){
		h1cdf += h1[j];
		h2cdf += h2[j];
		emd += fabs(h1cdf - h2cdf);
	}

	//euclidean distance
	//for( int i=0; i<bins; i++ )
	//	emd += (h1[i] - h2[i])*(h1[i] - h2[i])*(h1[i] - h2[i])*(h1[i] - h2[i]);
	//emd = sqrt(sqrt( emd ));
		
	return emd;
}



#endif