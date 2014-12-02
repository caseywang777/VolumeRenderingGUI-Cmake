#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include "stdio.h"
#include "math.h"
//#include "conio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector
using namespace std;

typedef struct _SelfInfo{int bin;float seflInfo;}SelfInfo;
bool selfInfoSort (SelfInfo i, SelfInfo j) { return (i.seflInfo>j.seflInfo); }

using namespace std;

//for sampling data from histogram
typedef struct _SamplingPoint{
	int order;
	float randNumber;
	float sampleValue;
}SamplingPoint;

bool samplingVectorSortByRandNumber( SamplingPoint a, SamplingPoint b )
{
	return (a.randNumber < b.randNumber );
}

bool samplingVectorSortByOrder( SamplingPoint a, SamplingPoint b )
{
	return (a.order < b.order );
}

class Histogram
{
public:
	//constructor: create a empty(all 0)histogram
	Histogram( int bins )
	{
		m_sum = 0;
		m_bins = bins;
		m_hist = (float*)malloc( sizeof(float)* m_bins );
		memset( m_hist, 0, sizeof(float)*m_bins );
	}

	//constructor: give a histogram by a float array
	Histogram( int bins, float* hist )
	{
		m_sum = 0;
		m_bins = bins;
		m_hist = (float*)malloc( sizeof(float)* m_bins );
		memset( m_hist, 0, sizeof(float)*m_bins );

		for( int i=0; i<=m_bins; i++ ){
			m_hist[i] = hist[i];
		}
	}

	//constructor: given sample points and create the histogram
	//bins: bin number
	//sample: number of samples
	//lowRange highRange: data range
	//data: 1D float array: data points
	Histogram( int bins, int samples, float lowRange, float highRange, float* data )
	{
		m_sum = 0;
		m_bins = bins;
		m_hist = (float*)malloc( sizeof(float)* m_bins );
		memset( m_hist, 0, sizeof(float)*m_bins );

		float range = highRange - lowRange;
		for( int i=0; i<samples; i++ ){
			int idx = (int)( ((data[i]-lowRange)/range) * bins);
			if( idx == bins )idx--;

			m_hist[idx] ++;
		}
	}
	
	//constructor: given a Histogram object, copy it and create this new one
	Histogram( Histogram *a )
	{
		m_sum = a->m_sum;
		m_bins = a->m_bins;
		m_hist = (float*)malloc( sizeof(float)* m_bins );
		memcpy( m_hist, a->m_hist, sizeof(float)*m_bins );
	}

	//destructor
	~Histogram()
	{
        //printf("freeeeeeeeeee  %f %f\n", m_hist[0], m_hist[1]);
		if( m_hist != NULL )
			free( m_hist );
	}

	//reset: reset all bins to zero
	//b: indicate a bin
	void reset( )
	{
		for( int i=0; i<m_bins; i++ )
			m_hist[i] = 0;
	}


	//getBin: get a bin value
	//b: indicate a bin
	float getBin( int b )
	{
		return m_hist[b];
	}

	//setBin: set a bin value
	//b: indicate a bin
	//v: bin value
	void setBin( int b, float v )
	{
		m_hist[b] = v;
	}

	//getBins: get bins in a row (include s and t position)
	//hist: return histogram the size should be t-s+1
	//s: start index
	//t: till index
	void getBins( float* hist, int s, int t )
	{
		for( int i=s; i<=t; i++ ){
			hist[i-s] = m_hist[i];
		}
	}

	//setBins: set bins in a row (include s and t position)
	//s: start index
	//t: till index
	//hist: input bin value, the size should be t-s+1
	void setBins( int s, int t, float* hist )
	{
		for( int i=s; i<=t; i++ ){
			m_hist[i] = hist[i-s];
		}
	}

	//getNorm: return Norm of this histogram(vector)
	float getNorm()
	{
		float norm = 0;
		for( int i=0; i< m_bins; i++ )
			norm += m_hist[i]*m_hist[i];
		return sqrt( norm );
	}

	//given s and t: bin index
	//return number of samples between s and t bin(include s and t )
	float getNSamples( int s, int t )
	{
		if( s > t ) return 0;

		float sum = 0;
		for( int i=s; i<=t; i++ )
			sum += m_hist[i];
		return sum;
	}

	//normalize: normalize this histogram(vector)
	void normalize()
	{
		m_sum = 0;

		for( int i=0; i< m_bins; i++ )
			m_sum += m_hist[i];

		for( int i=0; i< m_bins; i++ )
			m_hist[i] /= m_sum;
	}

	//getMean: return mean of this histogram
	float getMean()
	{
		float sum = 0;
		float cnt = 0;
		for( int i=0; i<m_bins; i++ ){
			cnt += m_hist[i];
			sum += m_hist[i]*(i/(float)m_bins);
		}
        
        float mean = sum/cnt;
        //if( mean != mean )return 0;
		//else return sum/cnt;
        return mean;
    }

	//getMean: return vairnace of this historgram
	float getVariance()
	{
		float sum = 0;
		float sum2 = 0;
		float cnt = 0;
		for( int i=0; i<m_bins; i++ ){
			cnt += m_hist[i];
			sum += m_hist[i]*(i/(float)m_bins);
			sum2 += m_hist[i]*(i/(float)m_bins)*(i/(float)m_bins);
		}
        
        float var = sum2/cnt - (sum/cnt)*(sum/cnt);
        //if( mean != mean )return 0;
		//else return sum/cnt;
        return var;
    }

	//getMean: return vairnace of this historgram
	//input: nomalized distribution
	float getEntropy()
	{
		float entropy = 0;
		for( int i=0; i<m_bins; i++ ){
			if( m_hist[i]==0 )continue;
			if( m_hist[i] > 0.9999999 )entropy = 0;//consider if( >1.0) , numerical problem
			else entropy += (-1*m_hist[i]*(log(m_hist[i])/log(2.0)) );
		}
        
        return entropy;
    }

	//getSelfInformationBingOrder: return a bin(every element is a bin(ing) array sort from 
	//largest information bin to smallest bin
	void getSeflInformationBinOrder( int* bins, float* selfInfo )
	{
		float ttFeq=0;
		vector<SelfInfo> sis;
		SelfInfo si;
		for( int i=0; i<m_bins; i++ ){
			ttFeq += m_hist[i];
			si.bin = i;
			si.seflInfo = 0;
			sis.push_back(si);
		}
		for( int i=0; i<m_bins; i++ ){
			float p = m_hist[i]/ttFeq;
			if( p==0 )sis[i].seflInfo=0;
			else{
				sis[i].seflInfo = -log( p );
				if( sis[i].seflInfo == 0 )sis[i].seflInfo==0.000001;//make the bin with frequency, it can have self information at least
			}
		}
		std::sort( sis.begin(), sis.end(), selfInfoSort );
		for( int i=0; i<m_bins; i++ ){
			bins[i] = sis[i].bin;
			selfInfo[i] = sis[i].seflInfo;
		}
    }

	//printHistogram: use printf output histogram info.
	void printHistogram()
	{
		for( int i=0; i<m_bins; i++ )
		{
			printf("%10.6f\n", m_hist[i]);
			//printf("%f ", m_hist[i]);
            if( i== 127)break;
		}
		printf("\n");
	}
    
	//output this histogram to file
    void writeToFile( FILE* fp )
    {
        fwrite( m_hist, sizeof(float), m_bins, fp);
    }
    
	//getBinNum: get m_bins, bin number of this histogram
	//getgetget test commend
    int getBinNum()
    {
        return m_bins;
    }

	int findMaxNonZeroBin()
	{
		int b = m_bins-1;
		for( ; b>=0 && m_hist[b]==0; b-- );
		return b;
	}

	int findMinNonZeroBin()
	{
		int b = 0;
		for( ; b<m_bins && m_hist[b]==0; b++ );
		if( b==m_bins)b==-1;//all zero
		return b;
	}

	int findMinFreqNonZeroBin()
	{
		float minFreq = 1000000;
		int minB = 0;
		for( int b=0; b<m_bins ; b++ ){
			if( m_hist[b] != 0 && minFreq>m_hist[b]){
				minB = b;
				minFreq = m_hist[b];
			}
		}
		return minB;
	}

	int findMaxFreqNonZeroBin()
	{
		float maxFreq = 0;
		int maxB = 0;
		for( int b=0; b<m_bins ; b++ ){
			if( m_hist[b] != 0 && maxFreq<m_hist[b]){
				maxB = b;
				maxFreq = m_hist[b];
			}
		}
		return maxB;
	}

	//sample a data point from histogram by the probability
	float sampleDataPoint()
	{
		int samples = 0;
		for( int b=0; b<m_bins; b++ )
			samples += m_hist[b];
		float randfloat = rand()/(float)RAND_MAX;
		int smpIdx = randfloat*samples;
		samples = 0;
		int b=0;
		for( b=0; b<m_bins; b++ ){
			samples += m_hist[b];
			if( samples >= smpIdx )
				break;
		}

		return b/(float)m_bins;
	}

	//sample a data point from out of 't' sigma data
	float sampleOutlierDataPoint( float t, float dis )
	{
		int bDis = dis*(float)m_bins;
		float mean = getMean();
		float sd = sqrt(getVariance());

		int bMean = mean*(float)m_bins;
		
		int bSt;
		if( mean - t*sd < 0 ){
			bSt = 0;
		}else{
			bSt = (mean - t*sd)*(float)m_bins;
			if( bSt > bMean - bDis )
				bSt = bMean - bDis;
			if( bSt < 0 )bSt = 0;
		}

		int bEd;
		if( mean + t*sd >= 1 ){
			bEd = m_bins-1;
		}else{
			bEd = (mean + t*sd)*(float)m_bins;
			if( bEd < bMean + bDis )
				bEd = bMean + bDis;
			if( bEd >= m_bins )bEd = m_bins-1;
		}



		int samples = 0;
		for( int b=0; b<m_bins; b++ ){
			if( b>=bSt && b<=bEd )continue;
			samples += m_hist[b];
		}
		if( samples == 0 )
			return -1;//no outlier
		
		float randfloat = rand()/(float)RAND_MAX;
		int smpIdx = randfloat*samples;
		samples = 0;
		int b=0;
		for( b=0; b<m_bins; b++ ){
			if( b>=bSt && b<=bEd )continue;
			samples += m_hist[b];
			if( samples >= smpIdx )
				break;
		}
		//printf("%d %f %d %d   %d\n", bMean, sd, bSt, bEd,b );getchar();

		return b/(float)m_bins;
	}

	//sample a bunch of point from this histogram
	//outputSamples is the return float array
	//numSample is the number of samples you need
	void samplePoints( float* outputSamples, int numSample )
	{
		vector<SamplingPoint> samplingVector;
		samplingVector.clear();

		for( int i = 0; i < numSample; i ++ ){
			float randNumber = rand() / (float) RAND_MAX;
			SamplingPoint sp;
			sp.order = i;
			sp.randNumber = randNumber;
			samplingVector.push_back( sp );
		}

		std::sort( samplingVector.begin(), samplingVector.end(), samplingVectorSortByRandNumber );

		float currentCDF = 0;
		float currentIndex = 0;
		for( int i = 0; i < m_bins; i ++ ){
			currentCDF += m_hist[ i ];
			while( currentIndex < samplingVector.size() && currentCDF > samplingVector[ currentIndex ].randNumber ){
				samplingVector[ currentIndex ].sampleValue = i / (float) m_bins;
				currentIndex ++;
			}
		}

		std::sort( samplingVector.begin(), samplingVector.end(), samplingVectorSortByOrder );

		for( int i = 0; i < numSample; i ++ ){
			outputSamples[ i ] = samplingVector[ i ].sampleValue;
		}
	}

	//assign: h1 = h2
	//h1: assigned histogram
	//h2: input histogram
	static void assign( Histogram* h1, Histogram* h2 )
	{
		memcpy( h1->m_hist, h2->m_hist, sizeof( float )*h1->m_bins );
	}

	//addHistogram: hist = a + b
	static void addHistogram( Histogram* hist, Histogram* a, Histogram* b )
	{
		for( int i=0; i<hist->m_bins; i++ ){
			hist->m_hist[i] = a->m_hist[i] + b->m_hist[i]; 
		}
	}

	//subHistogram: hist = b-a
	static void subHistogram( Histogram* hist, Histogram* a, Histogram* b )
	{
		for( int i=0; i<hist->m_bins; i++ ){
			hist->m_hist[i] = a->m_hist[i] - b->m_hist[i];
		}
	}

	//subSubHistogram: subtract sub-histogram: hist = b-a
	//include: s and t, only calculate s-t bins
	//hist: return histogram: the m_bins for this histogram shoule be t-s+1
	//s: start bin
	//t: till bin
	static void subSubHistogram( Histogram* hist, Histogram* a, Histogram* b, int s, int t )
	{
		for( int i=s; i<=t; i++ ){
			hist->m_hist[i] = a->m_hist[i] - b->m_hist[i];
		}
	}

	//calculateEuclideanDistance: cal Eu distance between h1 and h2
	static float calEuclideanDistance( Histogram* h1, Histogram* h2 )
	{
		float error = 0;	//calculate error by Eu distance
		for( int j=0; j<h1->m_bins; j++ )error += (h2->m_hist[j] - h1->m_hist[j]) * (h2->m_hist[j] - h1->m_hist[j]);
		
		return sqrt(error);
	}

	//calBhattacharyyaDistance
	//the return value is different from definition of BC in Bhattacharyya Distance
	//BC is similarily, 1 in BC means same distribution
	//we return 1.0 - BC, that means error
	static float calBhattacharyyaDistance( Histogram* h1, Histogram* h2 )
	{
		float bc = 0;	//calculate error by Eu distance
		for( int j=0; j<h1->m_bins; j++ )bc += sqrt(h2->m_hist[j] * h1->m_hist[j]);
		
		if( bc > 1.0 )bc = 1.0;
		if( bc < 0.0 )bc = 0.0;

		return 1.0 - bc;
	}

	//calEMD
	//the return value is different from definition of EMD in Bhattacharyya Distance
	//we normalize the result
	//0 is the same, 1 is the biggest error
	//h1 h2 is the normalized distribution
	static float calEMD( Histogram* h1, Histogram* h2 )
	{
		float emd = 0;
		float h1cdf = 0;
		float h2cdf = 0;
		for( int j=0; j<h1->m_bins; j++ ){
			h1cdf += h1->m_hist[j];
			h2cdf += h2->m_hist[j];
			emd += fabs(h1cdf - h2cdf);
		}
		
		return emd/float(h1->m_bins);
	}

	//h1:leftup h2:rightup h3:leftdown: h4:rightdown
	//ab is left and right propotion
	//cd is up and down proportion
	static void linearInterpolation2D( Histogram* hist, Histogram* h1, Histogram* h2,Histogram* h3, Histogram* h4,
										float a, float b, float c, float d )
	{
		Histogram tmp1(hist->m_bins);
		Histogram tmp2(hist->m_bins);
		linearInterpolation1D( &tmp1, h1, h2, a, b );
		linearInterpolation1D( &tmp2, h3, h4, a, b );
		linearInterpolation1D( hist, &tmp1, &tmp2, c, d );
	}

	static void linearInterpolation1D( Histogram* hist, Histogram* inh1, Histogram* inh2, float a, float b )
	{
		//make integral by bin order(cdf histogram)
		//and normalize to 0-1 (last bin is 1 )
		Histogram* h1 = new Histogram( inh1);
		Histogram* h2 = new Histogram( inh2);
		int bins = h1->getBinNum();
		float accH1 = 0;
		float accH2 = 0;
		for( int i=0; i<bins; i++ ){
			accH1 += h1->getBin( i );
			accH2 += h2->getBin( i );
			h1->setBin( i, accH1 );
			h2->setBin( i, accH2 );
		}
		for( int i=0; i<bins; i++ ){
			float h1b = h1->getBin( i );
			float h2b = h2->getBin( i );
			h1->setBin( i, h1b/accH1 );
			h2->setBin( i, h2b/accH2 );
		}
		//h1->printHistogram();
		//h2->printHistogram();
		

		int baseBin1 = 0;
		int baseBin2 = 0;
		float binWidth = 1.0/(float)(hist->getBinNum());
		float targetProb = 0;
		int lastBin = 0;
		float lastCdf = 0;

				while( targetProb<=1.0 ){
					float x1 = getXbyCDF( h1, targetProb, baseBin1 );
					float x2 = getXbyCDF( h2, targetProb, baseBin2 );
					float x = a*x1 + b*x2;

						int nowBin = (int)(x*bins);
						for( int j=lastBin+1; j<nowBin; j++ )hist->setBin( j, lastCdf );
						hist->setBin( nowBin, targetProb );
						lastBin = nowBin;
						lastCdf = targetProb;

					targetProb += 0.01;
				}

		for( int i=lastBin+1; i<bins; i++ ){
			hist->setBin( i, lastCdf );
		}

		//get a pdf after above process
		for( int i=hist->getBinNum()-1; i>0; i-- ){
			float thisBin = hist->getBin( i );
			float prevBin = hist->getBin( i-1 );
			hist->setBin(i, thisBin-prevBin );
		}
		hist->normalize();

		//this should have the following two lines code to avoid memory leak, but some unknown problem, I can not call it, fix it later
		//h1->~Histogram();
		//h2->~Histogram();
	}

	//baseBin: assure the cdf we want within this bin
	static float getXbyCDF( Histogram* hist, float cdf, int& baseBin )
	{
		float binWidth = 1.0/(float)(hist->getBinNum());
		float thisCdf = hist->getBin( baseBin );
		while( thisCdf < cdf ){
			baseBin++;
			thisCdf = hist->getBin( baseBin );
		}
		float baseCdf=0;
		float lastBin=0;
		float baseX=0;
		if( baseBin!=0){
			baseCdf= hist->getBin( baseBin-1 );
			lastBin = hist->getBin(baseBin-1);
			baseX = (baseBin)*binWidth;
		}
		
		float slope = (hist->getBin(baseBin) - lastBin)/binWidth;
		float offsetCdf = cdf - baseCdf;

		if( slope == 0 ){
			return baseX;
		}
		else{
			return (offsetCdf / slope)+baseX;
		}
	}

	//hist is a cdf histogram (normalized to 0-1 and is cdf)
	//x is the location we want
	//return the cdf value at x (0-1)
	static float getCDFfromCDFHistogram( Histogram* hist, float x )
	{
		int bin = (int)(x*hist->getBinNum());
		int baseBin = bin -1;
		float binWidth = 1.0/(float)(hist->getBinNum());
		float slope = (hist->getBin(bin) - hist->getBin(baseBin))/binWidth;
		float offset = x - binWidth*bin;
		
		float baseCDF = hist->getBin(baseBin);

		return baseCDF + offset*slope;
	}

	


private:
	int m_bins;
	float* m_hist;	//histogram data
	float m_sum;	//sum of bins before normalize
};

#endif