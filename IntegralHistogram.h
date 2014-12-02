#ifndef _INTEGRALHISTOGRAM_H_
#define _INTEGRALHISTOGRAM_H_

#include <vector>
using namespace std;
#include "Histogram.h"

class IntegralHistogram
{
public:
	//constructor: construct a empty IH
	IntegralHistogram()
	{
		m_hists.clear();
	}

	//deconstructor
	~IntegralHistogram()
	{
		for( int i=0; i<m_hists.size(); i++ ){
			delete m_hists[i];
		}
		m_hists.clear();
	}

	//pushHistogram: given a histogram and add it into this IH
	//it means the new location for this histogram in this IH
	//is this histogram+previous_integralHistogram
	void pushHistogram( Histogram *a)
	{
		Histogram* h = new Histogram(a);
		if( m_hists.size() != 0 )
			Histogram::addHistogram( h, h, m_hists[m_hists.size()-1] );
		m_hists.push_back( h );
    }

    int getbins(){
        return m_hists[0]->getBinNum();
    }

	//getHistogram: get histogram between s and t(include s and t )
	//hist: return histogram
	void getHistogram( Histogram* hist, int s, int t )
	{
		if( s == 0 )
			Histogram::assign( hist, m_hists[t] );
		else
			Histogram::subHistogram( hist, m_hists[t], m_hists[s-1] );
	}

	//getSubBinHistogram: the same as getHistogram but only calculate bin between binS and binT
	//hist: return histogram: bin number should be binT-binS+1
	//sliceS, sliceT: get histogram between sliceS and sliceT (include S and T )
	//binS and binT: calculate bin between binS and binT (include binS and binT)
	void getSubBinHistogram( Histogram* hist, int sliceS, int sliceT, int binS, int binT )
	{
		if( sliceS == 0 ){
			float* f = (float*) malloc( sizeof(float) * (binT-binS+1) );
			m_hists[sliceT]->getBins( f, binS, binT );
			hist->setBins( binS, binT, f );
			free(f);
		}else{
			Histogram::subSubHistogram( hist, m_hists[sliceT], m_hists[sliceS-1], binS, binT );
		}
	}

	//printIntegralHistogram: print this integral histogram
	void printIntegralHistogram()
	{
		for( int i=0; i<m_hists.size(); i++ ){
			printf( "%02d : ",i );
			m_hists[i]->printHistogram();
		}
	}
    
    //write integral histogram to file(not do any substraction)
    void writeToFile(FILE* fp)
    {
        for( int i=0; i<m_hists.size(); i++ ){
			m_hists[i]->writeToFile(fp);
		}
    }
    
	//getSize: return number of histogram in this integral histogram
    int getSize()
    {
        return (int)(m_hists.size());
    }

private:
	vector<Histogram*> m_hists;
	//histograms: every histogram is integraled( m_hists[i]-m_hist[i-1] can get raw histogram at i location)
};

#endif