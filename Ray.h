#ifndef _RAY_H_
#define _RAY_H_

#include "IntegralHistogram.h"

class Ray
{
public:
	//Ray: construct a ray by given many sampls
	//raw: a 3D array store the samples
	//dim1 dim2 dim3: dimension of raw array, dim1 is depth
	//bins: bin number you want
	Ray( float*** raw, int dim1, int dim2, int dim3, int bins )
	{
		//for a ray, dim1 is depth dimension
		//dim2 dim3 form a slice
		//data domain has to be 0-1
		m_dim1 = dim1;
		m_dim2 = dim2;
		m_dim3 = dim3;
		m_bins = bins;
		m_slcSample = m_dim2*m_dim3;

		float* data = (float*)malloc( sizeof(float)*m_dim2*m_dim3 );
		m_dnsIH = new IntegralHistogram();
        m_rlsDnsIH = false;
		m_idmIH = new IntegralHistogram();		
		
		for( int i = 0; i< m_dim1; i++ ){
			int cnt = 0;
			for( int j = 0; j<m_dim2; j++ ){
				for( int k = 0; k<m_dim3; k++ ){
					data[cnt++] = raw[i][j][k];
				}
			}
			Histogram h( m_bins, m_slcSample, 0, 1, data );
			m_dnsIH->pushHistogram( &h );
		}
		free( data );
	}

		Ray( float* histogram, int step, int bins, int slcSample )
	{
		//for a ray, dim1 is depth dimension
		//dim2 dim3 form a slice
		//data domain has to be 0-1
		m_dim1 = step;
		m_bins = bins;
		m_slcSample = slcSample;

		m_dnsIH = new IntegralHistogram();
        m_rlsDnsIH = false;
		m_idmIH = new IntegralHistogram();		
		
		int cnt = 0;
		for( int i = 0; i< m_dim1; i++ ){
			Histogram h( m_bins );

			for( int i=0; i< m_bins; i++ ){
				h.setBin( i, histogram[cnt++] );
			}

			m_dnsIH->pushHistogram( &h );
		}
	}
    
	//Ray: construct a ray 
	//just construct it and use loadFromFile to fill remain data later
    Ray( int bins, int slcSample )
    {
        m_dim1 = 0;
        m_dim2 = 0;
        m_dim3 = 0;
        m_idmIH = new IntegralHistogram();
        m_bins = bins;
        m_slcSample = slcSample;
        m_rlsDnsIH = true;
    }
    
    ~Ray()
    {
		
        if( !m_rlsDnsIH )//if call releaseDnsIH before, do release again
            delete m_dnsIH;
        delete m_idmIH;
    }

	void allocIdmColor()
	{
		m_r.resize( m_idmIH->getSize() );
		m_g.resize( m_idmIH->getSize() );
		m_b.resize( m_idmIH->getSize() );
		m_a.resize( m_idmIH->getSize() );
	}

	int getIdmIHSize()
	{
		return m_idmIH->getSize();
	}

	void getIdmColor( float &r, float &g, float &b, float &a, int idx )
	{
		r = m_r[idx];
		g = m_g[idx];
		b = m_b[idx];
		a = m_a[idx];
	}

	void setIdmColor( float r, float g, float b, float a, int idx )
	{
		m_r[idx] = r;
		m_g[idx] = g;
		m_b[idx] = b;
		m_a[idx] = a;
	}

	int getSingleIdmHistogram( Histogram* hist, int idx )
	{
		m_idmIH->getHistogram(hist, idx, idx );
		return m_idmTrend[idx];
	}
    
	//releaseDnsIH: for save memory
    void releaseDnsIH()
    {
        delete m_dnsIH;
    }

	//produceMonoRepresentation: from m_dnsIH to construct increase decrease representation
	//result store in m_idmIH, m_idmIdx, m_idmTrend
	//treshold: the error tollerance, between 0-1
	//boxSize: the evluated slices size
	int produceMonoRepresentation( float threshold, int boxSize, float& rayError, float& rayCnt )
	{
		float lastCnt = 0;
		float lastError = 0;
		int segments = 0;
		int from = 0;
		int status = 0;
		int k;

		rayCnt = 0;
		rayError = 0;
		for( k=1; k<m_dim1; k++ ){
			float cnt=0;
			float maxError = 0;
			float err = estMonoError( from, k, status, boxSize, cnt, maxError );//the error here is mean squred error
			//if( err/cnt > threshold ){//Mean squared error > threshold
			if( maxError > threshold ){//any histogram in this region, its error>threshold
				//have to build the histogram here
				Histogram h( m_bins );
				m_dnsIH->getHistogram( &h, from, k-1 );
				m_idmIH->pushHistogram( &h );
				m_idmIdx.push_back( k-1 );
				m_idmTrend.push_back( status );

				from = k;
				status = 0;
				segments++;

				rayCnt += lastCnt;
				rayError += lastError;

				lastCnt = 0;
				lastError = 0;
			}else{
				lastCnt = cnt;
				lastError = err;
			}
		}

		Histogram h( m_bins );
		m_dnsIH->getHistogram( &h, from, m_dim1-1 );
		m_idmIH->pushHistogram( &h );
		m_idmIdx.push_back( m_dim1-1 );
		m_idmTrend.push_back( 1 );//not determine inc or dec, set inc(randomly)
		segments++;
		rayCnt += lastCnt;
		rayError += lastError;


		//m_idmIH->printIntegralHistogram();

		return segments;
	}

	//return the start locaiton of this idm histogram
	int getIdmStartLocation( int index )
	{
		if( index == 0 )
			return 0;
		else
			return m_idmIdx[index-1]+1;
	}

	int getIdmEndLocation( int index )
	{
		return m_idmIdx[index];
	}

	int getIdmLength( int index )
	{
		return getIdmEndLocation( index ) - getIdmStartLocation( index ) + 1;
	}

	void printRayData()
	{
		m_dnsIH->printIntegralHistogram();
	}
    
    void printIdmData()
    {
		QString str = "";
        for( int i=0; i<m_idmIdx.size(); i++ ){
            //printf("%d:   Idx: %d,   Trend: %d\n", i, m_idmIdx[i], m_idmTrend[i]);
			str += QString::number(i) + "   Idx:" + QString::number(m_idmIdx[i]) + "    Trend:" + QString::number(m_idmTrend[i] ) + "\n";
        }
		qDebug() << str;
        //m_idmIH->printIntegralHistogram();
    }
    
    //getHistByRaw: (m_dnsIH)get the histogram between s and t(included s & t)
    void getHistByRaw( Histogram* hist, int s, int t )
    {
        m_dnsIH->getHistogram(hist, s, t);
    }

    
    //getHistByIdm: get the histogram between s and t(included s & t)
    void getHistByIdm( Histogram* hist, int s, int t )
    {
        Histogram histS(m_bins);
        Histogram histT(m_bins);
        if( s!= 0 )
            getIHByIdm(&histS, s-1);
        getIHByIdm(&histT, t);
        Histogram::subHistogram(hist, &histT, &histS);
    }
    
	//find Isovalue: use binary search in m_idmIH to find a iso value location
	//return location of the isovalue in m_d domain
	//vxl: the isovalue you want(between 0-1)
	//curIdmIdx: return index of m_idmIH, can search next voxel from next index of m_idmIH
	//thick: the thickness of this isovalue: in m_d domain
	//min: in domain of index of m_idmIH, if last search return curIdmIdx is i
	// set min as curIdmIdx+1 , then this search will start at the location of last result
    int findIsovalue( float vxl, int &curIdmIdx, int &thick, int min = -1 )
    {
        //if 3 is the largest one you donot need, min must be 3
		if( min >= m_idmIH->getSize() ){
			return -1;
		}

        int bin = (int)(vxl*m_bins);
        int max = m_idmIH->getSize()-1;
        int mid = (int)((min+max)/2.0);
        Histogram hist(m_bins);
        m_idmIH->getSubBinHistogram(&hist, min+1, max, bin, bin);
        if( hist.getBin(bin) == 0){
            thick = 0;
            curIdmIdx = -1;
            return -1;//no this voxel
        }
        
        do
        {
            m_idmIH->getSubBinHistogram(&hist, min+1, mid, bin, bin);
			bool hasVxl = false;
            if( hist.getBin(bin) != 0)hasVxl = true;
            if( hasVxl ){
                max = mid;
            }else{
                min = mid;
            }
            mid = (int)((min+max)/2.0);
        }while(max-min > 1);
        
        m_idmIH->getHistogram(&hist, max, max);
        int baseIdx;
        if( max == 0 )baseIdx = 0;
        else baseIdx = m_idmIdx[max-1];
        float belowSamples = 0;
        if( m_idmTrend[max] == 1 ){
            for( int i=0; i< bin; i++ )
                belowSamples += hist.getBin(i);
        }else{//== -1 (decrease)
            for( int i=m_bins-1; i>bin; i-- )
                belowSamples += hist.getBin(i);
        }
        int offsetIdx = belowSamples/(float)m_slcSample;
        curIdmIdx = max;
        
        return baseIdx+offsetIdx;
    }

	//find the isovalue from the end of the ray
	//this one is different findIsovalue, it can only find the first isovalue from the end of the ray
	int findLastIsovalue( float vxl )
    {
		int min = -1;
        int bin = (int)(vxl*m_bins);
        int max = m_idmIH->getSize()-1;
        int mid = (int)((min+max)/2.0);
        Histogram hist(m_bins);
        m_idmIH->getSubBinHistogram(&hist, min+1, max, bin, bin);
        if( hist.getBin(bin) == 0){
            return -1;//no this voxel
        }
        
        do
        {
            m_idmIH->getSubBinHistogram(&hist, mid+1, max, bin, bin);
			bool hasVxl = false;
            if( hist.getBin(bin) != 0)hasVxl = true;
            if( hasVxl ){
                min = mid;
            }else{
                max = mid;
            }
            mid = (int)((min+max)/2.0);
        }while(max-min > 1);
        
        m_idmIH->getHistogram(&hist, max, max);

        int baseIdx;
        if( max == 0 )baseIdx = 0;
        else baseIdx = m_idmIdx[max-1];
        float belowSamples = 0;
        if( m_idmTrend[max] == 1 ){
            for( int i=0; i< bin; i++ )
                belowSamples += hist.getBin(i);
        }else{//== -1 (decrease)
            for( int i=m_bins-1; i>bin; i-- )
                belowSamples += hist.getBin(i);
        }
        int offsetIdx = belowSamples/(float)m_slcSample;
        
        return baseIdx+offsetIdx;
    }

	//find Isovalue: use binary search in m_idmIH to find a iso value location
	//return location of the isovalue in m_d domain
	//vxl: the isovalue you want(between 0-1)
	//curIdmIdx: return index of m_idmIH, can search next voxel from next index of m_idmIH
	//thick: the thickness of this isovalue: in m_d domain
	//min: in domain of index of m_idmIH, if last search return curIdmIdx is i
	// set min as curIdmIdx+1 , then this search will start at the location of last result
    int findIsovalueRange( float vxl1, float vxl2, int &curIdmIdx, int &thick, int min = -1 )
    {
        //if 3 is the largest one you donot need, min must be 3
		if( min >= m_idmIH->getSize() ){
			return -1;
		}

        int bin1 = (int)(vxl1*m_bins);
		int bin2 = (int)(vxl2*m_bins);
        int max = m_idmIH->getSize()-1;
        int mid = (int)((min+max)/2.0);
        Histogram hist(m_bins);
        m_idmIH->getSubBinHistogram(&hist, min+1, max, bin1, bin2);
		bool hasVxl = false;
		for( int i=bin1; i<=bin2; i++ ){
			if( hist.getBin(i) != 0 ){
				hasVxl = true;
				break;
			}
		}
        if( !hasVxl ){
            thick = 0;
            curIdmIdx = -1;
            return -1;//no this voxel
        }
        
        do
        {
            m_idmIH->getSubBinHistogram(&hist, min+1, mid, bin1, bin2);
			bool hasVxl = false;
			for( int i=bin1; i<=bin2; i++ ){
				if( hist.getBin(i) != 0 ){
					hasVxl = true;
					break;
				}
			}
            if( hasVxl ){
                max = mid;
            }else{
                min = mid;
            }
            mid = (int)((min+max)/2.0);
        }while(max-min > 1);
        
        curIdmIdx = max;
        
        return 1;
    }

	//finSingleMaterial: give transfer funtion and starting search point find the first material block and return it using IDM
	//material: transfer 
	int findSingleMaterialByIdm( float* material, int numMaterial, int depth, int min = -1 )
    {
		float singleMaterialThreshold = 0.95;

		float* vote = (float*)malloc( sizeof(float) * numMaterial );
		//vote[i] indicate material[i-1]<data<=matrieal[i]
		int* voteIdx = (int*)malloc( sizeof(int) * m_bins );
		int vi = 0;
		for( int i=0; i<m_bins; i++ ){
			float vxl = i/(float)m_bins;

			int j;
			for( j=vi; vxl > material[j]; j++ );
			vi = j;
			voteIdx[i] = vi;
		}


        //if 3 is the largest one you donot need, min must be 3
		if( min >= depth ){
			return -1;
		}

        int max = depth-1;
        int mid = (int)((min+max)/2.0);
        Histogram hist(m_bins);
		//m_idmIH->getSubBinHistogram(&hist, min+1, max, 0, m_bins);
        int initMin = min;
        do
        {
            //m_dnsIH->getSubBinHistogram(&hist, min+1, mid, 0, m_bins);	
			//m_dnsIH->getSubBinHistogram(&hist, initMin+1, mid, 0, m_bins);
			getHistByIdm(&hist, initMin+1, mid );
			
			memset( vote, 0, sizeof(float)*numMaterial );
			for( int i=0; i<m_bins; i++ ){
				vote[ voteIdx[i] ] += hist.getBin( i );
			}
			int sumVote = 0;
			int largeVote = 0;
			int largeVoteIdx = 0;
			for( int i=0; i<numMaterial; i++ ){
				sumVote += vote[i];
				if( vote[i] > largeVote ){
					largeVote = vote[i];
					largeVoteIdx = i;
				}
			}
			
			bool singleMaterial = false;
            if( largeVote/(float)sumVote > singleMaterialThreshold )singleMaterial = true;

            if( singleMaterial ){
                min = mid;
            }else{
                max = mid;
            }
            mid = (int)((min+max)/2.0);
			
        }while(max-min > 1);
        
		free( vote );
		free( voteIdx );

        return max;
    }

	//finSingleMaterial: give transfer funtion and starting search point find the first material block and return it
	//material: transfer 
	int findSingleMaterial( float* material, int numMaterial, int min = -1 )
    {
		float singleMaterialThreshold = 0.95;

		float* vote = (float*)malloc( sizeof(float) * numMaterial );
		//vote[i] indicate material[i-1]<data<=matrieal[i]
		int* voteIdx = (int*)malloc( sizeof(int) * m_bins );
		int vi = 0;
		for( int i=0; i<m_bins; i++ ){
			float vxl = i/(float)m_bins;

			int j;
			for( j=vi; vxl > material[j]; j++ );
			vi = j;
			voteIdx[i] = vi;
		}


        //if 3 is the largest one you donot need, min must be 3
		if( min >= m_dnsIH->getSize() ){
			return -1;
		}

        int max = m_dnsIH->getSize()-1;
        int mid = (int)((min+max)/2.0);
        Histogram hist(m_bins);
		//m_idmIH->getSubBinHistogram(&hist, min+1, max, 0, m_bins);
        int initMin = min;
        do
        {
            //m_dnsIH->getSubBinHistogram(&hist, min+1, mid, 0, m_bins);
			m_dnsIH->getSubBinHistogram(&hist, initMin+1, mid, 0, m_bins);
			
			memset( vote, 0, sizeof(float)*numMaterial );
			for( int i=0; i<m_bins; i++ ){
				vote[ voteIdx[i] ] += hist.getBin( i );
			}
			int sumVote = 0;
			int largeVote = 0;
			int largeVoteIdx = 0;
			for( int i=0; i<numMaterial; i++ ){
				sumVote += vote[i];
				if( vote[i] > largeVote ){
					largeVote = vote[i];
					largeVoteIdx = i;
				}
			}
			
			bool singleMaterial = false;
            if( largeVote/(float)sumVote > singleMaterialThreshold )singleMaterial = true;

            if( singleMaterial ){
                min = mid;
            }else{
                max = mid;
            }
            mid = (int)((min+max)/2.0);
			
        }while(max-min > 1);
        
		free( vote );
		free( voteIdx );

        return max;
    }
    
    //test only
    void testRawIHandIDMonly( int idx )
    {
        printf("Groud truth: %d\n", idx);
        Histogram hist(m_bins);
        m_dnsIH->getHistogram(&hist, 0, idx);
        hist.printHistogram();
        
        printf("IDM result:\n");
        getIHByIdm(&hist, idx);
        hist.printHistogram();
    }
    
	//loadFromFile: load ray(histogram) from file
    void loadFromFile( FILE* fp )
    {
        float* tmp = new float[m_bins];
        fread(tmp, sizeof(float), m_bins, fp);
        Histogram *hist = new Histogram(m_bins, tmp);
        m_idmIH->pushHistogram(hist);
        //Histogram hist(m_bins, tmp);
        //m_idmIH->pushHistogram(&hist);
        
        //tep[0] is index, tmp[1] is trend
        fread(tmp, sizeof(float), 2, fp);
        m_idmIdx.push_back((int)tmp[0]);
        m_idmTrend.push_back((int)tmp[1]);
        delete tmp;
       // delete hist;
    }
    
	//writeToFile: write histogram(ray) to file
    void writeToFile( FILE* fp )
    {
        for( int i=0; i<m_idmIdx.size(); i++ ){
            Histogram hist(m_bins);
            m_idmIH->getHistogram(&hist, i, i);
            hist.writeToFile(fp);
            float idmIdx = m_idmIdx[i];
            fwrite( &idmIdx, sizeof(float), 1, fp);
            float idmTrend = m_idmTrend[i];
            fwrite( &idmTrend, sizeof(float), 1, fp);
        }
    }

	int getIdmFlag( int idx  )
	{
		int i=0;
		int flag = 0;
		for( i=0; i<m_idmIdx.size(); i++ ){
			if( idx <=i ){
				flag = m_idmTrend[i];
				break;
			}
		}
		return flag;
	}

private:
	//estMonoError: given m_dnsIH, from-to location, status of this section(inc or dec) and 
	//step(size of the evaluated block)
	//return the root mean square error
	float estMonoError( int from, int to, int &status, int step, float& returnCnt, float& maxError )
	{
		maxError = 0;
		float e=0;//return error
		float cnt=0;
		Histogram hist(m_bins); //histogram between from and to by ihRay
		m_dnsIH->getHistogram( &hist, from, to );

		if( status == 0 ){//mean shift to determine increase or decrease
			Histogram phist(m_bins); //previous histogram, no include 'to'
			//m_dnsIH->getHistogram( &phist, from, to-1 );
			m_dnsIH->getHistogram( &phist, from, from );

			Histogram thist(m_bins);
			m_dnsIH->getHistogram( &thist, to, to );

			float mean = thist.getMean();
			float pmean = phist.getMean();
			//printf("%f %f   %d  %d\n", mean, pmean, from , to);
			if( fabs(mean- pmean ) < 0.0001 )status = 0;//if mean and pmean very close, cannot detemine increase or decrease
			else if( mean > pmean )status = 1;
			else if( mean < pmean ) status = -1;

			e = 0;
			cnt = to - (from+(step-1)) +1;
		}
		else{//increase or decrease by status
			//printf("inc or dec\n");
			Histogram h(m_bins);//h: store ground truth
			Histogram mh(m_bins);//mh: store histogram from monoIncDec model
			for( int i=from+(step-1); i<=to; i++ ){
				m_dnsIH->getHistogram( &h, i-step+1, i );//get histogram at i-step+1~i slice , ground truth
				m_dnsIH->getHistogram( &mh, from, to );
				getHistByMonoAssump( &mh, i-from-(step-1), i-from, to-from+1, status );

				h.normalize();
				mh.normalize();

				//float error = Histogram::calEuclideanDistance( &mh, &h );	//calculate error by Eu distance
				//float error = Histogram::calBhattacharyyaDistance( &mh, &h );	//calculate error by 1.0- BC in BD
				float error = Histogram::calEMD( &mh, &h );	//calculate error by EMD (between 0-1)
				e+=(error*error);	//just sqaure sum now
				cnt++;
				if( error > maxError ){
					maxError = error;
				}
			}			
		}

		returnCnt = cnt;
		return e;	//just return saure sum now
	}

	//getHistByMonoAssump: given a histogram(hist)
	//and location of this histogram(from-to), 0 is smallest index
	//slices is the number of slices of this hitogram
	//mono is increase or decrease
	//return the histogram of the from-to region in "hist"
	void getHistByMonoAssump( Histogram* hist, int from , int to, int slices, int mono )
	{
		float cutBelow, cutUpper;
		if( mono == 1 ){
			cutBelow = from* m_slcSample;
			cutUpper = (slices - to - 1)*m_slcSample;
		}else{//mono ==-1
			cutUpper = from* m_slcSample;
			cutBelow = (slices - to - 1)*m_slcSample;
		}

		for( int i=0; i<m_bins && cutBelow>0; i++ ){
			float tmp = hist->getBin(i) - cutBelow;
			if( tmp < 0 ){
				hist->setBin(i, 0 );
				cutBelow = -1*tmp;
			}
			else{
				hist->setBin(i, tmp );
				cutBelow = 0;
			}
		}

		for( int i=m_bins-1; i>=0 && cutUpper>0; i-- ){
			float tmp = hist->getBin(i) - cutUpper;
			if( tmp < 0 ){
				hist->setBin(i, 0);
				cutUpper = -1*tmp;
			}
			else{
				hist->setBin(i, tmp);
				cutUpper = 0;
			}
		}
	}

    //getIHByIdm: get integral histogram to idx ( IH of 0-idx of the ray)
	void getIHByIdm( Histogram* hist, int idx )
	{
		Histogram pHist1(m_bins); //part histogram
		Histogram pHist2(m_bins);//complete histogram, lower
		int highIdx; //index of m_idmIdx
        int highIHIdx; //index of m_idxIH
		int lowIHIdx; //index of m_idxIH
        int trend;//inc or dec
        
		getHighIndex( highIdx, idx );//highIdx is the index in m_idmIH
        //if only use the first segment
		if( highIdx == 0 ){
			highIHIdx = m_idmIdx[highIdx];
            lowIHIdx = 0;
            m_idmIH->getHistogram( &pHist1, highIdx, highIdx );
            trend = m_idmTrend[highIdx];
            //make sure pHist2 is all 0, after construct it, it is all 0
		}else{
            highIHIdx = m_idmIdx[highIdx];
            lowIHIdx = m_idmIdx[highIdx-1]+1;
			m_idmIH->getHistogram( &pHist1, highIdx, highIdx );
            trend = m_idmTrend[highIdx];
            m_idmIH->getHistogram( &pHist2, 0, highIdx-1 );
		}
		getHistByMonoAssump( &pHist1, 0 , idx - lowIHIdx, highIHIdx - lowIHIdx + 1, trend );
        Histogram::addHistogram(hist, &pHist1, &pHist2);
	}


	//getHighIndex: used by getIHByIdm function
	void getHighIndex( int &h, int idx )
	{
		h = -1;
		//given idx, find and return the two indcs bound idx in idmIdx
		if( idx <= m_idmIdx[0] || m_idmIdx.size() == 1){
			h = 0;
			return;
		}else{
			int imax = m_idmIdx.size()-1;
			int imin = 0;
			while (imax >= imin)
			{   
				// calculate the midpoint for roughly equal partition
				int imid = floor( ( imax + imin )/2.0 );
				if(m_idmIdx[imid] < idx && idx<=m_idmIdx[imid+1]){
					h = imid+1;
					return;
				}
				else if (m_idmIdx[imid] < idx)
					imin = imid + 1;
				else         
					imax = imid;
			}
		}
	}

	IntegralHistogram* m_dnsIH;//raw integral histogram, almost ground truth
	IntegralHistogram* m_idmIH;//ih of increase decrease model 
	vector<int> m_idmIdx;	//idm model idx, record the localtion of every histogram
	vector<int> m_idmTrend;	//idm model, increase or decrease//1: increase,  -1:decrease
	vector<float> m_r;
	vector<float> m_g;
	vector<float> m_b;
	vector<float> m_a;
	int m_dim1, m_dim2, m_dim3;		//dimension of this data
	//*******m_dim1 is the depth of this ray****(importatnt)

	int m_bins;	//sin number
	int m_slcSample;//samples of a slice
    bool m_rlsDnsIH;//if true, m_dnsIH in memory
};

#endif 