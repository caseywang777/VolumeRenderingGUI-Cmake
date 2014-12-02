//
//  Image.h
//  TestRayVolumeRendering
//
//  Created by Casey on 2013/12/25.
//  Copyright (c) 2013年 Casey. All rights reserved.
//

#ifndef _IMAGEBASED_H_
#define _IMAGEBASED_H_

#include <omp.h>
#include "RegularScalarData.h"
#include "myUtility.h"
#include "math.h"
#include <vector>
#include "time.h"
#include "stdlib.h"
#include "omp.h"
#include <algorithm>

using namespace std;

typedef struct _Segment{
	int pxl;
	int from;
	int len;
	float metric;//entropy or variance ...
	float a;
	float b;
}Segment;

bool segmentSort( Segment a, Segment b )
{
	return (a.metric<b.metric);
}

bool segmentSortBigToSmall( Segment a, Segment b )
{
	return (a.metric>b.metric);
}

bool segmentSortByFrom( Segment a, Segment b )
{
	return (a.from<b.from);
}

class Image
{
public:
    //u and v is the image dimension (make sure what you want)
    Image( RegularScalarData* rsd, int u, int v )
    {
        m_rsd = rsd;
        m_u = u;
        m_v = v;
        m_color = (unsigned char*) malloc( m_u*m_v*3*sizeof(unsigned char) );
		m_color2 = (unsigned char*) malloc( m_u*m_v*3*sizeof(unsigned char) );
		m_color3 = (unsigned char*) malloc( m_u*m_v*3*sizeof(unsigned char) );
        memset(m_color, 0, m_u*m_v*4*sizeof(unsigned char) );
		memset(m_color2, 0, m_u*m_v*4*sizeof(unsigned char) );
		memset(m_color3, 0, m_u*m_v*4*sizeof(unsigned char) );
		m_depth = (float*) malloc( m_u*m_v*sizeof(float) );
		m_mrrV.resize(m_u*m_v);
		m_mrrE.resize(m_u*m_v);
		m_mrrL.resize(m_u*m_v);
    }
    

	//Image: input a preprocessing ray based file
	//pathName: path and name of the file
	//u,v : image dimension
	//d depth
	//bins: bin number
	//slcSample: how many sample per slice of a ray
	//rsd: is raw data
    Image( char* pathName, int u, int v, int d, int bins, int slcSample,  RegularScalarData* rsd=NULL )
    {
		if( rsd !=NULL ){
			m_rsd = rsd;
		}

        m_u = u;
        m_v = v;
        m_d = d;
		m_slcSample = slcSample;
        m_color = (unsigned char*) malloc( m_u*m_v*4*sizeof(unsigned char) );
		m_colorF = (float*) malloc( m_u*m_v*4*sizeof(float) );
		m_color2 = (unsigned char*) malloc( m_u*m_v*4*sizeof(unsigned char) );
		m_color2F = (float*) malloc( m_u*m_v*4*sizeof(float) );
		m_color3 = (unsigned char*) malloc( m_u*m_v*4*sizeof(unsigned char) );
        memset(m_color, 0, m_u*m_v*4*sizeof(unsigned char) );
		memset(m_colorF, 0, m_u*m_v*4*sizeof(float) );
		memset(m_color2, 0, m_u*m_v*4*sizeof(unsigned char) );
		memset(m_color3, 0, m_u*m_v*4*sizeof(unsigned char) );
		m_depth = (float*) malloc( m_u*m_v*sizeof(float) );
		m_mrrV.resize(m_u*m_v);
		m_mrrE.resize(m_u*m_v);
		m_mrrL.resize(m_u*m_v);
        m_bins = bins;
		m_ddaX = (int*) malloc( sizeof(int)*m_d );
		m_ddaY = (int*) malloc( sizeof(int)*m_d );
		m_ddaCnt = 0;

		//initialize entropy vector of every pixel
		entropyQueue = new vector<Segment>*[m_u*m_v];
		for( int i =0 ; i < m_u*m_v ; i++ )
			entropyQueue[i] = new std::vector<Segment>;
        
        FILE* fp = fopen( pathName, "rb");
        float check=0;
        fread(&check, sizeof(float), 1, fp);
        Ray* ray = new Ray( bins, slcSample );
		int i=0;
		size_t result;
        do{
			
            if(check >=0 ){//same ray
                fseek(fp, -1*sizeof(float), SEEK_CUR);
                ray->loadFromFile(fp);
            }else if( check == -1 ){//next ray
                m_rays.push_back(ray);
                ray = new Ray( bins, slcSample );
            }
            result = fread(&check, sizeof(float), 1, fp);
        }while( check!=-2 && result==1 );
        m_rays.push_back(ray);

        fclose( fp );
    }

	//****************procude IDM from different axis*******************************************
    
    //produceRayByZAsix:xy paralle to image plan
    //usf: upsampling factor
    //bins: number of bins
	//uncomplted function
    void produceRaysByZAxis(  float usf, int bins, float error, float block  )
    {
    	float evMSError = 0;
		float evCount = 0;
		float nonZeroBins = 0;
		Histogram histDns(bins);
		Histogram histIdm(bins);
        int allSegs = 0;
        int cnt=0;
        float*** data = alloc3DMatrix( m_rsd->getDim3()*usf, usf, usf);

        for( int d2 = 0; d2< m_rsd->getDim2()-1; d2+=1 ){
            for( int d1 = 0; d1< m_rsd->getDim1()-1; d1+=1 ){
                //prepare upsampling ray
                for( int k=0; k<(m_rsd->getDim3()-1)*usf; k++ ){//depth
                    for( int i=0; i<usf; i++ ){
                        for( int j=0; j<usf; j++ ){
                            data[k][i][j] = m_rsd->getUpsamplingVoxelValue(d1*usf+j, d2*usf+i, ((m_rsd->getDim3()-1)*((int)usf) -1) - k, usf);
                        }
                    }
                }
                Ray* r = new Ray( data, m_rsd->getDim3()*usf, usf, usf, bins );
                m_rays.push_back(r);
				float rayError = 0;
				float rayCnt = 0;
                int seg = r->produceMonoRepresentation( error, block, rayError, rayCnt );
                allSegs += seg;
				evCount += rayCnt;
				evMSError += rayError;

                r->releaseDnsIH();
				printf("u:%d v:%d allSegments:%d percetageOfSegs: %f  MSError: %10.9f\n", d2, d1, allSegs, allSegs/((d2*(m_rsd->getDim1()-1))+(d1+1))/((m_rsd->getDim3()-1)*usf), evMSError/evCount);
				//get nonZeroBins
				for( int i=0; i<r->getIdmIHSize(); i++ ){
					r->getSingleIdmHistogram( &histIdm, i );
					float nonZero = 0;
					for( int j=0; j<bins; j++ ){
						if( histIdm.getBin( j ) >=0.5 )nonZero++;
					}
					nonZeroBins += nonZero;
				}
				//end nonZeroBins
            }
        }

		printf("MSError : %f/%f= %10.9f      nonZeroBins: %f\n", evMSError, evCount, evMSError/evCount, nonZeroBins);
		printf("Done!\n");
    }
    
    //produceRayByXAsix:yz paralle to image plan
    //usf: upsampling factor
    //bins: number of bins
	//error: error tollerance: between 0-1
	//block: evaluated block size
	//best version right now
    void produceRaysByXAxis( float usf, int bins, float error, float block )
    {
		float evMSError = 0;
		float evCount = 0;
		float nonZeroBins = 0;
		Histogram histDns(bins);
		Histogram histIdm(bins);
        int allSegs = 0;
        int cnt=0;

		int uPixelSamples = ((m_rsd->getDim3()-1)*usf)/m_u;
		int vPixelSamples = ((m_rsd->getDim2()-1)*usf)/m_v;
		printf("UV-PixeSamples: %d %d  %f %f %d, %d\n", vPixelSamples, uPixelSamples, m_rsd->getDim1()*usf, usf, m_v, m_rsd->getDim1());
		
		float*** data = alloc3DMatrix( (int)(m_rsd->getDim1()*usf), vPixelSamples, uPixelSamples);
		
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
                //prepare upsampling ray
                for( int k=0; k<(m_rsd->getDim1()-1)*usf; k++ ){//depth
                    for( int i=0; i<vPixelSamples; i++ ){
                        for( int j=0; j<uPixelSamples; j++ ){
                            data[k][i][j] = m_rsd->getUpsamplingVoxelValue(k, v*vPixelSamples+i, u*uPixelSamples+j, usf);
                        }
                    }
                }
                Ray* r = new Ray( data, m_rsd->getDim1()*usf, vPixelSamples, uPixelSamples, bins );
                m_rays.push_back(r);
				float rayError = 0;
				float rayCnt = 0;
                int seg = r->produceMonoRepresentation( error, block, rayError, rayCnt );
                allSegs += seg;
				evCount += rayCnt;
				evMSError += rayError;

                r->releaseDnsIH();
				printf("u:%d v:%d allSegments:%d percetageOfSegs: %f  MSError: %10.9f\n", v, u, allSegs, allSegs/(float)((v*m_u+u)*((m_rsd->getDim1()-1)*usf)), evMSError/evCount);
				//get nonZeroBins
				for( int i=0; i<r->getIdmIHSize(); i++ ){
					r->getSingleIdmHistogram( &histIdm, i );
					float nonZero = 0;
					for( int j=0; j<bins; j++ ){
						if( histIdm.getBin( j ) >=0.5 )nonZero++;
					}
					nonZeroBins += nonZero;
				}
				//end nonZeroBins
            }
        }

		printf("MSError : %f/%f= %10.9f      nonZeroBins: %f\n", evMSError, evCount, evMSError/evCount, nonZeroBins);
		printf("Done!\n");
    }
	//****************procude IDM from different axis -END*******************************************
    

	//***********************DIRECT VOLUME REDNERING************************************

	//colorRendering: render color from m_rays and store to m_color
	//transfer function is from m_tf, so call makeTrasferfunction in advance is necesary
	//step: a voxel block(slices)
    void colorRenderingBySegmentMean( int step, float clrAdapt, float bgColor, float startDepth, float order, float lighting )
    {
        int rayCnt=0;
        int cnt=0;

		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram befHist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
				Ray* ray = m_rays[v*m_u+u];
                float cr = 0, cg = 0, cb = 0, ca = 0;

				
				for( int d = 0; d< m_d-step; d+=(step/8.0) ){
                    //tranditional implementation 
					ray->getHistByIdm(&hist, d, d+step-1);
                    float vxl = hist.getMean();

                    float r, g, b, a;
                    getRGBAbyTF(vxl, r, g, b, a, 1);
					//if( vxl < 0.3 ){
					//	//if( hist.getBin( 6 ) == 0 ){
					//	if( vxl < 0.0468 || vxl > 0.054 ){
					//		a = 0;
					//	}else{
					//		getRGBAbyTF(vxl, r, g, b, a, 1);
					//		a = 1.0;
					//	}
					//}else{
					//	a = 0.5;
					//}
					a*=clrAdapt;

					for( int i=0; i< step/8.0; i++ ){
						cr += (r*a) * (1.0 - ca);
						cg += (g*a) * (1.0 - ca);
						cb += (b*a) * (1.0 - ca);
						ca += a * (1.0 - ca);
					}

                    if( ca > 0.95 )break;
                }

				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				m_color[(v*m_u+u)*4+0] =  optR;
                m_color[(v*m_u+u)*4+1] =  optG;
                m_color[(v*m_u+u)*4+2] =  optB;
				m_color[(v*m_u+u)*4+3] =  255;//A
            }
        }
    }

	void colorRenderingByHistogramInfomation( )
    {
		float bgColor = m_sldBar8/(float)100.0;
		float colorModulation = m_lineEdit0;
		int step = m_lineEdit3;

		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram befHist(m_bins);
			Histogram rHist(m_bins);//use 16 too cluster visualization 
			Histogram gHist(m_bins);
			Histogram bHist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
				Ray* ray = m_rays[v*m_u+u];
                float cr = 0, cg = 0, cb = 0, ca = 0;
				float ucr = 0, ucg = 0, ucb = 0, uca = 0;	//uncertaiinty visualization color

				for( int d = 0; d< m_d-step; d+=step ){
					float frontSegTransparency, correctOpacity;

					if( d>0 )
						ray->getHistByIdm(&befHist, 0, d-1);
					else
						befHist.reset();
					ray->getHistByIdm(&hist, d, d+step-1);

					float rr, gg, bb;
					calculateSegmentColorContribution( &befHist, &hist, colorModulation, rr, gg, bb, frontSegTransparency, correctOpacity );

					cr += rr;
					cg += gg;
					cb += bb;
					ca = ( 1 - frontSegTransparency ) + correctOpacity * frontSegTransparency; //not necessary

					hist.reset();
					ray->getHistByIdm( &hist, d, d + step );

					//get weight histogram
					float sumWeight = 0;
					for( int i = 0; i < m_bins; i ++ ){
						float r, g, b, a;
						getRGBAbyTF( i / (float)m_bins, r, g, b, a, 1);
						a *= colorModulation;
						float w = hist.getBin( i ) * a;
						sumWeight += w;
						hist.setBin( i, w );
					}
					if( sumWeight< 0.001 ) continue;
					hist.normalize();

					//calculate weighted mean
					float rmean = 0, gmean = 0, bmean = 0;
					for( int i = 0; i < m_bins; i ++ ){
						float r, g, b, a;
						getRGBAbyTF( i / (float)m_bins, r, g, b, a, 1);
						float w = hist.getBin( i );
						rmean += r * w;
						gmean += g * w;
						bmean += b * w;
					}

					//calculate weighted variance
					float rvariance = 0, gvariance = 0, bvariance = 0;
					for( int i = 0; i < m_bins; i ++ ){
						float r, g, b, a;
						getRGBAbyTF( i / (float)m_bins, r, g, b, a, 1);
						float w = hist.getBin( i );
						rvariance += ( r - rmean ) * ( r - rmean ) * w;
						gvariance += ( g - gmean ) * ( g - gmean ) * w;
						bvariance += ( b - bmean ) * ( b - bmean ) * w;
					}
					rvariance = sqrt( rvariance / (float)m_bins );
					gvariance = sqrt( gvariance / (float)m_bins );
					bvariance = sqrt( bvariance / (float)m_bins );

					//pick up max varaince
					float variance = sqrt( ( ( rvariance * rvariance ) + ( gvariance * gvariance ) + ( bvariance * bvariance ) ) / 3.0 );
					
					////uncertainty color
					variance = ( variance )/ 0.05;	//rescale to certain value, otherwise the varaince is too small
					variance = pow( (double) variance, 1/1.0 );//also rescale
					float r, g, b, a;
					getRGBAbyTF( variance, r, g, b, a, 1);
					ucr += r * correctOpacity * frontSegTransparency;
					ucg += g * correctOpacity * frontSegTransparency;
					ucb += b * correctOpacity * frontSegTransparency;
                }
				
				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				m_color[(v*m_u+u)*4+0] =  optR;
                m_color[(v*m_u+u)*4+1] =  optG;
                m_color[(v*m_u+u)*4+2] =  optB;
				m_color[(v*m_u+u)*4+3] =  255;//A

				//uncertainty visualization
				unsigned char optUR = (unsigned char)(int)floor(ucr * 255.0);
                unsigned char optUG = (unsigned char)(int)floor(ucg * 255.0);
                unsigned char optUB = (unsigned char)(int)floor(ucb * 255.0);
				m_color2[(v*m_u+u)*4+0] =  optUR;
                m_color2[(v*m_u+u)*4+1] =  optUG;
                m_color2[(v*m_u+u)*4+2] =  optUB;
				m_color2[(v*m_u+u)*4+3] =  255;//A
            }
        }
    }

	void renderRenderingWorkload( )
    {
		//start to count time
		clock_t start_time, end_time;
		float total_time = 0;
		start_time = clock(); /* mircosecond */

		float colorModulation = m_lineEdit0;
		int step = m_lineEdit3;

		float transparantThreshold = 0.000001;
		float earlyTerminationThreshold = 0.9;

		int* renderFrom;	//every ray render from
		int* renderTo;		//every ray render until
		int* dispatchID;	//render job of a ray belong to whick core
		float sumWorkload = 0;	//sum of ( renderFrom - renderTo ) of all rays

		renderFrom = (int*) malloc( m_u * m_v * sizeof( int ) );
		renderTo = (int*) malloc( m_u * m_v * sizeof( int ) );
		dispatchID = (int*) malloc( m_u * m_v * sizeof( int ) );

		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);

            for( int u = 0; u< m_u; u+=1 ){
				int rayIdx = v*m_u+u;

				Ray* ray = m_rays[ rayIdx ];
                float cr = 0, cg = 0, cb = 0, ca = 0;
				int startRenderPosition;	//before this position, no visible data
				int endRenderPosition;		//rendering tail
				int endRenderLeapingPosition;	//after this position, no visible data
				int earlyTerminationPosition;

				//find startRenderPosition
				ray->getHistByIdm( &hist, 0, m_d - 1 );
				if( calculateSegmentOpacity( &hist, colorModulation ) < transparantThreshold ){
					startRenderPosition = m_d - 1;	//whole ray is transparancy
				}else{
					int min = 0;
					int max = m_d - 1;
					int mid = (int)((min+max)/2.0);
					int lastOpaquePosition = m_d - 1;
					bool flag = true; //true means it has opacity
			
					//bineary search -like algorithm
					//find the position, before it no visible data
					do
					{
						ray->getHistByIdm( &hist, min, mid );
						float opacity = calculateSegmentOpacity( &hist, colorModulation );
						if( opacity < transparantThreshold ){//min to mid is totally transparancy
							min = mid;
							flag = false;
						}else{
							lastOpaquePosition = mid;
							max = mid;
							flag = true;
						}
						mid = ( int )( ( min + max ) / 2.0 );
						if( max - min < step / 2.0 )break;
					}while( lastOpaquePosition - min > step || flag == true );
					startRenderPosition = min;
				}

				//find tail space leaping
				ray->getHistByIdm( &hist, 0, m_d - 1 );
				if( calculateSegmentOpacity( &hist, colorModulation ) < transparantThreshold ){
					endRenderLeapingPosition = 0;	//whole ray is transparancy
				}else{
					int min = 0;
					int max = m_d - 1;
					int mid = (int)((min+max)/2.0);
					int lastOpaquePosition = m_d - 1;
					bool flag = true; //true means it has opacity
			
					//bineary search -like algorithm
					//find the position, after it no visible data
					do
					{
						ray->getHistByIdm( &hist, mid, max );
						float opacity = calculateSegmentOpacity( &hist, colorModulation );
						if( opacity < transparantThreshold ){//min to mid is totally transparancy
							max = mid;
							flag = false;
						}else{
							lastOpaquePosition = mid;
							min = mid;
							flag = true;
						}
						mid = ( int )( ( min + max ) / 2.0 );
						if( max - min < step / 2.0 )break;
					}while( max - lastOpaquePosition > step || flag == true );
					endRenderLeapingPosition = max;
				}

				//find lastRenderPosition - early termination point
				ray->getHistByIdm( &hist, 0, m_d - 1 );
				if( calculateSegmentOpacity( &hist, colorModulation ) < earlyTerminationThreshold ){
					earlyTerminationPosition = m_d - 1;	//whole ray cannot reach 0.95 opacity
				}else{
					int min = 0;
					int max = m_d - 1;
					int mid = (int)((min+max)/2.0);
					int lastNotTerminatePoint = 0;
					bool flag = true; //true means it has opacity
			
					//bineary search -like algorithm
					//find the position, that is early termination point
					do
					{
						ray->getHistByIdm( &hist, 0, mid );
						float opacity = calculateSegmentOpacity( &hist, colorModulation );
						if( opacity < earlyTerminationThreshold ){//this mid point cannot be terminated( opacity < earlyTerminationThreshold )
							lastNotTerminatePoint = mid;
							min = mid;
							flag = false;
						}else{
							max = mid;
							flag = true;
						}
						mid = ( int )( ( min + max ) / 2.0 );
						if( max - min < step / 2.0 )break;
					}while( max - lastNotTerminatePoint  > step || flag == false );
					earlyTerminationPosition = max;
				}

				//check earlyTerminationPosition and endRenderLeapoingPosition, which one is more close to the tail
				if( earlyTerminationPosition < endRenderLeapingPosition )
					endRenderPosition = earlyTerminationPosition;
				else
					endRenderPosition = endRenderLeapingPosition;

				//store the render from and to information
				if( endRenderPosition < startRenderPosition )
					endRenderPosition = startRenderPosition;
				renderFrom[ rayIdx ] = startRenderPosition;
				renderTo[ rayIdx ] = endRenderPosition;
				#pragma omp critical
				{
					sumWorkload += ( renderTo[ rayIdx ] - renderFrom[ rayIdx ] );
				}

				cr = ( renderTo[ rayIdx ] - renderFrom[ rayIdx ] ) / (float) m_d;
                cg = ( renderTo[ rayIdx ] - renderFrom[ rayIdx ] ) / (float) m_d;
                cb = ( renderTo[ rayIdx ] - renderFrom[ rayIdx ] ) / (float) m_d;
				/*cr = ( renderTo[ rayIdx ] ) / (float) m_d;
                cg = ( renderTo[ rayIdx ] ) / (float) m_d;
                cb = ( renderTo[ rayIdx ] ) / (float) m_d;*/
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				m_color[(v*m_u+u)*4+0] =  optR;
                m_color[(v*m_u+u)*4+1] =  optG;
                m_color[(v*m_u+u)*4+2] =  optB;
				m_color[(v*m_u+u)*4+3] =  255;//A
            }//u loop
        }//v loop

		//dispatch the job and store the dispatch result in dispatchID
		//sumWorkload = m_u * m_v;	//for load imbalance
		float maxThreads = omp_get_max_threads();
		float avgWorkload = sumWorkload / maxThreads;

		int currentID = 0;
		float workload = 0;
		for( int v = 0; v< m_v; v+=1 ){
		    for( int u = 0; u< m_u; u+=1 ){
				int rayIdx = v * m_u + u;
				dispatchID[ rayIdx ] = currentID;
				workload += ( renderTo[ rayIdx ] - renderFrom[ rayIdx ] );	//for load balance
				//workload ++;	//for load imbalane
				if( workload >= avgWorkload ){
					workload = 0;
					currentID ++;
				}
			}
		}

		#pragma omp parallel for
		for( int p = 0; p < (int) maxThreads ; p ++ ){
			int threadID = omp_get_thread_num();
			
			clock_t start_time, end_time;
			float total_time = 0;
			start_time = clock(); /* mircosecond */

			int step = 8;
			for( int v = 0; v< m_v; v+=1 ){
				Histogram hist(m_bins);
				Histogram befHist(m_bins);
				for( int u = 0; u< m_u; u+=1 ){
					int rayIdx = v * m_u + u;
					if( threadID != dispatchID[ rayIdx ] )continue;

					Ray* ray = m_rays[ rayIdx ];
					float cr = 0, cg = 0, cb = 0, ca = 0;

					for( int d = renderFrom[ rayIdx ]; d< renderTo[ rayIdx ] - step; d += step ){
					//for( int d = 0; d < m_d - step; d += step ){
						float frontSegTransparency, correctOpacity;

						if( d>0 )
							ray->getHistByIdm(&befHist, 0, d-1);
						else
							befHist.reset();
						ray->getHistByIdm(&hist, d, d+step-1);

						float rr, gg, bb;
						calculateSegmentColorContribution( &befHist, &hist, colorModulation, rr, gg, bb, frontSegTransparency, correctOpacity );

						cr += rr;
						cg += gg;
						cb += bb;
						ca = ( 1 - frontSegTransparency ) + correctOpacity * frontSegTransparency; //not necessary
					}

					//for dispatch visualization
					/*float r, g, b, a;
					getRGBAbyTF( p / 8.0 , r, g, b, a, 1);
					cr = r;
					cg = g;
					cb = b;*/
					unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
					unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
					unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
					m_color[(v*m_u+u)*4+0] =  optR;
					m_color[(v*m_u+u)*4+1] =  optG;
					m_color[(v*m_u+u)*4+2] =  optB;
					m_color[(v*m_u+u)*4+3] =  255;//A
				}//u loop
			}//v loop	
			end_time = clock();
			char fn[500];
			sprintf( fn, "time%d.txt", p );
			FILE* fp = fopen( fn, "w" );
			total_time = (float)(end_time - start_time)/CLOCKS_PER_SEC;
			fprintf(fp, "Time : %f sec \n", total_time);
			fclose( fp );

		}//p loop

		free( dispatchID );
		free( renderTo );
		free( renderFrom );

		end_time = clock();
		FILE* fp = fopen( "tmp.txt", "w" );
		total_time = (float)(end_time - start_time)/CLOCKS_PER_SEC;
		fprintf(fp, "Time : %f sec \n", total_time);
		fclose( fp );
    }

	void calculateSegmentColorContribution( Histogram* histFront, Histogram* hist, float colorModulation, 
											float &rr, float &rg, float &rb, float &frontSegTransparency, float& correctOpacity )
	{
		
		//calculate the transparency of segment before current target segment
		frontSegTransparency = 1.0;
		frontSegTransparency = 1.0 - calculateSegmentOpacity( histFront, colorModulation );

		//calculate correct opacity and sum of opacity of target segment
		correctOpacity = 1;
		float sumOpacity = 0;
		correctOpacity = calculateSegmentOpacity( hist, colorModulation );

		//calculate color of target segment, but this is not correct
		float tmpR = 0, tmpG = 0, tmpB = 0;
		for( int i=0; i< m_bins; i++ ){
			float r, g, b, a;
			getRGBAbyTF( (i/(float)m_bins) , r, g, b, a, 1);
			a *= colorModulation;
			float freq = hist->getBin(i);
			tmpR += r * a * freq;
			tmpG += g * a * freq;
			tmpB += b * a * freq;
			sumOpacity += a * freq;
		}

		//calculate correct color and add into pixel color
		float opacityReachCamera = 0.0;
		if( sumOpacity > 0.0 ) 
			opacityReachCamera = frontSegTransparency * ( correctOpacity / sumOpacity );

		rr = tmpR * opacityReachCamera; 
		rg = tmpG * opacityReachCamera; 
		rb = tmpB * opacityReachCamera; 

		//the following is for mean color rendering, remove this part will back to histogram rendering directly
		Histogram h( hist );
		h.normalize();
		float mean = h.getMean();
		float r, g, b, a;
		getRGBAbyTF( mean , r, g, b, a, 1);
		rr = r * frontSegTransparency * correctOpacity;
		rg = g * frontSegTransparency * correctOpacity;
		rb = b * frontSegTransparency * correctOpacity;
	}

	bool entropyBasedRefineRendering()
	{  
		//initialize the UI inputing parameters
		float bgColor = m_sldBar8/(100.0);
		float colorModulation = m_lineEdit0;
		int step = m_lineEdit3;
		bool continueFlag = false;

		if( entropyQueue[0]->size() == 0 ){
			#pragma omp parallel for
			for( int v = 0; v< m_v; v+=1 ){
				Histogram hist(m_bins);
				Histogram befHist(m_bins);
				for( int u = 0; u< m_u; u+=1 ){
					int rayIdx = v*m_u+u;
					Ray* ray = m_rays[rayIdx];
					float cr = 0, cg = 0, cb = 0, ca = 0;	//color of a pixel
					float sumOfEntropy = 0;

					for( int d = 0; d< m_d-step; d+=step ){
						float frontSegTransparency, correctOpacity;

						if( d>0 )
							ray->getHistByIdm(&befHist, 0, d-1);
						else
							befHist.reset();
						ray->getHistByIdm(&hist, d, d+step-1);

						float rr, gg, bb;
						calculateSegmentColorContribution( &befHist, &hist, colorModulation, rr, gg, bb, frontSegTransparency, correctOpacity );

						cr += rr;
						cg += gg;
						cb += bb;
						ca = ( 1 - frontSegTransparency ) + correctOpacity * frontSegTransparency; //not necessary

						//calculate entropy and put in entropyQueue
						hist.normalize();
						float entropy = hist.getEntropy();
						Segment seg;
						seg.from = d;
						seg.len = step;
						seg.metric = frontSegTransparency * correctOpacity * entropy;
						seg.a = frontSegTransparency; seg.b = correctOpacity;
						std::vector<Segment>::iterator position;
						position = std::lower_bound (entropyQueue[rayIdx]->begin(), entropyQueue[rayIdx]->end(), seg, segmentSort );
						entropyQueue[rayIdx]->insert( position, seg );
						sumOfEntropy += seg.metric;
					}

					cr += (bgColor) * (1.0 - ca);
					cg += (bgColor) * (1.0 - ca);
					cb += (bgColor) * (1.0 - ca);
					unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
					unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
					unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
					m_color[(v*m_u+u)*4+0] =  optR;
					m_color[(v*m_u+u)*4+1] =  optG;
					m_color[(v*m_u+u)*4+2] =  optB;
					m_color[(v*m_u+u)*4+3] =  255;//A
					m_colorF[(v*m_u+u)*4+0] =  cr;
					m_colorF[(v*m_u+u)*4+1] =  cg;
					m_colorF[(v*m_u+u)*4+2] =  cb;

					//output entropy image
					float r, g, b, a;
					getRGBAbyTF( sumOfEntropy/7.0 , r, g, b, a, 1);	//7.0 is the maximum entropy if the bin number is 128
					m_color2[(v*m_u+u)*4+0] =  (unsigned char)(int)floor(r * 255.0);
					m_color2[(v*m_u+u)*4+1] =  (unsigned char)(int)floor(g * 255.0);
					m_color2[(v*m_u+u)*4+2] =  (unsigned char)(int)floor(b * 255.0);
					m_color2F[(v*m_u+u)*4+0] =  sumOfEntropy;
				}//u
			}//v
			continueFlag = true;
		}else{
			//not first run
			#pragma omp parallel for
			for( int v = 0; v< m_v; v+=1 ){
				Histogram hist(m_bins);
				Histogram frontHist(m_bins);
				for( int u = 0; u< m_u; u+=1 ){
					int rayIdx = v*m_u+u;
					Ray* ray = m_rays[rayIdx];
					vector<Segment>::iterator it = ( entropyQueue[rayIdx]->end() - 1 );
					int from = it->from;
					int len = it->len;
					int nextFrom = from + (int)( len / 2.0 ) ;
					entropyQueue[rayIdx]->pop_back();

					if( len >= 2 ){ 
						float sumOfEntropy = m_color2F[(v*m_u+u)*4+0];

						//old segment information
						float oldR, oldG, oldB;
						float oldFrontSegTransparency, oldCorrectOpacity;
						if( from > 0 )
							ray->getHistByIdm( &frontHist, 0, from - 1 );
						else
							frontHist.reset();
						ray->getHistByIdm( &hist, from, from + len - 1 );
						calculateSegmentColorContribution( &frontHist, &hist, colorModulation, oldR, oldG, oldB, oldFrontSegTransparency, oldCorrectOpacity );
						hist.normalize();
						float entropy = hist.getEntropy();
						sumOfEntropy -= ( oldFrontSegTransparency * oldCorrectOpacity * entropy );

						//first segment
						float r1, g1, b1;
						float frontSegTransparency1, correctOpacity1;
						if( from > 0 )
							ray->getHistByIdm( &frontHist, 0, from - 1 );
						else
							frontHist.reset();
						ray->getHistByIdm( &hist, from, nextFrom - 1 );
						calculateSegmentColorContribution( &frontHist, &hist, colorModulation, r1, g1, b1, frontSegTransparency1, correctOpacity1 );

						hist.normalize();
						float entropy1 = hist.getEntropy();
						Segment seg1;
						seg1.from = from;
						seg1.len = nextFrom - from;
						seg1.metric = frontSegTransparency1 * correctOpacity1 * entropy1;
						seg1.a = frontSegTransparency1; seg1.b = correctOpacity1;
						std::vector<Segment>::iterator position1;
						position1 = std::lower_bound (entropyQueue[rayIdx]->begin(), entropyQueue[rayIdx]->end(), seg1, segmentSort );
						entropyQueue[rayIdx]->insert( position1, seg1 );
						sumOfEntropy += seg1.metric;

						//second segment
						float r2, g2, b2;
						float frontSegTransparency2, correctOpacity2;
						if( nextFrom > 0 )
							ray->getHistByIdm( &frontHist, 0, nextFrom - 1 );
						else
							frontHist.reset();
						ray->getHistByIdm( &hist, nextFrom, from + len - 1 );
						calculateSegmentColorContribution( &frontHist, &hist, colorModulation, r2, g2, b2, frontSegTransparency2, correctOpacity2 );
						hist.normalize();
						float entropy2 = hist.getEntropy();
						Segment seg2;
						seg2.from = nextFrom;
						seg2.len = from + len - nextFrom;
						seg2.metric = frontSegTransparency2 * correctOpacity2 * entropy2;
						seg2.a = frontSegTransparency2; seg2.b = correctOpacity2;
						std::vector<Segment>::iterator position2;
						position2 = std::lower_bound (entropyQueue[rayIdx]->begin(), entropyQueue[rayIdx]->end(), seg2, segmentSort );
						entropyQueue[rayIdx]->insert( position2, seg2 );
						sumOfEntropy += seg2.metric;

						//update color
						float changeR = ( - oldR + r1 + r2 );
						float changeG = ( - oldG + g1 + g2 );
						float changeB = ( - oldB + b1 + b2 );
						m_colorF[(v*m_u+u)*4+0] =  m_colorF[(v*m_u+u)*4+0] + changeR;
						m_colorF[(v*m_u+u)*4+1] =  m_colorF[(v*m_u+u)*4+1] + changeG;
						m_colorF[(v*m_u+u)*4+2] =  m_colorF[(v*m_u+u)*4+2] + changeB;
						unsigned char optR = (unsigned char)(int)floor(m_colorF[(v*m_u+u)*4+0] * 255.0);
						unsigned char optG = (unsigned char)(int)floor(m_colorF[(v*m_u+u)*4+1] * 255.0);
						unsigned char optB = (unsigned char)(int)floor(m_colorF[(v*m_u+u)*4+2] * 255.0);
						m_color[(v*m_u+u)*4+0] =  optR;
						m_color[(v*m_u+u)*4+1] =  optG;
						m_color[(v*m_u+u)*4+2] =  optB;
						m_color[(v*m_u+u)*4+3] =  255;//A

						//stop or not
						if( changeR > 0.075 || changeG > 0.075 || changeB > 0.075 ){
							#pragma omp critical
							{
								continueFlag = true;
							}
						}

						//output entropy image
						float r, g, b, a;
						getRGBAbyTF( sumOfEntropy / 7.0 , r, g, b, a, 1);	//7.0 is the maximum entropy if the bin number is 128
						m_color2[(v*m_u+u)*4+0] =  (unsigned char)(int)floor(r * 255.0);
						m_color2[(v*m_u+u)*4+1] =  (unsigned char)(int)floor(g * 255.0);
						m_color2[(v*m_u+u)*4+2] =  (unsigned char)(int)floor(b * 255.0);
					}
				}
			}
		}

		return continueFlag;
   }


	

	void emphasizeIsovalueInBoxRegion( Histogram *reHist, int mouseXDown, int mouseYDown, int mouseXNow, int mouseYNow, int binSelect, 
															int rangeS, int rangeT, int rangeRatio, float enhanceS, float enhanceT, float bgColor )
    {
		printf("start to run by boundary search\n");
		//not real transfrer function analysis
		int mtrlBoundaryCnt = 9;
		float mtrlBoundary[9] = {0.125*0.0, 0.125*1.0, 0.125*2.0,
								0.125*3.0, 0.125*4.0, 0.125*5.0,
								0.125*6.0, 0.125*7.0, 0.125*8.0};
	    int rayCnt=0;
        int cnt=0;
        Histogram hist(m_bins);

		//clock_t start_time, end_time;
		//float total_time = 0;
		//start_time = clock();

		float* sumHist = (float*)malloc( sizeof(float)*m_bins );
		memset( sumHist, 0, sizeof(float)*m_bins );
        
		#pragma omp parallel
		{
			float* privateHist = (float*) malloc( sizeof(float)*m_bins );
			memset( privateHist, 0, sizeof(float)*m_bins );

			#pragma omp for
			for( int v = 0; v< m_v; v+=1 ){
				Histogram hist(m_bins);
				for( int u = 0; u< m_u; u+=1 ){
					Ray* r = m_rays[v*m_u+u];
				
					int loc = -1;
					float cr = 0, cg = 0, cb = 0, ca = 0;

					if( mouseXDown>=u || u>=mouseXNow || mouseYDown>=v || v>=mouseYNow ){
						//out of selection region, do nothging
					}else{//in the select region, render and collect histgraom
						int curIdx, thick;
						//****************************************************
						loc = -1;
						int lastLoc= -1;
						int step=8;
						for( int d = 0; d< m_d-step; d+=step ){
							r->getHistByIdm(&hist, d, d+step);
							if( hist.getBin(binSelect) != 0 ){
								loc = d;
								break;
							}					
						}
						for( int d = m_d-20; d>=0 ; d-=step ){
							r->getHistByIdm(&hist, d, d+step);
							if( hist.getBin(binSelect) != 0 ){
								lastLoc = d;
								break;
							}					
						}

						//****************************************************

						if( loc != -1 ){
							//add histogram to private hist
							if( loc < lastLoc ){
								r->getHistByIdm(&hist, loc, lastLoc);
								for( int i=0; i<m_bins; i++ ){
									privateHist[i] += hist.getBin(i);
								}
							
								for( int d = loc; d< lastLoc; d+=8 ){
									r->getHistByIdm(&hist, d, d+8-1);
									float vxl = hist.getMean();

									float r, g, b, a;
									getRGBAbyTF(vxl, r, g, b, a, 1);
									a*=0.1;

									cr += (r*a) * (1.0 - ca);
									cg += (g*a) * (1.0 - ca);
									cb += (b*a) * (1.0 - ca);
									ca += a * (1.0 - ca);

									if( ca > 0.95 )break;
								}
							}
						}else{//in box , but no isovalue we want
							cr = 1;
							cg = 0;
							cb = 0;
						}

						

						cr += (bgColor) * (1.0 - ca);
						cg += (bgColor) * (1.0 - ca);
						cb += (bgColor) * (1.0 - ca);
						unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
						unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
						unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
						m_color[(v*m_u+u)*4+0] = optR;
						m_color[(v*m_u+u)*4+1] = optG;
						m_color[(v*m_u+u)*4+2] = optB;
						m_color[(v*m_u+u)*4+3] = 255;
					}
				}//u loop
				printf("rendering V: %d\n" ,v);
			}//v loop

			#pragma omp critical
			{
				for( int i=0; i<m_bins; i++ ){
					sumHist[i] += privateHist[i];
				}
			}

			free(privateHist );
		}//omp


		//Histogram outHist(m_bins, sumHist );
		//outHist.normalize();

		for( int i=0; i<m_bins; i++ )
			reHist->setBin( i, sumHist[i] );
		reHist->normalize();

		//output file
		/*FILE *fp; 
		 do{
			fp =  fopen( "C:\\GravityLabDataSet\\DistributionQuery2.txt", "w" );
		}while( fp == NULL );
		 for( int i=0; i<m_bins; i++ ){
			 fprintf( fp, "%f\n", outHist.getBin(i) );
		 }
		 fclose(fp);*/
		free(sumHist);

		//end_time = clock();
		//total_time += (float)(end_time - start_time)/CLOCKS_PER_SEC;
		//printf("Time: %f sec\n", total_time);
    }

	//uneven multiresolution, high entropy with high resolution
	void multiResolutionEntropyEstimationRendering( float bgColor )
    {
		//not real transfrer function analysis
		int step =8;
		float clrAdp = 0.15;
	    int rayCnt=0;
        int cnt=0;
        Histogram hist(m_bins);

		//clock_t start_time, end_time;
		//float total_time = 0;
		//start_time = clock();

		//cluster transfer function
		int* group = (int*) malloc( sizeof(int)*m_bins);
		int gid = 0;
		float mr=0, mg=0, mb=0, ma=0;
		float mc=0;
		int si=0;
		float tre = 0.1;
		for( int i=0; i<m_bins; i++ ){
			float vxl = i/(float)m_bins;
			float r, g, b, a;
			getRGBAbyTF( vxl, r, g, b, a, 1 );

			mc++;
			mr += r;
			mg += g;
			mb += b;
			ma += a;
			bool flagR = true, flagG = true, flagB = true, flagA = true;
			for( int j=si; j<=i; j++ ){
				float rr, gg, bb, aa;
				float v = j/(float)m_bins;
				getRGBAbyTF( v, rr, gg, bb, aa, 1 );
				if( fabs( mr/mc - rr ) >tre )flagR = false;
				if( fabs( mg/mc - gg ) >tre )flagG = false;
				if( fabs( mb/mc - bb ) >tre )flagB = false;
				if( fabs( ma/mc - aa ) >tre )flagA = false;
			}
			
			if( flagR == false || flagG == false || flagB == false || flagA == false  ){
				gid++;
				group[i] = gid;
				mc = 0;
				mr= 0; mg= 0; mb= 0; ma= 0;
				si = i;
			}else{
				group[i] = gid;
			}
			qDebug() << "GROUP" + QString::number(i) + "  " +QString::number(group[i]);
		}

		#pragma omp parallel for
		for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram clr(gid+1);
            for( int u = 0; u< m_u; u+=1 ){

				float segs=0;
				float ttEtp = 0;
				
				float cr = 0, cg = 0, cb = 0, ca = 0;
				int idx = v*m_u+u;
				Ray* r = m_rays[idx];

				////output all voxel
				/*if( u!=347 || v!=66)continue;
				FILE* fp = fopen( "voxelOutput.txt", "w" );
				for( int d=0; d<504-1; d+=1 ){
					r->getHistByIdm( &hist, d, d+1 );
					fprintf(fp, "%f\n", hist.getMean());
				}
				fclose(fp);*/

				if( m_mrrL[idx].size() == 0 ){
					r->getHistByIdm( &hist, 0, m_d-1 );
					clr.reset();
					for( int i=0; i<m_bins; i++ ){
						clr.setBin( group[i], clr.getBin(group[i]) + hist.getBin(i) );
					}
					hist.normalize();
					clr.normalize();
					float vxl = hist.getMean();
					//float etp = hist.getEntropy();
					float etp = clr.getEntropy();

					m_mrrV[idx].push_back( vxl );
					m_mrrE[idx].push_back( etp );
					m_mrrL[idx].push_back( m_d );
					float r, g, b, a;
					float sr, sg, sb, sa;
					getRGBAbyTF( vxl, r, g, b, a, 1 );
					a*=clrAdp;
					renderSameVoxel( sr, sg, sb, sa, r, g, b, a, (m_d-1)/step );
					cr = sr;
					cg = sg;
					cb = sb;

					//ray statistic information
					segs++;
					ttEtp += etp;
				}else{
					int s = 0;
					//qDebug()<<"Not first" + QString::number( m_mrrV[idx].size() ) + " " +QString::number(m_mrrE[idx][0]);
					for( int i=0; i<m_mrrV[idx].size(); i++ ){
						if( m_mrrE[idx][i]  > 1 ){//entropy threshold and need to separate
							int t = m_mrrL[idx][i];

							//erase data
							m_mrrV[idx].erase( m_mrrV[idx].begin()+i, m_mrrV[idx].begin()+i+1 ); 
							m_mrrE[idx].erase( m_mrrE[idx].begin()+i, m_mrrE[idx].begin()+i+1 );
							m_mrrL[idx].erase( m_mrrL[idx].begin()+i, m_mrrL[idx].begin()+i+1 );

							//calculate and set first segment
							r->getHistByIdm( &hist, s, s + (t-s)/2 - 1 );
							clr.reset();
							for( int i=0; i<m_bins; i++ ){
								clr.setBin( group[i], clr.getBin(group[i]) + hist.getBin(i) );
							}
							hist.normalize();
							clr.normalize();
							float vxl1 = hist.getMean();
							//float etp1 = hist.getEntropy();
							float etp1 = clr.getEntropy();

							float loc1 = s + (t-s)/2;
							m_mrrV[idx].insert( m_mrrV[idx].begin()+i, vxl1 );
							m_mrrE[idx].insert( m_mrrE[idx].begin()+i, etp1 );
							m_mrrL[idx].insert( m_mrrL[idx].begin()+i, loc1 );

							//calculate and set second segment
							r->getHistByIdm( &hist, s + (t-s)/2, t -1 );
							clr.reset();
							for( int i=0; i<m_bins; i++ ){
								clr.setBin( group[i], clr.getBin(group[i]) + hist.getBin(i) );
							}
							hist.normalize();
							clr.normalize();
							float vxl2 = hist.getMean();
							//float etp2 = hist.getEntropy();
							float etp2 = clr.getEntropy();
							float loc2 = t;
							m_mrrV[idx].insert( m_mrrV[idx].begin()+i+1, vxl2 );
							m_mrrE[idx].insert( m_mrrE[idx].begin()+i+1, etp2 );
							m_mrrL[idx].insert( m_mrrL[idx].begin()+i+1, loc2 );

							//acculuate color
							//calculate the color for the current segment
							float r, g, b, a;
							float sr, sg, sb, sa;
							getRGBAbyTF( vxl1, r, g, b, a, 1 );
							a*=clrAdp;
							renderSameVoxel( sr, sg, sb, sa, r, g, b, a, (loc1-s)/step );
							cr += sr*sa*(1-ca);
							cg += sg*sa*(1-ca);
							cb += sb*sa*(1-ca);
							ca += sa * (1-ca);
							segs++;
							ttEtp+=etp1;

							getRGBAbyTF( vxl2, r, g, b, a, 1 );
							a*=clrAdp;
							renderSameVoxel( sr, sg, sb, sa, r, g, b, a, (loc2-loc1)/step );
							cr += sr*sa*(1-ca);
							cg += sg*sa*(1-ca);
							cb += sb*sa*(1-ca);
							ca += sa * (1-ca);
							segs++;
							ttEtp+=etp1;
							
							//skip second segment
							i++;

							//set s
							s = loc2;
							if( ca >= 0.95 )break;
						}else{
							float r, g, b, a;
							float sr, sg, sb, sa;
							float vxl = m_mrrV[idx][i];
							getRGBAbyTF( vxl, r, g, b, a, 1 );
							a*=clrAdp;
							renderSameVoxel( sr, sg, sb, sa, r, g, b, a, (m_mrrL[idx][i]-s)/step);
							cr += sr*sa*(1-ca);
							cg += sg*sa*(1-ca);
							cb += sb*sa*(1-ca);
							ca += sa * (1-ca);
							segs++;
							ttEtp += m_mrrE[idx][i];
							s = m_mrrL[idx][i];
							if( ca >= 0.95 )break;
						}//else if( m_mrrE[idx][i]  > 0.01 )					
					}//for
				}//else for if( m_mrrL[idx].size() == 0 )

				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
				unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
				unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
				unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(idx)*4+0] = optR;
				m_color[(idx)*4+1] = optG;
				m_color[(idx)*4+2] = optB;
				m_color[(idx)*4+3] = 255;

				//entropy
				/*m_color2[(idx)*4+0] = (unsigned char)(int)(ttEtp/segs)*50;
				m_color2[(idx)*4+1] = (unsigned char)(int)(ttEtp/segs)*50;
				m_color2[(idx)*4+2] = (unsigned char)(int)(ttEtp/segs)*50;
				m_color2[(idx)*4+3] = 255;*/
				//segments
				float rr, gg, bb, aa;
				getRGBAbyTF( segs/40.0, rr, gg, bb, aa, 1 );
				m_color2[(idx)*4+0] = (unsigned char)(int)((rr)*255.0);
				m_color2[(idx)*4+1] = (unsigned char)(int)((gg)*255.0);
				m_color2[(idx)*4+2] = (unsigned char)(int)((bb)*255.0);
				m_color2[(idx)*4+3] = 255;

				//qDebug()<<QString::number(ttEtp)+" "+QString::number(segs) + " "+QString::number(ttEtp/segs);
            }//u loop
            printf("rendering U: %d\n" ,v);
        }//v loop

		float max = 0;
		float sum = 0;
		for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
				int idx = v*m_u+u;
				sum+=m_color2[(idx)*4+0]/255.0*100.0;
				if( max < m_color2[(idx)*4+0]/255.0*100.0 )max= m_color2[(idx)*4+0]/255.0*100.0;
			}
		}
		qDebug()<<QString::number( max );

		free(group);
		//end_time = clock();
		//total_time += (float)(end_time - start_time)/CLOCKS_PER_SEC;
		//printf("Time: %f sec\n", total_time);
    }

	//very basic multiresolution and take mean
	void multiResolutionColorRendering(  int step )
    {
		int rayCnt=0;
        int cnt=0;
        Histogram hist(m_bins);
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[rayCnt++];
                float cr = 0, cg = 0, cb = 0, ca = 0;
               
				float lastVxl = 0;
				float lastPos = 0;
				for( int d = 0; d< m_d-step; d+=step ){
                    ray->getHistByIdm(&hist, d, d+step-1);
                    float vxl = hist.getMean();

					float endPos = (d-step/2.0);
					for( int i=lastPos; i<endPos; i+=8 ){
						float r, g, b,a;
						float v = lastVxl + ( (vxl-lastVxl) * ((i-lastPos)/(float)(endPos-lastPos)));
						//float v= vxl;
						getRGBAbyTF(v, r, g, b, a, 1);
						a*=0.1;

						cr += (r*a) * (1.0 - ca);
						cg += (g*a) * (1.0 - ca);
						cb += (b*a) * (1.0 - ca);
						ca += a * (1.0 - ca);
					}

					lastVxl = vxl;
					lastPos = endPos;
                    if( ca > 0.95 )break;
                }

				cr+= (1-ca)*1;
				cg+= (1-ca)*1;
				cb+= (1-ca)*1;
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
                m_color[cnt++] =  optR;
                m_color[cnt++] =  optG;
                m_color[cnt++] =  optB;
				m_color[cnt++] =  255;
                //printf("%d %d %d\n", optR, optG, optB);
            }
            printf("rendering U: %d\n" ,v);
        }
    }
	//***********************DIRECT VOLUME REDNERING - END************************************

	//take sample from a distribution of a segment
	void sampleRendering( int step, int startDepth, float clrAdapt, float bgColor )
    {
		startDepth = (startDepth/99.0)*(m_d-step);
        int rayCnt=0;
        int cnt=0;
		int cnt2=0;
        
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

                int init, iterate;
				if( m_order == 0 ){
					init = 0+startDepth;
					iterate = step;
				}else{
					init = m_d-step-1-startDepth;
					iterate = -step;
				}
				for( int d = init; d>=0 && d< m_d-step; d+=iterate ){
                    hist.reset();
					ray->getHistByIdm( &hist, d, d+16 );
					
					//float vxl = hist.getMean();
					float vxl = hist.sampleDataPoint();
					//float vxl = hist.getVariance();


					float r, g, b, a;
                    getRGBAbyTF(vxl, r, g, b, a, 1);
					a*=clrAdapt;
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
                    if( ca > 0.95 )break;
                }
				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
            printf("rendering U: %d\n" ,v);
        }
    }

	void entropyBaseOpacityRendering( int step, int startDepth, float clrAdapt, float bgColor, float maxE, float curve = 1 )
    {
		float maxEntropy = 0;
        
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 0; d< m_d-step; d+=8 ){
					int from, to;
					from = (int)(d - step/2.0);
					to = (int)(d + step/2.0);
					if( from < 0 )from = 0;
					if( to > m_d-1 )to = m_d-1;


                    hist.reset();
					ray->getHistByIdm( &hist, d, d);
					hist.normalize();
					float vxl = hist.getMean();

					hist.reset();
					ray->getHistByIdm( &hist, from, to);
					hist.normalize();
					float entropy = hist.getEntropy();
					

					#pragma omp critical
					{
						if( maxEntropy < entropy )maxEntropy = entropy;
					}

					entropy /= maxE;
					entropy = pow( (double)entropy, (double)curve );
					float r, g, b, a;
					//getRGBAbyTF(entropy, r, g, b, a, 1);
                    getRGBAbyTF(vxl, r, g, b, a, 1);
					//a*=clrAdapt;
					a*=entropy*clrAdapt;
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
                    //if( ca > 0.95 )break;
                }
				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
        }

		FILE* fp = fopen( "tmp.txt", "w" );
		fprintf(fp, "%f\n", maxEntropy);
		fclose(fp );
    }

	void entropyBaseSampleRendering(  )
    {
		int step = m_lineEdit3;
		float clrAdapt = m_lineEdit0;
		float curve = m_lineEdit5;
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram befHist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 1; d< m_d-step; d+=step ){
					hist.reset();
					befHist.reset();
					ray->getHistByIdm( &befHist, 0, d - 1 );
					ray->getHistByIdm( &hist, d, d + step );
					//ray->getHistByIdm( &hist, d, d + 1 );
					
					hist.normalize();
					//float rr, gg, bb;
					//float frontSegTransparency, correctOpacity;
					//calculateSegmentColorContribution( &befHist, &hist, clrAdapt, rr, gg, bb, frontSegTransparency, correctOpacity );
					//float entropy = hist.getEntropy();
					float mean = hist.getMean();

					//for( int i=0; i<m_bins; i++ ){
					//	float r, g, b, a;
					//	getRGBAbyTF( i/(float)m_bins, r, g, b, a, 1);
					//	hist.setBin( i, hist.getBin( i ) * a );
					//}
					hist.normalize();
					
					//random sample
					//float opa = ( correctOpacity * frontSegTransparency );
					//if( opa > 0.0001 ){
					int numSamples = 20;
					if( numSamples > 1 ){
							float* samples = (float*) malloc( sizeof( float ) * numSamples );
							hist.samplePoints( samples, numSamples );
							float aa = 0;
						
							for( int i = 0; i < numSamples; i ++ ){
								float r, g, b, a;
								getRGBAbyTF( samples[ i ], r, g, b, a, 1);
								a *= clrAdapt;
								cr += (r*a) * ( 1.0 - ca );
								cg += (g*a) * ( 1.0 - ca );
								cb += (b*a) * ( 1.0 - ca );
								//ca += a * ( 1.0 - ca );
								float addCa = a * ( 1.0 - ca );
								ca += addCa;	
								aa += addCa;
					//			if( aa > opa )break;
							}

							free( samples );
					}
					//}

					//mean
					//float r, g, b, a;
					//getRGBAbyTF( mean, r, g, b, a, 1);
					//a *= clrAdapt;
					//cr += (r*a) * (1.0 - ca);
					//cg += (g*a) * (1.0 - ca);
					//cb += (b*a) * (1.0 - ca);
					//ca += a * (1.0 - ca);
     //               
                    if( ca > 0.95 )break;
                }
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
        }
    }

	void visibilityBaseSampleRendering(  )
    {
		int step = m_lineEdit3;
		float clrAdapt = m_lineEdit0;
		float curve = m_lineEdit5;
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 0; d< m_d-step; d+=step ){
					hist.reset();
					ray->getHistByIdm( &hist, d, d + step );
					hist.normalize();
					float entropy = hist.getEntropy();
					

					//entropy = pow( (double)entropy, (double)curve );
					//int numSamples = (int) ( entropy * 3.0 );
					int numSamples = 1000.0 * ( step / (float) m_d ) ;
					if( numSamples > 1 ){
						float* samples = (float*) malloc( sizeof( float ) * numSamples );
						hist.samplePoints( samples, numSamples );

						for( int i = 0; i < numSamples; i ++ ){
							float r, g, b, a;
							getRGBAbyTF( samples[ i ], r, g, b, a, 1);
							a *= clrAdapt;
							cr += (r*a) * (1.0 - ca);
							cg += (g*a) * (1.0 - ca);
							cb += (b*a) * (1.0 - ca);
							ca += a * (1.0 - ca);
						}

						free( samples );
					}
                    
                    if( ca > 0.95 )break;
                }
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
        }
    }

	void uncertaintyRendering(  )
    {
		int step = m_lineEdit3;
		float clrAdapt = m_lineEdit0;
		float curve = m_lineEdit5;
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 0; d< m_d-step; d+=step ){
					hist.reset();
					ray->getHistByIdm( &hist, d, d + step );

					//get weight histogram
					for( int i = 0; i < m_bins; i ++ ){
						float r, g, b, a;
						getRGBAbyTF( i / (float)m_bins, r, g, b, a, 1);
						a *= clrAdapt;
						hist.setBin( i, hist.getBin( i ) * a );
					}
					hist.normalize();

					//calculate weighted mean
					float rmean = 0, gmean = 0, bmean = 0;
					for( int i = 0; i < m_bins; i ++ ){
						float r, g, b, a;
						getRGBAbyTF( i / (float)m_bins, r, g, b, a, 1);
						float w = hist.getBin( i );
						rmean += r * w;
						gmean += g * w;
						bmean += b * w;
					}

					//calculate weighted variance
					float rvariance = 0, gvariance = 0, bvariance = 0;
					for( int i = 0; i < m_bins; i ++ ){
						float r, g, b, a;
						getRGBAbyTF( i / (float)m_bins, r, g, b, a, 1);
						float w = hist.getBin( i );
						rvariance += ( r - rmean ) * ( r - rmean ) * w;
						gvariance += ( g - gmean ) * ( g - gmean ) * w;
						bvariance += ( b - bmean ) * ( b - bmean ) * w;
					}
					rvariance = sqrt( rvariance / (float)m_bins );
					gvariance = sqrt( gvariance / (float)m_bins );
					bvariance = sqrt( bvariance / (float)m_bins );

					//pick up max varaince
					float variance = rvariance;
					if( gvariance > variance ) variance = gvariance;
					if( bvariance > variance ) variance = bvariance;
					
					if( ca > 0.95 )break;
                }
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
        }
    }

	void selfInformationBaseSampleRendering(  )
    {
		int step = m_lineEdit3;
		float clrAdapt = m_lineEdit0;
		float curve = m_lineEdit5;
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram selfInfoHist( m_bins );
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 0; d< m_d-step; d+=step ){
					hist.reset();
					ray->getHistByIdm( &hist, d, d + step );
					hist.normalize();
					
					//create self informatoin histogram
					selfInfoHist.reset();
					float log2 = log( 2.0 );
					for( int i = 0; i < m_bins; i ++ ){
						if( hist.getBin( i ) > 0.000001 )
							selfInfoHist.setBin( i, -1 * ( log( hist.getBin( i ) ) / log2 ) );
					}
					selfInfoHist.normalize();

					//entropy = pow( (double)entropy, (double)curve );
					//int numSamples = (int) ( entropy * 3.0 );
					int numSamples = 1000.0 * ( step / (float) m_d ) ;
					if( numSamples > 1 ){
						float* samples = (float*) malloc( sizeof( float ) * numSamples );
						selfInfoHist.samplePoints( samples, numSamples );

						for( int i = 0; i < numSamples; i ++ ){
							float r, g, b, a;
							getRGBAbyTF( samples[ i ], r, g, b, a, 1);
							a *= clrAdapt;
							cr += (r*a) * (1.0 - ca);
							cg += (g*a) * (1.0 - ca);
							cb += (b*a) * (1.0 - ca);
							ca += a * (1.0 - ca);
						}

						free( samples );
					}
                    
                    if( ca > 0.95 )break;
                }
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
        }
    }

	//very basic volume rendering
	void tranditionalVolumeRendering( int step, int startDepth, float clrAdapt, float bgColor, float maxE )
    {
		startDepth = (startDepth/99.0)*(m_d-step);
        int rayCnt=0;
        int cnt=0;
		int cnt2=0;

		//float maxE = 2;
		float maxEntropy = 0;
        
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 0; d< m_d-step; d+=8 ){
                    hist.reset();
					ray->getHistByIdm( &hist, d, d+step);
					hist.normalize();
					float vxl = hist.getMean();
					
					//if( hist.getBin(80)!=0 )vxl = 0.625;

					float r, g, b, a;
					getRGBAbyTF(vxl, r, g, b, a, 1);
					a*=clrAdapt;
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
                    //if( ca > 0.95 )break;
                }
				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
        }

		FILE* fp = fopen( "tmp.txt", "w" );
		fprintf(fp, "%f\n", maxEntropy);
		fclose(fp );
    }

	//summerize flicker affect code
	void summerizeSampleRendering( int step, int startDepth, float clrAdapt, float bgColor )
    {
		startDepth = (startDepth/99.0)*(m_d-step);
        int rayCnt=0;
        int cnt=0;
		int cnt2=0;

		float *temp = (float*) malloc( sizeof(float)*(m_v*m_u ) );
        
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float rAry[200];
				float gAry[200];
				float bAry[200];
				int times = 100;
				for( int i=0; i<times; i++ ){
					float cr = 0, cg = 0, cb = 0, ca = 0;
					for( int d = 0; d< m_d-step; d+=step ){
						hist.reset();
						ray->getHistByIdm( &hist, d, d+step );
					
						//float vxl = hist.getMean();
						float vxl = hist.sampleDataPoint();
						//float vxl = hist.getVariance();


						float r, g, b, a;
						getRGBAbyTF(vxl, r, g, b, a, 1);
						a*=clrAdapt;
                    
						cr += (r*a) * (1.0 - ca);
						cg += (g*a) * (1.0 - ca);
						cb += (b*a) * (1.0 - ca);
						ca += a * (1.0 - ca);
						if( ca > 0.95 )break;
					}
					cr += (bgColor) * (1.0 - ca);
					cg += (bgColor) * (1.0 - ca);
					cb += (bgColor) * (1.0 - ca);

					rAry[i] = cr;
					gAry[i] = cg;
					bAry[i] = cb;
				}
				//calculate statistics
				float rMean=0, gMean=0, bMean=0;
				float rVariance=0, gVariance=0, bVariance=0;
				for( int i=0; i<times; i++ ){
					rMean += rAry[i];
					gMean += gAry[i];
					bMean += bAry[i];
				}
				rMean/=(float)times;
				gMean/=(float)times;
				bMean/=(float)times;

				for( int i=0; i<times; i++ ){
					rVariance += pow( (double)(rAry[i]-rMean), 2.0 );
					gVariance += pow( (double)(rAry[i]-rMean), 2.0 );
					bVariance += pow( (double)(rAry[i]-rMean), 2.0 );
				}
				rVariance /= (float)times;
				gVariance /= (float)times;
				bVariance /= (float)times;

				float var = (rVariance + gVariance + bVariance )/3.0;

				temp[v*m_u+u] = var;

                unsigned char optR = (unsigned char)(int)floor(var * 255.0);
                unsigned char optG = (unsigned char)(int)floor(var * 255.0);
                unsigned char optB = (unsigned char)(int)floor(var * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
            printf("rendering U: %d\n" ,v);
        }

		float maxVar = 0;
		for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
				if( maxVar < temp[v*m_u+u] )maxVar =  temp[v*m_u+u];
			}
		}

		for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
				float vxl = temp[v*m_u+u]/maxVar;
				vxl = log(vxl);
				float base = -5;
				if( vxl <base ) vxl = base;
				vxl = (vxl-base)/(0-base);
				if( vxl > 0.99999 )vxl=1;
				if( vxl<0)vxl = 0;

				float r, g, b, a;
                getRGBAbyTF(vxl, r, g, b, a, 1);

				unsigned char optR = (unsigned char)(int)floor(r * 255.0);
                unsigned char optG = (unsigned char)(int)floor(g * 255.0);
                unsigned char optB = (unsigned char)(int)floor(b * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
			}
		}

		free(temp);
    }

	float uncertaintyRendering( int step, int startDepth, float clrAdapt, float bgColor, float normFactor=0 )
    {
		startDepth = (startDepth/99.0)*(m_d-step);
        int rayCnt=0;
        int cnt=0;
		int cnt2=0;

		float maxEntropy = 0;
		float returnMaxEntropy;
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				for( int d = 0; d< m_d-step; d+=step ){
                    hist.reset();
					ray->getHistByIdm( &hist, d, d+step );
					hist.normalize();
					float vxl = hist.getEntropy();

					#pragma omp critical
					{
						if( maxEntropy < vxl )
							maxEntropy = vxl;
					}
				}
            }
        }
		returnMaxEntropy = maxEntropy;
		if( normFactor!=0 )maxEntropy=normFactor;

        
		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];
				float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 0; d< m_d-step; d+=step ){
                    hist.reset();
					ray->getHistByIdm( &hist, d, d+step );
					
					//float vxl = hist.getMean();
					//float vxl = hist.sampleDataPoint();
					//float vxl = hist.getVariance();
					hist.normalize();
					float vxl = hist.getEntropy();

					vxl /= maxEntropy;
					vxl = log(vxl);
					float base = -2.5;
					if( vxl <base ) vxl = base;
					vxl = (vxl-base)/(0-base);
					if( vxl > 0.99999 )vxl=1;
					if( vxl<0)vxl = 0;

					float r, g, b, a;
                    getRGBAbyTF(vxl, r, g, b, a, 1);
					a*=clrAdapt;
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
                    if( ca > 0.95 )break;
                }
				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
            printf("rendering U: %d\n" ,v);
        }

		return returnMaxEntropy;
    }

	
    
	

	
	//*********************statistic rendering - variance**********************

	//statisticRenderingEntropy: render variance by step
	//transfer function is from m_tf, so call makeTrasferfunction in advance is necesary
	//step: a voxel block(slices)

    void statisticRenderingEntropy( int step )
    {
        int rayCnt=0;
        int cnt=0;
		float maxEn = 0;
		float minEn = 10000;
		float maxEnV = 0;
		float minEnV = 10000;
        
		FILE* fp0 = fopen("tmpOut.txt", "w" );
		#pragma omp for
		for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram histSI(m_bins);
			Histogram histNow(m_bins);
			Histogram histVis(m_bins);

			for( int u = 0; u< m_u; u+=1 ){
				Ray* ray = m_rays[v*m_u+u];
				ray->getHistByIdm(&hist, 0, m_d-1);
				hist.normalize();
				float e = hist.getEntropy();
				/*#pragma omp critical
				{
					if( maxEn < e )maxEn = e;
					if( minEn > e )minEn = e;
				}*/

				//calcualte SI
				histSI.reset();
				for( int i=0; i<m_bins; i++ ){
					float v;
					if( hist.getBin(i) < 0.000001 ) v = 0;
					else if( hist.getBin(i) >0.99999 ) v= 0;
					else v = (-1*(log(hist.getBin(i))/log(2.0)));
					histSI.setBin( i, v );
				}

				//visibility
				histVis.reset();

				float cr = 0, cg = 0, cb = 0, ca = 0;
				for( int d=0; d<m_u-step; d+=step ){
					ray->getHistByIdm(&histNow, d, d+step-1);
					float vxl = histNow.getMean();
					float r, g, b, a;
					getRGBAbyTF(vxl, r, g, b, a, 1);
					a*=0.1;
					ca += a * (1.0 - ca);
					int bn = (int)(vxl*(float)m_bins);
					histVis.setBin(bn, histVis.getBin(bn)+ ca );
				}
				histVis.normalize();
				float ev = histVis.getEntropy();
			
				float m1=0, m2=0;
				for( int i=0; i<m_bins; i++ ){
					m1 += histSI.getBin( i );
					m2 += histVis.getBin( i );
				}
				m1/=(float)m_bins;
				m2/=(float)m_bins;
				for( int i=0; i<m_bins; i++ ){
					histSI.setBin( i, histSI.getBin( i )-m1 );
					histVis.setBin( i, histVis.getBin( i )-m2 );
				}
				float n1 = histSI.getNorm();
				float n2 = histVis.getNorm();
				float n=0;
				for( int i=0; i<m_bins; i++ ){
					n += histSI.getBin( i )* histVis.getBin( i );
				}

				/*#pragma omp critical
				{
					if( maxEnV < e )maxEnV = e;
					if( minEnV > e )minEnV = e;
				}*/
				#pragma omp critical
				{
					if( maxEn < n )maxEn = n;
					if( minEn > n )minEn = n;
				}
			}
		}


		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram histSI(m_bins);
			Histogram histNow(m_bins);
			Histogram histVis(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[v*m_u+u];

				ray->getHistByIdm(&hist, 0, m_d-1);
				hist.normalize();
				float e = hist.getEntropy();
				//e = (e-minEn)/(maxEn-minEn);
				//float r, g, b, a;
                //getRGBAbyTF(e, r, g, b, a, 1);

				//calcualte SI
				histSI.reset();
				for( int i=0; i<m_bins; i++ ){
					float v;
					if( hist.getBin(i) < 0.000001 ) v = 0;
					else if( hist.getBin(i) >0.99999 ) v= 0;
					else v = (-1*(log(hist.getBin(i))/log(2.0)));
					histSI.setBin( i, v );
				}

				//visibility
				histVis.reset();
				float cr = 0, cg = 0, cb = 0, ca = 0;
				for( int d=0; d<m_u-step; d+=step ){
					ray->getHistByIdm(&histNow, d, d+step-1);
					float vxl = histNow.getMean();
					float r, g, b, a;
					getRGBAbyTF(vxl, r, g, b, a, 1);
					a*=0.1;
					ca += a * (1.0 - ca);
					int bn = (int)(vxl*(float)m_bins);
					histVis.setBin(bn, histVis.getBin(bn)+ ca );
				}
				histVis.normalize();
				float ev = histVis.getEntropy();
				//ev = (e-minEnV)/(maxEnV-minEnV);
				//float rv, gv, bv, av;
                //getRGBAbyTF(e, rv, gv, bv, av, 1);

				float m1=0, m2=0;
				for( int i=0; i<m_bins; i++ ){
					m1 += histSI.getBin( i );
					m2 += histVis.getBin( i );
				}
				m1/=(float)m_bins;
				m2/=(float)m_bins;
				for( int i=0; i<m_bins; i++ ){
					histSI.setBin( i, histSI.getBin( i )-m1 );
					histVis.setBin( i, histVis.getBin( i )-m2 );
				}
				float n1 = histSI.getNorm();
				float n2 = histVis.getNorm();
				float n=0;
				for( int i=0; i<m_bins; i++ ){
					n += histSI.getBin( i )* histVis.getBin( i );
				}

				float fr=0.5, fg=0, fb=0, fa=0;
				float ef = n;
				ef = (ef-minEn)/(maxEn-minEn);
				getRGBAbyTF(ef, fr, fg, fb, fa, 1);

                unsigned char optR = (unsigned char)(int)floor(fr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(fg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(fb * 255.0);
                m_color[(v*m_u+u)*4+0] =  optR;
                m_color[(v*m_u+u)*4+1] =  optG;
                m_color[(v*m_u+u)*4+2] =  optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
        }
    }


	void adaptiveTrasferFunction( int step, float isovalue, float visStandard, float opaAdapt, float msDownX, float msDownY, float msRlsX, float msRlsY, float disapear )
    {
		if( disapear == 0 ){
			//reset m_depth
			memset( m_depth, 0, m_u * m_v * sizeof( float ) );

			//openmp paprallize: every column pixel for a thread
			#pragma omp parallel for

			//v is vertical axis of image
			for( int v = 0; v < m_v; v += 1 ){
				Histogram hist( m_bins );
				Histogram befHist( m_bins );

				//u is horizontal zxis of image
				for( int u = 0; u < m_u; u += 1 ){
					Ray* ray = m_rays[ v * m_u + u ];
					float cr = 0, cg = 0, cb = 0, ca = 0;	//rgba of this pixel
					float cbr = 0, cbg = 0, cbb = 0, cba = 0;	//rgba behind the position in m_depth

					//d is depth of a pixelr ay
					for( int d = 0; d < m_d - step; d += step ){
						float transp = 1.0;	//remained tranparency before 'd'

						//evaluate remianed transparency in front of 'd'
						//d == 0 dont evaluate
						if( d > 0 ){
							ray->getHistByIdm( &befHist, 0, d - 1 );
							for( int i = 0; i < m_bins; i ++ ){

								//evaluate by (1-a)^x * (1-b)^y * (1-c)^z........
								float rr, gg, bb, aa;
								getRGBAbyTF( ( i / (float) m_bins ), rr, gg, bb, aa, 1);
								aa *= opaAdapt;
								float freq = befHist.getBin( i );
								transp *= pow( (double) (1.0 - aa ), (double) freq );
							}
						}

						//get histogram of current block
						ray->getHistByIdm( &hist, d, d + step - 1 );

						//current block has desired block or not
						//and set depth to m_depth
						if( hist.getBin( (int) ( isovalue * m_bins ) ) > 0 && m_depth[ v * m_u + u ] == 0 ){
						//if( d > 200 && ( u > 195 && u < 370 ) && ( v > 40 && v < 80 ) && m_depth[ v * m_u + u ] == 0 ){
							m_depth[ v * m_u + u ] = d;

							//m_color3 is the color before the data we want to emphasize
							m_color3[ ( v * m_u + u ) * 4 + 0 ] =  (unsigned char) (int) floor( cbr * 255.0 );
							m_color3[ ( v * m_u + u ) * 4 + 1 ] =  (unsigned char) (int) floor( cbg * 255.0 );
							m_color3[ ( v * m_u + u ) * 4 + 2 ] =  (unsigned char) (int) floor( cbb * 255.0 );
							cbr = 0, cbg = 0, cbb = 0, cba = 0;
						}
                    
						float r = 0, g = 0, b = 0, a = 0;
						float vxlCnt = 0;
						float alpha = 1;	//the opacity this opacity can use
						float sumAlpha = 0;	//sum of opcity of this block
						for( int i = 0; i < m_bins; i ++ ){
							float rr, gg, bb, aa;
							getRGBAbyTF( ( i / (float) m_bins ) , rr, gg, bb, aa, 1);
							aa *= opaAdapt;
							float freq = hist.getBin( i );
							alpha *= pow( (double) ( 1.0 - aa ), (double) freq );
							sumAlpha += aa * freq;
						}
						alpha = 1.0 - alpha;
						float aRatio = 1;	//color adaptive of this block
						if( sumAlpha > 0 )aRatio = alpha / sumAlpha;

						//evalute the color of current block
						for( int i = 0; i < m_bins; i ++ ){
							float rr, gg, bb, aa;
							getRGBAbyTF( ( i / (float) m_bins ) , rr, gg, bb, aa, 1);
							aa *= opaAdapt;
							float freq = hist.getBin( i );
							cr += transp * aa * aRatio * rr * freq;
							cg += transp * aa * aRatio * gg * freq;
							cb += transp * aa * aRatio * bb * freq;
							cbr += transp * aa * aRatio * rr * freq;
							cbg += transp * aa * aRatio * gg * freq;
							cbb += transp * aa * aRatio * bb * freq;
						}
						ca = ( 1 - transp ) + alpha * transp;
						cba = ( 1 - transp ) + alpha * transp;
					}

					//set image color to image buffer
					float bgColor = 0;	//back groud color
					cr += ( bgColor ) * ( 1.0 - ca );
					cg += ( bgColor ) * ( 1.0 - ca );
					cb += ( bgColor ) * ( 1.0 - ca );
					unsigned char optR = (unsigned char) (int) floor( cr * 255.0 );
					unsigned char optG = (unsigned char) (int) floor( cg * 255.0 );
					unsigned char optB = (unsigned char) (int) floor( cb * 255.0 );
					m_color[ ( v * m_u + u ) * 4 + 0 ] =  optR;
					m_color[ ( v * m_u + u ) * 4 + 1 ] =  optG;
					m_color[ ( v * m_u + u ) * 4 + 2 ] =  optB;
					m_color[ ( v * m_u + u ) * 4 + 3 ] =  255;//A

					cbr += ( bgColor ) * ( 1.0 - cba );
					cbg += ( bgColor ) * ( 1.0 - cba );
					cbb += ( bgColor ) * ( 1.0 - cba );
					optR = (unsigned char) (int) floor( cbr * 255.0 );
					optG = (unsigned char) (int) floor( cbg * 255.0 );
					optB = (unsigned char) (int) floor( cbb * 255.0 );
					m_color2[ ( v * m_u + u ) * 4 + 0 ] =  optR;
					m_color2[ ( v * m_u + u ) * 4 + 1 ] =  optG;
					m_color2[ ( v * m_u + u ) * 4 + 2 ] =  optB;
					m_color2[ ( v * m_u + u ) * 4 + 3 ] =  255;//A
				}
			}
		}

		//adaptive rendering
		float clrAdapt = disapear;
		//openmp paprallize: every column pixel for a thread
		#pragma omp parallel for

		//v is vertical axis of image
        for( int v = 0; v < m_v; v += 1 ){
			Histogram hist( m_bins );

			//u is horizontal zxis of image
            for( int u = 0; u < m_u; u += 1 ){
				Ray* ray = m_rays[ v * m_u + u ];
				if( m_depth[ v * m_u + u ] == 0 ){
					m_color[ ( v * m_u + u ) * 4 + 0 ] =  m_color2[ ( v * m_u + u ) * 4 + 0 ];
					m_color[ ( v * m_u + u ) * 4 + 1 ] =  m_color2[ ( v * m_u + u ) * 4 + 1 ];
					m_color[ ( v * m_u + u ) * 4 + 2 ] =  m_color2[ ( v * m_u + u ) * 4 + 2 ];
				}else{
					ray->getHistByIdm( &hist, 0,  m_depth[ v * m_u + u ] - 1 );
					float transp = 1.0;
					for( int i = 0; i < m_bins; i ++ ){
						//evaluate by (1-a)^x * (1-b)^y * (1-c)^z........
						float rr, gg, bb, aa;
						getRGBAbyTF( ( i / (float) m_bins ), rr, gg, bb, aa, 1);
						aa *= ( clrAdapt * opaAdapt );
						float freq = hist.getBin( i );
						transp *= pow( (double) (1.0 - aa ), (double) freq );
					}

					float transpOri = 1.0;
					for( int i = 0; i < m_bins; i ++ ){
						//evaluate by (1-a)^x * (1-b)^y * (1-c)^z........
						float rr, gg, bb, aa;
						getRGBAbyTF( ( i / (float) m_bins ), rr, gg, bb, aa, 1);
						aa *= ( opaAdapt );
						float freq = hist.getBin( i );
						transpOri *= pow( (double) (1.0 - aa ), (double) freq );
					}

					float opa = ( 1 - transp ) / ( 1 - transpOri );

					//set image color to image buffer
					m_color[ ( v * m_u + u ) * 4 + 0 ] =  m_color3[ ( v * m_u + u ) * 4 + 0 ] * opa + ( m_color2[ ( v * m_u + u ) * 4 + 0 ] / transpOri ) * transp;
					m_color[ ( v * m_u + u ) * 4 + 1 ] =  m_color3[ ( v * m_u + u ) * 4 + 1 ] * opa + ( m_color2[ ( v * m_u + u ) * 4 + 1 ] / transpOri ) * transp;
					m_color[ ( v * m_u + u ) * 4 + 2 ] =  m_color3[ ( v * m_u + u ) * 4 + 2 ] * opa + ( m_color2[ ( v * m_u + u ) * 4 + 2 ] / transpOri ) * transp;
				}//else



            }//u loop
        }//v loop

    }

	void emphasize2()
	{
		int msDownX = -1;
		int msRlsX, msDownY, msRlsY;
		float isovalue = m_lineEdit1;
		float visStandard = 0.6;
		float opaAdapt = 0.05;
		int step = 8;


		Histogram hist(m_bins);
		Histogram flag(m_bins);
		flag.reset();

		if( msDownX == -1 ){
			msDownX = -1;
			msRlsX = m_u;
			msDownY= -1;
			msRlsY = m_v;
		}

		int *isoFlag = (int*) malloc( sizeof(int)* m_u* m_v );
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
                Ray* r = m_rays[v*m_u+u];
                
				if( msDownX<u && u<msRlsX && msDownY<v && v<msRlsY ){
					int thick;
					int curIdx;
					int loc = r->findIsovalue(isovalue, curIdx, thick, -1 );
					float color = loc/(float)(m_d);
				
					isoFlag[v*m_u+u] = loc;
				}else{
					isoFlag[v*m_u+u] = -1;
				}
            }
        }

		FILE* fp = fopen( "tmp.txt", "w" );
		float visibility = 0;
		float minVisibility = 1;
		float visCount = 0;
		float adtTransparent = 1;
		float width = 0.5;
		int bin = (int)(isovalue*(float)m_bins);
		do{
			visibility = 0;
			visCount= 0;
			minVisibility = 1;
			adtTransparent -= 0.01;
			width -= 1/(float)m_bins;
			for( int v = 0; v< m_v; v+=1 ){
				for( int u = 0; u< m_u; u+=1 ){
					if( isoFlag[v*m_u+u] == -1 )continue;

					Ray* r = m_rays[v*m_u+u];
                
					r->getHistByIdm(&hist, 0, isoFlag[v*m_u+u]-1);

					float tmpVis = 1;
					for( int i=0; i<m_bins; i++ ){
						if( i == bin || hist.getBin(i)<1 )continue;
						float r, g, b, a;
						getRGBAbyTF( (i/(float)m_bins), r, g, b, a, 1);
						a*=adtTransparent*opaAdapt;
						//a*= (1 - (fabs( (i/(float)(m_bins)) - isovalue)/width))*opaAdapt;
						tmpVis *= pow( (double)(1.0- a ), (double) hist.getBin(i)/(m_slcSample*8.0) );
						flag.setBin( i, 1 );
					}

					float rr, gg, bb, aa;
					getRGBAbyTF( isovalue, rr, gg, bb, aa, 1);
					visibility += tmpVis*aa*opaAdapt;
					visCount++;
					if( minVisibility > tmpVis ) minVisibility = tmpVis;
				}
			}
			fprintf(fp, "%f %f %f %f\n", adtTransparent, visibility, visCount, visibility/visCount );
		}while( minVisibility < visStandard );//while( visibility/visCount <= visStandard  );//;/
		
		for( int i=0; i<m_bins; i++ ){
			fprintf(fp, "%d %f\n", i, flag.getBin(i));
		}

		fclose(fp );
		free( isoFlag );

		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
				Ray* ray = m_rays[v*m_u+u];
                float cr = 0, cg = 0, cb = 0, ca = 0;

				for( int d = 0; d< m_d-step; d+=step ){
					
                    ray->getHistByIdm(&hist, d, d+step-1);
                    float vxl = hist.getMean();
					if( hist.getBin(bin) > 1 && ( msDownX<u && u<msRlsX && msDownY<v && v<msRlsY ) ) vxl = isovalue;

                    float r, g, b, a;
                    getRGBAbyTF(vxl, r, g, b, a, 1);
					if( flag.getBin( (int)(vxl*(float)m_bins) ) == 1 && ( msDownX<u && u<msRlsX && msDownY<v && v<msRlsY )  ){
						a*=adtTransparent*opaAdapt;
						//a*= (1 - (fabs( vxl - isovalue)/width))*opaAdapt;
					}
					else
						a*=opaAdapt;

					for( int i=0; i< step/8.0; i++ ){
						cr += (r*a) * (1.0 - ca);
						cg += (g*a) * (1.0 - ca);
						cb += (b*a) * (1.0 - ca);
						ca += a * (1.0 - ca);
					}

                    if( ca > 0.95 )break;
                }

				cr += (1.0) * (1.0 - ca);
                cg += (1.0) * (1.0 - ca);
                cb += (1.0) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				m_color[(v*m_u+u)*4+0] =  optR;
                m_color[(v*m_u+u)*4+1] =  optG;
                m_color[(v*m_u+u)*4+2] =  optB;
				m_color[(v*m_u+u)*4+3] =  255;//A
            }
            printf("rendering U: %d\n" ,v);
        }
	}

	//statisticRenderingVariance: render variance by step
	//transfer function is from m_tf, so call makeTrasferfunction in advance is necesary
	//step: a voxel block(slices)
    void statisticRenderingVariance( int step )
    {
        int rayCnt=0;
        int cnt=0;
        Histogram hist(m_bins);
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
                Ray* ray = m_rays[rayCnt++];
                float cr = 0, cg = 0, cb = 0, ca = 0;
                for( int d = 0; d< m_d; d+=step ){
                    ray->getHistByIdm(&hist, d, d+step-1);
					float vxl = hist.getVariance();
					vxl *= 250;
					if( vxl>=1 ){
						printf("%f\n", vxl);
						vxl = 1.0;
					}
                    float r, g, b, a;
                    getRGBAbyTF(vxl, r, g, b, a, 1);
                    a*=0.05;
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
                    if( ca > 0.95 )break;
                }
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
                m_color[cnt++] =  optR;
                m_color[cnt++] =  optG;
                m_color[cnt++] =  optB;
                //printf("%d %d %d\n", optR, optG, optB);
            }
            printf("rendering U: %d\n" ,v);
        }
    }

	//statisticRenderingVariance: render variance by step
	//transfer function is from m_tf, so call makeTrasferfunction in advance is necesary
	//step: a voxel block(slices)
    void statisticRenderingVarianceByIDR( int step, int bins, float usf )
    {
		clock_t start_time, end_time;
		float total_time = 0;

		float*** data = alloc3DMatrix( m_rsd->getDim3()*usf, usf, usf);

        int rayCnt=0;
        int cnt=0;
        Histogram hist(m_bins);
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
				//prepare upsampling ray
				//for Plume
                //for( int k=0; k<(m_rsd->getDim1()-1)*usf; k++ ){//depth
                //    for( int i=0; i<usf; i++ ){
                //        for( int j=0; j<usf; j++ ){
                //            data[k][i][j] = m_rsd->getUpsamplingVoxelValue(k, v*usf+i, u*usf+j, usf);
                //        }
                //    }
                //}
				//Ray* r = new Ray( data, m_rsd->getDim1()*usf, usf, usf, bins );

				//for Isabell
				for( int k=0; k<(m_rsd->getDim3()-1)*usf; k++ ){//depth
                    for( int i=0; i<usf; i++ ){
                        for( int j=0; j<usf; j++ ){
                            data[k][i][j] = m_rsd->getUpsamplingVoxelValue(u*usf+j, v*usf+i, k, usf);
                        }
                    }
                }

				Ray* r = new Ray( data, m_rsd->getDim3()*usf, usf, usf, bins );

                start_time = clock();
                float cr = 0, cg = 0, cb = 0, ca = 0;
                for( int d = 0; d< m_d-step; d+=step ){
					hist.reset();
					r->getHistByRaw(&hist, d, d+step-1);

					float vxl = hist.getVariance();
					vxl *= 250;
					if( vxl>=1 ){
						//printf("%f\n", vxl);
						vxl = 1.0;
					}
                    float r, g, b, a;
                    getRGBAbyTF(vxl, r, g, b, a, 1);
                    a*=0.05;
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
                    if( ca > 0.95 )break;
                }//d loop

                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
                m_color[cnt++] =  optR;
                m_color[cnt++] =  optG;
                m_color[cnt++] =  optB;
                
				end_time = clock();
				total_time += (float)(end_time - start_time)/CLOCKS_PER_SEC;

				r->releaseDnsIH();
            }//u loop
            printf("rendering U: %d\n" ,v);
        }//v loop

		printf("Time: %f sec\n", total_time);
    }

	//statisticRenderingVariance: render variance by step
	//transfer function is from m_tf, so call makeTrasferfunction in advance is necesary
	//step: a voxel block(slices)
    void statisticRenderingVarianceByRawdata( int step, int bins, float usf )
    {
		clock_t start_time, end_time;
		float total_time = 0;

        int rayCnt=0;
        int cnt=0;
        Histogram hist(m_bins);
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
                //Ray* ray = m_rays[rayCnt++];

				start_time = clock();
                float cr = 0, cg = 0, cb = 0, ca = 0;
                for( int d = 0; d< m_d-step; d+=step ){
                    //ray->getHistByIdm(&hist, d, d+step-1);
					hist.reset();
					//getUpsamplingHistFromRawData( hist, d, d+step-1, v*usf, (v+1)*usf-1, u*usf, (u+1)*usf-1, usf, bins); //for plume data
					getUpsamplingHistFromRawData( hist, u*usf, (u+1)*usf-1, v*usf, (v+1)*usf-1, d, d+step-1, usf, bins );//for Isabel data
					
					//printf("%d %d %d\n", d, (int)(v*usf),(int)(u*usf) );
					float vxl = hist.getVariance();
					vxl *= 250;
					if( vxl>=1 ){
						//printf("%f\n", vxl);
						vxl = 1.0;
					}
                    float r, g, b, a;
                    getRGBAbyTF(vxl, r, g, b, a, 1);
                    a*=0.05;
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
                    if( ca > 0.95 )break;
                }

                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
                m_color[cnt++] =  optR;
                m_color[cnt++] =  optG;
                m_color[cnt++] =  optB;

				end_time = clock();
				total_time += (float)(end_time - start_time)/CLOCKS_PER_SEC;
                //printf("%d %d %d\n", optR, optG, optB);
            }
            printf("rendering U: %d\n" ,v);
        }
		printf("Time: %f sec\n", total_time);
    }

	void fitStorageMultiResolution()
	{
		int segLimit = m_u*m_v*4;

		FILE* fp = fopen( "tmpOut.txt", "w" );

		vector<Segment> segs;
		Histogram hist(m_bins);
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
				Segment seg;
				seg.pxl = v*m_u+u;
				seg.from = 0;
				seg.len = m_d - seg.from;

				Ray* ray = m_rays[v*m_u+u];
                hist.reset();
				ray->getHistByIdm( &hist, seg.from, seg.from+seg.len-1 );
				hist.normalize();
				//seg.metric = hist.getEntropy();
				seg.metric = hist.getVariance();

				segs.push_back( seg );
			}
		}
		std::sort( segs.begin(), segs.end(), segmentSort );
		
		int numSeg = m_u*m_v;
		int minLen = 8;
		while( segLimit > numSeg ){
			int spIdx = segs.size()-1;
			while( segs[spIdx].len <= minLen )spIdx--;
			Segment seg = segs[spIdx];
			segs.erase( segs.begin()+spIdx );
			//here, miss a check, if the segment with largest entropy is too small to divide,
			//we should choose another one, but this should not happen

			std::vector<Segment>::iterator position;
			Segment seg1;
			Segment seg2;
			int l = (int)(seg.len/2.0);

			seg1.pxl = seg.pxl;
			seg1.from = seg.from;
			seg1.len = l;
			hist.reset();
			m_rays[seg1.pxl]->getHistByIdm( &hist, seg1.from, seg1.from+seg1.len-1 );
			hist.normalize();
			//seg1.metric = hist.getEntropy();
			seg1.metric = hist.getVariance();
			position = std::lower_bound (segs.begin(), segs.end(), seg1, segmentSort );
			segs.insert( position, seg1 );

			seg2.pxl = seg.pxl;
			seg2.from = seg.from+l;
			seg2.len = seg.len-l;
			hist.reset();
			m_rays[seg2.pxl]->getHistByIdm( &hist, seg2.from, seg2.from+seg2.len-1 );
			hist.normalize();
			//seg2.metric = hist.getEntropy();
			seg2.metric = hist.getVariance();
			position = std::lower_bound (segs.begin(), segs.end(), seg2, segmentSort );
			segs.insert( position, seg2 );

			numSeg++;
		}

		//reorgnize the data structure
		vector<Segment> **img;
		img = new vector<Segment>*[m_u*m_v];
		for( int i =0 ; i < m_u*m_v ; i++ )
			img[i] = new std::vector<Segment>;
		std::vector<Segment>::iterator position;
		for( int i=0; i<segs.size(); i++ ){
			int pxl = segs[i].pxl;
			position = std::lower_bound ( img[pxl]->begin(), img[pxl]->end(), segs[i], segmentSortByFrom );
			img[pxl]->insert( position, segs[i] );
		}

		float maxSegs = 0;
		for( int i=0; i<m_u*m_v; i++ ){
			if( maxSegs < img[i]->size() )
				maxSegs = img[i]->size();
		}

		//rendering
		#pragma omp parallel for
		for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
            for( int u = 0; u< m_u; u+=1 ){
				float cr = 0, cg = 0, cb = 0, ca = 0;
				int idx = v*m_u+u;
				Ray* r = m_rays[idx];

				for( int d=0; d < img[idx]->size(); d++ ){
					vector<Segment>::iterator t = img[idx]->begin()+d;
					int from  = (*t).from;
					int len = (*t).len;

					r->getHistByIdm( &hist, from, from+len-1 );
					hist.normalize();
					float vxl = hist.getMean();

					float r, g, b, a;
					float sr, sg, sb, sa;
					getRGBAbyTF( vxl, r, g, b, a, 1 );
					a*=0.1;
					/*renderSameVoxel( sr, sg, sb, sa, r, g, b, a, len/8.0);
					cr += sr*sa*(1-ca);
					cg += sg*sa*(1-ca);
					cb += sb*sa*(1-ca);
					ca += sa * (1-ca);*/
					for( int i=0; i<len/8.0; i++ ){
						cr += (r*a) * (1.0 - ca);
						cg += (g*a) * (1.0 - ca);
						cb += (b*a) * (1.0 - ca);
						ca += a * (1.0 - ca);
					}
					if( ca >= 0.95 )break;
				}

				float bgColor = 1.0;
				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
				unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
				unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
				unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(idx)*4+0] = optR;
				m_color[(idx)*4+1] = optG;
				m_color[(idx)*4+2] = optB;
				m_color[(idx)*4+3] = 255;

				float rr, gg, bb, aa;
				getRGBAbyTF( img[idx]->size()/maxSegs, rr, gg, bb, aa, 1 );
				m_color2[(idx)*4+0] = (unsigned char)(int)((rr)*255.0);
				m_color2[(idx)*4+1] = (unsigned char)(int)((gg)*255.0);
				m_color2[(idx)*4+2] = (unsigned char)(int)((bb)*255.0);
				m_color2[(idx)*4+3] = 255;

            }//u loop
            printf("rendering U: %d\n" ,v);
        }//v loop


		fprintf(fp, "%f\n", maxSegs);
		for( int i=0; i<m_u*m_v; i++ ){
			for( int j=0; j<img[i]->size(); j++ ){
				vector<Segment>::iterator t = img[i]->begin()+j;
				fprintf( fp, "%d %d %d %f\n", (*t).pxl, (*t).from, (*t).len, (*t).metric  );
			}
			fprintf(fp, "\n");
		}


		//for( int i=0; i<segs.size(); i++ )
		//	fprintf( fp, "%f %d %d\n", segs[i].metric, segs[i].len, segs[i].from  );
		fclose(fp);
	}


	//*********************statistic rendering - variance -END**********************

	void calculateVisibilityHistogram( Histogram* optVisHist, float colorModulation, int step )
	{
		optVisHist->reset();
		//not first run
		#pragma omp parallel for
		for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram rayVisHist(m_bins);	//store vis histogram of a local segment
			Histogram rayVisHist2(m_bins);
			for( int u = 0; u< m_u; u+=1 ){
				int rayIdx = v*m_u+u;
				Ray* ray = m_rays[rayIdx];
				rayVisHist.reset();
				rayVisHist2.reset();

				float alpha = 0;	//ray alpha at certain point in d loop
				for( int d = 0; d < m_d - step; d += step ) {
					ray->getHistByIdm( &hist, d, d + step - 1 );
					

					float transparancy = 1.0;	//transparancy of a segment
					float sumOpacity = 0.0;	//sum of opacity of a segment
					
					//calculate transparrancy and sumPacity of a segment
					for( int i=0; i<m_bins; i++ ){
						float r, g, b, a;
						getRGBAbyTF( (i/(float)m_bins), r, g, b, a, 1);
						a *= colorModulation;
						float freq = hist.getBin(i);
						//freq /= 64.0;	//****if the sample number is too large, freq need to be reduce, otherwise the transparancy calculation will overflow
						transparancy *= pow( (double)(1.0 - a ), (double) freq );
						sumOpacity += ( a * freq );
						rayVisHist.setBin( i, a * freq );//visibility before modulate by normalizing to (1.0-transparancy)
					}

					//if this segment is totally transparancy, do run this block
					if( sumOpacity > 0.00000001 ){
						float visModulation = ( ( 1.0 - transparancy ) / sumOpacity ) * ( 1.0 - alpha );	//modulate visibility histogram
						alpha += ( 1.0 - alpha ) * ( 1.0 - transparancy );	//accumlate alpha along this ray

						//modulate visibility histogram
						for( int i = 0; i < m_bins; i++ ){
							//rayVisHist.setBin( i, rayVisHist.getBin( i ) * visModulation+1 );
							rayVisHist2.setBin( i, rayVisHist2.getBin( i ) + rayVisHist.getBin( i ) * visModulation );
						}
					}
					

				}//for d

	
				//sum the visHistogram to optVisHistogram
				#pragma omp critical
				{
					for( int i=0; i<m_bins; i++ ){
						optVisHist->setBin( i, optVisHist->getBin( i ) + rayVisHist2.getBin( i ) );
					}
				}

			}//u loop
		}//v loop
		FILE* fp = fopen( "visHistOutput.txt", "w" );
		for( int i=0; i<128; i++ )
			fprintf( fp, "%f\n", optVisHist->getBin( i ) );
		fclose( fp );


		//**********************************************************
		optVisHist->reset();
		//not first run
		#pragma omp parallel for
		for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram rayVisHist(m_bins);	//store vis histogram of a local segment
			Histogram rayVisHist2(m_bins);
			for( int u = 0; u< m_u; u+=1 ){
				int rayIdx = v*m_u+u;
				Ray* ray = m_rays[rayIdx];
				rayVisHist.reset();
				rayVisHist2.reset();

				float alpha = 0;	//ray alpha at certain point in d loop
				for( int d = 0; d < ray->getIdmIHSize(); d ++ ) {
					int trend = ray->getSingleIdmHistogram( &hist, d );
														
					if( trend == 1 ){
						for( int i=0; i<m_bins ;i++ ){
							float r, g, b, a;
							getRGBAbyTF( (i/(float)m_bins), r, g, b, a, 1);
							a *= colorModulation;
							float freq = hist.getBin(i);

							for( int c = 0; c< freq; c++ ){
								float vis = ( 1.0 - alpha ) * a;
								alpha += vis;
								rayVisHist.setBin( i, rayVisHist.getBin(i) + vis );
							}
						}
					}else{
						for( int i=m_bins-1; i>=0 ;i-- ){
							float r, g, b, a;
							getRGBAbyTF( (i/(float)m_bins), r, g, b, a, 1);
							a *= colorModulation;
							float freq = hist.getBin(i);

							for( int c = 0; c< freq; c++ ){
								float vis = ( 1.0 - alpha ) * a;
								alpha += vis;
								rayVisHist.setBin( i, rayVisHist.getBin(i) + vis );
							}
						}
					}

				}//for d

	
				//sum the visHistogram to optVisHistogram
				#pragma omp critical
				{
					for( int i=0; i<m_bins; i++ ){
						optVisHist->setBin( i, optVisHist->getBin( i ) + rayVisHist.getBin( i ) );
					}
				}

			}//u loop
		}//v loop
		fp = fopen( "visHistOutput2.txt", "w" );
		for( int i=0; i<128; i++ )
			fprintf( fp, "%f\n", optVisHist->getBin( i ) );
		fclose( fp );
	}

	void depthVisibilityVisualization()
	{   
		//initialize the UI inputing parameters
		float colorModulation = m_lineEdit0;
		int depth = ( m_sldBar9 / 100.0 ) * m_d;	//render visibility from this depth
		int step = m_lineEdit3;		//render visiblity of this thickness

		//not first run
		#pragma omp parallel for
		for( int v = 0; v< m_v; v+=1 ){
			Histogram hist(m_bins);
			Histogram befHist( m_bins );
			for( int u = 0; u< m_u; u+=1 ){
				int rayIdx = v*m_u+u;
				Ray* ray = m_rays[rayIdx];

				//calculate transparrancy and sumPacity of a segment
				float cr = 0, cg = 0, cb = 0, ca = 0;	//color of a pixel
				float sumOfEntropy = 0;

				float frontSegTransparency, correctOpacity;

				//calculate the front histogrtam and histogram of this semgent
				float d = depth;
				if( d>0 )
					ray->getHistByIdm( &befHist, 0, d-1);
				else
					befHist.reset();
				ray->getHistByIdm( &hist, d, d + step - 1 );

				//use the two historgram to calcualte the real color and opacity and visibility contribution
				float rr, gg, bb;
				calculateSegmentColorContribution( &befHist, &hist, colorModulation, rr, gg, bb, frontSegTransparency, correctOpacity );

				ca = ( 1 - frontSegTransparency ) + correctOpacity * frontSegTransparency; //not necessary

				//render color contribution of this depth segment
				cr = rr * correctOpacity * frontSegTransparency;
				cg = gg * correctOpacity * frontSegTransparency;
				cb = bb * correctOpacity * frontSegTransparency;
				unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
				unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
				unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				m_color[(v*m_u+u)*4+0] =  optR;
				m_color[(v*m_u+u)*4+1] =  optG;
				m_color[(v*m_u+u)*4+2] =  optB;
				m_color[(v*m_u+u)*4+3] =  255;//A
			}//u loop
		}//v loop
		

    }

	//***********************ABOUT TRANSFER FUCNTION************************

    //makeTransferFunction: given control point make the transfer function to m_tf
    //nCtrlPoint: number of control point
    //ctrlPoint: array of contrlpoint. (if nCtrlPoint is n,
    //size of ctrlPoint is 5(rgba)*n
    //every 5-tuple is a (control point) rgbaVoxel
    //and first voxel in contraol point must be 0
    //last voxel in control point array must 1
    //this function assume, voxels of any two neighbor contoal point > 1/1024
	//base : base==1 use m_bsTF    0 use m_updTF
    void makeTransferFunction( int nCtrlPoint, float* ctrlPoint, int base )
    {
        //1024 interval so 1025 points
		if( base == 1 ){
			int cpNow = 0;
			for( int i=0; i<1025; i++ ){
				float vxl = i/1024.0;
				//current vxl > vxl in current control point, move to next control point
				if( vxl > ctrlPoint[5*(cpNow+1)+4] )cpNow++;
				m_bsTF[i][0] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
													ctrlPoint[5*cpNow+0], ctrlPoint[5*(cpNow+1)+0]);//r
				m_bsTF[i][1] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
												ctrlPoint[5*cpNow+1], ctrlPoint[5*(cpNow+1)+1]);//g
				m_bsTF[i][2] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
											 ctrlPoint[5*cpNow+2], ctrlPoint[5*(cpNow+1)+2]);//b
				m_bsTF[i][3] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
											 ctrlPoint[5*cpNow+3], ctrlPoint[5*(cpNow+1)+3]);//a
			}
		}else{
			int cpNow = 0;
			for( int i=0; i<1025; i++ ){
				float vxl = i/1024.0;
				//current vxl > vxl in current control point, move to next control point
				if( vxl > ctrlPoint[5*(cpNow+1)+4] )cpNow++;
				m_updTF[i][0] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
													ctrlPoint[5*cpNow+0], ctrlPoint[5*(cpNow+1)+0]);//r
				m_updTF[i][1] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
												ctrlPoint[5*cpNow+1], ctrlPoint[5*(cpNow+1)+1]);//g
				m_updTF[i][2] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
											 ctrlPoint[5*cpNow+2], ctrlPoint[5*(cpNow+1)+2]);//b
				m_updTF[i][3] = interpolation1D(vxl, ctrlPoint[5*cpNow+4], ctrlPoint[5*(cpNow+1)+4],
											 ctrlPoint[5*cpNow+3], ctrlPoint[5*(cpNow+1)+3]);//a
			}
		}
    }
    
	//makeUniformTransferFunction: given control point make the transfer function to m_tf, 
	//but the color and opacity of any interval is uniform
    //nCtrlPoint: number of control point
    //ctrlPoint: array of contrlpoint. (if nCtrlPoint is n,
    //size of ctrlPoint is 5(rgba)*n
    //every 5-tuple is a (control point) rgbaVoxel
    //and first voxel in contraol point must be 0
    //last voxel in control point array must 1
    //this function assume, voxels of any two neighbor contoal point > 1/1024
	//base : base==1 use m_bsTF    0 use m_updTF
    void makeUniformTransferFunction( int nCtrlPoint, float* ctrlPoint, int base )
    {
        //1024 interval so 1025 points
		if( base == 1 ){
			int cpNow = 0;
			for( int i=0; i<1025; i++ ){
				float vxl = i/1024.0;
				//current vxl > vxl in current control point, move to next control point
				if( vxl > ctrlPoint[5*(cpNow+1)+4] )cpNow++;
				m_bsTF[i][0] = ctrlPoint[5*cpNow+0];//r
				m_bsTF[i][1] = ctrlPoint[5*cpNow+1];//g
				m_bsTF[i][2] = ctrlPoint[5*cpNow+2];//b
				m_bsTF[i][3] = ctrlPoint[5*cpNow+3];//a
				if( i<200)printf( "%d  %f %f %f %f\n",i, m_bsTF[i][0], m_bsTF[i][1], m_bsTF[i][2], m_bsTF[i][3] );
			}
		}else{
			int cpNow = 0;
			for( int i=0; i<1025; i++ ){
				float vxl = i/1024.0;
				//current vxl > vxl in current control point, move to next control point
				if( vxl > ctrlPoint[5*(cpNow+1)+4] )cpNow++;
				m_updTF[i][0] = ctrlPoint[5*cpNow+0];//r
				m_updTF[i][1] = ctrlPoint[5*cpNow+1];//g
				m_updTF[i][2] = ctrlPoint[5*cpNow+2];//b
				m_updTF[i][3] = ctrlPoint[5*cpNow+3];//a
			}
		}
    }

	//getRGBAbyTF: after make transfer fucntion: get RGBA by voxel
    //vxl: given voxel ( must between 0-1)
    //rgba: return value of RGBA
	//base : base==1 use m_bsTF    0 use m_updTF
    void getRGBAbyTF( float vxl, float &r, float &g, float &b, float& a, int base  )
    {
        int idx = (int)(vxl*1024.0);
		if( base == 1 ){
			r = m_bsTF[idx][0];
			g = m_bsTF[idx][1];
			b = m_bsTF[idx][2];
			a = m_bsTF[idx][3];
		}else{
			r = m_updTF[idx][0];
			g = m_updTF[idx][1];
			b = m_updTF[idx][2];
			a = m_updTF[idx][3];
		}
    }

	void getRGBAbyHistogram( Histogram* hist, int trend, float slcSample, float &hr, float &hg, float &hb, float &ha, int tf )
	{
		hr = 0;
		hg = 0;
		hb = 0;
		ha = 0;
		int start, iter;
		if( trend == -1){
			start = m_bins-1;
			iter = -1;
		}else{
			start = 0;
			iter = 1;	
		}

		for( int i=start; i>=0 && i<m_bins; i+=iter ){
			float x = ceil((hist->getBin(i)/slcSample));//+0.5 is for round operation
			if( x!= 0 ){
				float r, g, b,a;
				getRGBAbyTF( i/(float)m_bins, r, g, b, a, tf);
				float br = 0, bg = 0, bb = 0, ba = 0;
				float aa=a;
				a*=0.005;
				if( a == 0 ) 
					ba = 0;
				else
					ba = (a - a*pow((1-a),x) )/(1-(1-a));

				br = r*ba;
				bg = g*ba;
				bb = b*ba;

				hr+= (1-ha)*br;
				hg+= (1-ha)*bg;
				hb+= (1-ha)*bb;
				ha += (1-ha)*ba;
			}
		}
	}
	//***********************ABOUT TRANSFER FUCNTION -END************************


	//********************OTHERS************************************

	//wirteToFile: output m_rays(histogram to a file)
    void writeToFile( char* pathName )
    {
        FILE*fp = fopen( pathName, "wb");
        for( int i=0; i<m_rays.size(); i++ ){
            m_rays[i]->writeToFile(fp);
            float tmp = -1;
            fwrite( &tmp, sizeof(float), 1, fp);
        }
        float tmp = -2;
        fwrite( &tmp, sizeof(float), 1, fp);
        fclose(fp);
    }

	void outputPPMImage( char* pathName )
    {
		unsigned char* data = (unsigned char*) malloc( sizeof(unsigned char) * m_u * m_v * 3 );
		for( int i = 0; i < m_u*m_v; i ++ ){
			data[ i * 3 + 0 ] = m_color[ i * 4 + 0 ];
			data[ i * 3 + 1 ] = m_color[ i * 4 + 1 ];
			data[ i * 3 + 2 ] = m_color[ i * 4 + 2 ];
		}

        FILE* f0 = fopen( pathName, "wb" );
        fprintf( f0, "P6\n%d %d\n255\n", m_u,m_v);
        fwrite( data, sizeof(unsigned char), m_u*m_v* 3, f0);
        fclose( f0 );

		free( data );
    }

	void outputPPMImage2( char* pathName )
    {
		unsigned char* data = (unsigned char*) malloc( sizeof(unsigned char) * m_u * m_v * 3 );
		for( int i = 0; i < m_u*m_v; i ++ ){
			data[ i * 3 + 0 ] = m_color2[ i * 4 + 0 ];
			data[ i * 3 + 1 ] = m_color2[ i * 4 + 1 ];
			data[ i * 3 + 2 ] = m_color2[ i * 4 + 2 ];
		}

        FILE* f0 = fopen( pathName, "wb" );
        fprintf( f0, "P6\n%d %d\n255\n", m_u,m_v);
        fwrite( data, sizeof(unsigned char), m_u*m_v* 3, f0);
        fclose( f0 );

		free( data );
    }

	unsigned char* getImgColorPointer()
	{
		return m_color;
	}

	unsigned char* getImgColorPointer2()
	{
		return m_color2;
	}

	//
	void getUpsamplingHistFromRawData( Histogram &hist, int x1, int x2, int y1, int y2, int z1, int z2,
										float usf, int bins )
	{
		for( int x = x1; x<=x2; x++ ){
			for( int y = y1; y<=y2; y++ ){
				for( int z = z1; z<=z2; z++ ){
					float vxl = m_rsd->getUpsamplingVoxelValue( x, y, z, usf );
					int b = (int)(vxl*bins);
					hist.setBin( b, hist.getBin(b)+1 );
				}
			}
		}
	}

	//n is the number of samples on this ray
	void renderSameVoxel(   float& ar, float &ag, float &ab, float &aa,
							float r, float g, float b, float a, float n )
	{
		if( a == 0 ){
			ar = 0;
			ag = 0;
			ab = 0;
			aa = 0;
		}else{
			ar = r;
			ag = g;
			ab = b;
			aa = (a - a*pow((1-a),n) )/(1-(1-a));
		}
	}

	void setClrAdapt( float clrAdapt )
	{
		m_clrAdapt = clrAdapt;
	}

	void setOrder(float order )
	{
		m_order = order;
	}

	void setBackgroundColor( float color )
	{
		m_bgColor = color;
	}

	void calculateDistributionQuerySinglePixel( int luu, int luv, float s, float t)
	{ 
		int rayCnt=0;
        int cnt=0;
		float maxEn = 0;
		float minEn = 10000;

		Ray* ray = m_rays[luv*m_u+luu];
		Histogram hist(m_bins);
				
		int start = (int)(m_d-1)*s;
		int til = (int)(m_d-1)*t;
		
		ray->getHistByIdm(&hist, start, til);

		hist.normalize();

		//output file
		FILE *fp; 
		 do{
			fp =  fopen( "C:\\GravityLabDataSet\\DistributionQuery2.txt", "w" );
		}while( fp == NULL );
		 for( int i=0; i<m_bins; i++ ){
			 fprintf( fp, "%f\n", hist.getBin(i) );
		 }
		 fclose(fp);
	}

	//lux, luy, rub, rby are left upper xy,and right bottom xy, s t is the depth
	void calculateDistributionQuery(  Histogram* reHist, int luu, int luv, int rbu, int rbv, float s, float t)
	{
		int rayCnt=0;
        int cnt=0;
		float maxEn = 0;
		float minEn = 10000;
		float* sumHist = (float*)malloc( sizeof(float)*m_bins );
		memset( sumHist, 0, sizeof(float)*m_bins );

		//for( int v = 0; v< m_v; v+=1 ){
		//	Histogram hist(m_bins);
		//	for( int u = 0; u< m_u; u+=1 ){
		//		if( u<luu || u>rbu || v<luv || v>rbv )continue;
		//		Ray* ray = m_rays[v*m_u+u];
		//	
		//		int start = (int)(m_d-1)*s;
		//		int til = (int)(m_d-1)*t;
		//
		//		//printf("%d %d\n", start, til);
		//		ray->getHistByIdm(&hist, start, til);
		//		for( int i=0; i<m_bins; i++ ){
		//			sumHist[i] += hist.getBin(i);
		//		}
		//	}
		//}

        #pragma omp parallel
		{
			float* privateHist = (float*) malloc( sizeof(float)*m_bins );
			memset( privateHist, 0, sizeof(float)*m_bins );

			#pragma omp for
			for( int v = 0; v< m_v; v+=1 ){
				Histogram hist(m_bins);
				for( int u = 0; u< m_u; u+=1 ){
					if( u<luu || u>rbu || v<luv || v>rbv )continue;
					Ray* ray = m_rays[v*m_u+u];
				
					int start = (int)(m_d-1)*s;
					int til = (int)(m_d-1)*t;

					if( m_order == 1 ){//reverse order
						start = (m_d-1)-start;
						til = (m_d-1)-til;
						int tmp;
						tmp = start;
						start = til;
						til = tmp;
					}
					//printf("%d %d\n", start, til);
					ray->getHistByIdm(&hist, start, til);
					for( int i=0; i<m_bins; i++ ){
						privateHist[i] += hist.getBin(i);
					}
				}
				printf("V %d\n", v);
			}

			#pragma omp critical
			{
				for( int i=0; i<m_bins; i++ ){
					sumHist[i] += privateHist[i];
				}
			}

			free(privateHist );
		}
		
		
		for( int i=0; i<m_bins; i++ )
			reHist->setBin( i, sumHist[i] );
		reHist->normalize();

		free(sumHist);
	}

	void getRayTrend(  Histogram* reHist, int bins, int luu, int luv )
	{
		int rayCnt=0;
		float maxEn = 0;
		float minEn = 10000;
		float* sumHist = (float*)malloc( sizeof(float)*bins );
		memset( sumHist, 0, sizeof(float)*bins );

		int step = 4;
		Ray* ray = m_rays[luv*m_u+luu];
		Histogram hist(m_bins);
		int cnt=0;
		for( int d = 0; d< m_d-step; d+=step ){
			ray->getHistByIdm(&hist, d, d+step-1);
            float vxl = hist.getMean();

			sumHist[cnt++]=vxl;
        }
		
		ray->printIdmData();

		for( int i=0; i<bins; i++ )
			reHist->setBin( i, sumHist[i] );

		free(sumHist);
	}

	void setUIParameters( int sldBar0, int sldBar1, int sldBar2, int sldBar3, int sldBar4, 
						  int sldBar5, int sldBar6, int sldBar7, int sldBar8, int sldBar9,
						  float lineEdit0, float lineEdit1, float lineEdit2, float lineEdit3, float lineEdit4, 
						  float lineEdit5, float lineEdit6, float lineEdit7, float lineEdit8, float lineEdit9 )
	{
		m_sldBar0 = sldBar0;  m_sldBar1 = sldBar1;  m_sldBar2 = sldBar2;  m_sldBar3 = sldBar3;  m_sldBar4 = sldBar4;
		m_sldBar5 = sldBar5;  m_sldBar6 = sldBar6;  m_sldBar7 = sldBar7;  m_sldBar8 = sldBar8;  m_sldBar9 = sldBar9;
		m_lineEdit0 = lineEdit0;  m_lineEdit1 = lineEdit1;  m_lineEdit2 = lineEdit2;  m_lineEdit3 = lineEdit3;  m_lineEdit4 = lineEdit4;  
		m_lineEdit5 = lineEdit5;  m_lineEdit6 = lineEdit6;  m_lineEdit7 = lineEdit7;  m_lineEdit8 = lineEdit8;  m_lineEdit9 = lineEdit9;
	}

	//calculate the energy absorbed by this histogram of segment (return value is normalized to 0-1)
	float calculateSegmentOpacity( Histogram* hist, float colorModulation )
	{
		float transparancy = 1.0;
		for( int i=0; i<m_bins; i++ ){
			float r, g, b, a;
			getRGBAbyTF( (i/(float)m_bins), r, g, b, a, 1);
			a *= colorModulation;
			float freq = hist->getBin(i);
			transparancy *= pow( (double)(1.0- a ), (double) freq );
		}
		return ( 1.0 - transparancy );
	}
	//********************OTHERS - END************************************

private:
    //interpolation 1D; do 1d intepolation
    //return value
    //p: point we want, ps: point start, pt: point util
    //vs: value start, vt: value til
    //for example:  interpolation( 4, 0 10, 0, 100 ) return 40
    float interpolation1D( float p, float ps, float pt, float vs, float vt)
    {
        float ratio = (p-ps)/(pt-ps);
        return ((vt-vs)*ratio + vs);
    }

	void crossProduct( float vR[], float v1[], float v2[] ) 
	{
		vR[0] =   ( (v1[1] * v2[2]) - (v1[2] * v2[1]) );
		vR[1] = - ( (v1[0] * v2[2]) - (v1[2] * v2[0]) );
		vR[2] =   ( (v1[0] * v2[1]) - (v1[1] * v2[0]) );
	}

	void normalize(float vR[], float v1[]) {
		float fMag;

		fMag = sqrt( pow(v1[0], 2) +
				   pow(v1[1], 2) +
				   pow(v1[2], 2)
				);

		vR[0] = v1[0] / fMag;
		vR[1] = v1[1] / fMag;
		vR[2] = v1[2] / fMag;
	}
    
    RegularScalarData* m_rsd;//regular scalar data object
    vector<Ray*> m_rays;	//rays
    float m_bsTF[1025][4];//transfer function table old
	float m_updTF[1025][4];//transfer function table new
    unsigned char* m_color;	//image color buffer
	float* m_colorF;	//image color buffer
	unsigned char* m_color2;	//image color buffer for draw
	float* m_color2F;	//image color buffer
	unsigned char* m_color3;	//image color buffer for draw
	float* m_depth;	//temp area for any need
    int m_u, m_v, m_d;//image dimension;
    int m_bins;	//number of bins
	int m_slcSample;
	int* m_ddaX;//x component of dda algorithm
	int* m_ddaY;
	int m_ddaCnt;
	float m_clrAdapt;
	float m_order;//0 is order(plume), 1 is reverse order for isabell
	float m_bgColor;	//back ground color 0-1
	vector<vector<float>> m_mrrV; //v:voxel, E: entropy, L, location
	vector<vector<float>> m_mrrE;
	vector<vector<float>> m_mrrL;
	int m_sldBar0;
	int m_sldBar1, m_sldBar2, m_sldBar3, m_sldBar4, 
		m_sldBar5, m_sldBar6, m_sldBar7, m_sldBar8, m_sldBar9;
	float m_lineEdit0, m_lineEdit1, m_lineEdit2, m_lineEdit3, m_lineEdit4, 
		  m_lineEdit5, m_lineEdit6, m_lineEdit7, m_lineEdit8, m_lineEdit9;
	vector<Segment> **entropyQueue;
};


#endif

