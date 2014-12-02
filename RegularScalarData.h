#ifndef _REGULARSCALARDATA_H_
#define _REGULARSCALARDATA_H_

#define E 2.7182818284

#include "stdio.h"
#include "math.h"
//#include "conio.h"
#include "stdlib.h"


#include "Histogram.h"
#include "Ray.h"

enum VolumeDataType
{
	RAW_DATA_UNSIGNEDCHAR,
	//raw data, char format, only tested for Bucky
	RAW_DATA_FLOAT,
	//raw data, float format, only tested for Isabala data
	VECTOR_MAGNITUDE,	
	//vector data, but create magnitude as scalar in class, mainly and only tested for plume .vec
};

class RegularScalarData{
public:
	RegularScalarData( char* fnPath, int dataType, int dim1, int dim2, int dim3 )
	{	
		m_dim1 = dim1;
		m_dim2 = dim2;
		m_dim3 = dim3;
		m_size = m_dim1* m_dim2 * m_dim3;
		m_dataType = dataType;
		
		if( dataType == VolumeDataType::RAW_DATA_UNSIGNEDCHAR ){
			loadRawDataChar(fnPath);
		}else if( dataType == VolumeDataType::RAW_DATA_FLOAT ){
			loadRawDataFloat(fnPath);
		}else if( dataType == VolumeDataType::VECTOR_MAGNITUDE ){
			loadVectorData( fnPath );
		}
	}

	~RegularScalarData()
	{
	
	}

	//***************************************
	// testRendering
	// only provide a simeple view to test the loading is correct or not
	// attFactor: attenuation factor, less indication more transparent
	//***************************************
	void testRendering( float attFactor = 0.2 )
	{
		//first view and image
		unsigned char* pxlBuffer = (unsigned char*) malloc( m_dim1*m_dim2*3*sizeof(unsigned char) );
		int cnt=0;
		for( int d1 = 0; d1< m_dim1; d1++ ){
			for( int d2 = 0; d2< m_dim2; d2++ ){
				float color = 0;
				float alpha = 0;
				for( int d3 = 0; d3< m_dim3; d3++ ){
					int idx = idx3dToIdx1d(d1, d2, d3);
					float vxl = m_data[idx];
					float opq = vxl * attFactor;
					color += (vxl*opq) * (1.0 - alpha);
					alpha += opq * (1.0 - alpha); 
				}
				unsigned char c = (unsigned char)(int)floor(color * 255.0);
				pxlBuffer[cnt++] =  c;
				pxlBuffer[cnt++] =  c;
				pxlBuffer[cnt++] =  c;
			}
		}
		FILE *f0 = fopen( "Dim1Dim2View.ppm", "wb" );    
		fprintf( f0, "P6\n%d %d\n255\n", m_dim2, m_dim1);
		fwrite( pxlBuffer, sizeof(unsigned char), m_dim1*m_dim2* 3, f0);
		fclose( f0 );
		free( pxlBuffer);
		
		//second view and image
		pxlBuffer = (unsigned char*) malloc( m_dim3*m_dim2*3*sizeof(unsigned char) );
		cnt=0;
		for( int d2 = 0; d2< m_dim2; d2++ ){
			for( int d3 = 0; d3< m_dim3; d3++ ){
				float color = 0;
				float alpha = 0;
				for( int d1 = 0; d1< m_dim1; d1++ ){
					int idx = idx3dToIdx1d(d1, d2, d3);
					float vxl = m_data[idx];
					float opq = vxl * attFactor;
					color += (vxl*opq) * (1.0 - alpha);
					alpha += opq * (1.0 - alpha); 
				}
				unsigned char c = (unsigned char)(int)floor(color * 255.0);
				pxlBuffer[cnt++] =  c;
				pxlBuffer[cnt++] =  c;
				pxlBuffer[cnt++] =  c;
			}
		}		
		f0 = fopen( "Dim2Dim3View.ppm", "wb" );  
		fprintf( f0, "P6\n%d %d\n255\n", m_dim3, m_dim2);
		fwrite( pxlBuffer, sizeof(unsigned char), m_dim3*m_dim2* 3, f0);
		fclose( f0 );
		free( pxlBuffer);
		
		//third view and image
		pxlBuffer = (unsigned char*) malloc( m_dim1*m_dim3*3*sizeof(unsigned char) );
		cnt=0;
		for( int d3 = 0; d3< m_dim3; d3++ ){
			for( int d1 = 0; d1< m_dim1; d1++ ){
				float color = 0;
				float alpha = 0;
				for( int d2 = 0; d2< m_dim2; d2++ ){
					int idx = idx3dToIdx1d(d1, d2, d3);
					float vxl = m_data[idx];
					//printf("%d\n", idx);
					float opq = vxl * attFactor;
					color += (vxl*opq) * (1.0 - alpha);
					alpha += opq * (1.0 - alpha); 
				}
				unsigned char c = (unsigned char)(int)floor(color * 255.0);
				pxlBuffer[cnt++] =  c;
				pxlBuffer[cnt++] =  c;
				pxlBuffer[cnt++] =  c;
			}
		}		
		f0 = fopen( "Dim3Dim1View.ppm", "wb" );
		fprintf( f0, "P6\n%d %d\n255\n", m_dim1, m_dim3);
		fwrite( pxlBuffer, sizeof(unsigned char), m_dim1*m_dim3* 3, f0);
		fclose( f0 );
		free( pxlBuffer);
	}



	float getVoxelValue(int idx1, int idx2, int idx3 )
	{
		return m_data[ idx3dToIdx1d( idx1, idx2, idx3 )];
	}
    
    //getUpscaleingVoxelValue: get vovel value after upsscaleing
    //this method do interpolation for every requesting, so inefficient
    //idx123: the index after upsampling
    //upsamplescacle: upsample scale
    //do linear interpolation
    float getUpsamplingVoxelValue( int idx1, int idx2, int idx3, int upSampleScale )
    {
        int i1 = idx1/upSampleScale;
        int i2 = idx2/upSampleScale;
        int i3 = idx3/upSampleScale;
        
        float c1 = (idx1 - i1*upSampleScale )/(float)upSampleScale;
        float c2 = (idx2 - i2*upSampleScale )/(float)upSampleScale;
        float c3 = (idx3 - i3*upSampleScale )/(float)upSampleScale;
        return trilinear(c1, c2, c3,
                         getVoxelValue(i1, i2, i3), getVoxelValue(i1, i2, i3+1), getVoxelValue(i1, i2+1, i3), getVoxelValue(i1, i2+1, i3+1),
                         getVoxelValue(i1+1, i2, i3), getVoxelValue(i1+1, i2, i3+1), getVoxelValue(i1+1, i2+1, i3), getVoxelValue(i1+1, i2+1, i3+1) );
    }
    
    int getDim1()
    {
        return m_dim1;
    }
    
    int getDim2()
    {
        return m_dim2;
    }
    
    int getDim3()
    {
        return m_dim3;
    }


private:
	//***************************************
	// loadRawDataChar
	// binary data and element should be float
	//***************************************
	void loadRawDataChar( char* fnPath )
	{
		FILE *fp = fopen(fnPath, "rb");
		if (!fp)
		{
			fprintf(stderr, "Error opening file '%s'\n", fnPath);
		}else{
			m_data = (float*)malloc(m_size*sizeof(float) );
			unsigned char* chData = (unsigned char*)malloc(m_size*sizeof(unsigned char) );
			fread((void*)chData, 1, m_size*sizeof(char), fp);
			fclose(fp);

			float min = 100000;
			float max = -100000;
			
			for( int i=0; i<m_size; i++ ){
				m_data[i] = (float)(int)chData[i];
				if( min > m_data[i] )min = m_data[i];
				if( max < m_data[i] )max = m_data[i];
			}
			printf("Raw data: max(%f)    min(%f)\n" , max, min);

			for( int i=0; i<m_size; i++ ){
				m_data[i] = ( m_data[i] - min )/ (max - min);
			}

			free( chData );
		}
	}

	//***************************************
	// loadRawDataFloat 
	// binary data and element should be float
	//***************************************
	void loadRawDataFloat( char* fnPath )
	{
		FILE *fp = fopen(fnPath, "rb");
		if (!fp)
		{
			fprintf(stderr, "Error opening file '%s'\n", fnPath);
		}else{
			m_data = (float*)malloc(m_size*sizeof(float) );
			fread((void*)m_data, 1, m_size*sizeof(float), fp);
			fclose(fp);

			float min = 100000;
			float max = -100000;
			
			for( int i=0; i<m_size; i++ ){
				if( min > m_data[i] )min = m_data[i];
				if( max < m_data[i] )max = m_data[i];
			}
			printf("Raw data: max(%f)    min(%f)\n" , max, min);

			for( int i=0; i<m_size; i++ ){
				m_data[i] = ( m_data[i] - min )/ (max - min);
			}

		}
	}

	//***************************************
	// loadVectorData
	// 
	//***************************************
	void loadVectorData( char* fnPath )
	{	
		//load data
		FILE * fIn; 
		int totalNum;
		float *pData; 
		int dimension[3];
	
		fIn = fopen(fnPath, "rb");
		fread(dimension, sizeof(int), 3, fIn);
		pData = (float*)malloc( m_size *3 * sizeof(float) );
		fread(pData, sizeof(float), m_size*3, fIn);
		fclose(fIn);

		m_data = (float*)malloc( sizeof(float) * m_size );

		//get the max and min value in this data for normalize later   
		float min=100000;
		float max=-100000;
		for( int i=0; i<m_size; i++ ){
			int idx = i*3;
			float value = pData[idx]*pData[idx] + pData[idx+1]*pData[idx+1] + pData[idx+2]*pData[idx+2];
			value = sqrt(value);
			if( value< min )min=value;
			if( value>max )max = value;
			m_data[i] = value;
		}
		//print max and min voxel in this data
		printf("This Data: max(%f)    min(%f)\n" , max, min);
	
		//normalize to 0-1
		for( int i=0; i<m_size; i++ ){	
			m_data[i] = ( m_data[i] - min )/ (max - min);
		}

		free( pData );
	}

	int idx3dToIdx1d( int idx1, int idx2, int idx3 )
	{
		//idx = idx1*m_dim2*m_dim3 + idx2*m_dim3 + idx3;
		int idx = 0;
		if( m_dataType == VolumeDataType::RAW_DATA_UNSIGNEDCHAR ){
			idx = idx3*m_dim2*m_dim1 + idx2*m_dim1 + idx1;
		}else if( m_dataType == VolumeDataType::RAW_DATA_FLOAT ){
			idx = idx3*m_dim2*m_dim1 + idx2*m_dim1 + idx1;
		}else if( m_dataType == VolumeDataType::VECTOR_MAGNITUDE ){
			idx = idx3*m_dim2*m_dim1 + idx2*m_dim1 + idx1;
		}

		return idx;
	}
    
    float trilinear( float x, float y, float z,
                    float v000, float v001, float v010, float v011,
                    float v100, float v101, float v110, float v111)
	{
		float x0 = 0, x1 = 1;
		float y0 = 0, y1 = 1;
		float z0 = 0, z1 = 1;
		float xd = (x-x0)/(x1-x0);
		float yd = (y-y0)/(y1-y0);
		float zd = (z-z0)/(z1-z0);
        
		float c00 = v000* (1-xd) + v100*xd;
		float c10 = v010* (1-xd) + v110*xd;
		float c01 = v001* (1-xd) + v101*xd;
		float c11 = v011* (1-xd) + v111*xd;
        
		float c0 = c00*(1-yd) + c10*yd;
		float c1 = c01*(1-yd) + c11*yd;
        
		float c = c0*(1-zd) + c1*zd;
        
		return c;
	}

	int m_dim1, m_dim2, m_dim3;
	float* m_data;
	int m_size;
	int m_dataType;
	float m_hRadius;
	int m_imgResU, m_imgResV, m_imgResD, m_imgResB;
	float m_gau1DTable[10];
	float m_gau2DTable[10][10];
	//float m_histo[512][126][64][64];
};


#endif
