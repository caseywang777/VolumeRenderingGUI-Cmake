// CUDA runtime
#include <cuda_runtime.h>

#include"device_launch_parameters.h"

#include <stdio.h>

__global__ void cudaKernel(float* visHist, int* rawHistogramRay, float* devOtf, int K, int B, int D)

{
	//int i = blockDim.x * blockIdx.x + threadIdx.x;
	int id = blockIdx.x;

	float localVisHist[128];

    int base = id * D * B;
	int slcSample = 0;
	for( int d = 0; d < D; d++ ){
		for( int b = 0; b < B; b++ ){
			localVisHist[b] += rawHistogramRay[base + d * D + b];
		}
	}

}

 
extern "C" double runCudaKernel( float* visHist, int K, int D, int B, float* otf, int* rawHistogramRays )

{

		FILE* fp = fopen( "output.txt", "w" );

         float* devVisHist = 0;

         float* devOtf = 0;

         int* devRawHistogramRay = 0;

         cudaError_t cudaStatus;

 

         //Choose which GPU to run on, change this on a multi-GPU system.

         cudaStatus= cudaSetDevice(0);

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaSetDevice failed!  Do you havea CUDA-capable GPU installed?");

                   goto Error;

         }

 

         //Allocate GPU buffers for three vectors (two input, one output)    .

         cudaStatus= cudaMalloc((void**)&devVisHist, B * sizeof(float));

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaMalloc failed!");

                   goto Error;

         }

 

         cudaStatus= cudaMalloc((void**)&devOtf, B * sizeof(float));

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaMalloc failed!");

                   goto Error;

         }

 

         cudaStatus= cudaMalloc((void**)&devRawHistogramRay, K * B * D * sizeof(int));

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaMalloc failed!");

                   goto Error;

         }

 

         //Copy input vectors from host memory to GPU buffers.

         cudaStatus= cudaMemcpy( devRawHistogramRay, rawHistogramRays, K * B * D * sizeof(int), cudaMemcpyHostToDevice);

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaMemcpy failed!");

                   goto Error;

         }

 

         cudaStatus= cudaMemcpy(devOtf, otf, B * sizeof(float), cudaMemcpyHostToDevice);

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaMemcpy failed!");

                   goto Error;

         }

 

         //Launch a kernel on the GPU with one thread for each element.
		 fprintf( fp, "before kernel\n " );fflush( fp );
         cudaKernel<<<K,1>>>(devVisHist, devRawHistogramRay, devOtf, K, B, D);
		 fprintf( fp, "after kernel\n " );fflush( fp );
 

         //Check for any errors launching the kernel

         cudaStatus= cudaGetLastError();

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));

                   goto Error;

         }

 

         //cudaDeviceSynchronize waits for the kernel to finish, and returns

         //any errors encountered during the launch.

         cudaStatus= cudaDeviceSynchronize();

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaDeviceSynchronize returned error code %d after launchingaddKernel!\n", cudaStatus);

                   goto Error;

         }

 

         //Copy output vector from GPU buffer to host memory.

         cudaStatus= cudaMemcpy(visHist, devVisHist, B * sizeof(float), cudaMemcpyDeviceToHost);

         if(cudaStatus != cudaSuccess) {

                   fprintf(fp,"cudaMemcpy failed!");

                   goto Error;

         }

		 //output result
		 fprintf( fp, "before output\n " );fflush( fp );
		 for( int i=0; i<B; i++ ){
			 fprintf(fp, "%d : %f \n", i, visHist[i]);
		 }
		 fprintf( fp, "after output\n " );fflush( fp );
 

Error:

         cudaFree(devRawHistogramRay);

         cudaFree(devVisHist);

         cudaFree(devOtf);

 


         if(cudaStatus != cudaSuccess) {

                   fprintf(stderr,"addWithCuda failed!");

                   return 1;

         }

 

}