#include <iostream>
#include<cstdio>
#include<fstream>
#include<sstream>
#include <cuda.h>
#include<bits/stdc++.h>
#include<vector>
#include<cmath>
using namespace std;
const int N=29;
#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
// Macro for timing kernel runs
#define START_METER {\
    cudaEvent_t start, stop;\
    float elapsedTime;\
    cudaEventCreate(&start);\
    cudaEventRecord(start, 0);
#define STOP_METER cudaEventCreate(&stop);\
    cudaEventRecord(stop, 0);\
    cudaEventSynchronize(stop);\
    cudaEventElapsedTime(&elapsedTime, start, stop);\
    printf("Elapsed time : %f ms\n", elapsedTime);\
                }
bool sortinrev(const pair < double,vector<double> > &a,  
               const pair < double,vector<double> > &b) 
{ 
       return (a.first > b.first); 
}
__global__ 
void mean_col(double *data,double *mean )
{
    int index = threadIdx.x + blockIdx.x*blockDim.x;
    if(index<29)
    {
    for(int j = 0;j<24998;j++)
        {
        mean[index] = mean[index] + data[29*j+index];
        }
    mean[index]=double(mean[index]/24998.0);
    }
}

__global__ 
void std_col(double *data,double *mean,double *std){
        int index = threadIdx.x;
        double t = 0;
        for(int j = 0;j<24998;j++)
        {
            double t2 = data[29*j+index]- mean[index];
            t = t + t2*t2;
        }
        t= t/24998;
        std[index] = sqrt(t);
}
__global__
void data_normalize(double *data,double *mean,double *std){
    int index = threadIdx.x;
    for(int i = 0;i<24998;i++)
    {
            data[i*29+index] = (data[i*29+index]- mean[index])/std[index];
    }
}
__global__
void covariance(double *data,double *r){
    int index = threadIdx.x;
    for(int j=0;j<29;j++)
    {
        for(int k=0;k<24998;k++)
        {
            r[index*29+j] +=  (data[k*29+index])*data[k*29+j];
        }
    }
}
__global__ void gpu_matrix_transpose(double* mat_in, double* mat_out) 
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

    if (idx < 29 && idy < 24998) 
    {
        unsigned int pos = idy * 29 + idx;
        unsigned int trans_pos = idx * 24998 + idy;
        mat_out[trans_pos] = mat_in[pos];
    }
}



int main(void)
{
    time_t start, end; 
    start = clock();
    ios_base::sync_with_stdio(false);

    int col = 29;
    int rw = 24998;
    int row = 29;

    double *r1,*data,*d_data,*d_mean,*transpose,*d_t;

    r1 = (double*)calloc(row*N,sizeof(double));

    data = (double*)malloc(rw*col*sizeof(double));
    transpose = (double*)malloc(rw*col*sizeof(double));

    cudaMalloc((void**)&d_data, rw*col*sizeof(double));
    cudaMalloc((void**)&d_t, rw*col*sizeof(double));

    //********** file read ***********
    std::ifstream file("intrusion.csv");
    
    std::string lin;
    std::getline(file, lin);
    
    for(int rows = 0; rows < rw; ++rows)
    {

        string line;
        getline(file, line);
        if ( !file.good() )
            break;
        std::stringstream iss(line);
        for (int cols = 0; cols < col; ++cols)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> data[rows*col + cols];
        }
        ////cout<<endl;
    }

    //******** file read end *********




    double *mean;
    mean = (double*)calloc(col,sizeof(double));
    cudaMalloc((void**)&d_mean, col*sizeof(double)); 
    cudaMemcpy(d_data, data, rw*col*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mean, mean, col*sizeof(double), cudaMemcpyHostToDevice);
    mean_col <<<1, 29>>>(d_data,d_mean);
    //mean_col <<<1, 29>>>(d_mean);
    cudaDeviceSynchronize();
    cudaMemcpy(mean, d_mean, col*sizeof(double), cudaMemcpyDeviceToHost);
    
    double *std,*d_std;
    std = (double*)malloc(col*sizeof(double));
    cudaMalloc((void**)&d_std, col*sizeof(double)); 
    //cudaMemcpy(data, d_data, rw*col*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(d_std, std, col*sizeof(double), cudaMemcpyHostToDevice);
    std_col <<< 1,29 >>>(d_data,d_mean,d_std);
    cudaDeviceSynchronize();
    cudaMemcpy(data, d_data, rw*col*sizeof(double), cudaMemcpyDeviceToHost);
    data_normalize <<< 1,29 >>>(d_data,d_mean,d_std);
    cudaDeviceSynchronize();
    
    cudaMemcpy(data, d_data, rw*col*sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_std);
    cudaFree(d_mean);
    free(mean);
    free(std);
    double *d_r;
    
    cudaMalloc((void**)&d_r, N*col*sizeof(double));
    cudaMemcpy(d_r, r1, N*col*sizeof(double), cudaMemcpyHostToDevice);
    
    covariance<<< 1, 29 >>>(d_data,d_r);
    cudaDeviceSynchronize();
    cudaMemcpy(r1, d_r, N*col*sizeof(double), cudaMemcpyDeviceToHost);
    ////cout<<"w"<<endl;

    
    
    double w[row];
    double v[row][N] ={0};
    double r[row][N]={0};
    for(int i=0;i<row;i++){
        for(int j=0;j<N;j++)
            r[i][j] = r1[i*N+j];
    }

    //cusolver *****************************
    cusolverDnHandle_t handle;
    cusolverDnCreate(&handle);
    int lwork;
    int rows = 29,cols=29;
    cusolverDnDgesvd_bufferSize(
        handle,
        rows,
        cols,
        &lwork);

    double *d_A;
    cudaMalloc(&d_A, sizeof(double)*N*cols);
    cudaMemcpy(d_A, r1, sizeof(double)*N*cols, cudaMemcpyHostToDevice);

    double *d_S;
    cudaMalloc(&d_S, sizeof(double)*rows);

    double *d_U;
    cudaMalloc(&d_U, sizeof(double)*rows*rows);

    double *d_VT,*VT,*S;

    VT = (double*)calloc(N*N,sizeof(double));
    S = (double*)calloc(N,sizeof(double));

    cudaMalloc(&d_VT, sizeof(double)*rows*rows);

    double *d_work;
    cudaMalloc(&d_work, sizeof(double)*lwork);

    double *d_rwork;
    cudaMalloc(&d_rwork, sizeof(double)*(rows - 1));

    int *devInfo;
    cudaMalloc(&devInfo, sizeof(int));

    for (int t = 0; t < 10; t++)
    {
        signed char jobu = 'A';
        signed char jobvt = 'A';
            cusolverDnDgesvd(
            handle,
            jobu,
            jobvt,
            rows,
            cols,
            d_A,
            rows,
            d_S,
            d_U,
            rows,
            d_VT,
            rows,
            d_work,
            lwork,
            d_rwork,
            devInfo);
    }

    cudaMemcpy(VT, d_VT, N*N*sizeof(double), cudaMemcpyDeviceToHost);
    
    cudaMemcpy(S, d_S, N*sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_rwork);
    cudaFree(d_S);
    cudaFree(d_U);
    cudaFree(d_VT);
    cudaFree(d_work);

   /* for(int i=0;i<N;i++){
        std::cout<<S[i]<<"\n";
    }
*/
    //

    //******
    //sort(w.begin(),w.end());
    double eigen_sum  = 0;
    for(int i =0;i<row;i++)
    {    
        eigen_sum  += fabs(S[i]);
        //cout<<w[i]<<" ";   
    }
    //cout<<endl;

    double thres= 0.9;

    vector< vector<double> > v1(row, vector<double> (col,0));

    for(int i =0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            v1[i][j] = VT[i+29*j];
        }
    }

    vector< pair < double,vector<double> > > vc;

    for(int i =0;i<row;i++)
    {
        vc.push_back( make_pair(S[i],v1[i]) );

    }
    //sort(vc.begin(), vc.end(),sortinrev);

    double sum = 0;
    int dim = 0;
    for(int i =0;i<row;i++)
    {
        sum += fabs(vc[i].first);
        dim++;
        double temp = sum/eigen_sum;
        if(temp>0.9)
        {
            break;
        }

    }


    vector< vector<double> > red_data(rw, vector<double> (dim,0));


    for(int i =0;i<rw;i++)
    {
        for(int j=0;j<dim;j++)
        {
            for(int k=0;k<col;k++)
            {
                red_data[i][j] = red_data[i][j] + data[i*29+k]*(vc[j].second[k]);
            }
        }
    }

    end = clock();  
  
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "Time taken by program is : " << fixed  
         << time_taken << setprecision(5); 
    cout << " sec " << endl; 

    ofstream outfile;
    outfile.open ("reduced_data.csv");
    for (int i = 0;i<dim-1;i++)
    {
        outfile<<"col"<<i<<",";
    }
    outfile<<"col"<<dim-1<<endl;

    for(int i =0;i<rw;i++)
    {
        for(int j=0;j<dim-1;j++)
        {
            outfile<<red_data[i][j]<<",";
        }
        outfile<<red_data[i][dim-1]<<endl;
    }

    outfile.close();
    return 0;
}



