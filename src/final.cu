#include <iostream>
#include <stdio.h>      
#include <math.h>   
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 
#include <cmath>
#include <vector>
#include <ctime>
#include <omp.h>
#include <fstream>

using namespace std;

#define Dirichlet
#define Gaussian_Pulse
#define zero_time_derivative
#define source_term0

//Denotes the no of points the axis is divided into
const int ARRAY_SIZE = 50; 

// Mesh_info >> x0,x1,y0,y1,t0,t1,Arraysize, X_resolution, Y_resolution , Time_resolution
__global__ void solver_kernel(const float* mesh_info,const float* X_MESH, const float* Y_MESH, const float* TIME_MESH, float* d_out ){

  int Px = blockIdx.x * blockDim.x + threadIdx.x;
  int Py = blockIdx.y * blockDim.y + threadIdx.y;
  int px = threadIdx.x;
  int py = threadIdx.y;

  //shared memory for mesh_info and then it is copied;
  __shared__ float mesh_sharedinfo[10];
  __shared__ float r1;
  __shared__ float r2;
  __shared__ float r3;
  __shared__ float r4;
  r1=r2=r3=r4=0;
  for(int i=0;i<10;i++){
    mesh_sharedinfo[i]=mesh_info[i];
  } 

  

  float Cx = mesh_sharedinfo[9]/mesh_sharedinfo[7], Cy = mesh_sharedinfo[9]/mesh_sharedinfo[7];

  //extern __shared__ float u1[(blockDim.x +2)*(blockDim.y +2)];
  extern __shared__ float u1[];
 

  ///////INITIAL CONDITION SETTING////////
  // U at time = T0;
  #ifdef Gaussian_Pulse
    //u=a*e^(-(x-b0)^2/(2*c0^2)-(y-b1)^2/(2*c1^2))
  float a= 2, b0 = 5, b1 =5, c0=1 , c1 =1 ;
  u1[(blockDim.x+2)*(px+1)+(py+1)] = a*exp((-1*(X_MESH[Px]-b0)*(X_MESH[Px]-b0)/(2*c0*c0)) +  (-1*(Y_MESH[Py]-b1)*(Y_MESH[Py]-b1)/(2*c1*c1)) );
  __syncthreads();
  #endif

  #ifdef sine_pulse
  float a= 1, b0 = 5, b1 =5, c0=1 , c1 =1 ;
  u1[(blockDim.x+2)*(px+1)+(py+1)] = a*sin(c0*(X_MESH[Px]-b0)+ c1*(Y_MESH[Py]-b1));
  __syncthreads();
  #endif
  extern __shared__ float u2[];
  extern __shared__ float u[];
  //////BOUNDARY CONDITION ///////
   #ifdef Dirichlet

   if(blockIdx.x == 0){
      u1[(blockDim.x+2)*(1)+(py+1)]=0;
   }
   if(blockIdx.x == gridDim.x-1){ 
  u1[(blockDim.x+2)*(blockDim.x)+(py+1)]=0;
   }

   if(blockIdx.y == gridDim.y-1){
  u1[(blockDim.x+2)*(px+1)+(blockDim.y)]=0;
   }  
   if(blockIdx.y == 0){
    u1[(blockDim.x+2)*(px+1)+(1)]=0;
   }  



  #endif

  #ifdef Von_Neumann

  #endif
  u2[(blockDim.x+2)*(px+1)+(py+1)]=u1[(blockDim.x+2)*(px+1)+(py+1)];
  d_out[Px*blockDim.x*gridDim.x + Py] =u1[(blockDim.x+2)*(px+1)+(py+1)];

for(int i=0;i<ARRAY_SIZE;i++){


  //int i=0;
   if(r1==0){
      u2[(blockDim.x+2)*(0)+(py+1)]= d_out[(ARRAY_SIZE*ARRAY_SIZE-1)*i + Px*blockDim.x*gridDim.x + Py];//d_out[(ARRAY_SIZE*ARRAY_SIZE-1)*i + ((blockIdx.x-1) * blockDim.x + (blockDim.x-1))*(blockDim.x)*gridDim.x + ((blockIdx.y) * blockDim.y + threadIdx.y) ];
      
    }

    if(r2==0){
      u2[(blockDim.x+2)*(px+1)+(0)]= d_out[Px*blockDim.x*gridDim.x + Py];//d_out[(ARRAY_SIZE*ARRAY_SIZE-1)*i+(blockIdx.x * blockDim.x + threadIdx.x)*(blockDim.x)*gridDim.x + ((blockIdx.y-1) * blockDim.y + (blockDim.y-1)) ];
      
    }

    if(r3==0){
      u2[(blockDim.x+2)*(blockDim.x+1)+(py+1)]= d_out[Px*blockDim.x*gridDim.x + Py];//d_out[(ARRAY_SIZE*ARRAY_SIZE-1)*i+((blockIdx.x+1) * (blockDim.x) + 0)*(blockDim.x)*gridDim.x + ((blockIdx.y) * blockDim.y + threadIdx.y) ];
      
    }

    if(r4==0){
      u2[(blockDim.x+2)*(px+1)+(blockDim.y+1)]= d_out[Px*blockDim.x*gridDim.x + Py];//d_out[(ARRAY_SIZE*ARRAY_SIZE-1)*i+(blockIdx.x * blockDim.x + threadIdx.x)*(blockDim.x)*gridDim.x + ((blockIdx.y+1) * blockDim.y + 0) ];
      
    }
     __syncthreads();
    //u[(blockDim.x+2)*(px+1)+(py+1)] = -u1[(blockDim.x+2)*(px+1)+(py+1)]+ 2*u2[(blockDim.x+2)*(px+1)+(py+1)];


    u[(blockDim.x+2)*(px+1)+(py+1)] = -u1[(blockDim.x+2)*(px+1)+(py+1)]+ 2*u2[(blockDim.x+2)*(px+1)+(py+1)] + Cx*Cx*(u1[(blockDim.x+2)*(px+2)+(py+1)]-2*u1[(blockDim.x+2)*(px+1)+(py+1)] + u1[(blockDim.x+2)*(px)+(py+1)]) + Cy*Cy*(u1[(blockDim.x+2)*(px+1)+(py+2)]-2*u1[(blockDim.x+2)*(px+1)+(py+1)] + u1[(blockDim.x+2)*(px+1)+(py)]);
    if(i==0){
     u[(blockDim.x+2)*(px+1)+(py+1)] = u[(blockDim.x+2)*(px+1)+(py+1)]/2;  
    }

  //////BOUNDARY CONDITION ///////
   #ifdef Dirichlet

   if(blockIdx.x == 0){
      u[(blockDim.x+2)*(1)+(py+1)]=0;
   }
   if(blockIdx.x == gridDim.x-1){ 
  u[(blockDim.x+2)*(blockDim.x)+(py+1)]=0;
   }

   if(blockIdx.y == gridDim.y-1){
  u[(blockDim.x+2)*(px+1)+(blockDim.y)]=0;
   }  
   if(blockIdx.y == 0){
    u[(blockDim.x+2)*(px+1)+(1)]=0;
   }  

  #endif

    d_out[(ARRAY_SIZE*ARRAY_SIZE -1)*(i+1) + Px*blockDim.x*gridDim.x + Py] = u[(blockDim.x+2)*(px+1)+(py+1)];
    __syncthreads();
    u1[(blockDim.x+2)*(px+1)+(py+1)] = u2[(blockDim.x+2)*(px+1)+(py+1)];
    __syncthreads();
    u2[(blockDim.x+2)*(px+1)+(py+1)] = u[(blockDim.x+2)*(px+1)+(py+1)];
    //__syncthreads();
   // u[(blockDim.x+2)*(px+1)+(py+1)]=0;


}
}


int main(int argc, char ** argv) {

ofstream myfile;
myfile.open ("solution3.txt");

// SURFACE DIMENSIONS

int X1 = 0, X2=10;    //X1<X<X2
int Y1=0,Y2=10;   //0<Y<10
int TIME = 10;   //0<TIME<1


///// Generate MESH for the SURFACE

float* x_mesh = new float[ARRAY_SIZE+1];
float* y_mesh= new float[ARRAY_SIZE+1];
float* time_mesh= new float[ARRAY_SIZE+1];
float X_res = (X2-X1)*1.0f/ARRAY_SIZE;
float Y_res = (Y2-Y1)*1.0f/ARRAY_SIZE;
float TIME_res = TIME*1.0f/ARRAY_SIZE;
float Cx = TIME_res/X_res, Cy = TIME_res/Y_res;

for(int i=0; i<=ARRAY_SIZE ;i++){
  x_mesh[i] = X1 + i*X_res;
  y_mesh[i] = Y1 + i*Y_res;
  time_mesh[i] = i*TIME_res;
  cout<<x_mesh[i]<<endl;
}

float mesh_info_cpu[10];
mesh_info_cpu[0]=X1,mesh_info_cpu[1]=X2,mesh_info_cpu[2]=Y1,mesh_info_cpu[3]=Y2;
mesh_info_cpu[4]=0,mesh_info_cpu[5]=TIME,mesh_info_cpu[6]=ARRAY_SIZE,mesh_info_cpu[7]=X_res;
mesh_info_cpu[8]=Y_res,mesh_info_cpu[9]=TIME_res;

float h_out[ARRAY_SIZE*ARRAY_SIZE*ARRAY_SIZE];
for(int i=0;i<2500*50;i++){
  h_out[i]=100;
}

float * d_out;
float * X_MESH;
float * Y_MESH;
float * TIME_MESH;
float * mesh_info;

cudaMalloc((void**) &mesh_info, 10*sizeof(float));
cudaMalloc((void**) &X_MESH, ARRAY_SIZE*sizeof(float));
cudaMalloc((void**) &Y_MESH, ARRAY_SIZE*sizeof(float));
cudaMalloc((void**) &TIME_MESH, ARRAY_SIZE*sizeof(float));
cudaMalloc((void**) &d_out,ARRAY_SIZE* ARRAY_SIZE*ARRAY_SIZE*sizeof(float));

cudaMemcpy(mesh_info, mesh_info_cpu, 10*sizeof(float), cudaMemcpyHostToDevice);
cudaMemcpy(X_MESH, x_mesh, ARRAY_SIZE*sizeof(float), cudaMemcpyHostToDevice);
cudaMemcpy(Y_MESH, y_mesh, ARRAY_SIZE*sizeof(float), cudaMemcpyHostToDevice);
cudaMemcpy(TIME_MESH, time_mesh, ARRAY_SIZE*sizeof(float), cudaMemcpyHostToDevice);

solver_kernel<<<dim3(2,2,1),dim3(25,25,1), 4*27*27*sizeof(float)>>>(mesh_info,X_MESH,Y_MESH,TIME_MESH,d_out);

cudaMemcpy(h_out, d_out,ARRAY_SIZE* ARRAY_SIZE*ARRAY_SIZE*sizeof(float) , cudaMemcpyDeviceToHost);

cudaFree(X_MESH);
cudaFree(d_out);
cudaFree(TIME_MESH);
cudaFree(Y_MESH);
cudaFree(mesh_info);

for(int t=0;t<ARRAY_SIZE;t++){
for(int i=0;i<ARRAY_SIZE;i++){
  for(int j=0;j<ARRAY_SIZE;j++){

  myfile<<h_out[(ARRAY_SIZE*ARRAY_SIZE-1)*t+ i*ARRAY_SIZE +j]<<" ";
}
myfile<<"\n";
}
myfile<<"TIMHdnvjsbsjbgs"<<"\n";

}
}
