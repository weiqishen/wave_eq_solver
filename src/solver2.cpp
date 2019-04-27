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

const int ARRAY_SIZE = 10;



///////// Boundary Conditions Function  //////////////


void boundary_setter(float** matrix){
	#ifdef Dirichlet
	// u = 0 at boundary 
	for(int i=0;i<ARRAY_SIZE+1;i++){
		matrix[0][i]=0;
		matrix[i][0]=0;
		matrix[ARRAY_SIZE][i]=0;
		matrix[i][ARRAY_SIZE]=0;
	}
	#endif

	#ifdef Von_Neumann
	// u = 0 at boundary 
	for( int i=0; i<ARRAY_SIZE+1;i++  ){

	}


	#endif


}


/// for ut(x, y, 0) = V (x, y), (x, y) ∈ Ω, (115)
float V(int x , int j ){
	#ifdef zero_time_derivative
	return 0;
	#endif
}


int main(int argc, char *argv[])
{


ofstream myfile;
myfile.open ("solution2.txt");

// SURFACE DIMENSIONS

int X = 10;    //0<X<10
int Y=10;		//0<Y<10
int TIME = 10;   //0<TIME<1


///// Generate MESH for the SURFACE

float* x_mesh = new float[ARRAY_SIZE+1];
float* y_mesh= new float[ARRAY_SIZE+1];
float* time_mesh= new float[ARRAY_SIZE+1];
float X_res = X*1.0f/ARRAY_SIZE;
float Y_res = Y*1.0f/ARRAY_SIZE;
float TIME_res = TIME*1.0f/ARRAY_SIZE;
float Cx = TIME_res/X_res, Cy = TIME_res/Y_res;

for(int i=0; i<=ARRAY_SIZE ;i++){
	x_mesh[i] = i*X_res;
	y_mesh[i] = i*Y_res;
	time_mesh[i] = i*TIME_res;
	cout<<x_mesh[i]<<endl;
}



//////////     PROBLEM TO SOLVE : DtDtu = DxDxu + DyDyu + f;
//////////      Source term is taken to be zero to be zero for now 
/////////		Non dimensionlising 
///				u* = u/u 0 , t∗ = t/(L/c) , x∗ = x/L , y∗ = y/L and f ∗ = f /(u0 c2 /L2 ) .
////////// 					non-dimensionalization
//								u_ref 1
//								c_ref 1
//								L_ref 1
//								t_ref 1



 float** u = new float*[ARRAY_SIZE+1];
for(int i = 0; i < ARRAY_SIZE+1; ++i){
  u[i] = new float[ARRAY_SIZE+1]; /// solution of wave equation at Time t;
}

 float** u_t_1 = new float*[ARRAY_SIZE+1];
for(int i = 0; i < ARRAY_SIZE+1; ++i){
  u_t_1[i] = new float[ARRAY_SIZE+1]; /// solution of wave equation at Time t;
} /// Sloution at time t-1 ;
 float** u_t_2 = new float*[ARRAY_SIZE+1];
for(int i = 0; i < ARRAY_SIZE+1; ++i){
  u_t_2[i] = new float[ARRAY_SIZE+1]; /// solution of wave equation at Time t;
} /// Solution at time t-2 ;


//float Wave_solution[ARRAY_SIZE+1][ARRAY_SIZE+1][ARRAY_SIZE+1];


/////////   Initial conditions //////////


#ifdef Gaussian_Pulse
	//u=a*e^(-(x-b0)^2/(2*c0^2)-(y-b1)^2/(2*c1^2))
		float a= 2, b0 = 5, b1 =5, c0=1 , c1 =1 ;
	for(int j=0;j<ARRAY_SIZE+1;j++){
		for(int i=0;i<ARRAY_SIZE+1;i++){
			u_t_2[j][i] = a*exp((-1*(x_mesh[j]-b0)*(x_mesh[j]-b0)/(2*c0*c0)) +  (-1*(y_mesh[i]-b1)*(y_mesh[i]-b1)/(2*c1*c1)) );
		}
	}


#endif

#ifdef sine_pulse
	//sin wave: u_t=a*sin(b0*(x-c0))sin(b1*(y-c1))
	float a= 1, b0 = 5, b1 =5, c0=1 , c1 =1 ;
	for(int j=0;j<ARRAY_SIZE+1;j++){
		for(int i=0;i<ARRAY_SIZE+1;i++){
			u_t_2[j][i] = a*sin(c0*(x_mesh[j]-b0)+ c1*(y_mesh[i]-b1));
		}	
	}


#endif




///// Computation of U_t_1 for t= 2*Time_res;

for(int i =1;i<ARRAY_SIZE;i++){
	for(int j =1; j<ARRAY_SIZE;j++){
		u_t_1[i][j]= u_t_2[i][j] + TIME_res*V(i,j) + 0.5*Cx*Cx*(u_t_2[i+1][j]-2*u_t_2[i][j] + u_t_2[i-1][j] ) + 0.5*Cy*Cy*(u_t_2[i][j+1]-2*u_t_2[i][j] + u_t_2[i][j-1]  );
	}
}

/// Setting Boundary Condition for Time t;
boundary_setter(u_t_1);
boundary_setter(u_t_2);

/// Storing Solution in Wave_solution variable 

myfile<<"Time 0\n";
for(int i =0 ; i<ARRAY_SIZE+1;i++){
	for(int j =0;j<ARRAY_SIZE+1;j++){
		myfile<<u_t_2[i][j]<<" ";
	}
	myfile<<"\n";
	
}

myfile<<"Time 1\n";
for(int i =0 ; i<ARRAY_SIZE+1;i++){
	for(int j =0;j<ARRAY_SIZE+1;j++){
		myfile<<u_t_1[i][j]<<" ";
	}
	myfile<<"\n";
	
}


for(int n=1;n<ARRAY_SIZE+1;n++){
	//myfile<<"Time "<<n+1<<"\n";
	for(int i =1;i<ARRAY_SIZE;i++){
		for(int j =1; j<ARRAY_SIZE;j++){
			u[i][j]= -u_t_2[i][j] +  2*u_t_1[i][j]  + 0.5*Cx*Cx*(u_t_1[i+1][j]-2*u_t_1[i][j] + u_t_1[i-1][j] ) + 0.5*Cy*Cy*(u_t_1[i][j+1]-2*u_t_1[i][j] + u_t_1[i][j-1]  );
			
		}
	}
	boundary_setter(u);
	myfile<<"Time "<<n+1<<"\n";
	for(int i =0 ; i<ARRAY_SIZE+1;i++){
		for(int j =0;j<ARRAY_SIZE+1;j++){
			myfile<<u[i][j]<<" ";
		}
	myfile<<"\n";

}

	
	for(int q=0;q<ARRAY_SIZE+1;q++){
		for(int w=0;w<ARRAY_SIZE+1;w++){
			float m= u_t_2[q][w];
			float e = u_t_1[q][q];
			float r = u[q][w];
			u_t_2[q][w]=e;
			u_t_1[q][q]=r;
		}
	}



}

for(int i = 0; i < ARRAY_SIZE+1; ++i){
delete [] u[i];}
delete [] u;

for(int i = 0; i < ARRAY_SIZE+1; ++i){
delete [] u_t_1[i];}
delete [] u_t_1;

for(int i = 0; i < ARRAY_SIZE+1; ++i){
delete [] u_t_2[i];}
delete [] u_t_2;
delete [] x_mesh;
delete [] y_mesh;
delete [] time_mesh;
return 0;

}

