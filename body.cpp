#pragma once
#include <limits> 
#include "head.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
#include <windows.h>
#include <locale>
#include <sstream>
using namespace std;

//Matrix-operations....
int MatrixOp::InitiateParam(int Vic,int k1,int k2){
	if (Vic==0){
		MatrixOp();
		d=0;
		this->Count1 = 0;
		this->Count2 = 0;
		this->a1 = 0;
		this->a2 = 0;
		this->a3 = 0;
		this->a4 = 0;
	}
	else{
		this->Count1 = Vic;
		d=1;
		this->Count2 = 0;
		this->a1 = k1;
		this->a2 = k2;
	}
	return 1;
}
vector<vector<double>> MatrixOp::defineMatrix(int i,int j){
	vector<vector<double>> Mat;
	Mat.resize(i);
	for (int I = 0; I<i; I++)
		Mat[I].resize(j);
	return Mat;
}
vector<double> MatrixOp::defineVector(int n){
	vector<double> v;
	v.resize(n);
	return v;
}
int MatrixOp::PrintMat(vector<vector<double>>& A,int n,int m){
	cout << endl;
	for (int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				if (j==0){
					cout <<"| "<< A[i][j] << "  ";
				}
				if (j!=0 && j<m-1){
					cout << A[i][j] << "  ";
				}
				if(j==(m-1)){
					cout << A[i][j] << " | ";
				}
			}
			cout << endl;
		}
		cout << endl;
	return 1;
}
int MatrixOp::PrintVec(vector<double>& u,int n){
	cout << endl;
	for (int h=0;h<n;h++){
		cout << "|"<< u[h] <<"|"<< endl;
	}
	cout << endl;
	return 1;
}
int MatrixOp::PrintVecint(vector<int>& u,int n){
	cout << endl;
	for (int h=0;h<n;h++){
		cout << "|"<< u[h] <<"|"<< endl;
	}
	cout << endl;
	return 1;
}
double MatrixOp::tr(vector<vector<double>>& A,int n,int m,double tr){
	if (n==m){
		for(int i=0;i<n;i++){
		tr = tr+A[i][i];
		}
		cout << " tr= " <<tr;
		return tr;
	}
	cout << "you need to squere the matrix to calculate a Trace.";
	return 1;
}
bool MatrixOp::InsertationmethodaMat(MatrixOp& mo,vector<vector<double>>& A,bool Rnd,bool Normalized,int i, int j){
	Unit = false;
	double Arg=0;
	double tr=0;
	double Check =0;
	double I = -1;
	std::srand(std::time(0));
	double s = 1;
	if (Rnd==true){
		for(int a=0;a<i;a++){
			for(int b=0;b<j;b++){
				int rand = std::rand();
				A[a][b] = pow(I,s+1)*rand;
				s=s+1;
			}
		}
	}
	if (Unit){
		for(int a=0;a<i;a++){
			for(int b=0;b<j;b++){
					A[a][b] = b;}
			}
	}
		if(Normalized==true){
			for(int a=0;a<i;a++){
				for(int b=0;b<j;b++){
					Arg = Arg+pow(A[a][b],2);
				}
			}
			for(int a=0;a<i;a++){
				for(int b=0;b<j;b++){
				A[a][b] = A[a][b]/sqrt(Arg);
				Check = Check +pow(A[a][b],2);
				}
			}
		}

	if (!Rnd&&!Unit){
		for(int a=0;a<i;a++){
			for(int b=0;b<j;b++){
				cout << endl << "please enter the Value of A[" << a <<"][" << b << "]" << endl;
				cin >> A[a][b];
				cout << endl;
			}
		}
		if(Normalized==true){
			for(int a=0;a<i;a++){
				for(int b=0;b<j;b++){
					Arg = Arg+pow(A[a][b],2);
				}
			}
			for(int a=0;a<i;a++){
				for(int b=0;b<j;b++){
				A[a][b] = A[a][b]/sqrt(Arg);
				Check = Check +pow(A[a][b],2);
				}
			}
		}
	}
	if (Count1==0){
		G.resize(i);
	for (int I = 0; I<i; I++){
		G[I].resize(j);
	}
	for(int f=0;f<i; f++){
		for(int z=0;z<i; z++){
			G[f][z] = A[f][z];
			}
		}
	}
	if (Count1!=0){
		P.resize(i);
	for (int I = 0; I<i; I++){
		P[I].resize(j);
	}
	for(int f=0;f<i; f++){
		if ( Unit ){
			for(int z=0;z<i; z++){
				if ( f==z ){
					P[f][z] = 0;
				}
				else { P[f][z] = 0; }
			}
		}
		else{
		for(int z=0;z<i; z++){
			P[f][z] = A[f][z];
			}
		}}
	}
	if (Count1==0){
		this->a1=i;
		this->a2=j;
	}
	else{
		//cout << " G is: " << endl;
		//cout << " P is: " << endl;
	}
	this->Count1++;
	return false;
}
vector<double> MatrixOp::Insertationmethodavec(MatrixOp& mo, vector<double>& u,bool Rnd,bool Normalized,int i){
	double Arg=0;
	double tr=0;
	double Check =0;
	double I = -1;
	std::srand(std::time(0));
	double s = 1;
	if (Rnd==true){
		for(int a=0;a<i;a++){
			int rand = pow(I,s+1)*std::rand();
			u[a] = rand;
			s=s+1;
			
		}
		if(Normalized==true){
			for(int a=0;a<i;a++){
				Arg = Arg+pow(u[a],2);
			}
			for(int a=0;a<i;a++){
				u[a] = u[a]/sqrt(Arg);
				Check = Check +pow(u[a],2);
			}
		}
	}
	else{
		for(int a=0;a<i;a++){
				cout << endl << "please enter the Value of A[" << a <<"]" << endl;
				cin >> u[a];
				cout << endl;
		}
		if(Normalized==true){
			for(int a=0;a<i;a++){
				Arg = Arg+pow(u[a],2);
			}
			for(int a=0;a<i;a++){
				u[a] = u[a]/sqrt(Arg);
				Check = Check +pow(u[a],2);
			}
		}
	}
	return u;
}
vector<vector<double>> MatrixOp::Primitivemulti(MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int i, int j){
	vector<vector<double>> Multip;
	Multip.resize(i);
	for (int I = 0; I<i; I++)
		Multip[I].resize(j);
	for (int x=0;x<i;x++){
		for (int y=0;y<j;y++){
			for(int k =0;k<i;k++){
				Multip[x][y] = Multip[x][y]+A[x][k]*B[k][y];
				}
			}
		}
	return Multip;
}
vector<vector<double>> MatrixOp::AddMat(MatrixOp& mo, vector<vector<double>>& O,vector<vector<double>>& M,int i, int j){
	vector<vector<double>> AddPro;
	AddPro.resize(i);
	for (int I = 0; I<i; I++){
		AddPro[I].resize(j);
	}
	for (int x=0;x<i;x++){
		for(int y=0;y<j;y++){
			AddPro[x][y] = O[x][y] + M[x][y];
		}
	}
	return AddPro;
}
vector<vector<double>> MatrixOp::SwAlmulti(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int N){
	if (N==1){
		Div = Div/2;
		lim = lim - Div;
	}
	if(N == 1){
		if (Div-lim==0){
			for (int I = 0; I<=3; I++){
				if(I<=1 && I==0){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II][III];
						}
					}
				}
				if(I<=1 && I==1){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II][III+Div];
						}
					}
				}
				if (I>=2 && I==2){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II+Div][III];
						}
					}
				}
				if (I>=2 && I==3){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II+Div][III+Div];
						}
					}
				}
			}
		}
		if (Div-lim!=0){
			for (int I = 0; I<=3; I++){
				if(I<=1 && I==0){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II][III];
						}
					}
				}
				if(I<=1 && I==1){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim-1;III++){
								Sw[I][II][III] = Sw[I][II][III]+G[II][III+lim];
						}
					}
				}
				if (I>=2 && I==2){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<lim;III++){
								Sw[I][II][III] = Sw[I][II][III]+G[II+lim][III];
						}
					}
				}
				if (I>=2 && I==3){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){
								Sw[I][II][III] = Sw[I][II][III]+G[II+lim][III+lim];
						}
					}
				}
			}
		}
	}
	if(N == 2){
		if (Div-lim==0){
			for (int I = 4; I<=7; I++){
				if(I<=7 && I==4){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II][III];
						}
					}
				}
				if(I<=7 && I==5){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II][III+Div];
						}
					}
				}
				if (I>=4 && I==6){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II+Div][III];
						}
					}
				}
				if (I>=4 && I==7){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II+Div][III+Div];
						}
					}
				}
			}
		}
		if (Div-lim!=0){
			for (int I = 4; I<=7; I++){
				if(I<=7 && I==4){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II][III];
						}
					}
				}
				if(I<=7 && I==5){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim-1;III++){
								Sw[I][II][III] = Sw[I][II][III]+P[II][III+lim];
						}
					}
				}
				if (I>=4 && I==6){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<lim;III++){
								Sw[I][II][III] = Sw[I][II][III]+P[II+lim][III];
						}
					}
				}
				if (I>=4 && I==7){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){
								Sw[I][II][III] = Sw[I][II][III]+P[II+lim][III+lim];
						}
					}
				}
			}
		}
	}
	if (N==1){
		N++;
		mo.SwAlmulti(Sw,mo,A,B,Div,lim,N);
		return A;
	}
	if (N==2){
		int az=8;
		SeqSwAl(Sw,mo,A,B,Div,lim,az);
	}
	return A;
}
vector<vector<double>> MatrixOp::initiateSwAl(MatrixOp& mo, vector<vector<double>>& A, int i, int j){
	int n,m;
	n=i/2;
	m=i/2;
	if (!Recursive){ 
		int N;
		N=1;
		if (i-j==0){
			Sw.resize(19);
			for (int I = 0; I<=18; I++){
				Sw[I].resize(n);
				for (int II = 0; II<n; II++){
					Sw[I][II].resize(m);
				}
			}
		}
		if(n-m!=0){
			Sw.resize(25);
			for (int I = 0; I<=24; I++){
				Sw[I].resize(n);
				for (int II = 0; II<n; II++){
					Sw[I][II].resize(m);
				}
			}
		}
		mo.SwAlmulti(Sw,mo,A,A,i,j,N);
	}
	else{
		if(Start){
			cout << "Start";
			OriginI = i;
			OriginJ = j;
			Return = 1;
			for (i;50<i;i = i-50){
				Return++;
			}
			for (int q=1;q<Return;q++){
				j = j-50;
			}
		}
		int N;
		N=1;
		Sw = this->Sw;
		if (i-j==0){
			Sw.resize(30);
			for (int I = 0; I<=29; I++){
				Sw[I].resize(i);
				for (int II = 0; II<i; II++){
					Sw[I][II].resize(j);
				}
			}
		}
		if(i-j!=0){
			Sw.resize(25);
			for (int I = 0; I<=24; I++){
				Sw[I].resize(i);
				for (int II = 0; II<j; II++){
					Sw[I][II].resize(j);
				}
			}
		}
		mo.SwAlmultiRec(Sw,mo,A,A,i,i,N);
	}
	return A;
}
vector<vector<double>> MatrixOp::SeqSwAl(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int I){
	bool Rec=true;
	MatrixOp Mo;
	if (Div-lim==0){
		if (Rec){
			Sw[I] = mo.initiateSwAlbranch(mo,AddMat(mo,Sw[0],Sw[3],Div,Div),Mo.AddMat(mo,Sw[4],Sw[7],Div,Div),Div,Div);
			I++;
			d++;
			Sw[I] = mo.initiateSwAlbranch(mo,mo.AddMat(mo,Sw[2],Sw[3],Div,Div),Sw[4],Div,Div);
			I++;
			d++;
			Sw[I] = mo.initiateSwAlbranch(mo,Sw[0],mo.MinusMat(mo,Sw[5],Sw[7],Div,Div),Div,Div);
			I++;
			d++;
			Sw[I] = mo.initiateSwAlbranch(mo,Sw[3],mo.MinusMat(mo,Sw[6],Sw[4],Div,Div),Div,Div);
			I++;
			d++;
			Sw[I] = mo.initiateSwAlbranch(mo,mo.AddMat(mo,Sw[0],Sw[1],Div,Div),Sw[7],Div,Div);
			I++;
			d++;
			Sw[I] = mo.initiateSwAlbranch(mo,mo.MinusMat(mo,Sw[2],Sw[0],Div,Div),mo.AddMat(mo,Sw[4],Sw[5],Div,Div),Div,Div);
			I++;
			d++;
			Sw[I] = mo.initiateSwAlbranch(mo,mo.MinusMat(mo,Sw[1],Sw[3],Div,Div),mo.AddMat(mo,Sw[6],Sw[7],Div,Div),Div,Div);
			I++;
			d++;
			Sw[I] = mo.AddMat(mo,AddMat(mo,Sw[8],Sw[14],Div,Div),MinusMat(mo,Sw[11],Sw[12],Div,Div),Div,Div);
			I++;
			Sw[I] = mo.AddMat(mo,Sw[10],Sw[12],Div,Div);
			I++;
			Sw[I] = mo.AddMat(mo,Sw[9],Sw[11],Div,Div);
			I++;
			Sw[I] = mo.AddMat(mo,AddMat(mo,Sw[10],Sw[13],Div,Div),MinusMat(mo,Sw[8],Sw[9],Div,Div),Div,Div);
			I++;
			D = mo.EndSeqSw(Sw,mo,Div,Div);
		}
		else{
			Sw[I] = mo.Primitivemulti(mo,mo.AddMat(mo,Sw[0],Sw[3],Div,Div),mo.AddMat(mo,Sw[4],Sw[7],Div,Div),Div,Div);
			I++;
			Sw[I] = mo.Primitivemulti(mo,mo.AddMat(mo,Sw[2],Sw[3],Div,Div),Sw[4],Div,Div);
			I++;
			Sw[I] = mo.Primitivemulti(mo,Sw[0],mo.MinusMat(mo,Sw[5],Sw[7],Div,Div),Div,Div);
			I++;
			Sw[I] = mo.Primitivemulti(mo,Sw[3],mo.MinusMat(mo,Sw[6],Sw[4],Div,Div),Div,Div);
			I++;
			Sw[I] = mo.Primitivemulti(mo,mo.AddMat(mo,Sw[0],Sw[1],Div,Div),Sw[7],Div,Div);
			I++;
			Sw[I] = mo.Primitivemulti(mo,mo.MinusMat(mo,Sw[2],Sw[0],Div,Div),mo.AddMat(mo,Sw[4],Sw[5],Div,Div),Div,Div);
			I++;
			Sw[I] = mo.Primitivemulti(mo,mo.MinusMat(mo,Sw[1],Sw[3],Div,Div),mo.AddMat(mo,Sw[6],Sw[7],Div,Div),Div,Div);
			I++;
			Sw[I] = mo.AddMat(mo,AddMat(mo,Sw[8],Sw[14],Div,Div),MinusMat(mo,Sw[11],Sw[12],Div,Div),Div,Div);
			I++;
			Sw[I] = mo.AddMat(mo,Sw[10],Sw[12],Div,Div);
			I++;
			Sw[I] = mo.AddMat(mo,Sw[9],Sw[11],Div,Div);
			I++;
			Sw[I] = mo.AddMat(mo,AddMat(mo,Sw[10],Sw[13],Div,Div),MinusMat(mo,Sw[8],Sw[9],Div,Div),Div,Div);
			I++;
			D = mo.EndSeqSw(Sw,mo,Div,Div);
		}
	}
	if (Div-lim!=0){
		int I = 8;
		Sw[I] = mo.Insertationmethodasw(Sw,mo,mo.Primitivemulti(mo,mo.AddMat(mo,Sw[0],Sw[3],lim,lim),mo.AddMat(mo,Sw[4],Sw[7],lim,lim),lim,lim),I,lim,lim);
	}
	return A;
}
vector<vector<double>> MatrixOp::Insertationmethodasw(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A ,int R,int C,int I){
	vector<vector<double>> X;
	X.resize(R);
	for (int c=0;c<C;c++){
		X[c].resize(C);
	}
	for (int i=0;i<(R-1);i++){
		for (int j=0;j<(C-1);j++){
			X[i][j] = X[i][j]+A[i][j];
		}
	}
	I++;
	return X;
}
vector<vector<double>> MatrixOp::MinusMat(MatrixOp& mo, vector<vector<double>>& O,vector<vector<double>>& M,int i, int j){
	vector<vector<double>> AddPro;
	AddPro.resize(i);
	for (int I = 0; I<i; I++){
		AddPro[I].resize(j);
	}
	for (int x=0;x<i;x++){
		for(int y=0;y<j;y++){
			AddPro[x][y] = O[x][y] - M[x][y];
		}
	}
	return AddPro;
}
vector<vector<double>> MatrixOp::EndSeqSw(vector<vector<vector<double>>>& Sw, MatrixOp& mo, int i, int j){
	for (int f=0;f<=3;f++){
		if(f==0){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					G[u][c] = Sw[15][u][c];
				}
			}
		}
		if(f==1){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					G[u][c+j] = Sw[16][u][c];
				}
			}
		}
		if(f==2){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					G[u+i][c] = Sw[17][u][c];
				}
			}
		}
		if(f==3){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					G[u+i][c+j] = Sw[18][u][c];
					}
				}
			}
		}
	cout << "The Result is.." << endl;
	for ( int k=0;k<800000000;k++ ){if(k==60000000){cout<<".";}if(k==200000000){cout<<".";}if(k==400000000){cout<<".";}}
	mo.PrintMat(G,2*i,2*i);
	return G;
}
vector<vector<double>> MatrixOp::Recursion(MatrixOp& mo, vector<vector<double>>& A, int i, int j){
	Recursive=false;
	Start=false;
	//Unit = true;
	Stop=0;
	Stop2=0;
	MoveRaw=0;
	MoveCol=0;
	fractal=0;
	initiateSwAl(mo,A,i,j);
	return A;
}
vector<vector<double>> MatrixOp::SwAlmultiRec(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int N){
	cout << endl << Return << endl;
	if (N==1){
		if(Start){
			Div = Div/2;
			lim = lim - Div;
			cout << Div;
			cout << lim;}}
	if(N == 1){
		if (Div-lim==0){
			for (int I = 0; I<=3; I++){
				if(I<=1 && I==0){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II+Div*Stop][III+Div*Stop2];
						}
					}
				}
				if(I<=1 && I==1){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II+Div*Stop][III+Div+Div*Stop2];
						}
					}
				}
				if (I>=2 && I==2){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II+Div+Div*Stop][III+Div*Stop2];
						}
					}
				}
				if (I>=2 && I==3){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II+Div+Div*Stop][III+Div+Div*Stop2];
						}
					}
				}
			}
		}
		if (Div-lim!=0){
			for (int I = 0; I<=3; I++){
				if(I<=1 && I==0){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+G[II+Div*Stop][III+lim*Stop2];
						}
					}
				}
				if(I<=1 && I==1){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim-1;III++){
								Sw[I][II][III] = Sw[I][II][III]+G[II+Div*Stop][III+lim+lim*Stop2];
						}
					}
				}
				if (I>=2 && I==2){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<lim;III++){
								Sw[I][II][III] = Sw[I][II][III]+G[II+lim+lim*Stop][III+lim*Stop2];
						}
					}
				}
				if (I>=2 && I==3){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){
								Sw[I][II][III] = Sw[I][II][III]+G[II+lim+lim*Stop][III+lim+lim*Stop2];
						}
					}
				}
			}
		}
	}
	if(N == 2){
		if (Div-lim==0){
			for (int I = 4; I<=7; I++){
				if(I<=7 && I==4){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II+Div*Stop][III+Div*Stop2];
						}
					}
				}
				if(I<=7 && I==5){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II+Div*Stop][III+Div+Div*Stop2];
						}
					}
				}
				if (I>=4 && I==6){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II+Div+Div*Stop][III+Div*Stop2];
						}
					}
				}
				if (I>=4 && I==7){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II+Div+Div*Stop][III+Div+Div*Stop2];
						}
					}
				}
			}
		}
		if (Div-lim!=0){
			for (int I = 4; I<=7; I++){
				if(I<=7 && I==4){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim;III++){ 
							Sw[I][II][III] = Sw[I][II][III]+P[II+Div*Stop][III+Div*Stop2];
						}
					}
				}
				if(I<=7 && I==5){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim-1;III++){
								Sw[I][II][III] = Sw[I][II][III]+P[II+Div*Stop][III+lim+lim*Stop2];
						}
					}
				}
				if (I>=4 && I==6){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<lim;III++){
								Sw[I][II][III] = Sw[I][II][III]+P[II+lim+lim*Stop][III+lim*Stop2];
						}
					}
				}
				if (I>=4 && I==7){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){
								Sw[I][II][III] = Sw[I][II][III]+P[II+lim+lim*Stop][III+lim+lim*Stop];
						}
					}
				}
			}
		}
	}
	if (N==1){
		cout << endl;
		N++;
		mo.SwAlmulti(Sw,mo,A,B,Div,lim,N);
		return A;
	}
	if (N==2){
		cout << endl;
		int az=8;
		cout << "!";
		SeqSwAl(Sw,mo,A,B,Div,lim,az);
	}
	return A;
}
MatrixOp::MatrixOp(void){
}
MatrixOp::~MatrixOp(void){
}
vector<vector<double>> MatrixOp::NotRecursive(MatrixOp& mo, vector<vector<double>>& A, int i, int j){
	MatrixOp Mo;
	Recursive=false;
	fractal=0;
	initiateSwAl(mo,A,i,j);
	return A;
}
vector<vector<double>> MatrixOp::SwAlmultibranch(vector<vector<vector<double>>> Swbranch, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int N){
	if (N==1){
		Div = Div/2;
		lim = lim - Div;
	}
	if(N == 1){
		if (Div-lim==0){
			for (int I = 0; I<=3; I++){
				if(I<=1 && I==0){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+A[II][III];
						}
					}
				}
				if(I<=1 && I==1){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+A[II][III+Div];
						}
					}
				}
				if (I>=2 && I==2){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+A[II+Div][III];
						}
					}
				}
				if (I>=2 && I==3){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+A[II+Div][III+Div];
						}
					}
				}
			}
		}
		if (Div-lim!=0){
			for (int I = 0; I<=3; I++){
				if(I<=1 && I==0){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+A[II][III];
						}
					}
				}
				if(I<=1 && I==1){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim-1;III++){
								Swbranch[I][II][III] = Swbranch[I][II][III]+A[II][III+lim];
						}
					}
				}
				if (I>=2 && I==2){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<lim;III++){
								Swbranch[I][II][III] = Swbranch[I][II][III]+A[II+lim][III];
						}
					}
				}
				if (I>=2 && I==3){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){
								Swbranch[I][II][III] = Swbranch[I][II][III]+A[II+lim][III+lim];
						}
					}
				}
			}
		}
	}
	if(N == 2){
		if (Div-lim==0){
			for (int I = 4; I<=7; I++){
				if(I<=7 && I==4){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+B[II][III];
						}
					}
				}
				if(I<=7 && I==5){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+B[II][III+Div];
						}
					}
				}
				if (I>=4 && I==6){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+B[II+Div][III];
						}
					}
				}
				if (I>=4 && I==7){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+B[II+Div][III+Div];
						}
					}
				}
			}
		}
		if (Div-lim!=0){
			for (int I = 4; I<=7; I++){
				if(I<=7 && I==4){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim;III++){ 
							Swbranch[I][II][III] = Swbranch[I][II][III]+B[II][III];
						}
					}
				}
				if(I<=7 && I==5){
					for (int II = 0; II<lim; II++){
						for(int III =0;III<lim-1;III++){
								Swbranch[I][II][III] = Swbranch[I][II][III]+B[II][III+lim];
						}
					}
				}
				if (I>=4 && I==6){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<lim;III++){
								Swbranch[I][II][III] = Swbranch[I][II][III]+B[II+lim][III];
						}
					}
				}
				if (I>=4 && I==7){
					for (int II = 0; II<Div; II++){
						for(int III =0;III<Div;III++){
								Swbranch[I][II][III] = Swbranch[I][II][III]+B[II+lim][III+lim];
						}
					}
				}
			}
		}
	}
	if (N==1){
		N++;
		mo.SwAlmultibranch(Swbranch,mo,A,B,Div,lim,N);
		return Dbranch[d];
	}
	if (N==2){
		int az=8;
		SeqSwAlbranch(Swbranch,mo,A,B,Div,lim,az);
	}
	return Dbranch[d];
}
vector<vector<double>> MatrixOp::initiateSwAlbranch(MatrixOp& mo,vector<vector<double>>& A, vector<vector<double>>& B, int i, int j){
	int n,m;
	int N;
	N=1;
	n = i/2;
	m = j/2;
	if (!Recursive){ 
			if (i-j==0){
				Swbranch.resize(19);
				for (int I = 0; I<=18; I++){
					Swbranch[I].resize(n);
					for (int II = 0; II<n; II++){
						Swbranch[I][II].resize(m);
						}
					}
					Dbranch.resize(9);
					for (int I = 0; I<=8; I++){
						Dbranch[I].resize(n*2);
					for (int II = 0; II<n*2; II++){
						Dbranch[I][II].resize(m*2);
					}
				}
			}
		}
		if(n-m!=0){
			Swbranch.resize(25);
			for (int I = 0; I<=24; I++){
				Swbranch[I].resize(n);
				for (int II = 0; II<n; II++){
					Swbranch[I][II].resize(m);
				}
			}
		}
	mo.SwAlmultibranch(Swbranch,mo,A,B,i,j,N);
	return Dbranch[d];
}
vector<vector<double>> MatrixOp::SeqSwAlbranch(vector<vector<vector<double>>> Swbranch, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int I){
	if (fractal==0){
	if (Div-lim==0){
		Swbranch[I] = mo.Primitivemulti(mo,mo.AddMat(mo,Swbranch[0],Swbranch[3],Div,Div),mo.AddMat(mo,Swbranch[4],Swbranch[7],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.Primitivemulti(mo,mo.AddMat(mo,Swbranch[2],Swbranch[3],Div,Div),Swbranch[4],Div,Div);
		I++;
		Swbranch[I] = mo.Primitivemulti(mo,Swbranch[0],mo.MinusMat(mo,Swbranch[5],Swbranch[7],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.Primitivemulti(mo,Swbranch[3],mo.MinusMat(mo,Swbranch[6],Swbranch[4],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.Primitivemulti(mo,mo.AddMat(mo,Swbranch[0],Swbranch[1],Div,Div),Swbranch[7],Div,Div);
		I++;
		Swbranch[I] = mo.Primitivemulti(mo,mo.MinusMat(mo,Swbranch[2],Swbranch[0],Div,Div),mo.AddMat(mo,Swbranch[4],Swbranch[5],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.Primitivemulti(mo,mo.MinusMat(mo,Swbranch[1],Swbranch[3],Div,Div),mo.AddMat(mo,Swbranch[6],Swbranch[7],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,mo.AddMat(mo,Swbranch[8],Swbranch[14],Div,Div),mo.MinusMat(mo,Swbranch[11],Swbranch[12],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,Swbranch[10],Swbranch[12],Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,Swbranch[9],Swbranch[11],Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,mo.AddMat(mo,Swbranch[10],Swbranch[13],Div,Div),mo.MinusMat(mo,Swbranch[8],Swbranch[9],Div,Div),Div,Div);
		I++;
		Dbranch[d] = mo.EndSeqSwbranch(Swbranch,mo,Div,Div);
	}
	}
	else{
		fractal--;
		Swbranch[I] = mo.initiateSwAlbranch(mo,mo.AddMat(mo,Swbranch[0],Swbranch[3],Div,Div),mo.AddMat(mo,Swbranch[4],Swbranch[7],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.initiateSwAlbranch(mo,mo.AddMat(mo,Swbranch[2],Swbranch[3],Div,Div),Swbranch[4],Div,Div);
		I++;
		Swbranch[I] = mo.initiateSwAlbranch(mo,Swbranch[0],mo.MinusMat(mo,Swbranch[5],Swbranch[7],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.initiateSwAlbranch(mo,Swbranch[3],mo.MinusMat(mo,Swbranch[6],Swbranch[4],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.initiateSwAlbranch(mo,mo.AddMat(mo,Swbranch[0],Swbranch[1],Div,Div),Swbranch[7],Div,Div);
		I++;
		Swbranch[I] = mo.initiateSwAlbranch(mo,mo.MinusMat(mo,Swbranch[2],Swbranch[0],Div,Div),mo.AddMat(mo,Swbranch[4],Swbranch[5],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.initiateSwAlbranch(mo,mo.MinusMat(mo,Swbranch[1],Swbranch[3],Div,Div),mo.AddMat(mo,Swbranch[6],Swbranch[7],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,mo.AddMat(mo,Swbranch[8],Swbranch[14],Div,Div),mo.MinusMat(mo,Swbranch[11],Swbranch[12],Div,Div),Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,Swbranch[10],Swbranch[12],Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,Swbranch[9],Swbranch[11],Div,Div);
		I++;
		Swbranch[I] = mo.AddMat(mo,mo.AddMat(mo,Swbranch[10],Swbranch[13],Div,Div),mo.MinusMat(mo,Swbranch[8],Swbranch[9],Div,Div),Div,Div);
		I++;
		Dbranch[d] = mo.EndSeqSwbranch(Swbranch,mo,Div,Div);
		fractal++;
	}
	return Dbranch[d];
}
vector<vector<double>> MatrixOp::EndSeqSwbranch(vector<vector<vector<double>>> Swbranch, MatrixOp& mo, int i, int j){
	if (Recursive){cout << "Recursive";}
	vector<vector<double>> J;
	J.resize(2*i);
	for (int x=0;x<2*j;x++){
			J[x].resize(2*j);
		}
	for (int f=0;f<=3;f++){
		if(f==0){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					J[u][c] = J[u][c]+Swbranch[15][u][c];
				}
			}
		}
		if(f==1){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					J[u][c+j] = J[u][c+j]+Swbranch[16][u][c];
				}
			}
		}
		if(f==2){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					J[u+i][c] = J[u+i][c]+Swbranch[17][u][c];
				}
			}
		}
		if(f==3){
			for(int u=0;u<i;u++){
				for(int c=0;c<j;c++){
					J[u+i][c+j] = J[u+i][c+j]+Swbranch[18][u][c];
				}
			}
		}
	}
	if(Recursive){
		cout << endl;
		cout << "Is. . .."<< endl;
		int z;
		cout << Return;
		cout << MoveCol;
		cin >> z;
		if(Return!=Stop){
			cout << "#";
			if(MoveCol<Return){
				cout << "___";
				Stop2++;
				Start=false;
				initiateSwAl(mo,D,i,j);
			}
			else{
				if(MoveRaw<Return){
					Stop++;
					Stop2=0;
					initiateSwAl(mo,D,i,j);
				}
			}
		}
	}
	return J;
}
