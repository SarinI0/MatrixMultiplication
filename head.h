#pragma once
	#include <vector>
	#include <map>
	#include <iostream>
	#include <math.h>
	#include "head.h"
	#include <vector>
	#include <iostream>
	#include <ctime>
	#include <cstdlib>
	#include <fstream>
	#include <string>
	#include <iomanip>
	#include <windows.h>
	#include <locale>
	#include <sstream>
using namespace std;
//MatrixOp
static bool Recursive, Start, Unit;
static int OriginI;
static int OriginJ;
static int Return,Stop,Stop2,MoveRaw,MoveCol,d,fractal;
static vector<vector<double>> G,P,D;
class MatrixOp{
private:

public:
	int MatrixOp::PrintVecint(vector<int>& u,int n);
	vector<vector<double>> defineMatrix(int i, int j);
	vector<double> defineVector(int n);
	int PrintMat(vector<vector<double>>& A,int n,int m);
	int PrintVec(vector<double>& u,int n);
	double tr(vector<vector<double>>& A,int n,int m,double tr);
	bool MatrixOp::InsertationmethodaMat(MatrixOp& mo,vector<vector<double>>& A,bool Rnd,bool Normalized,int i, int j);
	vector<vector<double>> MatrixOp::AddMat(MatrixOp& mo,vector<vector<double>>& O,vector<vector<double>>& M,int i, int j);
	vector<vector<double>> MatrixOp::AddMatsw(MatrixOp& mo,vector<vector<double>>& O,vector<vector<double>>& M,int i, int j);
	vector<double> MatrixOp::Insertationmethodavec(MatrixOp& mo,vector<double>& u,bool Rnd,bool Normalized,int i);
	int MatrixOp::Count1,Count2,a1,a2,a3,a4;
	int MatrixOp::InitiateParam(int Vic,int k1, int k2);
	vector<vector<double>> MatrixOp::Primitivemulti(MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int i, int j);
	vector<vector<double>> MatrixOp::EndSeqSw(vector<vector<vector<double>>>& Sw, MatrixOp& mo, int i, int j);
	vector<vector<double>> MatrixOp::SwAlmulti(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int N);
	vector<vector<double>> MatrixOp::initiateSwAl(MatrixOp& mo, vector<vector<double>>& A, int i, int j);
	vector<vector<double>> MatrixOp::Recursion(MatrixOp& mo, vector<vector<double>>& A, int i, int j);
	vector<vector<double>> MatrixOp::SeqSwAl(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int N);
	vector<vector<vector<double>>> Sw;
	vector<vector<vector<double>>> Swbranch;
	vector<vector<vector<double>>> Dbranch;
	vector<vector<double>> MatrixOp::Insertationmethodasw(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A ,int I,int R,int C);
	vector<vector<double>> MatrixOp::MinusMat(MatrixOp& mo, vector<vector<double>>& O,vector<vector<double>>& M,int i, int j);
	vector<vector<double>> MatrixOp::NotRecursive(MatrixOp& mo, vector<vector<double>>& A, int i, int j);
	vector<vector<double>> MatrixOp::SwAlmultiRec(vector<vector<vector<double>>>& Sw, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int N);
	vector<vector<double>> MatrixOp::initiateSwAlbranch(MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int i, int j);
	vector<vector<double>> MatrixOp::SeqSwAlbranch(vector<vector<vector<double>>> Swbranch, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int I);
	vector<vector<double>> MatrixOp::EndSeqSwbranch(vector<vector<vector<double>>> Swbranch, MatrixOp& mo, int i, int j);
	vector<vector<double>> MatrixOp::SwAlmultibranch(vector<vector<vector<double>>> Swbranch, MatrixOp& mo, vector<vector<double>>& A, vector<vector<double>>& B, int Div, int lim,int N);
	MatrixOp();
	~MatrixOp();
};
