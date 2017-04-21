#pragma once
#include "head.h"
#include <vector>
#include <iostream>
using namespace std;
int Sub(vector<vector<double>> A, int g,int o,int a1,int a2);
MatrixOp MO;
int n,m,j,Input;
bool Base;
bool Rnd=true, Normalized=true, Rec=false;
void main(){
	j = 0;
	int g = 0;
	int o=0;
	int a1 =0;
	int a2 =0;
	int t=0;
	string k="N";
	double tr=0;
	if ( j == 0 ){
	cout << "This Simple Class execute Matrix multiplication recursivly" << endl << "https://en.wikipedia.org/wiki/Strassen_algorithm" << endl;
	cout << "Please enter The metrix dimention;" << endl << "|1|-(100,100), |2| -(200,200)" << endl;
	cout << ":";}
	while ( Input != 1 || Input != 0 ){
		cin >> Input;
		if ( Input == 1 ) {n = 100; m = 100;break;}
		else if ( Input == 2 ) {n = 200; m = 200;break;}
		else if ( cin.fail()){
			cin.clear();
			cin.ignore();
			if ( t == 0 ){cout<< "wrong input please try again:"<<endl;t++;}
		}
	}
	cout << endl;
	vector<vector<double>> A = MO.defineMatrix(n,m);
	Sub(A,g,o,a1,a2);
}
int Sub(vector<vector<double>> A, int g, int o,int a1,int a2){
	if (g==0){
		MO.InitiateParam(g,a1,a2);
		MO.InsertationmethodaMat(MO, A,Rnd,Normalized,n,m);
		g++;
		o++;
	}
	if(g!=0){
		MO.InitiateParam(g,a1,a2);
		if (o>=2){
			int N=0;
			if (!Rec){
				MO.NotRecursive(MO,A,n,m);}
			if(Rec){
				MO.Recursion(MO,G,n,m);
			}
		}
		MO.InsertationmethodaMat(MO, A,Rnd,Normalized,n,m);
		g++;
		o++;
	}
	if (g<3){Sub(A,g,o,n,m);}
	return 1;
}

