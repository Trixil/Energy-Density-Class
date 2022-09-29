#define __STDCPP_WANT_MATH_SPEC_FUNCS__
#define _USE_MATH_DEFINES

#include <Eigenvalues>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "energyDensity.h"
using namespace std;
int main()
{
	energyDensity test;
	test.setConst(0.2, 0, 2, 0.5);
	test.setGrid(50, 0.05, 50, 0.05);
	//test.setPart(0, 0);
	test.setPart(0,0);
	/*
	test.setPart(-11, 6);
	test.setPart(1, -10);
	test.setPart(1, -10);
	test.setPart(1, -10);
	test.setPart(3, -3);
	test.setPart(5, 0);
	test.setPart(5, 0);
	test.setPart(1, 6);
	test.setPart(3, 6);
	test.setPart(3, 6);
	test.setPart(-11, 8);
	test.setPart(7, 5);
	*/
	//test.EDGrid(15);
	test.getParams();
	cout << "u_mu[2] is " << test.getu_mu(0, 3, 5)[2] << endl;
	cout << "j_mu[2] is " << test.getj_mu(0, 3, 5)[2] << endl;
	//test.FVSliceX(0.5, 0);
	test.FVSliceX(0.5, 0);
	/*for (double i = 0, percent = 0; i < 30; i+= 0.075, percent++)
	{
		test.EDGrid(i);
		cout << (percent / 400) * 100 << "% for tau = " << i << endl;
	}*/
	//test.FVEvolution();
}