#include "energyDensity.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
using namespace std;

fstream ede;
ofstream edg;
fstream fve;

void energyDensity::setConst(double sig, double tN, double tF, double tS)
{
	sigma = sig;
	tau0 = tN;
	tauFinal = tF;
	tauStep = tS;
}
void energyDensity::setGrid(double xR, double xS, double yR, double yS)
{
	xRange = xR;
	xScale = xS;
	yRange = yR;
	yScale = yS;
}
void energyDensity::setPart(double xN, double yN)
{
	partPos.push_back(make_pair(xN, yN));
	partNum = partPos.size();
}
void energyDensity::getParams()
{
	cout << "sigma = " << sigma << endl;
	cout << "tau0 = " << tau0 << " and final tau = " << tauFinal << " with a step of " << tauStep << endl;
	cout << "number of particles = " << partNum << endl;
	cout << "x range = " << xRange << " with a scale of " << xScale << endl;
	cout << "y range = " << yRange << " with a scale of " << yScale << endl;
	cout << "the particle position(s) are:\n";
	for (int part = 0; part < partNum; part++)
	{
		cout << "(" << partPos[part].first << ", " << partPos[part].second << ")\n";
	}

}

Eigen::MatrixXd energyDensity::getEMTensor(double x, double y, double tau)
{
	Eigen::Matrix<double, 4, 4> A;
	Eigen::Matrix<double, 4, 4> B;
	Eigen::Matrix<double, 4, 4> C;
	Eigen::Matrix<double, 4, 4> T;
	B(0, 0) = 1;
	B(1, 1) = B(2, 2) = B(3, 3) = -1;
	B(0, 1) = B(0, 2) = B(0, 3) = 
		B(1, 0) = B(1, 2) = B(1, 3) = 
		B(2, 0) = B(2, 1) = B(2, 3) = 
		B(3, 0) = B(3, 1) = B(3, 2) = 0;
	T(0, 0) = T(0, 1) = T(0, 2) = T(0, 3) =
		T(1, 0) = T(1, 1) = T(1, 2) = T(1, 3) =
		T(2, 0) = T(2, 1) = T(2, 2) = T(2, 3) =
		T(3, 0) = T(3, 1) = T(3, 2) = T(3, 3) = 0;

	for (int point = 0; point < partNum; point++)
	{
		double x0 = partPos[point].first;
		double y0 = partPos[point].second;
		double X = ((x - x0) * (tau - tau0)) / (sigma * sigma);
		double Y = ((y - y0) * (tau - tau0)) / (sigma * sigma);
		double Z = sqrt((X * X) + (Y * Y));
		double sin1 = X / (Z + 1e-16);
		double cos1 = Y / (Z + 1e-16);
		double sin2 = 2 * sin1 * cos1;
		double cos2 = (cos1 * cos1) - (sin1 * sin1);
		double I0 = std::cyl_bessel_i(0, Z);
		double I1 = std::cyl_bessel_i(1, Z);
		double I2 = std::cyl_bessel_i(2, Z);
		double ex = exp(((-1 * pow(x - x0, 2)) - pow(y - y0, 2) - pow(tau - tau0, 2)) / (2 * pow(sigma, 2)));
	
		A(0, 0) = ex * I0;
		A(0, 1) = ex * I1 * sin1;
		A(1, 0) = A(0, 1);
		A(0, 2) = ex * I1 * cos1;
		A(2, 0) = A(0, 2);
		A(0, 3) = 0;
		A(3, 0) = 0;
		A(1, 3) = 0;
		A(3, 1) = 0;
		A(2, 3) = 0;
		A(3, 2) = 0;
		A(3, 3) = 0;
		A(1, 1) = ex * ((I0 * sin1 * sin1) + ((I1 * cos2) / (Z + 1e-16)));
		A(1, 2) = ex * 0.5 * I2 * sin2;
		A(2, 1) = A(1, 2);
		A(2, 2) = ex * ((I0 * cos1 * cos1) - ((I1 * cos2) / (Z + 1e-16)));

		C = A * B;
		T = T + C;
	}
	return T;
}
double energyDensity::getEnergyDensity(double x, double y, double tau)
{
	Eigen::Matrix<double, 4, 4> D;
	D = getEMTensor(x, y, tau);
	Eigen::EigenSolver<Eigen::Matrix<double, 4, 4> > s(D);  // the instance s(A) includes the eigensystem

	return real(s.eigenvalues()[0]);
}
pair <double, double> energyDensity::getFlowVelocityXY(double x, double y, double tau)
{
	Eigen::Matrix<double, 4, 4> T;
	T = getEMTensor(x, y, tau);
	Eigen::EigenSolver<Eigen::Matrix<double, 4, 4> > s(T);
	return pair <double, double> (real(s.eigenvectors()(1, 0)), real(s.eigenvectors()(2, 0)));
}
double energyDensity::getFlowVelocityX(double x, double y, double tau)
{
	return getFlowVelocityXY(x, y, tau).first;
}
double energyDensity::getFlowVelocityY(double x, double y, double tau)
{
	return getFlowVelocityXY(x, y, tau).second;
}
void energyDensity::EDGrid(double tau)
{
	edg.open("energyDensityGrid_" + to_string(tau) + ".txt");
	double percent = 0;
	for (double x = -1 * xScale * xRange; x < (xRange + 1) * xScale; x += xScale)
	{
		for (double y = -1 * yScale * yRange; y < (yRange + 1) * yScale; y += yScale)
		{
			edg << x << " " << y << " " << real(getEnergyDensity(x, y, tau)) << endl;
			cout << (percent / ((((xRange * 2) + 1))*((yRange * 2) + 1))) * 100 << "% for tau = " << tau << endl;
			percent++;
		}
	}
	edg.close();
}
void energyDensity::EDEvolution(double x, double y)
{
	ede.open("energyDensityEvolution.txt");
	Eigen::Matrix<double, 4, 4> A;
	Eigen::Matrix<double, 4, 4> B;
	Eigen::Matrix<double, 4, 4> C;
	Eigen::Matrix<double, 4, 4> T;
	B(0, 0) = 1;
	B(1, 1) = B(2, 2) = B(3, 3) = -1;
	B(0, 1) = B(0, 2) = B(0, 3) = B(1, 0) = B(1, 2) = B(1, 3) = B(2, 0) = B(2, 1) = B(2, 3) = B(3, 0) = B(3, 1) = B(3, 2) = 0;
	ede << "(" << x << ", " << y << ")";
	for (double tau = tau0; tau <= tauFinal; tau++)
	{
		ede << " " << real(getEnergyDensity(x, y, tau));
	}
	ede << endl;
	ede.close();
}
void energyDensity::EDEvolution()
{
	ede.open("energyDensityEvolution.txt");
	Eigen::Matrix<double, 4, 4> A;
	Eigen::Matrix<double, 4, 4> B;
	Eigen::Matrix<double, 4, 4> C;
	Eigen::Matrix<double, 4, 4> T;
	B(0, 0) = 1;
	B(1, 1) = B(2, 2) = B(3, 3) = -1;
	B(0, 1) = B(0, 2) = B(0, 3) = B(1, 0) = B(1, 2) = B(1, 3) = B(2, 0) = B(2, 1) = B(2, 3) = B(3, 0) = B(3, 1) = B(3, 2) = 0;
	for (double x = -1 * xScale * xRange; x < (xRange + 1) * xScale; x += xScale)
	{
		for (double y = -1 * yScale * yRange; y < (yRange + 1) * yScale; y += yScale)
		{
			ede << "(" << x << ", " << y << ")";
			for (double tau = tau0; tau <= tauFinal; tau++)
			{
				ede << " " << real(getEnergyDensity(x, y, tau));
			}
			ede << endl;
		}
	}
	ede << endl;
	ede.close();
}
void energyDensity::FVEvolution()
{
	fve.open("flowVelocityEvolution.txt");
	for (double x = -1 * xScale * xRange; x < (xRange + 1) * xScale; x += xScale)
	{
		for (double y = -1 * yScale * yRange; y < (yRange + 1) * yScale; y += yScale)
		{
			fve << "(" << x << ", " << y << ") ";
			for (int tau = 0; tau < 30; tau++)
			{
				fve << "(" << getFlowVelocityX(x, y, tau) << ", " << getFlowVelocityY(x, y, tau) << ") ";
			}
			fve << endl;
		}
	}
	fve.close();
}
energyDensity::energyDensity()
{
	sigma = 1;
	tau0 = 0;
	tauFinal = 1;
	tauStep = 1;
	xRange = yRange = xScale = yScale = 1;
	partNum = 0;
}
energyDensity::energyDensity(double sig,
	double tN, double tF, double tS,
	double xR, double xS, double yR, double yS)
{
	partNum = 0;
	sigma = sig;
	tau0 = tN;
	tauFinal = tF;
	tauStep = tS;
	xRange = xR;
	xScale = xS;
	yRange = yR;
	yScale = yS;
}