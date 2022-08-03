#include <vector>
#include <Eigenvalues>
#include <utility>
using namespace std;

class energyDensity
{
private:
	int partNum;
	double sigma;
	double tau0;
	double tauFinal;
	double tauStep;
	double xRange;
	double xScale;
	double yRange;
	double yScale;
	vector<pair<double, double>> partPos;
public:
	void setConst(double sig, double tN, double tF, double tS);
	void setGrid(double xR, double xS, double yR, double yS);
	void setPart(double xN, double yN);
	void getParams();
	Eigen::MatrixXd getEMTensor(double x, double y, double tau);
	double getEnergyDensity(double x, double y, double tau);
	pair <double, double> getFlowVelocityXY(double x, double y, double tau);
	double getFlowVelocityX(double x, double y, double tau);
	double getFlowVelocityY(double x, double y, double tau);
	void EDGrid(double tau);
	void EDEvolution(double x, double y);
	void EDEvolution();
	void FVEvolution();
	energyDensity();
	energyDensity(double sig,
		double tN, double tF, double tS,
		double xR, double xS, double yR, double yS);
};