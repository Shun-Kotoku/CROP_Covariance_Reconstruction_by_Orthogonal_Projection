#pragma once
#include <string>
#include <armadillo>
#include <fstream> 

using namespace std;

struct Data {	
	string folderName; 
	string inputfileName; 
	int inputLength;

	int startTargetRow;
	int nodeNum;
	string filenameBase;
	string statedatafileName;
	string directIPCSavefile; 
	string CROPIPCSavefile; 

	vector <unsigned int> stepsForMemory;
	unsigned int maxDegree;

	arma::Mat<double> stateMatrix;

	arma::Col<double> trainingTargetFromInputFile;


	double dataLength;
	int trainLength;

	arma::Mat<double> stateMatrix_Train;
	arma::Mat<double> stateMatrix_Test;

	arma::Col<double> targetTrain;
	arma::Col<double> targetTest;
	arma::Col<double> reconstTrain;
	arma::Col<double> reconstTest;

	arma::Mat<double> mumu;
	arma::Mat<double> mumu_Test;

	arma::Mat<double> muylmuyl;
	arma::Mat<double> muylmuyl_Test;

	vector<int> startSteps;

	arma::Mat<double> pMatrix;
	arma::Mat<double> pTest;

	bool debug = false;
	Data();
};

Data::Data()
{
}
