#pragma once
#include "Data.h"
#include <functional>

double legendreFunction(int n, double x) {
    if (n == 0)
        return 1;
    else if (n == 1)
        return x;
    else
        return ((2.0 * n - 1.0) * x * legendreFunction(n - 1, x) - (n - 1.0) * legendreFunction(n - 2, x)) / n;
}

vector<vector<int>> getCombinations(int n, int k) {   
    if (k == 1) {
        return { {n} };
    }
    if (k == 2) {
        vector<vector<int>> tmp = vector<vector<int>>(n - 1, vector<int>(2));
        for (int i = 0; i < n - 1; i++) {
            tmp[i][0] = n - i - 1;
            tmp[i][1] = i + 1;
        }
        return tmp;
    }
    else {
        vector<vector<int>> tmp;
        vector<vector<int>> tmp2;
        vector<int> tmp3;
        for (int i = 1; i <= n - (k - 1); i++) {
            tmp2 = getCombinations(n - i, k - 1);
            tmp3 = { i };
            for (unsigned int j = 0; j < tmp2.size(); j++) {
                tmp2[j].insert(tmp2[j].begin(), tmp3.begin(), tmp3.end());
            }
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
        }
        return tmp;
    }
}

void loopForOneLegendrePolynomial(Data& data, vector<int>& powerList, int size, vector<int>& steps,
                                   int startStepBefore, int power, const std::function<void(vector<int>&)>& leafProcessor) {
    unsigned int startStep = startStepBefore + 1;
    for (unsigned int i = startStep; i < data.stepsForMemory[power - 1]; i++) {
        steps[powerList.size() - size] = i;
        data.startSteps[powerList.size() - size] = startStep;
        if (size > 1) {
            loopForOneLegendrePolynomial(data, powerList, size - 1, steps, i, power, leafProcessor);
        }
        else {
            leafProcessor(steps);
        }
    }
}

void writeResultLine(Data& data, ofstream& dataFile, vector<int>& steps, vector<int>& powerList, double C_T) {
    for (unsigned int j = 0; j < steps.size(); j++) {
        dataFile << steps[j] << "\t";
    }
    for (unsigned int j = 1; j <= data.maxDegree - steps.size(); j++) {
        dataFile << 0 << "\t";
    }
    for (unsigned int j = 0; j < powerList.size(); j++) {
        dataFile << powerList[j] << "\t";
    }
    for (unsigned int j = 1; j <= data.maxDegree - steps.size(); j++) {
        dataFile << 0 << "\t";
    }
    dataFile << C_T << endl;
}

void calcErrors(Data& data, vector<int>& powerList, ofstream& dataFile, arma::Mat<double>& ssti, int power) {
    vector<int> steps = vector<int>(powerList.size());
    data.startSteps = vector<int>(powerList.size());

    auto leafProcessor = [&](vector<int>& steps) {
        for (unsigned int k = 0; k < data.stateMatrix_Train.n_rows; k++) {
            data.reconstTrain(k) = 1.0;
            for (unsigned int j = 0; j < powerList.size(); j++) {
                data.reconstTrain(k) *= legendreFunction(powerList[j], data.trainingTargetFromInputFile[k + data.startTargetRow - steps[j]]) * sqrt(2.0 * powerList[j] + 1.0); 
            }
        }
        for (unsigned int k = 0; k < data.stateMatrix_Test.n_rows; k++) {
            data.reconstTest(k) = 1.0;
            for (unsigned int j = 0; j < powerList.size(); j++) {
                data.reconstTest(k) *= legendreFunction(powerList[j], data.trainingTargetFromInputFile[k + data.trainLength + data.startTargetRow - steps[j]]) * sqrt(2.0 * powerList[j] + 1.0); 
            }
        }
        
        data.pMatrix = data.stateMatrix_Train.t() * data.reconstTrain;  
        data.pTest = data.stateMatrix_Test.t() * data.reconstTest; 

        arma::Col<double> w = ssti * data.pMatrix;
        arma::Col<double> yPred = data.stateMatrix_Test * w;

        double mse = pow(norm(yPred - data.reconstTest), 2);
        double varZ = dot(data.reconstTest, data.reconstTest);

        double C_T = 1.0 - (mse / varZ);

        data.muylmuyl = data.pMatrix * data.pMatrix.t() / data.trainLength;
        data.muylmuyl_Test = data.pTest * data.pTest.t() / (data.dataLength - data.trainLength);

        if (C_T > 0) {
            data.mumu += data.muylmuyl;
            data.mumu_Test += data.muylmuyl_Test;
        }
        writeResultLine(data, dataFile, steps, powerList, C_T);
    };
    loopForOneLegendrePolynomial(data, powerList, powerList.size(), steps, -1, power, leafProcessor);
}

void calcestimation(Data& data, vector<int>& powerList, ofstream& dataFile, int power) {
    vector<int> steps = vector<int>(powerList.size());
    data.startSteps = vector<int>(powerList.size());

    auto leafProcessor = [&](vector<int>& steps) {
        for (unsigned int k = 0; k < data.stateMatrix_Train.n_rows; k++) {
            data.targetTrain(k) = 1.0;
            for (unsigned int j = 0; j < powerList.size(); j++) {
                data.targetTrain(k) *= legendreFunction(powerList[j], data.trainingTargetFromInputFile[k + data.startTargetRow - steps[j]]) * sqrt(2.0 * powerList[j] + 1.0); 
            }
        }

        for (unsigned int k = 0; k < data.stateMatrix_Test.n_rows; k++) {
            data.targetTest(k) = 1.0;
            for (unsigned int j = 0; j < powerList.size(); j++) {
                data.targetTest(k) *= legendreFunction(powerList[j], data.trainingTargetFromInputFile[k + data.trainLength + data.startTargetRow - steps[j]]) * sqrt(2.0 * powerList[j] + 1.0);
            }
        }

        data.pMatrix = data.stateMatrix_Train.t() * data.targetTrain / data.trainLength;  
        data.pTest = data.stateMatrix_Test.t() * data.targetTest / (data.dataLength - data.trainLength);  

        double varZ = dot(data.targetTest, data.targetTest) / (data.dataLength - data.trainLength);

        double C_T = arma::as_scalar((2 * data.pMatrix.t() * data.mumu * data.pTest - data.pMatrix.t() * data.mumu * data.mumu_Test * data.mumu * data.pMatrix) / varZ);

        writeResultLine(data, dataFile, steps, powerList, C_T);
    };

    loopForOneLegendrePolynomial(data, powerList, powerList.size(), steps, -1, power, leafProcessor);
}

void calcNonlinearMemory(Data& data) {
    data.dataLength = data.stateMatrix.n_rows;
    data.trainLength = (1 * data.stateMatrix.n_rows) / 2;

    data.stateMatrix_Train = data.stateMatrix.rows(0, data.trainLength - 1);
    data.stateMatrix_Test = data.stateMatrix.rows(data.trainLength, data.stateMatrix.n_rows - 1);

    data.targetTrain.set_size(data.stateMatrix_Train.n_rows);
    data.targetTest.set_size(data.stateMatrix_Test.n_rows);
    data.reconstTrain.set_size(data.stateMatrix_Train.n_rows);
    data.reconstTest.set_size(data.stateMatrix_Test.n_rows);

    arma::Mat<double> sst;
    arma::Mat<double> ssttest;

    sst = data.stateMatrix_Train.t() * data.stateMatrix_Train;
    ssttest = data.stateMatrix_Test.t() * data.stateMatrix_Test;

    arma::Mat<double> ssti = arma::pinv(sst);

    data.mumu = arma::zeros<arma::Mat<double>>(sst.n_rows, sst.n_cols);
    data.mumu_Test = arma::zeros<arma::Mat<double>>(ssttest.n_rows, ssttest.n_cols);

    ofstream data_file_nonlinear_memory;
    data_file_nonlinear_memory.open(data.directIPCSavefile);

    ofstream data_file_estimate_noiseless;
    data_file_estimate_noiseless.open(data.CROPIPCSavefile);

    for (unsigned int power = 1; power <= data.maxDegree; power++) {
        cout << "Power: " << power << endl;

        for (unsigned int variables = 1; variables <= power; variables++) {
            if (data.debug)
                cout << "Variables: " << variables << endl;

            vector<vector<int>> powerLists = getCombinations(power, variables);

            for (unsigned int combination = 0; combination < powerLists.size(); combination++) {
                calcErrors(data, powerLists[combination], data_file_nonlinear_memory, ssti, power);
            }
        }
    }

    data.mumu += sst.col(0) * sst.row(0) / data.trainLength;
    data.mumu /= data.trainLength;
    data.mumu = arma::pinv(data.mumu);

    data.mumu_Test += ssttest.col(0) * ssttest.row(0) / (data.dataLength - data.trainLength);
    data.mumu_Test /= (data.dataLength - data.trainLength);

    for (unsigned int power = 1; power <= data.maxDegree; power++) {
        cout << "Power: " << power << endl;
        for (unsigned int variables = 1; variables <= power; variables++) {
            vector<vector<int>> powerLists = getCombinations(power, variables);
            for (unsigned int combination = 0; combination < powerLists.size(); combination++) {
                calcestimation(data, powerLists[combination], data_file_estimate_noiseless, power);
            }
        }
    }

    data_file_nonlinear_memory.close();
    data_file_estimate_noiseless.close();
}
