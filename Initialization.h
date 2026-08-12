#pragma once
#include "Data.h"
#include <iomanip>

void setOutputFile(Data& data) {
    data.folderName = "";
    data.inputfileName = data.folderName + "input_1.01e6.bin";
    data.inputLength = 11000;

    data.startTargetRow = 10000;
    data.nodeNum = 100; 
    data.filenameBase = "quadsys_0.0005_rep=1_leng=1e3";
    data.statedatafileName = data.folderName + "stateData_" + data.filenameBase + ".bin";
    data.directIPCSavefile = data.folderName + "output_IPC_" + data.filenameBase + ".txt";
    data.CROPIPCSavefile = data.folderName + "output_est_IPC_" + data.filenameBase + ".txt";

    data.stepsForMemory = { 21, 20, 19, 16, 11, 6 }; 
    data.maxDegree = data.stepsForMemory.size();
    }

void loadData(Data& data) {
    data.stateMatrix.load(data.statedatafileName);
    data.stateMatrix.reshape(data.inputLength - data.startTargetRow, data.nodeNum + 1);
    cout << data.directIPCSavefile << endl;
    cout << data.CROPIPCSavefile << endl;

    data.trainingTargetFromInputFile.load(data.inputfileName);
    if (data.trainingTargetFromInputFile.n_elem > data.inputLength) {
        data.trainingTargetFromInputFile.resize(data.inputLength);
    }
    std::cout << "stateMatrix size: " << data.stateMatrix.n_rows << " x " << data.stateMatrix.n_cols << std::endl;
}

void initializeVariables(Data& data) {
    setOutputFile(data);
    loadData(data);
}

