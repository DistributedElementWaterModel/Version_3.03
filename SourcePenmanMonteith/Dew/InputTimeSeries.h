#pragma once

class DateTime;
class ParametersGeneral;

#include "stdafx.h"

using namespace std;

class InputTimeSeries
{
public:
    InputTimeSeries(int numberRows, int numberColumns, DateTime firstTime, DateTime lastTime, int secondsPerTimeStep);
    ~InputTimeSeries();
    int GetNumberTimeSteps() const;
    int GetNumberInputSeries() const;
    void SetGeneralPar(ParametersGeneral *parObj);
    ParametersGeneral *GetGeneralPar() const;
    DateTime GetDateTime(int i) const;
    void SetInput(ifstream &fin);
    double GetInput(int i, int j) const;
    void WriteInput();

private:
    ParametersGeneral *commonPar;
    int timeSteps;
    int numberSeries;
    DateTime * datetime;
    double * inputArray;
};
