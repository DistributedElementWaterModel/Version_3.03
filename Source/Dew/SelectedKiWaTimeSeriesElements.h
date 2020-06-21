#pragma once

#include "stdafx.h"

using namespace std;

class SelectedKiWaTimeSeriesElements
{
public:
    SelectedKiWaTimeSeriesElements();
    ~SelectedKiWaTimeSeriesElements();
    void SetNumberElements(int value);
    int GetNumberElements() const;
    void SetKiWaTimeSeriesElement(int index, int value);
    int GetKiWaTimeSeriesElement(int index) const;
    void SetLengthFractionOne(int index, double value);
    double GetLengthFractionOne(int index) const;
    void SetLengthFractionTwo(int index, double value);
    double GetLengthFractionTwo(int index) const;
    void SelectedKiWaTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout);

private:
    int numberElements;
    int * timeSeriesElements;
    double * lengthFractionOne;
    double * lengthFractionTwo;
};
