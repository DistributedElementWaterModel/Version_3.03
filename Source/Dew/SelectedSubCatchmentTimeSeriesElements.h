#pragma once

#include "stdafx.h"

using namespace std;

class SelectedSubCatchmentTimeSeriesElements
{
public:
    SelectedSubCatchmentTimeSeriesElements();
    ~SelectedSubCatchmentTimeSeriesElements();
    void SetNumberElements(int value);
    int GetNumberElements() const;
    void SetSubCatchmentTimeSeriesElement(int index, int value);
    int GetSubCatchmentTimeSeriesElement(int index) const;
    void SelectedSubCatchmentTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout);

private:
    int numberElements;
    int * timeSeriesElements;
};
