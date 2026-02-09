#pragma once

#include "stdafx.h"

using namespace std;

class SelectedHbvTimeSeriesElements
{
public:
    SelectedHbvTimeSeriesElements();
    ~SelectedHbvTimeSeriesElements();
    void SetNumberElements(int value);
    int GetNumberElements() const;
    void SetHbvTimeSeriesElement(int index, int value);
    int GetHbvTimeSeriesElement(int index) const;
    void SetEffPor(int index, double value);
    double GetEffPor(int index) const;
    void SetGroundWaterRef(int index, double value);
    double GetGroundWaterRef(int index) const;
    void SelectedHbvTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout);

private:
    int numberElements;
    int * timeSeriesElements;
    double * effectivePorosity;
    double * groundWaterReference;
};
