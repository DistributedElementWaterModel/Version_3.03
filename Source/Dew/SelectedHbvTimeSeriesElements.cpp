#include "SelectedHbvTimeSeriesElements.h"
#include "stdafx.h"

using namespace std;

SelectedHbvTimeSeriesElements::SelectedHbvTimeSeriesElements() :
    numberElements(0)
{
}

SelectedHbvTimeSeriesElements::~SelectedHbvTimeSeriesElements()
{
}

void  SelectedHbvTimeSeriesElements::SetNumberElements(int value)
{
    numberElements = value;
}
int  SelectedHbvTimeSeriesElements::GetNumberElements() const
{
    return numberElements;
}
void  SelectedHbvTimeSeriesElements::SetHbvTimeSeriesElement(int index, int value)
{
    timeSeriesElements[index] = value;
}
int  SelectedHbvTimeSeriesElements::GetHbvTimeSeriesElement(int index) const
{
    return timeSeriesElements[index];
}
void  SelectedHbvTimeSeriesElements::SetEffPor(int index, double value)
{
    effectivePorosity[index] = value;
}
double  SelectedHbvTimeSeriesElements::GetEffPor(int index) const
{
    return effectivePorosity[index];
}
void  SelectedHbvTimeSeriesElements::SetGroundWaterRef(int index, double value)
{
    groundWaterReference[index] = value;
}
double  SelectedHbvTimeSeriesElements::GetGroundWaterRef(int index) const
{
    return groundWaterReference[index];
}

void SelectedHbvTimeSeriesElements::SelectedHbvTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout)
{
    char fileName[80], buffer[1024];
    int i, numElements, elementId;
    double effPor, gwRef;

    /*  cout << " File with landscape elements selected for HBV time series output: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream fileTimeSeries(fileName);  // Open for reading
    if (!fileTimeSeries.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        //exit(1);
    }
    fileTimeSeries.ignore(100, ':');
    fileTimeSeries >> numElements;
    SetNumberElements(numElements);
    timeSeriesElements = new int[numElements];
    effectivePorosity = new double[numElements];
    groundWaterReference = new double[numElements];
    fileTimeSeries.ignore(1024, '\n');
    fileTimeSeries.getline(buffer, 1024);
    for (i = 0; i < numElements; i++)
    {
        fileTimeSeries >> elementId >> gwRef >> effPor;
        SetHbvTimeSeriesElement(i, elementId);
        SetGroundWaterRef(i, gwRef);
        SetEffPor(i, effPor);
        fileTimeSeries.ignore(1024, '\n');
    }
    fout << "Number of landscape elements selected for HBV time series output: " << GetNumberElements() << endl;
    for (i = 0; i < GetNumberElements(); i++)
    {
        fout << GetHbvTimeSeriesElement(i) << "  " << GetGroundWaterRef(i) << "  " << GetEffPor(i) << endl;
    }
    fout << endl;
}
