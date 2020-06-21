#include "SelectedKiWaTimeSeriesElements.h"

#include "stdafx.h"

using namespace std;

SelectedKiWaTimeSeriesElements::SelectedKiWaTimeSeriesElements() :
    numberElements(0)
{
}

SelectedKiWaTimeSeriesElements::~SelectedKiWaTimeSeriesElements()
{
}

void  SelectedKiWaTimeSeriesElements::SetNumberElements(int value)
{
    numberElements = value;
}
int  SelectedKiWaTimeSeriesElements::GetNumberElements() const
{
    return numberElements;
}
void  SelectedKiWaTimeSeriesElements::SetKiWaTimeSeriesElement(int index, int value)
{
    timeSeriesElements[index] = value;
}
int  SelectedKiWaTimeSeriesElements::GetKiWaTimeSeriesElement(int index) const
{
    return timeSeriesElements[index];
}
void  SelectedKiWaTimeSeriesElements::SetLengthFractionOne(int index, double value)
{
    lengthFractionOne[index] = value;
}
double  SelectedKiWaTimeSeriesElements::GetLengthFractionOne(int index) const
{
    return lengthFractionOne[index];
}
void  SelectedKiWaTimeSeriesElements::SetLengthFractionTwo(int index, double value)
{
    lengthFractionTwo[index] = value;
}
double  SelectedKiWaTimeSeriesElements::GetLengthFractionTwo(int index) const
{
    return lengthFractionTwo[index];
}

void SelectedKiWaTimeSeriesElements::SelectedKiWaTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout)
{
    char fileName[80], buffer[1024];
    int i, numElements, elementId;
    double lthFracOne, lthFracTwo;

    /*  cout << " File with landscape elements selected for KiWa time series output: ";
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
    lengthFractionOne = new double[numElements];
    lengthFractionTwo = new double[numElements];
    fileTimeSeries.ignore(1024, '\n');
    fileTimeSeries.getline(buffer, 1024);
    for (i = 0; i < numElements; i++)
    {
        fileTimeSeries >> elementId >> lthFracOne >> lthFracTwo;
        if (lthFracOne < 0.0 || lthFracOne > 1.0)
        {
            lthFracOne = 0.5;
        }
        if (lthFracTwo < 0.0 || lthFracTwo > 1.0)
        {
            lthFracTwo = 1.0;
        }
        SetKiWaTimeSeriesElement(i, elementId);
        SetLengthFractionOne(i, lthFracOne);
        SetLengthFractionTwo(i, lthFracTwo);
        fileTimeSeries.ignore(1024, '\n');
    }
    fout << endl << "Number of landscape elements selected for KiWa time series output: " << GetNumberElements() << endl;
    for (i = 0; i < GetNumberElements(); i++)
    {
        fout << GetKiWaTimeSeriesElement(i) << "  " << GetLengthFractionOne(i) << "  " << GetLengthFractionTwo(i) << endl;
    }
    fout << endl;
}
