#include "SelectedSubCatchmentTimeSeriesElements.h"


SelectedSubCatchmentTimeSeriesElements::SelectedSubCatchmentTimeSeriesElements() :
    numberElements(0)
{
}

SelectedSubCatchmentTimeSeriesElements::~SelectedSubCatchmentTimeSeriesElements()
{
}

void  SelectedSubCatchmentTimeSeriesElements::SetNumberElements(int value)
{
    numberElements = value;
}
int  SelectedSubCatchmentTimeSeriesElements::GetNumberElements() const
{
    return numberElements;
}
void  SelectedSubCatchmentTimeSeriesElements::SetSubCatchmentTimeSeriesElement(int index, int value)
{
    timeSeriesElements[index] = value;
}
int  SelectedSubCatchmentTimeSeriesElements::GetSubCatchmentTimeSeriesElement(int index) const
{
    return timeSeriesElements[index];
}

void  SelectedSubCatchmentTimeSeriesElements::SelectedSubCatchmentTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout)
{
    char fileName[80], buffer[1024];
    char ch;
    int i, j, numElements, elementId;

    /*  cout << " File with sub-catchment elements selected for time series output: ";
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
    fileTimeSeries.ignore(1024, '\n');
    for (i = 0; i < numElements; i++)
    {
        fileTimeSeries >> j >> ch >> elementId;
        if (j != i)
        {
            cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << j << endl;
            //exit(1);
        }
        SetSubCatchmentTimeSeriesElement(i, elementId);
        fileTimeSeries.ignore(1024, '\n');
    }
    fout << "Number of sub-catchment elements selected for time series output: " << GetNumberElements() << endl;
    for (i = 0; i < GetNumberElements(); i++)
    {
        fout << GetSubCatchmentTimeSeriesElement(i) << endl;
    }
    fout << endl;
}
