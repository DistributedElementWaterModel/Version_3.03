#include "MeteorologicalStations.h"


MeteorologicalStations::MeteorologicalStations() :
    numberPrecStations(0),
    numberTempStations(0)
{
}

MeteorologicalStations::~MeteorologicalStations()
{
}

void  MeteorologicalStations::SetNumPrecStations(int value)
{
    numberPrecStations = value;
}
int  MeteorologicalStations::GetNumPrecStations()
{
    return numberPrecStations;
}
void  MeteorologicalStations::SetNumTempStations(int value)
{
    numberTempStations = value;
}
int  MeteorologicalStations::GetNumTempStations()
{
    return numberTempStations;
}
int  MeteorologicalStations::GetStationNumber(int k)
{
    return stationNumber[k];
}
double  MeteorologicalStations::GetStationCoordX(int k)
{
    return stationCoordX[k];
}
double  MeteorologicalStations::GetStationCoordY(int k)
{
    return stationCoordY[k];
}
double  MeteorologicalStations::GetStationAltitude(int k)
{
    return stationAltitude[k];
}
void SetMeteorologicalStations(ifstream &fileControl, ofstream &fout);

void MeteorologicalStations::SetMeteorologicalStations(ifstream &fileControl, ofstream &fout)
{
    int i, numberPrecStations, numberTempStations;
    char stationType;
    char fileName[80];
    /*  cout << " File with meteorological stations: ";
    cin >> fileName;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    //  cout << fileName << endl;
    ifstream finMetSta(fileName);  // Open for reading
    if (!finMetSta.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        //exit(1);
    }
    finMetSta.ignore(100, ':');
    finMetSta >> numberPrecStations;
    finMetSta.ignore(100, ':');
    finMetSta >> numberTempStations;
    SetNumPrecStations(numberPrecStations);
    SetNumTempStations(numberTempStations);
    stationNumber = new int[numberPrecStations + numberTempStations];
    stationCoordX = new double[numberPrecStations + numberTempStations];
    stationCoordY = new double[numberPrecStations + numberTempStations];
    stationAltitude = new double[numberPrecStations + numberTempStations];
    for (i = 0; i < numberPrecStations; i++)
    {
        finMetSta >> stationType >> stationNumber[i] >> stationCoordX[i] >> stationCoordY[i] >> stationAltitude[i];
        finMetSta.ignore(1024, '\n');
        if (stationType != 'P' && stationType != 'p')
        {
            cout << endl << " Not legal station type : " << stationType << " , Station number = " << i << endl << endl;
            //exit(1);
        }
    }
    for (i = 0; i < numberTempStations; i++)
    {
        finMetSta >> stationType >> stationNumber[numberPrecStations + i] >> stationCoordX[numberPrecStations + i]
                  >> stationCoordY[numberPrecStations + i] >> stationAltitude[numberPrecStations + i];
        finMetSta.ignore(1024, '\n');
        if (stationType != 'T' && stationType != 't')
        {
            cout << endl << " Not legal station type : " << stationType
                 << " station number = " << numberPrecStations + i << endl << endl;
            //exit(1);
        }
    }
    finMetSta.close();
    fout << endl << "Meteorological stations: \n";
    fout << GetNumPrecStations() << endl;
    fout << GetNumTempStations() << endl;
    fout << "Prec.    " << endl;
    for (i = 0; i < GetNumPrecStations(); i++)
    {
        fout << GetStationNumber(i) << "  ";
        fout << GetStationCoordX(i) << "  ";
        fout << GetStationCoordY(i) << "  ";
        fout << GetStationAltitude(i) << endl;
    }
    fout << "Temp.    " << endl;
    for (i = 0; i < GetNumTempStations(); i++)
    {
        fout << GetStationNumber(GetNumPrecStations() + i) << "  ";
        fout << GetStationCoordX(GetNumPrecStations() + i) << "  ";
        fout << GetStationCoordY(GetNumPrecStations() + i) << "  ";
        fout << GetStationAltitude(GetNumPrecStations() + i) << endl;
    }
    fout << endl;
}
