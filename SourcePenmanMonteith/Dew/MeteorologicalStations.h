#pragma once

#include "stdafx.h"

using namespace std;

class MeteorologicalStations
{
public:
    MeteorologicalStations();
    ~MeteorologicalStations();
    void SetNumPrecStations(int value);
    int GetNumPrecStations();
    void SetNumTempStations(int value);
    int GetNumTempStations();
    int GetStationNumber(int k) ;
    double GetStationCoordX(int k) ;
    double GetStationCoordY(int k) ;
    double GetStationAltitude(int k) ;
    void SetMeteorologicalStations(ifstream &fileControl, ofstream &fout);

private:
    int numberPrecStations;
    int numberTempStations;
    int * stationNumber;
    double * stationCoordX;
    double * stationCoordY;
    double * stationAltitude;
};
