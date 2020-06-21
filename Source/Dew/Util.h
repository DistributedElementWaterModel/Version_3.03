#pragma once
double power(double base, double exponent);
int leapYear(int year);
int dayNumber(int year, int month, int day);
void dayNo2Date(int dayNo, int year, int * month, int * day);
double potentialEvapTemperatureIndex(double temp, double epotPar);
double HBVTranspSoilEvapTemperatureIndex(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel);
double KiWaTranspSoilEvapTemperatureIndex(double teta, double temp, double epotPar, double eactPar, double tSat0, double wiltPoint);
double potentialEvapLongTermMean(double temp, double potentialEvaporation);
double HBVTranspSoilEvapLongTermMean(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel, double potentialEvaporation);
double KiWaTranspSoilEvapLongTermMean(double teta, double temp, double epotPar, double eactPar, double tSat0, double wiltPoint, double potentialEvaporation);
double absValue(double);


