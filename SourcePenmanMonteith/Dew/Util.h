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
double HBVTranspSoilEvap(double soilMoist, double potev, double fieldCapacity, double fcDel);
double KiWaTranspSoilEvapLongTermMean(double teta, double temp, double epotPar, double eactPar, double tSat0, double wiltPoint, double potentialEvaporation);
double absValue(double);
double potentialEvap(double temp, double tempMax, double tempMin, double radiationS, double vp, 
	double wind1, int dayofyear, double elemLatitute, double elemElevation, double windHeight, double humidityHeight,
	double treeHeight, double treeLai, double treeLaiCorr, double treeAlbedo, double bulkResist,double albedo_snow,double deciduous, 
        double snowCoverFraction, double tree_wind_H, double Topen_min, double Tclose_min, double VPDclose,
		     double VPDopen, double gh, double cl, double z0g,
		     double gsmax, double CR, double D50, double Q50);
double potentialEvapLake(double temp, double tempMax, double tempMin, double radiationS, double vp, 
	double wind1, int dayofyear, double elemLatitute, double elemElevation,double windHeight);


