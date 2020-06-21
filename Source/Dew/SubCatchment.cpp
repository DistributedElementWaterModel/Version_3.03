#include "SubCatchment.h"
#include "DateTime.h"
#include "Dew.h"
#include "stdafx.h"

SubCatchment::SubCatchment() :
    lakeNumber(-9999),
    numUpStream(0),
    numLandScape(0),
    maxBas(1.0),
    correction(0.0),
    manning(-9999.0),
    riverSlope(0.0),
    width(0.0)
{
    UpStream = 0;
    SetLandScapeElement(0);
    SetSelectedSubCatchmentTimeSeriesElements(0);
}

SubCatchment::~SubCatchment()
{
}

void  SubCatchment::SetSubCatchmentIndex(int value)
{
    subCatchmentIndex = value;
}
int  SubCatchment::GetSubCatchmentIndex() const
{
    return subCatchmentIndex;
}
void  SubCatchment::SetIdentifier(int value)
{
    identifier = value;
}
int  SubCatchment::GetIdentifier() const
{
    return identifier;
}
int  SubCatchment::GetNumUpStream() const
{
    return numUpStream;
}
void  SubCatchment::SetUpStream(int k, SubCatchment *theSubCatchment)
{
    UpStream[k] = theSubCatchment;
}
SubCatchment * SubCatchment::GetUpStream(int k) const
{
    return UpStream[k];
}
void  SubCatchment::SetNumLandScape(int value)
{
    numLandScape = value;
}
int  SubCatchment::GetNumLandScape() const
{
    return numLandScape;
}
void  SubCatchment::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * SubCatchment::GetLandScapeElement() const
{
    return landScapeElement;
}
void  SubCatchment::SetLakeNumber(int value)
{
    lakeNumber = value;
}
int  SubCatchment::GetLakeNumber() const
{
    return lakeNumber;
}
void  SubCatchment::SetSelectedSubCatchmentTimeSeriesElements(SelectedSubCatchmentTimeSeriesElements *object)
{
    selectedSubCatchmentElements = object;
}
SelectedSubCatchmentTimeSeriesElements * SubCatchment::GetSelectedSubCatchmentTimeSeriesElements() const
{
    return selectedSubCatchmentElements;
}
void  SubCatchment::SetMaxBas(double value)
{
    maxBas = value;
}
double  SubCatchment::GetMaxBas() const
{
    return maxBas;
}
void  SubCatchment::SetCorrection(double value)
{
    correction = value;
}
double  SubCatchment::GetCorrection() const
{
    return correction;
}
void  SubCatchment::SetManning(double value)
{
    manning = value;
}
double  SubCatchment::GetManning() const
{
    return manning;
}
double  SubCatchment::GetRiverSlope() const
{
    return riverSlope;
}
void  SubCatchment::SetWidth(double value)
{
    width = value;
}
double  SubCatchment::GetWidth() const
{
    return width;
}
void  SubCatchment::SetAccumulatedRouteInput(double value)
{
    accumulatedRouteInput = value;
};
double  SubCatchment::GetAccumulatedRouteInput() const
{
    return accumulatedRouteInput;
}
void  SubCatchment::SetAccumulatedDischarge(int index, double value)
{
    accumulatedDischarge[index] = value;
}
double  SubCatchment::GetAccumulatedDischarge(int index) const
{
    return accumulatedDischarge[index];
}
void  SubCatchment::SetAccumulatedUpperDischarge(int index, double value)
{
    accumulatedUpperDischarge[index] = value;
}
double  SubCatchment::GetAccumulatedUpperDischarge(int index) const
{
    return accumulatedUpperDischarge[index];
}
void  SubCatchment::SetAccumulatedLowerDischarge(int index, double value)
{
    accumulatedLowerDischarge[index] = value;
}
double  SubCatchment::GetAccumulatedLowerDischarge(int index) const
{
    return accumulatedLowerDischarge[index];
}
void  SubCatchment::SetAccumulatedInFlow(int index, double value)
{
    accumulatedInFlow[index] = value;
}
double  SubCatchment::GetAccumulatedInFlow(int index) const
{
    return accumulatedInFlow[index];
}
//  void SetAccumulatedUpperInFlow(int index, double value) { accumulatedUpperInFlow[index] = value; }
//  double GetAccumulatedUpperInFlow(int index) const { return accumulatedUpperInFlow[index]; }
//  void SetAccumulatedLowerInFlow(int index, double value) { accumulatedLowerInFlow[index] = value; }
//  double GetAccumulatedLowerInFlow(int index) const { return accumulatedLowerInFlow[index]; }
void  SubCatchment::SetAccumulatedPrecipitation(int index, double value)
{
    accumulatedPrecipitation[index] = value;
}
double  SubCatchment::GetAccumulatedPrecipitation(int index) const
{
    return accumulatedPrecipitation[index];
}
void  SubCatchment::SetAccumulatedTemperature(int index, double value)
{
    accumulatedTemperature[index] = value;
}
double  SubCatchment::GetAccumulatedTemperature(int index) const
{
    return accumulatedTemperature[index];
}
void  SubCatchment::SetAccumulatedLakeStorage(int index, double value)
{
    accumulatedLakeStorage[index] = value;
}
double  SubCatchment::GetAccumulatedLakeStorage(int index) const
{
    return accumulatedLakeStorage[index];
}
void  SubCatchment::SetAccumulatedSnowStore(int index, double value)
{
    accumulatedSnowStore[index] = value;
}
double  SubCatchment::GetAccumulatedSnowStore(int index) const
{
    return accumulatedSnowStore[index];
}
void  SubCatchment::SetAccumulatedSnowCoverFraction(int index, double value)
{
    accumulatedSnowCoverFraction[index] = value;
}
double  SubCatchment::GetAccumulatedSnowCoverFraction(int index) const
{
    return accumulatedSnowCoverFraction[index];
}
void  SubCatchment::SetAccumulatedMeltWater(int index, double value)
{
    accumulatedMeltWater[index] = value;
}
double  SubCatchment::GetAccumulatedMeltWater(int index) const
{
    return accumulatedMeltWater[index];
}
void  SubCatchment::SetAccumulatedSnowWaterEquivalentChange(int index, double value)
{
    accumulatedSnowWaterEquivalentChange[index] = value;
}
double  SubCatchment::GetAccumulatedSnowWaterEquivalentChange(int index) const
{
    return accumulatedSnowWaterEquivalentChange[index];
}
void  SubCatchment::SetAccumulatedWaterOutput(int index, double value)
{
    accumulatedWaterOutput[index] = value;
}
double  SubCatchment::GetAccumulatedWaterOutput(int index) const
{
    return accumulatedWaterOutput[index];
}
void  SubCatchment::SetAccumulatedGlacierSnowStore(int index, double value)
{
    accumulatedGlacierSnowStore[index] = value;
}
double  SubCatchment::GetAccumulatedGlacierSnowStore(int index) const
{
    return accumulatedGlacierSnowStore[index];
}
void  SubCatchment::SetAccumulatedGlacierSnowMeltWater(int index, double value)
{
    accumulatedGlacierSnowMeltWater[index] = value;
}
double  SubCatchment::GetAccumulatedGlacierSnowMeltWater(int index) const
{
    return accumulatedGlacierSnowMeltWater[index];
}
void  SubCatchment::SetAccumulatedGlacierMassBalance(int index, double value)
{
    accumulatedGlacierMassBalance[index] = value;
}
double  SubCatchment::GetAccumulatedGlacierMassBalance(int index) const
{
    return accumulatedGlacierMassBalance[index];
}
void  SubCatchment::SetAccumulatedAreaGlacierMassBalance(int index, double value)
{
    accumulatedAreaGlacierMassBalance[index] = value;
}
double  SubCatchment::GetAccumulatedAreaGlacierMassBalance(int index) const
{
    return accumulatedAreaGlacierMassBalance[index];
}
void  SubCatchment::SetAccumulatedGlacierIceMelt(int index, double value)
{
    accumulatedGlacierIceMelt[index] = value;
}
double  SubCatchment::GetAccumulatedGlacierIceMelt(int index) const
{
    return accumulatedGlacierIceMelt[index];
}
void  SubCatchment::SetAccumulatedAnnualMassBalance(int index, double value)
{
    accumulatedAnnualMassBalance[index] = value;
}
double  SubCatchment::GetAccumulatedAnnualMassBalance(int index) const
{
    return accumulatedAnnualMassBalance[index];
}
void  SubCatchment::SetAccumulatedGlacierIceVolume(int index, double value)
{
    accumulatedGlacierIceVolume[index] = value;
}
double  SubCatchment::GetAccumulatedGlacierIceVolume(int index) const
{
    return accumulatedGlacierIceVolume[index];
}
void  SubCatchment::SetAccumulatedEvapotranspiration(int index, double value)
{
    accumulatedEvapotranspiration[index] = value;
}
double  SubCatchment::GetAccumulatedEvapotranspiration(int index) const
{
    return accumulatedEvapotranspiration[index];
}
void  SubCatchment::SetAccumulatedRunoff(int index, double value)
{
    accumulatedRunoff[index] = value;
}
double  SubCatchment::GetAccumulatedRunoff(int index) const
{
    return accumulatedRunoff[index];
}
void  SubCatchment::SetAccumulatedUpperRunoff(int index, double value)
{
    accumulatedUpperRunoff[index] = value;
}
double  SubCatchment::GetAccumulatedUpperRunoff(int index) const
{
    return accumulatedUpperRunoff[index];
}
void  SubCatchment::SetAccumulatedLowerRunoff(int index, double value)
{
    accumulatedLowerRunoff[index] = value;
}
double  SubCatchment::GetAccumulatedLowerRunoff(int index) const
{
    return accumulatedLowerRunoff[index];
}
void  SubCatchment::SetAccumulatedHbvSoilMoisture(int index, double value)
{
    accumulatedHbvSoilMoisture[index] = value;
}
double  SubCatchment::GetAccumulatedHbvSoilMoisture(int index) const
{
    return accumulatedHbvSoilMoisture[index];
}
void  SubCatchment::SetAccumulatedHbvSoilMoistureDeficit(int index, double value)
{
    accumulatedHbvSoilMoistureDeficit[index] = value;
}
double  SubCatchment::GetAccumulatedHbvSoilMoistureDeficit(int index) const
{
    return accumulatedHbvSoilMoistureDeficit[index];
}
void  SubCatchment::SetAccumulatedHbvPercSoilUpper(int index, double value)
{
    accumulatedHbvPercSoilUpper[index] = value;
}
double  SubCatchment::GetAccumulatedHbvPercSoilUpper(int index) const
{
    return accumulatedHbvPercSoilUpper[index];
}
void  SubCatchment::SetAccumulatedHbvUpperZone(int index, double value)
{
    accumulatedHbvUpperZone[index] = value;
}
double  SubCatchment::GetAccumulatedHbvUpperZone(int index) const
{
    return accumulatedHbvUpperZone[index];
}
void  SubCatchment::SetAccumulatedHbvLowerZone(int index, double value)
{
    accumulatedHbvLowerZone[index] = value;
}
double  SubCatchment::GetAccumulatedHbvLowerZone(int index) const
{
    return accumulatedHbvLowerZone[index];
}
void  SubCatchment::SetAccumulatedSum(int index, double value)
{
    accumulatedSum[index] = value;
}
double  SubCatchment::GetAccumulatedSum(int index) const
{
    return accumulatedSum[index];
}
void  SubCatchment::SetAccumulatedSumLake(int index, double value)
{
    accumulatedSumLake[index] = value;
}
double  SubCatchment::GetAccumulatedSumLake(int index) const
{
    return accumulatedSumLake[index];
}
void  SubCatchment::SetAccumulatedSumSnow(int index, double value)
{
    accumulatedSumSnow[index] = value;
}
double  SubCatchment::GetAccumulatedSumSnow(int index) const
{
    return accumulatedSumSnow[index];
}
void  SubCatchment::SetAccumulatedSumGlacier(int index, double value)
{
    accumulatedSumGlacier[index] = value;
}
double  SubCatchment::GetAccumulatedSumGlacier(int index) const
{
    return accumulatedSumGlacier[index];
}
void  SubCatchment::SetAccumulatedSumHbv(int index, double value)
{
    accumulatedSumHbv[index] = value;
}
double  SubCatchment::GetAccumulatedSumHbv(int index) const
{
    return accumulatedSumHbv[index];
}
void  SubCatchment::SetSubCatchmentPrecipitation(int index, double value)
{
    subCatchmentPrecipitation[index] = value;
}
double  SubCatchment::GetSubCatchmentPrecipitation(int index) const
{
    return subCatchmentPrecipitation[index];
}
void  SubCatchment::SetSubCatchmentTemperature(int index, double value)
{
    subCatchmentTemperature[index] = value;
}
double  SubCatchment::GetSubCatchmentTemperature(int index) const
{
    return subCatchmentTemperature[index];
}
void  SubCatchment::SetSubCatchmentLakeStorage(int index, double value)
{
    subCatchmentLakeStorage[index] = value;
}
double  SubCatchment::GetSubCatchmentLakeStorage(int index) const
{
    return subCatchmentLakeStorage[index];
}
void  SubCatchment::SetSubCatchmentSnowStore(int index, double value)
{
    subCatchmentSnowStore[index] = value;
}
double  SubCatchment::GetSubCatchmentSnowStore(int index) const
{
    return subCatchmentSnowStore[index];
}
void  SubCatchment::SetSubCatchmentSnowCoverFraction(int index, double value)
{
    subCatchmentSnowCoverFraction[index] = value;
}
double  SubCatchment::GetSubCatchmentSnowCoverFraction(int index) const
{
    return subCatchmentSnowCoverFraction[index];
}
void  SubCatchment::SetSubCatchmentMeltWater(int index, double value)
{
    subCatchmentMeltWater[index] = value;
}
double  SubCatchment::GetSubCatchmentMeltWater(int index) const
{
    return subCatchmentMeltWater[index];
}
void  SubCatchment::SetSubCatchmentSnowWaterEquivalentChange(int index, double value)
{
    subCatchmentSnowWaterEquivalentChange[index] = value;
}
double  SubCatchment::GetSubCatchmentSnowWaterEquivalentChange(int index) const
{
    return subCatchmentSnowWaterEquivalentChange[index];
}
void  SubCatchment::SetSubCatchmentWaterOutput(int index, double value)
{
    subCatchmentWaterOutput[index] = value;
}
double  SubCatchment::GetSubCatchmentWaterOutput(int index) const
{
    return subCatchmentWaterOutput[index];
}
void  SubCatchment::SetSubCatchmentGlacierSnowStore(int index, double value)
{
    subCatchmentGlacierSnowStore[index] = value;
}
double  SubCatchment::GetSubCatchmentGlacierSnowStore(int index) const
{
    return subCatchmentGlacierSnowStore[index];
}
void  SubCatchment::SetSubCatchmentGlacierSnowMeltWater(int index, double value)
{
    subCatchmentGlacierSnowMeltWater[index] = value;
}
double  SubCatchment::GetSubCatchmentGlacierSnowMeltWater(int index) const
{
    return subCatchmentGlacierSnowMeltWater[index];
}
void  SubCatchment::SetSubCatchmentGlacierMassBalance(int index, double value)
{
    subCatchmentGlacierMassBalance[index] = value;
}
double  SubCatchment::GetSubCatchmentGlacierMassBalance(int index) const
{
    return subCatchmentGlacierMassBalance[index];
}
void  SubCatchment::SetSubCatchmentAreaGlacier(int index, double value)
{
    subCatchmentAreaGlacier[index] = value;
}
double  SubCatchment::GetSubCatchmentAreaGlacier(int index) const
{
    return subCatchmentAreaGlacier[index];
}
void  SubCatchment::SetSubCatchmentAreaGlacierMassBalance(int index, double value)
{
    subCatchmentAreaGlacierMassBalance[index] = value;
}
double  SubCatchment::GetSubCatchmentAreaGlacierMassBalance(int index) const
{
    return subCatchmentAreaGlacierMassBalance[index];
}
void  SubCatchment::SetSubCatchmentGlacierIceMelt(int index, double value)
{
    subCatchmentGlacierIceMelt[index] = value;
}
double  SubCatchment::GetSubCatchmentGlacierIceMelt(int index) const
{
    return subCatchmentGlacierIceMelt[index];
}
void  SubCatchment::SetSubCatchmentAnnualMassBalance(int index, double value)
{
    subCatchmentAnnualMassBalance[index] = value;
}
double  SubCatchment::GetSubCatchmentAnnualMassBalance(int index) const
{
    return subCatchmentAnnualMassBalance[index];
}
void  SubCatchment::SetSubCatchmentGlacierIceVolume(int index, double value)
{
    subCatchmentGlacierIceVolume[index] = value;
}
double  SubCatchment::GetSubCatchmentGlacierIceVolume(int index) const
{
    return subCatchmentGlacierIceVolume[index];
}
void  SubCatchment::SetSubCatchmentEvapotranspiration(int index, double value)
{
    subCatchmentEvapotranspiration[index] = value;
}
double  SubCatchment::GetSubCatchmentEvapotranspiration(int index) const
{
    return subCatchmentEvapotranspiration[index];
}
void  SubCatchment::SetSubCatchmentRunoff(int index, double value)
{
    subCatchmentRunoff[index] = value;
}
double  SubCatchment::GetSubCatchmentRunoff(int index) const
{
    return subCatchmentRunoff[index];
}
void  SubCatchment::SetSubCatchmentUpperRunoff(int index, double value)
{
    subCatchmentUpperRunoff[index] = value;
}
double  SubCatchment::GetSubCatchmentUpperRunoff(int index) const
{
    return subCatchmentUpperRunoff[index];
}
void  SubCatchment::SetSubCatchmentLowerRunoff(int index, double value)
{
    subCatchmentLowerRunoff[index] = value;
}
double  SubCatchment::GetSubCatchmentLowerRunoff(int index) const
{
    return subCatchmentLowerRunoff[index];
}
void  SubCatchment::SetSubCatchmentHbvSoilMoisture(int index, double value)
{
    subCatchmentHbvSoilMoisture[index] = value;
}
double  SubCatchment::GetSubCatchmentHbvSoilMoisture(int index) const
{
    return subCatchmentHbvSoilMoisture[index];
}
void  SubCatchment::SetSubCatchmentHbvSoilMoistureDeficit(int index, double value)
{
    subCatchmentHbvSoilMoistureDeficit[index] = value;
}
double  SubCatchment::GetSubCatchmentHbvSoilMoistureDeficit(int index) const
{
    return subCatchmentHbvSoilMoistureDeficit[index];
}
void  SubCatchment::SetSubCatchmentHbvPercSoilUpper(int index, double value)
{
    subCatchmentHbvPercSoilUpper[index] = value;
}
double  SubCatchment::GetSubCatchmentHbvPercSoilUpper(int index) const
{
    return subCatchmentHbvPercSoilUpper[index];
}
void  SubCatchment::SetSubCatchmentHbvUpperZone(int index, double value)
{
    subCatchmentHbvUpperZone[index] = value;
}
double  SubCatchment::GetSubCatchmentHbvUpperZone(int index) const
{
    return subCatchmentHbvUpperZone[index];
}
void  SubCatchment::SetSubCatchmentHbvLowerZone(int index, double value)
{
    subCatchmentHbvLowerZone[index] = value;
}
double  SubCatchment::GetSubCatchmentHbvLowerZone(int index) const
{
    return subCatchmentHbvLowerZone[index];
}
double  SubCatchment::GetObsData(int index) const
{
    return observedData[index];
}

void SubCatchment::SetNumUpStream(int value)
{
    numUpStream = value;
    SubCatchment **upStr = new SubCatchment *[value];
    UpStream = upStr;
}

void SubCatchment::SetRiverSlope(double value)
{
    if (value < 0.01)
    {
        value = 0.01;
    }
    riverSlope = tan(value * acos(-1.0) / 180.0);
}

void SubCatchment::AllocateAccumulatedDischarge(int numberTimeSteps)
{
    accumulatedDischarge = new double[numberTimeSteps];
    accumulatedUpperDischarge = new double[numberTimeSteps];
    accumulatedLowerDischarge = new double[numberTimeSteps];
}

void SubCatchment::AllocateAccumulatedInFlow(int numberTimeSteps)
{
    accumulatedInFlow = new double[numberTimeSteps];
    //  accumulatedUpperInFlow = new double [numberTimeSteps];
    //  accumulatedLowerInFlow = new double [numberTimeSteps];
}

void SubCatchment::AllocateAccumulatedWaterBalance(int numberTimeSteps)
{
    accumulatedPrecipitation = new double[numberTimeSteps];
    accumulatedTemperature = new double[numberTimeSteps];
    accumulatedLakeStorage = new double[numberTimeSteps];
    accumulatedSnowStore = new double[numberTimeSteps];
    accumulatedSnowCoverFraction = new double[numberTimeSteps];
    accumulatedMeltWater = new double[numberTimeSteps];
    accumulatedSnowWaterEquivalentChange = new double[numberTimeSteps];
    accumulatedWaterOutput = new double[numberTimeSteps];
    accumulatedGlacierSnowStore = new double [numberTimeSteps];
    accumulatedGlacierSnowMeltWater = new double [numberTimeSteps];
    accumulatedGlacierMassBalance = new double[numberTimeSteps];
    accumulatedAreaGlacierMassBalance = new double[numberTimeSteps];
    accumulatedGlacierIceMelt = new double[numberTimeSteps];
    accumulatedAnnualMassBalance = new double[numberTimeSteps];
    accumulatedGlacierIceVolume = new double[numberTimeSteps];
    accumulatedEvapotranspiration = new double[numberTimeSteps];
    accumulatedRunoff = new double[numberTimeSteps];
    accumulatedUpperRunoff = new double[numberTimeSteps];
    accumulatedLowerRunoff = new double[numberTimeSteps];
    accumulatedHbvSoilMoisture = new double[numberTimeSteps];
    accumulatedHbvSoilMoistureDeficit = new double[numberTimeSteps];
    accumulatedHbvPercSoilUpper = new double[numberTimeSteps];
    accumulatedHbvUpperZone = new double[numberTimeSteps];
    accumulatedHbvLowerZone = new double[numberTimeSteps];
    accumulatedSum = new double[numberTimeSteps];
    accumulatedSumLake = new double[numberTimeSteps];
    accumulatedSumSnow = new double[numberTimeSteps];
    accumulatedSumGlacier = new double[numberTimeSteps];
    accumulatedSumHbv = new double[numberTimeSteps];
}

void SubCatchment::AllocateWaterBalance(int numberTimeSteps)
{
    subCatchmentPrecipitation = new double[numberTimeSteps];
    subCatchmentTemperature = new double[numberTimeSteps];
    subCatchmentLakeStorage = new double[numberTimeSteps];
    subCatchmentSnowStore = new double[numberTimeSteps];
    subCatchmentSnowCoverFraction = new double[numberTimeSteps];
    subCatchmentMeltWater = new double[numberTimeSteps];
    subCatchmentSnowWaterEquivalentChange = new double[numberTimeSteps];
    subCatchmentWaterOutput = new double[numberTimeSteps];
    subCatchmentGlacierSnowStore = new double [numberTimeSteps];
    subCatchmentGlacierSnowMeltWater = new double [numberTimeSteps];
    subCatchmentGlacierMassBalance = new double[numberTimeSteps];
    subCatchmentAreaGlacier = new double[numberTimeSteps];
    subCatchmentAreaGlacierMassBalance = new double[numberTimeSteps];
    subCatchmentGlacierIceMelt = new double[numberTimeSteps];
    subCatchmentAnnualMassBalance = new double[numberTimeSteps];
    subCatchmentGlacierIceVolume = new double[numberTimeSteps];
    subCatchmentEvapotranspiration = new double[numberTimeSteps];
    subCatchmentRunoff = new double[numberTimeSteps];
    subCatchmentUpperRunoff = new double[numberTimeSteps];
    subCatchmentLowerRunoff = new double[numberTimeSteps];
    subCatchmentHbvSoilMoisture = new double[numberTimeSteps];
    subCatchmentHbvSoilMoistureDeficit = new double[numberTimeSteps];
    subCatchmentHbvPercSoilUpper = new double[numberTimeSteps];
    subCatchmentHbvUpperZone = new double[numberTimeSteps];
    subCatchmentHbvLowerZone = new double[numberTimeSteps];
}

void SubCatchment::ObsDataInput(char * fileObsStreamflow, DateTime firstTime, DateTime lastTime, int numberTimeSteps, int secondsPerTimeStep)
{
    char c, buffer[1024];
    int i, subCatchmentId, year, month, day, hour, min;
    double value;
    bool catchmentFound = false;
    DateTime * datetime = new DateTime[numberTimeSteps];
    datetime[0] = firstTime;
    for (i = 1; i < numberTimeSteps; i++) 
    {
        datetime[i] = datetime[i-1] + secondsPerTimeStep; 
    }
    /*    for (i = 0; i < numberTimeSteps; i++)
	  {
	  datetime[i] = firstTime + i * secondsPerTimeStep;
	  }*/
    if (datetime[0] != firstTime || datetime[numberTimeSteps - 1] != lastTime)
    {
        cout << endl << " DateTime error during initialisation of observed data array: " <<
             datetime[0] << "  " << firstTime << "  " << datetime[numberTimeSteps - 1] << "  " << lastTime << endl << endl;
        exit(1);
    }
    observedData = new double[numberTimeSteps];
    for (i = 0; i < numberTimeSteps; i++)
    {
        observedData[i] = missingData;
    }
    ifstream fileIn(fileObsStreamflow);
    if (!fileIn.is_open())
    {
        cout << endl << " Error opening file " << fileObsStreamflow << endl << endl;
        //    exit (1);
    }
    while (fileIn.getline(buffer, 1024) && !catchmentFound)
    {
        if (buffer[0] == '#')
        {
            sscanf(buffer, "%c %c %c %c %c %d ", &c, &c, &c, &c, &c, &subCatchmentId);
            if (subCatchmentId == GetIdentifier())
            {
                i = 0;
                while (fileIn.peek() != '#' && fileIn.peek() != EOF)
                {
                    fileIn >> buffer >> value;
                    sscanf(buffer, "%4d %2d %2d %c %2d %2d ", &year, &month, &day, &c, &hour, &min);
                    while ((datetime[i].getYear() < year ||
                            (datetime[i].getYear() == year && datetime[i].getMonth() < month) ||
                            (datetime[i].getYear() == year && datetime[i].getMonth() == month && datetime[i].getDay() < day) ||
                            (datetime[i].getYear() == year && datetime[i].getMonth() == month && datetime[i].getDay() == day &&
                             datetime[i].getHour() < hour) ||
                            (datetime[i].getYear() == year && datetime[i].getMonth() == month && datetime[i].getDay() == day &&
                             datetime[i].getHour() == hour && datetime[i].getMinute() < min)) &&
                            i < numberTimeSteps - 1)
                    {
                        i++;
                        //      cout << " i = " << i << endl;
                    }
                    //    cout << datetime[i].getYear() << " " << datetime[i].getMonth() << " " << datetime[i].getDay() << " " <<
                    //      datetime[i].getHour() << " " << datetime[i].getMinute() << endl;
                    if (datetime[i].getYear() == year && datetime[i].getMonth() == month && datetime[i].getDay() == day &&
                            datetime[i].getHour() == hour && datetime[i].getMinute() == min)
                    {
                        //      cout << " data found  " << i << " " << year << " " << month << " " << day << " " << hour << " " << min << " " << endl;
                        if (value <= -9.0)
                        {
                            observedData[i] = missingData;
                        }
                        else if (observedData[i] > -9.0 && observedData[i] < 0.0)
                        {
                            observedData[i] = 0.0;
                        }
                        else
                        {
                            observedData[i] = value;
                        }
                        //    cout << subCatchmentId << " " << i << " " << year << " " << month << " ";
                        //    cout << day << " " << hour << " " << min << " " << " " << value << endl;
                    }
                    fileIn.ignore(1024, '\n');
                }
            }
            else
            {
                while (fileIn.peek() != '#' && fileIn.peek() != EOF)
                {
                    //      fileIn.getline(buffer, 1024);
                    fileIn.ignore(1024, '\n');
                }
            }
        }
    }
    fileIn.close();
}
