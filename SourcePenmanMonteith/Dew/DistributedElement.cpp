#include "DistributedElement.h"
#include "stdafx.h"
#include "DateTime.h"
#include "Glacier.h"
#include "HbvAquifer.h"
#include "KwaAquifer.h"
#include "ParametersGeneral.h"
#include "TotalReservoirStorage.h"
#include "Lake.h"

DistributedElement::DistributedElement() :
    xCoord(0.0),
    yCoord(0.0),
    lakeNumber(-9999),
    numUpLand(0),
    flowDirection(0),
    accumulatedSum(0.0),
    accumulatedSumLake(0.0),
    accumulatedSumSnow(0.0),
    accumulatedSumGlacier(0.0),
    accumulatedSumHbv(0.0),
    latitude(0.0),
    area(0.0),
    elevation(0.0),
    timefalling(0.0),
    slopeLength(0.0),
    slopeAngle(0.0),
    aspect(0.0),
    travelTime(-9999.0),
    Pcorr(0.0),
    precipitationCorrection(1.0),
    temperatureCorrection(0.0),
    accumulatedDischarge(0.0),
    accumulatedLowerDischarge(0.0),
    accumulatedUpperDischarge(0.0),
    accumulatedPrecipitation(0.0),
    accumulatedTemperature(0.0),
    accumulatedLakeStorage(0.0),
    accumulatedSnowStore(0.0),
    accumulatedSnowCoverFraction(0.0),
    accumulatedMeltWater(0.0),
    accumulatedSnowWaterEquivalentChange(0.0),
    accumulatedWaterOutput(0.0),
    accumulatedGlacierSnowStore(0.0),
    accumulatedGlacierSnowMeltWater(0.0),
    accumulatedGlacierMassBalance(0.0),
    accumulatedAreaGlacierMassBalance(0.0),
    accumulatedGlacierIceMelt(0.0),
    accumulatedAnnualMassBalance(0.0),
    accumulatedGlacierIceVolume(0.0),
    accumulatedEvapotranspiration(0.0),
    accumulatedRunoff(0.0),
    accumulatedUpperRunoff(0.0),
    accumulatedLowerRunoff(0.0),
    accumulatedHbvSoilMoisture(0.0),
    accumulatedHbvSoilMoistureDeficit(0.0),
    accumulatedHbvPercSoilUpper(0.0),
    accumulatedHbvUpperZone(0.0),
    accumulatedHbvLowerZone(0.0),
    precStationsWeightedElevation(0.0),
    tempStationsWeightedElevation(0.0),
    numberSum(0),
    initialStorage(0),
    finalStorage(0),
    sumPrecipitation(0.0),
    sumEvapotranspiration(0.0),
    sumRunoff(0.0),
    sumGlacierIceMelt(0.0)
{
    upLandFlow = 0;
    SetSubCatchmentElement(0);
    SetNextElement(0);
    SetNextGlacierElement(0);
    SetLake(0);
    SetGlacier(0);
    SetHbvAquifer(0);
    SetKwaAquifer(0);
    SetSelectedKiWaTimeSeriesElements(0);
    SetSelectedHbvTimeSeriesElements(0);
    SetGeneralPar(0);
    SetTotalReservoir(0);
    SetModelControlObj(0);
}

DistributedElement::~DistributedElement()
{
}

void DistributedElement::SetGeoIndex(int value)
{
    geoIndex = value;
}
int DistributedElement::GetGeoIndex() const
{
    return geoIndex;
}
void DistributedElement::SetLandIndex(int value)
{
    landIndex = value;
}
int DistributedElement::GetLandIndex() const
{
    return landIndex;
}
void DistributedElement::SetXCoord(double value)
{
    xCoord = value;
}
double DistributedElement::GetXCoord() const
{
    return xCoord;
}
void DistributedElement::SetYCoord(double value)
{
    yCoord = value;
}
double DistributedElement::GetYCoord() const
{
    return yCoord;
}
void DistributedElement::SetLakeNumber(int value)
{
    lakeNumber = value;
}
int DistributedElement::GetLakeNumber() const
{
    return lakeNumber;
}
void AllocateInputValues(int value);
int DistributedElement::GetNumberInputValues() const
{
    return numberInputValues;
}
void DistributedElement::SetInputValue(int i, double value)
{
    inputArray[i] = value;
}
double DistributedElement::GetInputValue(int i) const
{
    return inputArray[i];
}
void DistributedElement::SetArea(double value)
{
    area = value;
}
double DistributedElement::GetArea() const
{
    return area;
}
void DistributedElement::SetElevation(double value)
{
    elevation = value;
}
double DistributedElement::GetElevation() const
{
    return elevation;
}
void DistributedElement::SetSlopeLength(double value)
{
    slopeLength = value;
}
double DistributedElement::GetSlopeLength() const
{
    return slopeLength;
}
double DistributedElement::GetSlopeAngle() const
{
    return slopeAngle;
}
void DistributedElement::SetAspect(double value)
{
    aspect = value;
}
double DistributedElement::GetAspect() const
{
    return aspect;
}
void DistributedElement::SetFlowDirection(int value)
{
    flowDirection = value;
}
int DistributedElement::GetFlowDirection() const
{
    return flowDirection;
}
void DistributedElement::SetTravelTime(double value)
{
    travelTime = value;
}
double DistributedElement::GetTravelTime() const
{
    return travelTime;
}
void DistributedElement::SetSelectedKiWaTimeSeriesElements(SelectedKiWaTimeSeriesElements *object)
{
    selectedKiWaElements = object;
}
SelectedKiWaTimeSeriesElements *DistributedElement::GetSelectedKiWaTimeSeriesElements() const
{
    return selectedKiWaElements;
}
void DistributedElement::SetSelectedHbvTimeSeriesElements(SelectedHbvTimeSeriesElements *object)
{
    selectedHbvElements = object;
}
SelectedHbvTimeSeriesElements *DistributedElement::GetSelectedHbvTimeSeriesElements() const
{
    return selectedHbvElements;
}
void DistributedElement::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral *DistributedElement::GetGeneralPar() const
{
    return commonPar;
}
void DistributedElement::SetTotalReservoir(TotalReservoirStorage *resObj)
{
    totalReservoir = resObj;
}
TotalReservoirStorage *DistributedElement::GetTotalReservoir() const
{
    return totalReservoir;
}
void DistributedElement::SetModelControlObj(ModelControl *controlObj)
{
    modelControlObject = controlObj;
}
ModelControl *DistributedElement::GetModelControlObj() const
{
    return modelControlObject;
}
void DistributedElement::SetPrecipitationCorrection(double value)
{
    precipitationCorrection = value;
}
double DistributedElement::GetPrecipitationCorrection() const
{
    return precipitationCorrection;
}
void DistributedElement::SetTemperatureCorrection(double value)
{
    temperatureCorrection = value;
}
double DistributedElement::GetTemperatureCorrection() const
{
    return temperatureCorrection;
}
int DistributedElement::GetNumUpLand() const
{
    return numUpLand;
}
void DistributedElement::SetUpLandFlow(int k, DistributedElement *theElement)
{
    upLandFlow[k] = theElement;
}
DistributedElement *DistributedElement::GetUpLandFlow(int k) const
{
    return upLandFlow[k];
}
void DistributedElement::SetNextElement(DistributedElement *theElement)
{
    nextElement = theElement;
}
DistributedElement *DistributedElement::GetNextElement() const
{
    return nextElement;
}
void DistributedElement::SetNextGlacierElement(DistributedElement *theElement)
{
    nextGlacierElement = theElement;
}
DistributedElement *DistributedElement::GetNextGlacierElement() const
{
    return nextGlacierElement;
}
void DistributedElement::SetSubCatchmentElement(SubCatchment *theElement)
{
    subCatchmentElement = theElement;
}
SubCatchment *DistributedElement::GetSubCatchmentElement() const
{
    return subCatchmentElement;
}
void DistributedElement::SetLake(Lake *theLake)
{
    lake = theLake;
}
Lake *DistributedElement::GetLake() const
{
    return lake;
}
void DistributedElement::SetGlacier(Glacier *theGlacier)
{
    glacier = theGlacier;
}
Glacier *DistributedElement::GetGlacier() const
{
    return glacier;
}
void DistributedElement::SetHbvAquifer(HbvAquifer *theHbvAquifer)
{
    hbvAquifer = theHbvAquifer;
}
HbvAquifer *DistributedElement::GetHbvAquifer() const
{
    return hbvAquifer;
}
void DistributedElement::SetKwaAquifer(KwaAquifer *theKwaAquifer)
{
    kwaAquifer = theKwaAquifer;
}
KwaAquifer *DistributedElement::GetKwaAquifer() const
{
    return kwaAquifer;
}
void DistributedElement::SetPrecStationsWeightedElevation(double value)
{
    precStationsWeightedElevation = value;
}
double DistributedElement::GetPrecStationsWeightedElevation()
{
    return precStationsWeightedElevation;
}
void DistributedElement::SetTempStationsWeightedElevation(double value)
{
    tempStationsWeightedElevation = value;
}
double DistributedElement::GetTempStationsWeightedElevation()
{
    return tempStationsWeightedElevation;
}
void AllocateMetSeries(int numberPrecSeries, int numberTempSeries);
void DistributedElement::SetMetSeriesNumber(int k, int number)
{
    metSeriesNumber[k] = number;
}
void DistributedElement::SetMetSeriesWeight(int k, double weight)
{
    metSeriesWeight[k] = weight;
}
int DistributedElement::GetMetSeriesNumber(int k)
{
    return metSeriesNumber[k];
}
double DistributedElement::GetMetSeriesWeight(int k)
{
    return metSeriesWeight[k];
}
void AllocateKiWaWaterBalance(int numberTimeSteps);
void AllocateHbvWaterBalance(int numberTimeSteps);
void DistributedElement::SetDistributedElementPrecipitation(int index, double value)
{
    distributedElementPrecipitation[index] = value;
}
double DistributedElement::GetDistributedElementPrecipitation(int index) const
{
    return distributedElementPrecipitation[index];
}
void DistributedElement::SetDistributedElementTemperature(int index, double value)
{
    distributedElementTemperature[index] = value;
}
double DistributedElement::GetDistributedElementTemperature(int index) const
{
    return distributedElementTemperature[index];
}
void DistributedElement::SetDistributedElementLakeStorage(int index, double value)
{
    distributedElementLakeStorage[index] = value;
}
double DistributedElement::GetDistributedElementLakeStorage(int index) const
{
    return distributedElementLakeStorage[index];
}
void DistributedElement::SetDistributedElementSnowStore(int index, double value)
{
    distributedElementSnowStore[index] = value;
}
double DistributedElement::GetDistributedElementSnowStore(int index) const
{
    return distributedElementSnowStore[index];
}
void DistributedElement::SetDistributedElementSnowCoverFraction(int index, double value)
{
    distributedElementSnowCoverFraction[index] = value;
}
double DistributedElement::GetDistributedElementSnowCoverFraction(int index) const
{
    return distributedElementSnowCoverFraction[index];
}
void DistributedElement::SetDistributedElementMeltWater(int index, double value)
{
    distributedElementMeltWater[index] = value;
}
double DistributedElement::GetDistributedElementMeltWater(int index) const
{
    return distributedElementMeltWater[index];
}
void DistributedElement::SetDistributedElementSnowWaterEquivalentChange(int index, double value)
{
    distributedElementSnowWaterEquivalentChange[index] = value;
}
double DistributedElement::GetDistributedElementSnowWaterEquivalentChange(int index) const
{
    return distributedElementSnowWaterEquivalentChange[index];
}
void DistributedElement::SetDistributedElementWaterOutput(int index, double value)
{
    distributedElementWaterOutput[index] = value;
}
double DistributedElement::GetDistributedElementWaterOutput(int index) const
{
    return distributedElementWaterOutput[index];
}
void DistributedElement::SetDistributedElementGlacierSnowStore(int index, double value)
{
    distributedElementGlacierSnowStore[index] = value;
}
double DistributedElement::GetDistributedElementGlacierSnowStore(int index) const
{
    return distributedElementGlacierSnowStore[index];
}
void DistributedElement::SetDistributedElementGlacierSnowMeltWater(int index, double value)
{
    distributedElementGlacierSnowMeltWater[index] = value;
}
double DistributedElement::GetDistributedElementGlacierSnowMeltWater(int index) const
{
    return distributedElementGlacierSnowMeltWater[index];
}
void DistributedElement::SetDistributedElementGlacierMassBalance(int index, double value)
{
    distributedElementGlacierMassBalance[index] = value;
}
double DistributedElement::GetDistributedElementGlacierMassBalance(int index) const
{
    return distributedElementGlacierMassBalance[index];
}
void DistributedElement::SetDistributedElementAreaGlacierMassBalance(int index, double value)
{
    distributedElementAreaGlacierMassBalance[index] = value;
}
double DistributedElement::GetDistributedElementAreaGlacierMassBalance(int index) const
{
    return distributedElementAreaGlacierMassBalance[index];
}
void DistributedElement::SetDistributedElementGlacierIceMelt(int index, double value)
{
    distributedElementGlacierIceMelt[index] = value;
}
double DistributedElement::GetDistributedElementGlacierIceMelt(int index) const
{
    return distributedElementGlacierIceMelt[index];
}
void DistributedElement::SetDistributedElementAnnualMassBalance(int index, double value)
{
    distributedElementAnnualMassBalance[index] = value;
}
double DistributedElement::GetDistributedElementAnnualMassBalance(int index) const
{
    return distributedElementAnnualMassBalance[index];
}
void DistributedElement::SetDistributedElementGlacierIceVolume(int index, double value)
{
    distributedElementGlacierIceVolume[index] = value;
}
double DistributedElement::GetDistributedElementGlacierIceVolume(int index) const
{
    return distributedElementGlacierIceVolume[index];
}
void DistributedElement::SetDistributedElementEvapotranspiration(int index, double value)
{
    distributedElementEvapotranspiration[index] = value;
}
double DistributedElement::GetDistributedElementEvapotranspiration(int index) const
{
    return distributedElementEvapotranspiration[index];
}
void DistributedElement::SetDistributedElementRunoff(int index, double value)
{
    distributedElementRunoff[index] = value;
}
double DistributedElement::GetDistributedElementRunoff(int index) const
{
    return distributedElementRunoff[index];
}
void DistributedElement::SetDistributedElementUpperRunoff(int index, double value)
{
    distributedElementUpperRunoff[index] = value;
}
double DistributedElement::GetDistributedElementUpperRunoff(int index) const
{
    return distributedElementUpperRunoff[index];
}
void DistributedElement::SetDistributedElementLowerRunoff(int index, double value)
{
    distributedElementLowerRunoff[index] = value;
}
double DistributedElement::GetDistributedElementLowerRunoff(int index) const
{
    return distributedElementLowerRunoff[index];
}
void DistributedElement::SetDistributedElementKiWaSoilMoistureOne(int index, double value)
{
    distributedElementKiWaSoilMoistureOne[index] = value;
}
double DistributedElement::GetDistributedElementKiWaSoilMoistureOne(int index) const
{
    return distributedElementKiWaSoilMoistureOne[index];
}
void DistributedElement::SetDistributedElementKiWaSoilMoistureTwo(int index, double value)
{
    distributedElementKiWaSoilMoistureTwo[index] = value;
}
double DistributedElement::GetDistributedElementKiWaSoilMoistureTwo(int index) const
{
    return distributedElementKiWaSoilMoistureTwo[index];
}
void DistributedElement::SetDistributedElementKiWaGroundWaterDepthOne(int index, double value)
{
    distributedElementKiWaGroundWaterDepthOne[index] = value;
}
double DistributedElement::GetDistributedElementKiWaGroundWaterDepthOne(int index) const
{
    return distributedElementKiWaGroundWaterDepthOne[index];
}
void DistributedElement::SetDistributedElementKiWaGroundWaterDepthTwo(int index, double value)
{
    distributedElementKiWaGroundWaterDepthTwo[index] = value;
}
double DistributedElement::GetDistributedElementKiWaGroundWaterDepthTwo(int index) const
{
    return distributedElementKiWaGroundWaterDepthTwo[index];
}
void DistributedElement::SetDistributedElementHbvSoilMoisture(int index, double value)
{
    distributedElementHbvSoilMoisture[index] = value;
}
double DistributedElement::GetDistributedElementHbvSoilMoisture(int index) const
{
    return distributedElementHbvSoilMoisture[index];
}
void DistributedElement::SetDistributedElementHbvSoilMoistureDeficit(int index, double value)
{
    distributedElementHbvSoilMoistureDeficit[index] = value;
}
double DistributedElement::GetDistributedElementHbvSoilMoistureDeficit(int index) const
{
    return distributedElementHbvSoilMoistureDeficit[index];
}
void DistributedElement::SetDistributedElementHbvPercSoilUpper(int index, double value)
{
    distributedElementHbvPercSoilUpper[index] = value;
}
double DistributedElement::GetDistributedElementHbvPercSoilUpper(int index) const
{
    return distributedElementHbvPercSoilUpper[index];
}
void DistributedElement::SetDistributedElementHbvUpperZone(int index, double value)
{
    distributedElementHbvUpperZone[index] = value;
}
double DistributedElement::GetDistributedElementHbvUpperZone(int index) const
{
    return distributedElementHbvUpperZone[index];
}
void DistributedElement::SetDistributedElementHbvLowerZone(int index, double value)
{
    distributedElementHbvLowerZone[index] = value;
}
double DistributedElement::GetDistributedElementHbvLowerZone(int index) const
{
    return distributedElementHbvLowerZone[index];
}
double DistributedElement::GetAccumulatedDischarge() const
{
    return accumulatedDischarge;
}
double DistributedElement::GetAccumulatedLowerDischarge() const
{
    return accumulatedLowerDischarge;
}
double DistributedElement::GetAccumulatedUpperDischarge() const
{
    return accumulatedUpperDischarge;
}
void DistributedElement::SetAccumulatedSum(double value)
{
    accumulatedSum = value;
}
double DistributedElement::GetAccumulatedSum() const
{
    return accumulatedSum;
}
void DistributedElement::SetAccumulatedSumLake(double value)
{
    accumulatedSumLake = value;
}
double DistributedElement::GetAccumulatedSumLake() const
{
    return accumulatedSumLake;
}
void DistributedElement::SetAccumulatedSumSnow(double value)
{
    accumulatedSumSnow = value;
}
double DistributedElement::GetAccumulatedSumSnow() const
{
    return accumulatedSumSnow;
}
void DistributedElement::SetAccumulatedSumGlacier(double value)
{
    accumulatedSumGlacier = value;
}
double DistributedElement::GetAccumulatedSumGlacier() const
{
    return accumulatedSumGlacier;
}
void DistributedElement::SetAccumulatedSumHbv(double value)
{
    accumulatedSumHbv = value;
}
double DistributedElement::GetAccumulatedSumHbv() const
{
    return accumulatedSumHbv;
}
void DistributedElement::SetAccumulatedPrecipitation(double value)
{
    accumulatedPrecipitation = value;
}
double DistributedElement::GetAccumulatedPrecipitation() const
{
    return accumulatedPrecipitation;
}
void DistributedElement::SetAccumulatedTemperature(double value)
{
    accumulatedTemperature = value;
}
double DistributedElement::GetAccumulatedTemperature() const
{
    return accumulatedTemperature;
}
void DistributedElement::SetAccumulatedLakeStorage(double value)
{
    accumulatedLakeStorage = value;
}
double DistributedElement::GetAccumulatedLakeStorage() const
{
    return accumulatedLakeStorage;
}
void DistributedElement::SetAccumulatedSnowStore(double value)
{
    accumulatedSnowStore = value;
}
double DistributedElement::GetAccumulatedSnowStore() const
{
    return accumulatedSnowStore;
}
void DistributedElement::SetAccumulatedSnowCoverFraction(double value)
{
    accumulatedSnowCoverFraction = value;
}
double DistributedElement::GetAccumulatedSnowCoverFraction() const
{
    return accumulatedSnowCoverFraction;
}
void DistributedElement::SetAccumulatedMeltWater(double value)
{
    accumulatedMeltWater = value;
}
double DistributedElement::GetAccumulatedMeltWater() const
{
    return accumulatedMeltWater;
}
void DistributedElement::SetAccumulatedSnowWaterEquivalentChange(double value)
{
    accumulatedSnowWaterEquivalentChange = value;
}
double DistributedElement::GetAccumulatedSnowWaterEquivalentChange() const
{
    return accumulatedSnowWaterEquivalentChange;
}
void DistributedElement::SetAccumulatedWaterOutput(double value)
{
    accumulatedWaterOutput = value;
}
double DistributedElement::GetAccumulatedWaterOutput() const
{
    return accumulatedWaterOutput;
}
void DistributedElement::SetAccumulatedGlacierMassBalance(double value)
{
    accumulatedGlacierMassBalance = value;
}
double DistributedElement::GetAccumulatedGlacierMassBalance() const
{
    return accumulatedGlacierMassBalance;
}
void DistributedElement::SetAccumulatedGlacierSnowStore(double value)
{
    accumulatedGlacierSnowStore = value;
}
double DistributedElement::GetAccumulatedGlacierSnowStore() const
{
    return accumulatedGlacierSnowStore;
}
void DistributedElement::SetAccumulatedGlacierSnowMeltWater(double value)
{
    accumulatedGlacierSnowMeltWater = value;
}
double DistributedElement::GetAccumulatedGlacierSnowMeltWater() const
{
    return accumulatedGlacierSnowMeltWater;
}
void DistributedElement::SetAccumulatedAreaGlacierMassBalance(double value)
{
    accumulatedAreaGlacierMassBalance = value;
}
double DistributedElement::GetAccumulatedAreaGlacierMassBalance() const
{
    return accumulatedAreaGlacierMassBalance;
}
void DistributedElement::SetAccumulatedGlacierIceMelt(double value)
{
    accumulatedGlacierIceMelt = value;
}
double DistributedElement::GetAccumulatedGlacierIceMelt() const
{
    return accumulatedGlacierIceMelt;
}
void DistributedElement::SetAccumulatedAnnualMassBalance(double value)
{
    accumulatedAnnualMassBalance = value;
}
double DistributedElement::GetAccumulatedAnnualMassBalance() const
{
    return accumulatedAnnualMassBalance;
}
void DistributedElement::SetAccumulatedGlacierIceVolume(double value)
{
    accumulatedGlacierIceVolume = value;
}
double DistributedElement::GetAccumulatedGlacierIceVolume() const
{
    return accumulatedGlacierIceVolume;
}
void DistributedElement::SetAccumulatedEvapotranspiration(double value)
{
    accumulatedEvapotranspiration = value;
}
double DistributedElement::GetAccumulatedEvapotranspiration() const
{
    return accumulatedEvapotranspiration;
}
void DistributedElement::SetAccumulatedRunoff(double value)
{
    accumulatedRunoff = value;
}
double DistributedElement::GetAccumulatedRunoff() const
{
    return accumulatedRunoff;
}
void DistributedElement::SetAccumulatedUpperRunoff(double value)
{
    accumulatedUpperRunoff = value;
}
double DistributedElement::GetAccumulatedUpperRunoff() const
{
    return accumulatedUpperRunoff;
}
void DistributedElement::SetAccumulatedLowerRunoff(double value)
{
    accumulatedLowerRunoff = value;
}
double DistributedElement::GetAccumulatedLowerRunoff() const
{
    return accumulatedLowerRunoff;
}
void DistributedElement::SetAccumulatedHbvSoilMoisture(double value)
{
    accumulatedHbvSoilMoisture = value;
}
double DistributedElement::GetAccumulatedHbvSoilMoisture() const
{
    return accumulatedHbvSoilMoisture;
}
void DistributedElement::SetAccumulatedHbvSoilMoistureDeficit(double value)
{
    accumulatedHbvSoilMoistureDeficit = value;
}
double DistributedElement::GetAccumulatedHbvSoilMoistureDeficit() const
{
    return accumulatedHbvSoilMoistureDeficit;
}
void DistributedElement::SetAccumulatedHbvPercSoilUpper(double value)
{
    accumulatedHbvPercSoilUpper = value;
}
double DistributedElement::GetAccumulatedHbvPercSoilUpper() const
{
    return accumulatedHbvPercSoilUpper;
}
void DistributedElement::SetAccumulatedHbvUpperZone(double value)
{
    accumulatedHbvUpperZone = value;
}
double DistributedElement::GetAccumulatedHbvUpperZone() const
{
    return accumulatedHbvUpperZone;
}
void DistributedElement::SetAccumulatedHbvLowerZone(double value)
{
    accumulatedHbvLowerZone = value;
}
double DistributedElement::GetAccumulatedHbvLowerZone() const
{
    return accumulatedHbvLowerZone;
}

void DistributedElement::SetSlopeAngle(double value)
{
    if (value < 0.01)
    {
        value = 0.01;
    }
    slopeAngle = value;
}

void DistributedElement::SetNumUpLand(int value)
{
    numUpLand = value;
    DistributedElement **upLand = new DistributedElement *[value];
    upLandFlow = upLand;
}

void DistributedElement::AllocateInputValues(int value)
{
    int i;
    numberInputValues = value;
    inputArray = new double[numberInputValues];
    for (i = 0; i < numberInputValues; i++)
    {
        inputArray[i] = missingData;
    }
}

void DistributedElement::AllocateMetSeries(int numberPrecSeries, int numberTempSeries)
{
    metSeriesNumber = new int[numberPrecSeries + numberTempSeries];
    metSeriesWeight = new double[numberPrecSeries + numberTempSeries];
}

void DistributedElement::AllocateKiWaWaterBalance(int numberTimeSteps)
{
    distributedElementPrecipitation = new double[numberTimeSteps];
    distributedElementTemperature = new double[numberTimeSteps];
    distributedElementLakeStorage = new double[numberTimeSteps];
    distributedElementSnowStore = new double[numberTimeSteps];
    distributedElementSnowCoverFraction = new double[numberTimeSteps];
    distributedElementMeltWater = new double[numberTimeSteps];
    distributedElementSnowWaterEquivalentChange = new double[numberTimeSteps];
    distributedElementWaterOutput = new double[numberTimeSteps];
    distributedElementGlacierMassBalance = new double[numberTimeSteps];
    distributedElementGlacierSnowStore = new double [numberTimeSteps];
    distributedElementGlacierSnowMeltWater = new double [numberTimeSteps];
    distributedElementAreaGlacierMassBalance = new double[numberTimeSteps];
    distributedElementGlacierIceMelt = new double[numberTimeSteps];
    distributedElementAnnualMassBalance = new double[numberTimeSteps];
    distributedElementGlacierIceVolume = new double[numberTimeSteps];
    distributedElementEvapotranspiration = new double[numberTimeSteps];
    distributedElementRunoff = new double[numberTimeSteps];
    distributedElementUpperRunoff = new double[numberTimeSteps];
    distributedElementLowerRunoff = new double[numberTimeSteps];
    distributedElementKiWaSoilMoistureOne = new double[numberTimeSteps];
    distributedElementKiWaSoilMoistureTwo = new double[numberTimeSteps];
    distributedElementKiWaGroundWaterDepthOne = new double[numberTimeSteps];
    distributedElementKiWaGroundWaterDepthTwo = new double[numberTimeSteps];
}

void DistributedElement::AllocateHbvWaterBalance(int numberTimeSteps)
{
    distributedElementPrecipitation = new double[numberTimeSteps];
    distributedElementTemperature = new double[numberTimeSteps];
    distributedElementLakeStorage = new double[numberTimeSteps];
    distributedElementSnowStore = new double[numberTimeSteps];
    distributedElementSnowCoverFraction = new double[numberTimeSteps];
    distributedElementMeltWater = new double[numberTimeSteps];
    distributedElementSnowWaterEquivalentChange = new double[numberTimeSteps];
    distributedElementWaterOutput = new double[numberTimeSteps];
    distributedElementGlacierSnowStore = new double [numberTimeSteps];
    distributedElementGlacierSnowMeltWater = new double [numberTimeSteps];
    distributedElementGlacierMassBalance = new double[numberTimeSteps];
    distributedElementAreaGlacierMassBalance = new double[numberTimeSteps];
    distributedElementGlacierIceMelt = new double[numberTimeSteps];
    distributedElementAnnualMassBalance = new double[numberTimeSteps];
    distributedElementGlacierIceVolume = new double[numberTimeSteps];
    distributedElementEvapotranspiration = new double[numberTimeSteps];
    distributedElementRunoff = new double[numberTimeSteps];
    distributedElementUpperRunoff = new double[numberTimeSteps];
    distributedElementLowerRunoff = new double[numberTimeSteps];
    distributedElementHbvSoilMoisture = new double[numberTimeSteps];
    distributedElementHbvSoilMoistureDeficit = new double[numberTimeSteps];
    distributedElementHbvPercSoilUpper = new double[numberTimeSteps];
    distributedElementHbvUpperZone = new double[numberTimeSteps];
    distributedElementHbvLowerZone = new double[numberTimeSteps];
}

void DistributedElement::SetAccumulatedDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect)
{
    // Algorithm to be performed in case: no input to landscape element from upstream elements
    //  accumulatedDischarge = upLandValue+localValue;
    // Algorithm to be performed in case: input to landscape element from upstream elements
    //  if (flowHierarchy && GetHbvAquifer())
    if (flowHierarchy && !forceDirect)
    {
        accumulatedDischarge = localValue;
    }
    else
    {
        accumulatedDischarge = upLandValue + localValue;
    }
}

void DistributedElement::SetAccumulatedLowerDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect)
{
    if (flowHierarchy && !forceDirect)
    {
        accumulatedLowerDischarge = localValue;
    }
    else
    {
        accumulatedLowerDischarge = upLandValue + localValue;
    }
}

void DistributedElement::SetAccumulatedUpperDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect)
{
    if (flowHierarchy && !forceDirect)
    {
        accumulatedUpperDischarge = localValue;
    }
    else
    {
        accumulatedUpperDischarge = upLandValue + localValue;
    }
}

void DistributedElement::SetThisYearAnnualGlacierValues(DateTime datetime)
{
    if (GetGlacier())
    {
        GetGlacier()->SetThisYearAnnualGlacierValues(datetime);
    }
}

void DistributedElement::SnowToGlacierIce()
{
    if (GetGlacier())
    {
        GetGlacier()->SnowToGlacierIce();
    }
}

void DistributedElement::RemoveSnowOnGlacierIce()
{
    if (GetGlacier())
    {
        GetGlacier()->RemoveSnowOnGlacierIce();
    }
}

void DistributedElement::SetAnnualMassBalance(double value)
{
    if (GetGlacier())
    {
        GetGlacier()->SetAnnualMassBalance(value);
    }
}

void DistributedElement::SetSnowStore(double value)
{
    //  cout << "\n * SetSnowStore " << GetGeoIndex() << endl;
    if (GetGlacier())
    {
        GetGlacier()->SetSnowStore(value);
    }
    if (GetHbvAquifer())
    {
        GetHbvAquifer()->SetSnowStore(value);
    }
    if (GetKwaAquifer())
    {
        GetKwaAquifer()->SetSnowStore(value);
    }
}

void DistributedElement::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
    if (GetGlacier())
        if (GetGlacier()->GetHBV())
        {
            GetGlacier()->SetSubSurfaceHbvStore(sm, uz, lz);
        }
    if (GetHbvAquifer())
    {
        GetHbvAquifer()->SetSubSurfaceHbvStore(sm, uz, lz);
    }
}

double DistributedElement::GetLandArea() const
{
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetLandArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    if (elementFound)
    {
        return sumAreaFraction * GetArea();
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetLakeArea() const
{
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetLakeArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    if (elementFound)
    {
        return sumAreaFraction * GetArea();
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierArea() const
{
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetGlacierArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    if (elementFound)
    {
        return sumAreaFraction * GetArea();
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierIceArea() const
{
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetGlacierIceArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    if (elementFound)
    {
        return sumAreaFraction * GetArea();
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetKiWaArea() const
{
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetKinematicWave())
        {
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetKwaAquifer())
    {
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetKiWaArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    if (elementFound)
    {
        return sumAreaFraction * GetArea();
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetHbvArea() const
{
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetHBV())
        {
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetHbvAquifer())
    {
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetHbvArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    if (elementFound)
    {
        return sumAreaFraction * GetArea();
    }
    else
    {
        return missingData;
    }
}

// Algorithm to be performed in case: no input to landscape element from upstream elements
//void DistributedElement::WaterBalance(int timeStep) const
// Algorithm to be performed in case: input to landscape element from upstream elements
void DistributedElement::WaterBalance(int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge, int dayofyear) const
{
    HbvAquifer *lastHbvAquifer;
    KwaAquifer *lastKwaAquifer;
    //  cout << "start WaterBalance " << GetLandIndex() << endl;
    if (GetLake())
    {
        GetLake()->WaterBalance(timeStep, datetime, upLandAccumulatedLowerDischarge + upLandAccumulatedUpperDischarge, dayofyear);
    }
    if (GetGlacier())
    {
        GetGlacier()->WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge, dayofyear);
    }
    if (GetHbvAquifer())
    {
        lastHbvAquifer = GetHbvAquifer();
        while (lastHbvAquifer)
        {
            // Algorithm to be performed in case: no input to landscape element from upstream elements
            //      lastHbvAquifer->WaterBalance(timeStep);
            // Algorithm to be performed in case: input to landscape element from upstream elements
    	    lastHbvAquifer->WaterBalance(timeStep, datetime, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge, dayofyear);
            lastHbvAquifer = lastHbvAquifer->GetNextHbvAquifer();
        }
    }
    if (GetKwaAquifer())
    {
        lastKwaAquifer = GetKwaAquifer();
        while (lastKwaAquifer)
        {
            lastKwaAquifer->WaterBalance(timeStep, datetime, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge, dayofyear);
            lastKwaAquifer = lastKwaAquifer->GetNextKwaAquifer();
        }
    }
    //  cout << "  end WaterBalance\n";
}

double DistributedElement::GetPrecipitation() const
{
    double precipitation = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        precipitation = GetLake()->GetPrecipitation();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetGlacier())
    {
        precipitation = precipitation + GetGlacier()->GetPrecipitation();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        precipitation = precipitation + GetHbvAquifer()->GetPrecipitation();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        precipitation = precipitation + GetKwaAquifer()->GetPrecipitation();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetPrecipitation " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return precipitation/sumAreaFraction;
    if (elementFound)
    {
        return precipitation/sumAreaFraction;
	//        return precipitation;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetTemperature() const
{
    double temperature = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        temperature = GetLake()->GetTemperature();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetGlacier())
    {
        temperature = temperature + GetGlacier()->GetTemperature();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        temperature = temperature + GetHbvAquifer()->GetTemperature();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        temperature = temperature + GetKwaAquifer()->GetTemperature();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetTemperature " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return temperature/sumAreaFraction;
    if (elementFound)
    {
        return temperature/sumAreaFraction;
	//        return temperature;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetLakeStorage() const
{
    double lakeStorage = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        lakeStorage = GetLake()->GetLakeStorage();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetLakeStorage " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return lakeStorage/sumAreaFraction;
    if (elementFound)
    {
        return lakeStorage/sumAreaFraction;
	//        return lakeStorage;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetLakeStorageChange() const
{
    double lakeStorageChange = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        lakeStorageChange = GetLake()->GetLakeStorageChange();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetLakeStorageChange " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return lakeStorageChange/sumAreaFraction;
    if (elementFound)
    {
        return lakeStorageChange/sumAreaFraction;
	//        return lakeStorageChange;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetSnowStore() const
{
    double snowStore = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        snowStore = GetGlacier()->GetSnowStore();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        snowStore = snowStore + GetHbvAquifer()->GetSnowStore();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        snowStore = snowStore + GetKwaAquifer()->GetSnowStore();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetSnowStore " << elementFound << "  " << snowStore;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetSnowStore " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return snowStore/sumAreaFraction;
    if (elementFound)
    {
        return snowStore/sumAreaFraction;
	//        return snowStore;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetSnowCoverFraction() const
{
    double snowCoverFraction = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        snowCoverFraction = GetGlacier()->GetSnowCoverFraction();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        snowCoverFraction = snowCoverFraction + GetHbvAquifer()->GetSnowCoverFraction();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        snowCoverFraction = snowCoverFraction + GetKwaAquifer()->GetSnowCoverFraction();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetSnowCoverFraction " << elementFound << "  " << snowCoverFraction;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetSnowCoverFraction " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return snowCoverFraction/sumAreaFraction;
    if (elementFound)
    {
        return snowCoverFraction/sumAreaFraction;
	//        return snowCoverFraction;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetMeltWater() const
{
    double meltWater = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        meltWater = GetGlacier()->GetMeltWater();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        meltWater = meltWater + GetHbvAquifer()->GetMeltWater();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        meltWater = meltWater + GetKwaAquifer()->GetMeltWater();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetMeltWater " << elementFound << "  " << meltWater;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetMeltWater " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return meltWater/sumAreaFraction;
    if (elementFound)
    {
        return meltWater/sumAreaFraction;
	//        return meltWater;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetSnowWaterEquivalentChange() const
{
    double snowWaterEquivalentChange = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        snowWaterEquivalentChange = GetGlacier()->GetSnowWaterEquivalentChange();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        snowWaterEquivalentChange = snowWaterEquivalentChange + GetHbvAquifer()->GetSnowWaterEquivalentChange();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        snowWaterEquivalentChange = snowWaterEquivalentChange + GetKwaAquifer()->GetSnowWaterEquivalentChange();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetSnowWaterEquivalentChange " << elementFound << "  " << snowWaterEquivalentChange;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetSnowWaterEquivalentChange " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return snowWaterEquivalentChange/sumAreaFraction;
    if (elementFound)
    {
        return snowWaterEquivalentChange/sumAreaFraction;
	//        return snowWaterEquivalentChange;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetWaterOutput() const
{
    double waterOutput = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        waterOutput = GetGlacier()->GetWaterOutput();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        waterOutput = waterOutput + GetHbvAquifer()->GetWaterOutput();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        waterOutput = waterOutput + GetKwaAquifer()->GetWaterOutput();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetWaterOutput " << elementFound << "  " << waterOutput;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetWaterOutput " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return waterOutput/sumAreaFraction;
    if (elementFound)
    {
        return waterOutput/sumAreaFraction;
	//        return waterOutput;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierSnowStore() const
{
  double snowStore=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0) {
    snowStore=GetGlacier()->GetSnowStore();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetGlacierSnowStore " << elementFound << "  " << snowStore;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetGlacierSnowStore " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  //  if (elementFound) return snowStore/sumAreaFraction;
  if (elementFound) 
  {
      return snowStore/sumAreaFraction;
      //        return snowStore;
  }
  else 
  {
      return missingData;
  }
}

double DistributedElement::GetGlacierSnowMeltWater() const
{
  double meltWater=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0) {
    meltWater=GetGlacier()->GetMeltWater();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetGlacierSnowMeltWater " << elementFound << "  " << meltWater;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetGlacierSnowMeltWater " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  //  if (elementFound) return meltWater/sumAreaFraction;
  if (elementFound) 
  {
      return meltWater/sumAreaFraction;
      //        return meltWater;
  }
  else 
  {
      return missingData;
  }
}

double DistributedElement::GetAreaGlacierMassBalance() const
{
    double areaGlacierMassBalance = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    //  cout << "start GetGlacierMassBalance\n";
    if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0)
    {
        areaGlacierMassBalance = GetGlacier()->GetGlacierMassBalance();
        //    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetGlacierMassBalance " << elementFound << "  " << areaGlacierMassBalance;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetGlacierMassBalance " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  cout << "  end GetGlacierMassBalance\n";
    //  if (elementFound) return 1.0;
    //  if (elementFound) return areaGlacierMassBalance/sumAreaFraction;
    if (elementFound)
    {
        return areaGlacierMassBalance/sumAreaFraction;
	//        return areaGlacierMassBalance;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierMassBalance() const
{
    double glacierMassBalance = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    //  cout << "start GetGlacierMassBalance\n";
    if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0)
    {
        glacierMassBalance = GetGlacier()->GetGlacierMassBalance();
        //    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetGlacierMassBalance " << elementFound << "  " << glacierMassBalance;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetGlacierMassBalance " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  cout << "  end GetGlacierMassBalance\n";
    //  if (elementFound) return 1.0;
    //  if (elementFound) return glacierMassBalance/sumAreaFraction;
    if (elementFound)
    {
        return glacierMassBalance/sumAreaFraction;
	//        return glacierMassBalance;
    }
    //  if (elementFound) return glacierMassBalance;
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierIceMelt() const
{
    double glacierIceMelt = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    //  cout << "start GetGlacierIceMelt\n";
    if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0)
    {
        glacierIceMelt = GetGlacier()->GetGlacierIceMelt();
        //  sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetGlacierIceMelt " << elementFound << "  " << glacierIceMelt;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetGlacierIceMelt " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  cout << "  end GetGlacierIceMelt\n";
    //  if (elementFound) return 1.0;
    //  if (elementFound) return glacierIceMelt/sumAreaFraction;
    if (elementFound)
    {
        return glacierIceMelt/sumAreaFraction;
	//        return glacierIceMelt;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetAnnualMassBalance() const
{
    double annualMassBalance = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    //  cout << "start GetAnnualMassBalance\n";
    if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0 && GetGlacier()->GetAnnualMassBalance() != missingData)
    {
        annualMassBalance = GetGlacier()->GetAnnualMassBalance();
        //    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetAnnualMassBalance " << elementFound << "  " << annualMassBalance;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetAnnualMassBalance " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  cout << "  end GetAnnualMassBalance\n";
    //  if (elementFound) return 1.0;
    //  if (elementFound) return annualMassBalance/sumAreaFraction;
    if (elementFound)
    {
        return annualMassBalance/sumAreaFraction;
	//        return annualMassBalance;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierIceThickness() const
{
    if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0)
    {
        return GetGlacier()->GetGlacierIceThickness();
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierIceVolume() const
{
    if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0)
    {
        return GetGlacier()->GetGlacierIceVolume();
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetGlacierSurfaceElevation() const
{
    if (GetGlacier() && GetGlacier()->GetGlacierIceAreaFraction() > 0.0)
    {
        return GetGlacier()->GetGlacierSurfaceElevation();
    }
    else
    {
        return missingData;
    }
}

// Algorithm to be performed in case calibration using glacier surface elevation
// Missing data cells should not be output in sourceDew.cpp: WriteAsciiGridAnnualGlacierValues
//double DistributedElement::GetGlacierSurfaceElevation() const
//{
//  if (GetGlacier()) {
//    if (GetGlacier()->GetGlacierIceAreaFraction() > 0.0)
//      return GetGlacier()->GetGlacierSurfaceElevation();
//    else
//      return GetElevation();
//  }
//  else {
//    return missingData;
//  }
//}

double DistributedElement::GetKiWaSoilMoisture(double lengthFraction) const
{
    double soilMoisture = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetKinematicWave())
        {
            soilMoisture = GetGlacier()->GetKiWaSoilMoisture(lengthFraction);
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetKwaAquifer())
    {
        soilMoisture = soilMoisture + GetKwaAquifer()->GetSoilMoisture(lengthFraction);
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetKiWaSoilMoisture " << elementFound << "  " << soilMoisture;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetKiWaSoilMoisture " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return soilMoisture/sumAreaFraction;
    if (elementFound)
    {
        return soilMoisture/sumAreaFraction;
	//        return soilMoisture;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetKiWaGroundWaterDepth(double lengthFraction) const
{
    double groundWaterDepth = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetKinematicWave())
        {
            groundWaterDepth = GetGlacier()->GetKiWaGroundWaterDepth(lengthFraction);
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetKwaAquifer())
    {
        groundWaterDepth = groundWaterDepth + GetKwaAquifer()->GetGroundWaterDepth(lengthFraction);
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetKiWaGroundWaterDepth " << elementFound << "  " << groundWaterDepth;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetKiWaGroundWaterDepth " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return groundWaterDepth/sumAreaFraction;
    if (elementFound)
    {
        return groundWaterDepth/sumAreaFraction;
	//        return groundWaterDepth;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetHbvSoilMoisture() const
{
    double soilMoisture = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetHBV())
        {
            soilMoisture = GetGlacier()->GetHbvSoilMoisture();
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetHbvAquifer())
    {
        soilMoisture = soilMoisture + GetHbvAquifer()->GetSoilMoisture();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetSoilMoisture " << elementFound << "  " << soilMoisture;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetSoilMoisture " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return soilMoisture/sumAreaFraction;
    if (elementFound)
    {
        return soilMoisture/sumAreaFraction;
	//        return soilMoisture;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetHbvSoilMoistureDeficit() const
{
    double soilMoistureDeficit = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetHBV())
        {
            soilMoistureDeficit = GetGlacier()->GetHbvSoilMoistureDeficit();
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetHbvAquifer())
    {
        soilMoistureDeficit = soilMoistureDeficit + GetHbvAquifer()->GetSoilMoistureDeficit();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetSoilMoistureDeficit " << elementFound << "  " << soilMoistureDeficit;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetSoilMoistureDeficit " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return soilMoistureDeficit/sumAreaFraction;
    if (elementFound)
    {
        return soilMoistureDeficit/sumAreaFraction;
	//        return soilMoistureDeficit;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetHbvPercSoilUpper() const
{
    double percSoilUpper = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetHBV())
        {
            percSoilUpper = GetGlacier()->GetHbvPercSoilUpper();
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetHbvAquifer())
    {
        percSoilUpper = percSoilUpper + GetHbvAquifer()->GetPercSoilUpper();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetPercSoilUpper " << elementFound << "  " << percSoilUpper;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetPercSoilUpper " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return percSoilUpper/sumAreaFraction;
    if (elementFound)
    {
        return percSoilUpper/sumAreaFraction;
	//        return percSoilUpper;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetHbvUpperZone() const
{
    double upperZone = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetHBV())
        {
            upperZone = GetGlacier()->GetHbvUpperZone();
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetHbvAquifer())
    {
        upperZone = upperZone + GetHbvAquifer()->GetUpperZone();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetHbvUpperZone " << elementFound << "  " << upperZone;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetHbvUpperZone " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return upperZone/sumAreaFraction;
    if (elementFound)
    {
        return upperZone/sumAreaFraction;
	//        return upperZone;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetHbvLowerZone() const
{
    double lowerZone = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        if (GetGlacier()->GetHBV())
        {
            lowerZone = GetGlacier()->GetHbvLowerZone();
            sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
            elementFound = true;
        }
    }
    if (GetHbvAquifer())
    {
        lowerZone = lowerZone + GetHbvAquifer()->GetLowerZone();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetHbvLowerZone " << elementFound << "  " << lowerZone;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetHbvLowerZone " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return lowerZone/sumAreaFraction;
    if (elementFound)
    {
        return lowerZone/sumAreaFraction;
	//        return lowerZone;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetEvapotranspiration() const
{
    double evapotranspiration = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        evapotranspiration = GetGlacier()->GetInterceptionLoss() + GetGlacier()->GetTranspSoilEvap();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        evapotranspiration = evapotranspiration + GetHbvAquifer()->GetInterceptionLoss() + GetHbvAquifer()->GetTranspSoilEvap();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        evapotranspiration = evapotranspiration + GetKwaAquifer()->GetInterceptionLoss() + GetKwaAquifer()->GetTranspSoilEvap();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetLake())
    {
        evapotranspiration = evapotranspiration + GetLake()->GetLakeEvap();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetInterceptionLoss " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return evapotranspiration/sumAreaFraction;
    if (elementFound)
    {
        return evapotranspiration/sumAreaFraction;
	//        return evapotranspiration;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetPET() const
{
	double pet = 0.0;
	double sumAreaFraction = 0.0;
	bool elementFound = false;
	if (GetGlacier()) {
		sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
		elementFound = true;
	}
	if (GetHbvAquifer()) {
		pet = pet + GetHbvAquifer()->GetPET();
		sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
		elementFound = true;
	}
	if (GetLake()) {
		pet = pet + GetLake()->GetLakeEvap();
		sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
		elementFound = true;
	}
	if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction>1.0 + epsilon)) {
		cout << " sumAreaFraction GetInterceptionLoss " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
		exit(1);
	}
	//  if (elementFound) return evapotranspiration/sumAreaFraction;
	if (elementFound) return pet;
	else return missingData;
}

double DistributedElement::GetInterceptionLoss() const
{
    double interceptionLoss = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        interceptionLoss = GetGlacier()->GetInterceptionLoss();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        interceptionLoss = interceptionLoss + GetHbvAquifer()->GetInterceptionLoss();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        interceptionLoss = interceptionLoss + GetKwaAquifer()->GetInterceptionLoss();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetInterceptionLoss " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return interceptionLoss/sumAreaFraction;
    if (elementFound)
    {
        return interceptionLoss/sumAreaFraction;
	//        return interceptionLoss;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetTranspSoilEvap() const
{
    double transpSoilEvap = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetGlacier())
    {
        transpSoilEvap = GetGlacier()->GetTranspSoilEvap();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        transpSoilEvap = transpSoilEvap + GetHbvAquifer()->GetTranspSoilEvap();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        transpSoilEvap = transpSoilEvap + GetKwaAquifer()->GetTranspSoilEvap();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetTranspSoilEvap " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return transpSoilEvap/sumAreaFraction;
    if (elementFound)
    {
        return transpSoilEvap/sumAreaFraction;
	//        return transpSoilEvap;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetLakeEvap() const
{
    double lakeEvap = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        lakeEvap = GetLake()->GetLakeEvap();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetLakeEvap " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return lakeEvap/sumAreaFraction;
    if (elementFound)
    {
        return lakeEvap/sumAreaFraction;
	//        return lakeEvap;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetRunoff() const
{
    double runoff = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        runoff = GetLake()->GetRunoff();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetGlacier())
    {
        runoff = runoff + GetGlacier()->GetRunoff();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        runoff = runoff + GetHbvAquifer()->GetRunoff();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        runoff = runoff + GetKwaAquifer()->GetRunoff();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetRunoff " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return runoff/sumAreaFraction;
    if (elementFound)
    {
        return runoff/sumAreaFraction;
	//        return runoff;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetLowerRunoff() const
{
    double lowerRunoff = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        lowerRunoff = GetLake()->GetRunoff();
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetGlacier())
    {
        lowerRunoff = lowerRunoff + GetGlacier()->GetLowerRunoff();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        lowerRunoff = lowerRunoff + GetHbvAquifer()->GetLowerRunoff();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        lowerRunoff = lowerRunoff + GetKwaAquifer()->GetLowerRunoff();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetLowerRunoff " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return lowerRunoff/sumAreaFraction;
    if (elementFound)
    {
        return lowerRunoff/sumAreaFraction;
	//        return lowerRunoff;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetUpperRunoff() const
{
    double upperRunoff = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    if (GetLake())
    {
        upperRunoff = 0.0;
        sumAreaFraction = sumAreaFraction + GetLake()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetGlacier())
    {
        upperRunoff = upperRunoff + GetGlacier()->GetUpperRunoff();
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetHbvAquifer())
    {
        upperRunoff = upperRunoff + GetHbvAquifer()->GetUpperRunoff();
        sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
        elementFound = true;
    }
    if (GetKwaAquifer())
    {
        upperRunoff = upperRunoff + GetKwaAquifer()->GetUpperRunoff();
        sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
        elementFound = true;
    }
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetUpperRunoff " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  if (elementFound) return upperRunoff/sumAreaFraction;
    if (elementFound)
    {
        return upperRunoff/sumAreaFraction;
	//        return upperRunoff;
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetDischarge() const
{
    if (GetRunoff() != missingData)
    {
        return GetRunoff() * area / commonPar->GetSECONDS_TIMESTEP();    /*  runoff (m) -> discharge (m3/s)  */
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetLowerDischarge() const
{
    if (GetLowerRunoff() != missingData)
    {
        return GetLowerRunoff() * area / commonPar->GetSECONDS_TIMESTEP();    /*  runoff (m) -> discharge (m3/s)  */
    }
    else
    {
        return missingData;
    }
}

double DistributedElement::GetUpperDischarge() const
{
    if (GetUpperRunoff() != missingData)
    {
        return GetUpperRunoff() * area / commonPar->GetSECONDS_TIMESTEP();    /*  runoff (m) -> discharge (m3/s)  */
    }
    else
    {
        return missingData;
    }
}

int DistributedElement::GetNumberSum() const
{
    return numberSum;
}

void DistributedElement::SetSumWaterBalance()
{
  if (GetRunoff() != missingData) {
    sumPrecipitation = sumPrecipitation + GetPrecipitation();
    sumEvapotranspiration = sumEvapotranspiration + GetEvapotranspiration();
    sumRunoff = sumRunoff + GetRunoff();
    if(sumRunoff<0) printf("classdist iha %f\n",sumRunoff);
    sumGlacierIceMelt = sumGlacierIceMelt + GetGlacierIceMelt();
    numberSum++;
  }
}

void DistributedElement::SetInitialStorage()
{
  initialStorage = GetSnowStore() + GetMeltWater();
  if (GetHbvAquifer()) {
    initialStorage = initialStorage + GetHbvSoilMoisture() + GetHbvUpperZone() + GetHbvLowerZone() + GetLakeStorage();
  }
  if (GetLake()) {
    initialStorage = initialStorage + GetLakeStorage();
  }
  if (GetGlacier()) {
    initialStorage = initialStorage - GetSumGlacierIceMelt();
  }
}

double DistributedElement::GetInitialStorage() const
{
    return initialStorage;
}

void DistributedElement::SetFinalStorage()
{
  finalStorage = GetSnowStore() + GetMeltWater();
  if (GetHbvAquifer()) {
    finalStorage = finalStorage + GetHbvSoilMoisture() + GetHbvUpperZone() + GetHbvLowerZone() + GetLakeStorage();
  }
  if (GetLake()) {
    finalStorage = finalStorage + GetLakeStorage();
  }
  if (GetGlacier()) {
    finalStorage = finalStorage - GetSumGlacierIceMelt();
  }
}

double DistributedElement::GetFinalStorage() const
{
    return finalStorage;
}

double DistributedElement::GetSumPrecipitation() const
{
    return sumPrecipitation;
}

double DistributedElement::GetSumEvapotranspiration() const
{
    return sumEvapotranspiration;
}

double DistributedElement::GetSumRunoff() const
{
    return sumRunoff;
}

double DistributedElement::GetSumGlacierIceMelt() const
{
    return sumGlacierIceMelt;
}

void DistributedElement::SetTotalReservoirStorage(int timeStep, int initialTimeSteps, int numberTimeSteps, bool * inputDataFound, bool * firstTotal)
{
    double totalStorage;
    // Initial value
    if (*firstTotal)
    {
        if (timeStep > 0)
        {
            totalStorage = totalReservoir->GetTotalReservoirStorage(timeStep - 1);
        }
        else
        {
            totalStorage = totalReservoir->GetInitialTotalReservoirStorage();
        }
    }
    else
    {
        totalStorage = totalReservoir->GetTotalReservoirStorage(timeStep);
    }
    if (*inputDataFound)
    {
        //      if (GetLake()) totalStorage = totalStorage + GetLake()->GetLakeStorageChange()*1000.0;
        if (GetLake())
        {
            totalStorage = totalStorage + GetLakeStorageChange() * GetArea();
        }
    }
    totalReservoir->SetTotalReservoirStorage(timeStep, totalStorage);
}

/*void DistributedElement::SetTotalReservoirStorage(int timeStep, int initialTimeSteps, int numberTimeSteps, bool * inputDataFound, bool * firstTotal)
{
double totalStorage;
// Initial value
if (*firstTotal) {
if (*inputDataFound && GetLake()) {
totalStorage = totalStorage + GetLakeStorage()*GetArea();
}
else {
totalStorage = 0.0;
}
}
else {
totalStorage = totalReservoir->GetTotalReservoirStorage(timeStep);
if (*inputDataFound && GetLake()) {
totalStorage = totalStorage + GetLakeStorage()*GetArea();
}
}
totalReservoir->SetTotalReservoirStorage(timeStep, totalStorage);
}*/

/*void DistributedElement::SetTotalReservoirStorage(int timeStep, int initialTimeSteps, int numberTimeSteps, bool * inputDataFound, bool * firstTotal)
{
double totalStorage;
if (*firstTotal) {
if (*inputDataFound && GetLake()) {
totalStorage = totalStorage + GetLake()->GetLakeStorage()*area;
}
else {
totalStorage = 0.0;
}
}
else {
totalStorage = totalReservoir->GetTotalReservoirStorage(timeStep);
if (*inputDataFound && GetLake()) {
totalStorage = totalStorage + GetLake()->GetLakeStorage()*area;
}
}
totalReservoir->SetTotalReservoirStorage(timeStep, totalStorage);
}*/

