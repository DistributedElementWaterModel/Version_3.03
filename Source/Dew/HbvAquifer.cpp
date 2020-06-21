#include "HbvAquifer.h"
#include "Snow.h"
#include "HBV.h"
#include "Vegetation.h"
#include "stdafx.h"
#include "DateTime.h"

HbvAquifer::HbvAquifer()
    //    HbvAquifer::HbvAquifer() :
    //    surfaceAlbedo(0.0)
{
    SetLandScapeElement(0);
    //    SetLandAtmosphereInterfaceElement(0);
    //    SetAlbedoElement(0);
    //    SetAtmosphereElement(0);
    //    SetPenmanMonteithElement(0);
    SetAreaFraction(0.0);
    SetNextHbvAquifer(0);
    SetHBV(0);
}

HbvAquifer::~HbvAquifer()
{
}

void  HbvAquifer::SetNextHbvAquifer(HbvAquifer *theHbvAquifer)
{
    nextHbvAquifer = theHbvAquifer;
}
HbvAquifer * HbvAquifer::GetNextHbvAquifer() const
{
    return nextHbvAquifer;
}
void  HbvAquifer::SetVegetation(Vegetation *theVegetation)
{
    ptrVeg = theVegetation;
}
Vegetation * HbvAquifer::GetVegetation() const
{
    return ptrVeg;
}
void  HbvAquifer::SetSnow(Snow *theSnow)
{
    ptrSnow = theSnow;
}
Snow * HbvAquifer::GetSnow() const
{
    return ptrSnow;
}
void  HbvAquifer::SetHBV(HBV *theHBV)
{
    ptrHbv = theHBV;
}
HBV * HbvAquifer::GetHBV() const
{
    return ptrHbv;
}
void  HbvAquifer::SetAreaFraction(double value)
{
    areaFraction = value;
}
double  HbvAquifer::GetAreaFraction() const
{
    return areaFraction;
}
void  HbvAquifer::SetLandSurfaceType(LANDSURFACE landSurf)
{
    landSurfaceType = landSurf;
}
LANDSURFACE  HbvAquifer::GetLandSurfaceType() const
{
    return landSurfaceType;
}
void  HbvAquifer::SetSoilType(SOIL soil)
{
    soilType = soil;
}
SOIL  HbvAquifer::GetSoilType() const
{
    return soilType;
}
/*void  HbvAquifer::SetSurfaceAlbedo(double value)
{
    surfaceAlbedo = value;
}
double  HbvAquifer::GetSurfaceAlbedo() const
{
    return surfaceAlbedo;
    }*/
void  HbvAquifer::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * HbvAquifer::GetLandScapeElement() const
{
    return landScapeElement;
}
/*void  HbvAquifer::SetLandAtmosphereInterfaceElement(LandAtmosphereInterface *theElement)
{
    landAtmosphereInterfaceElement = theElement;
}
LandAtmosphereInterface * HbvAquifer::GetLandAtmosphereInterfaceElement() const
{
    return landAtmosphereInterfaceElement;
}
void  HbvAquifer::SetAlbedoElement(Albedo *theElement)
{
    albedoElement = theElement;
}
Albedo * HbvAquifer::GetAlbedoElement() const
{
    return albedoElement;
}
void  HbvAquifer::SetAtmosphereElement(Atmosphere *theElement)
{
    atmosphereElement = theElement;
}
Atmosphere * HbvAquifer::GetAtmosphereElement() const
{
    return atmosphereElement;
}
void  HbvAquifer::SetPenmanMonteithElement(PenmanMonteith *theElement)
{
    penmanMonteithElement = theElement;
}
PenmanMonteith * HbvAquifer::GetPenmanMonteithElement() const
{
    return penmanMonteithElement;
}*/

void HbvAquifer::SetSnowStore(double value)
{
    if (GetNextHbvAquifer())
    {
        GetNextHbvAquifer()->SetSnowStore(value);
    }
    GetSnow()->SetSnowStore(value);
}

void HbvAquifer::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
    if (GetNextHbvAquifer())
    {
        GetNextHbvAquifer()->SetSubSurfaceHbvStore(sm, uz, lz);
    }
    GetHBV()->SetSubSurfaceHbvStore(sm, uz, lz);
}

// Algorithm to be performed in case: no input to landscape element from upstream elements
/*void HbvAquifer::WaterBalance(int timeStep)
{
GetVegetation()->WaterBalance(timeStep);
GetSnow()->WaterBalance(timeStep, GetVegetation()->GetThroughFall());
GetHBV()->WaterBalance(timeStep, GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), GetVegetation()->GetDryPeriod());
}*/
// End algorithm to be performed in case: no input to landscape element from upstream elements

// Algorithm to be performed in case: input to landscape element from upstream elements
void HbvAquifer::WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge)
{
    GetVegetation()->WaterBalance(timeStep, datetime, missingData, missingData);
    GetSnow()->WaterBalance(timeStep, datetime, GetVegetation()->GetThroughFall(), GetVegetation()->GetTemperature());
    GetHBV()->WaterBalance(timeStep, datetime, GetSnow()->GetWaterOutput(), GetVegetation()->GetTemperature(), GetSnow()->GetSnowCoverFraction(), GetVegetation()->GetDryPeriod(), upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
}
// End algorithm to be performed in case: input to landscape element from upstream elements

double HbvAquifer::GetTotalHbvAreaFraction() const
{
    double totalHbvAreaFraction = 0.0;
    if (GetNextHbvAquifer())
    {
        totalHbvAreaFraction = totalHbvAreaFraction + GetNextHbvAquifer()->GetTotalHbvAreaFraction();
    }
    totalHbvAreaFraction = totalHbvAreaFraction + GetAreaFraction();
    return totalHbvAreaFraction;
}

double HbvAquifer::GetPrecipitation() const
{
    double precipitation = 0.0;
    if (GetNextHbvAquifer())
    {
        precipitation = precipitation + GetNextHbvAquifer()->GetPrecipitation();
    }
    precipitation = precipitation + (GetVegetation()->GetPrecipitation() * GetAreaFraction() / 100.0);
    return precipitation;
}

double HbvAquifer::GetTemperature() const
{
    double temperature = 0.0;
    if (GetNextHbvAquifer())
    {
        temperature = temperature + GetNextHbvAquifer()->GetTemperature();
    }
    temperature = temperature + (GetVegetation()->GetTemperature() * GetAreaFraction() / 100.0);
    return temperature;
}

double HbvAquifer::GetSnowCoverFraction() const
{
    double snowCoverFraction = 0.0;
    if (GetNextHbvAquifer())
    {
        snowCoverFraction = snowCoverFraction + GetNextHbvAquifer()->GetSnowCoverFraction();
    }
    snowCoverFraction = snowCoverFraction + (GetSnow()->GetSnowCoverFraction() * GetAreaFraction() / 100.0);
    return snowCoverFraction;
}

double HbvAquifer::GetSnowStore() const
{
    double snowStore = 0.0;
    if (GetNextHbvAquifer())
    {
        snowStore = snowStore + GetNextHbvAquifer()->GetSnowStore();
    }
    snowStore = snowStore + (GetSnow()->GetSnowStore() * GetAreaFraction() / 100.0);
    return snowStore;
}

double HbvAquifer::GetMeltWater() const
{
    double meltWater = 0.0;
    if (GetNextHbvAquifer())
    {
        meltWater = meltWater + GetNextHbvAquifer()->GetMeltWater();
    }
    meltWater = meltWater + (GetSnow()->GetMeltWater() * GetAreaFraction() / 100.0);
    return meltWater;
}

double HbvAquifer::GetSnowWaterEquivalentChange() const
{
    double snowWaterEquivalentChange = 0.0;
    if (GetNextHbvAquifer())
    {
        snowWaterEquivalentChange = snowWaterEquivalentChange + GetNextHbvAquifer()->GetSnowWaterEquivalentChange();
    }
    snowWaterEquivalentChange = snowWaterEquivalentChange + (GetSnow()->GetSnowWaterEquivalentChange() * GetAreaFraction() / 100.0);
    return snowWaterEquivalentChange;
}

double HbvAquifer::GetWaterOutput() const
{
    double waterOutput = 0.0;
    if (GetNextHbvAquifer())
    {
        waterOutput = waterOutput + GetNextHbvAquifer()->GetWaterOutput();
    }
    waterOutput = waterOutput + (GetSnow()->GetWaterOutput() * GetAreaFraction() / 100.0);
    return waterOutput;
}

double HbvAquifer::GetSoilMoisture() const
{
    double soilMoisture = 0.0;
    if (GetNextHbvAquifer())
    {
        soilMoisture = soilMoisture + GetNextHbvAquifer()->GetSoilMoisture();
    }
    soilMoisture = soilMoisture + (GetHBV()->GetSoilMoisture() * GetAreaFraction() / 100.0);
    return soilMoisture;
}

double HbvAquifer::GetSoilMoistureDeficit() const
{
    double soilMoistureDeficit = 0.0;
    if (GetNextHbvAquifer())
    {
        soilMoistureDeficit = soilMoistureDeficit + GetNextHbvAquifer()->GetSoilMoistureDeficit();
    }
    soilMoistureDeficit = soilMoistureDeficit + (GetHBV()->GetSoilMoistureDeficit() * GetAreaFraction() / 100.0);
    return soilMoistureDeficit;
}

double HbvAquifer::GetPercSoilUpper() const
{
    double percSoilUpper = 0.0;
    if (GetNextHbvAquifer())
    {
        percSoilUpper = percSoilUpper + GetNextHbvAquifer()->GetPercSoilUpper();
    }
    percSoilUpper = percSoilUpper + (GetHBV()->GetPercSoilUpper() * GetAreaFraction() / 100.0);
    return percSoilUpper;
}

double HbvAquifer::GetUpperZone() const
{
    double upperZone = 0.0;
    if (GetNextHbvAquifer())
    {
        upperZone = upperZone + GetNextHbvAquifer()->GetUpperZone();
    }
    upperZone = upperZone + (GetHBV()->GetUpperZone() * GetAreaFraction() / 100.0);
    return upperZone;
}

double HbvAquifer::GetLowerZone() const
{
    double lowerZone = 0.0;
    if (GetNextHbvAquifer())
    {
        lowerZone = lowerZone + GetNextHbvAquifer()->GetLowerZone();
    }
    lowerZone = lowerZone + (GetHBV()->GetLowerZone() * GetAreaFraction() / 100.0);
    return lowerZone;
}

double HbvAquifer::GetInterceptionLoss() const
{
    double interceptionLoss = 0.0;
    if (GetNextHbvAquifer())
    {
        interceptionLoss = interceptionLoss + GetNextHbvAquifer()->GetInterceptionLoss();
    }
    interceptionLoss = interceptionLoss + (GetVegetation()->GetInterceptionLoss() * GetAreaFraction() / 100.0);
    return interceptionLoss;
}

double HbvAquifer::GetTranspSoilEvap() const
{
    double transpSoilEvap = 0.0;
    if (GetNextHbvAquifer())
    {
        transpSoilEvap = transpSoilEvap + GetNextHbvAquifer()->GetTranspSoilEvap();
    }
    transpSoilEvap = transpSoilEvap + (GetHBV()->GetTranspSoilEvap() * GetAreaFraction() / 100.0);
    return transpSoilEvap;
}

double HbvAquifer::GetRunoff() const
{
    double runoff = 0.0;
    if (GetNextHbvAquifer())
    {
        runoff = runoff + GetNextHbvAquifer()->GetRunoff();
    }
    runoff = runoff + (GetHBV()->GetRunoff() * GetAreaFraction() / 100.0);
    return runoff;
}

double HbvAquifer::GetLowerRunoff() const
{
    double lowerRunoff = 0.0;
    if (GetNextHbvAquifer())
    {
        lowerRunoff = lowerRunoff + GetNextHbvAquifer()->GetLowerRunoff();
    }
    lowerRunoff = lowerRunoff + (GetHBV()->GetLowerRunoff() * GetAreaFraction() / 100.0);
    return lowerRunoff;
}

double HbvAquifer::GetUpperRunoff() const
{
    double upperRunoff = 0.0;
    if (GetNextHbvAquifer())
    {
        upperRunoff = upperRunoff + GetNextHbvAquifer()->GetUpperRunoff();
    }
    upperRunoff = upperRunoff + (GetHBV()->GetUpperRunoff() * GetAreaFraction() / 100.0);
    return upperRunoff;
}
