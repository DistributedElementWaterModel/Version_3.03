#include "KwaAquifer.h"

#include "Vegetation.h"
#include "Snow.h"
#include "KinematicWave.h"
#include "stdafx.h"
#include "DateTime.h"

KwaAquifer::KwaAquifer()
    //    KwaAquifer::KwaAquifer() :
    //    surfaceAlbedo(0.0)
{
    SetLandScapeElement(0);
    //    SetLandAtmosphereInterfaceElement(0);
    //    SetAlbedoElement(0);
    //    SetAtmosphereElement(0);
    //    SetPenmanMonteithElement(0);
    SetAreaFraction(0);
    SetNextKwaAquifer(0);
    SetKinematicWave(0);
}

KwaAquifer::~KwaAquifer()
{
}

void  KwaAquifer::SetNextKwaAquifer(KwaAquifer *theKwaAquifer)
{
    nextKwaAquifer = theKwaAquifer;
}
KwaAquifer * KwaAquifer::GetNextKwaAquifer() const
{
    return nextKwaAquifer;
}
void  KwaAquifer::SetVegetation(Vegetation *theVegetation)
{
    ptrVeg = theVegetation;
}
Vegetation * KwaAquifer::GetVegetation() const
{
    return ptrVeg;
}
void  KwaAquifer::SetSnow(Snow *theSnow)
{
    ptrSnow = theSnow;
}
Snow * KwaAquifer::GetSnow() const
{
    return ptrSnow;
}
void  KwaAquifer::SetKinematicWave(KinematicWave *theKinematicWave)
{
    ptrKwa = theKinematicWave;
}
KinematicWave * KwaAquifer::GetKinematicWave() const
{
    return ptrKwa;
}
void  KwaAquifer::SetAreaFraction(double value)
{
    areaFraction = value;
}
double  KwaAquifer::GetAreaFraction() const
{
    return areaFraction;
}
void  KwaAquifer::SetLandSurfaceType(LANDSURFACE landSurf)
{
    landSurfaceType = landSurf;
}
LANDSURFACE  KwaAquifer::GetLandSurfaceType() const
{
    return landSurfaceType;
}
void  KwaAquifer::SetSoilType(SOIL soil)
{
    soilType = soil;
}
SOIL  KwaAquifer::GetSoilType() const
{
    return soilType;
}
/*void  KwaAquifer::SetSurfaceAlbedo(double value)
{
    surfaceAlbedo = value;
}
double  KwaAquifer::GetSurfaceAlbedo() const
{
    return surfaceAlbedo;
    }*/
void  KwaAquifer::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * KwaAquifer::GetLandScapeElement() const
{
    return landScapeElement;
}
/*void  KwaAquifer::SetLandAtmosphereInterfaceElement(LandAtmosphereInterface *theElement)
{
    landAtmosphereInterfaceElement = theElement;
}
LandAtmosphereInterface * KwaAquifer::GetLandAtmosphereInterfaceElement() const
{
    return landAtmosphereInterfaceElement;
}
void  KwaAquifer::SetAlbedoElement(Albedo *theElement)
{
    albedoElement = theElement;
}
Albedo * KwaAquifer::GetAlbedoElement() const
{
    return albedoElement;
}
void  KwaAquifer::SetAtmosphereElement(Atmosphere *theElement)
{
    atmosphereElement = theElement;
}
Atmosphere * KwaAquifer::GetAtmosphereElement() const
{
    return atmosphereElement;
}
void  KwaAquifer::SetPenmanMonteithElement(PenmanMonteith *theElement)
{
    penmanMonteithElement = theElement;
}
PenmanMonteith * KwaAquifer::GetPenmanMonteithElement() const
{
    return penmanMonteithElement;
}*/

void KwaAquifer::SetSnowStore(double value)
{
    if (GetNextKwaAquifer())
    {
        GetNextKwaAquifer()->SetSnowStore(value);
    }
    GetSnow()->SetSnowStore(value);
}

void KwaAquifer::WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge, int dayofyear)
{
    //  cout << "kwaAquifer\n";
    GetVegetation()->WaterBalance(timeStep, datetime, missingData, missingData, dayofyear, GetSnow()->GetSnowCoverFraction());
    GetSnow()->WaterBalance(timeStep, datetime, GetVegetation()->GetThroughFall(), GetVegetation()->GetTemperature());
    GetKinematicWave()->WaterBalance(timeStep, datetime, GetSnow()->GetWaterOutput(), GetVegetation()->GetTemperature(), GetSnow()->GetSnowCoverFraction(), GetVegetation()->GetDryPeriod(), upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
}

double KwaAquifer::GetTotalKwaAreaFraction() const
{
    double totalKwaAreaFraction = 0.0;
    if (GetNextKwaAquifer())
    {
        totalKwaAreaFraction = totalKwaAreaFraction + GetNextKwaAquifer()->GetTotalKwaAreaFraction();
    }
    totalKwaAreaFraction = totalKwaAreaFraction + GetAreaFraction();
    return totalKwaAreaFraction;
}

double KwaAquifer::GetPrecipitation() const
{
    double precipitation = 0.0;
    if (GetNextKwaAquifer())
    {
        precipitation = precipitation + GetNextKwaAquifer()->GetPrecipitation();
    }
    precipitation = precipitation + (GetVegetation()->GetPrecipitation() * GetAreaFraction() / 100.0);
    return precipitation;
}

double KwaAquifer::GetTemperature() const
{
    double temperature = 0.0;
    if (GetNextKwaAquifer())
    {
        temperature = temperature + GetNextKwaAquifer()->GetTemperature();
    }
    temperature = temperature + (GetVegetation()->GetTemperature() * GetAreaFraction() / 100.0);
    return temperature;
}

double KwaAquifer::GetSnowCoverFraction() const
{
    double snowCoverFraction = 0.0;
    if (GetNextKwaAquifer())
    {
        snowCoverFraction = snowCoverFraction + GetNextKwaAquifer()->GetSnowCoverFraction();
    }
    snowCoverFraction = snowCoverFraction + (GetSnow()->GetSnowCoverFraction() * GetAreaFraction() / 100.0);
    return snowCoverFraction;
}

double KwaAquifer::GetSnowStore() const
{
    double snowStore = 0.0;
    if (GetNextKwaAquifer())
    {
        snowStore = snowStore + GetNextKwaAquifer()->GetSnowStore();
    }
    snowStore = snowStore + (GetSnow()->GetSnowStore() * GetAreaFraction() / 100.0);
    return snowStore;
}

double KwaAquifer::GetMeltWater() const
{
    double meltWater = 0.0;
    if (GetNextKwaAquifer())
    {
        meltWater = meltWater + GetNextKwaAquifer()->GetMeltWater();
    }
    meltWater = meltWater + (GetSnow()->GetMeltWater() * GetAreaFraction() / 100.0);
    return meltWater;
}

double KwaAquifer::GetSnowWaterEquivalentChange() const
{
    double snowWaterEquivalentChange = 0.0;
    if (GetNextKwaAquifer())
    {
        snowWaterEquivalentChange = snowWaterEquivalentChange + GetNextKwaAquifer()->GetSnowWaterEquivalentChange();
    }
    snowWaterEquivalentChange = snowWaterEquivalentChange + (GetSnow()->GetSnowWaterEquivalentChange() * GetAreaFraction() / 100.0);
    return snowWaterEquivalentChange;
}

double KwaAquifer::GetWaterOutput() const
{
    double waterOutput = 0.0;
    if (GetNextKwaAquifer())
    {
        waterOutput = waterOutput + GetNextKwaAquifer()->GetWaterOutput();
    }
    waterOutput = waterOutput + (GetSnow()->GetWaterOutput() * GetAreaFraction() / 100.0);
    return waterOutput;
}

double KwaAquifer::GetSoilMoisture(double lengthFraction) const
{
    double soilMoisture = 0.0;
    if (GetNextKwaAquifer())
    {
        soilMoisture = soilMoisture + GetNextKwaAquifer()->GetSoilMoisture(lengthFraction);
    }
    soilMoisture = soilMoisture + (GetKinematicWave()->GetSoilMoisture(lengthFraction) * GetAreaFraction() / 100.0);
    return soilMoisture;
}

double KwaAquifer::GetGroundWaterDepth(double lengthFraction) const
{
    double groundWaterDepth = 0.0;
    if (GetNextKwaAquifer())
    {
        groundWaterDepth = groundWaterDepth + GetNextKwaAquifer()->GetGroundWaterDepth(lengthFraction);
    }
    groundWaterDepth = groundWaterDepth + (GetKinematicWave()->GetGroundWaterDepth(lengthFraction) * GetAreaFraction() / 100.0);
    return groundWaterDepth;
}

double KwaAquifer::GetInterceptionLoss() const
{
    double interceptionLoss = 0.0;
    if (GetNextKwaAquifer())
    {
        interceptionLoss = interceptionLoss + GetNextKwaAquifer()->GetInterceptionLoss();
    }
    interceptionLoss = interceptionLoss + (GetVegetation()->GetInterceptionLoss() * GetAreaFraction() / 100.0);
    return interceptionLoss;
}

double KwaAquifer::GetTranspSoilEvap() const
{
    double transpSoilEvap = 0.0;
    if (GetNextKwaAquifer())
    {
        transpSoilEvap = transpSoilEvap + GetNextKwaAquifer()->GetTranspSoilEvap();
    }
    transpSoilEvap = transpSoilEvap + (GetKinematicWave()->GetTranspSoilEvap() * GetAreaFraction() / 100.0);
    return transpSoilEvap;
}

double KwaAquifer::GetRunoff() const
{
    double runoff = 0;
    if (GetNextKwaAquifer())
    {
        runoff = runoff + GetNextKwaAquifer()->GetRunoff();
    }
    runoff = runoff + GetKinematicWave()->GetRunoff() * GetAreaFraction() / 100.0;
    return runoff;
}

double KwaAquifer::GetLowerRunoff() const
{
    double lowerRunoff = 0;
    if (GetNextKwaAquifer())
    {
        lowerRunoff = lowerRunoff + GetNextKwaAquifer()->GetLowerRunoff();
    }
    lowerRunoff = lowerRunoff + GetKinematicWave()->GetLowerRunoff() * GetAreaFraction() / 100.0;
    return lowerRunoff;
}

double KwaAquifer::GetUpperRunoff() const
{
    double upperRunoff = 0;
    if (GetNextKwaAquifer())
    {
        upperRunoff = upperRunoff + GetNextKwaAquifer()->GetUpperRunoff();
    }
    upperRunoff = upperRunoff + GetKinematicWave()->GetUpperRunoff() * GetAreaFraction() / 100.0;
    return upperRunoff;
}
