#include "Glacier.h"
#include "GlacierIce.h"
#include "DateTime.h"
#include "Dew.h"
#include "HBV.h"
#include "Snow.h"
#include "KinematicWave.h"
#include "Vegetation.h"


Glacier::Glacier(DateTime startModelTime, DateTime startSimulationTime, DateTime endSimulationTime) :
  //    surfaceAlbedo(0.0),
    glacierMassBalance(0.0),
    previousYearMassBalance(0.0),
    annualMassBalance(missingData)
{
    int i;
    int numberYears = endSimulationTime.getYear() - startModelTime.getYear() + 1;
    SetLandScapeElement(0);
    //    SetLandAtmosphereInterfaceElement(0);
    //    SetAlbedoElement(0);
    //    SetAtmosphereElement(0);
    //    SetPenmanMonteithElement(0);
    SetAreaFraction(0.0);
    SetGlacierIceAreaFraction(0.0);
    SetHBV(0);
    SetKinematicWave(0);
    SetVegetation(0);
}

Glacier::~Glacier()
{
}

void Glacier::SetSnow(Snow *theSnow)
{
    ptrSnow = theSnow;
}
Snow *Glacier::GetSnow() const
{
    return ptrSnow;
}
//  void Glacier::SetGlacierSurface(GlacierSurface *theGlacierSurface) { ptrSurface = theGlacierSurface; }
//  GlacierSurface *Glacier::GetGlacierSurface() const { return ptrSurface; }
void Glacier::SetGlacierIce(GlacierIce *theGlacierIce)
{
    ptrIce = theGlacierIce;
}
GlacierIce *Glacier::GetGlacierIce() const
{
    return ptrIce;
}
void Glacier::SetHBV(HBV *theHBV)
{
    ptrHbv = theHBV;
}
HBV *Glacier::GetHBV() const
{
    return ptrHbv;
}
void Glacier::SetKinematicWave(KinematicWave *theKinematicWave)
{
    ptrKwa = theKinematicWave;
}
KinematicWave *Glacier::GetKinematicWave() const
{
    return ptrKwa;
}
void Glacier::SetVegetation(Vegetation *theVegetation)
{
    ptrVeg = theVegetation;
}
Vegetation *Glacier::GetVegetation() const
{
    return ptrVeg;
}
void Glacier::SetAreaFraction(double value)
{
    areaFraction = value;
}
double Glacier::GetAreaFraction() const
{
    return areaFraction;
}
void Glacier::SetLandSurfaceType(LANDSURFACE landSurf)
{
    landSurfaceType = landSurf;
}
LANDSURFACE Glacier::GetLandSurfaceType() const
{
    return landSurfaceType;
}
void Glacier::SetSoilType(SOIL soil)
{
    soilType = soil;
}
SOIL Glacier::GetSoilType() const
{
    return soilType;
}
void Glacier::SetGlacierIceAreaFraction(double value)
{
    glacierIceAreaFraction = value;
}
double Glacier::GetGlacierIceAreaFraction() const
{
    return glacierIceAreaFraction;
}

/*void Glacier::SetSurfaceAlbedo(double value)
{
    surfaceAlbedo = value;
}
double Glacier::GetSurfaceAlbedo() const
{
    return surfaceAlbedo;
    }*/
void Glacier::SetPreviousYearMassBalance(double value)
{
    previousYearMassBalance = value;
}
double Glacier::GetPreviousYearMassBalance() const
{
    return previousYearMassBalance;
}
void Glacier::SetAnnualMassBalance(double value)
{
    annualMassBalance = value;
}
void Glacier::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement *Glacier::GetLandScapeElement() const
{
    return landScapeElement;
}
/*void Glacier::SetLandAtmosphereInterfaceElement(LandAtmosphereInterface *theElement)
{
    landAtmosphereInterfaceElement = theElement;
}
LandAtmosphereInterface *Glacier::GetLandAtmosphereInterfaceElement() const
{
    return landAtmosphereInterfaceElement;
}
void Glacier::SetAlbedoElement(Albedo *theElement)
{
    albedoElement = theElement;
}
Albedo *Glacier::GetAlbedoElement() const
{
    return albedoElement;
}
void Glacier::SetAtmosphereElement(Atmosphere *theElement)
{
    atmosphereElement = theElement;
}
Atmosphere *Glacier::GetAtmosphereElement() const
{
    return atmosphereElement;
}
void Glacier::SetPenmanMonteithElement(PenmanMonteith *theElement)
{
    penmanMonteithElement = theElement;
}
PenmanMonteith *Glacier::GetPenmanMonteithElement() const
{
    return penmanMonteithElement;
    }*/


void Glacier::SetGlacierSurfaceElevation(double value)
{
    GetGlacierIce()->SetGlacierSurfaceElevation(value);
}

double Glacier::GetGlacierSurfaceElevation() const
{
    double glacierSurfaceElevation = GetGlacierIce()->GetGlacierSurfaceElevation();
    return glacierSurfaceElevation;
}

void Glacier::SetGlacierSurfaceElevationChangeNormalized(double value)
{
    GetGlacierIce()->SetGlacierSurfaceElevationChangeNormalized(value);
}

double Glacier::GetGlacierSurfaceElevationChangeNormalized() const
{
    double glacierSurfaceElevationChangeNormalized = GetGlacierIce()->GetGlacierSurfaceElevationChangeNormalized();
    return glacierSurfaceElevationChangeNormalized;
}

void Glacier::SetRestrictedElevationChange(double value)
{
    GetGlacierIce()->SetRestrictedElevationChange(value);
}

double Glacier::GetRestrictedElevationChange() const
{
    double restrictedElevationChange = GetGlacierIce()->GetRestrictedElevationChange();
    return restrictedElevationChange;
}

void Glacier::SetGlacierIceThickness(double value)
{
    GetGlacierIce()->SetGlacierIceThickness(value);
}

double Glacier::GetGlacierIceThickness() const
{
    double glacierIceThickness = GetGlacierIce()->GetGlacierIceThickness();
    return glacierIceThickness;
}

void Glacier::SetGlacierIceVolume(double value)
{
    GetGlacierIce()->SetGlacierIceVolume(value);
}

double Glacier::GetGlacierIceVolume() const
{
    double glacierIceVolume = GetGlacierIce()->GetGlacierIceVolume();
    return glacierIceVolume;
}

void Glacier::SetThisYearAnnualGlacierValues(DateTime datetime)
{
    double annualMassBal = glacierMassBalance - GetPreviousYearMassBalance();
    //  cout << "2 SetThisYearAnnualGlacierValues " << "  glacierMassBalance " << glacierMassBalance << "  GetPreviousYearMassBalance() " << GetPreviousYearMassBalance() << "           annualMassBal "  << annualMassBal << endl;
    if (GetGlacierIceAreaFraction() > 0.0)
    {
        SetAnnualMassBalance(annualMassBal);
        SetPreviousYearMassBalance(glacierMassBalance);
    }
    else
    {
        SetAnnualMassBalance(missingData);
        SetPreviousYearMassBalance(0.0);
    }
    //  cout << "3 SetThisYearAnnualGlacierValues " << "  glacierMassBalance " << glacierMassBalance << "  GetPreviousYearMassBalance() " << GetPreviousYearMassBalance() << "  annualMassBalance "  << annualMassBalance << endl;
}

void Glacier::SnowToGlacierIce()
{
    //  cout << " Glacier SnowToGlacierIce \n";
    GetGlacierIce()->SnowToGlacierIce(GetSnow()->GetSnowStore());
}

void Glacier::RemoveSnowOnGlacierIce()
{
    //  cout << " Glacier RemoveSnowOnGlacierIce \n";
    if (GetGlacierIceAreaFraction() > 0.0)
    {
        GetSnow()->SetSnowStore(0.0);
    }
}

void Glacier::SetSnowStore(double value)
{
    GetSnow()->SetSnowStore(value);
}

void Glacier::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
    GetHBV()->SetSubSurfaceHbvStore(sm, uz, lz);
}

/*void Glacier::SetInitialKwaValues()
{
GetKinematicWave()->SetInitialKwaValues();
}*/

void Glacier::WaterBalance(int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge, int dayofyear)
{
    double dryPeriod, snowIceCoverFraction;
    GetGlacierIce()->WaterBalance(timeStep, datetime, GetGlacierIceAreaFraction());
    if (GetGlacierIceAreaFraction() == GetAreaFraction())
    {
        GetSnow()->WaterBalance(timeStep, datetime, GetGlacierIce()->GetPrecipitation(), GetGlacierIce()->GetTemperature());
        dryPeriod = missingData;
    }
    else
    {
        GetVegetation()->WaterBalance(timeStep, datetime, GetGlacierIce()->GetPrecipitation(), GetGlacierIce()->GetTemperature(), dayofyear, GetSnow()->GetSnowCoverFraction());
        GetSnow()->WaterBalance(timeStep, datetime, GetVegetation()->GetThroughFall() * (1.0 - (GetGlacierIceAreaFraction() / GetAreaFraction())), GetVegetation()->GetTemperature());
        dryPeriod = GetVegetation()->GetDryPeriod();
    }
    //  GetGlacierIce()->WaterBalance(timeStep, GetGlacierSurface()->GetWaterOutput());
    //  if (timeStep >= initialTimeSteps) {
    if (timeStep > initialTimeSteps)
    {
        glacierMassBalance = glacierMassBalance + GetSnow()->GetSnowWaterEquivalentChange() - GetGlacierIce()->GetIceMelt() * GetGlacierIceAreaFraction() / GetAreaFraction() * (1.0 - GetSnow()->GetSnowCoverFraction());
    }
    else
    {
        glacierMassBalance = 0.0;
    }
    // Algorithm to be performed in case: no input to landscape element from upstream elements
    /*  GetHBV()->WaterBalance(timeStep, GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction())
    + GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), missingData);*/
    // Algorithm to be performed in case: input to landscape element from upstream elements
    /*  if (GetGlacierIceAreaFraction() == GetAreaFraction())
    dryPeriod = missingData;
    else
    dryPeriod = 1.0;*/
    if (GetGlacierIceAreaFraction() / GetAreaFraction() > GetSnow()->GetSnowCoverFraction())
    {
        snowIceCoverFraction = GetGlacierIceAreaFraction() / GetAreaFraction();
    }
    else
    {
        snowIceCoverFraction = GetSnow()->GetSnowCoverFraction();
    }
    if (GetHBV())
    {
        //    GetHBV()->WaterBalance(timeStep, GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction())+GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), dryPeriod, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
        GetHBV()->WaterBalance(timeStep, datetime, GetGlacierIce()->GetIceMelt()*GetGlacierIceAreaFraction() / GetAreaFraction() * (1.0 - GetSnow()->GetSnowCoverFraction()) + GetSnow()->GetWaterOutput(), GetGlacierIce()->GetTemperature(), snowIceCoverFraction, dryPeriod, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge, dayofyear);
    }
    if (GetKinematicWave())
    {
        //    GetKinematicWave()->WaterBalance(timeStep, GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction())+GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), dryPeriod, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
        GetKinematicWave()->WaterBalance(timeStep, datetime, GetGlacierIce()->GetIceMelt()*GetGlacierIceAreaFraction() / GetAreaFraction() * (1.0 - GetSnow()->GetSnowCoverFraction()) + GetSnow()->GetWaterOutput(), GetGlacierIce()->GetTemperature(), snowIceCoverFraction, dryPeriod, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
    }
    //  cout << "  end GlacierWaterBalance\n";
}

double Glacier::GetPrecipitation() const
{
    double precipitation = GetGlacierIce()->GetPrecipitation();
    return precipitation * GetAreaFraction() / 100.0;
}

double Glacier::GetTemperature() const
{
    double temperature = GetGlacierIce()->GetTemperature();
    return temperature * GetAreaFraction() / 100.0;
}

double Glacier::GetSnowCoverFraction() const
{
    double snowCoverFraction = GetSnow()->GetSnowCoverFraction();
    return snowCoverFraction * GetAreaFraction() / 100.0;
}

double Glacier::GetSnowStore() const
{
    double snowStore = GetSnow()->GetSnowStore();
    return snowStore * GetAreaFraction() / 100.0;
}

double Glacier::GetMeltWater() const
{
    double meltWater = GetSnow()->GetMeltWater();
    return meltWater * GetAreaFraction() / 100.0;
}

double Glacier::GetSnowWaterEquivalentChange() const
{
    double snowWaterEquivalentChange = GetSnow()->GetSnowWaterEquivalentChange();
    return snowWaterEquivalentChange * GetAreaFraction() / 100.0;
}

double Glacier::GetWaterOutput() const
{
    //  double waterOutput=GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction())+GetSnow()->GetWaterOutput();
    double waterOutput = GetGlacierIce()->GetIceMelt() * GetGlacierIceAreaFraction() / GetAreaFraction() * (1.0 - GetSnow()->GetSnowCoverFraction()) + GetSnow()->GetWaterOutput();
    return waterOutput * GetAreaFraction() / 100.0;
}

double Glacier::GetGlacierMassBalance() const
{
    //  return glacierMassBalance*GetAreaFraction()/100.0;
    return glacierMassBalance * GetGlacierIceAreaFraction() / 100.0;
}

double Glacier::GetGlacierIceMelt() const
{
    double iceMelt = GetGlacierIce()->GetIceMelt() * GetGlacierIceAreaFraction() / GetAreaFraction() * (1.0 - GetSnow()->GetSnowCoverFraction());
    //  return iceMelt*GetAreaFraction()/100.0;
    return iceMelt * GetGlacierIceAreaFraction() / 100.0;
}

double Glacier::GetAnnualMassBalance() const
{
    if (annualMassBalance != missingData)
    {
        return annualMassBalance * GetGlacierIceAreaFraction() / 100.0;
    }
    else
    {
        return missingData;
    }
}

double Glacier::GetKiWaSoilMoisture(double lengthFraction) const
{
    double soilMoisture = GetKinematicWave()->GetSoilMoisture(lengthFraction);
    return soilMoisture * GetAreaFraction() / 100.0;
}

double Glacier::GetKiWaGroundWaterDepth(double lengthFraction) const
{
    double groundWaterDepth = GetKinematicWave()->GetGroundWaterDepth(lengthFraction);
    return groundWaterDepth * GetAreaFraction() / 100.0;
}

double Glacier::GetHbvSoilMoisture() const
{
    double soilMoisture = GetHBV()->GetSoilMoisture();
    return soilMoisture * GetAreaFraction() / 100.0;
}

double Glacier::GetHbvSoilMoistureDeficit() const
{
    double soilMoistureDeficit = GetHBV()->GetSoilMoistureDeficit();
    return soilMoistureDeficit * GetAreaFraction() / 100.0;
}

double Glacier::GetHbvPercSoilUpper() const
{
    double percSoilUpper = GetHBV()->GetPercSoilUpper();
    return percSoilUpper * GetAreaFraction() / 100.0;
}

double Glacier::GetHbvUpperZone() const
{
    double upperZone = GetHBV()->GetUpperZone();
    return upperZone * GetAreaFraction() / 100.0;
}

double Glacier::GetHbvLowerZone() const
{
    double lowerZone = GetHBV()->GetLowerZone();
    return lowerZone * GetAreaFraction() / 100.0;
}

double Glacier::GetInterceptionLoss() const
{
    double interceptionLoss = 0.0;
    interceptionLoss = GetVegetation()->GetInterceptionLoss() * (1.0 - (GetGlacierIceAreaFraction() / GetAreaFraction())) * GetAreaFraction() / 100.0;
    return interceptionLoss;
}

double Glacier::GetTranspSoilEvap() const
{
    double transpSoilEvap = 0.0;
    if (GetHBV())
    {
        transpSoilEvap = transpSoilEvap + GetHBV()->GetTranspSoilEvap();
    }
    if (GetKinematicWave())
    {
        transpSoilEvap = transpSoilEvap + GetKinematicWave()->GetTranspSoilEvap();
    }
    return transpSoilEvap * GetAreaFraction() / 100.0;
}

double Glacier::GetRunoff() const
{
    double runoff = 0.0;
    if (GetHBV())
    {
        runoff = runoff + GetHBV()->GetRunoff();
    }
    if (GetKinematicWave())
    {
        runoff = runoff + GetKinematicWave()->GetRunoff();
    }
    //  cout << "glacier runoff " << runoff << "  " << GetGlacierSurface()->GetIceMelt() + GetSnow()->GetWaterOutput() << endl;
    return runoff * GetAreaFraction() / 100.0;
}

double Glacier::GetLowerRunoff() const
{
    double lowerRunoff = 0.0;
    if (GetHBV())
    {
        lowerRunoff = lowerRunoff + GetHBV()->GetLowerRunoff();
    }
    if (GetKinematicWave())
    {
        lowerRunoff = lowerRunoff + GetKinematicWave()->GetLowerRunoff();
    }
    return lowerRunoff * GetAreaFraction() / 100.0;
}

double Glacier::GetUpperRunoff() const
{
    double upperRunoff = 0.0;
    if (GetHBV())
    {
        upperRunoff = upperRunoff + GetHBV()->GetUpperRunoff();
    }
    if (GetKinematicWave())
    {
        upperRunoff = upperRunoff + GetKinematicWave()->GetUpperRunoff();
    }
    return upperRunoff * GetAreaFraction() / 100.0;
}
