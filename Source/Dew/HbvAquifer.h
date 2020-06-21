#pragma once

using namespace std;
class DateTime;
class Vegetation;
class Snow;
class HBV;
class DistributedElement;
//class LandAtmosphereInterface;
//class Albedo;
//class Atmosphere;
//class PenmanMonteith;

#include "Dew.h"

class HbvAquifer
{
public:
    HbvAquifer();
    ~HbvAquifer();
    void SetNextHbvAquifer(HbvAquifer *theHbvAquifer);
    HbvAquifer *GetNextHbvAquifer() const;
    void SetVegetation(Vegetation *theVegetation);
    Vegetation *GetVegetation() const;
    void SetSnow(Snow *theSnow);
    Snow *GetSnow() const;
    void SetHBV(HBV *theHBV);
    HBV *GetHBV() const;
    void SetAreaFraction(double value);
    double GetAreaFraction() const;
    void SetLandSurfaceType(LANDSURFACE landSurf);
    LANDSURFACE GetLandSurfaceType() const;
    void SetSoilType(SOIL soil);
    SOIL GetSoilType() const;
    void SetSnowStore(double value);
    void SetSubSurfaceHbvStore(double sm, double uz, double lz);
    // Algorithm to be performed in case: no input to landscape element from upstream elements
    //  void WaterBalance(int timeStep);
    // Algorithm to be performed in case: input to landscape element from upstream elements
    void WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge);
    //    void SetSurfaceAlbedo(double value);
    //    double GetSurfaceAlbedo() const;
    double GetTotalHbvAreaFraction() const;
    double GetPrecipitation() const;
    double GetTemperature() const;
    double GetSnowCoverFraction() const;
    double GetSnowStore() const;
    double GetMeltWater() const;
    double GetSnowWaterEquivalentChange() const;
    double GetWaterOutput() const;
    double GetSoilMoisture() const;
    double GetSoilMoistureDeficit() const;
    double GetPercSoilUpper() const;
    double GetUpperZone() const;
    double GetLowerZone() const;
    double GetInterceptionLoss() const;
    double GetTranspSoilEvap() const;
    double GetLowerRunoff() const;
    double GetUpperRunoff() const;
    double GetRunoff() const;
    void SetLandScapeElement(DistributedElement *theElement);
    DistributedElement *GetLandScapeElement() const;
    //    void SetLandAtmosphereInterfaceElement(LandAtmosphereInterface *theElement);
    //    LandAtmosphereInterface *GetLandAtmosphereInterfaceElement() const;
    //    void SetAlbedoElement(Albedo *theElement);
    //    Albedo *GetAlbedoElement() const;
    //    void SetAtmosphereElement(Atmosphere *theElement);
    //    Atmosphere *GetAtmosphereElement() const;
    //    void SetPenmanMonteithElement(PenmanMonteith *theElement);
    //    PenmanMonteith *GetPenmanMonteithElement() const;

private:
    //    double surfaceAlbedo;                      /*  Albedo of land surface  */
    Vegetation *ptrVeg;
    Snow *ptrSnow;
    HBV *ptrHbv;
    HbvAquifer *nextHbvAquifer;
    DistributedElement *landScapeElement;
    //    LandAtmosphereInterface * landAtmosphereInterfaceElement;
    //    Albedo * albedoElement;
    //    Atmosphere * atmosphereElement;
    //    PenmanMonteith * penmanMonteithElement;
    double areaFraction;
    LANDSURFACE landSurfaceType;
    SOIL soilType;
};
