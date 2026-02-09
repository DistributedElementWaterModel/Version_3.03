#pragma once

class Snow;
class GlacierIce;
class HBV;
class KinematicWave;
class Vegetation;
class DistributedElement;
//class LandAtmosphereInterface;
//class Albedo;
//class Atmosphere;
//class PenmanMonteith;
class DateTime;
#include "Dew.h"

class Glacier
{
public:
    Glacier(DateTime startModelTime, DateTime startSimulationTime, DateTime endSimulationTime);
    ~Glacier();
    void SetSnow(Snow *theSnow);
    Snow *GetSnow() const;
    //  void SetGlacierSurface(GlacierSurface *theGlacierSurface);
    //  GlacierSurface *GetGlacierSurface() const;
    void SetGlacierIce(GlacierIce *theGlacierIce);
    GlacierIce *GetGlacierIce() const;
    void SetHBV(HBV *theHBV);
    HBV *GetHBV() const;
    void SetKinematicWave(KinematicWave *theKinematicWave);
    KinematicWave *GetKinematicWave() const;
    void SetVegetation(Vegetation *theVegetation);
    Vegetation *GetVegetation() const;
    void SetAreaFraction(double value);
    double GetAreaFraction() const;
    void SetLandSurfaceType(LANDSURFACE landSurf);
    LANDSURFACE GetLandSurfaceType() const;
    void SetSoilType(SOIL soil);
    SOIL GetSoilType() const;
    void SetGlacierIceAreaFraction(double value);
    double GetGlacierIceAreaFraction() const;
    void SetGlacierSurfaceElevation(double value);
    double GetGlacierSurfaceElevation() const;
    void SetGlacierSurfaceElevationChangeNormalized(double value);
    double GetGlacierSurfaceElevationChangeNormalized() const;
    void SetRestrictedElevationChange(double value);
    double GetRestrictedElevationChange() const;
    void SetGlacierIceThickness(double value);
    double GetGlacierIceThickness() const;
    void SetGlacierIceVolume(double value);
    double GetGlacierIceVolume() const;
    void SnowToGlacierIce();
    void RemoveSnowOnGlacierIce();
    void SetSnowStore(double value);
    void SetSubSurfaceHbvStore(double sm, double uz, double lz);
    //  void SetInitialKwaValues();
    void WaterBalance(int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge, int dayofyear);
    //    void SetSurfaceAlbedo(double value);
    //    double GetSurfaceAlbedo() const;
    void SetPreviousYearMassBalance(double value);
    double GetPreviousYearMassBalance() const;
    void SetThisYearAnnualGlacierValues(DateTime datetime);
    void SetAnnualMassBalance(double value);
    double GetAnnualMassBalance() const;
    double GetPrecipitation() const;
    double GetTemperature() const;
    double GetSnowCoverFraction() const;
    double GetSnowStore() const;
    double GetMeltWater() const;
    double GetSnowWaterEquivalentChange() const;
    double GetWaterOutput() const;
    double GetGlacierMassBalance() const;
    double GetGlacierIceMelt() const;
    double GetKiWaSoilMoisture(double lengthFraction) const;
    double GetKiWaGroundWaterDepth(double lengthFraction) const;
    double GetHbvSoilMoisture() const;
    double GetHbvSoilMoistureDeficit() const;
    double GetHbvPercSoilUpper() const;
    double GetHbvUpperZone() const;
    double GetHbvLowerZone() const;
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
    Snow *ptrSnow;
    //  GlacierSurface *ptrSurface;
    GlacierIce *ptrIce;
    HBV *ptrHbv;
    KinematicWave *ptrKwa;
    Vegetation *ptrVeg;
    DistributedElement *landScapeElement;
    //    LandAtmosphereInterface * landAtmosphereInterfaceElement;
    //    Albedo * albedoElement;
    //    Atmosphere * atmosphereElement;
    //    PenmanMonteith * penmanMonteithElement;
    double areaFraction;
    double glacierIceAreaFraction;
    double glacierMassBalance;
    double previousYearMassBalance;
    double annualMassBalance;
    LANDSURFACE landSurfaceType;
    SOIL soilType;
};

