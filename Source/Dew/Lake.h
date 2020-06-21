#pragma once

#include "Dew.h"

using namespace std;
class DateTime;
class LakeWaterBalance;
class DistributedElement;
//class LandAtmosphereInterface;
//class Albedo;
//class Atmosphere;
//class PenmanMonteith;

class Lake
{
public:
    Lake();
    ~Lake();
    void SetLakeWaterBalance(LakeWaterBalance *theLakeWaterBalance);
    LakeWaterBalance *GetLakeWaterBalance() const;
    void SetAreaFraction(double value);
    double GetAreaFraction() const;
    void SetLandSurfaceType(LANDSURFACE landSurf);
    LANDSURFACE GetLandSurfaceType() const;
    void SetSoilType(SOIL soil);
    SOIL GetSoilType() const;
    void WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedDischarge);
    //    void SetSurfaceAlbedo(double value);
    //    double GetSurfaceAlbedo() const;
    double GetPrecipitation() const;
    double GetTemperature() const;
    double GetLakeEvap() const;
    double GetRunoff() const;
    double GetLakeStorage() const;
    double GetLakeStorageChange() const;
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
    LakeWaterBalance *ptrLakeWaterBalance;
    DistributedElement *landScapeElement;
    //    LandAtmosphereInterface * landAtmosphereInterfaceElement;
    //    Albedo * albedoElement;
    //    Atmosphere * atmosphereElement;
    //    PenmanMonteith * penmanMonteithElement;
    double areaFraction;
    LANDSURFACE landSurfaceType;
    SOIL soilType;
};
