#pragma once

#include "Dew.h"

using namespace std;
class DateTime;
class Vegetation;
class Snow;
class KinematicWave;
class DistributedElement;
//class LandAtmosphereInterface;
//class Albedo;
//class Atmosphere;
//class PenmanMonteith;

class KwaAquifer
{
public:
    KwaAquifer();
    ~KwaAquifer();
    void SetNextKwaAquifer(KwaAquifer *theKwaAquifer);
    KwaAquifer *GetNextKwaAquifer() const;
    void SetVegetation(Vegetation *theVegetation);
    Vegetation *GetVegetation() const;
    void SetSnow(Snow *theSnow);
    Snow *GetSnow() const;
    void SetKinematicWave(KinematicWave *theKinematicWave);
    KinematicWave *GetKinematicWave() const;
    void SetAreaFraction(double value);
    double GetAreaFraction() const;
    void SetLandSurfaceType(LANDSURFACE landSurf);
    LANDSURFACE GetLandSurfaceType() const;
    void SetSoilType(SOIL soil);
    SOIL GetSoilType() const;
    void SetSnowStore(double value);
    void WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge, int dayofyear);
    //    void SetSurfaceAlbedo(double value);
    //    double GetSurfaceAlbedo() const;
    double GetTotalKwaAreaFraction() const;
    double GetPrecipitation() const;
    double GetTemperature() const;
    double GetSnowCoverFraction() const;
    double GetSnowStore() const;
    double GetMeltWater() const;
    double GetSnowWaterEquivalentChange() const;
    double GetWaterOutput() const;
    double GetSoilMoisture(double lengthFraction) const;
    double GetGroundWaterDepth(double lengthFraction) const;
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
    KinematicWave *ptrKwa;
    KwaAquifer *nextKwaAquifer;
    DistributedElement *landScapeElement;
    //    LandAtmosphereInterface * landAtmosphereInterfaceElement;
    //    Albedo * albedoElement;
    //    Atmosphere * atmosphereElement;
    //    PenmanMonteith * penmanMonteithElement;
    double areaFraction;
    LANDSURFACE landSurfaceType;
    SOIL soilType;
};
