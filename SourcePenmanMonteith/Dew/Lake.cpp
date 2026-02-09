#include "Lake.h"
#include "LakeWaterBalance.h"
#include "DistributedElement.h"
#include "InputElement.h"
#include "ModelControl.h"
#include "ParametersGeneral.h"
#include "stdafx.h"
#include "DateTime.h"

Lake::Lake()
    //    Lake::Lake() :
    //    surfaceAlbedo(0.0)
{
    SetLandScapeElement(0);
    //    SetLandAtmosphereInterfaceElement(0);
    //    SetAlbedoElement(0);
    //    SetAtmosphereElement(0);
    //    SetPenmanMonteithElement(0);
    SetAreaFraction(0.0);
}

Lake::~Lake()
{
}

void  Lake::SetLakeWaterBalance(LakeWaterBalance *theLakeWaterBalance)
{
    ptrLakeWaterBalance = theLakeWaterBalance;
}
LakeWaterBalance * Lake::GetLakeWaterBalance() const
{
    return ptrLakeWaterBalance;
}
void  Lake::SetAreaFraction(double value)
{
    areaFraction = value;
}
double  Lake::GetAreaFraction() const
{
    return areaFraction;
}
void  Lake::SetLandSurfaceType(LANDSURFACE landSurf)
{
    landSurfaceType = landSurf;
}
LANDSURFACE  Lake::GetLandSurfaceType() const
{
    return landSurfaceType;
}
void  Lake::SetSoilType(SOIL soil)
{
    soilType = soil;
}
SOIL  Lake::GetSoilType() const
{
    return soilType;
}
/*void  Lake::SetSurfaceAlbedo(double value)
{
    surfaceAlbedo = value;
}
double  Lake::GetSurfaceAlbedo() const
{
    return surfaceAlbedo;
    }*/
void  Lake::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * Lake::GetLandScapeElement() const
{
    return landScapeElement;
}
/*void  Lake::SetLandAtmosphereInterfaceElement(LandAtmosphereInterface *theElement)
{
    landAtmosphereInterfaceElement = theElement;
}
LandAtmosphereInterface * Lake::GetLandAtmosphereInterfaceElement() const
{
    return landAtmosphereInterfaceElement;
}
void  Lake::SetAlbedoElement(Albedo *theElement)
{
    albedoElement = theElement;
}
Albedo * Lake::GetAlbedoElement() const
{
    return albedoElement;
}
void  Lake::SetAtmosphereElement(Atmosphere *theElement)
{
    atmosphereElement = theElement;
}
Atmosphere * Lake::GetAtmosphereElement() const
{
    return atmosphereElement;
}
void  Lake::SetPenmanMonteithElement(PenmanMonteith *theElement)
{
    penmanMonteithElement = theElement;
}
PenmanMonteith * Lake::GetPenmanMonteithElement() const
{
    return penmanMonteithElement;
}*/

void Lake::WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedDischarge, int dayofyear_PM2)
{
    //  GetAlbedoElement()->SurfaceAlbedo(this, 0, 0, 0);
    GetLakeWaterBalance()->WaterBalance(timeStep, datetime, upLandAccumulatedDischarge, dayofyear_PM2);
}

double Lake::GetPrecipitation() const
{
    double precipitation = GetLakeWaterBalance()->GetPrecipitation();
    return precipitation * GetAreaFraction() / 100.0;
}

double Lake::GetTemperature() const
{
    double temperature = GetLakeWaterBalance()->GetTemperature();
    return temperature * GetAreaFraction() / 100.0;
}

double Lake::GetLakeEvap() const
{
    double lakeEvap = GetLakeWaterBalance()->GetLakeEvap();
    return lakeEvap * GetAreaFraction() / 100.0;
}

double Lake::GetRunoff() const
{
    double runoff = GetLakeWaterBalance()->GetRunoff();
    return runoff * GetAreaFraction() / 100.0;
}

double Lake::GetLakeStorage() const
{
    double lakeStorage = GetLakeWaterBalance()->GetWaterLevel();
    return lakeStorage * GetAreaFraction() / 100.0;
}

double Lake::GetLakeStorageChange() const
{
    double storageChange = GetLakeWaterBalance()->GetWaterLevelChange();
    return storageChange * GetAreaFraction() / 100.0;
}
