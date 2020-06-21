#pragma once

using namespace std;
class DateTime;
class ParametersGeneral;
class ParametersLandSurface;
class InputElement;
class DistributedElement;
class ModelControl;

class LakeWaterBalance
{
public:
    LakeWaterBalance();
    ~LakeWaterBalance();
    void WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedDischarge);
    void SetLakeValues(double temperature, double level);
    double GetPrecipitation() const;
    double GetTemperature() const;
    double GetLakeEvap() const;
    double GetRunoff() const;
    double GetWaterLevel() const;
    double GetWaterLevelChange() const;
    void SetGeneralPar(ParametersGeneral *parObj);
    ParametersGeneral *GetGeneralPar() const;
    void SetLandSurfacePar(ParametersLandSurface *parObj);
    ParametersLandSurface *GetLandSurfacePar();
    //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj);
    //  InputTimeSeries *GetInputTimeSeries() const;
    void SetInputElement(InputElement *inElementObj);
    InputElement *GetInputElement() const;
    void SetLandScapeElement(DistributedElement *theElement);
    DistributedElement *GetLandScapeElement() const;

private:
    ParametersGeneral *commonPar;
    ParametersLandSurface *landSurfacePar;
    //  InputTimeSeries *inTimeSeries;
    InputElement *inElement;
    DistributedElement *landScapeElement;
    double precipitation;                       /*  Precipitation (m/timestep)  */
    double temp;                                /*  Air temperature (deg. C)   */
    double lakeTemp;                            /*  Lake temperature (deg. C)  */
    double lakeEvaporation;                     /*  Evaporation (m/timestep)  */
    double waterLevel;                          /*  Lake water level (m)  */
    double waterLevelChange;                    /*  Lake water level change (m)  */
    double runoff;                              /*  Runoff (m/timestep)  */
    double discharge;                           /*  Lake outflow (m3/s)  */
};
