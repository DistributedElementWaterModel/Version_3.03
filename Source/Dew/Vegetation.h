#pragma once

using namespace std;
class DateTime;
class ParametersGeneral;
class ParametersLandSurface;
class InputElement;
class DistributedElement;
class ModelControl;

class Vegetation
{
public:
    Vegetation();
    ~Vegetation();
    void WaterBalance(int timeStep, DateTime datetime, double inputPrecipitation, double inputTemperature);
    double GetPrecipitation() const;
    double GetTemperature() const;
    double GetInterceptionLoss() const;
    double GetThroughFall() const;
    void SetInterceptionStore(double value);
    double GetInterceptionStore() const;
    void SetLeafAreaIndex(double value);
    double GetLeafAreaIndex() const;
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
    double GetDryPeriod();

private:
    ParametersGeneral *commonPar;
    ParametersLandSurface *landSurfacePar;
    //  InputTimeSeries *inTimeSeries;
    InputElement *inElement;
    DistributedElement *landScapeElement;
    double precipitation;                    /*  Precipitation (m/timestep)  */
    double temp;                             /*  Air temperature (deg. C)  */
    double potev;                            /*  Potential evapotranspiration (m/timestep)  */
    double prevInterception;                 /*  Interception from previous time step (m)  */
    double interceptionStore;                /*  Interception store (m)  */
    double interceptionLoss;                 /*  Interception loss from vegetation (m)  */
    double throughFall;                      /*  Throughfall (m/timestep)  */
    double wetPeriod, dryPeriod;             /*  Length of wet and dry period during evapotranspiration (fraction of timestep)  */
    double leafAreaIndex;                    /*  Leaf area index  */
};
