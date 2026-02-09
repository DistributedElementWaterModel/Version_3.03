#pragma once

using namespace std;
class DateTime;
class ParametersGeneral;
class ParametersLandSurface;
class ParametersGlacierRetreat;
class InputElement;
class DistributedElement;

class GlacierIce
{
public:
    GlacierIce();
    ~GlacierIce();
    void WaterBalance(int timeStep, DateTime datetime, double glacierIceAreaFraction);
    double GetPrecipitation() const;
    double GetTemperature() const;
    double GetIceMelt() const;
    void SetGlacierIceThickness(double value);
    double GetGlacierIceThickness() const;
    void SetGlacierIceVolume(double value);
    double GetGlacierIceVolume() const;
    void SetGlacierSurfaceElevation(double value);
    double GetGlacierSurfaceElevation() const;
    void SetGlacierSurfaceElevationChangeNormalized(double value);
    double GetGlacierSurfaceElevationChangeNormalized() const;
    void SetRestrictedElevationChange(double value);
    double GetRestrictedElevationChange();
    void SnowToGlacierIce(double snowStore);
    void SetGeneralPar(ParametersGeneral *parObj);
    ParametersGeneral *GetGeneralPar() const;
    void SetLandSurfacePar(ParametersLandSurface *parObj);
    ParametersLandSurface *GetLandSurfacePar();
    void SetGlacierRetreatPar(ParametersGlacierRetreat *parObj);
    ParametersGlacierRetreat *GetGlacierRetreatPar() const;
    //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj);
    //  InputTimeSeries *GetInputTimeSeries() const;
    void SetInputElement(InputElement *inElementObj);
    InputElement *GetInputElement() const;
    void SetLandScapeElement(DistributedElement *theElement);
    DistributedElement *GetLandScapeElement() const;

private:
    ParametersGeneral *commonPar;
    ParametersLandSurface *landSurfacePar;
    ParametersGlacierRetreat *glacRetPar;
    //  InputTimeSeries *inTimeSeries;
    InputElement *inElement;
    DistributedElement *landScapeElement;
    double precipitation;                    /*  Precipitation (m/timestep)  */
    double temp;                             /*  Air temperature (deg. C)  */
    double iceMelt;                          /*  Meltwater from ice (m/timestep)  */
    double glacierIceThickness;              /*  Thickness of glacier ice  */
    double glacierIceVolume;                 /*  Volume of glacier ice  */
    double glacierSurfaceElevation;          /*  Elevation of glacier surface  */
    double glacierSurfaceElevationChangeNormalized;    /*  Elevation change of glacier surface  */
    double restrictedGlacierSurfaceElevationChange;    /*  Elevation change restricted by glacier tongue mass balance  */
};
