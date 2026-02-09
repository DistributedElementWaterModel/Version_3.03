#pragma once

using namespace std;
class DateTime;
class ParametersGeneral;
class ParametersLandSurface;
class ParametersSubSurfaceHbv;
class InputElement;
class DistributedElement;
class ModelControl;

class HBV
{
public:
    HBV();
    ~HBV();
    void SetInitialHbvValues();
    void SetSubSurfaceHbvStore(double sm, double uz, double lz);
    // Algorithm to be performed in case: no input to landscape element from upstream elements
    //  void WaterBalance(int timeStep, double waterInput, double snowCoverFraction, double dryPeriod);
    // Algorithm to be performed in case: input to landscape element from upstream elements
    void WaterBalance(int timeStep, DateTime datetime, double waterInput, double temperature, double snowCoverFraction, double dryPeriod, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge, double potev);
    double GetSoilMoisture() const;
    double GetSoilMoistureDeficit() const;
    double GetPercSoilUpper() const;
    double GetUpperZone() const;
    double GetLowerZone() const;
    double GetTranspSoilEvap() const;
    double GetLowerRunoff() const;
    double GetUpperRunoff() const;
    double GetRunoff() const;
    void SetGeneralPar(ParametersGeneral *parObj) ;
    ParametersGeneral *GetGeneralPar() const ;
    void SetLandSurfacePar(ParametersLandSurface *parObj) ;
    ParametersLandSurface *GetLandSurfacePar() ;
    void SetSubSurfaceHbvPar(ParametersSubSurfaceHbv *parObj) ;
    ParametersSubSurfaceHbv *GetSubSurfaceHbvPar() const ;
    //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) ;
    //  InputTimeSeries *GetInputTimeSeries() const ;
    void SetInputElement(InputElement *inElementObj) ;
    InputElement *GetInputElement() const ;
    void SetLandScapeElement(DistributedElement *theElement) ;
    DistributedElement *GetLandScapeElement() const ;

private:
    ParametersGeneral *commonPar;
    ParametersSubSurfaceHbv *subSurfacePar;
    ParametersLandSurface *landSurfacePar;
    //  InputTimeSeries *inTimeSeries;
    InputElement *inElement;
    DistributedElement *landScapeElement;
    //  double temp;                                /*  Temperature (deg. C)  */
    double soilMoisture;                        /*  Soil moisture content (m)  */
    double percSoilUpper;                       /*  Percolation from soil moisture zone to upper zone (m/timestep)  */
    double upperZone;                           /*  Upper groundwater zone water content (m)  */
    double lowerZone;                           /*  Lower groundwater zone water content (m)  */
    double lowerRunoff;                         /*  Runoff from lower layer (m/timestep)  */
    double upperRunoff;                         /*  Runoff from upper layer (m/timestep)  */
    double transpSoilEvap;                      /*  Water lost from subsurface by evapotranspiration (m)  */
    double runoff;                              /*  Runoff (m/timestep)  */
};
