#pragma once

using namespace std;
class DateTime;
class ParametersGeneral;
class ParametersLandSurface;
class InputElement;
class DistributedElement;

#include "Dew.h"

class Snow
{
public:
    Snow();
    ~Snow();
    void WaterBalance(int timeStep, DateTime datetime, double waterInput, double temperature);
    //    void SetSnowAlbedo(double value);
    //    double GetSnowAlbedo() const;
    void SetSnowStore(double value);
    double GetSnowStore() const;
    double GetMeltWater() const;
    double GetSnowWaterEquivalentChange() const;
    double GetWaterOutput() const;
    double GetSnowCoverFraction() const;
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
    //  double temp;                               /*  Air temperature (deg. C)  */
    //    double snowAlbedo;                      /*  Albedo of snow surface  */
    double snowStore;                          /*  Snow store (m)  */
    double meltWater;                          /*  Meltwater in snow (m)  */
    double snowWaterEquivalentChange;          /*  Change of snow water equivalent (m)  */
    double waterOutput;                        /*  Output of meltwater from snow store (m/timestep)  */
    double snowCoverFraction;                  /*  Fraction of area covered by snow  */
    double distSnowStore[numberSnowClasses];   /*  Distributed snow store (m)  */
    double distMeltWater[numberSnowClasses];   /*  Distributed meltwater in snow (m)  */
    double distWaterOutput[numberSnowClasses]; /*  Output of meltwater from snow store (m/timestep)  */
};
