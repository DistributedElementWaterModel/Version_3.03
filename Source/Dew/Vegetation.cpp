#include "Vegetation.h"
#include "DistributedElement.h"
#include "InputElement.h"
#include "ModelControl.h"
#include "ParametersLandSurface.h"
#include "Util.h"
#include "stdafx.h"
#include "DateTime.h"


Vegetation::Vegetation() :
    precipitation(0.0),
    temp(0.0),
    potev(0.0),
    prevInterception(0.0),
    interceptionStore(0.0),
    interceptionLoss(0.0),
    throughFall(0.0),
    wetPeriod(0.0),
    dryPeriod(0.0),
    leafAreaIndex(0.0)
{
    SetGeneralPar(0);
    SetLandSurfacePar(0);
    //  SetInputTimeSeries(0);
    SetInputElement(0);
    SetLandScapeElement(0);
}

Vegetation::~Vegetation()
{
}

void  Vegetation::SetInterceptionStore(double value)
{
    interceptionStore = value;
}
double  Vegetation::GetInterceptionStore() const
{
    return interceptionStore;
}
void  Vegetation::SetLeafAreaIndex(double value)
{
    leafAreaIndex = value;
}
double  Vegetation::GetLeafAreaIndex() const
{
    return leafAreaIndex;
}
void  Vegetation::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral * Vegetation::GetGeneralPar() const
{
    return commonPar;
}
void  Vegetation::SetLandSurfacePar(ParametersLandSurface *parObj)
{
    landSurfacePar = parObj;
}
ParametersLandSurface * Vegetation::GetLandSurfacePar()
{
    return landSurfacePar;
}
//  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
//  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
void  Vegetation::SetInputElement(InputElement *inElementObj)
{
    inElement = inElementObj;
}
InputElement * Vegetation::GetInputElement() const
{
    return inElement;
}
void  Vegetation::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * Vegetation::GetLandScapeElement() const
{
    return landScapeElement;
}
double  Vegetation::GetDryPeriod()
{
    return dryPeriod;
}

double Vegetation::GetPrecipitation() const
{
    return precipitation;
}

double Vegetation::GetTemperature() const
{
    return temp;
}

double Vegetation::GetInterceptionLoss() const
{
    return interceptionLoss;
}

double Vegetation::GetThroughFall() const
{
    return throughFall;
}


void Vegetation::WaterBalance(int timeStep, DateTime datetime, double inputPrecipitation, double inputTemperature)
{
    double timeResolution = 1.0;
    int i;
    /*  if (inputPrecipitation > missingData && inputTemperature > missingData) {
    precipitation = inputPrecipitation;
    temp = inputTemperature;
    }
    else {
    precipitation=GetInputElement()->GetInput(0);
    temp = GetInputElement()->GetInput(1);
    }*/
    if (inputPrecipitation > missingData)
    {
        precipitation = inputPrecipitation;
    }
    else
    {
        precipitation = GetInputElement()->GetInput(0);
    }
    if (inputTemperature > missingData)
    {
        temp = inputTemperature;
    }
    else
    {
        temp = GetInputElement()->GetInput(1);
    }
    //  cout << "    timeStep " << timeStep << "    Vegetation precipitation " << precipitation;
    //  cout << "    Temperature " << temp << endl;

    /* Potential evapotranspiration */
    wetPeriod = 0.0;
    dryPeriod = timeResolution;
    throughFall = 0.0;
    interceptionLoss = 0.0;
    if (GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'M' || 
        GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'm') {
      i = datetime.getMonth()-1;
      potev = potentialEvapLongTermMean(temp, GetLandScapeElement()->GetModelControlObj()->GetEvaporationArray(i)/1000.0);
      /*for (i=0; i<numberPotentialEvaporationValuesPerYear; i++) {
        cout << " " << i << "    " << GetLandScapeElement()->GetModelControlObj()->GetEvaporationArray(i) << "\n";
        }*/
    }
    else {
	potev = potentialEvapTemperatureIndex(temp, landSurfacePar->GetEPOT_PAR());
    }

    /*    if (GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'T' || GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 't')
	  {
	  potev = potentialEvapTemperatureIndex(temp, landSurfacePar->GetEPOT_PAR());
	  }*/

    /* Water input from precipitation and/or snowmelt > 0 or interception remaining from previous time step */
    if (precipitation > 0.0 || prevInterception > 0.0)
    {
        interceptionStore = prevInterception + precipitation;
        if (potev > 0.0)
        {
            wetPeriod = interceptionStore / potev;
        }
        else
        {
            wetPeriod = 0.0;
        }
        //    cout << interceptionStore << " " << potev  << " " << timeResolution  << " " << wetPeriod << " " << dryPeriod << endl;
        if (wetPeriod < timeResolution)
        {
            interceptionLoss = potev * wetPeriod;
            dryPeriod = timeResolution - (wetPeriod * landSurfacePar->GetWET_PER_CORR());
        }
        else
        {
            interceptionLoss = potev * timeResolution;
            dryPeriod = 0.0;
        }
        interceptionStore = interceptionStore - interceptionLoss;

        /* If interception store > INTER_MAX, surplus water is infiltrated through the soil surface
        Soil moistured deficit and percolation is calculated */
        //    cout << interceptionStore << " " << potev  << " " << timeResolution  << " " << wetPeriod << " " << dryPeriod << endl;
        if (interceptionStore > landSurfacePar->GetINTER_MAX())
        {
            throughFall = interceptionStore - landSurfacePar->GetINTER_MAX();
            interceptionStore = landSurfacePar->GetINTER_MAX();
        }

        /* If interception store < 0 following evaporation at potential rate then dry period > 0,
        transpiration and soil evaporation will be calculated for dryPeriod */
        if (interceptionStore < 0.0 - epsilon * epsilon)
        {
            printf("    timeStep = %d      interceptionStore = %f\n", timeStep, interceptionStore);
            /*        exit(1);*/
        }
        if (wetPeriod < 0.0)
        {
            printf("    timeStep = %d      wetPeriod = %f\n", timeStep, wetPeriod);
            /*        exit(1);*/
        }
        if (dryPeriod < 0.0 || dryPeriod > timeResolution)
        {
            printf("    timeStep = %d      dryPeriod = %f\n", timeStep, dryPeriod);
            /*        exit(1);*/
        }

        /* If 0 <= interception store <= INTER_MAX, no action is necessary
        Infiltration through soil surface = 0, soil moisture deficit and depth of saturated zone is unchanged */
    }
    prevInterception = interceptionStore;
    //  cout << "throughFall " << throughFall << endl;
}

