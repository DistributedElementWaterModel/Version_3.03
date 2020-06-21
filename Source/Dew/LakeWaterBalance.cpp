#include "LakeWaterBalance.h"
#include "Lake.h"
#include "DistributedElement.h"
#include "InputElement.h"
#include "ModelControl.h"
#include "ParametersGeneral.h"
#include "Util.h"
#include "stdafx.h"
#include "DateTime.h"


LakeWaterBalance::LakeWaterBalance() :
    precipitation(0.0),
    temp(0.0),
    lakeEvaporation(0.0),
    waterLevel(0.0),
    waterLevelChange(0.0),
    runoff(0.0),
    discharge(0.0)
{
    SetGeneralPar(0);
    SetLandSurfacePar(0);
    //  SetInputTimeSeries(0);
    SetInputElement(0);
    SetLandScapeElement(0);
}

LakeWaterBalance::~LakeWaterBalance()
{
}

void  LakeWaterBalance::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral * LakeWaterBalance::GetGeneralPar() const
{
    return commonPar;
}
void  LakeWaterBalance::SetLandSurfacePar(ParametersLandSurface *parObj)
{
    landSurfacePar = parObj;
}
ParametersLandSurface * LakeWaterBalance::GetLandSurfacePar()
{
    return landSurfacePar;
}
//  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
//  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
void  LakeWaterBalance::SetInputElement(InputElement *inElementObj)
{
    inElement = inElementObj;
}
InputElement * LakeWaterBalance::GetInputElement() const
{
    return inElement;
}
void  LakeWaterBalance::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * LakeWaterBalance::GetLandScapeElement() const
{
    return landScapeElement;
}

void LakeWaterBalance::SetLakeValues(double temperature, double level)
{
    lakeTemp = temperature;
    waterLevel = level;
}

double LakeWaterBalance::GetPrecipitation() const
{
    return precipitation;
}

double LakeWaterBalance::GetTemperature() const
{
    return temp;
}

double LakeWaterBalance::GetLakeEvap() const
{
    return lakeEvaporation;
}

double LakeWaterBalance::GetRunoff() const
{
    return runoff;
}

double LakeWaterBalance::GetWaterLevel() const
{
    return waterLevel;
}

double LakeWaterBalance::GetWaterLevelChange() const
{
    return waterLevelChange;
}


// ** Algorithm to be performed in case: no input to landscape element from upstream elements
//void LakeWaterBalance::WaterBalance(int timeStep)
// ** Algorithm to be performed in case: input to landscape element from upstream elements
void LakeWaterBalance::WaterBalance(int timeStep, DateTime datetime, double upLandAccumulatedDischarge)
{
    double tempMemory;
    double stage, newStage;
    double kLake;
    double deltaLevel;
    double nLake;
    double maximumLevel;
    double waterInput;
    int i;
    precipitation = GetInputElement()->GetInput(0);
    temp = GetInputElement()->GetInput(1);
    /*  cout << "    timeStep " << timeStep << "          Lake precipitation " << precipitation;
    cout << "    Temperature " << temp << endl;*/

    // ** Algorithm to be performed in case: no input to landscape element from upstream elements
    //  waterInput = precipitation;
    // ** Algorithm to be performed in case: input to landscape element from upstream elements
    /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
    waterInput = precipitation + upLandAccumulatedDischarge * commonPar->GetSECONDS_TIMESTEP() / GetLandScapeElement()->GetArea();

    /* Maximum lake storage */
    maximumLevel = commonPar->GetMAXIMUM_LEVEL();

    /*  Lake temperature and evaporation  */
    tempMemory = commonPar->GetDAY_TEMP_MEMORY() * 86400.0/commonPar->GetSECONDS_TIMESTEP();
    lakeTemp = lakeTemp*(1.0-(1.0/tempMemory))+temp/tempMemory;
    if (GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'M' || 
        GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'm') {
      i = datetime.getMonth()-1;
      lakeEvaporation = potentialEvapLongTermMean(lakeTemp, GetLandScapeElement()->GetModelControlObj()->GetEvaporationArray(i)/1000.0);
      /*for (i=0; i<numberPotentialEvaporationValuesPerYear; i++) {
        cout << " " << i << "    " << GetLandScapeElement()->GetModelControlObj()->GetEvaporationArray(i) << "\n";
        }*/
    }
    else {
      lakeEvaporation = potentialEvapTemperatureIndex(lakeTemp, commonPar->GetLAKE_EPOT_PAR());
    }
    //  else {
    //    PenmanMonteith();
    //  }

    /*  Lake rating curve parameters  */
    if (GetLandScapeElement()->GetLakeNumber() < 0)
    {
        deltaLevel = commonPar->GetDELTA_LEVEL();
        nLake = commonPar->GetNLAKE();
        kLake = commonPar->GetKLAKE();
    }
    else
    {
        deltaLevel = 0.0;
        nLake = 1.0;
        kLake = 1.0;
    }

    /*  Lake water level minimum  */
    if (waterLevel < (-1)*deltaLevel) {
      lakeEvaporation = 0.0;
    }

    /*  Initial lake water level  */
    stage = waterLevel + waterInput - lakeEvaporation;

    /* Water exceeding lake maximum level to runoff */
    if (stage > maximumLevel)
    {
        runoff = stage - maximumLevel;
        stage = maximumLevel;
    }
    else
    {
        runoff = 0.0;
    }

    /*  Lake outflow  */
    if (stage < (-1)*deltaLevel)
    {
        runoff = 0.0;
    }
    else
    {
        runoff = runoff + kLake * pow((stage + deltaLevel), nLake);
    }

    /*  Final runoff  */
    newStage = stage - runoff;
    //    if (runoff > 0.0)                                    // New test added in order to allow lake evaporation to draw water below -deltaLevel
    //    {
    //        if (newStage + deltaLevel < 0.0)
    //        {
    //            runoff = runoff + newStage + deltaLevel;
    //            if (runoff < 0.0)
    //            {
    //                runoff = 0.0;
    //            }
    //            newStage = (-1) * deltaLevel;
    //        }
    //    }

    /*  Final lake water level and lake water level change  */
    waterLevelChange = newStage - waterLevel;
    waterLevel = newStage;
    discharge = runoff * (GetLandScapeElement()->GetArea() * GetLandScapeElement()->GetLake()->GetAreaFraction() / 100.0) / commonPar->GetSECONDS_TIMESTEP();

}
