#include "GlacierIce.h"
#include "ParametersGeneral.h"
#include "DistributedElement.h"
#include "InputElement.h"
#include "ParametersLandSurface.h"
#include "stdafx.h"
#include "DateTime.h"

GlacierIce::GlacierIce() :
    precipitation(0.0),
    temp(0.0),
    iceMelt(0.0),
    glacierIceThickness(0.0),
    glacierIceVolume(0.0),
    glacierSurfaceElevation(0.0),
    glacierSurfaceElevationChangeNormalized(0.0),
    restrictedGlacierSurfaceElevationChange(1.0e-6)
{
    SetGeneralPar(0);
    SetLandSurfacePar(0);
    SetGlacierRetreatPar(0);
    //  SetInputTimeSeries(0);
    SetInputElement(0);
    SetLandScapeElement(0);
}

GlacierIce::~GlacierIce()
{
}

void GlacierIce::SetGlacierIceThickness(double value)
{
    glacierIceThickness = value;
}
double GlacierIce::GetGlacierIceThickness() const
{
    return glacierIceThickness;
}
void GlacierIce::SetGlacierIceVolume(double value)
{
    glacierIceVolume = value;
}
double GlacierIce::GetGlacierIceVolume() const
{
    return glacierIceVolume;
}
void GlacierIce::SetGlacierSurfaceElevation(double value)
{
    glacierSurfaceElevation = value;
}
double GlacierIce::GetGlacierSurfaceElevation() const
{
    return glacierSurfaceElevation;
}
void GlacierIce::SetGlacierSurfaceElevationChangeNormalized(double value)
{
    glacierSurfaceElevationChangeNormalized = value;
}
double GlacierIce::GetGlacierSurfaceElevationChangeNormalized() const
{
    return glacierSurfaceElevationChangeNormalized;
}
void GlacierIce::SetRestrictedElevationChange(double value)
{
    restrictedGlacierSurfaceElevationChange = value;
}
double GlacierIce::GetRestrictedElevationChange()
{
    return restrictedGlacierSurfaceElevationChange;
}
void GlacierIce::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral *GlacierIce::GetGeneralPar() const
{
    return commonPar;
}
void GlacierIce::SetLandSurfacePar(ParametersLandSurface *parObj)
{
    landSurfacePar = parObj;
}
ParametersLandSurface *GlacierIce::GetLandSurfacePar()
{
    return landSurfacePar;
}
void GlacierIce::SetGlacierRetreatPar(ParametersGlacierRetreat *parObj)
{
    glacRetPar = parObj;
}
ParametersGlacierRetreat *GlacierIce::GetGlacierRetreatPar() const
{
    return glacRetPar;
}
//  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
//  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
void GlacierIce::SetInputElement(InputElement *inElementObj)
{
    inElement = inElementObj;
}
InputElement *GlacierIce::GetInputElement() const
{
    return inElement;
}
void GlacierIce::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement *GlacierIce::GetLandScapeElement() const
{
    return landScapeElement;
}

double GlacierIce::GetPrecipitation() const
{
    return precipitation;
}

double GlacierIce::GetTemperature() const
{
    return temp;
}

double GlacierIce::GetIceMelt() const
{
    return iceMelt;
}

void GlacierIce::SnowToGlacierIce(double snowStore)
{
    //  cout << " SnowToGlacierIce  " << snowStore << "  " << commonPar->GetDENSITY_ICE() << "  " << snowStore/commonPar->GetDENSITY_ICE() << "\n" ;
    //  cout << " SnowToGlacierIce  glacierIceThickness " << glacierIceThickness << "  glacierSurfaceElevation " << glacierSurfaceElevation << "\n" ;
    glacierIceThickness = glacierIceThickness + snowStore / commonPar->GetDENSITY_ICE();
    glacierSurfaceElevation = glacierSurfaceElevation + snowStore / commonPar->GetDENSITY_ICE();
    //  cout << " SnowToGlacierIce  glacierIceThickness " << glacierIceThickness << "  glacierSurfaceElevation " << glacierSurfaceElevation << "\n" ;
}


void GlacierIce::WaterBalance(int timeStep, DateTime datetime, double glacierIceAreaFraction)
{
    //  double precGradLow, precGradHigh, lapseRate, lapseRateDry, lapseRateWet, elementElevation, gradElevation;
    double elementElevation;
    double weight = 1.0;

    //  cout << "start GlacierIceWaterBalance " << endl;
    precipitation = GetInputElement()->GetInput(0);
    temp = GetInputElement()->GetInput(1);
    elementElevation = GetLandScapeElement()->GetElevation();
    /*  cout << " 1  timeStep " << timeStep << "       Glacier precipitation " << precipitation;
    cout << "    Temperature " << temp << endl;*/

    /*  Precipitation and temperature input correction caused by change in glacier surface elevation  */
    //  cout << "   elementElevation  " << elementElevation << "    glacierSurfaceElevation  " << glacierSurfaceElevation << "\n";
    if (glacierSurfaceElevation != elementElevation)
    {
        //    gradElevation = commonPar->GetGRAD_CHANGE_ALT();
        //    precGradLow = commonPar->GetPREC_GRAD_LOW();
        //    precGradHigh = commonPar->GetPREC_GRAD_HIGH();
        //    lapseRateDry = commonPar->GetLAPSE_DRY();
        //    lapseRateWet = commonPar->GetLAPSE_WET();
        // Precipitation data
        if (precipitation >= 0.0)
        {
            /*      if (gradElevation==0.0 || (glacierSurfaceElevation<=gradElevation && elementElevation<=gradElevation))
            precipitation = (precipitation)*pow(precGradLow,(glacierSurfaceElevation-elementElevation)/100.0);
            else if (glacierSurfaceElevation>gradElevation && elementElevation<gradElevation)
            precipitation = (precipitation)*pow(precGradLow,(gradElevation-elementElevation)/100.0)*
            pow(precGradHigh,(glacierSurfaceElevation-gradElevation)/100.0);
            else if (glacierSurfaceElevation<gradElevation && elementElevation>gradElevation)
            precipitation = (precipitation)*pow(precGradLow,(glacierSurfaceElevation-gradElevation)/100.0)*
            pow(precGradHigh,(gradElevation-elementElevation)/100.0);
            else
            precipitation = (precipitation)*pow(precGradHigh,(glacierSurfaceElevation-elementElevation)/100.0);*/
            precipitation = GetInputElement()->PrecipitationElevationCorrected(GetGeneralPar(), precipitation, elementElevation, glacierSurfaceElevation, weight);
        }
        // Temperature data
        //  if (precipitation == 0.0)
        //    lapseRate = lapseRateDry;
        //  else
        //    lapseRate = lapseRateWet;
        if (temp > -99.0)
        {
            //    temp = (temp + lapseRate*(glacierSurfaceElevation-elementElevation)/100.0);
            temp = GetInputElement()->TemperatureElevationCorrected(GetGeneralPar(), temp, precipitation, elementElevation, glacierSurfaceElevation, weight);
        }
    }
    /*  cout << " 2  timeStep " << timeStep << "       Glacier precipitation " << precipitation;
    cout << "    Temperature " << temp << endl;*/

    /*  Glacier ice melt  */
    //  if (temp > landSurfacePar->GetMELT_TEMP() && GetLandScapeElement()->GetGlacier()->GetSnowStore() == 0.0)
    if (temp > landSurfacePar->GetMELT_TEMP() && glacierIceAreaFraction > 0.0)
    {
        iceMelt = landSurfacePar->GetICE_MELT_RATE() * landSurfacePar->GetSNOW_MELT_RATE() * (temp - landSurfacePar->GetMELT_TEMP());
    }
    else
    {
        iceMelt = 0.0;
    }
    //  cout << temp << "  " << landSurfacePar->GetMELT_TEMP() << "  " << landSurfacePar->GetICE_MELT_RATE() << "  " << GetLandScapeElement()->GetSnowStore() << endl;
    //  cout << "  end GlacierIceWaterBalance " << endl;
}
