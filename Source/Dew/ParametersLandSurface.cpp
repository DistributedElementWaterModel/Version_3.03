#include "ParametersLandSurface.h"


ParametersLandSurface::ParametersLandSurface() :
    INTER_MAX(0.0),
    EPOT_PAR(0.0),
    WET_PER_CORR(0.0),
    ACC_TEMP(0.0),
    MELT_TEMP(0.0),
    SNOW_MELT_RATE(0.0),
    ICE_MELT_RATE(0.0),
    FREEZE_EFF(0.0),
    MAX_REL(0.0),
    ALBEDO(0.0),
    CV_SNOW(0.0)
{
    int i;
    for (i = 0; i < numberSnowClasses; i++)
    {
        SNOW_WEIGHT[i] = 0.0;
    }
}

ParametersLandSurface::~ParametersLandSurface()
{
}

void  ParametersLandSurface::SetINTER_MAX(double value)
{
    INTER_MAX = value;
}
double  ParametersLandSurface::GetINTER_MAX() const
{
    return INTER_MAX;
}
void  ParametersLandSurface::SetEPOT_PAR(double value)
{
    EPOT_PAR = value;
}
double  ParametersLandSurface::GetEPOT_PAR() const
{
    return EPOT_PAR;
}
void  ParametersLandSurface::SetWET_PER_CORR(double value)
{
    WET_PER_CORR = value;
}
double  ParametersLandSurface::GetWET_PER_CORR() const
{
    return WET_PER_CORR;
}
void  ParametersLandSurface::SetACC_TEMP(double value)
{
    ACC_TEMP = value;
}
double  ParametersLandSurface::GetACC_TEMP() const
{
    return ACC_TEMP;
}
void  ParametersLandSurface::SetMELT_TEMP(double value)
{
    MELT_TEMP = value;
}
double  ParametersLandSurface::GetMELT_TEMP() const
{
    return MELT_TEMP;
}
void  ParametersLandSurface::SetSNOW_MELT_RATE(double value)
{
    SNOW_MELT_RATE = value;
}
double  ParametersLandSurface::GetSNOW_MELT_RATE() const
{
    return SNOW_MELT_RATE;
}
void  ParametersLandSurface::SetICE_MELT_RATE(double value)
{
    ICE_MELT_RATE = value;
}
double  ParametersLandSurface::GetICE_MELT_RATE() const
{
    return ICE_MELT_RATE;
}
void  ParametersLandSurface::SetFREEZE_EFF(double value)
{
    FREEZE_EFF = value;
}
double  ParametersLandSurface::GetFREEZE_EFF() const
{
    return FREEZE_EFF;
}
void  ParametersLandSurface::SetMAX_REL(double value)
{
    MAX_REL = value;
}
double  ParametersLandSurface::GetMAX_REL() const
{
    return MAX_REL;
}
void  ParametersLandSurface::SetALBEDO(double value)
{
    ALBEDO = value;
}
double  ParametersLandSurface::GetALBEDO() const
{
    return ALBEDO;
}
void  ParametersLandSurface::SetCV_SNOW(double value)
{
    CV_SNOW = value;
}
double  ParametersLandSurface::GetCV_SNOW() const
{
    return CV_SNOW;
}
void  ParametersLandSurface::SetSNOW_WEIGHT(int k, double value)
{
    SNOW_WEIGHT[k] = value;
}
double  ParametersLandSurface::GetSNOW_WEIGHT(int k) const
{
    return SNOW_WEIGHT[k];
}
