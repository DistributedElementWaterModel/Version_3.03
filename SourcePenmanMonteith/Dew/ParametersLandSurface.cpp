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
    CV_SNOW(0.0),
    TREE_HEIGHT(0.0),
    TREE_LAI(0.0),
    TREE_LAI_CORR(0.0),
    BULK_RESISTANCE(0.0),
    ALBEDO_SNOW(0.0),
    DECIDUOUS_SHARE(0.0),
    WIND_H(0.0),
    Topen_min(0.0),
    Tclose_min(0.0),
    VPDclose(0.0),
    VPDopen(0.0),
    gh(0.0),
    cl(0.0),
    z0g(0.0),
    gsmax(0.0),
    CR(0.0),
    D50(0.0),
    Q50(0.0)

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
void  ParametersLandSurface::SetTREE_HEIGHT(double value)
{
    TREE_HEIGHT = value;
}
double  ParametersLandSurface::GetTREE_HEIGHT() const
{
    return TREE_HEIGHT;
}
void  ParametersLandSurface::SetTREE_LAI(double value)
{
    TREE_LAI = value;
}
double  ParametersLandSurface::GetTREE_LAI() const
{
    return TREE_LAI;
}
void  ParametersLandSurface::SetTREE_LAI_CORR(double value)
{
    TREE_LAI_CORR = value;
}
double  ParametersLandSurface::GetTREE_LAI_CORR() const
{
    return TREE_LAI_CORR;
}
void  ParametersLandSurface::SetBULK_RESISTANCE(double value)
{
    BULK_RESISTANCE = value;
}
double  ParametersLandSurface::GetBULK_RESISTANCE() const
{
    return BULK_RESISTANCE;
}
void  ParametersLandSurface::SetALBEDO_SNOW(double value)
{
    ALBEDO_SNOW = value;
}
double  ParametersLandSurface::GetALBEDO_SNOW() const
{
    return ALBEDO_SNOW;
}
void  ParametersLandSurface::SetDECIDUOUS_SHARE(double value)
{
    DECIDUOUS_SHARE = value;
}
double  ParametersLandSurface::GetDECIDUOUS_SHARE() const
{
    return DECIDUOUS_SHARE;
}
void  ParametersLandSurface::SetWind_H(double value)
{
    WIND_H = value;
}
double  ParametersLandSurface::GetWind_H() const
{
    return WIND_H;
}
void  ParametersLandSurface::SetTopen_min(double value)
{
    Topen_min = value;
}
double  ParametersLandSurface::GetTopen_min() const
{
    return Topen_min;
}
void  ParametersLandSurface::SetTclose_min(double value)
{
    Tclose_min = value;
}
double  ParametersLandSurface::GetTclose_min() const
{
    return Tclose_min;
}
void  ParametersLandSurface::SetVPDclose(double value)
{
    VPDclose = value;
}
double  ParametersLandSurface::GetVPDclose() const
{
    return VPDclose;
}
void  ParametersLandSurface::SetVPDopen(double value)
{
    VPDopen = value;
}
double  ParametersLandSurface::GetVPDopen() const
{
    return VPDopen;
}
void  ParametersLandSurface::SetGH(double value)
{
    gh = value;
}
double  ParametersLandSurface::GetGH() const
{
    return gh;
}
void  ParametersLandSurface::SetCL(double value)
{
    cl = value;
}
double  ParametersLandSurface::GetCL() const
{
    return cl;
}
void  ParametersLandSurface::SetZ0G(double value)
{
    z0g = value;
}
double  ParametersLandSurface::GetZ0G() const
{
    return z0g;
}
void  ParametersLandSurface::SetGSMAX(double value)
{
    gsmax = value;
}
double  ParametersLandSurface::GetGSMAX() const
{
    return gsmax;
}
void  ParametersLandSurface::SetCR(double value)
{
    CR = value;
}
double  ParametersLandSurface::GetCR() const
{
    return CR;
}
void  ParametersLandSurface::SetD50(double value)
{
    D50 = value;
}
double  ParametersLandSurface::GetD50() const
{
    return D50;
}
void  ParametersLandSurface::SetQ50(double value)
{
    Q50 = value;
}
double  ParametersLandSurface::GetQ50() const
{
    return Q50;
}
