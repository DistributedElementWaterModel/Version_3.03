#include "ParametersGeneral.h"



ParametersGeneral::ParametersGeneral() :
    SECONDS_TIMESTEP(0),
    NUM_PREC_SERIES(0),
    NUM_TEMP_SERIES(0),
    DAY_SNOW_ZERO(0),
    DAY_ANNUAL_GLACIER(0),
    PREC_GRAD_LOW(0.0),
    PREC_GRAD_HIGH(0.0),
    GRAD_CHANGE_ALT(0.0),
    PREC_CORR_RAIN(0.0),
    PREC_CORR_SNOW(0.0),
    LAPSE_DRY(0.0),
    LAPSE_WET(0.0),
    DAY_TEMP_MEMORY(0.0),
    LAKE_EPOT_PAR(0.0),
    KLAKE(0.0),
    DELTA_LEVEL(0.0),
    NLAKE(0.0),
    MAXIMUM_LEVEL(0.0),
    DENSITY_ICE(0.0),
    INITIAL_SOIL_MOISTURE(0.0),
    INITIAL_UPPER_ZONE(0.0),
    INITIAL_LOWER_ZONE(0.0),
    INITIAL_SATURATED_ONE(0.0),
    INITIAL_SATURATED_TWO(0.0),
    INITIAL_LAKE_TEMP(0.0),
    INITIAL_LAKE_LEVEL(0.0),
    INITIAL_SNOW(0.0),
    INITIAL_TOTAL_RESERVOIR(0.0)
{
}

ParametersGeneral::~ParametersGeneral()
{
}

void   ParametersGeneral::SetSECONDS_TIMESTEP(int value)
{
    SECONDS_TIMESTEP = value;
}
int   ParametersGeneral::GetSECONDS_TIMESTEP() const
{
    return SECONDS_TIMESTEP;
}
void   ParametersGeneral::SetNUM_PREC_SERIES(int value)
{
    NUM_PREC_SERIES = value;
}
int  ParametersGeneral::GetNUM_PREC_SERIES() const
{
    return NUM_PREC_SERIES;
}
void  ParametersGeneral::SetNUM_TEMP_SERIES(int value)
{
    NUM_TEMP_SERIES = value;
}
int  ParametersGeneral::GetNUM_TEMP_SERIES() const
{
    return NUM_TEMP_SERIES;
}
void  ParametersGeneral::SetPREC_GRAD_LOW(double value)
{
    PREC_GRAD_LOW = value;
}
double  ParametersGeneral::GetPREC_GRAD_LOW() const
{
    return PREC_GRAD_LOW;
}
void  ParametersGeneral::SetPREC_GRAD_HIGH(double value)
{
    PREC_GRAD_HIGH = value;
}
double  ParametersGeneral::GetPREC_GRAD_HIGH() const
{
    return PREC_GRAD_HIGH;
}
void  ParametersGeneral::SetGRAD_CHANGE_ALT(double value)
{
    GRAD_CHANGE_ALT = value;
}
double  ParametersGeneral::GetGRAD_CHANGE_ALT() const
{
    return GRAD_CHANGE_ALT;
}
void  ParametersGeneral::SetPREC_CORR_RAIN(double value)
{
    PREC_CORR_RAIN = value;
}
double  ParametersGeneral::GetPREC_CORR_RAIN() const
{
    return PREC_CORR_RAIN;
}
void  ParametersGeneral::SetPREC_CORR_SNOW(double value)
{
    PREC_CORR_SNOW = value;
}
double  ParametersGeneral::GetPREC_CORR_SNOW() const
{
    return PREC_CORR_SNOW;
}
void  ParametersGeneral::SetLAPSE_DRY(double value)
{
    LAPSE_DRY = value;
}
double  ParametersGeneral::GetLAPSE_DRY() const
{
    return LAPSE_DRY;
}
void  ParametersGeneral::SetLAPSE_WET(double value)
{
    LAPSE_WET = value;
}
double  ParametersGeneral::GetLAPSE_WET() const
{
    return LAPSE_WET;
}
void  ParametersGeneral::SetDAY_TEMP_MEMORY(double value)
{
    DAY_TEMP_MEMORY = value;
}
double  ParametersGeneral::GetDAY_TEMP_MEMORY() const
{
    return DAY_TEMP_MEMORY;
}
void  ParametersGeneral::SetLAKE_EPOT_PAR(double value)
{
    LAKE_EPOT_PAR = value;
}
double  ParametersGeneral::GetLAKE_EPOT_PAR() const
{
    return LAKE_EPOT_PAR;
}
void  ParametersGeneral::SetKLAKE(double value)
{
    KLAKE = value;
}
double  ParametersGeneral::GetKLAKE() const
{
    return KLAKE;
}
void  ParametersGeneral::SetDELTA_LEVEL(double value)
{
    DELTA_LEVEL = value;
}
double  ParametersGeneral::GetDELTA_LEVEL() const
{
    return DELTA_LEVEL;
}
void  ParametersGeneral::SetNLAKE(double value)
{
    NLAKE = value;
}
double  ParametersGeneral::GetNLAKE() const
{
    return NLAKE;
}
void  ParametersGeneral::SetMAXIMUM_LEVEL(double value)
{
    MAXIMUM_LEVEL = value;
}
double  ParametersGeneral::GetMAXIMUM_LEVEL() const
{
    return MAXIMUM_LEVEL;
}
void  ParametersGeneral::SetDENSITY_ICE(double value)
{
    DENSITY_ICE = value;
}
double  ParametersGeneral::GetDENSITY_ICE() const
{
    return DENSITY_ICE;
}
void  ParametersGeneral::SetINITIAL_SOIL_MOISTURE(double value)
{
    INITIAL_SOIL_MOISTURE = value;
}
double  ParametersGeneral::GetINITIAL_SOIL_MOISTURE() const
{
    return INITIAL_SOIL_MOISTURE;
}
void  ParametersGeneral::SetINITIAL_UPPER_ZONE(double value)
{
    INITIAL_UPPER_ZONE = value;
}
double  ParametersGeneral::GetINITIAL_UPPER_ZONE() const
{
    return INITIAL_UPPER_ZONE;
}
void  ParametersGeneral::SetINITIAL_LOWER_ZONE(double value)
{
    INITIAL_LOWER_ZONE = value;
}
double  ParametersGeneral::GetINITIAL_LOWER_ZONE() const
{
    return INITIAL_LOWER_ZONE;
}
void  ParametersGeneral::SetINITIAL_SATURATED_ONE(double value)
{
    INITIAL_SATURATED_ONE = value;
}
double  ParametersGeneral::GetINITIAL_SATURATED_ONE() const
{
    return INITIAL_SATURATED_ONE;
}
void  ParametersGeneral::SetINITIAL_SATURATED_TWO(double value)
{
    INITIAL_SATURATED_TWO = value;
}
double  ParametersGeneral::GetINITIAL_SATURATED_TWO() const
{
    return INITIAL_SATURATED_TWO;
}
void  ParametersGeneral::SetINITIAL_LAKE_TEMP(double value)
{
    INITIAL_LAKE_TEMP = value;
}
double  ParametersGeneral::GetINITIAL_LAKE_TEMP() const
{
    return INITIAL_LAKE_TEMP;
}
void  ParametersGeneral::SetINITIAL_LAKE_LEVEL(double value)
{
    INITIAL_LAKE_LEVEL = value;
}
double  ParametersGeneral::GetINITIAL_LAKE_LEVEL() const
{
    return INITIAL_LAKE_LEVEL;
}
void  ParametersGeneral::SetINITIAL_SNOW(double value)
{
    INITIAL_SNOW = value;
}
double  ParametersGeneral::GetINITIAL_SNOW() const
{
    return INITIAL_SNOW;
}
void  ParametersGeneral::SetINITIAL_TOTAL_RESERVOIR(double value)
{
    INITIAL_TOTAL_RESERVOIR = value;
}
double  ParametersGeneral::GetINITIAL_TOTAL_RESERVOIR()
{
    return INITIAL_TOTAL_RESERVOIR;
}
void  ParametersGeneral::SetDAY_SNOW_ZERO(int value)
{
    DAY_SNOW_ZERO = value;
}
int  ParametersGeneral::GetDAY_SNOW_ZERO() const
{
    return DAY_SNOW_ZERO;
}
void  ParametersGeneral::SetDAY_ANNUAL_GLACIER(int value)
{
    DAY_ANNUAL_GLACIER = value;
}
int  ParametersGeneral::GetDAY_ANNUAL_GLACIER() const
{
    return DAY_ANNUAL_GLACIER;
}
