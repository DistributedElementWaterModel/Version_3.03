#pragma once

#include "Dew.h"

class ParametersLandSurface
{
public:
    ParametersLandSurface();
    ~ParametersLandSurface();
    void SetINTER_MAX(double value);
    double GetINTER_MAX() const;
    void SetEPOT_PAR(double value);
    double GetEPOT_PAR() const;
    void SetWET_PER_CORR(double value);
    double GetWET_PER_CORR() const;
    void SetACC_TEMP(double value);
    double GetACC_TEMP() const;
    void SetMELT_TEMP(double value);
    double GetMELT_TEMP() const;
    void SetSNOW_MELT_RATE(double value);
    double GetSNOW_MELT_RATE() const;
    void SetICE_MELT_RATE(double value);
    double GetICE_MELT_RATE() const;
    void SetFREEZE_EFF(double value);
    double GetFREEZE_EFF() const;
    void SetMAX_REL(double value);
    double GetMAX_REL() const;
    void SetALBEDO(double value);
    double GetALBEDO() const;
    void SetCV_SNOW(double value);
    double GetCV_SNOW() const;
    void SetSNOW_WEIGHT(int k, double value);
    double GetSNOW_WEIGHT(int k) const;

private:
    double INTER_MAX;
    double EPOT_PAR;
    double WET_PER_CORR;
    double ACC_TEMP;
    double MELT_TEMP;
    double SNOW_MELT_RATE;
    double ICE_MELT_RATE;
    double FREEZE_EFF;
    double MAX_REL;
    double ALBEDO;
    double CV_SNOW;
    double SNOW_WEIGHT[numberSnowClasses];
};
