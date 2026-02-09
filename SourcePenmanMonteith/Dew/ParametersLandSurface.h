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
    //    double GetICE_MELT_RATE() const { return ICE_MELT_RATE; }
    //    void SetFREEZE_EFF(double value) { FREEZE_EFF = value; }
    //    double GetFREEZE_EFF() const { return FREEZE_EFF; }
    //    void SetMAX_REL(double value) { MAX_REL = value; }
    //    double GetMAX_REL() const { return MAX_REL; }
    //    void SetALBEDO(double value) { ALBEDO = value; }
    //    double GetALBEDO() const { return ALBEDO; }
    //    void SetCV_SNOW(double value) { CV_SNOW = value; }
    //    double GetCV_SNOW() const { return CV_SNOW; }
    //    void SetSNOW_WEIGHT(int k, double value) { SNOW_WEIGHT[k]= value; }
    //    double GetSNOW_WEIGHT(int k) const { return SNOW_WEIGHT[k]; }
    void SetTREE_HEIGHT(double value) { TREE_HEIGHT = value; }
    double GetTREE_HEIGHT() const { return TREE_HEIGHT; }
    void SetTREE_LAI(double value) { TREE_LAI = value; }
    double GetTREE_LAI() const { return TREE_LAI; }
    void SetTREE_LAI_CORR(double value) { TREE_LAI_CORR = value; }
    double GetTREE_LAI_CORR() const { return TREE_LAI_CORR; }
    void SetBULK_RESISTANCE(double value) { BULK_RESISTANCE = value; }
    double GetBULK_RESISTANCE() const { return BULK_RESISTANCE; }
    void SetALBEDO_SNOW(double value) { ALBEDO_SNOW = value; }
    double GetALBEDO_SNOW() const { return ALBEDO_SNOW; }
    void SetDECIDUOUS_SHARE(double value) { DECIDUOUS_SHARE = value; }
    double GetDECIDUOUS_SHARE() const { return DECIDUOUS_SHARE; }
    void SetWind_H(double value) { WIND_H = value; }
    double GetWind_H() const { return WIND_H; }
    void SetTopen_min(double value) { Topen_min = value; }
    double GetTopen_min() const { return Topen_min; }
    void SetTclose_min(double value) { Tclose_min = value; }
    double GetTclose_min() const { return Tclose_min; }
    void SetVPDclose(double value) { VPDclose = value; }
    double GetVPDclose() const { return VPDclose; }
    void SetVPDopen(double value) { VPDopen = value; }
    double GetVPDopen() const { return VPDopen; }
    void SetGH(double value) { gh = value; }
    double GetGH() const { return gh; }
    void SetCL(double value) { cl = value; }
    double GetCL() const { return cl; }
    void SetZ0G(double value) { z0g = value; }
    double GetZ0G() const { return z0g; }
    void SetGSMAX(double value) { gsmax = value; }
    double GetGSMAX() const { return gsmax; }
    void SetCR(double value) { CR = value; }
    double GetCR() const { return CR; }
    void SetD50(double value) { D50 = value; }
    double GetD50() const { return D50; }
    void SetQ50(double value) { Q50 = value; }
    double GetQ50() const { return Q50; }

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
    double TREE_HEIGHT;
    double TREE_LAI;
    double TREE_LAI_CORR;
    double BULK_RESISTANCE;
    double ALBEDO_SNOW;
    double DECIDUOUS_SHARE;
    double WIND_H;
    double Topen_min;
    double Tclose_min;
    double VPDclose;
    double VPDopen;
    double gh;
    double cl;
    double z0g;
    double gsmax;
    double CR;
    double D50;
    double Q50;
};
