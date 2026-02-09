#include "ParametersKiWa.h"


ParametersKiWa::ParametersKiWa() :
//  SLOPE_LENGTH(0),
    SOIL_DEPTH(0),
    OV_PAR_1(0),
    OV_PAR_2(0),
    TSAT_0(0),
    EFF_POR(0),
    KSAT_0(0),
    A(0),
    DELTA(0),
    LAMBDA_KW(0),
    ROOT_DEPTH(0),
    WILT_POINT(0),
    EACT_PAR(0.0)
{
}

ParametersKiWa::~ParametersKiWa()
{
}

void  ParametersKiWa::SetSOIL_DEPTH(double value)
{
    SOIL_DEPTH = value;
}
double  ParametersKiWa::GetSOIL_DEPTH() const
{
    return SOIL_DEPTH;
}
void  ParametersKiWa::SetOV_PAR_1(double value)
{
    OV_PAR_1 = value;
}
double  ParametersKiWa::GetOV_PAR_1() const
{
    return OV_PAR_1;
}
void  ParametersKiWa::SetOV_PAR_2(double value)
{
    OV_PAR_2 = value;
}
double  ParametersKiWa::GetOV_PAR_2() const
{
    return OV_PAR_2;
}
void  ParametersKiWa::SetTSAT_0(double value)
{
    TSAT_0 = value;
}
double  ParametersKiWa::GetTSAT_0() const
{
    return TSAT_0;
}
void  ParametersKiWa::SetEFF_POR(double value)
{
    EFF_POR = value;
}
double  ParametersKiWa::GetEFF_POR() const
{
    return EFF_POR;
}
void  ParametersKiWa::SetKSAT_0(double value)
{
    KSAT_0 = value;
}
double  ParametersKiWa::GetKSAT_0() const
{
    return KSAT_0;
}
void  ParametersKiWa::SetA(double value)
{
    A = value;
}
double  ParametersKiWa::GetA() const
{
    return A;
}
void  ParametersKiWa::SetDELTA(double value)
{
    DELTA = value;
}
double  ParametersKiWa::GetDELTA() const
{
    return DELTA;
}
void  ParametersKiWa::SetLAMBDA_KW(double value)
{
    LAMBDA_KW = value;
}
double  ParametersKiWa::GetLAMBDA_KW() const
{
    return LAMBDA_KW;
}
void  ParametersKiWa::SetROOT_DEPTH(double value)
{
    ROOT_DEPTH = value;
}
double  ParametersKiWa::GetROOT_DEPTH() const
{
    return ROOT_DEPTH;
}
void  ParametersKiWa::SetWILT_POINT(double value)
{
    WILT_POINT = value;
}
double  ParametersKiWa::GetWILT_POINT() const
{
    return WILT_POINT;
}
void  ParametersKiWa::SetEACT_PAR(double value)
{
    EACT_PAR = value;
}
double  ParametersKiWa::GetEACT_PAR() const
{
    return EACT_PAR;
}
