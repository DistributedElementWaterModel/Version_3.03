#include "ParametersSubSurfaceHbv.h"

ParametersSubSurfaceHbv::ParametersSubSurfaceHbv() :
    FC(0.0),
    FCDEL(0.0),
    BETA(0.0),
    INFMAX(0.0),
    KUZ(0.0),
    ALFA(0.0),
    PERC(0.0),
    KLZ(0.0),
    DRAW(0.0)
{
}

ParametersSubSurfaceHbv::~ParametersSubSurfaceHbv()
{
}

void  ParametersSubSurfaceHbv::SetFC(double value)
{
    FC = value;
}
double  ParametersSubSurfaceHbv::GetFC() const
{
    return FC;
}
void  ParametersSubSurfaceHbv::SetFCDEL(double value)
{
    FCDEL = value;
}
double  ParametersSubSurfaceHbv::GetFCDEL() const
{
    return FCDEL;
}
void  ParametersSubSurfaceHbv::SetBETA(double value)
{
    BETA = value;
}
double  ParametersSubSurfaceHbv::GetBETA() const
{
    return BETA;
}
void  ParametersSubSurfaceHbv::SetINFMAX(double value)
{
    INFMAX = value;
}
double  ParametersSubSurfaceHbv::GetINFMAX() const
{
    return INFMAX;
}
void  ParametersSubSurfaceHbv::SetKUZ(double value)
{
    KUZ = value;
}
double  ParametersSubSurfaceHbv::GetKUZ() const
{
    return KUZ;
}
void  ParametersSubSurfaceHbv::SetALFA(double value)
{
    ALFA = value;
}
double  ParametersSubSurfaceHbv::GetALFA() const
{
    return ALFA;
}
void  ParametersSubSurfaceHbv::SetPERC(double value)
{
    PERC = value;
}
double  ParametersSubSurfaceHbv::GetPERC() const
{
    return PERC;
}
void  ParametersSubSurfaceHbv::SetKLZ(double value)
{
    KLZ = value;
}
double  ParametersSubSurfaceHbv::GetKLZ() const
{
    return KLZ;
}
void  ParametersSubSurfaceHbv::SetDRAW(double value)
{
    DRAW = value;
}
double  ParametersSubSurfaceHbv::GetDRAW() const
{
    return DRAW;
}
