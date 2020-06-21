#include "ParametersGlacierRetreat.h"


ParametersGlacierRetreat::ParametersGlacierRetreat() :
    A(0.0),
    B(0.0),
    C(0.0),
    GAMMA(0.0)
{
}

ParametersGlacierRetreat::~ParametersGlacierRetreat()
{
}

void  ParametersGlacierRetreat::SetA(double value)
{
    A = value;
}
double  ParametersGlacierRetreat::GetA() const
{
    return A;
}
void  ParametersGlacierRetreat::SetB(double value)
{
    B = value;
}
double  ParametersGlacierRetreat::GetB() const
{
    return B;
}
void  ParametersGlacierRetreat::SetC(double value)
{
    C = value;
}
double  ParametersGlacierRetreat::GetC() const
{
    return C;
}
void  ParametersGlacierRetreat::SetGAMMA(double value)
{
    GAMMA = value;
}
double  ParametersGlacierRetreat::GetGAMMA() const
{
    return GAMMA;
}
void  ParametersGlacierRetreat::SetINCREASE_THRESH(double value)
{
    INCREASE_THRESH = value;
}
double  ParametersGlacierRetreat::GetINCREASE_THRESH() const
{
    return INCREASE_THRESH;
}
void  ParametersGlacierRetreat::SetNUMBER_ADVANCE(int value)
{
    NUMBER_ADVANCE = value;
}
int  ParametersGlacierRetreat::GetNUMBER_ADVANCE() const
{
    return NUMBER_ADVANCE;
}
