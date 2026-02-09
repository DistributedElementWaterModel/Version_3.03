#pragma once
class ParametersGlacierRetreat
{
public:
    ParametersGlacierRetreat();
    ~ParametersGlacierRetreat();
    void SetA(double value);
    double GetA() const;
    void SetB(double value);
    double GetB() const;
    void SetC(double value);
    double GetC() const;
    void SetGAMMA(double value);
    double GetGAMMA() const;
    void SetINCREASE_THRESH(double value);
    double GetINCREASE_THRESH() const;
    void SetNUMBER_ADVANCE(int value);
    int GetNUMBER_ADVANCE() const;

private:
    double A;
    double B;
    double C;
    double GAMMA;
    double INCREASE_THRESH;
    int NUMBER_ADVANCE;
};

