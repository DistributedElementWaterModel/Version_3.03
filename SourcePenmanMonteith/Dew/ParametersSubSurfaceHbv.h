#pragma once
class ParametersSubSurfaceHbv
{
public:
    ParametersSubSurfaceHbv();
    ~ParametersSubSurfaceHbv();
    void SetFC(double value);
    double GetFC() const;
    void SetFCDEL(double value);
    double GetFCDEL() const;
    void SetBETA(double value);
    double GetBETA() const;
    void SetINFMAX(double value);
    double GetINFMAX() const;
    void SetKUZ(double value);
    double GetKUZ() const;
    void SetALFA(double value);
    double GetALFA() const;
    void SetPERC(double value);
    double GetPERC() const;
    void SetKLZ(double value);
    double GetKLZ() const;
    void SetDRAW(double value);
    double GetDRAW() const;

private:
    double FC;
    double FCDEL;
    double BETA;
    double INFMAX;
    double KUZ;
    double ALFA;
    double PERC;
    double KLZ;
    double DRAW;
};
