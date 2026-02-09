#pragma once
class ParametersKiWa
{
public:
    ParametersKiWa();
    ~ParametersKiWa();
    //  void SetSLOPE_LENGTH(double value) { SLOPE_LENGTH = value; }
    //  double GetSLOPE_LENGTH() const { return SLOPE_LENGTH; }
    void SetSOIL_DEPTH(double value);
    double GetSOIL_DEPTH() const;
    void SetOV_PAR_1(double value);
    double GetOV_PAR_1() const;
    void SetOV_PAR_2(double value);
    double GetOV_PAR_2() const;
    void SetTSAT_0(double value);
    double GetTSAT_0() const;
    void SetEFF_POR(double value);
    double GetEFF_POR() const;
    void SetKSAT_0(double value);
    double GetKSAT_0() const;
    void SetA(double value);
    double GetA() const;
    void SetDELTA(double value);
    double GetDELTA() const;
    void SetLAMBDA_KW(double value);
    double GetLAMBDA_KW() const;
    void SetROOT_DEPTH(double value);
    double GetROOT_DEPTH() const;
    void SetWILT_POINT(double value);
    double GetWILT_POINT() const;
    void SetEACT_PAR(double value);
    double GetEACT_PAR() const;

private:
    //  double SLOPE_LENGTH;
    double SOIL_DEPTH;
    double OV_PAR_1;
    double OV_PAR_2;
    double TSAT_0;
    double EFF_POR;
    double KSAT_0;
    double A;
    double DELTA;
    double LAMBDA_KW;
    double ROOT_DEPTH;
    double WILT_POINT;
    double EACT_PAR;
};
