#pragma once

class ParametersGeneral;

class TotalReservoirStorage
{
public:
    TotalReservoirStorage();
    ~TotalReservoirStorage();
    void SetGeneralPar(ParametersGeneral *parObj);
    ParametersGeneral *GetGeneralPar() const;
    void AllocateTotalReservoirStorage(int numberTimeSteps);
    void SetInitialTotalReservoirStorage();
    double GetInitialTotalReservoirStorage() const;;
    void SetTotalReservoirStorage(int index, double value);
    double GetTotalReservoirStorage(int index) const;

private:
    ParametersGeneral *commonPar;
    double initialTotalReservoirStorage;
    double *totalReservoirStorage;
};


