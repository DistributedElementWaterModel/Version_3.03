#include "TotalReservoirStorage.h"
#include "ParametersGeneral.h"
#include "Dew.h"

TotalReservoirStorage::TotalReservoirStorage() :
    initialTotalReservoirStorage(0)
{
    SetGeneralPar(0);
}

TotalReservoirStorage::~TotalReservoirStorage()
{
}

void  TotalReservoirStorage::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral * TotalReservoirStorage::GetGeneralPar() const
{
    return commonPar;
}
double  TotalReservoirStorage::GetInitialTotalReservoirStorage() const
{
    return initialTotalReservoirStorage;
};
void  TotalReservoirStorage::SetTotalReservoirStorage(int index, double value)
{
    totalReservoirStorage[index] = value;
}
double  TotalReservoirStorage::GetTotalReservoirStorage(int index) const
{
    return totalReservoirStorage[index];
}

void TotalReservoirStorage::SetInitialTotalReservoirStorage()
{
    initialTotalReservoirStorage = commonPar->GetINITIAL_TOTAL_RESERVOIR();
}

void TotalReservoirStorage::AllocateTotalReservoirStorage(int numberTimeSteps)
{
    int i;
    totalReservoirStorage = new double[numberTimeSteps];
    for (i = 0; i < numberTimeSteps; i++)
    {
        totalReservoirStorage[i] = missingData;
    }
}
