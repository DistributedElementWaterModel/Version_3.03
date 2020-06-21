#include "Parameters.h"
#include "ParametersLandSurface.h"
#include "ParametersKiWa.h"
#include "ParametersGeneral.h"
#include "ParametersGlacierRetreat.h"
#include "ParametersSubSurfaceHbv.h"

using namespace std;

void SetLandSurfaceParameters(ParametersLandSurface * const ParLandSurfaceStore, ifstream &fileControl, ofstream &fout)
{
    char fileName[80];
    char buffer[1024];
    int i, j, k;
    double interMax, epotPar, wetPerCorr;
    double accTemp, meltTemp, snowMeltRate, iceMeltRate, freezeEff, maxRel, albedo, cvSnow;

    /*  cout << " File with landsurface parameters: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finLandSurfacePar(fileName);  // Open for reading
    if (!finLandSurfacePar.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        //exit(1);
    }
    finLandSurfacePar.getline(buffer, 1024);
    for (i = 0; i < numberLandSurfaceClasses; i++)
    {
        finLandSurfacePar >> buffer >> j >> interMax >> epotPar >> wetPerCorr >> accTemp >> meltTemp;
        finLandSurfacePar >> snowMeltRate >> iceMeltRate >> freezeEff >> maxRel >> albedo >> cvSnow;
        if (i != j)
        {
            cout << endl << " Error in land surface parameter file, parameter no. "
                 << i << endl << endl;
            //exit(1);
        }
        ParLandSurfaceStore[i].SetINTER_MAX(interMax);
        ParLandSurfaceStore[i].SetEPOT_PAR(epotPar);
        ParLandSurfaceStore[i].SetWET_PER_CORR(wetPerCorr);
        ParLandSurfaceStore[i].SetACC_TEMP(accTemp);
        ParLandSurfaceStore[i].SetMELT_TEMP(meltTemp);
        ParLandSurfaceStore[i].SetSNOW_MELT_RATE(snowMeltRate);
        ParLandSurfaceStore[i].SetICE_MELT_RATE(iceMeltRate);
        ParLandSurfaceStore[i].SetFREEZE_EFF(freezeEff);
        ParLandSurfaceStore[i].SetMAX_REL(maxRel);
        ParLandSurfaceStore[i].SetALBEDO(albedo);
        ParLandSurfaceStore[i].SetCV_SNOW(cvSnow);
        // Set snow distribution parameters
        SetSnowDistribution(&ParLandSurfaceStore[i], cvSnow);
    }
    finLandSurfacePar.close();

    fout << "Land surface parameters: \n";
    for (i = 0; i < numberLandSurfaceClasses; i++)
    {
        fout << ParLandSurfaceStore[i].GetINTER_MAX() << "    ";
        fout << ParLandSurfaceStore[i].GetEPOT_PAR() << "    ";
        fout << ParLandSurfaceStore[i].GetWET_PER_CORR() << "    ";
        fout << ParLandSurfaceStore[i].GetACC_TEMP() << "    ";
        fout << ParLandSurfaceStore[i].GetMELT_TEMP() << "    ";
        fout << ParLandSurfaceStore[i].GetSNOW_MELT_RATE() << "    ";
        fout << ParLandSurfaceStore[i].GetICE_MELT_RATE() << "    ";
        fout << ParLandSurfaceStore[i].GetFREEZE_EFF() << "    ";
        fout << ParLandSurfaceStore[i].GetMAX_REL() << "    ";
        fout << ParLandSurfaceStore[i].GetALBEDO() << "    ";
        fout << ParLandSurfaceStore[i].GetCV_SNOW() << endl;
        for (k = 0; k < numberSnowClasses; k++)
        {
            fout << ParLandSurfaceStore[i].GetSNOW_WEIGHT(k) << "  ";
        }
        fout << endl;
    }
    fout << endl;
}


void SetSnowDistribution(ParametersLandSurface * thisParLandSurface, double cvSnow)
{
    int k;
    double stdDevNorm, meanNorm, sumNorm, sumNorm2;
    double stdNormVar[numberSnowClasses], logNormWeight[numberSnowClasses];
    stdNormVar[0] = -2.326347;
    stdNormVar[1] = -1.644853476;
    stdNormVar[2] = -1.036433474;
    stdNormVar[3] = -0.385320604;
    stdNormVar[4] = 0.385320604;
    stdNormVar[5] = 1.036433474;
    stdNormVar[6] = 1.644853476;
    stdNormVar[7] = 2.326347;
    stdNormVar[8] = 3.71909027;
    //  stdNormVar[8]=4.265043367;
    meanNorm = 0.5 * log(1.0 / (1.0 + cvSnow * cvSnow));
    stdDevNorm = sqrt(log(1.0 + cvSnow * cvSnow));
    sumNorm = 0.0;
    for (k = 0; k < numberSnowClasses; k++)
    {
        logNormWeight[k] = exp(stdNormVar[k] * stdDevNorm + meanNorm);
        sumNorm = sumNorm + logNormWeight[k] * probNorm[k];
    }
    sumNorm2 = 0.0;
    for (k = 0; k < numberSnowClasses; k++)
    {
        logNormWeight[k] = logNormWeight[k] / sumNorm;
        sumNorm2 = sumNorm2 + logNormWeight[k] * probNorm[k];
        thisParLandSurface->SetSNOW_WEIGHT(k, logNormWeight[k]);
        //    cout << k << "  " << logNormWeight[k] << endl;
    }
    if (sumNorm2 < 1.0 - epsilon || sumNorm2 > 1.0 + epsilon)
    {
        cout << endl << " Sum of snow distribution weights = " << sumNorm2 << endl << endl;
        //exit(1);
    }
}


void SetSubSurfaceHbvParameters(ParametersSubSurfaceHbv * const ParSubSurfaceHbvStore, ifstream &fileControl, ofstream &fout)
{
    char fileName[80];
    char buffer[1024];
    int i, j;
    double fc, fcdel, beta, infmax;
    double kuz, alfa, perc, klz, draw;

    /*  cout << " File with HBV subsurface parameters: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finSubSurface(fileName);  // Open for reading
    if (!finSubSurface.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        //exit(1);
    }
    finSubSurface.getline(buffer, 1024);
    for (i = 0; i < numberSoilClasses; i++)
    {
        finSubSurface >> buffer >> j >> fc >> fcdel >> beta >> infmax;
        finSubSurface >> kuz >> alfa >> perc >> klz >> draw;
        if (i != j)
        {
            cout << endl << " Error in HBV subsurface parameter file, parameter no. "
                 << i << endl << endl;
            //exit(1);
        }
        ParSubSurfaceHbvStore[i].SetFC(fc);
        ParSubSurfaceHbvStore[i].SetFCDEL(fcdel);
        ParSubSurfaceHbvStore[i].SetBETA(beta);
        ParSubSurfaceHbvStore[i].SetINFMAX(infmax);
        ParSubSurfaceHbvStore[i].SetKUZ(kuz);
        ParSubSurfaceHbvStore[i].SetALFA(alfa);
        ParSubSurfaceHbvStore[i].SetPERC(perc);
        ParSubSurfaceHbvStore[i].SetKLZ(klz);
        ParSubSurfaceHbvStore[i].SetDRAW(draw);
    }
    finSubSurface.close();
    fout << "HBV subsurface parameters: \n";
    for (i = 0; i < numberSoilClasses; i++)
    {
        fout << ParSubSurfaceHbvStore[i].GetFC() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetFCDEL() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetBETA() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetINFMAX() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetKUZ() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetALFA() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetPERC() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetKLZ() << "    ";
        fout << ParSubSurfaceHbvStore[i].GetDRAW() << "\n";
    }
    fout << endl;
}


void SetKiWaParameters(ParametersKiWa * const ParKiWaStore, ifstream &fileControl, ofstream &fout)
{
    char fileName[80];
    char buffer[1024];
    int i, j;
    //  double slopeLength, soilDepth, ovPar1, ovPar2;
    double soilDepth, ovPar1, ovPar2;
    double tSat0, effPor, kSat0, a, delta, lambdaKw, rootDepth, wiltPoint, eactPar;

    /*  cout << " File with KiWa subsurface parameters: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finKiWa(fileName);  // Open for reading
    if (!finKiWa.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        //exit(1);
    }
    finKiWa.getline(buffer, 1024);
    for (i = 0; i < numberSoilClasses; i++)
    {
        //    finKiWa >> buffer >> j >> slopeLength >> soilDepth >> ovPar1 >> ovPar2 >> tSat0 >> effPor;
        finKiWa >> buffer >> j >> soilDepth >> ovPar1 >> ovPar2 >> tSat0 >> effPor;
        finKiWa >> kSat0 >> a >> delta >> lambdaKw >> rootDepth >> wiltPoint >> eactPar;
        if (i != j)
        {
            cout << endl << " Error in KiWa subsurface parameter file, parameter no. "
                 << i << endl << endl;
            //exit(1);
        }
        //    ParKiWaStore[i].SetSLOPE_LENGTH(slopeLength);
        ParKiWaStore[i].SetSOIL_DEPTH(soilDepth);
        ParKiWaStore[i].SetOV_PAR_1(ovPar1);
        ParKiWaStore[i].SetOV_PAR_2(ovPar2);
        ParKiWaStore[i].SetTSAT_0(tSat0);
        ParKiWaStore[i].SetEFF_POR(effPor);
        ParKiWaStore[i].SetKSAT_0(kSat0);
        ParKiWaStore[i].SetA(a);
        ParKiWaStore[i].SetDELTA(delta);
        ParKiWaStore[i].SetLAMBDA_KW(lambdaKw);
        ParKiWaStore[i].SetROOT_DEPTH(rootDepth);
        ParKiWaStore[i].SetWILT_POINT(wiltPoint);
        ParKiWaStore[i].SetEACT_PAR(eactPar);
    }
    finKiWa.close();
    fout << "KiWa subsurface parameters: \n";
    for (i = 0; i < numberSoilClasses; i++)
    {
        //    fout << ParKiWaStore[i].GetSLOPE_LENGTH() << "    ";
        fout << ParKiWaStore[i].GetSOIL_DEPTH() << "    ";
        fout << ParKiWaStore[i].GetOV_PAR_1() << "    ";
        fout << ParKiWaStore[i].GetOV_PAR_2() << "    ";
        fout << ParKiWaStore[i].GetTSAT_0() << "    ";
        fout << ParKiWaStore[i].GetEFF_POR() << "    ";
        fout << ParKiWaStore[i].GetKSAT_0() << "    ";
        fout << ParKiWaStore[i].GetA() << "    ";
        fout << ParKiWaStore[i].GetDELTA() << "    ";
        fout << ParKiWaStore[i].GetLAMBDA_KW() << "    ";
        fout << ParKiWaStore[i].GetROOT_DEPTH() << "    ";
        fout << ParKiWaStore[i].GetWILT_POINT() << "    ";
        fout << ParKiWaStore[i].GetEACT_PAR() << "\n";
    }
    fout << endl;
}


void SetGlacierRetreatParameters(ParametersGlacierRetreat * const ParGlacRetStore, ifstream &fileControl, ofstream &fout)
{
    char fileName[80];
    char buffer[1024];
    int i, j;
    int numberAdvance;
    double a, b, c, gamma, increaseThresh;

    /*  cout << " File with glacier retreat parameters: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finGlacRet(fileName);  // Open for reading
    if (!finGlacRet.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        //exit(1);
    }
    finGlacRet.getline(buffer, 1024);
    for (i = 0; i < numberGlacierClasses; i++)
    {
        //        finGlacRet >> buffer >> j >> a >> b >> c >> gamma >> increaseThresh >> numberAdvance;
        finGlacRet >> buffer >> j >> a >> b >> c >> gamma;
        if (i != j)
        {
            cout << endl << " Error in glacier retreat parameter file, parameter no. "
                 << i << endl << endl;
            //exit(1);
        }
        ParGlacRetStore[i].SetA(a);
        ParGlacRetStore[i].SetB(b);
        ParGlacRetStore[i].SetC(c);
        ParGlacRetStore[i].SetGAMMA(gamma);
	//        ParGlacRetStore[i].SetINCREASE_THRESH(increaseThresh);
	//        ParGlacRetStore[i].SetNUMBER_ADVANCE(numberAdvance);
    }
    finGlacRet.close();
    fout << "Glacier retreat parameters: \n";
    for (i = 0; i < numberGlacierClasses; i++)
    {
        fout << ParGlacRetStore[i].GetA() << "    ";
        fout << ParGlacRetStore[i].GetB() << "    ";
        fout << ParGlacRetStore[i].GetC() << "    ";
        fout << ParGlacRetStore[i].GetGAMMA() << "    ";
	//        fout << ParGlacRetStore[i].GetINCREASE_THRESH() << "    ";
	//        fout << ParGlacRetStore[i].GetNUMBER_ADVANCE() << "\n";
    }
    fout << endl;
}


void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout)
{
    char fileName[80];
    double precGradLow, precGradHigh, gradChangeAltitude, lapseDry, lapseWet;
    double precCorrRain, precCorrSnow;
    double dayTempMem, lakeEpotPar, kLake, deltaLevel, nLake, maximumLevel, densityIce;
    double initialSoilMoisture, initialUpperZone, initialLowerZone;
    double saturatedFractionOne, saturatedFractionTwo;
    double initialLakeTemp, initialLakeLevel, initialSnow;
    double initialTotalReservoir;
    int secondsTimestep, numPrec, numTemp, daySnowZero, dayAnnualGlacier;

    /*  cout << " File with common parameters: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finGeneralPar(fileName);  // Open for reading
    if (!finGeneralPar.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        //exit(1);
    }
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> secondsTimestep;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> numPrec;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> numTemp;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> precGradLow;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> precGradHigh;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> gradChangeAltitude;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> precCorrRain;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> precCorrSnow;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> lapseDry;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> lapseWet;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> dayTempMem;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> lakeEpotPar;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> kLake;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> deltaLevel;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> nLake;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> maximumLevel;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> densityIce;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> initialSoilMoisture;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> initialUpperZone;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> initialLowerZone;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> saturatedFractionOne;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> saturatedFractionTwo;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> initialLakeTemp;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> initialLakeLevel;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> initialSnow;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> initialTotalReservoir;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> daySnowZero;
    finGeneralPar.ignore(100, ':');
    finGeneralPar >> dayAnnualGlacier;
    ParGeneralStore->SetSECONDS_TIMESTEP(secondsTimestep);
    ParGeneralStore->SetNUM_PREC_SERIES(numPrec);
    ParGeneralStore->SetNUM_TEMP_SERIES(numTemp);
    ParGeneralStore->SetPREC_GRAD_LOW(precGradLow);
    ParGeneralStore->SetPREC_GRAD_HIGH(precGradHigh);
    ParGeneralStore->SetGRAD_CHANGE_ALT(gradChangeAltitude);
    ParGeneralStore->SetPREC_CORR_RAIN(precCorrRain);
    ParGeneralStore->SetPREC_CORR_SNOW(precCorrSnow);
    ParGeneralStore->SetLAPSE_DRY(lapseDry);
    ParGeneralStore->SetLAPSE_WET(lapseWet);
    ParGeneralStore->SetDAY_TEMP_MEMORY(dayTempMem);
    ParGeneralStore->SetLAKE_EPOT_PAR(lakeEpotPar);
    ParGeneralStore->SetKLAKE(kLake);
    ParGeneralStore->SetDELTA_LEVEL(deltaLevel);
    ParGeneralStore->SetNLAKE(nLake);
    ParGeneralStore->SetMAXIMUM_LEVEL(maximumLevel);
    ParGeneralStore->SetDENSITY_ICE(densityIce);
    ParGeneralStore->SetINITIAL_SOIL_MOISTURE(initialSoilMoisture);
    ParGeneralStore->SetINITIAL_UPPER_ZONE(initialUpperZone);
    ParGeneralStore->SetINITIAL_LOWER_ZONE(initialLowerZone);
    ParGeneralStore->SetINITIAL_SATURATED_ONE(saturatedFractionOne);
    ParGeneralStore->SetINITIAL_SATURATED_TWO(saturatedFractionTwo);
    ParGeneralStore->SetINITIAL_LAKE_TEMP(initialLakeTemp);
    ParGeneralStore->SetINITIAL_LAKE_LEVEL(initialLakeLevel);
    ParGeneralStore->SetINITIAL_SNOW(initialSnow);
    ParGeneralStore->SetINITIAL_TOTAL_RESERVOIR(initialTotalReservoir);
    ParGeneralStore->SetDAY_SNOW_ZERO(daySnowZero);
    ParGeneralStore->SetDAY_ANNUAL_GLACIER(dayAnnualGlacier);
    //  ParGeneralStore->SetNumStations(finGeneralPar, numPrec, numTemp);
    finGeneralPar.close();
    if (daySnowZero > 0 && dayAnnualGlacier > 0 && daySnowZero != dayAnnualGlacier)
    {
        cout << "\n" << " File " << fileName << " DAY_SNOW_ZERO " << daySnowZero << " DAY_ANNUAL_GLACIER " << dayAnnualGlacier << "\n\n";
        //exit(1);
    }

    fout << "Common parameters: \n";
    fout << ParGeneralStore->GetSECONDS_TIMESTEP() << endl;
    fout << ParGeneralStore->GetNUM_PREC_SERIES() << endl;
    fout << ParGeneralStore->GetNUM_TEMP_SERIES() << endl;
    fout << ParGeneralStore->GetPREC_GRAD_LOW() << endl;
    fout << ParGeneralStore->GetPREC_GRAD_HIGH() << endl;
    fout << ParGeneralStore->GetGRAD_CHANGE_ALT() << endl;
    fout << ParGeneralStore->GetPREC_CORR_RAIN() << endl;
    fout << ParGeneralStore->GetPREC_CORR_SNOW() << endl;
    fout << ParGeneralStore->GetLAPSE_DRY() << endl;
    fout << ParGeneralStore->GetLAPSE_WET() << endl;
    fout << ParGeneralStore->GetDAY_TEMP_MEMORY() << endl;
    fout << ParGeneralStore->GetLAKE_EPOT_PAR() << endl;
    fout << ParGeneralStore->GetKLAKE() << endl;
    fout << ParGeneralStore->GetDELTA_LEVEL() << endl;
    fout << ParGeneralStore->GetNLAKE() << endl;
    fout << ParGeneralStore->GetMAXIMUM_LEVEL() << endl;
    fout << ParGeneralStore->GetDENSITY_ICE() << endl;
    fout << ParGeneralStore->GetINITIAL_SOIL_MOISTURE() << endl;
    fout << ParGeneralStore->GetINITIAL_UPPER_ZONE() << endl;
    fout << ParGeneralStore->GetINITIAL_LOWER_ZONE() << endl;
    fout << ParGeneralStore->GetINITIAL_SATURATED_ONE() << endl;
    fout << ParGeneralStore->GetINITIAL_SATURATED_TWO() << endl;
    fout << ParGeneralStore->GetINITIAL_LAKE_TEMP() << endl;
    fout << ParGeneralStore->GetINITIAL_LAKE_LEVEL() << endl;
    fout << ParGeneralStore->GetINITIAL_SNOW() << endl;
    fout << ParGeneralStore->GetINITIAL_TOTAL_RESERVOIR() << endl;
    fout << ParGeneralStore->GetDAY_SNOW_ZERO() << endl;
    fout << ParGeneralStore->GetDAY_ANNUAL_GLACIER() << endl;
    /*  fout << " Prec.    ";
    for (i=0; i<ParGeneralStore->GetNUM_PREC_SERIES(); i++) {
    fout << ParGeneralStore->GetSTATION_ALTITUDE(i) << "  ";
    fout << ParGeneralStore->GetSTATION_WEIGHT(i) << "  ";
    }
    fout << endl<< " Temp.    ";
    for (i=0; i<ParGeneralStore->GetNUM_TEMP_SERIES(); i++) {
    fout << ParGeneralStore->GetSTATION_ALTITUDE(ParGeneralStore->GetNUM_PREC_SERIES()+i) << "  ";
    fout << ParGeneralStore->GetSTATION_WEIGHT(ParGeneralStore->GetNUM_PREC_SERIES()+i) << "  ";
    }*/
    fout << endl;
}


