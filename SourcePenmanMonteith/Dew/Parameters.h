#pragma once

#include "stdafx.h"

class ParametersLandSurface;
class ParametersSubSurfaceHbv;
class ParametersGeneral;
class ParametersKiWa;
class ParametersGlacierRetreat;

using namespace std;

void SetLandSurfaceParameters(ParametersLandSurface * const ParLandSurfaceStore, ifstream &fileControl, ofstream &fout);
void SetSnowDistribution(ParametersLandSurface * const thisParLandSurface, double cvSnow);
void SetSubSurfaceHbvParameters(ParametersSubSurfaceHbv * const ParSubSurfaceHbvStore, ifstream &fileControl, ofstream &fout);
void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout);
void SetKiWaParameters(ParametersKiWa * const ParKiWaStore, ifstream &fileControl, ofstream &fout);
void SetGlacierRetreatParameters(ParametersGlacierRetreat * const ParGlacierRetStore, ifstream &fileControl, ofstream &fout);
