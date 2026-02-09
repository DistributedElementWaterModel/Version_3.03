#include "HBV.h"
#include "ParametersGeneral.h"
#include "ParametersSubSurfaceHbv.h"
#include "DistributedElement.h"
#include "ParametersLandSurface.h"
#include "ModelControl.h"
#include "Util.h"
#include "stdafx.h"
#include "DateTime.h"

HBV::HBV() :
//  temp(0.0),
    soilMoisture(0.0),
    percSoilUpper(0.0),
    upperZone(0.0),
    lowerZone(0.0),
    transpSoilEvap(0.0),
    lowerRunoff(0.0),
    upperRunoff(0.0),
    runoff(0.0)
{
    SetGeneralPar(0);
    SetSubSurfaceHbvPar(0);
    SetLandSurfacePar(0);
    //  SetInputTimeSeries(0);
    SetInputElement(0);
    SetLandScapeElement(0);
}

HBV::~HBV()
{
}

void HBV::SetInitialHbvValues()
{
    soilMoisture = commonPar->GetINITIAL_SOIL_MOISTURE();
    upperZone = commonPar->GetINITIAL_UPPER_ZONE();
    lowerZone = commonPar->GetINITIAL_LOWER_ZONE();
}

void HBV::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
    soilMoisture = sm;
    upperZone = uz;
    lowerZone = lz;
    //  cout << " SetSubSurfaceHbvStore  " << soilMoisture << "  " << upperZone << "  " << lowerZone << endl;
}

void HBV::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral *HBV::GetGeneralPar() const
{
    return commonPar;
}
void HBV::SetLandSurfacePar(ParametersLandSurface *parObj)
{
    landSurfacePar = parObj;
}
ParametersLandSurface *HBV::GetLandSurfacePar()
{
    return landSurfacePar;
}
void HBV::SetSubSurfaceHbvPar(ParametersSubSurfaceHbv *parObj)
{
    subSurfacePar = parObj;
}
ParametersSubSurfaceHbv *HBV::GetSubSurfaceHbvPar() const
{
    return subSurfacePar;
}
//  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
//  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
void HBV::SetInputElement(InputElement *inElementObj)
{
    inElement = inElementObj;
}
InputElement *HBV::GetInputElement() const
{
    return inElement;
}
void HBV::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement *HBV::GetLandScapeElement() const
{
    return landScapeElement;
}

double HBV::GetSoilMoisture() const
{
    return soilMoisture;
}

double HBV::GetSoilMoistureDeficit() const
{
    return subSurfacePar->GetFC() - soilMoisture;
}

double HBV::GetPercSoilUpper() const
{
    return percSoilUpper;
}

double HBV::GetUpperZone() const
{
    return upperZone;
}

double HBV::GetLowerZone() const
{
    return lowerZone;
}

double HBV::GetTranspSoilEvap() const
{
    return transpSoilEvap;
}

double HBV::GetRunoff() const
{
    return runoff;
}

double HBV::GetLowerRunoff() const
{
    return lowerRunoff;
}

double HBV::GetUpperRunoff() const
{
    return upperRunoff;
}


// ** Algorithm to be performed in case: no input to landscape element from upstream elements
//void HBV::WaterBalance(int timeStep, double waterInput, double snowCoverFraction, double dryPeriod)
// ** Algorithm to be performed in case: input to landscape element from upstream elements
void HBV::WaterBalance(int timeStep, DateTime datetime, double waterInput, double temp, double snowCoverFraction, double dryPeriod, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge, double potev)
{
    int i, hoursPerTimeStep;
    double inSoilHour, inSoil, outSoil;
    double lowerPercolation, upperPercolation, upperTemporary;
    double upperZoneMax, lowerZoneMax, drawUp, directRunoff;
    //  temp = GetInputElement()->GetInput(1);
    //  cout << "timeStep " << timeStep << " " << soilMoisture << " " << upperZone << " " << lowerZone << " " << endl;
    //  cout << "\ntimeStep " << timeStep << "   subSurfacePar->GetFC()  " << subSurfacePar->GetFC() << endl;

    // ** Algorithm to be performed in case: input to landscape element from upstream elements through soil surface
    /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
    //  cout << timeStep << "  " << waterInput*1000 << "  " << temp << "\n  " ;
    waterInput = waterInput + (upLandAccumulatedLowerDischarge + upLandAccumulatedUpperDischarge) * commonPar->GetSECONDS_TIMESTEP() / GetLandScapeElement()->GetArea();
    // ** End algorithm to be performed in case: input to landscape element from upstream elements through soil surface

    // ** Algorithm to be performed in case: input to landscape element from upstream elements through subsurface
    /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
    //  upperZone = upperZone + upLandAccumulatedUpperDischarge*commonPar->GetSECONDS_TIMESTEP()/GetLandScapeElement()->GetArea();
    //  lowerZone = lowerZone + upLandAccumulatedLowerDischarge*commonPar->GetSECONDS_TIMESTEP()/GetLandScapeElement()->GetArea();
    // ** End algorithm to be performed in case: input to landscape element from upstream elements through subsurface

    /* Maximum upper zone storage */
    upperZoneMax = commonPar->GetMAXIMUM_LEVEL();
    //  upperZoneMax = pow(subSurfacePar->GetINFMAX()/subSurfacePar->GetKUZ(),(1/subSurfacePar->GetALFA()));

    /* Maximum lower zone storage */
    lowerZoneMax = subSurfacePar->GetPERC() / subSurfacePar->GetKLZ();

    /* Capillary rise from lower zone to soil moisture zone */
    drawUp = (2.0 * subSurfacePar->GetDRAW()) * (lowerZone / lowerZoneMax) * (subSurfacePar->GetFC() - soilMoisture) / (subSurfacePar->GetFC());
    lowerZone = lowerZone - drawUp;
    if (lowerZone < 0.0)
    {
        drawUp = drawUp + lowerZone;
        lowerZone = 0.0;
    }
    soilMoisture = soilMoisture + drawUp;

    /* Water exceeding lower zone maximum to upper zone */
    if (lowerZone > lowerZoneMax)
    {
        upperZone = upperZone + lowerZone - lowerZoneMax;
        lowerZone = lowerZoneMax;
    }

    /* Water exceeding upper zone maximum to runoff */
    if (upperZone > upperZoneMax)
    {
        directRunoff = upperZone - upperZoneMax;
        upperZone = upperZoneMax;
    }
    else
    {
        directRunoff = 0.0;
    }

    /* Soil moisture deficit in use for glacier-free areas */
    if (dryPeriod > missingData)
    {
        /* Infiltration to soil moisture zone */
        if (waterInput > subSurfacePar->GetINFMAX())
        {
            upperPercolation = waterInput - subSurfacePar->GetINFMAX();
            inSoil = subSurfacePar->GetINFMAX();
        }
        else
        {
            upperPercolation = 0.0;
            inSoil = waterInput;
        }

        //    cout << " subSurfaceHbv   snowCoverFraction   " << snowCoverFraction << endl;
        /* Transpiration and soil evaporation for snow free areas */
        if (dryPeriod > 0.0)
        {
            //    if (dryPeriod > 0.0 && GetLandScapeElement()->GetSnowStore() <= 0.0) {
            //    if (dryPeriod > 0.0 && GetLandScapeElement()->GetSnowStore() > 0.0) {
           if (GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'P' || 
	       GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'p') 
	   {
	     transpSoilEvap = dryPeriod * (1.0 - snowCoverFraction) * 
	       HBVTranspSoilEvap (soilMoisture,  potev, 
				  subSurfacePar->GetFC(), subSurfacePar->GetFCDEL());
	   }
           else if (GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'M' || 
	       GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'm') 
	   {
	       i = datetime.getMonth()-1;
	       transpSoilEvap = dryPeriod * (1.0 - snowCoverFraction) *
		 HBVTranspSoilEvapLongTermMean(soilMoisture, temp, landSurfacePar->GetEPOT_PAR(),
                                               subSurfacePar->GetFC(), subSurfacePar->GetFCDEL(),
					       GetLandScapeElement()->GetModelControlObj()->GetEvaporationArray(i)/1000.0);
	   }
	   else
	   {
	       transpSoilEvap = dryPeriod * (1.0 - snowCoverFraction) * 
		 HBVTranspSoilEvapTemperatureIndex(soilMoisture, temp, landSurfacePar->GetEPOT_PAR(), 
						   subSurfacePar->GetFC(), subSurfacePar->GetFCDEL());
	   }
	}
        else
        {
            transpSoilEvap = 0.0;
        }
        soilMoisture = soilMoisture - transpSoilEvap;
        if (soilMoisture < 0.0)
        {
            transpSoilEvap = transpSoilEvap + soilMoisture;
            soilMoisture = 0.0;
        }

        /* Water balance for soil moisture zone per hour, percolation to upper groundwater zone */
        if (inSoil > 0.0)
        {
            hoursPerTimeStep = commonPar->GetSECONDS_TIMESTEP() / minimumTimeStep;
            inSoilHour = inSoil / hoursPerTimeStep;
            //      cout << " hoursPerTimeStep : " << hoursPerTimeStep << " * " << minimumTimeStep << " = " << commonPar->GetSECONDS_TIMESTEP() << endl;
            //      cout << " inSoilHour : " << inSoilHour << " * " << hoursPerTimeStep << " = " << inSoil << endl;
            for (i = 0; i < hoursPerTimeStep; i++)
            {
                if (soilMoisture < subSurfacePar->GetFC())
                {
                    outSoil = inSoilHour * pow(soilMoisture / subSurfacePar->GetFC(), subSurfacePar->GetBETA());
                }
                else
                {
                    outSoil = inSoilHour;
                }
                soilMoisture = soilMoisture + inSoilHour - outSoil;
                if (soilMoisture > subSurfacePar->GetFC())
                {
                    outSoil = outSoil + soilMoisture - subSurfacePar->GetFC();
                    soilMoisture = subSurfacePar->GetFC();
                }
                upperPercolation = upperPercolation + outSoil;
            }
        }
        else
        {
            outSoil = 0.0;
        }
        //    cout << "upperPercolation not Glacier  " << upperPercolation << endl;
    }

    /* Soil moisture deficit not in use */
    else
    {
        /*    transpSoilEvap = 0.0;
        upperPercolation = waterInput;
        soilMoisture = subSurfacePar->GetFC();*/
        transpSoilEvap = 0.0;
        upperPercolation = 0.0;
        inSoil = waterInput;

        /* Water balance for soil moisture zone, percolation to upper groundwater zone */
        if (inSoil > 0.0)
        {
            if (soilMoisture < subSurfacePar->GetFC())
            {
                outSoil = inSoil * power(soilMoisture / subSurfacePar->GetFC(), subSurfacePar->GetBETA());
            }
            else
            {
                outSoil = inSoil;
            }
            soilMoisture = soilMoisture + inSoil - outSoil;
            if (soilMoisture > subSurfacePar->GetFC())
            {
                outSoil = outSoil + soilMoisture - subSurfacePar->GetFC();
                soilMoisture = subSurfacePar->GetFC();
            }
            upperPercolation = upperPercolation + outSoil;
        }
        else
        {
            outSoil = 0.0;
        }
        //    cout << "Glacier:   upperPercolation " << upperPercolation << "  ";
    }

    /* Water balance for upper groundwater zone, percolation to lower groundwater zone */
    lowerPercolation = subSurfacePar->GetPERC();
    //  if (dryPeriod == missingData) cout << " subSurfacePar->GetPERC() " << subSurfacePar->GetPERC() << "  ";
    //  if (dryPeriod == missingData) cout << " lowerPercolation " << lowerPercolation << "  ";
    upperTemporary = upperZone + 0.5 * (upperPercolation - lowerPercolation);
    //  if (dryPeriod == missingData) cout << " upperTemporary " << upperTemporary << "  ";
    if (lowerPercolation > upperZone + upperPercolation)
    {
        lowerPercolation = upperZone + upperPercolation;
        upperRunoff = 0.0;
    }
    else
    {
        upperRunoff = subSurfacePar->GetKUZ() * pow(upperTemporary, subSurfacePar->GetALFA());
    }
    upperZone = upperZone + upperPercolation - lowerPercolation - upperRunoff;
    if (upperZone < 0.0)
    {
        upperRunoff = upperRunoff + upperZone;
        if (upperRunoff < 0.0)
        {
            lowerPercolation = lowerPercolation + upperRunoff;
            upperRunoff = 0.0;
        }
        if (lowerPercolation < 0.0)
        {
            lowerPercolation = 0.0;
        }
        upperZone = 0.0;
    }
    //  if (dryPeriod == missingData) cout << " lowerPercolation " << lowerPercolation << "  " << endl;

    /* Water balance for lower groundwater zone */
    lowerZone = lowerZone + lowerPercolation;
    if (lowerZone < 0.0)
    {
        lowerZone = 0.0;
    }
    lowerRunoff = lowerZone * subSurfacePar->GetKLZ();
    lowerZone = lowerZone - lowerRunoff;

    /* Add direct runoff to upper runoff */
    upperRunoff = upperRunoff + directRunoff;

    /* Runoff from both zones */
    runoff = lowerRunoff + upperRunoff;

    /* Percolation from soil moisture zone to upper zone */
    percSoilUpper = upperPercolation;
}
