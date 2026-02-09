#include "KinematicWave.h"
#include "stdafx.h"
#include "DistributedElement.h"
#include "ParametersGeneral.h"
#include "ParametersKiWa.h"
#include "ParametersLandSurface.h"
#include "ModelControl.h"
#include "Util.h"
#include "stdafx.h"
#include "DateTime.h"

KinematicWave::KinematicWave() :
//  temp(0.0),
    vap_def(0.0),
    u_high(0),
    mean_perc(0.0),
    field_capacity(0.0),
    teta(0.0),
    gw_h(0.0),
    evaporated_volume(0.0),
    evol_upper(0.0),
    volume_root(0.0),
    def_par(0.0),
    actev(0.0),
    actev_loss(0.0),
    transpSoilEvap(0.0),
    lower_runoff(0.0),
    upper_runoff(0.0),
    runoff(0.0)
{
    SetKiWaPar(0);
    SetGeneralPar(0);
    SetLandSurfacePar(0);
    //  SetInputTimeSeries(0);
    SetInputElement(0);
    SetLandScapeElement(0);
}

KinematicWave::~KinematicWave()
{
}

void KinematicWave::SetInitialKwaValues()
{
}

void  KinematicWave::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral * KinematicWave::GetGeneralPar() const
{
    return commonPar;
}
void  KinematicWave::SetLandSurfacePar(ParametersLandSurface *parObj)
{
    landSurfacePar = parObj;
}
ParametersLandSurface * KinematicWave::GetLandSurfacePar()
{
    return landSurfacePar;
}
void  KinematicWave::SetKiWaPar(ParametersKiWa *parObj)
{
    kiWaPar = parObj;
}
ParametersKiWa * KinematicWave::GetKiWaPar() const
{
    return kiWaPar;
}
//  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
//  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
void  KinematicWave::SetInputElement(InputElement *inElementObj)
{
    inElement = inElementObj;
}
InputElement * KinematicWave::GetInputElement() const
{
    return inElement;
}
void  KinematicWave::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * KinematicWave::GetLandScapeElement() const
{
    return landScapeElement;
}

double KinematicWave::GetTranspSoilEvap() const
{
    return transpSoilEvap;
}

double KinematicWave::GetRunoff() const
{
    return runoff;
}

double KinematicWave::GetLowerRunoff() const
{
    return lower_runoff;
}

double KinematicWave::GetUpperRunoff() const
{
    return upper_runoff;
}


void KinematicWave::WaterBalance(int timeStep, DateTime datetime, double waterInput, double temp, double snowCoverFraction, double dry_period, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge)
{
  int i, j;
    double SLOPE_LENGTH, SLOPE_ANGLE;
    //  temp = GetInputElement()->GetInput(1);
    SLOPE_LENGTH = GetLandScapeElement()->GetSlopeLength();
    SLOPE_ANGLE = GetLandScapeElement()->GetSlopeAngle() * acos(-1.0) / 180.0;

    // ** Algorithm to be performed in case: input to landscape element from upstream elements
    /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
    waterInput = waterInput + (upLandAccumulatedLowerDischarge + upLandAccumulatedUpperDischarge) * commonPar->GetSECONDS_TIMESTEP() / GetLandScapeElement()->GetArea();
    // ** End algorithm to be performed in case: input to landscape element from upstream elements

    if (timeStep == 0)
    {
        for (i = 1; i <= numberCharacteristic; i++)
        {
            smdef[i] = 0;
            perc[i] = 0;
            len_coord[i] = ((double)(i - 1) / (double)(numberCharacteristic - 1)) * SLOPE_LENGTH;
            fixed_length[i] = len_coord[i];
            sat_depth[i] = kiWaPar->GetSOIL_DEPTH() * (commonPar->GetINITIAL_SATURATED_ONE() +
                           (commonPar->GetINITIAL_SATURATED_TWO() - commonPar->GetINITIAL_SATURATED_ONE()) *
                           (i - 1) / (numberCharacteristic - 1));
            fixed_sat[i] = 0;
            upp_tim[i] = 0;
            upp_dep[i] = 0;
            upp_len[i] = 0;
        }
    }
    actev = 0;
    actev_loss = 0;
    mean_perc = 0;
    lower_runoff = 0;
    upper_runoff = 0;
    //  cout << "timeStep " << timeStep << "    Kwa waterInput " << waterInput << endl;
    //  cout << "timeStep " << timeStep << "    i " << "50" << "    smdef " << smdef[50] << endl;
    if (waterInput > 0)
    {
        for (i = 1; i <= numberCharacteristic; i++)
        {
            smdef[i] = smdef[i] + waterInput;
            if (smdef[i] >= 0)
            {
                perc[i] = smdef[i];
                smdef[i] = 0;
            }
            else
            {
                perc[i] = 0;
            }
            mean_perc = mean_perc + perc[i];
        }
        mean_perc = mean_perc / numberCharacteristic;
    }
    //  cout << "timeStep " << timeStep << "    mean_perc " << mean_perc << endl;
    //  cout << "timeStep " << timeStep << "    sat_depth[50] " << sat_depth[50] << endl;

    if (mean_perc > 0)     /* Kinematic wave with lateral inflow */
        kinematic_wave_with_lateral_inflow(mean_perc, len_coord, sat_depth, &lower_runoff, &upper_runoff, upp_tim, upp_len, upp_dep, &u_high,
                                           SLOPE_LENGTH, SLOPE_ANGLE, kiWaPar->GetSOIL_DEPTH(),
                                           kiWaPar->GetOV_PAR_1(), kiWaPar->GetOV_PAR_2(), kiWaPar->GetEFF_POR(),
                                           kiWaPar->GetKSAT_0(), kiWaPar->GetA());
    else                   /* Kinematic wave without lateral inflow */
        kinematic_wave_without_lateral_inflow(len_coord, sat_depth, &lower_runoff, &upper_runoff, upp_tim, upp_len, upp_dep, &u_high,
                                              SLOPE_LENGTH, SLOPE_ANGLE, kiWaPar->GetSOIL_DEPTH(),
                                              kiWaPar->GetOV_PAR_1(), kiWaPar->GetOV_PAR_2(), kiWaPar->GetEFF_POR(),
                                              kiWaPar->GetKSAT_0(), kiWaPar->GetA());
    //  cout << "timeStep " << timeStep << "    lower_runoff " << lower_runoff << "    upper_runoff " << upper_runoff << endl;
    /*  Redistribution of saturated depth profile to fixed length coordinates  */
    i = numberCharacteristic;
    j = numberCharacteristic;
    while (len_coord[numberCharacteristic] <= fixed_length[j])
    {
        fixed_sat[j] = sat_depth[numberCharacteristic];
        j--;
    }
    while (j > 0)
    {
        while (len_coord[i] > fixed_length[j] && i > 1)
        {
            i--;
        }
        if (len_coord[i] > fixed_length[j])
        {
            fixed_sat[j] = sat_depth[i];
        }
        else
        {
            fixed_sat[j] = sat_depth[i] + (fixed_length[j] - len_coord[i]) * (sat_depth[i + 1] - sat_depth[i]) / (len_coord[i + 1] - len_coord[i]);
        }
        j--;
    }
    for (i = 1; i <= numberCharacteristic; i++)
    {
        len_coord[i] = fixed_length[i];
        if (fixed_sat[i] > 0)
        {
            sat_depth[i] = fixed_sat[i];
        }
        else
        {
            sat_depth[i] = 0;
        }
    }
    /* If dry period > 0, calculate volumetric water content at the soil surface, actual evapotranspiration,
    saturated depth and soil moisture deficit along each characteristic curve in lower layer.
    If upper layer detph > 0, the water demanded by evapotranspiration is extracted from this layer first */
    if (dry_period > 0)
    {
        /* Upper layer */
        /*   if (u_high > 0) {
        evol_upper = dry_period * (1.0 - snowCoverFraction) * potentialEvap (temp, landSurfacePar->GetEPOT_PAR());
        for (j = 1; j <= u_high; j++) {
        upp_dep[j] = upp_dep[j] - evol_upper;
        }
        }*/
        /* Lower layer */
        for (i = 1; i <= numberCharacteristic; i++)
        {
            gw_h = (kiWaPar->GetSOIL_DEPTH() - sat_depth[i]) / cos(SLOPE_ANGLE);
            field_capacity = kiWaPar->GetTSAT_0() * exp(kiWaPar->GetLAMBDA_KW() * gw_h);
            teta = field_capacity + smdef[i] / kiWaPar->GetROOT_DEPTH();
            /*        if (gw_h < 0-EPS_KW) {
            printf("    gw_h = %f < 0    i = %d    index = %d \n",gw_h,i,index);
            exit(1);
            }*/
            /*        if (teta < kiWaPar->GetWILT_POINT()-kiWaPar->GetEPS_KW() || teta > kiWaPar->GetTSAT_0()+kiWaPar->GetEPS_KW()) {
            if (teta < 0-kiWaPar->GetEPS_KW() || teta > kiWaPar->GetTSAT_0()+kiWaPar->GetEPS_KW()) {
            printf("\n\n    kiWaPar->GetWILT_POINT() = %f    teta = %f    gw_h = %f    fc = %f    smdef = %f    i = %d    index = %d \n\n",
            kiWaPar->GetWILT_POINT(),teta,gw_h,field_capacity,smdef[i],i,index);
            exit(1);
            }*/
            if (teta < 0)
            {
                teta = 0;
            }
            if (teta > kiWaPar->GetTSAT_0())
            {
                teta = kiWaPar->GetTSAT_0();
            }
           if (GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'M' || 
	       GetLandScapeElement()->GetModelControlObj()->GetEvaporationModelling() == 'm') 
	   {
	       j = datetime.getMonth()-1;
	       evaporated_volume = dry_period * (1.0 - snowCoverFraction) * KiWaTranspSoilEvapLongTermMean(teta, temp, landSurfacePar->GetEPOT_PAR(),
				   kiWaPar->GetEACT_PAR(), kiWaPar->GetTSAT_0(), kiWaPar->GetWILT_POINT(),
				   GetLandScapeElement()->GetModelControlObj()->GetEvaporationArray(i)/1000.0);
	   }
	   else
	   {
	       evaporated_volume = dry_period * (1.0 - snowCoverFraction) * KiWaTranspSoilEvapTemperatureIndex(teta, temp, landSurfacePar->GetEPOT_PAR(),
				   kiWaPar->GetEACT_PAR(), kiWaPar->GetTSAT_0(), kiWaPar->GetWILT_POINT());
	   }
           /* Find index of characteristic curve in upper layer with length coordinate <= length coordinate along characteristic curve in lower layer */
            /*     if (u_high > 0) {
            if (len_coord[i] >= upp_len[u_high] ) {
            j = u_high;
            if (j > 1) {
            while (len_coord[i] >= upp_len[j-1] && j > 1) j--;
            }
            if (upp_dep[j] < 0)
            evaporated_volume = (-1)*upp_dep[j];
            else
            evaporated_volume = 0;
            }
            }*/
            /* Saturated depth and soil moisture deficit in lower layer */
            if (evaporated_volume > 0)
            {
                if (gw_h > kiWaPar->GetROOT_DEPTH())
                {
                    def_par = 1.0 - exp(kiWaPar->GetDELTA() * (gw_h - kiWaPar->GetROOT_DEPTH()));
                    gw_h = gw_h + ((evaporated_volume * (1.0 - def_par)) / kiWaPar->GetEFF_POR());
                    smdef[i] = smdef[i] - evaporated_volume * def_par;
                }
                else
                {
                    volume_root = kiWaPar->GetEFF_POR() * (kiWaPar->GetROOT_DEPTH() - gw_h);
                    if (volume_root >= evaporated_volume)
                    {
                        gw_h = gw_h + (evaporated_volume / kiWaPar->GetEFF_POR());
                    }
                    else
                    {
                        gw_h = kiWaPar->GetROOT_DEPTH() + ((evaporated_volume - volume_root) / kiWaPar->GetEFF_POR());
                    }
                }
                if (gw_h > kiWaPar->GetSOIL_DEPTH() / cos(SLOPE_ANGLE))
                {
                    evaporated_volume = evaporated_volume - (gw_h - kiWaPar->GetSOIL_DEPTH() / cos(SLOPE_ANGLE)) * kiWaPar->GetEFF_POR();
                    gw_h = kiWaPar->GetSOIL_DEPTH() / cos(SLOPE_ANGLE);
                }
                sat_depth[i] = kiWaPar->GetSOIL_DEPTH() - (gw_h * cos(SLOPE_ANGLE));
            }
            /*      actev = actev + (1.0 - snowCoverFraction) * KiWaTranspSoilEvapTemperatureIndex(teta, temp, landSurfacePar->GetEPOT_PAR(),
            kiWaPar->GetEACT_PAR(), kiWaPar->GetTSAT_0(), kiWaPar->GetWILT_POINT());*/
            actev = actev + evaporated_volume;
        }
        actev = actev / (numberCharacteristic);
        actev_loss = actev * dry_period;
    }
    /* Soil moisture deficit must be less or equal to zero */
    for (i = 1; i <= numberCharacteristic; i++)
    {
        if (smdef[i] > 0 + epsilon)
        {
            printf("    smdef[i] = %f > 0    i = %d    index = %d \n", smdef[i], i, timeStep);
            /*        exit(1);*/
        }
    }
    /* New value of u_high in upper layer */
    if (u_high > 0)
    {
        j = u_high;
        while (upp_dep[j] <= 0 && j > 0)
        {
            j--;
        }
        u_high = j;
    }
    /* Control saturated depth profile, soil moisture deficit */
    for (i = 1; i <= numberCharacteristic; i++)
    {
        if (sat_depth[i] < 0)
        {
            sat_depth[i] = 0;
        }
        if (smdef[i] > 0)
        {
            smdef[i] = 0;
        }
    }
    /* Combine runoff from upper and lower layer */
    runoff = lower_runoff + upper_runoff;

    /* Evapotranspiration */
    transpSoilEvap = actev_loss;

    /* Output data */
    //  if (GetLandScapeElement()->GetLandIndex() == selectedLandIndex) cout << "evapotranspiration     "
    //              << (GetLandScapeElement()->GetSoilEvaporation() + GetLandScapeElement()->GetInterceptionLoss())*1000 << endl;
    /*  if (GetLandScapeElement()->GetLandIndex() == selectedLandIndex) {
    cout << "groundwater slope    ";
    cout.width(10);cout.precision(5);
    cout << -(kiWaPar->GetSOIL_DEPTH() - sat_depth[numberCharacteristic/2]) << endl;
    }*/
    //  cout << GetLandScapeElement()->GetLandIndex() << "   " << kiWaPar->GetNumberSelectedKiWaHillslopeElements() << endl;
    /*  if (GetLandScapeElement()->GetLandIndex() == selectedLandIndex) {
    //    if (timeStep == 943) {
    if (timeStep == 2159) {
    for (i = 1; i <= numberCharacteristic; i++) {
    gw_h = kiWaPar->GetSOIL_DEPTH() - sat_depth[i];
    field_capacity = kiWaPar->GetTSAT_0() * exp(kiWaPar->GetLAMBDA_KW() * gw_h);
    teta = field_capacity + smdef[i]/kiWaPar->GetROOT_DEPTH();
    cout << len_coord[i] << "    " << -gw_h << "    " << teta << endl;
    }
    }
    }*/
}


void KinematicWave::kinematic_wave_with_lateral_inflow(double inflow, double *len_coord, double *sat_depth, double *lower_runoff, double *upper_runoff,
        double *upp_tim, double *upp_len, double *upp_dep, int *u_high, double SLOPE_LENGTH, double SLOPE_ANGLE,
        double SOIL_DEPTH, double OV_PAR_1, double OV_PAR_2, double EFF_POR, double KSAT_0, double A)
{
    int i;
    int acc_num_new = 0;                   /*  Number of characteristic curves in lower layer that has reached the end of the hillslope during current time step  */
    int num_end = 0;                       /*  Number of characteristic curves in lower or upper layer layer that reaches the downslope end during current time step  */
    int new_high = 0;                      /*  Number of new characteristic curves in upper soil layer  */

    double sum_end = 0;                     /*  Accumulated discharge from characteristic curves in lower or upper layer layer that reaches
											the downslope end of hillslope during current time step  */
    double remain_distance;                 /*  Remaining distance to downslope end of hillslope  */
    double remain_time;                     /*  Remaining time to downslope end of hillslope  */
    double downslope_depth;                 /*  Saturated depth at downslope end of hillslope in lower or upper layer  */
    double increment_dist;                  /*  Distance increment between new characteristic curves in lower layer  */
    double increment_depth;                 /*  Saturated depth increment between new characteristic curves in lower layer  */
    double inter_time[numberCharacteristic + 1];          /*  Time from start of time step to intersection of characteristic curve in lower layer
														  with level SOIL_DEPTH  */
    double inter_length[numberCharacteristic + 1];        /*  Distance from initial position of characteristic curve in lower layer
														  in current time step to intersection with level SOIL_DEPTH  */
    double timeResolution = 1.0;

    for (i = 1; i <= numberCharacteristic; i++)
    {
        inter_time[i] = -999;
        inter_length[i] = -999;
    }

    /* Calculate time and length from start of time step to intersection of characteristic curve in lower layer with level SOIL_DEPTH  */
    i = numberCharacteristic + 1;
    while (i > 1 && (*u_high + new_high < numberCharacteristic))
    {
        i--;
        if (*u_high > 0 && len_coord[i] >= upp_len[*u_high])
        {
            inter_time[i] = 0;
            inter_length[i] = len_coord[i];
        }
        else if (sat_depth[i] >= SOIL_DEPTH)
        {
            new_high++;
            inter_time[i] = 0;
            inter_length[i] = len_coord[i];
            upp_tim[*u_high + new_high] = 0;
            upp_len[*u_high + new_high] = len_coord[i];
            upp_dep[*u_high + new_high] = 0;
        }
        else
        {
            inter_time[i] = (SOIL_DEPTH - sat_depth[i]) * EFF_POR / inflow;
            inter_length[i] = len_coord[i] + ((KSAT_0 * sin(SLOPE_ANGLE)) / (A * inflow)) * exp(A * (SOIL_DEPTH - sat_depth[i]))
                              * (1.0 - exp(((-1) * A * inflow / EFF_POR) * inter_time[i]));
            if (inter_time[i] < timeResolution && inter_length[i] < SLOPE_LENGTH)
            {
                new_high++;
                upp_tim[*u_high + new_high] = inter_time[i];
                upp_len[*u_high + new_high] = inter_length[i];
                upp_dep[*u_high + new_high] = 0;
            }
            else
            {
                inter_time[i] = -999;
                inter_length[i] = -999;
            }
        }
    }

    /* New values are assigned to index u_high */
    *u_high = *u_high + new_high;
    if (*u_high > numberCharacteristic)
    {
        printf("    u_high = %d    new_high = %d\n", *u_high, new_high);
        /*    exit(1);*/
    }

    /* Calculate time of arrival at downslope end for each characteristic curve in lower layer.
    If time of arrival is within current time step, calculate discharge at downslope end and assign
    the average of these values to upper layer discharge from hillslope for current time step. */
    for (i = numberCharacteristic; i >= 1; i--)
    {
        if (inter_time[i] >= 0)
        {
            remain_distance = SLOPE_LENGTH - inter_length[i];
            remain_time = inter_time[i] + EFF_POR * exp(A * (SOIL_DEPTH - SOIL_DEPTH)) * remain_distance / (KSAT_0 * sin(SLOPE_ANGLE));
            downslope_depth = SOIL_DEPTH;
        }
        else
        {
            remain_distance = SLOPE_LENGTH - len_coord[i];
            remain_time = EFF_POR * (1.0 / (A * inflow)) * (A * (SOIL_DEPTH - sat_depth[i]) -
                          log(exp(A * (SOIL_DEPTH - sat_depth[i])) - A * inflow * remain_distance / (KSAT_0 * sin(SLOPE_ANGLE))));
            downslope_depth = sat_depth[i] + inflow * remain_time / EFF_POR;
        }
        if (remain_time <= timeResolution)
        {
            num_end = num_end + 1;                                /* New characteristic curve in lower layer */
            sum_end = sum_end + KSAT_0 * sin(SLOPE_ANGLE) * (1.0 / A) * exp(A * SOIL_DEPTH) * (1 - exp((-1) * A * downslope_depth));
        }
    }
    if (num_end > 0)
    {
        *lower_runoff = sum_end / num_end;
    }
    else
    {
        *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0 / A) * exp(A * SOIL_DEPTH) * (1 - exp((-1) * A * sat_depth[numberCharacteristic]));
    }
    *lower_runoff = *lower_runoff / SLOPE_LENGTH;                           /* Convert m2/time_step to m/time_step */
    /*  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                           * Convert m2/time_step to mm/time_step */
    /*  *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);        * Convert m2/time_step to m3/second */

    /* Update number of new characteristic curves in lower layer to be started */
    acc_num_new = num_end;
    if (acc_num_new > numberCharacteristic)
    {
        printf("    ERROR    acc_num_new = %d    numberCharacteristic = %d\n", acc_num_new, numberCharacteristic);
        /*    exit(1);*/
    }

    /* Calculate length coordinates along all characteristic curves in lower layer that will not reach end of hillslope during current time step  */
    for (i = numberCharacteristic - acc_num_new; i >= 1; i--)
    {
        if (inter_time[i] >= 0)
        {
            len_coord[i] = inter_length[i] + KSAT_0 * sin(SLOPE_ANGLE) * (1.0 / EFF_POR) * exp(A * (SOIL_DEPTH - SOIL_DEPTH)) * (timeResolution - inter_time[i]);
            sat_depth[i] = SOIL_DEPTH;
        }
        else
        {
            len_coord[i] = len_coord[i] + KSAT_0 * sin(SLOPE_ANGLE) * (1.0 / (A * inflow)) * exp(A * (SOIL_DEPTH - sat_depth[i])) *
                           (1 - exp((-1) * A * inflow * timeResolution / EFF_POR));
            sat_depth[i] = sat_depth[i] + inflow * timeResolution / EFF_POR;
        }
    }

    if (acc_num_new > 0)
    {
        /* Update indices for each characteristic curve in lower layer that has not reached the downslope end */
        for (i = numberCharacteristic - acc_num_new; i >= 1; i--)
        {
            len_coord[i + acc_num_new] = len_coord[i];
            sat_depth[i + acc_num_new] = sat_depth[i];
            inter_time[i + acc_num_new] = inter_time[i];
            inter_length[i + acc_num_new] = inter_length[i];
        }

        /* For each characteristic curve in lower layer that has reached the downslope end, start a new characteristic curve at the upslope end */
        increment_dist = len_coord[1] / (acc_num_new + 1);
        increment_depth = sat_depth[1] / (acc_num_new + 1);
        for (i = 1; i <= acc_num_new; i++)
        {
            len_coord[i] = increment_dist * i;
            sat_depth[i] = increment_depth * i;
        }
        acc_num_new = 0;
    }

    /* Calculate discharge from lower layer */
    /*  i = numberCharacteristic;
    *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*sat_depth[i]));
    *  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step *
    *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */

    if (*u_high > 0)
    {
        /* Calculate time of arrival at downslope end for each characteristic curve in the upper layer.
        If time of arrival is within current time step, calculate discharge at downslope end and assign
        the average of these values to upper layer discharge from hillslope for current time step. */
        num_end = 0;
        sum_end = 0;
        for (i = 1; i <= *u_high; i++)
        {
            remain_distance = SLOPE_LENGTH - upp_len[i];
            remain_time = (1 / inflow) * power((inflow * remain_distance / OV_PAR_1 + power(upp_dep[i], OV_PAR_2)), (1 / OV_PAR_2)) - (upp_dep[i] / inflow);
            downslope_depth = upp_dep[i] + inflow * remain_time;
            if (remain_time <= timeResolution)
            {
                num_end = num_end + 1;
                sum_end = sum_end + OV_PAR_1 * power(downslope_depth, OV_PAR_2) / SLOPE_LENGTH;                       /* Convert m2/time_step to m/time_step */
                /*      sum_end =sum_end + OV_PAR_1 * power(downslope_depth,OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
                /*        sum_end =sum_end + OV_PAR_1 * power(downslope_depth,OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */
            }
        }
        if (num_end > 0)
        {
            *upper_runoff = sum_end / num_end;
        }
        else
        {
            *upper_runoff = OV_PAR_1 * power(upp_dep[1], OV_PAR_2) / SLOPE_LENGTH;    /* Convert m2/time_step to m/time_step */
        }
        /*      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
        /*      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */

        /* Update indices for each characteristic curve in the upper layer that has not reached the downslope end and index u_high */
        i = 1;
        while (i + num_end <= *u_high)
        {
            upp_tim[i] = upp_tim[i + num_end];
            upp_len[i] = upp_len[i + num_end];
            upp_dep[i] = upp_dep[i + num_end];
            i++;
        }
        *u_high = *u_high - num_end;
        if (*u_high < 0)
        {
            *u_high = 0;
        }

        /* Calculate length coordinates along all characteristic curves in the upper layer that
        will not reach the end of the hillslope during the current time step  */
        for (i = 1; i <= *u_high; i++)
        {
            upp_len[i] = upp_len[i] + (OV_PAR_1 / (inflow)) * (power((upp_dep[i] + inflow * (timeResolution - upp_tim[i])), OV_PAR_2) - power(upp_dep[i], OV_PAR_2));
            upp_dep[i] = upp_dep[i] + inflow * (timeResolution - upp_tim[i]);
            upp_tim[i] = 0;
            if (upp_len[i] > SLOPE_LENGTH)
            {
                /*      printf("\n    upp_len[i] = %f    i = %d\n\n",upp_len[i],i);*/
                upp_len[i] = SLOPE_LENGTH;
            }
        }
    }

    return;
}


void KinematicWave::kinematic_wave_without_lateral_inflow(double *len_coord, double *sat_depth, double *lower_runoff, double *upper_runoff,
        double *upp_tim, double *upp_len, double *upp_dep, int *u_high, double SLOPE_LENGTH, double SLOPE_ANGLE,
        double SOIL_DEPTH, double OV_PAR_1, double OV_PAR_2, double EFF_POR, double KSAT_0, double A)
{
    int i;
    int acc_num_new = 0;                   /*  Number of characteristic curves in lower layer that has reached the end of the hillslope during current time step  */
    int num_end = 0;                       /*  Number of characteristic curves in lower or upper layer that reaches the downslope end during current time step  */

    double sum_end = 0;                     /*  Accumulated discharge from characteristic curves in lower or upper layer layer that reaches
											the downslope end of hillslope during current time step  */
    double remain_distance;                 /*  Remaining distance to downslope end of hillslope  */
    double remain_time;                     /*  Remaining time to downslope end of hillslope  */
    double increment_dist;                  /*  Distance increment between new characteristic curves in lower layer  */
    double increment_depth;                 /*  Saturated depth increment between new characteristic curves in lower layer  */
    double timeResolution = 1.0;

    /* Calculate time of arrival at downslope end for each characteristic curve in lower layer.
    If time of arrival is within current time step, calculate discharge at downslope end and assign
    the average of these values to upper layer discharge from hillslope for current time step. */
    for (i = numberCharacteristic; i >= 1; i--)
    {
        remain_distance = SLOPE_LENGTH - len_coord[i];
        remain_time = EFF_POR * exp(A * (sat_depth[i] - SOIL_DEPTH)) * remain_distance / (KSAT_0 * sin(SLOPE_ANGLE));
        if (remain_time <= timeResolution)
        {
            num_end = num_end + 1;                                /* New characteristic curve in lower layer */
            sum_end = sum_end + KSAT_0 * sin(SLOPE_ANGLE) * (1.0 / A) * exp(A * SOIL_DEPTH) * (1 - exp((-1) * A * sat_depth[i]));
        }
    }
    if (num_end > 0)
    {
        *lower_runoff = sum_end / num_end;
    }
    else
    {
        *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0 / A) * exp(A * SOIL_DEPTH) * (1 - exp((-1) * A * sat_depth[numberCharacteristic]));
    }
    *lower_runoff = *lower_runoff / SLOPE_LENGTH;                           /* Convert m2/time_step to m/time_step */
    /*  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                           * Convert m2/time_step to mm/time_step */
    /*  *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);        * Convert m2/time_step to m3/second */

    /* Update number of new characteristic curves in lower layer to be started */
    acc_num_new = num_end;
    if (acc_num_new > numberCharacteristic)
    {
        printf("    ERROR    acc_num_new = %d    numberCharacteristic = %d\n", acc_num_new, numberCharacteristic);
        /*    exit(1);*/
    }

    /* Calculate length coordinates along all characteristic curves in lower layer that will not reach end of hillslope during current time step.
    Saturated thickness remains unchanged in the absence of lateral inflow */
    for (i = numberCharacteristic - acc_num_new; i >= 1; i--)
    {
        len_coord[i] = len_coord[i] + KSAT_0 * sin(SLOPE_ANGLE) * (1.0 / EFF_POR) * exp(A * (SOIL_DEPTH - sat_depth[i])) * timeResolution;
    }

    if (acc_num_new > 0)
    {
        /* Update indices for each characteristic curve in lower layer that has not reached the downslope end */
        for (i = numberCharacteristic - acc_num_new; i >= 1; i--)
        {
            len_coord[i + acc_num_new] = len_coord[i];
            sat_depth[i + acc_num_new] = sat_depth[i];
        }

        /* For each characteristic curve in lower layer that has reached the downslope end, start a new characteristic curve at the upslope end */
        increment_dist = len_coord[1] / (acc_num_new + 1);
        increment_depth = sat_depth[1] / (acc_num_new + 1);
        for (i = 1; i <= acc_num_new; i++)
        {
            len_coord[i] = increment_dist * i;
            sat_depth[i] = increment_depth * i;
        }
        acc_num_new = 0;
    }

    /* Calculate discharge from lower layer */
    /*  i = numberCharacteristic;
    *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*sat_depth[i]));
    *  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                           * Convert m2/time_step to mm/time_step *
    *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);        * Convert m2/time_step to m3/second */

    if (*u_high > 0)
    {
        /* Calculate time of arrival at downslope end for each characteristic curve in the upper layer.
        If time of arrival is within current time step, calculate discharge at downslope end and assign
        the average of these values to upper layer discharge from hillslope for current time step. */
        num_end = 0;
        sum_end = 0;
        for (i = 1; i <= *u_high; i++)
        {
            remain_distance = SLOPE_LENGTH - upp_len[i];
            remain_time = remain_distance / (OV_PAR_1 * OV_PAR_2 * power(upp_dep[i], OV_PAR_2 - 1.0));
            if (remain_time <= timeResolution)
            {
                num_end = num_end + 1;
                sum_end = sum_end + OV_PAR_1 * power(upp_dep[i], OV_PAR_2) / SLOPE_LENGTH;                       /* Convert m2/time_step to m/time_step */
                /*      sum_end =sum_end + OV_PAR_1 * power(upp_dep[i],OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
                /*        sum_end =sum_end + OV_PAR_1 * power(upp_dep[i],OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */
            }
        }
        if (num_end > 0)
        {
            *upper_runoff = sum_end / num_end;
        }
        else
        {
            *upper_runoff = OV_PAR_1 * power(upp_dep[1], OV_PAR_2) / SLOPE_LENGTH;    /* Convert m2/time_step to m/time_step */
        }
        /*      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
        /*      upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */

        /* Update indices for each characteristic curve in the upper layer that has not reached the downslope end and index u_high */
        i = 1;
        while (i + num_end <= *u_high)
        {
            upp_tim[i] = upp_tim[i + num_end];
            upp_len[i] = upp_len[i + num_end];
            upp_dep[i] = upp_dep[i + num_end];
            i++;
        }
        *u_high = *u_high - num_end;
        if (*u_high < 0)
        {
            *u_high = 0;
        }

        /* Calculate length coordinates along all characteristic curves in the upper layer that
        will not reach the end of the hillslope during the current time step
        Depth of upper layer remains unchanged in the absence of lateral inflow */
        for (i = 1; i <= *u_high; i++)
        {
            upp_tim[i] = 0;
            upp_len[i] = upp_len[i] + OV_PAR_1 * OV_PAR_2 * power(upp_dep[i], OV_PAR_2 - 1.0) * timeResolution;
            if (upp_len[i] > SLOPE_LENGTH)
            {
                /*      printf("\n    upp_len[i] = %f    i = %d\n\n",upp_len[i],i);*/
                upp_len[i] = SLOPE_LENGTH;
            }
        }
    }

    return;
}


/* Groundwater table depth at a fraction of total hillslope length (from top to bottom) */
double KinematicWave::GetGroundWaterDepth(double lengthFraction) const
{
    return -(kiWaPar->GetSOIL_DEPTH() - sat_depth[(int)(numberCharacteristic * lengthFraction)]);
}


/* Volumetric soil moisture content at a fraction of total hillslope length (from top to bottom) */
double KinematicWave::GetSoilMoisture(double lengthFraction) const
{
    double gw_h, field_capacity, teta;
    gw_h = (kiWaPar->GetSOIL_DEPTH() - sat_depth[(int)(numberCharacteristic * lengthFraction)]);
    field_capacity = kiWaPar->GetTSAT_0() * exp(kiWaPar->GetLAMBDA_KW() * gw_h);
    teta = field_capacity + smdef[(int)(numberCharacteristic * lengthFraction)] / kiWaPar->GetROOT_DEPTH();
    return teta;
}

