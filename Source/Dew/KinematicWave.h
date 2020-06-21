#pragma once

using namespace std;
class DateTime;
class ParametersGeneral;
class ParametersLandSurface;
class ParametersKiWa;
class InputElement;
class DistributedElement;
class ModelControl;

#include "Dew.h"

class KinematicWave
{
public:
    KinematicWave();
    ~KinematicWave();
    void SetInitialKwaValues();
    void WaterBalance(int timeStep, DateTime datetime, double waterInput, double temperature, double snowCoverFraction, double dry_period, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge);
    void kinematic_wave_with_lateral_inflow(double, double *, double *, double *, double *, double *, double *, double *, int *,
                                            double, double, double, double, double, double, double, double);
    void kinematic_wave_without_lateral_inflow(double *, double *, double *, double *, double *, double *, double *, int *,
            double, double, double, double, double, double, double, double);
    void KiWaGroundWaterTable(int elementId, int timeStep);
    void KiWaSoilMoisture(int elementId, int timeStep);
    double GetSoilMoisture(double lengthFraction) const;
    double GetGroundWaterDepth(double lengthFraction) const;
    double GetTranspSoilEvap() const;
    double GetLowerRunoff() const;
    double GetUpperRunoff() const;
    double GetRunoff() const;
    void SetGeneralPar(ParametersGeneral *parObj);
    ParametersGeneral *GetGeneralPar() const;
    void SetLandSurfacePar(ParametersLandSurface *parObj);
    ParametersLandSurface *GetLandSurfacePar();
    void SetKiWaPar(ParametersKiWa *parObj);
    ParametersKiWa *GetKiWaPar() const;
    //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj);
    //  InputTimeSeries *GetInputTimeSeries() const;
    void SetInputElement(InputElement *inElementObj);
    InputElement *GetInputElement() const;
    void SetLandScapeElement(DistributedElement *theElement);
    DistributedElement *GetLandScapeElement() const;

private:
    ParametersGeneral *commonPar;
    ParametersKiWa *kiWaPar;
    ParametersLandSurface *landSurfacePar;
    //  InputTimeSeries *inTimeSeries;
    InputElement *inElement;
    DistributedElement *landScapeElement;
    int u_high;                                 /*  Number of characteristic curves in upper layer  */
    //  double temp;                                /*  Temperature (deg. C)  */
    double vap_def;                             /*  Vapor pressure deficit (hPa)  */
    double mean_perc;                           /*  Mean percolation along hillslope (m/time_step)  */
    double field_capacity;                      /*  Field capacity  */
    double teta;                                /*  Volumetric water content */
    double gw_h;                                /*  Depth to groundwater table (m)  */
    double evaporated_volume;                   /*  Volume of water lost from lower layer by evapotranspiration (m)  */
    double evol_upper;                          /*  Volume of water lost from upper layer by evapotranspiration (m)  */
    double volume_root;                         /*  Volume of water removed from root zone by evapotranspiration (m)  */
    double def_par;                             /*  Fraction of actual evapotranspiration removed from soil moisture  */
    double smdef[numberCharacteristic + 1];       /*  Soil moisture deficit (m) ( <= 0 )  */
    double perc[numberCharacteristic + 1];        /*  Volume of water percolating to saturated zone (m)  */
    double len_coord[numberCharacteristic + 1];   /*  Length-coordinate along characteristic curves in lower layer (m)  */
    double fixed_length[numberCharacteristic + 1];/*  Fixed length-coordinate along hillslope (m)  */
    double sat_depth[numberCharacteristic + 1];   /*  Depth of saturated zone along characteristic curves in lower layer (m)  */
    double fixed_sat[numberCharacteristic + 1];   /*  Saturated depth at fixed length coordinates  */
    double upp_tim[numberCharacteristic + 1];     /*  Initial time within time step of characteristic curves in upper layer  */
    double upp_dep[numberCharacteristic + 1];     /*  Depth along characteristic curves in upper layer  */
    double upp_len[numberCharacteristic + 1];     /*  Length coordinate along characteristic curves in upper layer  */
    double actev;                               /*  Actual evapotranspiration (m/time_step)  */
    double actev_loss;                          /*  Actual evapotranspiration loss from soil (m)  */
    double transpSoilEvap;                      /*  Water lost from subsurface by evapotranspiration (m)  */
    double lower_runoff;                        /*  Runoff from lower layer (mm/time_step)  */
    double upper_runoff;                        /*  Runoff from upper layer (mm/time_step)  */
    double runoff;                              /*  Runoff (m/timestep)  */
};
