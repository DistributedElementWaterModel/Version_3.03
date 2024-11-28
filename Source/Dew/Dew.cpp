/****************************************************************************************
*                                                                                      *
*  Program for calculating water balance for distributed hydrological response units,  *
*  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
*  river networks, lakes and landscape elements.                                       *
*                                                                                      *
*  Author: Stein Beldring, Norway                                                      *
*                                                                                      *
****************************************************************************************/

#include "stdafx.h"

using namespace std;

#include "DateTime.h"
#include "Date.h"
#include "CTime.h"
#include "DistributedElement.h"
#include "SubCatchment.h"
#include "ParametersGeneral.h"
#include "MeteorologicalStations.h"
#include "GlacierElements.h"
#include "InputTimeSeries.h"
#include "InputElement.h"
#include "ModelControl.h"
#include "SelectedHbvTimeSeriesElements.h"
#include "SelectedKiWaTimeSeriesElements.h"
#include "SelectedSubCatchmentTimeSeriesElements.h"
#include "TotalReservoirStorage.h"
//#include "Albedo.h"
//#include "Atmosphere.h"
#include "Lake.h"
#include "Glacier.h"
#include "GlacierIce.h"
#include "HBV.h"
#include "HbvAquifer.h"
#include "KinematicWave.h"
#include "LakeWaterBalance.h"
#include "KwaAquifer.h"
#include "Snow.h"
#include "Vegetation.h"
#include "ParametersGlacierRetreat.h"
#include "ParametersKiWa.h"
#include "ParametersLandSurface.h"
#include "ParametersSubSurfaceHbv.h"
//#include "LandAtmosphereInterface.h"
//#include "PenmanMonteith.h"
#include "Util.h"
#include "Routing.h"
#include "Parameters.h"

void ReadSubCatchmentIdentifier(DistributedElement * const Dew, SubCatchment * const CatchmentElement,
                                int numWatc, bool flowHierarchy, ifstream &fileControl, ofstream &fout);
void ReadLandscapeHierarchy(DistributedElement * const Dew, ifstream &fileControl, ofstream &fout);
void SnowGlacierIceReDistribution(SubCatchment ** const Outlet, DistributedElement * const Dew, ParametersGeneral * ParGeneralStore, int initialTimeSteps, int numberTimeSteps,
                                  int numLand, int numWatcOut, int timeStep, DateTime datetime, char glacierModelling, GlacierElements * DewGlacierElements,
                                  int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, bool modelCalibration, ofstream &fout);
void WaterBalanceTimeSeries(DistributedElement * const Dew, ParametersGeneral * const ParGeneralStore,
                            MeteorologicalStations * const MetStations, int initialTimeSteps, int numberTimeSteps,
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore,
                            int numLand, int timeStep, DateTime datetime, bool * inputDataFound, bool * inputDataError, bool * firstTotal);
void SetTimeSeriesInput(DistributedElement * const Dew, ParametersGeneral * const ParGeneralStore,
                        MeteorologicalStations * const MetStations, int initialTimeSteps, int numberTimeSteps,
                        InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore,
                        int timeStep, bool * inputDataFound, bool * inputDataError);
void WaterBalanceGrid(DistributedElement * const Dew, ParametersGeneral * ParGeneralStore,
                      InputElement * InputElementStore, int initialTimeSteps, int numberTimeSteps,
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy,
                      bool * inputDataFound, bool * firstTotal);
void WaterBalanceGridVariable(DistributedElement * const Dew, ParametersGeneral * ParGeneralStore,
                      InputElement * InputElementStore, int initialTimeSteps, int numberTimeSteps,
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy,
                      bool * inputDataFound, bool * firstTotal);
void WaterBalanceElements(DistributedElement * const Dew, ParametersGeneral * ParGeneralStore,
                      InputElement * InputElementStore, int initialTimeSteps, int numberTimeSteps,
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy,
                      bool * inputDataFound, bool * firstTotal);
void SetElementCorrection(DistributedElement * const Dew, int numLand, int numberCorrectionCatchments, int * correctionCatchments, 
			  double * correctionPrecipitation, double * correctionTemperature, ofstream &fout);
void TraverseOutletRoutingSubCatchment(SubCatchment * const thisSubCatchment, int numWatc, int * outletRoutingSubCatchments);
void TraverseCorrectionSubCatchment(SubCatchment * const thisSubCatchment, int numberCorrectionCatchments,
                                    int * correctionCatchments, double * correctionPrecipitation,
                                    double * correctionTemperature);
void TraverseCorrectionLandScape(DistributedElement * const thisElement, double precCorr, double tempCorr);
void TraverseSubCatchment(SubCatchment * const thisSubCatchment, ParametersGeneral * const ParGeneralStore,
                          MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
                          InputElement * InputElementStore, int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, char inputFormat,
                          bool flowHierarchy, bool forceDirect, bool * inputDataFound, bool * inputDataError, bool * firstTotal, ofstream &fout);
void TraverseLandScape(DistributedElement * const thisElement, ParametersGeneral * const ParGeneralStore,
                       MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
                       InputElement * InputElementStore, int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, char inputFormat,
                       bool flowHierarchy, bool forceDirect, bool * inputDataFound, bool * inputDataError, bool * firstTotal, ofstream &fout);
void TraverseMissingDataSubCatchment(SubCatchment * const thisSubCatchment, int timeStep, ofstream &fout);
void TraverseMissingDataLandScape(DistributedElement * const thisElement, int timeStep, ofstream &fout);
void WriteSubCatchmentIdentifier(SubCatchment * const CatchmentElement, int numWatc, ofstream &fout);
void WriteSubCatchmentDischarge(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime,
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                int secondsPerTimeStep, bool modelCalibration, ofstream &fout);
void WriteSubCatchmentWaterBalance(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime,
                                   DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                   int secondsPerTimeStep);
void WriteOutletDischarge(SubCatchment * const thisSubCatchment, DateTime startSimulationTime,
                          DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                          int secondsPerTimeStep, bool modelCalibration, ofstream &fout);
void WriteOutletWaterBalance(SubCatchment * const thisSubCatchment, DateTime startSimulationTime,
                             DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                             int secondsPerTimeStep, bool flowHierarchy, bool forceDirect);
void WriteReducedBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows, int nCols,
                            int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void WriteBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows, int nCols,
                     int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void WriteAsciiGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows, int nCols,
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void WriteAsciiGridWaterBalance(DistributedElement * const Dew, DateTime startSimulationTime, DateTime endSimulationTime, int numLand, int nRows, int nCols, 
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void WriteAsciiGridAnnualGlacierValues(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,
                                       int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void ObjectiveCriteria(int numberTimeSteps, double *obs, double *sim,
                       double *ns, double *rmse, double *bias, double *pears);
void WriteDistributedElementTimeSeries(DistributedElement * const Dew, int numLand,
                                       DateTime startSimulationTime, DateTime endSimulationTime,
                                       int initialTimeSteps, int numberTimeSteps, int secondsPerTimeStep);
void WriteTotalReservoirStorage(TotalReservoirStorage * const TotalReservoirStore, DateTime startSimulationTime,
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                bool modelCalibration, int secondsPerTimeStep);
void GetRowCol(int elementNo, int nCols, int &row, int &col);

int main(int argc, char *argv[])
{
    std::cout << "\n\n Distributed Element Water Model \n\n";

    bool inputDataFound = false;
    bool inputDataError = false;
    bool modelCalibration = false;
    bool flowHierarchy = false;
    bool forceDirect = false;
    bool firstTotal = false;
    bool maxBasTransform=false;
    char fileName[100];
    char fileObsStreamflow[100];
    char fileNameInput[100];
    char fileNameSubCatchment[100];
    char fileNameRouting[100];
    char buffer[1024];
    char manningVegetationType[100];
    char modelRun = 'M';
    char hierarchy = 'H';
    char inputFormat = 'F';
    char evaporationModelling = 'T';
    char glacierModelling = 'G';
    char routingType = 'R';
    char ch;
    char *metPath;
    char *metMask;
    int g, h, i, j, k;
    int numLand;
    int geoIndex;
    int landIndex;
    int nRows, nCols, noData;
    int modStruct;
    int landSurf;
    int soil;
    int elementFlowDir;
    int numberTimeSteps, initialTimeSteps, timeStep;
    int numWatc, numWatcUp, numWatcOut;
    int subCatchmentId;
    int day, mth, year, hour, minute;
    //  int time, startModelTime, startSimulationTime, endSimulationTime;
    int numberCorrectionCatchments;
    int seriesNumber;
    int lakeNumber;
    int velocityIndex, numWatcRoute, currentCat, manningIndex;
    double seriesWeight;
    double sumWeight;
    double sumArea;
    double precStationsWeightedElevation, tempStationsWeightedElevation;
    double correction, maxBas;
    double travelTime, riverSlope, width;
    double xllCorner, yllCorner, cellSize;
    double elementArea, elementElevation, elementSlopeLength, elementSlopeAngle, elementAspect, lakePercent, glacierPercent, glacSurfElev, glacIceThick;
    double areaFraction[maximumNumberLandClasses];
    double manningRoughness[maximumNumberManningClasses];
    DateTime datetime;
    Lake *ptrLake;
    Glacier *ptrGlacier;
    HbvAquifer *ptrHbvAquifer;
    HbvAquifer *lastHbvAquifer;
    KwaAquifer *ptrKwaAquifer;
    KwaAquifer *lastKwaAquifer;
    LakeWaterBalance *ptrLakeWaterBalance;
    Vegetation *ptrVegetation;
    Snow *ptrSnow;
    //  GlacierSurface * ptrSurface;
    GlacierIce *ptrIce;
    HBV *ptrHbv;
    KinematicWave *ptrKwa;
    //    LandAtmosphereInterface *ptrLandAtmInt;
    //    Albedo *ptrAlbedo;
    //    Atmosphere *ptrAtm;
    //    PenmanMonteith *ptrPenMon;
    MODEL_STRUCTURE modelStructure;
    LANDSURFACE landSurfType[numberLandSurfaceClasses];
    SOIL soilType[numberSoilClasses];
    GLACIER_TYPE glacierType[numberGlacierClasses];

    if (argc != 2)
    {
        cout << " " << argv[0] << "  <control file name>\n\n";
        exit(1);
    }
    ifstream fileControl(argv[1]);
    if (!fileControl.is_open())
    {
        cout << " Error opening file " << argv[1] << endl << endl;
        exit(1);
    }

    /*  while (modelRun != 'S' && modelRun != 'C' && modelRun != 's' && modelRun != 'c') {
    cout << " Type of model run, simulation(S) or calibration(C): ";
    cin >> modelRun;
    cout << endl;
    }*/
    fileControl.ignore(100, ':');
    fileControl >> modelRun;
    fileControl.ignore(1024, '\n');
    if (modelRun != 'S' && modelRun != 'C' && modelRun != 's' && modelRun != 'c')
    {
        cout << "\n Type of model run, simulation(S) or calibration(C) \n\n";
        exit(1);
    }
    if (modelRun == 'C' || modelRun == 'c')
    {
        modelCalibration = true;
    }

    /*  while (hierarchy != 'N' && hierarchy != 'C' && hierarchy != 'n' && hierarchy != 'c') {
    cout << " Landscape elements hierarchy, flow direction network(N) or nested catchments(C): ";
    cin >> hierarchy;
    cout << endl;
    }*/
    fileControl.ignore(100, ':');
    fileControl >> hierarchy;
    fileControl.ignore(1024, '\n');
    if (hierarchy != 'N' && hierarchy != 'C' && hierarchy != 'n' && hierarchy != 'c')
    {
        cout << "\n Landscape elements hierarchy, flow direction network(N) or nested catchments(C) \n\n";
        exit(1);
    }
    if (hierarchy == 'N' || hierarchy == 'n')
    {
        flowHierarchy = true;
    }

    /*  while (inputFormat != 'G' && inputFormat != 'T' && inputFormat != 'g' && inputFormat != 't') {
    cout << " Input data format, grid files(G) or time series file(T): ";
    cin >> inputFormat;
    cout << endl;
    }*/
    fileControl.ignore(100, ':');
    fileControl >> inputFormat;
    fileControl.ignore(1024, '\n');
    if (inputFormat != 'G' && inputFormat != 'T' && inputFormat != 'g' && inputFormat != 't')
    {
        cout << "\n Input data format, grid files(G) or time series file(T) \n\n";
        exit(1);
    }

    /*  while (evaporationModelling != 'T' && evaporationModelling != 'M' && evaporationModelling != 't' && evaporationModelling != 't') {
        cout << " Potential evaporation, temperature index (T) or long-term mean monthly values (M): ";
        cin >> evaporationModelling;
        cout << endl;
        }*/
    fileControl.ignore(256,':');
    fileControl >> evaporationModelling;
    fileControl.ignore(1024,'\n');
    if (evaporationModelling != 'T' && evaporationModelling != 'M' && evaporationModelling != 't' && evaporationModelling != 'm') {
      cout << "\n Potential evaporation, temperature index (T) or long-term mean monthly values (M) \n\n";
      exit (1);
    }

    /*  while (evaporationModelling != 'T' && evaporationModelling != 'P' && evaporationModelling != 't' && evaporationModelling != 'p') {
    cout << " Evaporation modelling, temperature index (T) or Penman-Monteith parameterization (P): ";
    cin >> evaporationModelling;
    cout << endl;
    }*/
    /*  fileControl.ignore(100,':');
    fileControl >> evaporationModelling;
    fileControl.ignore(1024,'\n');
    if (evaporationModelling != 'T' && evaporationModelling != 'P' && evaporationModelling != 't' && evaporationModelling != 'p') {
    cout << "\n Evaporation modelling, temperature index (T) or Penman-Monteith parameterization (P) \n\n";
    exit (1);
    }*/

    /*  while (glacierModelling != 'S' && glacierModelling != 'E' && glacierModelling != 'M' && glacierModelling != 's' && glacierModelling != 'e' && glacierModelling != 'm') {
    cout << " Glacier modelling, static ice (S), elevation based parameterisation of ice dynamics (E) or surface elevation change based on mass balance (M): ";
    cin >> glacierModelling;
    cout << endl;
    }*/
    fileControl.ignore(100, ':');
    fileControl >> glacierModelling;
    fileControl.ignore(1024, '\n');
    if (glacierModelling != 'S' && glacierModelling != 'E' && glacierModelling != 'M' && glacierModelling != 's' && glacierModelling != 'e' && glacierModelling != 'm')
    {
        cout << "\n Glacier modelling, static ice (S), elevation based parameterisation of ice dynamics (E) or surface elevation change based on mass balance (M) \n\n";
        exit(1);
    }

    /*  cout << " Routing type, No routing(0), Muskingum-Cunge(M), Source-To-Sink(S) or Triangular distribution of weights(T): ";
    cin >> routingType;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> routingType;
    fileControl.ignore(1024, '\n');
    if (routingType != 'M' && routingType != 'S' && routingType != 'T' && routingType != 'm' && routingType != 's' && routingType != 't' && routingType != '0')
    {
        cout << "\n Routing type, No routing(0), Muskingum-Cunge(M) or Source-To-Sink(S) \n\n";
        exit(1);
    }
    if (routingType == 'S' || routingType == 's')
    {
        forceDirect = true;
    }
    if (routingType == 'T' || routingType == 't') 
    {
      maxBasTransform = true;
    }

    /*  cout << " Output file: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ofstream fout(fileName);  // Open for writing

    // Read dates for start model, start simulation, end simulation
    /*  cout << " Start model date and time (day, month, year, hour, minute): ";
    cin >> day >> mth >> year >> hour >> minute;*/
    fileControl.ignore(100, ':');
    fileControl >> day >> mth >> year >> hour >> minute;
    fileControl.ignore(1024, '\n');
    DateTime startModelTime(year, mth, day, hour, minute, 0);
    if (startModelTime.legal() != 1)
    {
        cout << endl << " Not legal start model date and time " << endl << endl;
        exit(1);
    }
    /*  cout << endl;
    cout << " Start simulation date and time (day, month, year, hour, minute): ";
    cin >> day >> mth >> year >> hour >> minute;*/
    fileControl.ignore(100, ':');
    fileControl >> day >> mth >> year >> hour >> minute;
    fileControl.ignore(1024, '\n');
    DateTime startSimulationTime(year, mth, day, hour, minute, 0);
    if (startSimulationTime.legal() != 1)
    {
        cout << endl << " Not legal start simulation date and time " << endl << endl;
        exit(1);
    }
    /*  cout << endl;
    cout << " End simulation date and time (day, month, year, hour, minute): ";
    cin >> day >> mth >> year >> hour >> minute;*/
    fileControl.ignore(100, ':');
    fileControl >> day >> mth >> year >> hour >> minute;
    fileControl.ignore(1024, '\n');
    DateTime endSimulationTime(year, mth, day, hour, minute, 0);
    if (endSimulationTime.legal() != 1)
    {
        cout << endl << " Not legal end simulation date and time " << endl << endl;
        exit(1);
    }
    cout << endl;
    cout << "  " << startModelTime << "  " << startSimulationTime << "  " << endSimulationTime << endl;
    if (startModelTime > startSimulationTime || startSimulationTime > endSimulationTime)
    {
        cout << endl << " DateTime error " << endl << endl;
        exit(1);
    }

    // Object for storing model control information
    ModelControl * ModelControlObject = new ModelControl;
    ModelControlObject->SetModelControl(modelRun, hierarchy, inputFormat, evaporationModelling, glacierModelling, routingType);

    // Object for storing meteorological station information
    MeteorologicalStations * MetStations = new MeteorologicalStations;
    MetStations->SetMeteorologicalStations(fileControl, fout);

    // Read common parameters file and set parameter values
    ParametersGeneral * ParGeneralStore = new ParametersGeneral;
    SetGeneralParameters(ParGeneralStore, fileControl, fout);
    if ((ParGeneralStore->GetSECONDS_TIMESTEP() % minimumTimeStep) != 0)
    {
        cout << endl << " Temporal resolution must be a multiple of " << minimumTimeStep << " seconds : " << ParGeneralStore->GetSECONDS_TIMESTEP() << endl << endl;
        exit(1);
    }
    // End read common parameters

    // Read landsurface parameters file and set parameter values
    ParametersLandSurface * ParLandSurfaceStore = new ParametersLandSurface[numberLandSurfaceClasses];
    SetLandSurfaceParameters(ParLandSurfaceStore, fileControl, fout);
    // End read landsurface parameters

    // Read subsurface parameters file and set parameter values
    ParametersSubSurfaceHbv * ParSubSurfaceHbvStore = new ParametersSubSurfaceHbv[numberSoilClasses];
    SetSubSurfaceHbvParameters(ParSubSurfaceHbvStore, fileControl, fout);
    // End read subsurface parameters

    // Read kinematic wave parameters file and set parameter values
    ParametersKiWa * ParKiWaStore = new ParametersKiWa[numberSoilClasses];
    SetKiWaParameters(ParKiWaStore, fileControl, fout);
    // End read kinematic wave parameters

    // Read glacier retreat parameters file and set parameter values
    ParametersGlacierRetreat * ParGlacRetStore = new ParametersGlacierRetreat[numberGlacierClasses];
    SetGlacierRetreatParameters(ParGlacRetStore, fileControl, fout);
    // End read glacier retreat parameters

    // Object for storing landscape elements selected for HBV time series output
    SelectedHbvTimeSeriesElements * SelectedHbvTimeSeriesElementsStore = new SelectedHbvTimeSeriesElements;
    SelectedHbvTimeSeriesElementsStore->SelectedHbvTimeSeriesElementsInput(fileControl, fout);

    // Object for storing landscape elements selected for kinematic wave time series output
    SelectedKiWaTimeSeriesElements * SelectedKiWaTimeSeriesElementsStore = new SelectedKiWaTimeSeriesElements;
    SelectedKiWaTimeSeriesElementsStore->SelectedKiWaTimeSeriesElementsInput(fileControl, fout);

    // Object for storing subcatchment elements selected for time series output
    SelectedSubCatchmentTimeSeriesElements * SelectedSubCatchmentTimeSeriesElementsStore = new SelectedSubCatchmentTimeSeriesElements;
    SelectedSubCatchmentTimeSeriesElementsStore->SelectedSubCatchmentTimeSeriesElementsInput(fileControl, fout);

    // Object for storing input data for each landscape element
    InputElement * InputElementStore = new InputElement(numberInputSeries);

    // Calculate no. of initial and simulation time steps
    initialTimeSteps = (startSimulationTime.date2jday() - startModelTime.date2jday()) * (int)(numberSecondsDay/ParGeneralStore->GetSECONDS_TIMESTEP());
    numberTimeSteps = (1 + endSimulationTime.date2jday() - startSimulationTime.date2jday()) * (int)(numberSecondsDay/ParGeneralStore->GetSECONDS_TIMESTEP());
    //    initialTimeSteps = (int)((startSimulationTime - startModelTime) / (double)ParGeneralStore->GetSECONDS_TIMESTEP());
    //    numberTimeSteps = 1 + (int)((endSimulationTime - startSimulationTime) / (double)ParGeneralStore->GetSECONDS_TIMESTEP());
    cout << "  initialTimeSteps " << initialTimeSteps << "    numberTimeSteps " << numberTimeSteps << endl;

    // Define total reservoir storage object
    TotalReservoirStorage * TotalReservoirStore = new TotalReservoirStorage;
    TotalReservoirStore->AllocateTotalReservoirStorage(initialTimeSteps + numberTimeSteps);
    TotalReservoirStore->SetGeneralPar(ParGeneralStore);
    TotalReservoirStore->SetInitialTotalReservoirStorage();
    // End define total reservoir storage object

    // Read landscape element file and generate landscape element objects
    /*  cout << " File with landscape element information: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finLandScape(fileName);  // Open for reading
    if (!finLandScape.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        exit(1);
    }
    finLandScape >> buffer >> nCols;
    finLandScape >> buffer >> nRows;
    finLandScape >> buffer >> xllCorner;
    finLandScape >> buffer >> yllCorner;
    finLandScape >> buffer >> cellSize;
    finLandScape >> buffer >> noData;
    finLandScape.ignore(100, ':');
    finLandScape >> numLand;
    DistributedElement * Dew = new DistributedElement[numLand];
    if (!Dew)
    {
        cout << endl << " Error during dynamic memory allocating for Dew[" << numLand << "]" << endl << endl;
        exit(1);
    }

    for (i = 0; i < numLand; i++)
    {
        finLandScape >> landIndex >> geoIndex >> modStruct >> elementArea >> elementElevation;
        finLandScape >> elementSlopeLength >> elementSlopeAngle >> elementAspect >> elementFlowDir >> lakePercent >> glacierPercent >> glacSurfElev >> glacIceThick;
        modelStructure = MODEL_STRUCTURE(modStruct);
        //    cout << landIndex << "  " << geoIndex << "  " << elementArea << "  " << elementElevation << "  " << elementSlopeLength << "  " << elementSlopeAngle << "  " << elementAspect << "  " << elementFlowDir << "  " << lakePercent << "  " << glacierPercent << "  " << endl;
        sumArea = lakePercent + glacierPercent;
        for (j = 0; j < maximumNumberLandClasses; j++)
        {
            finLandScape >> landSurf >> soil >> areaFraction[j];
            landSurfType[j] = LANDSURFACE(landSurf);
            soilType[j] = SOIL(soil);
            sumArea = sumArea + areaFraction[j];
        }
        //    for (j=0; j<maximumNumberLandClasses; j++) {
        //      cout << landSurfType[j] << "  " << soilType[j] << "  " << areaFraction[j] << "  ";
        //    }
        //    cout << endl;
        if (sumArea != 100.0)
        {
            lakePercent = lakePercent * 100.0 / sumArea;
            glacierPercent = glacierPercent * 100.0 / sumArea;
            for (j = 0; j < maximumNumberLandClasses; j++)
            {
                areaFraction[j] = areaFraction[j] * 100.0 / sumArea;
            }
        }
        //    cout << lakePercent << "  " << glacierPercent << "  ";
        //    for (j=0; j<maximumNumberLandClasses; j++) {
        //      cout << areaFraction[j] << "  ";
        //    }
        //    cout << endl;
        Dew[landIndex].SetModelControlObj(ModelControlObject);
        Dew[landIndex].SetGeoIndex(geoIndex);
        Dew[landIndex].SetLandIndex(landIndex);
        Dew[landIndex].SetArea(elementArea);
        Dew[landIndex].SetElevation(elementElevation);
        Dew[landIndex].SetAspect(elementAspect);
        Dew[landIndex].SetFlowDirection(elementFlowDir);
        Dew[landIndex].SetSelectedKiWaTimeSeriesElements(SelectedKiWaTimeSeriesElementsStore);
        Dew[landIndex].SetSelectedHbvTimeSeriesElements(SelectedHbvTimeSeriesElementsStore);
        Dew[landIndex].SetSlopeLength(elementSlopeLength);
        Dew[landIndex].SetSlopeAngle(elementSlopeAngle);
        Dew[landIndex].SetGeneralPar(ParGeneralStore);
        Dew[landIndex].SetTotalReservoir(TotalReservoirStore);
        Dew[landIndex].AllocateInputValues(InputElementStore->GetNumberValues());
        // Allocate space for water balance time series for landscape elements
        for (j = 0; j < Dew[landIndex].GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++)
        {
            if (Dew[landIndex].GetLandIndex() == Dew[landIndex].GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j))
            {
                Dew[landIndex].AllocateKiWaWaterBalance(initialTimeSteps + numberTimeSteps);
            }
        }
        for (j = 0; j < Dew[landIndex].GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++)
        {
            if (Dew[landIndex].GetLandIndex() == Dew[landIndex].GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j))
            {
                Dew[landIndex].AllocateHbvWaterBalance(initialTimeSteps + numberTimeSteps);
            }
        }
        // Time series input data format - read meteorological station information
        if (inputFormat == 'T' || inputFormat == 't')
        {
            sumWeight = 0.0;
            precStationsWeightedElevation = 0.0;
            Dew[landIndex].AllocateMetSeries(ParGeneralStore->GetNUM_PREC_SERIES(), ParGeneralStore->GetNUM_TEMP_SERIES());
            for (j = 0; j < ParGeneralStore->GetNUM_PREC_SERIES(); j++)
            {
                finLandScape >> seriesNumber >> seriesWeight;
                Dew[landIndex].SetMetSeriesNumber(j, seriesNumber);
                Dew[landIndex].SetMetSeriesWeight(j, seriesWeight);
                precStationsWeightedElevation = precStationsWeightedElevation + MetStations->GetStationAltitude(seriesNumber) * seriesWeight;
                sumWeight = sumWeight + seriesWeight;
            }
            if (sumWeight != 1.0)
            {
                for (j = 0; j < ParGeneralStore->GetNUM_PREC_SERIES(); j++)
                {
                    seriesWeight = Dew[landIndex].GetMetSeriesWeight(j) / sumWeight;
                    Dew[landIndex].SetMetSeriesWeight(j, seriesWeight);
                    precStationsWeightedElevation = precStationsWeightedElevation / sumWeight;
                }
            }
	    Dew[landIndex].SetPrecStationsWeightedElevation(precStationsWeightedElevation);
            //      cout << " precStationsWeightedElevation = " << precStationsWeightedElevation << endl;
            sumWeight = 0.0;
            tempStationsWeightedElevation = 0.0;
            for (j = 0; j < ParGeneralStore->GetNUM_TEMP_SERIES(); j++)
            {
                finLandScape >> seriesNumber >> seriesWeight;
                Dew[landIndex].SetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES() + j, seriesNumber);
                Dew[landIndex].SetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES() + j, seriesWeight);
                tempStationsWeightedElevation = tempStationsWeightedElevation + MetStations->GetStationAltitude(seriesNumber) * seriesWeight;
                sumWeight = sumWeight + seriesWeight;
            }
            if (sumWeight != 1.0)
            {
                for (j = 0; j < ParGeneralStore->GetNUM_TEMP_SERIES(); j++)
                {
                    seriesWeight = Dew[landIndex].GetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES() + j) / sumWeight;
                    Dew[landIndex].SetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES() + j, seriesWeight);
                    tempStationsWeightedElevation = tempStationsWeightedElevation / sumWeight;
                }
            }
	    Dew[landIndex].SetTempStationsWeightedElevation(tempStationsWeightedElevation);
            //      cout << " tempStationsWeightedElevation = " << tempStationsWeightedElevation << endl;

            /*      for (j=0; j<ParGeneralStore->GetNUM_PREC_SERIES(); j++) {
            cout << "P " << Dew[landIndex].GetMetSeriesNumber(j) << "  " << Dew[landIndex].GetMetSeriesWeight(j) << endl;
            }
            for (j=0; j<ParGeneralStore->GetNUM_TEMP_SERIES(); j++) {
            cout << "T " << Dew[landIndex].GetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES()+j) <<
            "  " << Dew[landIndex].GetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES()+j) << endl;
            }*/
        }
        else
        {
            //      for (j=0; j<ParGeneralStore->GetNUM_PREC_SERIES(); j++) {
            //        finLandScape >> seriesNumber >> seriesWeight;
            //      }
            //      for (j=0; j<ParGeneralStore->GetNUM_TEMP_SERIES(); j++) {
            //        finLandScape >> seriesNumber >> seriesWeight;
            //      }
            finLandScape.ignore(1024, '\n');
        }
        if (lakePercent > 0.0)
        {
            ptrLake = new Lake;
            ptrLake->SetLandScapeElement(&Dew[landIndex]);
            ptrLake->SetAreaFraction(lakePercent);
	    //            ptrLake->SetLandSurfaceType(M_WATER);
	    //            ptrLake->SetSoilType(M_WATER_WATER);
	    /*if (evaporationModelling == 'P' || evaporationModelling == 'p')
            {
	      //                ptrLandAtmInt = new LandAtmosphereInterface();
	      //                ptrLake->SetLandAtmosphereInterfaceElement(ptrLandAtmInt);
	      //                ptrAlbedo = new Albedo;
	      //                ptrLake->SetAlbedoElement(ptrAlbedo);
	      //                ptrAtm = new Atmosphere;
	      //                ptrLake->SetAtmosphereElement(ptrAtm);
	      //                ptrPenMon = new PenmanMonteith();
	      //                ptrLake->SetPenmanMonteithElement(ptrPenMon);
	      }*/
            ptrLakeWaterBalance = new LakeWaterBalance;
            ptrLakeWaterBalance->SetGeneralPar(ParGeneralStore);
	    //            ptrLakeWaterBalance->SetLandSurfacePar(&ParLandSurfaceStore[M_WATER]);
	    //            ptrLakeWaterBalance->SetInputTimeSeries(InputTimeSeriesStore);
            ptrLakeWaterBalance->SetInputElement(InputElementStore);
            ptrLakeWaterBalance->SetLandScapeElement(&Dew[landIndex]);
            ptrLakeWaterBalance->SetLakeValues(ParGeneralStore->GetINITIAL_LAKE_TEMP(), ParGeneralStore->GetINITIAL_LAKE_LEVEL());
            ptrLake->SetLakeWaterBalance(ptrLakeWaterBalance);
            Dew[landIndex].SetLake(ptrLake);
        }
        if (glacierPercent > 0.0)
        {
            ptrGlacier = new Glacier(startModelTime, startSimulationTime, endSimulationTime);
            ptrGlacier->SetLandScapeElement(&Dew[landIndex]);
            ptrGlacier->SetAreaFraction(glacierPercent);
            ptrGlacier->SetLandSurfaceType(GLACIER);
            ptrGlacier->SetSoilType(GLACIER_BED);
            ptrGlacier->SetGlacierIceAreaFraction(glacierPercent);
	    /*if (evaporationModelling == 'P' || evaporationModelling == 'p')
            {
	      //                ptrLandAtmInt = new LandAtmosphereInterface;
	      //                ptrGlacier->SetLandAtmosphereInterfaceElement(ptrLandAtmInt);
	      //                ptrAlbedo = new Albedo;
	      //                ptrGlacier->SetAlbedoElement(ptrAlbedo);
	      //                ptrAtm = new Atmosphere;
	      //                ptrGlacier->SetAtmosphereElement(ptrAtm);
	      //                ptrPenMon = new PenmanMonteith;
	      //                ptrGlacier->SetPenmanMonteithElement(ptrPenMon);
	      }*/
            //      ptrSurface = new GlacierSurface;
            //      ptrSurface->SetGlacierSurfaceElevation(glacSurfElev);
            //      ptrSurface->SetGeneralPar(ParGeneralStore);
            //      ptrSurface->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]);
            //      ptrSurface->SetInputTimeSeries(InputTimeSeriesStore);
            //      ptrSurface->SetInputElement(InputElementStore);
            //      ptrSurface->SetLandScapeElement(&Dew[landIndex]);
            //      ptrGlacier->SetGlacierSurface(ptrSurface);
            ptrSnow = new Snow;
            ptrSnow->SetGeneralPar(ParGeneralStore);
            ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]);
            //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
            ptrSnow->SetInputElement(InputElementStore);
            ptrSnow->SetLandScapeElement(&Dew[landIndex]);
            ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
            ptrGlacier->SetSnow(ptrSnow);
            ptrIce = new GlacierIce;
            ptrIce->SetGlacierIceThickness(glacIceThick);
            ptrIce->SetGlacierIceVolume(glacIceThick * elementArea * glacierPercent / 100.0);
            // Glacier ice surface elevation
            ptrIce->SetGlacierSurfaceElevation(glacSurfElev);
            // Glacier bedrock surface elevation
            //    ptrIce->SetGlacierSurfaceElevation(glacSurfElev+glacIceThick);
            ptrIce->SetGeneralPar(ParGeneralStore);
            ptrIce->SetGlacierRetreatPar(&ParGlacRetStore[VALLEY]);
            ptrIce->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]);
            //      ptrIce->SetInputTimeSeries(InputTimeSeriesStore);
            ptrIce->SetInputElement(InputElementStore);
            ptrIce->SetLandScapeElement(&Dew[landIndex]);
            ptrGlacier->SetGlacierIce(ptrIce);
            ptrVegetation = new Vegetation;
            ptrVegetation->SetGeneralPar(ParGeneralStore);
            ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]);
            //      ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
            ptrVegetation->SetInputElement(InputElementStore);
            ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
            ptrGlacier->SetVegetation(ptrVegetation);
            if (modelStructure == HBV_MODEL)
            {
                ptrHbv = new HBV;
                ptrHbv->SetGeneralPar(ParGeneralStore);
                ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[GLACIER_BED]);
                ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]);
                //      ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
                ptrHbv->SetInputElement(InputElementStore);
                ptrHbv->SetLandScapeElement(&Dew[landIndex]);
                ptrHbv->SetInitialHbvValues();
                ptrHbv->SetSubSurfaceHbvStore(ParGeneralStore->GetINITIAL_SOIL_MOISTURE(),
                                              ParGeneralStore->GetINITIAL_UPPER_ZONE(), ParGeneralStore->GetINITIAL_LOWER_ZONE());
                ptrGlacier->SetHBV(ptrHbv);
            }
            else if (modelStructure == KWA_MODEL)
            {
                ptrKwa = new KinematicWave;
                ptrKwa->SetGeneralPar(ParGeneralStore);
                ptrKwa->SetKiWaPar(&ParKiWaStore[GLACIER_BED]);
                ptrKwa->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]);
                //      ptrKwa->SetInputTimeSeries(InputTimeSeriesStore);
                ptrKwa->SetInputElement(InputElementStore);
                ptrKwa->SetLandScapeElement(&Dew[landIndex]);
                ptrKwa->SetInitialKwaValues();
                ptrGlacier->SetKinematicWave(ptrKwa);
            }
            Dew[landIndex].SetGlacier(ptrGlacier);
        }
        if (modelStructure == HBV_MODEL)
        {
            if (areaFraction[0] > 0.0)
            {
                ptrHbvAquifer = new HbvAquifer;
                ptrHbvAquifer->SetLandScapeElement(&Dew[landIndex]);
                ptrHbvAquifer->SetAreaFraction(areaFraction[0]);
                ptrHbvAquifer->SetLandSurfaceType(landSurfType[0]);
                ptrHbvAquifer->SetSoilType(soilType[0]);
		/*if (evaporationModelling == 'P' || evaporationModelling == 'p')
                {
		  //                    ptrLandAtmInt = new LandAtmosphereInterface;
		  //                    ptrHbvAquifer->SetLandAtmosphereInterfaceElement(ptrLandAtmInt);
		  //                    ptrAlbedo = new Albedo;
		  //                    ptrHbvAquifer->SetAlbedoElement(ptrAlbedo);
		  //                    ptrAtm = new Atmosphere;
		  //                    ptrHbvAquifer->SetAtmosphereElement(ptrAtm);
		  //                    ptrPenMon = new PenmanMonteith;
		  //                    ptrHbvAquifer->SetPenmanMonteithElement(ptrPenMon);
		  }*/
                ptrVegetation = new Vegetation;
                ptrVegetation->SetGeneralPar(ParGeneralStore);
                ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]);
                //      ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
                ptrVegetation->SetInputElement(InputElementStore);
                ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
                ptrHbvAquifer->SetVegetation(ptrVegetation);
                ptrSnow = new Snow;
                ptrSnow->SetGeneralPar(ParGeneralStore);
                ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]);
                //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
                ptrSnow->SetInputElement(InputElementStore);
                ptrSnow->SetLandScapeElement(&Dew[landIndex]);
                ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
                ptrHbvAquifer->SetSnow(ptrSnow);
                ptrHbv = new HBV;
                ptrHbv->SetGeneralPar(ParGeneralStore);
                ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[soilType[0]]);
                ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]);
                //      ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
                ptrHbv->SetInputElement(InputElementStore);
                ptrHbv->SetLandScapeElement(&Dew[landIndex]);
                ptrHbv->SetInitialHbvValues();
                ptrHbv->SetSubSurfaceHbvStore(ParGeneralStore->GetINITIAL_SOIL_MOISTURE(),
                                              ParGeneralStore->GetINITIAL_UPPER_ZONE(), ParGeneralStore->GetINITIAL_LOWER_ZONE());
                ptrHbvAquifer->SetHBV(ptrHbv);
                Dew[landIndex].SetHbvAquifer(ptrHbvAquifer);
                j = 1;
                while (j < maximumNumberLandClasses && areaFraction[j] > 0.0)
                {
                    lastHbvAquifer = ptrHbvAquifer;
                    ptrHbvAquifer = new HbvAquifer;
                    ptrHbvAquifer->SetLandScapeElement(&Dew[landIndex]);
                    ptrHbvAquifer->SetAreaFraction(areaFraction[j]);
                    ptrHbvAquifer->SetLandSurfaceType(landSurfType[j]);
                    ptrHbvAquifer->SetSoilType(soilType[j]);
		    /*if (evaporationModelling == 'P' || evaporationModelling == 'p')
                    {
		      //                        ptrLandAtmInt = new LandAtmosphereInterface;
		      //                        ptrHbvAquifer->SetLandAtmosphereInterfaceElement(ptrLandAtmInt);
		      //                        ptrAlbedo = new Albedo;
		      //                        ptrHbvAquifer->SetAlbedoElement(ptrAlbedo);
		      //                        ptrAtm = new Atmosphere;
		      //                        ptrHbvAquifer->SetAtmosphereElement(ptrAtm);
		      //                        ptrPenMon = new PenmanMonteith;
		      //                        ptrHbvAquifer->SetPenmanMonteithElement(ptrPenMon);
		      }*/
                    ptrVegetation = new Vegetation;
                    ptrVegetation->SetGeneralPar(ParGeneralStore);
                    ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]);
                    //        ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
                    ptrVegetation->SetInputElement(InputElementStore);
                    ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
                    ptrHbvAquifer->SetVegetation(ptrVegetation);
                    ptrSnow = new Snow;
                    ptrSnow->SetGeneralPar(ParGeneralStore);
                    ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]);
                    //        ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
                    ptrSnow->SetInputElement(InputElementStore);
                    ptrSnow->SetLandScapeElement(&Dew[landIndex]);
                    ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
                    ptrHbvAquifer->SetSnow(ptrSnow);
                    ptrHbv = new HBV;
                    ptrHbv->SetGeneralPar(ParGeneralStore);
                    ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[soilType[j]]);
                    ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]);
                    //        ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
                    ptrHbv->SetInputElement(InputElementStore);
                    ptrHbv->SetLandScapeElement(&Dew[landIndex]);
                    ptrHbv->SetInitialHbvValues();
                    ptrHbv->SetSubSurfaceHbvStore(ParGeneralStore->GetINITIAL_SOIL_MOISTURE(),
                                                  ParGeneralStore->GetINITIAL_UPPER_ZONE(), ParGeneralStore->GetINITIAL_LOWER_ZONE());
                    ptrHbvAquifer->SetHBV(ptrHbv);
                    lastHbvAquifer->SetNextHbvAquifer(ptrHbvAquifer);
                    j++;
                }
            }
        }
        else if (modelStructure == KWA_MODEL)
        {
            if (areaFraction[0] > 0.0)
            {
                ptrKwaAquifer = new KwaAquifer;
                ptrKwaAquifer->SetLandScapeElement(&Dew[landIndex]);
                ptrKwaAquifer->SetAreaFraction(areaFraction[0]);
                ptrKwaAquifer->SetLandSurfaceType(landSurfType[0]);
                ptrKwaAquifer->SetSoilType(soilType[0]);
		/*if (evaporationModelling == 'P' || evaporationModelling == 'p')
                {
		  //                    ptrLandAtmInt = new LandAtmosphereInterface;
		  //                    ptrKwaAquifer->SetLandAtmosphereInterfaceElement(ptrLandAtmInt);
		  //                    ptrAlbedo = new Albedo;
		  //                    ptrKwaAquifer->SetAlbedoElement(ptrAlbedo);
		  //                    ptrAtm = new Atmosphere;
		  //                    ptrKwaAquifer->SetAtmosphereElement(ptrAtm);
		  //                    ptrPenMon = new PenmanMonteith;
		  //                    ptrKwaAquifer->SetPenmanMonteithElement(ptrPenMon);
		  }*/
                ptrVegetation = new Vegetation;
                ptrVegetation->SetGeneralPar(ParGeneralStore);
                ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]);
                //      ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
                ptrVegetation->SetInputElement(InputElementStore);
                ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
                ptrKwaAquifer->SetVegetation(ptrVegetation);
                ptrSnow = new Snow;
                ptrSnow->SetGeneralPar(ParGeneralStore);
                ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]);
                //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
                ptrSnow->SetInputElement(InputElementStore);
                ptrSnow->SetLandScapeElement(&Dew[landIndex]);
                ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
                ptrKwaAquifer->SetSnow(ptrSnow);
                ptrKwa = new KinematicWave;
                ptrKwa->SetGeneralPar(ParGeneralStore);
                ptrKwa->SetKiWaPar(&ParKiWaStore[soilType[0]]);
                ptrKwa->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]);
                //      ptrKwa->SetInputTimeSeries(InputTimeSeriesStore);
                ptrKwa->SetInputElement(InputElementStore);
                ptrKwa->SetLandScapeElement(&Dew[landIndex]);
                ptrKwa->SetInitialKwaValues();
                ptrKwaAquifer->SetKinematicWave(ptrKwa);
                Dew[landIndex].SetKwaAquifer(ptrKwaAquifer);
                j = 1;
                while (j < maximumNumberLandClasses && areaFraction[j] > 0.0)
                {
                    lastKwaAquifer = ptrKwaAquifer;
                    ptrKwaAquifer = new KwaAquifer;
                    ptrKwaAquifer->SetLandScapeElement(&Dew[landIndex]);
                    ptrKwaAquifer->SetAreaFraction(areaFraction[j]);
                    ptrKwaAquifer->SetLandSurfaceType(landSurfType[j]);
                    ptrKwaAquifer->SetSoilType(soilType[j]);
		    /*if (evaporationModelling == 'P' || evaporationModelling == 'p')
                    {
		      //                        ptrLandAtmInt = new LandAtmosphereInterface;
		      //                        ptrKwaAquifer->SetLandAtmosphereInterfaceElement(ptrLandAtmInt);
		      //                        ptrAlbedo = new Albedo;
		      //                        ptrKwaAquifer->SetAlbedoElement(ptrAlbedo);
		      //                        ptrAtm = new Atmosphere;
		      //                        ptrKwaAquifer->SetAtmosphereElement(ptrAtm);
		      //                        ptrPenMon = new PenmanMonteith;
		      //                        ptrKwaAquifer->SetPenmanMonteithElement(ptrPenMon);
		      }*/
                    ptrVegetation = new Vegetation;
                    ptrVegetation->SetGeneralPar(ParGeneralStore);
                    ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]);
                    //        ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
                    ptrVegetation->SetInputElement(InputElementStore);
                    ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
                    ptrKwaAquifer->SetVegetation(ptrVegetation);
                    ptrSnow = new Snow;
                    ptrSnow->SetGeneralPar(ParGeneralStore);
                    ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]);
                    //        ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
                    ptrSnow->SetInputElement(InputElementStore);
                    ptrSnow->SetLandScapeElement(&Dew[landIndex]);
                    ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
                    ptrKwaAquifer->SetSnow(ptrSnow);
                    ptrKwa = new KinematicWave;
                    ptrKwa->SetGeneralPar(ParGeneralStore);
                    ptrKwa->SetKiWaPar(&ParKiWaStore[soilType[j]]);
                    ptrKwa->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]);
                    //        ptrKwa->SetInputTimeSeries(InputTimeSeriesStore);
                    ptrKwa->SetInputElement(InputElementStore);
                    ptrKwa->SetLandScapeElement(&Dew[landIndex]);
                    ptrKwa->SetInitialKwaValues();
                    ptrKwaAquifer->SetKinematicWave(ptrKwa);
                    lastKwaAquifer->SetNextKwaAquifer(ptrKwaAquifer);
                    j++;
                }
            }
        }
    }

    ptrLake = 0;
    ptrGlacier = 0;
    ptrHbvAquifer = 0;
    lastHbvAquifer = 0;
    ptrKwaAquifer = 0;
    lastKwaAquifer = 0;
    ptrVegetation = 0;
    ptrSnow = 0;
    ptrIce = 0;
    ptrHbv = 0;
    ptrKwa = 0;
    //    ptrLandAtmInt = 0;
    //    ptrAlbedo = 0;
    //    ptrAtm = 0;
    //    ptrPenMon = 0;
    finLandScape.close();
    // End read file with landscape elements

    // Landscape element coordinates
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
		  Dew[k].SetXCoord((double)(xllCorner + (j+0.5)*cellSize));
		  Dew[k].SetYCoord((double)(yllCorner + (nRows-i-0.5)*cellSize));
		  //		  cout << " k " << k << " x " << (double)(xllCorner + (j+0.5)*cellSize);
		  //		  cout << " y " << (double)(yllCorner + (nRows-i-0.5)*cellSize) << "\n";
		  //		  cout << " k " << k << " x " << Dew[k].GetXCoord();
		  //		  cout << " y " << Dew[k].GetYCoord() << "\n";
		  k++;
		}
	    }
	}
    }
    // End landscape element coordinates

    // Read name of file with observed streamflow data 
    /*  cout << " File with observed streamflow data: ";
    cin >> fileObsStreamflow;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileObsStreamflow;
    fileControl.ignore(1024, '\n');
    // End read name of file with observed streamflow data

    // Read sub-catchment information file and generate sub-catchment objects
    /*  cout << " File with sub-catchment hierarchy: ";
    cin >> fileNameSubCatchment;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileNameSubCatchment;
    fileControl.ignore(1024, '\n');
    ifstream fileWCo(fileNameSubCatchment);
    if (!fileWCo.is_open())
    {
        cout << endl << " Error opening file " << fileNameSubCatchment << endl << endl;
        exit(1);
    }
    // Sub-catchment elements
    fileWCo.ignore(100, ':');
    fileWCo >> numWatc;
    cout << "\n # Number of sub-catchments " << numWatc << endl;
    SubCatchment *CatchmentElement = new SubCatchment[numWatc];
    // Outlet lake routing sub-catchment
    int * outletRoutingSubCatchments = new int[numWatc];
    for (i = 0; i < numWatc; i++)
    {
        fileWCo >> j >> ch >> subCatchmentId >> lakeNumber >> maxBas >> correction;
        if (j != i)
        {
            cout << endl << "Error reading file " << fileNameSubCatchment << "\t" << i << "\t" << j << endl;
            exit(1);
        }
	//	if (maxBas != 1.0) maxBasTransform = true;
        fileWCo.ignore(1024, '\n');
        CatchmentElement[i].SetSubCatchmentIndex(i);
        CatchmentElement[i].SetIdentifier(subCatchmentId);
        CatchmentElement[i].SetLakeNumber(lakeNumber);
	CatchmentElement[i].SetMaxBas(maxBas);
        CatchmentElement[i].SetCorrection(correction);
        CatchmentElement[i].AllocateAccumulatedDischarge(initialTimeSteps + numberTimeSteps);
        CatchmentElement[i].AllocateAccumulatedInFlow(initialTimeSteps + numberTimeSteps);
        CatchmentElement[i].AllocateAccumulatedWaterBalance(initialTimeSteps + numberTimeSteps);
        CatchmentElement[i].AllocateWaterBalance(initialTimeSteps + numberTimeSteps);
        CatchmentElement[i].ObsDataInput(fileObsStreamflow, startSimulationTime, endSimulationTime, numberTimeSteps,
                                         ParGeneralStore->GetSECONDS_TIMESTEP());
        CatchmentElement[i].SetSelectedSubCatchmentTimeSeriesElements(SelectedSubCatchmentTimeSeriesElementsStore);
	outletRoutingSubCatchments[i] = lakeNumber;
        cout << " Sub-catchment index " << i << "  " << "Sub-catchment identifier " << CatchmentElement[i].GetIdentifier() << "  " << "Sub-catchment correction " << CatchmentElement[i].GetCorrection() << endl;
    }
    // Watercourse outlets
    fileWCo.ignore(100, ':');
    fileWCo >> numWatcOut;
    cout << "\n # Number of watercourse outlets " << numWatcOut << endl;
    SubCatchment ** Outlet = new SubCatchment *[numWatcOut];
    for (i = 0; i < numWatcOut; i++)
    {
        fileWCo >> j;
        Outlet[i] = &CatchmentElement[j];
        cout << " Outlet no. " << i << "\t" << " Sub-catchment no. " << j << "\t" << endl;
        fileWCo.ignore(1024, '\n');
    }
    // Hierarchy of sub-catchments
    fileWCo.getline(buffer, 1024);
    cout << "\n " << buffer << endl;
    while (fileWCo >> i)
    {
        fileWCo >> numWatcUp;
        CatchmentElement[i].SetNumUpStream(numWatcUp);
        fileWCo.ignore(100, ':');
        cout << " Downstream, sub-catchment no.  " << i << "    Identifier  " << CatchmentElement[i].GetIdentifier() << endl;
        cout << " No. of upstream sub-catchments " << numWatcUp << endl;
        k = 0;
        while (fileWCo.peek() != '\n')
        {
            fileWCo >> j;
            cout << "\t" << "Upstream, sub-catchment no. " << j;
            CatchmentElement[i].SetUpStream(k, &CatchmentElement[j]);
            cout << "\t" << "UpStream[" << k << "]" << "    Identifier  " << CatchmentElement[i].GetUpStream(k)->GetIdentifier() << endl;
            while (fileWCo.peek() == ' ')
            {
                fileWCo.ignore(1, ' ');
            }
            k++;
        }
        fileWCo.ignore(1024, '\n');
        if (numWatcUp != k)
        {
            cout << endl << " Error in number of upstream pointers for sub-catchment no. " << i << endl << endl;
            exit(1);
        }
    }
    fileWCo.close();
    // End read file with sub-catchment information

     // Find lake outlets 
    for (i = 0; i < numWatcOut; i++)
    {
        TraverseOutletRoutingSubCatchment(Outlet[i], numWatc, outletRoutingSubCatchments);
    }

    // Read information about sub-catchment elements and landscape elements
    ReadSubCatchmentIdentifier(Dew, CatchmentElement, numWatc, flowHierarchy, fileControl, fout);

    // Read hierarchy of landscape elements
    if (flowHierarchy)
    {
        ReadLandscapeHierarchy(Dew, fileControl, fout);
    }
    else
    {
        fileControl.getline(buffer, 1024);
        //    fileControl.ignore(1024,'\n');
        //    cout << "\n " << buffer << endl;
    }

    // Write information about sub-catchment outlets and landscape elements to file test_waterland.txt
    //  WriteSubCatchmentIdentifier(CatchmentElement, numWatc, fout);

    // Lists of glacier elements
    GlacierElements * DewGlacierElements = new GlacierElements(Outlet, numWatcOut);
    DewGlacierElements->BuildGlacierLists();

    // Precipitation and temperature correction for catchments
    int * correctionCatchments = new int[maximumCorrectionCatchments];
    double * correctionPrecipitation = new double[maximumCorrectionCatchments];
    double * correctionTemperature = new double[maximumCorrectionCatchments];
    /*  cout << " File with precipitation and temperature correction for catchments: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    //  cout << fileName << endl;
    ifstream finCorrection(fileName);  // Open for reading
    if (!finCorrection.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        exit(1);
    }
    finCorrection.getline(buffer, 1024);
    i = 0;
    fout << "Precipitation and temperature correction for catchments: \n";
    while (finCorrection >> correctionCatchments[i])
    {
        finCorrection >> correctionPrecipitation[i] >> correctionTemperature[i];
        fout << i << "  " << correctionCatchments[i] << "  " << correctionPrecipitation[i] << "  " << correctionTemperature[i] << endl;
        i++;
    }
    numberCorrectionCatchments = i;
    SetElementCorrection(Dew, numLand, numberCorrectionCatchments, correctionCatchments, correctionPrecipitation, correctionTemperature, fout);
    for (i = 0; i < numWatcOut; i++)
    {
        TraverseCorrectionSubCatchment(Outlet[i], numberCorrectionCatchments, correctionCatchments, correctionPrecipitation, correctionTemperature);
    }
    finCorrection.close();
    delete[] outletRoutingSubCatchments;
    delete[] correctionCatchments;
    delete[] correctionPrecipitation;
    delete[] correctionTemperature;
    // End precipitation and temperature correction for catchments

    if (inputFormat == 'T' || inputFormat == 't')
    {
        /*    cout << " File with meteorological input data: ";
        cin >> fileName;
        cout << endl;*/
        fileControl.ignore(100, ':');
        fileControl >> fileNameInput;
        fileControl.ignore(1024, '\n');
        //    cout << fileNameInput << endl;
    }
    else
    {
        fileControl.getline(buffer, 1024);
        //    fileControl.ignore(1024,'\n');
        //    cout << "\n " << buffer << endl;
    }

    // Source-To-Sink routing travel times for all landscape elements
    if (routingType == 'S' || routingType == 's')
    {
        // Water velocities for landscape elements
        double * waterVelocity = new double[maximumNumberWaterVelocities];
        /*    cout << " File with water velocities: ";
        cin >> fileName;
        cout << endl;*/
        fileControl.ignore(100, ':');
        fileControl >> fileName;
        fileControl.ignore(1024, '\n');
        ifstream finVelocity(fileName);  // Open for reading
        if (!finVelocity.is_open())
        {
            cout << endl << " Error opening file " << fileName << endl << endl;
            exit(1);
        }
        finVelocity.getline(buffer, 1024);
        i = 0;
	cout << endl << " # Water velocities for landscape elements " << endl << endl;
        while (finVelocity.getline(buffer, 1024))
        {
            //    cout << buffer << endl;
            sscanf(buffer, "%d %lf ", &velocityIndex, &waterVelocity[i]);
            if (velocityIndex != i)
            {
                cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << velocityIndex << endl;
                exit(1);
            }
            cout << velocityIndex << " \t" << waterVelocity[velocityIndex] << endl;
            i++;
        }
        finVelocity.close();
        travelTime = 0.0;
        for (i = 0; i < numWatcOut; i++)
        {
            TraverseWaterCourseTravelTime(Outlet[i], travelTime, waterVelocity, fout);
        }
        delete[] waterVelocity;
    }
    else
    {
        fileControl.getline(buffer, 1024);
        //    fileControl.ignore(1024,'\n');
        //    cout << "\n " << buffer << endl;
    }
    double * sourceToSinkDischarge = new double[initialTimeSteps + numberTimeSteps];
    for (i = 0; i < initialTimeSteps + numberTimeSteps; i++)
    {
        sourceToSinkDischarge[i] = 0.0;
    }
    // End Source-To-Sink routing travel times for all landscape elements

    // Time series format input data
    if (inputFormat == 'T' || inputFormat == 't')
    {
        // Object for storing input data for time series
        InputTimeSeries * InputTimeSeriesStore = new
        InputTimeSeries(initialTimeSteps + numberTimeSteps, MetStations->GetNumPrecStations() + MetStations->GetNumTempStations(),
                        startModelTime, endSimulationTime, ParGeneralStore->GetSECONDS_TIMESTEP());
        InputTimeSeriesStore->SetGeneralPar(ParGeneralStore);
        ifstream finInput(fileNameInput);  // Open for reading
        if (!finInput.is_open())
        {
            cout << endl << " Error opening file " << fileNameInput << endl << endl;
            exit(1);
        }
        InputTimeSeriesStore->SetInput(finInput);
        finInput.close();
	// Write input data to test file
	//	InputTimeSeriesStore->WriteInput();

        // Water balance for all elements and time steps for spin-up period and simulation period
        timeStep = 0;
        for (datetime = startModelTime; datetime <= endSimulationTime; datetime += ParGeneralStore->GetSECONDS_TIMESTEP())
        {
            //    cout << "  " << datetime.getYear() << "  " << datetime.getMonth() << "  " << datetime.getDay() << "  "
            //         << datetime.getHour() << "  " << datetime.getMinute() << "  " << datetime.getSecond() << endl;
            //      if (timeStep==0) WriteAsciiGridAnnualGlacierValues(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);

            // ** Algorithm to be performed in case: no input to landscape element from upstream elements
            if (!flowHierarchy)
            {
                inputDataFound = false;
                inputDataError = false;
                firstTotal = true;
                WaterBalanceTimeSeries(Dew, ParGeneralStore, MetStations, initialTimeSteps, numberTimeSteps,
                                       InputTimeSeriesStore, InputElementStore,
                                       numLand, timeStep, datetime, &inputDataFound, &inputDataError, &firstTotal);
                // Traverse sub-catchments and landscape elements
                if (!inputDataError)
                {
                    for (i = 0; i < numWatcOut; i++)
                    {
                        TraverseSubCatchment(Outlet[i], ParGeneralStore, MetStations, InputTimeSeriesStore, InputElementStore, timeStep, datetime, initialTimeSteps, numberTimeSteps,
                                             inputFormat, flowHierarchy, forceDirect, &inputDataFound, &inputDataError, &firstTotal, fout);
                    }
                }
                else
                {
                    for (i = 0; i < numWatcOut; i++)
                    {
                        TraverseMissingDataSubCatchment(Outlet[i], timeStep, fout);
                    }
                }
            }
            // ** End algorithm to be performed in case: no input to landscape element from upstream elements

            // ** Algorithm to be performed in case: input to landscape element from upstream elements
            else
            {
                inputDataFound = false;
                inputDataError = false;
                firstTotal = true;
                for (i = 0; i < numWatcOut; i++)
                {
                    TraverseSubCatchment(Outlet[i], ParGeneralStore, MetStations, InputTimeSeriesStore, InputElementStore, timeStep, datetime, initialTimeSteps, numberTimeSteps,
                                         inputFormat, flowHierarchy, forceDirect, &inputDataFound, &inputDataError, &firstTotal, fout);
                    if (inputDataError)
                    {
                        TraverseMissingDataSubCatchment(Outlet[i], timeStep, fout);
                    }
                }
            }
            // ** End algorithm to be performed in case: input to landscape element from upstream elements

            SnowGlacierIceReDistribution(Outlet, Dew, ParGeneralStore, initialTimeSteps, numberTimeSteps, numLand, numWatcOut, timeStep, datetime, glacierModelling, DewGlacierElements,
					 nRows, nCols, noData, xllCorner, yllCorner, cellSize, modelCalibration, fout); 
	    //            if (timeStep == (int)(initialTimeSteps) || timeStep == (int)(initialTimeSteps + numberTimeSteps / 2.0) || timeStep == initialTimeSteps + numberTimeSteps - 1)
	    //	    {
	    //	      if (!modelCalibration) 
	    //	      {
	    //		WriteAsciiGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
	    //	      }
	    //            }
            if (routingType == 'S' || routingType == 's')
            {
                RouteSourceToSinkDischarge(Dew, numLand, sourceToSinkDischarge, initialTimeSteps, numberTimeSteps, timeStep, ParGeneralStore->GetSECONDS_TIMESTEP());
            }
            timeStep++;
        }
        if (timeStep != initialTimeSteps + numberTimeSteps)
        {
            cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
            exit(1);
        }
        delete InputTimeSeriesStore;
    }

    // Grid file format input data
    else
    {
        // Dummy variable InputTimeSeriesStore
        InputTimeSeries * InputTimeSeriesStore = new
        InputTimeSeries(initialTimeSteps + numberTimeSteps, MetStations->GetNumPrecStations() + MetStations->GetNumTempStations(),
                        startModelTime, endSimulationTime, ParGeneralStore->GetSECONDS_TIMESTEP());

        // Path to grid files with meteorological input data
        metPath = getenv("METDATA");
        if (!metPath)
        {
            cout << " Environment variable METDATA not defined " << endl << endl;
            exit(1);
        }
        unsigned short int * precip10 = new unsigned short int [nRows*nCols];
        unsigned short int * temp10K = new unsigned short int [nRows*nCols];
        //    unsigned short int * precip10 = new unsigned short int[numRows1km * numCols1km];
        //    unsigned short int * temp10K = new unsigned short int[numRows1km * numCols1km];
	/*        metMask = getenv("METMASK");
        if (!metMask)
        {
            cout << " Environment variable METMASK not defined " << endl << endl;
            exit(1);
        }
        unsigned short int * precip10 = new unsigned short int [numLand];
        unsigned short int * temp10K = new unsigned short int [numLand];*/

        // Water balance for all elements and time steps for spin-up period and simulation period
        timeStep = 0;
        for (datetime = startModelTime; datetime <= endSimulationTime; datetime += ParGeneralStore->GetSECONDS_TIMESTEP())
        {
            //    cout << "  " << datetime.getYear() << "  " << datetime.getMonth() << "  " << datetime.getDay() << "  "
            //           << datetime.getHour() << "  " << datetime.getMinute() << "  " << datetime.getSecond() << endl;
            inputDataFound = false;
            inputDataError = false;
            firstTotal = true;
	    WaterBalanceGrid(Dew, ParGeneralStore, InputElementStore, initialTimeSteps, numberTimeSteps,
                             numLand, timeStep, nRows, nCols,
                             datetime, metPath, precip10, temp10K, flowHierarchy, &inputDataFound, &firstTotal);
	    /*            WaterBalanceGridVariable(Dew, ParGeneralStore, InputElementStore, initialTimeSteps, numberTimeSteps,
                             numLand, timeStep, nRows, nCols,
                             datetime, metPath, precip10, temp10K, flowHierarchy, &inputDataFound, &firstTotal);*/
	    /*            WaterBalanceElements(Dew, ParGeneralStore, InputElementStore, initialTimeSteps, numberTimeSteps,
                             numLand, timeStep, nRows, nCols,
                             datetime, metMask, precip10, temp10K, flowHierarchy, &inputDataFound, &firstTotal);*/
	    // Traverse sub-catchments and landscape elements
            if (inputDataFound)
            {
                for (i = 0; i < numWatcOut; i++)
                {
                    TraverseSubCatchment(Outlet[i], ParGeneralStore, MetStations, InputTimeSeriesStore, InputElementStore, timeStep, datetime, initialTimeSteps, numberTimeSteps,
                                         inputFormat, flowHierarchy, forceDirect, &inputDataFound, &inputDataError, &firstTotal, fout);
                }
            }
            else
            {
                for (i = 0; i < numWatcOut; i++)
                {
                    TraverseMissingDataSubCatchment(Outlet[i], timeStep, fout);
                }
            }
            SnowGlacierIceReDistribution(Outlet, Dew, ParGeneralStore, initialTimeSteps, numberTimeSteps, numLand, numWatcOut, timeStep, datetime, glacierModelling, DewGlacierElements,
                                         nRows, nCols, noData, xllCorner, yllCorner, cellSize, modelCalibration, fout);
            // Write state variables for all landscape elements
	    if (timeStep == (int)(initialTimeSteps) || timeStep == (int)(initialTimeSteps + numberTimeSteps / 2.0) || timeStep == initialTimeSteps + numberTimeSteps - 1)
	      {
	    //	      WriteReducedBinaryGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
	    //	      WriteBinaryGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
		//		if (!modelCalibration) 
		//		  {
		//		    WriteAsciiGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
		//		  }
	      }
            if (routingType == 'S' || routingType == 's')
            {
                RouteSourceToSinkDischarge(Dew, numLand, sourceToSinkDischarge, initialTimeSteps, numberTimeSteps, timeStep, ParGeneralStore->GetSECONDS_TIMESTEP());
            }
            timeStep++;
        }
        if (timeStep != initialTimeSteps + numberTimeSteps)
        {
            cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
            exit(1);
        }
        delete[] precip10;
        delete[] temp10K;
        delete InputTimeSeriesStore;
    }

    // Muskingum-Cunge routing
    if (routingType == 'M' || routingType == 'm')
    {

        // Read Manning roughness coefficients
        /*    cout << " File with Manning roughness coefficients: ";
        cin >> fileName;
        cout << endl;*/
        fileControl.ignore(100, ':');
        fileControl >> fileName;
        fileControl.ignore(1024, '\n');
        ifstream finManning(fileName);  // Open for reading
        if (!finManning.is_open())
        {
            cout << endl << " Error opening file " << fileName << endl << endl;
            exit(1);
        }
        finManning.getline(buffer, 1024);
        i = 0;
	cout << endl << " # Manning roughness coefficients " << endl << endl;
        while (finManning.getline(buffer, 1024))
        {
            //    cout << buffer << endl;
            sscanf(buffer, "%s %d %lf ", manningVegetationType, &manningIndex, &manningRoughness[i]);
            if (manningIndex != i)
            {
                cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << manningIndex << endl;
                exit(1);
            }
            cout << manningVegetationType << " \t" << manningIndex << " \t" << manningRoughness[manningIndex] << endl;
            i++;
        }
        finManning.close();

        // Read routing parameters of river and lake elements
        /*    cout << " File with routing parameters: ";
        cin >> fileNameRouting;
        cout << endl;*/
        fileControl.ignore(100, ':');
        fileControl >> fileNameRouting;
        fileControl.ignore(1024, '\n');
        ifstream fileRout(fileNameRouting);
        if (!fileRout.is_open())
        {
            cout << endl << "Error opening file " << fileNameRouting << endl << endl;
            exit(1);
        }
        // Sub-catchment elements with routing parameters
        fileRout.ignore(100, ':');
        fileRout >> numWatcRoute;
        cout << "\n # Number of sub-catchments " << numWatcRoute << endl;
        for (i = 0; i < numWatcRoute; i++)
        {
            fileRout >> j >> ch >> subCatchmentId >> manningIndex >> riverSlope >> width;
            if (j != i)
            {
                cout << endl << "Error reading file " << fileNameRouting << "\t" << i << "\t" << j << endl;
                exit(1);
            }
            /*    if (i != CatchmentElement[i].GetSubCatchmentIndex()) {
            cout << endl << "Error reading file " << fileNameRouting << "\t" << i << "\t" << CatchmentElement[i].GetSubCatchmentIndex() << endl;
            exit (1);
            }
            if (subCatchmentId != CatchmentElement[i].GetIdentifier()) {
            cout << endl << "Error reading file " << fileNameRouting << "\t" << i << "\t" << CatchmentElement[i].GetIdentifier() << endl;
            exit (1);
            }*/

            // Find the sub-catchment identification of river cell
            for (j = 0; j < numWatc; j++)
            {
                if (subCatchmentId == CatchmentElement[j].GetIdentifier())
                {
                    currentCat = j;
                    break;
                }
            }
            if (j == numWatc)
            {
                cout << endl << "Error reading file " << fileNameRouting << "\t Sub-catchment " << subCatchmentId << "  not found in CatchmentElement " << endl;
                exit(1);
            }

            fileRout.ignore(1024, '\n');
            CatchmentElement[currentCat].SetManning(manningRoughness[manningIndex]);
            CatchmentElement[currentCat].SetRiverSlope(riverSlope);
            CatchmentElement[currentCat].SetWidth(width);
        }
        fileRout.close();

        // Routing of water through network of river and lake elements
        timeStep = 0;
        for (datetime = startModelTime; datetime < startSimulationTime; datetime += ParGeneralStore->GetSECONDS_TIMESTEP())
        {
            for (i = 0; i < numWatcOut; i++)
            {
                TraverseRouteWaterCourse(Outlet[i], timeStep, ParGeneralStore->GetSECONDS_TIMESTEP(), fout);
            }
            timeStep++;
        }
        if (timeStep != initialTimeSteps)
        {
            cout << " timeStep != initialTimeSteps " << timeStep << "  " << initialTimeSteps << endl << endl;
            exit(1);
        }
        for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += ParGeneralStore->GetSECONDS_TIMESTEP())
        {
            for (i = 0; i < numWatcOut; i++)
            {
                TraverseRouteWaterCourse(Outlet[i], timeStep, ParGeneralStore->GetSECONDS_TIMESTEP(), fout);
            }
            timeStep++;
        }
        if (timeStep != initialTimeSteps + numberTimeSteps)
        {
            cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
            exit(1);
        }
    }
    else
    {
        fileControl.getline(buffer, 1024);
        //    fileControl.ignore(1024,'\n');
        //    cout << "\n " << buffer << endl;
        fileControl.getline(buffer, 1024);
        //    fileControl.ignore(1024,'\n');
        //    cout << "\n " << buffer << endl;
    }
    // End Muskingum-Cunge routing

    //Transform discharge from sub-catchments with triangular weight function
    if (maxBasTransform) {
      for (i=0; i<numWatcOut; i++) {
	TraverseMaxBas(Outlet[i], initialTimeSteps, numberTimeSteps, fout);
      }
    }
 
    // Write water balance grid
    WriteAsciiGridWaterBalance(Dew, startSimulationTime, endSimulationTime, numLand, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
 
    // Write discharge from all watercourse/sub-catchment elements in sub-catchment hierarchy to output files
    WriteSubCatchmentDischarge(CatchmentElement, numWatc, startSimulationTime, endSimulationTime, initialTimeSteps,
                               numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP(), modelCalibration, fout);

    // Write water balance from all sub-catchment elements in sub-catchment hierarchy to output files
    WriteSubCatchmentWaterBalance(CatchmentElement, numWatc, startSimulationTime, endSimulationTime,
                                  initialTimeSteps, numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP());

    // Write discharge from all watercourse outlets to output files
    for (i = 0; i < numWatcOut; i++)
        WriteOutletDischarge(Outlet[i], startSimulationTime, endSimulationTime, initialTimeSteps,
                             numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP(), false, fout);

    // Write water balance from all watercourse outlets to output files
    for (i = 0; i < numWatcOut; i++)
        WriteOutletWaterBalance(Outlet[i], startSimulationTime, endSimulationTime, initialTimeSteps,
                                numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP(), flowHierarchy, forceDirect);

    // Write state variable time series for landscape elements selected for output
    WriteDistributedElementTimeSeries(Dew, numLand, startSimulationTime, endSimulationTime,
                                      initialTimeSteps, numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP());

    // Write total reservoir storage for all time steps
    /*    WriteTotalReservoirStorage(TotalReservoirStore, startSimulationTime, endSimulationTime, initialTimeSteps,
	  numberTimeSteps, modelCalibration, ParGeneralStore->GetSECONDS_TIMESTEP());*/

    // Write Source-To-Sink discharge for catchment outlet
    if (routingType == 'S' || routingType == 's')
    {
        WriteSourceToSinkDischarge(Outlet[0], sourceToSinkDischarge, startSimulationTime, endSimulationTime, initialTimeSteps, numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP());
        delete[] sourceToSinkDischarge;
    }


    /*  k=0;
    fout << "\nPrecipitation correction grid:\n";
    for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
    if (k<numLand) {
    if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
    fout.width(10);
    fout << Dew[k].GetPrecipitationCorrection() << "  ";
    k++;
    }
    else {
    fout.width(10); fout << noData << "  ";
    }
    }
    else {
    fout.width(10); fout << noData << "  ";
    }
    }
    fout << endl;
    }
    fout << endl;
    k=0;
    fout << "\nTemperature correction grid:\n";
    for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
    if (k<numLand) {
    if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
    fout.width(10);
    fout << Dew[k].GetTemperatureCorrection() << "  ";
    k++;
    }
    else {
    fout.width(10); fout << noData << "  ";
    }
    }
    else {
    fout.width(10); fout << noData << "  ";
    }
    }
    fout << endl;
    }*/


    fout << endl;
    fout.close();
    fileControl.close();

    delete[] CatchmentElement;
    delete[] Outlet;
    delete[] Dew;
    delete ParGeneralStore;
    delete[] ParLandSurfaceStore;
    delete[] ParSubSurfaceHbvStore;
    delete InputElementStore;

    return 0;
}


void ReadSubCatchmentIdentifier(DistributedElement * const Dew, SubCatchment * const CatchmentElement, int numWatc, bool flowHierarchy, ifstream &fileControl, ofstream &fout)
{
    char fileName[80];
    char ch;
    int i, j;
    int subCatchmentId, numLandScape;
    int landIndex, geoIndex;
    DistributedElement * thisElement;

    // Read information about sub-catchment elements and landscape elements
    /*  cout << "\n File with information about sub-catchment elements and landscape elements: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finSubCatchment(fileName);  // Open for reading
    if (!finSubCatchment.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        exit(1);
    }
    // Connect landscape elements to sub-catchment elements
    for (i = 0; i < numWatc; i++)
    {
        finSubCatchment >> ch >> subCatchmentId >> ch >> numLandScape;
        finSubCatchment.ignore(1024, '\n');
        CatchmentElement[i].SetNumLandScape(numLandScape);
        if (subCatchmentId != CatchmentElement[i].GetIdentifier())
        {
            cout << endl << " Error reading file " << fileName << " for sub-catchment " << i << "\t"
                 << subCatchmentId << endl << endl;
            exit(1);
        }
        if (numLandScape > 0)
        {
            finSubCatchment >> landIndex >> geoIndex;
            thisElement = &Dew[landIndex];
            thisElement->SetLakeNumber(CatchmentElement[i].GetLakeNumber());
            CatchmentElement[i].SetLandScapeElement(thisElement);
            //    cout << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
            if (flowHierarchy)
            {
                thisElement->SetSubCatchmentElement(&CatchmentElement[i]);
            }
            for (j = 1; j < numLandScape; j++)
            {
                finSubCatchment >> landIndex >> geoIndex;
                thisElement->SetNextElement(&Dew[landIndex]);
                thisElement = &Dew[landIndex];
                //      cout << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
            }
        }
    }
}


void ReadLandscapeHierarchy(DistributedElement * const Dew, ifstream &fileControl, ofstream &fout)
{
    char fileName[80];
    int i, j, k;
    int numUp;

    // Read file with pointers between landscape elements
    /*  cout << " File with landscape element hierarchy: ";
    cin >> fileName;
    cout << endl;*/
    fileControl.ignore(100, ':');
    fileControl >> fileName;
    fileControl.ignore(1024, '\n');
    ifstream finLandUpFlow(fileName);  // Open for reading
    if (!finLandUpFlow.is_open())
    {
        cout << endl << " Error opening file " << fileName << endl << endl;
        exit(1);
    }
    while (finLandUpFlow >> i)
    {
        finLandUpFlow >> numUp;
        Dew[i].SetNumUpLand(numUp);
        finLandUpFlow.ignore(100, ':');
        k = 0;
        while (finLandUpFlow.peek() != '\n')
        {
            finLandUpFlow >> j;
            Dew[i].SetUpLandFlow(k, &Dew[j]);
            k++;
            while (finLandUpFlow.peek() == ' ')
            {
                finLandUpFlow.ignore(1, ' ');
            }
        }
        finLandUpFlow.ignore(1024, '\n');
        if (numUp != k)
        {
            cout << endl << " Error in number of upland pointers for landscape element no. "
                 << i << endl << endl;
            exit(1);
        }
    }
    finLandUpFlow.close();
}


void SnowGlacierIceReDistribution(SubCatchment ** const Outlet, DistributedElement * const Dew, ParametersGeneral * ParGeneralStore, int initialTimeSteps, int numberTimeSteps,
                                  int numLand, int numWatcOut, int timeStep, DateTime datetime, char glacierModelling, GlacierElements * DewGlacierElements,
                                  int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, bool modelCalibration, ofstream &fout)
{
    int i;
    // Snow store at glaciers is converted to ice at day no. DAY_ANNUAL_GLACIER
    if (ParGeneralStore->GetDAY_ANNUAL_GLACIER() > 0 &&
            dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay()) ==
            ParGeneralStore->GetDAY_ANNUAL_GLACIER() + leapYear(datetime.getYear()))
    {
        //    cout << "    Glacier snow " << datetime.getYear() << " " << datetime.getMonth() << " " << datetime.getDay() << " " << endl;
        for (i = 0; i < numLand; i++)
        {
            //    cout << "\n * Glacier SetSnowStore " << Dew[i].GetGeoIndex() << endl;
	    // Convert snow to glacier ice when changing elevation of glacier surface based on mass balance
            if ((timeStep > initialTimeSteps) && (glacierModelling == 'M' || glacierModelling == 'm')) 
            {
                Dew[i].SnowToGlacierIce();
            }
            Dew[i].RemoveSnowOnGlacierIce();
            Dew[i].SetAnnualMassBalance(missingData);
        }
    }
    // Snow store is removed at day no. DAY_SNOW_ZERO
    if (ParGeneralStore->GetDAY_SNOW_ZERO() > 0 &&
            dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay()) ==
            ParGeneralStore->GetDAY_SNOW_ZERO() + leapYear(datetime.getYear()))
    {
        //    cout << "            Snow " << datetime.getYear() << " " << datetime.getMonth() << " " << datetime.getDay() << " " << endl;
        for (i = 0; i < numLand; i++)
        {
            Dew[i].SetSnowStore(0.0);
        }
    }
    // Glacier surface elevation is redistributed at day no. DAY_ANNUAL_GLACIER
    if (ParGeneralStore->GetDAY_ANNUAL_GLACIER() > 0 &&
            dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay()) ==
            ParGeneralStore->GetDAY_ANNUAL_GLACIER() + leapYear(datetime.getYear()))
    {
        DewGlacierElements->SetThisYearAnnualGlacierValues(datetime);
        if (glacierModelling == 'E' || glacierModelling == 'e' || glacierModelling == 'M' || glacierModelling == 'm')
        {
            if (timeStep > initialTimeSteps)
            {
	        // Redistribution of glacier surface based on Huss
	        if (glacierModelling == 'E' || glacierModelling == 'e')
	        {
		    DewGlacierElements->GlacierSurfaceElevationReDistribution();
		}
	        // Change elevation of glacier surface based on mass balance
	        else if (glacierModelling == 'M' || glacierModelling == 'm')
	        {
		    DewGlacierElements->ChangeGlacierSurfaceElevation();
		}
                DewGlacierElements->RemoveElementsGlacierLists();
                DewGlacierElements->SortGlacierLists();
            }
	    //	    if (!modelCalibration) 
	    //            {
	    WriteAsciiGridAnnualGlacierValues(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
	    //	    }
            cout << " Glacier surface " << datetime.getYear() << " " << datetime.getMonth() << " " << datetime.getDay() << " " << endl;
        }
    }
}


void WaterBalanceTimeSeries(DistributedElement * const Dew, ParametersGeneral * const ParGeneralStore,
                            MeteorologicalStations * const MetStations, int initialTimeSteps, int numberTimeSteps,
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore,
                            int numLand, int timeStep, DateTime datetime, bool * inputDataFound, bool * inputDataError, bool * firstTotal)
{
    int i;
    for (i = 0; i < numLand; i++)
    {
        SetTimeSeriesInput(&Dew[i], ParGeneralStore, MetStations, initialTimeSteps, numberTimeSteps,
                           InputTimeSeriesStore, InputElementStore, timeStep, inputDataFound, inputDataError);
        if (*inputDataFound)
        {
            Dew[i].WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, 0.0, 0.0);
	    if (timeStep == initialTimeSteps-1) Dew[i].SetInitialStorage();
	    if (timeStep >= initialTimeSteps) Dew[i].SetSumWaterBalance();
	    if (timeStep == initialTimeSteps+numberTimeSteps-1) Dew[i].SetFinalStorage();
        }
        Dew[i].SetTotalReservoirStorage(timeStep, initialTimeSteps, numberTimeSteps, inputDataFound, firstTotal);
        if (*firstTotal)
        {
            *firstTotal = false;
        }
    }
}


void SetTimeSeriesInput(DistributedElement * const thisElement, ParametersGeneral * const ParGeneralStore,
                        MeteorologicalStations * const MetStations, int initialTimeSteps, int numberTimeSteps,
                        InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore,
                        int timeStep, bool * inputDataFound, bool * inputDataError)
{
    int indSta, indMet, indWeight, indexPrec, indexTemp;
    double precipitation, temperature, weight, elementPrecipitation, elementTemperature;
    //  double precGradLow, precGradHigh, lapseRate, lapseRateDry, lapseRateWet;
    double lapseRate, lapseRateDry, lapseRateWet, correctedTemperature;
    //  double elementElevation, gradElevation, metStationElevation, metStationWeight;
    double elementElevation, metStationElevation, metStationWeight;
    double rainSnowTemperature;
    //  gradElevation = ParGeneralStore->GetGRAD_CHANGE_ALT();
    //  precGradLow = ParGeneralStore->GetPREC_GRAD_LOW();
    //  precGradHigh = ParGeneralStore->GetPREC_GRAD_HIGH();
    lapseRateDry = ParGeneralStore->GetLAPSE_DRY();
    lapseRateWet = ParGeneralStore->GetLAPSE_WET();
    elementElevation = thisElement->GetElevation();
    // Precipitation data
    elementPrecipitation = 0.0;
    precipitation = 0.0;
    weight = 0.0;
    for (indexPrec = 0; indexPrec < ParGeneralStore->GetNUM_PREC_SERIES(); indexPrec++)
    {
        indSta = thisElement->GetMetSeriesNumber(indexPrec);
        indWeight = indexPrec;
        indMet = indSta;
        metStationElevation = MetStations->GetStationAltitude(indMet);
        metStationWeight = thisElement->GetMetSeriesWeight(indexPrec);
        precipitation = InputTimeSeriesStore->GetInput(timeStep, indMet);
        if (precipitation >= 0.0)
        {
            precipitation = precipitation / 1000.0;
            /*      if (gradElevation==0.0 || (elementElevation<=gradElevation && metStationElevation<=gradElevation))
            precipitation = (precipitation)*pow(precGradLow,(elementElevation-metStationElevation)/100.0)*metStationWeight;
            else if (elementElevation>gradElevation && metStationElevation<gradElevation)
            precipitation = (precipitation)*pow(precGradLow,(gradElevation-metStationElevation)/100.0)*
            pow(precGradHigh,(elementElevation-gradElevation)/100.0)*metStationWeight;
            else if (elementElevation<gradElevation && metStationElevation>gradElevation)
            precipitation = (precipitation)*pow(precGradLow,(elementElevation-gradElevation)/100.0)*
            pow(precGradHigh,(gradElevation-metStationElevation)/100.0)*metStationWeight;
            else
            precipitation = (precipitation)*pow(precGradHigh,(elementElevation-metStationElevation)/100.0)*metStationWeight;
            elementPrecipitation = elementPrecipitation + precipitation;*/
            elementPrecipitation = elementPrecipitation + InputElementStore->PrecipitationElevationCorrected(ParGeneralStore, precipitation, metStationElevation, elementElevation, metStationWeight);
            weight = weight + metStationWeight;
        }
        /*      cout << indexPrec << "  " << indSta << "  " << indMet << "  " << metStationWeight << "  " << weight << endl;
        if (indexPrec==ParGeneralStore->GetNUM_PREC_SERIES()-1) cout << endl;*/
    }
    if (weight > 0.0)
    {
        elementPrecipitation = elementPrecipitation / weight;
    }
    else
    {
        elementPrecipitation = missingData;
    }
    //    cout << indexPrec << "  " << metStationElevation << "  " << weight << "  " << elementPrecipitation << endl;
    /*    for (j=0; j<ParGeneralStore->GetNUM_PREC_SERIES(); j++) {
    cout << "P " << thisElement->GetMetSeriesNumber(j) << "  " << thisElement->GetMetSeriesWeight(j) << endl;
    }
    for (j=0; j<ParGeneralStore->GetNUM_TEMP_SERIES(); j++) {
    cout << "T " << thisElement->GetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES()+j) <<
    "  " << thisElement->GetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES()+j) << endl;
    }*/
    // Temperature data
    if (elementPrecipitation == 0.0)
    {
        lapseRate = lapseRateDry;
    }
    else
    {
        lapseRate = lapseRateWet;
    }
    elementTemperature = 0.0;
    temperature = 0.0;
    weight = 0.0;
    for (indexTemp = 0; indexTemp < ParGeneralStore->GetNUM_TEMP_SERIES(); indexTemp++)
    {
        indSta = thisElement->GetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES() + indexTemp);
        indWeight = ParGeneralStore->GetNUM_PREC_SERIES() + indexTemp;
        indMet = MetStations->GetNumPrecStations() + indSta;
        metStationElevation = MetStations->GetStationAltitude(indMet);
        metStationWeight = thisElement->GetMetSeriesWeight(indWeight);
        temperature = InputTimeSeriesStore->GetInput(timeStep, indMet);
        if (temperature > -99.0)
        {
            /*          temperature = (temperature/1.0 + lapseRate*(elementElevation-metStationElevation)/100.0)*metStationWeight;*/
            correctedTemperature = InputElementStore->TemperatureElevationCorrected(ParGeneralStore, temperature, elementPrecipitation, metStationElevation, elementElevation, metStationWeight);
            elementTemperature = elementTemperature + correctedTemperature;
            weight = weight + metStationWeight;
        }
        /*      cout << indexTemp << "  " << indSta << "  " << indMet << "  " << metStationWeight << "  " << weight << endl;
        if (indexTemp==ParGeneralStore->GetNUM_TEMP_SERIES()-1) cout << endl;*/
    }
    if (weight > 0.0)
    {
        elementTemperature = elementTemperature / weight;
    }
    else
    {
        elementTemperature = missingData;
    }
    if (elementPrecipitation > 0.0)
    {
        rainSnowTemperature = elementTemperature +
                              lapseRate * (elementElevation - thisElement->GetPrecStationsWeightedElevation()) / 100.0;
        //        lapseRate*(thisElement->GetTempStationsWeightedElevation()-thisElement->GetPrecStationsWeightedElevation())/100.0;
        if (rainSnowTemperature >= 0.0)
        {
            elementPrecipitation = elementPrecipitation * ParGeneralStore->GetPREC_CORR_RAIN();
        }
        else
        {
            elementPrecipitation = elementPrecipitation * ParGeneralStore->GetPREC_CORR_RAIN() * ParGeneralStore->GetPREC_CORR_SNOW();
        }
    }
    //    cout << indexPrec << "  "  << indexTemp << "  " << metStationElevation << "  " << weight << "  " << elementTemperature << endl;
    if (elementPrecipitation > missingData && elementTemperature > missingData)
    {
        *inputDataFound = true;
        InputElementStore->SetInput(0, elementPrecipitation * thisElement->GetPrecipitationCorrection());
        InputElementStore->SetInput(1, elementTemperature + thisElement->GetTemperatureCorrection());
        //      cout << elementPrecipitation*1000 << "  " << elementTemperature << endl;
        /*      if (dayNumber(InputTimeSeriesStore->GetYear(timeStep),InputTimeSeriesStore->GetMth(timeStep),
        InputTimeSeriesStore->GetDay(timeStep)) == ParGeneralStore->GetDAY_SNOW_ZERO()+
        leapYear(InputTimeSeriesStore->GetYear(timeStep))) {
        InputElementStore->SetInput(0,0.0);
        }*/
    }
    else
    {
        *inputDataFound = false;
        *inputDataError = true;
        InputElementStore->SetInput(0, missingData);
        InputElementStore->SetInput(1, missingData);
    }
}


void WaterBalanceGrid(DistributedElement * Dew, ParametersGeneral * ParGeneralStore, InputElement * InputElementStore,
                      int initialTimeSteps, int numberTimeSteps,
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy,
                      bool * inputDataFound, bool * firstTotal)
{
    ifstream filePrec, fileTemp;
    int i, j, k, jNew;
    unsigned short int metMissing = 10000;
    double precipitation, temperature;
    char fileName[100];
    char precFileName[100];
    char tempFileName[100];
    char hydYear[5];
    strcpy(precFileName, metPath);
    strcpy(tempFileName, metPath);
    strcat(precFileName, "/rr/");
    strcat(tempFileName, "/tm/");
    /*    if (datetime.getMonth() < 9)
          {*/
        sprintf(hydYear, "%04d", datetime.getYear());
    /*    }
          else
          {
          sprintf(hydYear, "%04d", datetime.getYear() + 1);
          }*/
    strcat(precFileName, hydYear);
    strcat(tempFileName, hydYear);
    sprintf(fileName, "/tm_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
    strcat(tempFileName, fileName);
    sprintf(fileName, "/rr_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
    strcat(precFileName, fileName);
    //  cout << precFileName << "  " << tempFileName << endl;

    filePrec.open(precFileName, ios::in | ios::binary);
    if (!filePrec.is_open())
    {
        cout << endl << "Error opening file " << precFileName << endl << endl;
        exit(1);
    }
    fileTemp.open(tempFileName, ios::in | ios::binary);
    if (!fileTemp.is_open())
    {
        cout << endl << "Error opening file " << tempFileName << endl << endl;
        exit(1);
    }
    //  filePrec.read((unsigned short int*) precip10, sizeof(unsigned short int)*nRows*nCols);
    //  fileTemp.read((unsigned short int*) temp10K, sizeof(unsigned short int)*nRows*nCols);
    streamoff newPosition;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    newPosition = ELEMENT(i, j) * (sizeof(unsigned short int));
                    filePrec.seekg(newPosition, ios::beg);
                    fileTemp.seekg(newPosition, ios::beg);
                    //          filePrec.read((unsigned short int *) &(precip10[ELEMENT(i,j)]), sizeof (unsigned short int));
                    //          fileTemp.read((unsigned short int *) &(temp10K[ELEMENT(i,j)]), sizeof (unsigned short int));
                    filePrec.read(reinterpret_cast<char *> (&(precip10[ELEMENT(i, j)])), sizeof (unsigned short int));
                    fileTemp.read(reinterpret_cast<char *> (&(temp10K[ELEMENT(i, j)])), sizeof (unsigned short int));
                    //      printf("%d  %hu  %hu\n",ELEMENT(i,j),precip10[ELEMENT(i,j)],temp10K[ELEMENT(i,j)]);
                    if (precip10[ELEMENT(i, j)] < metMissing && temp10K[ELEMENT(i, j)] < metMissing)
                    {
                        *inputDataFound = true;
                        precipitation = Dew[k].GetPrecipitationCorrection() * (double)precip10[ELEMENT(i, j)] / 10000.0;
                        temperature = Dew[k].GetTemperatureCorrection() + (double)(temp10K[ELEMENT(i, j)] - 2731) / 10.0;
                        //      printf("%d  %f  %f\n",ELEMENT(i,j),precipitation,temperature);
                        InputElementStore->SetInput(0, precipitation);
                        InputElementStore->SetInput(1, temperature);
                    }
                    else
                    {
                        jNew = j;
                        while ((precip10[ELEMENT(i, jNew)] >= metMissing || temp10K[ELEMENT(i, jNew)] >= metMissing) && jNew > 0)
                        {
                            jNew--;
                            newPosition = ELEMENT(i, jNew) * (sizeof(unsigned short int));
                            filePrec.seekg(newPosition, ios::beg);
                            fileTemp.seekg(newPosition, ios::beg);
                            filePrec.read(reinterpret_cast<char *> (&(precip10[ELEMENT(i, jNew)])), sizeof (unsigned short int));
                            fileTemp.read(reinterpret_cast<char *> (&(temp10K[ELEMENT(i, jNew)])), sizeof (unsigned short int));
                        }
                        if (precip10[ELEMENT(i, jNew)] < metMissing && temp10K[ELEMENT(i, jNew)] < metMissing)
                        {
                            *inputDataFound = true;
                            precipitation = Dew[k].GetPrecipitationCorrection() * (double)precip10[ELEMENT(i, jNew)] / 10000.0;
                            temperature = Dew[k].GetTemperatureCorrection() + (double)(temp10K[ELEMENT(i, jNew)] - 2731) / 10.0;
                            //      printf("%d  %f  %f\n",ELEMENT(i,jNew),precipitation,temperature);
                            InputElementStore->SetInput(0, precipitation);
                            InputElementStore->SetInput(1, temperature);
                        }
                        else
                        {
                            //      cout << endl << " Missing meterological data for: " << endl;
                            //      cout << precFileName << " or " << tempFileName << endl;
                            //      cout << " row = " << i << "  col = " << j << "  element no. = " << ELEMENT(i,j) << endl;
                            //      printf("  Precipitation %hu  Temperature %hu\n",precip10,temp10K);
                            //      cout << endl << endl;
                            //      exit (1);
                            //      *inputDataFound=false;
                            *inputDataFound = true;
                            InputElementStore->SetInput(0, 0.0);
                            InputElementStore->SetInput(1, 5.0);
                        }
                    }
                    // Water balance for landscape elements in case: no input to landscape element from upstream elements
                    if (!flowHierarchy)
                    {
                        Dew[k].WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, 0.0, 0.0);
			if (timeStep == initialTimeSteps-1) Dew[k].SetInitialStorage();
			if (timeStep >= initialTimeSteps) Dew[k].SetSumWaterBalance();
			if (timeStep == initialTimeSteps+numberTimeSteps-1) Dew[k].SetFinalStorage();
                        Dew[k].SetTotalReservoirStorage(timeStep, initialTimeSteps, numberTimeSteps, inputDataFound, firstTotal);
                        if (*firstTotal)
                        {
                            *firstTotal = false;
                        }
                    }
                    // End water balance for landscape elements in case: no input to landscape element from upstream elements
                    Dew[k].SetInputValue(0, InputElementStore->GetInput(0));
                    Dew[k].SetInputValue(1, InputElementStore->GetInput(1));
                    k++;
                }
            }
        }
    }
    filePrec.close();
    fileTemp.close();
}


void WaterBalanceGridVariable(DistributedElement * Dew, ParametersGeneral * ParGeneralStore, InputElement * InputElementStore,
                      int initialTimeSteps, int numberTimeSteps,
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy,
                      bool * inputDataFound, bool * firstTotal)
{
    ifstream filePrec, fileTemp;
    int g, h, i, j, k, jNew;
    unsigned short int metMissing = 10000;
    double precipitation, temperature;
    char fileName[100];
    char precFileName[100];
    char tempFileName[100];
    char hydYear[5];
    strcpy(precFileName, metPath);
    strcpy(tempFileName, metPath);
    strcat(precFileName, "/rr/");
    strcat(tempFileName, "/tm/");
    if (datetime.getMonth() < 9)
    {
        sprintf(hydYear, "%04d", datetime.getYear());
    }
    else
    {
        sprintf(hydYear, "%04d", datetime.getYear() + 1);
    }
    strcat(precFileName, hydYear);
    strcat(tempFileName, hydYear);
    sprintf(fileName, "/tm_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
    strcat(tempFileName, fileName);
    sprintf(fileName, "/rr_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
    strcat(precFileName, fileName);
    //  cout << precFileName << "  " << tempFileName << endl;

    filePrec.open(precFileName, ios::in | ios::binary);
    if (!filePrec.is_open())
    {
        cout << endl << "Error opening file " << precFileName << endl << endl;
        exit(1);
    }
    fileTemp.open(tempFileName, ios::in | ios::binary);
    if (!fileTemp.is_open())
    {
        cout << endl << "Error opening file " << tempFileName << endl << endl;
        exit(1);
    }
    //  filePrec.read((unsigned short int*) precip10, sizeof(unsigned short int)*nRows*nCols);
    //  fileTemp.read((unsigned short int*) temp10K, sizeof(unsigned short int)*nRows*nCols);
    for (i = 0; i < numRows1km; i++) 
      {
        for (j = 0; j < numCols1km; j++)
          {
            filePrec.read(reinterpret_cast<char *> (&(precip10[ELEMENT(i, j)])), sizeof (unsigned short int));
            fileTemp.read(reinterpret_cast<char *> (&(temp10K[ELEMENT(i, j)])), sizeof (unsigned short int));
          }
      }

    k = 0;
    //  for (i = 0; i < nRows; i++)
    for (g = 0; g < nRows; g++)
    {
        for (h = 0; h < nCols; h++)
	  //      for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(g, h))
                {
                    // Find row and column in SeNorge grid corresponding to landscape element

 		    i = (int)(numRows1km - 1 - floor((Dew[k].GetYCoord() - yllCorner1km)/cellSize1km));
		    j = (int)(floor((Dew[k].GetXCoord() - xllCorner1km)/cellSize1km));

                    if (precip10[ELEMENT(i, j)] < metMissing && temp10K[ELEMENT(i, j)] < metMissing)
                    {
                        *inputDataFound = true;
                        precipitation = Dew[k].GetPrecipitationCorrection() * (double)precip10[ELEMENT(i, j)] / 10000.0;
                        temperature = Dew[k].GetTemperatureCorrection() + (double)(temp10K[ELEMENT(i, j)] - 2731) / 10.0;
                        //      printf("%d  %f  %f\n",ELEMENT(i,j),precipitation,temperature);
                        InputElementStore->SetInput(0, precipitation);
                        InputElementStore->SetInput(1, temperature);
                    }
                    else
                    {
                        jNew = j;
                        while ((precip10[ELEMENT(i, jNew)] >= metMissing || temp10K[ELEMENT(i, jNew)] >= metMissing) && jNew > 0)
                        {
                            jNew--;
                        }
                        if (precip10[ELEMENT(i, jNew)] < metMissing && temp10K[ELEMENT(i, jNew)] < metMissing)
                        {
                            *inputDataFound = true;
                            precipitation = Dew[k].GetPrecipitationCorrection() * (double)precip10[ELEMENT(i, jNew)] / 10000.0;
                            temperature = Dew[k].GetTemperatureCorrection() + (double)(temp10K[ELEMENT(i, jNew)] - 2731) / 10.0;
                            //      printf("%d  %f  %f\n",ELEMENT(i,jNew),precipitation,temperature);
                            InputElementStore->SetInput(0, precipitation);
                            InputElementStore->SetInput(1, temperature);
                        }
                        else
                        {
                            //      cout << endl << " Missing meterological data for: " << endl;
                            //      cout << precFileName << " or " << tempFileName << endl;
                            //      cout << " row = " << i << "  col = " << j << "  element no. = " << ELEMENT(i,j) << endl;
                            //      printf("  Precipitation %hu  Temperature %hu\n",precip10,temp10K);
                            //      cout << endl << endl;
                            //      exit (1);
                            //      *inputDataFound=false;
                            *inputDataFound = true;
                            InputElementStore->SetInput(0, 0.0);
                            InputElementStore->SetInput(1, 5.0);
                        }
                    }
                    // Water balance for landscape elements in case: no input to landscape element from upstream elements
                    if (!flowHierarchy)
                    {
                        Dew[k].WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, 0.0, 0.0);
			if (timeStep == initialTimeSteps-1) Dew[k].SetInitialStorage();
			if (timeStep >= initialTimeSteps) Dew[k].SetSumWaterBalance();
			if (timeStep == initialTimeSteps+numberTimeSteps-1) Dew[k].SetFinalStorage();
                        Dew[k].SetTotalReservoirStorage(timeStep, initialTimeSteps, numberTimeSteps, inputDataFound, firstTotal);
                        if (*firstTotal)
                        {
                            *firstTotal = false;
                        }
                    }
                    // End water balance for landscape elements in case: no input to landscape element from upstream elements
                    Dew[k].SetInputValue(0, InputElementStore->GetInput(0));
                    Dew[k].SetInputValue(1, InputElementStore->GetInput(1));
                    k++;
                }
            }
        }
    }
    filePrec.close();
    fileTemp.close();
}


void WaterBalanceElements(DistributedElement * Dew, ParametersGeneral * ParGeneralStore, InputElement * InputElementStore,
                      int initialTimeSteps, int numberTimeSteps,
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metMask,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy,
                      bool * inputDataFound, bool * firstTotal)
{
    ifstream filePrec, fileTemp;
    int g, h, i, j, k, jNew;
    unsigned short int metMissing = 10000;
    double precipitation, temperature;
    char fileName[100];
    char precFileName[100];
    char tempFileName[100];
    char hydYear[5];

    strcpy(precFileName, metMask);
    strcpy(tempFileName, metMask);
    strcat(precFileName, "/Pre/");
    strcat(tempFileName, "/Tem/");
    sprintf(fileName, "/tem_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
    strcat(tempFileName, fileName);
    sprintf(fileName, "/pre_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
    strcat(precFileName, fileName);
    //  cout << precFileName << "  " << tempFileName << endl;

    filePrec.open(precFileName, ios::in | ios::binary);
    if (!filePrec.is_open())
    {
        cout << endl << "Error opening file " << precFileName << endl << endl;
        exit(1);
    }
    fileTemp.open(tempFileName, ios::in | ios::binary);
    if (!fileTemp.is_open())
    {
        cout << endl << "Error opening file " << tempFileName << endl << endl;
        exit(1);
    }
    for (k=0; k<numLand; k++)
      {
	filePrec.read(reinterpret_cast<char *> (&(precip10[k])), sizeof (unsigned short int));
	fileTemp.read(reinterpret_cast<char *> (&(temp10K[k])), sizeof (unsigned short int));
      }
    for (k=0; k<numLand; k++)
      {
	if (precip10[k] < metMissing && temp10K[k] < metMissing)
	  {
	    *inputDataFound = true;
	    precipitation = Dew[k].GetPrecipitationCorrection() * (double)precip10[k] / 10000.0;
	    temperature = Dew[k].GetTemperatureCorrection() + (double)(temp10K[k] - 2731) / 10.0;
	    //      printf("%d  %f  %f\n",ELEMENT(i,j),precipitation,temperature);
	    InputElementStore->SetInput(0, precipitation);
	    InputElementStore->SetInput(1, temperature);
	  }
	else
	  {
	    if (precip10[ELEMENT(i, j)] >= metMissing) 
	      {
		cout << "\n " << "Missing meteorological data " << precFileName << " index = " << k << "\n\n ";
		
	      }
	    if (temp10K[ELEMENT(i, j)] >= metMissing)
	      {
		cout << "\n " << "Missing meteorological data " << tempFileName << " index = " << k << "\n\n ";
		
	      }
	    //      exit (1);
	    //      *inputDataFound=false;
	    *inputDataFound = true;
	    InputElementStore->SetInput(0, 0.0);
	    InputElementStore->SetInput(1, 5.0);
	  }
	//      }
	// Water balance for landscape elements in case: no input to landscape element from upstream elements
	if (!flowHierarchy)
	  {
	    Dew[k].WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, 0.0, 0.0);
	    if (timeStep == initialTimeSteps-1) Dew[k].SetInitialStorage();
	    if (timeStep >= initialTimeSteps) Dew[k].SetSumWaterBalance();
	    if (timeStep == initialTimeSteps+numberTimeSteps-1) Dew[k].SetFinalStorage();
	    Dew[k].SetTotalReservoirStorage(timeStep, initialTimeSteps, numberTimeSteps, inputDataFound, firstTotal);
	    if (*firstTotal)
	      {
		*firstTotal = false;
	      }
	  }
	// End water balance for landscape elements in case: no input to landscape element from upstream elements
	Dew[k].SetInputValue(0, InputElementStore->GetInput(0));
	Dew[k].SetInputValue(1, InputElementStore->GetInput(1));
      }
    filePrec.close();
    fileTemp.close();
}


void GetRowCol(int elementNo, int nCols, int &row, int &col)
{
    row = elementNo / nCols;
    col = elementNo % nCols;
}


void TraverseOutletRoutingSubCatchment(SubCatchment * const thisSubCatchment, int numWatc, int * outletRoutingSubCatchments)
{
    bool outletRoutingFound=false;
    int i, j;

    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        outletRoutingFound=false;
        for (j = 0; j < numWatc; j++)
        {
	    if (!outletRoutingFound && outletRoutingSubCatchments[j] == thisSubCatchment->GetLakeNumber())
	    {
 	        outletRoutingFound=true;
		thisSubCatchment->SetLakeOutlet(thisSubCatchment->GetLakeNumber());
	    }
	    else
	    {
		outletRoutingSubCatchments[j] = -9999;
	    }
	}
	TraverseOutletRoutingSubCatchment(thisSubCatchment->GetUpStream(i), numWatc, outletRoutingSubCatchments);
    }
}


void SetElementCorrection(DistributedElement * const Dew, int numLand, int numberCorrectionCatchments, int * correctionCatchments, 
			  double * correctionPrecipitation, double * correctionTemperature, ofstream &fout)
{

  int i, k;
  double precCorr, tempCorr;

  if (numberCorrectionCatchments > 0) {
    precCorr = correctionPrecipitation[0];
    tempCorr = correctionTemperature[0];
  }
  else {
    precCorr = 1.0;
    tempCorr = 0.0;
  }

  for (k=0; k<numLand; k++) {
    Dew[k].SetPrecipitationCorrection(precCorr);
    Dew[k].SetTemperatureCorrection(tempCorr);
  }
}


void TraverseCorrectionSubCatchment(SubCatchment * const thisSubCatchment, int numberCorrectionCatchments,
                                    int * correctionCatchments, double * correctionPrecipitation, double * correctionTemperature)
{
    int i;
    double precCorr, tempCorr;
    DistributedElement * thisElement;

    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        TraverseCorrectionSubCatchment(thisSubCatchment->GetUpStream(i), numberCorrectionCatchments,
                                       correctionCatchments, correctionPrecipitation, correctionTemperature);
    }
    if (numberCorrectionCatchments > 0)
    {
        precCorr = correctionPrecipitation[0];
        tempCorr = correctionTemperature[0];
    }
    else
    {
        precCorr = 1.0;
        tempCorr = 0.0;
    }
    for (i = 0; i < numberCorrectionCatchments; i++)
    {
        if (thisSubCatchment->GetIdentifier() == correctionCatchments[i])
        {
            precCorr = correctionPrecipitation[i];
            tempCorr = correctionTemperature[i];
        }
    }
    thisElement = thisSubCatchment->GetLandScapeElement();
    while (thisElement)
    {
        TraverseCorrectionLandScape(thisElement, precCorr, tempCorr);
        thisElement = thisElement->GetNextElement();
    }
}


void TraverseCorrectionLandScape(DistributedElement * const thisElement, double precCorr, double tempCorr)
{
    int i;
    for (i = 0; i < thisElement->GetNumUpLand(); i++)
    {
        TraverseCorrectionLandScape(thisElement->GetUpLandFlow(i), precCorr, tempCorr);
    }
    thisElement->SetPrecipitationCorrection(precCorr);
    thisElement->SetTemperatureCorrection(tempCorr);
}


void TraverseSubCatchment(SubCatchment * const thisSubCatchment, ParametersGeneral * const ParGeneralStore,
                          MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
                          InputElement * InputElementStore, int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, char inputFormat,
                          bool flowHierarchy, bool forceDirect, bool * inputDataFound, bool * inputDataError, bool * firstTotal, ofstream &fout)
{
    int i;
    double accumulatedSum = 0.0;
    double accumulatedSumLake = 0.0;
    double accumulatedSumSnow = 0.0;
    double accumulatedSumGlacier = 0.0;
    double accumulatedSumHbv = 0.0;
    double subCatchmentSum = 0.0;
    double subCatchmentSumLake = 0.0;
    double subCatchmentSumSnow = 0.0;
    double subCatchmentSumGlacier = 0.0;
    double subCatchmentSumHbv = 0.0;
    double accumulatedDischarge = 0.0;
    double accumulatedUpperDischarge = 0.0;
    double accumulatedLowerDischarge = 0.0;
    double accumulatedInFlow = 0.0;
    double accumulatedUpperInFlow = 0.0;
    double accumulatedLowerInFlow = 0.0;
    double accumulatedPrecipitation = 0.0;
    double accumulatedTemperature = 0.0;
    double accumulatedLakeStorage = 0.0;
    double accumulatedSnowStore = 0.0;
    double accumulatedMeltWater = 0.0;
    double accumulatedSnowWaterEquivalentChange = 0.0;
    double accumulatedWaterOutput = 0.0;
    double accumulatedSnowCoverFraction = 0.0;
    double accumulatedGlacierSnowStore=0.0;
    double accumulatedGlacierSnowMeltWater=0.0;
    double accumulatedGlacierMassBalance = 0.0;
    double accumulatedAreaGlacierMassBalance = 0.0;
    double accumulatedGlacierIceMelt = 0.0;
    double accumulatedAnnualMassBalance = 0.0;
    double accumulatedGlacierIceVolume = 0.0;
    double accumulatedEvapotranspiration = 0.0;
    double accumulatedRunoff = 0.0;
    double accumulatedUpperRunoff = 0.0;
    double accumulatedLowerRunoff = 0.0;
    double accumulatedHbvSoilMoisture = 0.0;
    double accumulatedHbvSoilMoistureDeficit = 0.0;
    double accumulatedHbvPercSoilUpper = 0.0;
    double accumulatedHbvUpperZone = 0.0;
    double accumulatedHbvLowerZone = 0.0;
    //  double subCatchmentDischarge=0.0;
    double subCatchmentPrecipitation = 0.0;
    double subCatchmentTemperature = 0.0;
    double subCatchmentLakeStorage = 0.0;
    double subCatchmentSnowStore = 0.0;
    double subCatchmentMeltWater = 0.0;
    double subCatchmentSnowWaterEquivalentChange = 0.0;
    double subCatchmentWaterOutput = 0.0;
    double subCatchmentSnowCoverFraction = 0.0;
    double subCatchmentGlacierSnowStore=0.0;
    double subCatchmentGlacierSnowMeltWater=0.0;
    double subCatchmentGlacierMassBalance = 0.0;
    double subCatchmentAreaGlacierMassBalance = 0.0;
    double subCatchmentGlacierIceMelt = 0.0;
    double subCatchmentAnnualMassBalance = 0.0;
    double subCatchmentGlacierIceVolume = 0.0;
    double subCatchmentEvapotranspiration = 0.0;
    double subCatchmentRunoff = 0.0;
    double subCatchmentUpperRunoff = 0.0;
    double subCatchmentLowerRunoff = 0.0;
    double subCatchmentHbvSoilMoisture = 0.0;
    double subCatchmentHbvSoilMoistureDeficit = 0.0;
    double subCatchmentHbvPercSoilUpper = 0.0;
    double subCatchmentHbvUpperZone = 0.0;
    double subCatchmentHbvLowerZone = 0.0;
    DistributedElement * thisElement;
    //  double valueSub=0.0;
    //  double valueAcc=0.0;

    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        TraverseSubCatchment(thisSubCatchment->GetUpStream(i), ParGeneralStore, MetStations, InputTimeSeriesStore, InputElementStore,
                             timeStep, datetime, initialTimeSteps, numberTimeSteps, inputFormat, flowHierarchy, forceDirect, inputDataFound, inputDataError, firstTotal, fout);
        // Fluxes accumulated
        accumulatedDischarge = accumulatedDischarge + thisSubCatchment->GetUpStream(i)->GetAccumulatedDischarge(timeStep);
        accumulatedUpperDischarge = accumulatedUpperDischarge + thisSubCatchment->GetUpStream(i)->GetAccumulatedUpperDischarge(timeStep);
        accumulatedLowerDischarge = accumulatedLowerDischarge + thisSubCatchment->GetUpStream(i)->GetAccumulatedLowerDischarge(timeStep);
        accumulatedPrecipitation = accumulatedPrecipitation + thisSubCatchment->GetUpStream(i)->GetAccumulatedPrecipitation(timeStep);
        accumulatedTemperature = accumulatedTemperature + thisSubCatchment->GetUpStream(i)->GetAccumulatedTemperature(timeStep);
        accumulatedEvapotranspiration = accumulatedEvapotranspiration + thisSubCatchment->GetUpStream(i)->GetAccumulatedEvapotranspiration(timeStep);
        accumulatedRunoff = accumulatedRunoff + thisSubCatchment->GetUpStream(i)->GetAccumulatedRunoff(timeStep);
        accumulatedUpperRunoff = accumulatedUpperRunoff + thisSubCatchment->GetUpStream(i)->GetAccumulatedUpperRunoff(timeStep);
        accumulatedLowerRunoff = accumulatedLowerRunoff + thisSubCatchment->GetUpStream(i)->GetAccumulatedLowerRunoff(timeStep);
        accumulatedSum = accumulatedSum + thisSubCatchment->GetUpStream(i)->GetAccumulatedSum(timeStep);
        // Lake state variables accumulated
        accumulatedLakeStorage = accumulatedLakeStorage + thisSubCatchment->GetUpStream(i)->GetAccumulatedLakeStorage(timeStep);
        accumulatedSumLake = accumulatedSumLake + thisSubCatchment->GetUpStream(i)->GetAccumulatedSumLake(timeStep);
        // Snow state variables accumulated
        accumulatedSnowStore = accumulatedSnowStore + thisSubCatchment->GetUpStream(i)->GetAccumulatedSnowStore(timeStep);
        accumulatedMeltWater = accumulatedMeltWater + thisSubCatchment->GetUpStream(i)->GetAccumulatedMeltWater(timeStep);
        accumulatedSnowWaterEquivalentChange = accumulatedSnowWaterEquivalentChange + thisSubCatchment->GetUpStream(i)->GetAccumulatedSnowWaterEquivalentChange(timeStep);
        accumulatedWaterOutput = accumulatedWaterOutput + thisSubCatchment->GetUpStream(i)->GetAccumulatedWaterOutput(timeStep);
        accumulatedSnowCoverFraction = accumulatedSnowCoverFraction + thisSubCatchment->GetUpStream(i)->GetAccumulatedSnowCoverFraction(timeStep);
        accumulatedSumSnow = accumulatedSumSnow + thisSubCatchment->GetUpStream(i)->GetAccumulatedSumSnow(timeStep);
        // Glacier state variables accumulated
        //    cout << " 1 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << " " << thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierMassBalance(timeStep) << endl;
	accumulatedGlacierSnowStore=accumulatedGlacierSnowStore+thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierSnowStore(timeStep);
	accumulatedGlacierSnowMeltWater=accumulatedGlacierSnowMeltWater+thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierSnowMeltWater(timeStep);
        accumulatedGlacierMassBalance = accumulatedGlacierMassBalance + thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierMassBalance(timeStep);
        accumulatedAreaGlacierMassBalance = accumulatedAreaGlacierMassBalance + thisSubCatchment->GetUpStream(i)->GetAccumulatedAreaGlacierMassBalance(timeStep);
        accumulatedGlacierIceMelt = accumulatedGlacierIceMelt + thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierIceMelt(timeStep);
        accumulatedAnnualMassBalance = accumulatedAnnualMassBalance + thisSubCatchment->GetUpStream(i)->GetAccumulatedAnnualMassBalance(timeStep);
        accumulatedGlacierIceVolume = accumulatedGlacierIceVolume + thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierIceVolume(timeStep);
        accumulatedSumGlacier = accumulatedSumGlacier + thisSubCatchment->GetUpStream(i)->GetAccumulatedSumGlacier(timeStep);
        //    cout << " 1 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
        // Hbv state variables accumulated
        accumulatedHbvSoilMoisture = accumulatedHbvSoilMoisture + thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvSoilMoisture(timeStep);
        accumulatedHbvSoilMoistureDeficit = accumulatedHbvSoilMoistureDeficit + thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvSoilMoistureDeficit(timeStep);
        accumulatedHbvPercSoilUpper = accumulatedHbvPercSoilUpper + thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvPercSoilUpper(timeStep);
        accumulatedHbvUpperZone = accumulatedHbvUpperZone + thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvUpperZone(timeStep);
        accumulatedHbvLowerZone = accumulatedHbvLowerZone + thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvLowerZone(timeStep);
        accumulatedSumHbv = accumulatedSumHbv + thisSubCatchment->GetUpStream(i)->GetAccumulatedSumHbv(timeStep);
    }
    thisElement = thisSubCatchment->GetLandScapeElement();
    while (thisElement)
    {
        TraverseLandScape(thisElement, ParGeneralStore, MetStations, InputTimeSeriesStore, InputElementStore,
                          timeStep, datetime, initialTimeSteps, numberTimeSteps, inputFormat, flowHierarchy, forceDirect, inputDataFound, inputDataError, firstTotal, fout);
        // Fluxes accumulated
        accumulatedDischarge = accumulatedDischarge + thisElement->GetAccumulatedDischarge();
        accumulatedUpperDischarge = accumulatedUpperDischarge + thisElement->GetAccumulatedUpperDischarge();
        accumulatedLowerDischarge = accumulatedLowerDischarge + thisElement->GetAccumulatedLowerDischarge();
        accumulatedInFlow = accumulatedInFlow + thisElement->GetAccumulatedDischarge();
        accumulatedUpperInFlow = accumulatedUpperInFlow + thisElement->GetAccumulatedUpperDischarge();
        accumulatedLowerInFlow = accumulatedLowerInFlow + thisElement->GetAccumulatedLowerDischarge();
        accumulatedPrecipitation = accumulatedPrecipitation + thisElement->GetAccumulatedPrecipitation();
        accumulatedTemperature = accumulatedTemperature + thisElement->GetAccumulatedTemperature();
        accumulatedEvapotranspiration = accumulatedEvapotranspiration + thisElement->GetAccumulatedEvapotranspiration();
        accumulatedRunoff = accumulatedRunoff + thisElement->GetAccumulatedRunoff();
        accumulatedUpperRunoff = accumulatedUpperRunoff + thisElement->GetAccumulatedUpperRunoff();
        accumulatedLowerRunoff = accumulatedLowerRunoff + thisElement->GetAccumulatedLowerRunoff();
        accumulatedSum = accumulatedSum + thisElement->GetAccumulatedSum();
        //    subCatchmentDischarge=subCatchmentDischarge+thisElement->GetAccumulatedDischarge();
        subCatchmentPrecipitation = subCatchmentPrecipitation + thisElement->GetAccumulatedPrecipitation();
        subCatchmentTemperature = subCatchmentTemperature + thisElement->GetAccumulatedTemperature();
        subCatchmentEvapotranspiration = subCatchmentEvapotranspiration + thisElement->GetAccumulatedEvapotranspiration();
        subCatchmentRunoff = subCatchmentRunoff + thisElement->GetAccumulatedRunoff();
        subCatchmentUpperRunoff = subCatchmentUpperRunoff + thisElement->GetAccumulatedUpperRunoff();
        subCatchmentLowerRunoff = subCatchmentLowerRunoff + thisElement->GetAccumulatedLowerRunoff();
        subCatchmentSum = subCatchmentSum + thisElement->GetAccumulatedSum();
        // Lake state variables accumulated
        if (thisElement->GetAccumulatedSumLake() > 0)
        {
            accumulatedLakeStorage = accumulatedLakeStorage + thisElement->GetAccumulatedLakeStorage();
            accumulatedSumLake = accumulatedSumLake + thisElement->GetAccumulatedSumLake();
            subCatchmentLakeStorage = subCatchmentLakeStorage + thisElement->GetAccumulatedLakeStorage();
            subCatchmentSumLake = subCatchmentSumLake + thisElement->GetAccumulatedSumLake();
        }
        // Snow state variables accumulated
        if (thisElement->GetAccumulatedSumSnow() > 0)
        {
            accumulatedSnowStore = accumulatedSnowStore + thisElement->GetAccumulatedSnowStore();
            accumulatedMeltWater = accumulatedMeltWater + thisElement->GetAccumulatedMeltWater();
            accumulatedSnowWaterEquivalentChange = accumulatedSnowWaterEquivalentChange + thisElement->GetAccumulatedSnowWaterEquivalentChange();
            accumulatedWaterOutput = accumulatedWaterOutput + thisElement->GetAccumulatedWaterOutput();
            accumulatedSnowCoverFraction = accumulatedSnowCoverFraction + thisElement->GetAccumulatedSnowCoverFraction();
            accumulatedSumSnow = accumulatedSumSnow + thisElement->GetAccumulatedSumSnow();
            subCatchmentSnowStore = subCatchmentSnowStore + thisElement->GetAccumulatedSnowStore();
            subCatchmentMeltWater = subCatchmentMeltWater + thisElement->GetAccumulatedMeltWater();
            subCatchmentSnowWaterEquivalentChange = subCatchmentSnowWaterEquivalentChange + thisElement->GetAccumulatedSnowWaterEquivalentChange();
            subCatchmentWaterOutput = subCatchmentWaterOutput + thisElement->GetAccumulatedWaterOutput();
            subCatchmentSnowCoverFraction = subCatchmentSnowCoverFraction + thisElement->GetAccumulatedSnowCoverFraction();
            subCatchmentSumSnow = subCatchmentSumSnow + thisElement->GetAccumulatedSumSnow();
        }
        // Glacier state variables accumulated
        if (thisElement->GetAccumulatedSumGlacier() > 0)
        {
            //        cout << " 2 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
	    accumulatedGlacierSnowStore=accumulatedGlacierSnowStore+thisElement->GetAccumulatedGlacierSnowStore();
            accumulatedGlacierSnowMeltWater=accumulatedGlacierSnowMeltWater+thisElement->GetAccumulatedGlacierSnowMeltWater();
            accumulatedGlacierMassBalance = accumulatedGlacierMassBalance + thisElement->GetAccumulatedGlacierMassBalance();
            accumulatedAreaGlacierMassBalance = accumulatedAreaGlacierMassBalance + thisElement->GetAccumulatedAreaGlacierMassBalance();
            accumulatedGlacierIceMelt = accumulatedGlacierIceMelt + thisElement->GetAccumulatedGlacierIceMelt();
            accumulatedAnnualMassBalance = accumulatedAnnualMassBalance + thisElement->GetAccumulatedAnnualMassBalance();
            accumulatedGlacierIceVolume = accumulatedGlacierIceVolume + thisElement->GetAccumulatedGlacierIceVolume();
            accumulatedSumGlacier = accumulatedSumGlacier + thisElement->GetAccumulatedSumGlacier();
            //        cout << " 2 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
	    subCatchmentGlacierSnowStore=subCatchmentGlacierSnowStore+thisElement->GetAccumulatedGlacierSnowStore();
	    subCatchmentGlacierSnowMeltWater=subCatchmentGlacierSnowMeltWater+thisElement->GetAccumulatedGlacierSnowMeltWater();
            subCatchmentGlacierMassBalance = subCatchmentGlacierMassBalance + thisElement->GetAccumulatedGlacierMassBalance();
            subCatchmentAreaGlacierMassBalance = subCatchmentAreaGlacierMassBalance + thisElement->GetAccumulatedAreaGlacierMassBalance();
            subCatchmentGlacierIceMelt = subCatchmentGlacierIceMelt + thisElement->GetAccumulatedGlacierIceMelt();
            subCatchmentAnnualMassBalance = subCatchmentAnnualMassBalance + thisElement->GetAccumulatedAnnualMassBalance();
            subCatchmentGlacierIceVolume = subCatchmentGlacierIceVolume + thisElement->GetAccumulatedGlacierIceVolume();
            subCatchmentSumGlacier = subCatchmentSumGlacier + thisElement->GetAccumulatedSumGlacier();
        }
        // Hbv state variables accumulated
        if (thisElement->GetAccumulatedSumHbv() > 0)
        {
            accumulatedHbvSoilMoisture = accumulatedHbvSoilMoisture + thisElement->GetAccumulatedHbvSoilMoisture();
            accumulatedHbvSoilMoistureDeficit = accumulatedHbvSoilMoistureDeficit + thisElement->GetAccumulatedHbvSoilMoistureDeficit();
            accumulatedHbvPercSoilUpper = accumulatedHbvPercSoilUpper + thisElement->GetAccumulatedHbvPercSoilUpper();
            accumulatedHbvUpperZone = accumulatedHbvUpperZone + thisElement->GetAccumulatedHbvUpperZone();
            accumulatedHbvLowerZone = accumulatedHbvLowerZone + thisElement->GetAccumulatedHbvLowerZone();
            accumulatedSumHbv = accumulatedSumHbv + thisElement->GetAccumulatedSumHbv();
            subCatchmentHbvSoilMoisture = subCatchmentHbvSoilMoisture + thisElement->GetAccumulatedHbvSoilMoisture();
            subCatchmentHbvSoilMoistureDeficit = subCatchmentHbvSoilMoistureDeficit + thisElement->GetAccumulatedHbvSoilMoistureDeficit();
            subCatchmentHbvPercSoilUpper = subCatchmentHbvPercSoilUpper + thisElement->GetAccumulatedHbvPercSoilUpper();
            subCatchmentHbvUpperZone = subCatchmentHbvUpperZone + thisElement->GetAccumulatedHbvUpperZone();
            subCatchmentHbvLowerZone = subCatchmentHbvLowerZone + thisElement->GetAccumulatedHbvLowerZone();
            subCatchmentSumHbv = subCatchmentSumHbv + thisElement->GetAccumulatedSumHbv();
        }
        thisElement = thisElement->GetNextElement();
    }
    // Fluxes accumulated
    thisSubCatchment->SetAccumulatedDischarge(timeStep, accumulatedDischarge);
    thisSubCatchment->SetAccumulatedUpperDischarge(timeStep, accumulatedUpperDischarge);
    thisSubCatchment->SetAccumulatedLowerDischarge(timeStep, accumulatedLowerDischarge);
    thisSubCatchment->SetAccumulatedInFlow(timeStep, accumulatedInFlow);
    //  thisSubCatchment->SetAccumulatedUpperInFlow(timeStep, accumulatedUpperInFlow);
    //  thisSubCatchment->SetAccumulatedLowerInFlow(timeStep, accumulatedLowerInFlow);
    if (accumulatedInFlow > accumulatedDischarge)
    {
        cout << endl << thisSubCatchment->GetIdentifier() << " accumulatedInFlow > accumulatedDischarge  " << accumulatedInFlow << "  " << accumulatedDischarge << endl << endl;
        //    exit(1);
    }
    thisSubCatchment->SetAccumulatedPrecipitation(timeStep, accumulatedPrecipitation);
    thisSubCatchment->SetAccumulatedTemperature(timeStep, accumulatedTemperature);
    thisSubCatchment->SetAccumulatedEvapotranspiration(timeStep, accumulatedEvapotranspiration);
    thisSubCatchment->SetAccumulatedRunoff(timeStep, accumulatedRunoff);
    thisSubCatchment->SetAccumulatedUpperRunoff(timeStep, accumulatedUpperRunoff);
    thisSubCatchment->SetAccumulatedLowerRunoff(timeStep, accumulatedLowerRunoff);
    thisSubCatchment->SetAccumulatedSum(timeStep, accumulatedSum);
    thisSubCatchment->SetSubCatchmentPrecipitation(timeStep, subCatchmentPrecipitation / subCatchmentSum);
    thisSubCatchment->SetSubCatchmentTemperature(timeStep, subCatchmentTemperature / subCatchmentSum);
    thisSubCatchment->SetSubCatchmentEvapotranspiration(timeStep, subCatchmentEvapotranspiration / subCatchmentSum);
    if (flowHierarchy && !forceDirect)
    {
        thisSubCatchment->SetSubCatchmentRunoff(timeStep, accumulatedInFlow * ParGeneralStore->GetSECONDS_TIMESTEP() / subCatchmentSum); /*  discharge (m3/s) -> runoff (m)  */
        thisSubCatchment->SetSubCatchmentUpperRunoff(timeStep, accumulatedUpperInFlow * ParGeneralStore->GetSECONDS_TIMESTEP() / subCatchmentSum); /*  discharge (m3/s) -> runoff (m)  */
        thisSubCatchment->SetSubCatchmentLowerRunoff(timeStep, accumulatedLowerInFlow * ParGeneralStore->GetSECONDS_TIMESTEP() / subCatchmentSum); /*  discharge (m3/s) -> runoff (m)  */
    }
    else
    {
        thisSubCatchment->SetSubCatchmentRunoff(timeStep, subCatchmentRunoff / subCatchmentSum);
        thisSubCatchment->SetSubCatchmentUpperRunoff(timeStep, subCatchmentUpperRunoff / subCatchmentSum);
        thisSubCatchment->SetSubCatchmentLowerRunoff(timeStep, subCatchmentLowerRunoff / subCatchmentSum);
    }
    // Lake state variables accumulated
    if (accumulatedSumLake > 0)
    {
        thisSubCatchment->SetAccumulatedLakeStorage(timeStep, accumulatedLakeStorage);
        thisSubCatchment->SetAccumulatedSumLake(timeStep, accumulatedSumLake);
    }
    else
    {
        thisSubCatchment->SetAccumulatedLakeStorage(timeStep, accumulatedLakeStorage);
        thisSubCatchment->SetAccumulatedSumLake(timeStep, accumulatedSumLake);
    }
    if (subCatchmentSumLake > 0)
    {
        thisSubCatchment->SetSubCatchmentLakeStorage(timeStep, subCatchmentLakeStorage/subCatchmentSumLake);
        //    thisSubCatchment->SetSubCatchmentLakeStorage(timeStep, subCatchmentLakeStorage / subCatchmentSum);
    }
    else
    {
        thisSubCatchment->SetSubCatchmentLakeStorage(timeStep, missingData);
    }
    // Snow state variables accumulated
    if (accumulatedSumSnow > 0)
    {
        //    if (timeStep > 1) valueAcc=thisSubCatchment->GetAccumulatedSnowWaterEquivalentChange(timeStep-1);
        thisSubCatchment->SetAccumulatedSnowStore(timeStep, accumulatedSnowStore);
        thisSubCatchment->SetAccumulatedMeltWater(timeStep, accumulatedMeltWater);
        thisSubCatchment->SetAccumulatedSnowWaterEquivalentChange(timeStep, accumulatedSnowWaterEquivalentChange);
        //    thisSubCatchment->SetAccumulatedSnowWaterEquivalentChange(timeStep, valueAcc+accumulatedSnowWaterEquivalentChange);
        thisSubCatchment->SetAccumulatedWaterOutput(timeStep, accumulatedWaterOutput);
        thisSubCatchment->SetAccumulatedSnowCoverFraction(timeStep, accumulatedSnowCoverFraction);
        thisSubCatchment->SetAccumulatedSumSnow(timeStep, accumulatedSumSnow);
    }
    else
    {
        thisSubCatchment->SetAccumulatedSnowStore(timeStep, accumulatedSnowStore);
        thisSubCatchment->SetAccumulatedMeltWater(timeStep, accumulatedMeltWater);
        thisSubCatchment->SetAccumulatedSnowWaterEquivalentChange(timeStep, accumulatedSnowWaterEquivalentChange);
        thisSubCatchment->SetAccumulatedWaterOutput(timeStep, accumulatedWaterOutput);
        thisSubCatchment->SetAccumulatedSnowCoverFraction(timeStep, accumulatedSnowCoverFraction);
        thisSubCatchment->SetAccumulatedSumSnow(timeStep, accumulatedSumSnow);
    }
    if (subCatchmentSumSnow > 0)
    {
        //    if (timeStep > 1) valueSub=thisSubCatchment->GetSubCatchmentSnowWaterEquivalentChange(timeStep-1);
        thisSubCatchment->SetSubCatchmentSnowStore(timeStep, subCatchmentSnowStore/subCatchmentSumSnow);
        thisSubCatchment->SetSubCatchmentMeltWater(timeStep, subCatchmentMeltWater/subCatchmentSumSnow);
        thisSubCatchment->SetSubCatchmentSnowWaterEquivalentChange(timeStep, subCatchmentSnowWaterEquivalentChange/subCatchmentSumSnow);
        //    thisSubCatchment->SetSubCatchmentSnowWaterEquivalentChange(timeStep, valueSub+subCatchmentSnowWaterEquivalentChange/subCatchmentSumSnow);
        thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, subCatchmentWaterOutput/subCatchmentSumSnow);
        thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, subCatchmentSnowCoverFraction/subCatchmentSumSnow);
        //    thisSubCatchment->SetSubCatchmentSnowStore(timeStep, subCatchmentSnowStore / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentMeltWater(timeStep, subCatchmentMeltWater / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentSnowWaterEquivalentChange(timeStep, subCatchmentSnowWaterEquivalentChange / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, subCatchmentWaterOutput / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, subCatchmentSnowCoverFraction / subCatchmentSum);
    }
    else
    {
        thisSubCatchment->SetSubCatchmentSnowStore(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentMeltWater(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentSnowWaterEquivalentChange(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, missingData);
    }
    // Glacier state variables accumulated
    if (accumulatedSumGlacier > 0)
    {
        //    cout << " 3 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
        thisSubCatchment->SetAccumulatedGlacierSnowStore(timeStep, accumulatedGlacierSnowStore);
        thisSubCatchment->SetAccumulatedGlacierSnowMeltWater(timeStep, accumulatedGlacierSnowMeltWater);
        thisSubCatchment->SetAccumulatedGlacierMassBalance(timeStep, accumulatedGlacierMassBalance);
        thisSubCatchment->SetAccumulatedAreaGlacierMassBalance(timeStep, accumulatedAreaGlacierMassBalance);
        thisSubCatchment->SetAccumulatedGlacierIceMelt(timeStep, accumulatedGlacierIceMelt);
        thisSubCatchment->SetAccumulatedAnnualMassBalance(timeStep, accumulatedAnnualMassBalance);
        thisSubCatchment->SetAccumulatedGlacierIceVolume(timeStep, accumulatedGlacierIceVolume);
        thisSubCatchment->SetAccumulatedSumGlacier(timeStep, accumulatedSumGlacier);
        //    cout << " 3 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
    }
    else
    {
        //    cout << " 4 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
        thisSubCatchment->SetAccumulatedGlacierSnowStore(timeStep, accumulatedGlacierSnowStore);
        thisSubCatchment->SetAccumulatedGlacierSnowMeltWater(timeStep, accumulatedGlacierSnowMeltWater);
        thisSubCatchment->SetAccumulatedGlacierMassBalance(timeStep, accumulatedGlacierMassBalance);
        thisSubCatchment->SetAccumulatedAreaGlacierMassBalance(timeStep, accumulatedAreaGlacierMassBalance);
        thisSubCatchment->SetAccumulatedGlacierIceMelt(timeStep, accumulatedGlacierIceMelt);
        thisSubCatchment->SetAccumulatedAnnualMassBalance(timeStep, accumulatedAnnualMassBalance);
        thisSubCatchment->SetAccumulatedGlacierIceVolume(timeStep, accumulatedGlacierIceVolume);
        thisSubCatchment->SetAccumulatedSumGlacier(timeStep, accumulatedSumGlacier);
        //    cout << " 4 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
    }
    if (subCatchmentSumGlacier > 0)
    {
        thisSubCatchment->SetSubCatchmentAreaGlacier(timeStep, subCatchmentSumGlacier);
        thisSubCatchment->SetSubCatchmentGlacierSnowStore(timeStep, subCatchmentGlacierSnowStore/subCatchmentSumGlacier);
        thisSubCatchment->SetSubCatchmentGlacierSnowMeltWater(timeStep, subCatchmentGlacierSnowMeltWater/subCatchmentSumGlacier);
        thisSubCatchment->SetSubCatchmentAreaGlacierMassBalance(timeStep, subCatchmentGlacierMassBalance/subCatchmentSum);
        thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, subCatchmentGlacierMassBalance/subCatchmentSumGlacier);
        thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, subCatchmentGlacierIceMelt/subCatchmentSum);
        thisSubCatchment->SetSubCatchmentAnnualMassBalance(timeStep, subCatchmentAnnualMassBalance/subCatchmentSumGlacier);
        thisSubCatchment->SetSubCatchmentGlacierIceVolume(timeStep, subCatchmentGlacierIceVolume);
        //    thisSubCatchment->SetSubCatchmentGlacierSnowStore(timeStep, subCatchmentGlacierSnowStore/subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentGlacierSnowMeltWater(timeStep, subCatchmentGlacierSnowMeltWater/subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, subCatchmentGlacierMassBalance / subCatchmentSumGlacier);
        //    thisSubCatchment->SetSubCatchmentAreaGlacierMassBalance(timeStep, subCatchmentGlacierMassBalance / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, subCatchmentGlacierIceMelt / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentAnnualMassBalance(timeStep, subCatchmentAnnualMassBalance / subCatchmentSumGlacier);
        //    thisSubCatchment->SetSubCatchmentGlacierIceVolume(timeStep, subCatchmentGlacierIceVolume);
    }
    else
    {
        thisSubCatchment->SetSubCatchmentAreaGlacier(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentGlacierSnowStore(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentGlacierSnowMeltWater(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentAreaGlacierMassBalance(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentAnnualMassBalance(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentGlacierIceVolume(timeStep, missingData);
    }
    // Hbv state variables accumulated
    if (accumulatedSumHbv > 0)
    {
        thisSubCatchment->SetAccumulatedHbvSoilMoisture(timeStep, accumulatedHbvSoilMoisture);
        thisSubCatchment->SetAccumulatedHbvSoilMoistureDeficit(timeStep, accumulatedHbvSoilMoistureDeficit);
        thisSubCatchment->SetAccumulatedHbvPercSoilUpper(timeStep, accumulatedHbvPercSoilUpper);
        thisSubCatchment->SetAccumulatedHbvUpperZone(timeStep, accumulatedHbvUpperZone);
        thisSubCatchment->SetAccumulatedHbvLowerZone(timeStep, accumulatedHbvLowerZone);
        thisSubCatchment->SetAccumulatedSumHbv(timeStep, accumulatedSumHbv);
    }
    else
    {
        thisSubCatchment->SetAccumulatedHbvSoilMoisture(timeStep, accumulatedHbvSoilMoisture);
        thisSubCatchment->SetAccumulatedHbvSoilMoistureDeficit(timeStep, accumulatedHbvSoilMoistureDeficit);
        thisSubCatchment->SetAccumulatedHbvPercSoilUpper(timeStep, accumulatedHbvPercSoilUpper);
        thisSubCatchment->SetAccumulatedHbvUpperZone(timeStep, accumulatedHbvUpperZone);
        thisSubCatchment->SetAccumulatedHbvLowerZone(timeStep, accumulatedHbvLowerZone);
        thisSubCatchment->SetAccumulatedSumHbv(timeStep, accumulatedSumHbv);
    }
    if (subCatchmentSumHbv > 0)
    {
        thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, subCatchmentHbvSoilMoisture/subCatchmentSumHbv);
        thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, subCatchmentHbvSoilMoistureDeficit/subCatchmentSumHbv);
        thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, subCatchmentHbvPercSoilUpper/subCatchmentSumHbv);
        thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, subCatchmentHbvUpperZone/subCatchmentSumHbv);
        thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, subCatchmentHbvLowerZone/subCatchmentSumHbv);
        //    thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, subCatchmentHbvSoilMoisture / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, subCatchmentHbvSoilMoistureDeficit / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, subCatchmentHbvPercSoilUpper / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, subCatchmentHbvUpperZone / subCatchmentSum);
        //    thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, subCatchmentHbvLowerZone / subCatchmentSum);
    }
    else
    {
        thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, missingData);
        thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, missingData);
    }
    //  cout << thisSubCatchment->GetIdentifier() << "  " << timeStep << "  " << "  Discharge: "
    //       << accumulatedDischarge << endl;
    //       << accumulatedDischarge << "  " << thisSubCatchment->GetAccumulatedDischarge(timeStep) << endl;
}


void TraverseLandScape(DistributedElement * const thisElement, ParametersGeneral * const ParGeneralStore,
                       MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
                       InputElement * InputElementStore, int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, char inputFormat,
                       bool flowHierarchy, bool forceDirect, bool * inputDataFound, bool * inputDataError, bool * firstTotal, ofstream &fout)
{
    int i, j;
    double accumulatedSum = 0.0;
    double accumulatedSumLake = 0.0;
    double accumulatedSumSnow = 0.0;
    double accumulatedSumGlacier = 0.0;
    double accumulatedSumHbv = 0.0;
    double accumulatedDischarge = 0.0;
    double accumulatedLowerDischarge = 0.0;
    double accumulatedUpperDischarge = 0.0;
    double accumulatedPrecipitation = 0.0;
    double accumulatedTemperature = 0.0;
    double accumulatedLakeStorage = 0.0;
    double accumulatedSnowStore = 0.0;
    double accumulatedMeltWater = 0.0;
    double accumulatedSnowWaterEquivalentChange = 0.0;
    double accumulatedWaterOutput = 0.0;
    double accumulatedSnowCoverFraction = 0.0;
    double accumulatedGlacierSnowStore=0.0;
    double accumulatedGlacierSnowMeltWater=0.0;
    double accumulatedGlacierMassBalance = 0.0;
    double accumulatedAreaGlacierMassBalance = 0.0;
    double accumulatedGlacierIceMelt = 0.0;
    double accumulatedAnnualMassBalance = 0.0;
    double accumulatedGlacierIceVolume = 0.0;
    double accumulatedEvapotranspiration = 0.0;
    double accumulatedRunoff = 0.0;
    double accumulatedUpperRunoff = 0.0;
    double accumulatedLowerRunoff = 0.0;
    double accumulatedHbvSoilMoisture = 0.0;
    double accumulatedHbvSoilMoistureDeficit = 0.0;
    double accumulatedHbvPercSoilUpper = 0.0;
    double accumulatedHbvUpperZone = 0.0;
    double accumulatedHbvLowerZone = 0.0;
    // cout << "start TraverselandScape " << thisElement->GetLandIndex() << endl;
    for (i = 0; i < thisElement->GetNumUpLand(); i++)
    {
        TraverseLandScape(thisElement->GetUpLandFlow(i), ParGeneralStore, MetStations, InputTimeSeriesStore, InputElementStore,
                          timeStep, datetime, initialTimeSteps, numberTimeSteps, inputFormat, flowHierarchy, forceDirect, inputDataFound, inputDataError, firstTotal, fout);
        if (thisElement->GetUpLandFlow(i)->GetAccumulatedPrecipitation() != missingData)
        {
            // Fluxes accumulated
            accumulatedDischarge = accumulatedDischarge + thisElement->GetUpLandFlow(i)->GetAccumulatedDischarge();
            accumulatedLowerDischarge = accumulatedLowerDischarge + thisElement->GetUpLandFlow(i)->GetAccumulatedLowerDischarge();
            accumulatedUpperDischarge = accumulatedUpperDischarge + thisElement->GetUpLandFlow(i)->GetAccumulatedUpperDischarge();
            accumulatedPrecipitation = accumulatedPrecipitation + thisElement->GetUpLandFlow(i)->GetAccumulatedPrecipitation();
            accumulatedTemperature = accumulatedTemperature + thisElement->GetUpLandFlow(i)->GetAccumulatedTemperature();
            accumulatedEvapotranspiration = accumulatedEvapotranspiration + thisElement->GetUpLandFlow(i)->GetAccumulatedEvapotranspiration();
            accumulatedRunoff = accumulatedRunoff + thisElement->GetUpLandFlow(i)->GetAccumulatedRunoff();
            accumulatedUpperRunoff = accumulatedUpperRunoff + thisElement->GetUpLandFlow(i)->GetAccumulatedUpperRunoff();
            accumulatedLowerRunoff = accumulatedLowerRunoff + thisElement->GetUpLandFlow(i)->GetAccumulatedLowerRunoff();
            accumulatedSum = accumulatedSum + thisElement->GetUpLandFlow(i)->GetAccumulatedSum();
            // Lake state variables accumulated
            accumulatedLakeStorage = accumulatedLakeStorage + thisElement->GetUpLandFlow(i)->GetAccumulatedLakeStorage();
            accumulatedSumLake = accumulatedSumLake + thisElement->GetUpLandFlow(i)->GetAccumulatedSumLake();
            // Snow state variables accumulated
            accumulatedSnowStore = accumulatedSnowStore + thisElement->GetUpLandFlow(i)->GetAccumulatedSnowStore();
            accumulatedMeltWater = accumulatedMeltWater + thisElement->GetUpLandFlow(i)->GetAccumulatedMeltWater();
            accumulatedSnowWaterEquivalentChange = accumulatedSnowWaterEquivalentChange + thisElement->GetUpLandFlow(i)->GetAccumulatedSnowWaterEquivalentChange();
            accumulatedWaterOutput = accumulatedWaterOutput + thisElement->GetUpLandFlow(i)->GetAccumulatedWaterOutput();
            accumulatedSnowCoverFraction = accumulatedSnowCoverFraction + thisElement->GetUpLandFlow(i)->GetAccumulatedSnowCoverFraction();
            accumulatedSumSnow = accumulatedSumSnow + thisElement->GetUpLandFlow(i)->GetAccumulatedSumSnow();
            // Glacier state variables accumulated
            //      cout << " 5 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
            accumulatedGlacierSnowStore=accumulatedGlacierSnowStore+thisElement->GetUpLandFlow(i)->GetAccumulatedGlacierSnowStore();
            accumulatedGlacierSnowMeltWater=accumulatedGlacierSnowMeltWater+thisElement->GetUpLandFlow(i)->GetAccumulatedGlacierSnowMeltWater();
            accumulatedGlacierMassBalance = accumulatedGlacierMassBalance + thisElement->GetUpLandFlow(i)->GetAccumulatedGlacierMassBalance();
            accumulatedAreaGlacierMassBalance = accumulatedAreaGlacierMassBalance + thisElement->GetUpLandFlow(i)->GetAccumulatedAreaGlacierMassBalance();
            accumulatedGlacierIceMelt = accumulatedGlacierIceMelt + thisElement->GetUpLandFlow(i)->GetAccumulatedGlacierIceMelt();
            accumulatedAnnualMassBalance = accumulatedAnnualMassBalance + thisElement->GetUpLandFlow(i)->GetAccumulatedAnnualMassBalance();
            accumulatedGlacierIceVolume = accumulatedGlacierIceVolume + thisElement->GetUpLandFlow(i)->GetAccumulatedGlacierIceVolume();
            accumulatedSumGlacier = accumulatedSumGlacier + thisElement->GetUpLandFlow(i)->GetAccumulatedSumGlacier();
            //      cout << " 5 " << accumulatedSumGlacier+thisElement->GetUpLandFlow(i)->GetAccumulatedSumGlacier()
            //           << " " << accumulatedGlacierMassBalance+thisElement->GetUpLandFlow(i)->GetAccumulatedGlacierMassBalance() << endl;
            // Hbv state variables accumulated
            accumulatedHbvSoilMoisture = accumulatedHbvSoilMoisture + thisElement->GetUpLandFlow(i)->GetAccumulatedHbvSoilMoisture();
            accumulatedHbvSoilMoistureDeficit = accumulatedHbvSoilMoistureDeficit + thisElement->GetUpLandFlow(i)->GetAccumulatedHbvSoilMoistureDeficit();
            accumulatedHbvPercSoilUpper = accumulatedHbvPercSoilUpper + thisElement->GetUpLandFlow(i)->GetAccumulatedHbvPercSoilUpper();
            accumulatedHbvUpperZone = accumulatedHbvUpperZone + thisElement->GetUpLandFlow(i)->GetAccumulatedHbvUpperZone();
            accumulatedHbvLowerZone = accumulatedHbvLowerZone + thisElement->GetUpLandFlow(i)->GetAccumulatedHbvLowerZone();
            accumulatedSumHbv = accumulatedSumHbv + thisElement->GetUpLandFlow(i)->GetAccumulatedSumHbv();
        }
    }

    // ** Algorithm to be performed in case: input to landscape element from upstream elements
    if (flowHierarchy)
    {
        if (inputFormat == 'T' || inputFormat == 't')
        {
            *inputDataFound = false;
            SetTimeSeriesInput(thisElement, ParGeneralStore, MetStations, initialTimeSteps, numberTimeSteps,
                               InputTimeSeriesStore, InputElementStore, timeStep, inputDataFound, inputDataError);
            if (*inputDataFound)
            {
                if (!forceDirect)
                {
                    thisElement->WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, accumulatedLowerDischarge, accumulatedUpperDischarge);
                }
                else
                {
                    thisElement->WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, 0.0, 0.0);
                }
            }
            thisElement->SetTotalReservoirStorage(timeStep, initialTimeSteps, numberTimeSteps, inputDataFound, firstTotal);
            if (*firstTotal)
            {
                *firstTotal = false;
            }
        }
        else if (inputFormat == 'G' || inputFormat == 'g')
        {
            InputElementStore->SetInput(0, thisElement->GetInputValue(0));
            InputElementStore->SetInput(1, thisElement->GetInputValue(1));
            *inputDataFound = true;
            if (!forceDirect)
            {
                thisElement->WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, accumulatedLowerDischarge, accumulatedUpperDischarge);
            }
            else
            {
                thisElement->WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps, 0.0, 0.0);
            }
            thisElement->SetTotalReservoirStorage(timeStep, initialTimeSteps, numberTimeSteps, inputDataFound, firstTotal);
            if (*firstTotal)
            {
                *firstTotal = false;
            }
        }
	if (timeStep == initialTimeSteps-1) thisElement->SetInitialStorage();
	if (timeStep >= initialTimeSteps) thisElement->SetSumWaterBalance();
	if (timeStep == initialTimeSteps+numberTimeSteps-1) thisElement->SetFinalStorage();
    }
    // ** End algorithm to be performed in case: input to landscape element from upstream elements

    // Fluxes accumulated
    thisElement->SetAccumulatedDischarge(accumulatedDischarge, thisElement->GetDischarge(), flowHierarchy, forceDirect);
    thisElement->SetAccumulatedLowerDischarge(accumulatedLowerDischarge, thisElement->GetLowerDischarge(), flowHierarchy, forceDirect);
    thisElement->SetAccumulatedUpperDischarge(accumulatedUpperDischarge, thisElement->GetUpperDischarge(), flowHierarchy, forceDirect);
    thisElement->SetAccumulatedPrecipitation(accumulatedPrecipitation + thisElement->GetPrecipitation()*thisElement->GetArea());
    thisElement->SetAccumulatedTemperature(accumulatedTemperature + thisElement->GetTemperature()*thisElement->GetArea());
    thisElement->SetAccumulatedEvapotranspiration(accumulatedEvapotranspiration + thisElement->GetEvapotranspiration()*thisElement->GetArea());
    thisElement->SetAccumulatedRunoff(accumulatedRunoff + thisElement->GetRunoff()*thisElement->GetArea());
    thisElement->SetAccumulatedUpperRunoff(accumulatedUpperRunoff + thisElement->GetUpperRunoff()*thisElement->GetArea());
    thisElement->SetAccumulatedLowerRunoff(accumulatedLowerRunoff + thisElement->GetLowerRunoff()*thisElement->GetArea());
    thisElement->SetAccumulatedSum(accumulatedSum + thisElement->GetArea());
    //  cout << thisElement->GetAccumulatedSum() << endl;
    // Lake state variables accumulated
    if (thisElement->GetLakeStorage() != missingData)
    {
        thisElement->SetAccumulatedLakeStorage(accumulatedLakeStorage + thisElement->GetLakeStorage()*thisElement->GetLakeArea());
        thisElement->SetAccumulatedSumLake(accumulatedSumLake + thisElement->GetLakeArea());
    }
    else
    {
        thisElement->SetAccumulatedLakeStorage(accumulatedLakeStorage);
        thisElement->SetAccumulatedSumLake(accumulatedSumLake);
    }
    // Snow state variables accumulated
    if (thisElement->GetSnowStore() != missingData)
    {
        thisElement->SetAccumulatedSnowStore(accumulatedSnowStore + thisElement->GetSnowStore()*thisElement->GetLandArea());
        thisElement->SetAccumulatedMeltWater(accumulatedMeltWater + thisElement->GetMeltWater()*thisElement->GetLandArea());
        thisElement->SetAccumulatedSnowWaterEquivalentChange(accumulatedSnowWaterEquivalentChange + thisElement->GetSnowWaterEquivalentChange()*thisElement->GetLandArea());
        thisElement->SetAccumulatedWaterOutput(accumulatedWaterOutput + thisElement->GetWaterOutput()*thisElement->GetLandArea());
        thisElement->SetAccumulatedSnowCoverFraction(accumulatedSnowCoverFraction + thisElement->GetSnowCoverFraction()*thisElement->GetLandArea());
        thisElement->SetAccumulatedSumSnow(accumulatedSumSnow + thisElement->GetLandArea());
    }
    else
    {
        thisElement->SetAccumulatedSnowStore(accumulatedSnowStore);
        thisElement->SetAccumulatedMeltWater(accumulatedMeltWater);
        thisElement->SetAccumulatedSnowWaterEquivalentChange(accumulatedSnowWaterEquivalentChange);
        thisElement->SetAccumulatedWaterOutput(accumulatedWaterOutput);
        thisElement->SetAccumulatedSnowCoverFraction(accumulatedSnowCoverFraction);
        thisElement->SetAccumulatedSumSnow(accumulatedSumSnow);
    }
    // Glacier state variables accumulated
    if (thisElement->GetGlacierMassBalance() != missingData)
    {
        //    cout << " 6 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
        thisElement->SetAccumulatedGlacierSnowStore(accumulatedGlacierSnowStore+thisElement->GetGlacierSnowStore()*thisElement->GetGlacierIceArea());
	thisElement->SetAccumulatedGlacierSnowMeltWater(accumulatedGlacierSnowMeltWater+thisElement->GetGlacierSnowMeltWater()*thisElement->GetGlacierIceArea());
        thisElement->SetAccumulatedGlacierMassBalance(accumulatedGlacierMassBalance + thisElement->GetGlacierMassBalance()*thisElement->GetGlacierIceArea());
        thisElement->SetAccumulatedAreaGlacierMassBalance(accumulatedAreaGlacierMassBalance + thisElement->GetAreaGlacierMassBalance()*thisElement->GetGlacierIceArea());
        thisElement->SetAccumulatedGlacierIceMelt(accumulatedGlacierIceMelt + thisElement->GetGlacierIceMelt()*thisElement->GetGlacierIceArea());
        thisElement->SetAccumulatedGlacierIceVolume(accumulatedGlacierIceVolume + thisElement->GetGlacierIceVolume());
        thisElement->SetAccumulatedSumGlacier(accumulatedSumGlacier + thisElement->GetGlacierIceArea());
        //    cout << " 6 " << accumulatedSumGlacier+thisElement->GetGlacierIceArea()
        //   << " " << accumulatedGlacierMassBalance+thisElement->GetGlacierMassBalance()*thisElement->GetGlacierIceArea() << endl;
    }
    else
    {
        //    cout << " 7 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
        thisElement->SetAccumulatedGlacierSnowStore(accumulatedGlacierSnowStore);
        thisElement->SetAccumulatedGlacierSnowMeltWater(accumulatedGlacierSnowMeltWater);
        thisElement->SetAccumulatedGlacierMassBalance(accumulatedGlacierMassBalance);
        thisElement->SetAccumulatedAreaGlacierMassBalance(accumulatedAreaGlacierMassBalance);
        thisElement->SetAccumulatedGlacierIceMelt(accumulatedGlacierIceMelt);
        thisElement->SetAccumulatedGlacierIceVolume(accumulatedGlacierIceVolume);
        thisElement->SetAccumulatedSumGlacier(accumulatedSumGlacier);
        //    cout << " 7 " << accumulatedSumGlacier << " " << accumulatedGlacierMassBalance << endl;
    }
    if (thisElement->GetAnnualMassBalance() != missingData)
    {
        thisElement->SetAccumulatedAnnualMassBalance(accumulatedAnnualMassBalance + thisElement->GetAnnualMassBalance()*thisElement->GetGlacierIceArea());
    }
    else
    {
        thisElement->SetAccumulatedAnnualMassBalance(accumulatedAnnualMassBalance);
    }
    // Hbv state variables accumulated
    if (thisElement->GetHbvSoilMoisture() != missingData)
    {
        thisElement->SetAccumulatedHbvSoilMoisture(accumulatedHbvSoilMoisture + thisElement->GetHbvSoilMoisture()*thisElement->GetHbvArea());
        thisElement->SetAccumulatedSumHbv(accumulatedSumHbv + thisElement->GetHbvArea());
    }
    else
    {
        thisElement->SetAccumulatedHbvSoilMoisture(accumulatedHbvSoilMoisture);
        thisElement->SetAccumulatedSumHbv(accumulatedSumHbv);
    }
    if (thisElement->GetHbvSoilMoistureDeficit() != missingData)
    {
        thisElement->SetAccumulatedHbvSoilMoistureDeficit(accumulatedHbvSoilMoistureDeficit + thisElement->GetHbvSoilMoistureDeficit()*thisElement->GetHbvArea());
    }
    else
    {
        thisElement->SetAccumulatedHbvSoilMoistureDeficit(accumulatedHbvSoilMoistureDeficit);
    }
    if (thisElement->GetHbvPercSoilUpper() != missingData)
    {
        thisElement->SetAccumulatedHbvPercSoilUpper(accumulatedHbvPercSoilUpper + thisElement->GetHbvPercSoilUpper()*thisElement->GetHbvArea());
    }
    else
    {
        thisElement->SetAccumulatedHbvPercSoilUpper(accumulatedHbvPercSoilUpper);
    }
    if (thisElement->GetHbvUpperZone() != missingData)
    {
        thisElement->SetAccumulatedHbvUpperZone(accumulatedHbvUpperZone + thisElement->GetHbvUpperZone()*thisElement->GetHbvArea());
    }
    else
    {
        thisElement->SetAccumulatedHbvUpperZone(accumulatedHbvUpperZone);
    }
    if (thisElement->GetHbvLowerZone() != missingData)
    {
        thisElement->SetAccumulatedHbvLowerZone(accumulatedHbvLowerZone + thisElement->GetHbvLowerZone()*thisElement->GetHbvArea());
    }
    else
    {
        thisElement->SetAccumulatedHbvLowerZone(accumulatedHbvLowerZone);
    }
    // Time series for landscape elements
    for (j = 0; j < thisElement->GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++)
    {
        if (thisElement->GetLandIndex() == thisElement->GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j))
        {
            thisElement->SetDistributedElementPrecipitation(timeStep, thisElement->GetPrecipitation());
            thisElement->SetDistributedElementTemperature(timeStep, thisElement->GetTemperature());
            thisElement->SetDistributedElementEvapotranspiration(timeStep, thisElement->GetEvapotranspiration());
            thisElement->SetDistributedElementRunoff(timeStep, thisElement->GetRunoff());
            thisElement->SetDistributedElementUpperRunoff(timeStep, thisElement->GetUpperRunoff());
            thisElement->SetDistributedElementLowerRunoff(timeStep, thisElement->GetLowerRunoff());
            thisElement->SetDistributedElementLakeStorage(timeStep, thisElement->GetLakeStorage());
            thisElement->SetDistributedElementSnowStore(timeStep, thisElement->GetSnowStore());
            thisElement->SetDistributedElementSnowCoverFraction(timeStep, thisElement->GetSnowCoverFraction());
            thisElement->SetDistributedElementMeltWater(timeStep, thisElement->GetMeltWater());
            thisElement->SetDistributedElementSnowWaterEquivalentChange(timeStep, thisElement->GetSnowWaterEquivalentChange());
            thisElement->SetDistributedElementWaterOutput(timeStep, thisElement->GetWaterOutput());
            thisElement->SetDistributedElementGlacierSnowStore(timeStep,thisElement->GetGlacierSnowStore());
            thisElement->SetDistributedElementGlacierSnowMeltWater(timeStep,thisElement->GetGlacierSnowMeltWater());
            thisElement->SetDistributedElementGlacierMassBalance(timeStep, thisElement->GetGlacierMassBalance());
            thisElement->SetDistributedElementAreaGlacierMassBalance(timeStep, thisElement->GetAreaGlacierMassBalance());
            thisElement->SetDistributedElementGlacierIceMelt(timeStep, thisElement->GetGlacierIceMelt());
            thisElement->SetDistributedElementAnnualMassBalance(timeStep, thisElement->GetAnnualMassBalance());
            thisElement->SetDistributedElementGlacierIceVolume(timeStep, thisElement->GetGlacierIceVolume());
            thisElement->SetDistributedElementHbvSoilMoisture(timeStep, thisElement->GetHbvSoilMoisture());
            thisElement->SetDistributedElementHbvSoilMoistureDeficit(timeStep, thisElement->GetHbvSoilMoistureDeficit());
            thisElement->SetDistributedElementHbvPercSoilUpper(timeStep, thisElement->GetHbvPercSoilUpper());
            thisElement->SetDistributedElementHbvUpperZone(timeStep, thisElement->GetHbvUpperZone());
            thisElement->SetDistributedElementHbvLowerZone(timeStep, thisElement->GetHbvLowerZone());
        }
    }
    for (j = 0; j < thisElement->GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++)
    {
        if (thisElement->GetLandIndex() == thisElement->GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j))
        {
            thisElement->SetDistributedElementPrecipitation(timeStep, thisElement->GetPrecipitation());
            thisElement->SetDistributedElementTemperature(timeStep, thisElement->GetTemperature());
            thisElement->SetDistributedElementEvapotranspiration(timeStep, thisElement->GetEvapotranspiration());
            thisElement->SetDistributedElementRunoff(timeStep, thisElement->GetRunoff());
            thisElement->SetDistributedElementUpperRunoff(timeStep, thisElement->GetUpperRunoff());
            thisElement->SetDistributedElementLowerRunoff(timeStep, thisElement->GetLowerRunoff());
            thisElement->SetDistributedElementLakeStorage(timeStep, thisElement->GetLakeStorage());
            thisElement->SetDistributedElementSnowStore(timeStep, thisElement->GetSnowStore());
            thisElement->SetDistributedElementSnowCoverFraction(timeStep, thisElement->GetSnowCoverFraction());
            thisElement->SetDistributedElementMeltWater(timeStep, thisElement->GetMeltWater());
            thisElement->SetDistributedElementSnowWaterEquivalentChange(timeStep, thisElement->GetSnowWaterEquivalentChange());
            thisElement->SetDistributedElementWaterOutput(timeStep, thisElement->GetWaterOutput());
            thisElement->SetDistributedElementGlacierSnowStore(timeStep,thisElement->GetGlacierSnowStore());
            thisElement->SetDistributedElementGlacierSnowMeltWater(timeStep,thisElement->GetGlacierSnowMeltWater());
            thisElement->SetDistributedElementGlacierMassBalance(timeStep, thisElement->GetGlacierMassBalance());
            thisElement->SetDistributedElementAreaGlacierMassBalance(timeStep, thisElement->GetAreaGlacierMassBalance());
            thisElement->SetDistributedElementGlacierIceMelt(timeStep, thisElement->GetGlacierIceMelt());
            thisElement->SetDistributedElementAnnualMassBalance(timeStep, thisElement->GetAnnualMassBalance());
            thisElement->SetDistributedElementGlacierIceVolume(timeStep, thisElement->GetGlacierIceVolume());
            thisElement->SetDistributedElementKiWaSoilMoistureOne(timeStep, thisElement->GetKiWaSoilMoisture
                    (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionOne(j)));
            thisElement->SetDistributedElementKiWaSoilMoistureTwo(timeStep, thisElement->GetKiWaSoilMoisture
                    (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionTwo(j)));
            thisElement->SetDistributedElementKiWaGroundWaterDepthOne(timeStep, thisElement->GetKiWaGroundWaterDepth
                    (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionOne(j)));
            thisElement->SetDistributedElementKiWaGroundWaterDepthTwo(timeStep, thisElement->GetKiWaGroundWaterDepth
                    (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionTwo(j)));
        }
    }
}


void TraverseMissingDataSubCatchment(SubCatchment * const thisSubCatchment, int timeStep, ofstream &fout)
{
    int i;
    DistributedElement * thisElement;

    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        TraverseMissingDataSubCatchment(thisSubCatchment->GetUpStream(i), timeStep, fout);
    }
    thisElement = thisSubCatchment->GetLandScapeElement();
    while (thisElement)
    {
        TraverseMissingDataLandScape(thisElement, timeStep, fout);
        thisElement = thisElement->GetNextElement();
    }
    thisSubCatchment->SetAccumulatedDischarge(timeStep, missingData);
    thisSubCatchment->SetAccumulatedUpperDischarge(timeStep, missingData);
    thisSubCatchment->SetAccumulatedLowerDischarge(timeStep, missingData);
    thisSubCatchment->SetAccumulatedInFlow(timeStep, missingData);
    //  thisSubCatchment->SetAccumulatedUpperInFlow(timeStep, missingData);
    //  thisSubCatchment->SetAccumulatedLowerInFlow(timeStep, missingData);
    thisSubCatchment->SetAccumulatedPrecipitation(timeStep, missingData);
    thisSubCatchment->SetAccumulatedTemperature(timeStep, missingData);
    thisSubCatchment->SetAccumulatedLakeStorage(timeStep, missingData);
    thisSubCatchment->SetAccumulatedSnowStore(timeStep, missingData);
    thisSubCatchment->SetAccumulatedMeltWater(timeStep, missingData);
    thisSubCatchment->SetAccumulatedSnowWaterEquivalentChange(timeStep, missingData);
    thisSubCatchment->SetAccumulatedWaterOutput(timeStep, missingData);
    thisSubCatchment->SetAccumulatedSnowCoverFraction(timeStep, missingData);
    thisSubCatchment->SetAccumulatedGlacierSnowStore(timeStep, missingData);
    thisSubCatchment->SetAccumulatedGlacierSnowMeltWater(timeStep, missingData);
    thisSubCatchment->SetAccumulatedGlacierMassBalance(timeStep, missingData);
    thisSubCatchment->SetAccumulatedAreaGlacierMassBalance(timeStep, missingData);
    thisSubCatchment->SetAccumulatedGlacierIceMelt(timeStep, missingData);
    thisSubCatchment->SetAccumulatedAnnualMassBalance(timeStep, missingData);
    thisSubCatchment->SetAccumulatedGlacierIceVolume(timeStep, missingData);
    thisSubCatchment->SetAccumulatedEvapotranspiration(timeStep, missingData);
    thisSubCatchment->SetAccumulatedRunoff(timeStep, missingData);
    thisSubCatchment->SetAccumulatedUpperRunoff(timeStep, missingData);
    thisSubCatchment->SetAccumulatedLowerRunoff(timeStep, missingData);
    thisSubCatchment->SetAccumulatedHbvSoilMoisture(timeStep, missingData);
    thisSubCatchment->SetAccumulatedHbvSoilMoistureDeficit(timeStep, missingData);
    thisSubCatchment->SetAccumulatedHbvPercSoilUpper(timeStep, missingData);
    thisSubCatchment->SetAccumulatedHbvUpperZone(timeStep, missingData);
    thisSubCatchment->SetAccumulatedHbvLowerZone(timeStep, missingData);
    thisSubCatchment->SetAccumulatedSum(timeStep, missingData);
    thisSubCatchment->SetAccumulatedSumSnow(timeStep, missingData);
    thisSubCatchment->SetAccumulatedSumGlacier(timeStep, missingData);
    thisSubCatchment->SetAccumulatedSumHbv(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentPrecipitation(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentTemperature(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentLakeStorage(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentSnowStore(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentMeltWater(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentSnowWaterEquivalentChange(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentGlacierSnowStore(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentGlacierSnowMeltWater(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentAreaGlacier(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentAreaGlacierMassBalance(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentAnnualMassBalance(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentGlacierIceVolume(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentEvapotranspiration(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentRunoff(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentUpperRunoff(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentLowerRunoff(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, missingData);
}


void TraverseMissingDataLandScape(DistributedElement * const thisElement, int timeStep, ofstream &fout)
{
    int i, j;
    for (i = 0; i < thisElement->GetNumUpLand(); i++)
    {
        TraverseMissingDataLandScape(thisElement->GetUpLandFlow(i), timeStep, fout);
    }
    thisElement->SetAccumulatedDischarge(0.0, missingData, false, false);
    thisElement->SetAccumulatedLowerDischarge(0.0, missingData, false, false);
    thisElement->SetAccumulatedUpperDischarge(0.0, missingData, false, false);
    thisElement->SetAccumulatedPrecipitation(missingData);
    thisElement->SetAccumulatedTemperature(missingData);
    thisElement->SetAccumulatedLakeStorage(missingData);
    thisElement->SetAccumulatedSnowStore(missingData);
    thisElement->SetAccumulatedMeltWater(missingData);
    thisElement->SetAccumulatedSnowWaterEquivalentChange(missingData);
    thisElement->SetAccumulatedWaterOutput(missingData);
    thisElement->SetAccumulatedSnowCoverFraction(missingData);
    thisElement->SetAccumulatedGlacierSnowStore(missingData);
    thisElement->SetAccumulatedGlacierSnowMeltWater(missingData);
    thisElement->SetAccumulatedGlacierMassBalance(missingData);
    thisElement->SetAccumulatedAreaGlacierMassBalance(missingData);
    thisElement->SetAccumulatedGlacierIceMelt(missingData);
    thisElement->SetAccumulatedAnnualMassBalance(missingData);
    thisElement->SetAccumulatedGlacierIceVolume(missingData);
    thisElement->SetAccumulatedEvapotranspiration(missingData);
    thisElement->SetAccumulatedRunoff(missingData);
    thisElement->SetAccumulatedUpperRunoff(missingData);
    thisElement->SetAccumulatedLowerRunoff(missingData);
    thisElement->SetAccumulatedHbvSoilMoisture(missingData);
    thisElement->SetAccumulatedHbvSoilMoistureDeficit(missingData);
    thisElement->SetAccumulatedHbvPercSoilUpper(missingData);
    thisElement->SetAccumulatedHbvUpperZone(missingData);
    thisElement->SetAccumulatedHbvLowerZone(missingData);
    thisElement->SetAccumulatedSum(missingData);
    thisElement->SetAccumulatedSumSnow(missingData);
    thisElement->SetAccumulatedSumGlacier(missingData);
    thisElement->SetAccumulatedSumHbv(missingData);
    for (j = 0; j < thisElement->GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++)
    {
        if (thisElement->GetLandIndex() == thisElement->GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j))
        {
            thisElement->SetDistributedElementPrecipitation(timeStep, missingData);
            thisElement->SetDistributedElementTemperature(timeStep, missingData);
            thisElement->SetDistributedElementLakeStorage(timeStep, missingData);
            thisElement->SetDistributedElementSnowStore(timeStep, missingData);
            thisElement->SetDistributedElementSnowCoverFraction(timeStep, missingData);
            thisElement->SetDistributedElementMeltWater(timeStep, missingData);
            thisElement->SetDistributedElementSnowWaterEquivalentChange(timeStep, missingData);
            thisElement->SetDistributedElementWaterOutput(timeStep, missingData);
            thisElement->SetDistributedElementGlacierSnowStore(timeStep,missingData);
            thisElement->SetDistributedElementGlacierSnowMeltWater(timeStep,missingData);
            thisElement->SetDistributedElementGlacierMassBalance(timeStep, missingData);
            thisElement->SetDistributedElementAreaGlacierMassBalance(timeStep, missingData);
            thisElement->SetDistributedElementGlacierIceMelt(timeStep, missingData);
            thisElement->SetDistributedElementAnnualMassBalance(timeStep, missingData);
            thisElement->SetDistributedElementGlacierIceVolume(timeStep, missingData);
            thisElement->SetDistributedElementEvapotranspiration(timeStep, missingData);
            thisElement->SetDistributedElementRunoff(timeStep, missingData);
            thisElement->SetDistributedElementUpperRunoff(timeStep, missingData);
            thisElement->SetDistributedElementLowerRunoff(timeStep, missingData);
            thisElement->SetDistributedElementHbvSoilMoisture(timeStep, missingData);
            thisElement->SetDistributedElementHbvSoilMoistureDeficit(timeStep, missingData);
            thisElement->SetDistributedElementHbvPercSoilUpper(timeStep, missingData);
            thisElement->SetDistributedElementHbvUpperZone(timeStep, missingData);
            thisElement->SetDistributedElementHbvLowerZone(timeStep, missingData);
            //      cout << "HBV  " << timeStep << "  " << missingData << "  " << thisElement->GetDistributedElementHbvUpperZone(timeStep) << endl;
        }
    }
    for (j = 0; j < thisElement->GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++)
    {
        if (thisElement->GetLandIndex() == thisElement->GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j))
        {
            thisElement->SetDistributedElementPrecipitation(timeStep, missingData);
            thisElement->SetDistributedElementTemperature(timeStep, missingData);
            thisElement->SetDistributedElementLakeStorage(timeStep, missingData);
            thisElement->SetDistributedElementSnowStore(timeStep, missingData);
            thisElement->SetDistributedElementSnowCoverFraction(timeStep, missingData);
            thisElement->SetDistributedElementMeltWater(timeStep, missingData);
            thisElement->SetDistributedElementSnowWaterEquivalentChange(timeStep, missingData);
            thisElement->SetDistributedElementWaterOutput(timeStep, missingData);
            thisElement->SetDistributedElementGlacierSnowStore(timeStep,missingData);
            thisElement->SetDistributedElementGlacierSnowMeltWater(timeStep,missingData);
            thisElement->SetDistributedElementGlacierMassBalance(timeStep, missingData);
            thisElement->SetDistributedElementAreaGlacierMassBalance(timeStep, missingData);
            thisElement->SetDistributedElementGlacierIceMelt(timeStep, missingData);
            thisElement->SetDistributedElementAnnualMassBalance(timeStep, missingData);
            thisElement->SetDistributedElementGlacierIceVolume(timeStep, missingData);
            thisElement->SetDistributedElementEvapotranspiration(timeStep, missingData);
            thisElement->SetDistributedElementRunoff(timeStep, missingData);
            thisElement->SetDistributedElementUpperRunoff(timeStep, missingData);
            thisElement->SetDistributedElementLowerRunoff(timeStep, missingData);
            thisElement->SetDistributedElementKiWaSoilMoistureOne(timeStep, missingData);
            thisElement->SetDistributedElementKiWaSoilMoistureTwo(timeStep, missingData);
            thisElement->SetDistributedElementKiWaGroundWaterDepthOne(timeStep, missingData);
            thisElement->SetDistributedElementKiWaGroundWaterDepthTwo(timeStep, missingData);
            //      cout << "KiWa " << timeStep << "  " << missingData << "  " << thisElement->GetDistributedElementKiWaGroundWaterDepthOne(timeStep) << endl;
        }
    }
}


// For test purpose only
void WriteSubCatchmentIdentifier(SubCatchment * const CatchmentElement, int numWatc, ofstream &fout)
{
    int i;
    DistributedElement * thisElement;

    // Sub-Catchment outlets
    ofstream subCatchmentOut("test_waterland.txt");  // Open for writing"
    for (i = 0; i < numWatc; i++)
    {
        if (CatchmentElement[i].GetLandScapeElement())
        {
            subCatchmentOut << "#  " << CatchmentElement[i].GetIdentifier() << "  #  "
                            << CatchmentElement[i].GetNumLandScape() << endl;
            thisElement = CatchmentElement[i].GetLandScapeElement();
            while (thisElement)
            {
                subCatchmentOut << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
                thisElement = thisElement->GetNextElement();
            }
        }
    }
    subCatchmentOut << endl;
    subCatchmentOut.close();
}


void WriteSubCatchmentDischarge(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime,
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                int secondsPerTimeStep, bool modelCalibration, ofstream &fout)
{
    FILE *fpOut, *fpInf;
    char fileName[100];
    double sumValue;
    double nashSutcliffe;        /*  Nash-Sutcliffe efficiency criterion  */
    double rootMSError;          /*  Root mean square error criterion  */
    double biasVol;              /*  Bias (volume error) criterion  */
    double pearsCorr;            /*  Pearson's product-moment correlation coefficient  */
    int i, j, k, timeStep;
    DateTime datetime;
    double * observed = new double[numberTimeSteps];
    double * simulated = new double[numberTimeSteps];

    for (i = 0; i < numWatc; i++)
    {
        //    cout << CatchmentElement[i].GetSelectedSubCatchmentTimeSeriesElements()->GetNumberElements() << endl;
        for (k = 0; k < CatchmentElement[i].GetSelectedSubCatchmentTimeSeriesElements()->GetNumberElements(); k++)
        {
            //      cout << CatchmentElement[i].GetIdentifier() << "    " <<  CatchmentElement[i].GetSelectedSubCatchmentTimeSeriesElements()->GetSubCatchmentTimeSeriesElement(k) << endl;
            if (CatchmentElement[i].GetIdentifier() == CatchmentElement[i].GetSelectedSubCatchmentTimeSeriesElements()->GetSubCatchmentTimeSeriesElement(k))
            {
                sprintf(fileName, "dew_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpOut = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "inf_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpInf = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sumValue = 0.0;
                j = 0;
                timeStep = initialTimeSteps;
                for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
                {
                    if (CatchmentElement[i].GetAccumulatedDischarge(timeStep) > missingData)
                    {
                        fprintf(fpOut, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetAccumulatedDischarge(timeStep)*CatchmentElement[i].GetCorrection());
                        if (modelCalibration && CatchmentElement[i].GetObsData(j) > missingData)
                        {
                            sumValue = sumValue + CatchmentElement[i].GetAccumulatedDischarge(timeStep);
                        }
                    }
                    else
                    {
                        fprintf(fpOut, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    if (CatchmentElement[i].GetAccumulatedDischarge(timeStep) > missingData)
                    {
                        fprintf(fpInf, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetAccumulatedInFlow(timeStep));
                    }
                    else
                    {
                        fprintf(fpInf, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    observed[j] = CatchmentElement[i].GetObsData(j);
                    simulated[j] = CatchmentElement[i].GetAccumulatedDischarge(timeStep) * CatchmentElement[i].GetCorrection();
                    j++;
                    timeStep++;
                }
                if (modelCalibration)
                {
                    fprintf(fpOut, "%29.6f\n", sumValue * CatchmentElement[i].GetCorrection());
                }
                fclose(fpOut);
                fclose(fpInf);
                if (timeStep != initialTimeSteps + numberTimeSteps)
                {
                    cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
                    exit(1);
                }
                ObjectiveCriteria(numberTimeSteps, observed, simulated, &nashSutcliffe, &rootMSError, &biasVol, &pearsCorr);
                cout.precision(2);
                cout.setf(ios::fixed);
                cout.setf(ios::showpoint);
                cout << "\n    Sub-catchment :  " << CatchmentElement[i].GetIdentifier() << endl << endl;
                cout.width(10);
                cout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;;
                cout.width(10);
                cout << "    Root mean square error                              " << rootMSError << endl;
                cout.width(10);
                cout << "    Bias (volume error)                                 " << biasVol << endl;
                cout.width(10);
                cout << "    Volume error correction factor                      " << 1/(1+biasVol) << endl;
                cout.width(10);
                cout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
                cout << endl;
                fout.precision(2);
                fout.setf(ios::fixed);
                fout.setf(ios::showpoint);
                fout << "\n    Sub-catchment :  " << CatchmentElement[i].GetIdentifier() << endl << endl;
                fout.width(10);
                fout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;;
                fout.width(10);
                fout << "    Root mean square error                              " << rootMSError << endl;
                fout.width(10);
                fout << "    Bias (volume error)                                 " << biasVol << endl;
                fout.width(10);
                fout << "    Volume error correction factor                      " << 1/(1+biasVol) << endl;
                fout.width(10);
                fout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
                fout << endl;
            }
        }
    }

    delete[] observed;
    delete[] simulated;
}


void WriteSubCatchmentWaterBalance(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime,
                                   DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                   int secondsPerTimeStep)
{
    FILE *fpPre, *fpTem, *fpLak, *fpSwe, *fpSch, *fpIns, *fpGar, *fpGmb, *fpGmc, *fpGsw, *fpGim, *fpGiv, *fpGam, *fpEva,* fpRun,* fpRup,* fpRlo, *fpHsd, *fpHsm, *fpHpe, *fpHuz, *fpHlz, *fpHgw;
    char fileName[100];
    int i, k, timeStep;
    DateTime datetime;
    for (i = 0; i < numWatc; i++)
    {
        for (k = 0; k < CatchmentElement[i].GetSelectedSubCatchmentTimeSeriesElements()->GetNumberElements(); k++)
        {
            if (CatchmentElement[i].GetIdentifier() == CatchmentElement[i].GetSelectedSubCatchmentTimeSeriesElements()->GetSubCatchmentTimeSeriesElement(k))
            {
                sprintf(fileName, "pre_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpPre = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "tem_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpTem = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "lak_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpLak = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "swe_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpSwe = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "sch_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpSch = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "ins_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpIns = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "gar_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpGar = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "gmb_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpGmb = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "gmc_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpGmc = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
		/*		sprintf(fileName,"gsw_%08d.var",CatchmentElement[i].GetIdentifier());
		if ((fpGsw = fopen(fileName, "w")) == NULL ) 
		{
		    printf("\n Filen %s ikke funnet!\n\n",fileName);
		    exit(1);
		    }*/
                sprintf(fileName, "gim_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpGim = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "giv_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpGiv = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "gam_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpGam = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
		}
                sprintf(fileName, "eva_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpEva = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "run_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpRun = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "rup_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpRup = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "rlo_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpRlo = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "hsd_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpHsd = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "hsm_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpHsm = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "hpe_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpHpe = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "huz_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpHuz = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "hlz_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpHlz = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                sprintf(fileName, "hgw_%08d.var", CatchmentElement[i].GetIdentifier());
                if ((fpHgw = fopen(fileName, "w")) == NULL)
                {
                    printf("\n Filen %s ikke funnet!\n\n", fileName);
                    exit(1);
                }
                timeStep = initialTimeSteps;
                for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
                {
                    // Fluxes
                    if (CatchmentElement[i].GetSubCatchmentPrecipitation(timeStep) != missingData)
                    {
                        fprintf(fpPre, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentPrecipitation(timeStep) * 1000.0);
                        fprintf(fpTem, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentTemperature(timeStep));
                        fprintf(fpEva, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentEvapotranspiration(timeStep) * 1000.0);
                        fprintf(fpRun, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentRunoff(timeStep) * 1000.0);
                        fprintf(fpRup, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentUpperRunoff(timeStep) * 1000.0);
                        //                  CatchmentElement[i].GetSubCatchmentUpperRunoff(timeStep)*1000.0+CatchmentElement[i].GetSubCatchmentLowerRunoff(timeStep)*1000.0);
                        fprintf(fpRlo, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentLowerRunoff(timeStep) * 1000.0);
                    }
                    else
                    {
                        fprintf(fpPre, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpTem, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpEva, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpRun, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpRup, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpRlo, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    // Lake state variables
                    if (CatchmentElement[i].GetSubCatchmentLakeStorage(timeStep) != missingData)
                    {
                        fprintf(fpLak, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentLakeStorage(timeStep) * 1000.0);
                    }
                    else
                    {
                        fprintf(fpLak, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    // Snow state variables
                    if (CatchmentElement[i].GetSubCatchmentSnowStore(timeStep) != missingData)
                    {
                        fprintf(fpSwe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                (CatchmentElement[i].GetSubCatchmentSnowStore(timeStep) + CatchmentElement[i].GetSubCatchmentMeltWater(timeStep)) * 1000.0);
                        fprintf(fpSch, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentSnowWaterEquivalentChange(timeStep) * 1000.0);
                        fprintf(fpIns, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentWaterOutput(timeStep) * 1000.0);
                    }
                    else
                    {
                        fprintf(fpSwe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpSch, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpIns, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    // Glacier area, mass balance, ice melt, snow store and ice volume
                    if (CatchmentElement[i].GetSubCatchmentGlacierMassBalance(timeStep) != missingData)
                    {
                        fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentGlacierMassBalance(timeStep) * 1000.0);
                        //                    CatchmentElement[i].GetSubCatchmentGlacierMassBalance(timeStep));
                        fprintf(fpGim, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentGlacierIceMelt(timeStep) * 1000.0);
                        //                    CatchmentElement[i].GetSubCatchmentGlacierIceMelt(timeStep));
                        fprintf(fpGmc, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentAreaGlacierMassBalance(timeStep) * 1000.0);
			/*			fprintf(fpGsw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
				datetime.getDay(), datetime.getHour(), datetime.getMinute(),
				(CatchmentElement[i].GetSubCatchmentGlacierSnowStore(timeStep)+CatchmentElement[i].GetSubCatchmentGlacierSnowMeltWater(timeStep))*1000.0);*/
                        fprintf(fpGar, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentAreaGlacier(timeStep) / 1.0e6);
                    }
                    else
                    {
                        fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpGim, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpGmc, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
			/*			fprintf(fpGsw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
						datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);*/
                        fprintf(fpGar, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    if (CatchmentElement[i].GetSubCatchmentGlacierIceVolume(timeStep) != missingData)
                    {
                        fprintf(fpGiv, "%4d%02d%02d/%02d%02d%36.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentGlacierIceVolume(timeStep));
                    }
                    else
                    {
                        fprintf(fpGiv, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    if (CatchmentElement[i].GetSubCatchmentAnnualMassBalance(timeStep) != missingData)
                    {
                        fprintf(fpGam, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentAnnualMassBalance(timeStep) * 1000.0);
                        //                    CatchmentElement[i].GetSubCatchmentAnnualMassBalance(timeStep));
                    }
                    else
                    {
                        fprintf(fpGam, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
		    }
                    // Hbv state variables
                    if (CatchmentElement[i].GetSubCatchmentHbvSoilMoisture(timeStep) != missingData)
                    {
                        fprintf(fpHsd, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentHbvSoilMoistureDeficit(timeStep) * 1000.0);
                        fprintf(fpHsm, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentHbvSoilMoisture(timeStep) * 1000.0);
                        fprintf(fpHpe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentHbvPercSoilUpper(timeStep) * 1000.0);
                        fprintf(fpHuz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentHbvUpperZone(timeStep) * 1000.0);
                        fprintf(fpHlz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                CatchmentElement[i].GetSubCatchmentHbvLowerZone(timeStep) * 1000.0);
                        fprintf(fpHgw, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                (CatchmentElement[i].GetSubCatchmentHbvUpperZone(timeStep) + CatchmentElement[i].GetSubCatchmentHbvLowerZone(timeStep)) * 1000.0);
                    }
                    else
                    {
                        fprintf(fpHsd, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpHsm, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpHpe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpHuz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpHlz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                        fprintf(fpHgw, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    timeStep++;
                }
                fclose(fpPre);
                fclose(fpTem);
                fclose(fpLak);
                fclose(fpSwe);
                fclose(fpSch);
                fclose(fpIns);
                fclose(fpGar);
                fclose(fpGmb);
                fclose(fpGmc);
		//		fclose(fpGsw);
                fclose(fpGim);
                fclose(fpGiv);
                fclose(fpGam);
                fclose(fpEva);
                fclose(fpRun);
                fclose(fpRup);
                fclose(fpRlo);
                fclose(fpHsd);
                fclose(fpHsm);
                fclose(fpHpe);
                fclose(fpHuz);
                fclose(fpHlz);
                if (timeStep != initialTimeSteps + numberTimeSteps)
                {
                    cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
                    exit(1);
                }
            }
        }
    }
}


void WriteOutletDischarge(SubCatchment * const thisSubCatchment, DateTime startSimulationTime,
                          DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                          int secondsPerTimeStep, bool modelCalibration, ofstream &fout)
{
    FILE *fpOut;
    char fileName[100];
    double sumValue;
    double nashSutcliffe;        /*  Nash-Sutcliffe efficiency criterion  */
    double rootMSError;          /*  Root mean square error criterion  */
    double biasVol;              /*  Bias (volume error) criterion  */
    double pearsCorr;            /*  Pearson's product-moment correlation coefficient  */
    int i, j, k, timeStep;
    DateTime datetime;
    double * observed = new double[numberTimeSteps];
    double * simulated = new double[numberTimeSteps];

    sprintf(fileName, "dew_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpOut = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sumValue = 0.0;
    j = 0;
    timeStep = initialTimeSteps;
    for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
    {
        if (thisSubCatchment->GetAccumulatedDischarge(timeStep) > missingData)
        {
            fprintf(fpOut, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedDischarge(timeStep)*thisSubCatchment->GetCorrection());
            if (modelCalibration && thisSubCatchment->GetObsData(j) > missingData)
            {
                sumValue = sumValue + thisSubCatchment->GetAccumulatedDischarge(timeStep);
            }
        }
        else
        {
            fprintf(fpOut, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
        }
        observed[j] = thisSubCatchment->GetObsData(j);
        simulated[j] = thisSubCatchment->GetAccumulatedDischarge(timeStep) * thisSubCatchment->GetCorrection();
        j++;
        timeStep++;
    }
    if (modelCalibration)
    {
        fprintf(fpOut, "%29.6f\n", sumValue * thisSubCatchment->GetCorrection());
    }
    fclose(fpOut);
    if (timeStep != initialTimeSteps + numberTimeSteps)
    {
        cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
        exit(1);
    }
    ObjectiveCriteria(numberTimeSteps, observed, simulated, &nashSutcliffe, &rootMSError, &biasVol, &pearsCorr);
    cout.precision(2);
    cout.setf(ios::fixed);
    cout.setf(ios::showpoint);
    cout << "\n    Watercourse outlet :  " << thisSubCatchment->GetIdentifier() << endl << endl;
    cout.width(10);
    cout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;;
    cout.width(10);
    cout << "    Root mean square error                              " << rootMSError << endl;
    cout.width(10);
    cout << "    Bias (volume error)                                 " << biasVol << endl;
    cout.width(10);
    cout << "    Volume error correction factor                      " << 1/(1+biasVol) << endl;
    cout.width(10);
    cout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
    cout << endl;
    fout.precision(2);
    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout << "\n    Watercourse outlet :  " << thisSubCatchment->GetIdentifier() << endl << endl;
    fout.width(10);
    fout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;;
    fout.width(10);
    fout << "    Root mean square error                              " << rootMSError << endl;
    fout.width(10);
    fout << "    Bias (volume error)                                 " << biasVol << endl;
    fout.width(10);
    fout << "    Volume error correction factor                      " << 1/(1+biasVol) << endl;
    fout.width(10);
    fout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
    fout << endl;

    delete[] observed;
    delete[] simulated;
}


void WriteOutletWaterBalance(SubCatchment * const thisSubCatchment, DateTime startSimulationTime,
                             DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                             int secondsPerTimeStep, bool flowHierarchy, bool forceDirect)
{
    FILE *fpPre, *fpTem, *fpLak, *fpSwe, *fpSch, *fpIns, *fpGar, *fpGmb, *fpGmc, *fpGsw, *fpGim, *fpGiv, *fpGam, *fpEva, *fpRun, *fpRup, *fpRlo, *fpHsd, *fpHsm, *fpHpe, *fpHuz, *fpHlz, *fpHgw;
    char fileName[100];
    int i, k, timeStep;
    DateTime datetime;

    sprintf(fileName, "pre_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpPre = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "tem_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpTem = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    /*    sprintf(fileName, "lak_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpLak = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
	}*/
    sprintf(fileName, "swe_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpSwe = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "sch_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpSch = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "ins_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpIns = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "gar_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpGar = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "gmb_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpGmb = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "gmc_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpGmc = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    /*    sprintf(fileName,"gsw_%08d.out",thisSubCatchment->GetIdentifier());
    if ((fpGsw = fopen(fileName, "w")) == NULL ) 
    {
        printf("\n Filen %s ikke funnet!\n\n",fileName);
        exit(1);
	}*/
    /*    sprintf(fileName, "gim_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpGim = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
	}*/
    sprintf(fileName, "giv_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpGiv = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "gam_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpGam = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "eva_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpEva = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "run_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpRun = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "rup_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpRup = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "rlo_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpRlo = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "hsd_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpHsd = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "hsm_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpHsm = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "hpe_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpHpe = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "huz_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpHuz = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "hlz_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpHlz = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    sprintf(fileName, "hgw_%08d.out", thisSubCatchment->GetIdentifier());
    if ((fpHgw = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }

    timeStep = initialTimeSteps;
    for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
    {
        // Fluxes
        if (thisSubCatchment->GetAccumulatedPrecipitation(timeStep) != missingData && thisSubCatchment->GetAccumulatedSum(timeStep) > 0.0)
        {
            fprintf(fpPre, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedPrecipitation(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            fprintf(fpTem, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedTemperature(timeStep) / thisSubCatchment->GetAccumulatedSum(timeStep));
            fprintf(fpEva, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedEvapotranspiration(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            if (flowHierarchy && !forceDirect)
            {
                fprintf(fpRun, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                        datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                        thisSubCatchment->GetAccumulatedDischarge(timeStep)*secondsPerTimeStep * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
                fprintf(fpRup, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                        datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                        thisSubCatchment->GetAccumulatedUpperDischarge(timeStep)*secondsPerTimeStep * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
                fprintf(fpRlo, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                        datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                        thisSubCatchment->GetAccumulatedLowerDischarge(timeStep)*secondsPerTimeStep * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            }
            else
            {
                fprintf(fpRun, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                        datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                        thisSubCatchment->GetAccumulatedRunoff(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
                fprintf(fpRup, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                        datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                        thisSubCatchment->GetAccumulatedUpperRunoff(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
                fprintf(fpRlo, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                        datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                        thisSubCatchment->GetAccumulatedLowerRunoff(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            }
        }
        else
        {
            fprintf(fpPre, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpTem, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpEva, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpRun, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpRup, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpRlo, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
        }
        // Lake state variables
	/*        if (thisSubCatchment->GetAccumulatedLakeStorage(timeStep) != missingData && thisSubCatchment->GetAccumulatedSumLake(timeStep) > 0.0)
        {
            fprintf(fpLak, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedLakeStorage(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedLakeStorage(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumLake(timeStep));
        }
        else
        {
            fprintf(fpLak, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
		    }*/
        // Snow state variables
        if (thisSubCatchment->GetAccumulatedSnowStore(timeStep) != missingData && thisSubCatchment->GetAccumulatedSumSnow(timeStep) > 0.0)
        {
            fprintf(fpSwe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    (thisSubCatchment->GetAccumulatedSnowStore(timeStep) + thisSubCatchment->GetAccumulatedMeltWater(timeStep)) * 1000.0
                    / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              /thisSubCatchment->GetAccumulatedSumSnow(timeStep));
            fprintf(fpSch, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedSnowWaterEquivalentChange(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedSnowWaterEquivalentChange(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumSnow(timeStep));
            fprintf(fpIns, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedWaterOutput(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedWaterOutput(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumSnow(timeStep));
        }
        else
        {
            fprintf(fpSwe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpSch, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpIns, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
        }
	// Glacier area, mass balance and ice volume
        if (thisSubCatchment->GetAccumulatedGlacierMassBalance(timeStep) != missingData && thisSubCatchment->GetAccumulatedSumGlacier(timeStep) > 0.0)
        {
            fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedGlacierMassBalance(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSumGlacier(timeStep));
            //              thisSubCatchment->GetAccumulatedGlacierMassBalance(timeStep)/thisSubCatchment->GetAccumulatedSumGlacier(timeStep));
	    /*      fprintf(fpGim, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedGlacierIceMelt(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
		    //              thisSubCatchment->GetAccumulatedGlacierIceMelt(timeStep)/thisSubCatchment->GetAccumulatedSumGlacier(timeStep));*/
            fprintf(fpGmc, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedAreaGlacierMassBalance(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
	    /*	    fprintf(fpGsw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
		    (thisSubCatchment->GetAccumulatedGlacierSnowStore(timeStep)+thisSubCatchment->GetAccumulatedGlacierSnowMeltWater(timeStep))*1000.0
		    /thisSubCatchment->GetAccumulatedSum(timeStep));*/
            fprintf(fpGar, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedSumGlacier(timeStep) / 1.0e6);
        }
        else
        {
            fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            /*	    fprintf(fpGim, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);*/
            fprintf(fpGmc, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
	    /*	    fprintf(fpGsw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);*/
            fprintf(fpGar, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
        }
        if (thisSubCatchment->GetAccumulatedGlacierIceVolume(timeStep) != missingData)
        {
            fprintf(fpGiv, "%4d%02d%02d/%02d%02d%36.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedGlacierIceVolume(timeStep));
        }
        else
        {
            fprintf(fpGiv, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
        }
        if (thisSubCatchment->GetAccumulatedAnnualMassBalance(timeStep) != missingData && thisSubCatchment->GetAccumulatedSumGlacier(timeStep) > 0.0)
        {
            fprintf(fpGam, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedAnnualMassBalance(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSumGlacier(timeStep));
            //              thisSubCatchment->GetAccumulatedAnnualMassBalance(timeStep)/thisSubCatchment->GetAccumulatedSumGlacier(timeStep));
        }
        else
        {
            fprintf(fpGam, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
	}
        // Hbv state variables
        if (thisSubCatchment->GetAccumulatedHbvSoilMoisture(timeStep) != missingData && thisSubCatchment->GetAccumulatedSumHbv(timeStep) > 0.0)
        {
            fprintf(fpHsd, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedHbvSoilMoistureDeficit(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedHbvSoilMoistureDeficit(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumHbv(timeStep));
            fprintf(fpHsm, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedHbvSoilMoisture(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedHbvSoilMoisture(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumHbv(timeStep));
            fprintf(fpHpe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedHbvPercSoilUpper(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedHbvPercSoilUpper(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumHbv(timeStep));
            fprintf(fpHuz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedHbvUpperZone(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedHbvUpperZone(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumHbv(timeStep));
            fprintf(fpHlz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    thisSubCatchment->GetAccumulatedHbvLowerZone(timeStep) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
            //              thisSubCatchment->GetAccumulatedHbvLowerZone(timeStep)*1000.0/thisSubCatchment->GetAccumulatedSumHbv(timeStep));
            fprintf(fpHgw, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    //              (thisSubCatchment->GetAccumulatedHbvUpperZone(timeStep)+thisSubCatchment->GetAccumulatedHbvLowerZone(timeStep))*1000.0/thisSubCatchment->GetAccumulatedSumHbv(timeStep));
                    (thisSubCatchment->GetAccumulatedHbvUpperZone(timeStep) + thisSubCatchment->GetAccumulatedHbvLowerZone(timeStep)) * 1000.0 / thisSubCatchment->GetAccumulatedSum(timeStep));
        }
        else
        {
            fprintf(fpHsd, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpHsm, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpHpe, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpHuz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpHlz, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
            fprintf(fpHgw, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
        }
        timeStep++;
    }
    fclose(fpPre);
    fclose(fpTem);
    //    fclose(fpLak);
    fclose(fpSwe);
    fclose(fpSch);
    fclose(fpIns);
    fclose(fpGar);
    fclose(fpGmb);
    fclose(fpGmc);
    //    fclose(fpGsw);
    //    fclose(fpGim);
    fclose(fpGiv);
    fclose(fpGam);
    fclose(fpEva);
    fclose(fpRun);
    fclose(fpRup);
    fclose(fpRlo);
    fclose(fpHsd);
    fclose(fpHsm);
    fclose(fpHpe);
    fclose(fpHuz);
    fclose(fpHlz);
    if (timeStep != initialTimeSteps + numberTimeSteps)
    {
        cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
        exit(1);
    }
}


void WriteReducedBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,
                            int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
    int i, j, k;
    char fileName[100];
    char dirName[100];
    char timeName[100];
    //  unsigned short int noData16=65535;
    ofstream fileTemp, filePre, fileEva, fileSwe, fileScf, fileSmw, fileRun, fileDew, fileHsd, fileHsm, fileHgw, fileHlz;
    sprintf(dirName, "./%04d/", datetime.getYear());
    sprintf(timeName, "_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());

    //  Temperature in landscape elements
    float * temp = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "tem");
    strcat(fileName, timeName);
    fileTemp.open(fileName, ios::out | ios::binary);
    if (!fileTemp.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //    if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        //      if (Dew[k].GetAccumulatedTemperature() > missingData)
        if (Dew[k].GetTemperature() > missingData)
        {
            temp[k] = (float)((Dew[k].GetTemperature()));
        }
        else
        {
            temp[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    temp[ELEMENT(i,j)] = noData;
    }
    }
    else {
    temp[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileTemp.write(reinterpret_cast<char *>(temp), sizeof(float)*numLand);
    fileTemp.close();
    delete[] temp;

    //  Precipitation in landscape elements
    float * pre = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "pre");
    strcat(fileName, timeName);
    filePre.open(fileName, ios::out | ios::binary);
    if (!filePre.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //    if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        //      if (Dew[k].GetAccumulatedPrecipitation() > missingData)
        if (Dew[k].GetPrecipitation() > missingData)
        {
            pre[k] = (float)((Dew[k].GetPrecipitation()) * 1000.0);
        }
        else
        {
            pre[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    pre[ELEMENT(i,j)] = noData;
    }
    }
    else {
    pre[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    filePre.write(reinterpret_cast<char *>(pre), sizeof(float)*numLand);
    filePre.close();
    delete[] pre;

    //  Evapotranspiration from landscape elements
    float * eva = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "eva");
    strcat(fileName, timeName);
    fileEva.open(fileName, ios::out | ios::binary);
    if (!fileEva.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //      if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      //        if (Dew[k].GetInterceptionLoss() > missingData &&
      //                Dew[k].GetTranspSoilEvap() > missingData &&
      //                Dew[k].GetLakeEvap() > missingData)
            //      eva[k] = (float) ((Dew[k].GetInterceptionLoss() + Dew[k].GetTranspSoilEvap() +
            //                         Dew[k].GetLakeEvap())*1000.0);
        if (Dew[k].GetEvapotranspiration() > missingData)
        {
            eva[k] = (float)((Dew[k].GetEvapotranspiration()) * 1000.0);
        }
        else
        {
            eva[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    eva[ELEMENT(i,j)] = noData;
    }
    }
    else {
    eva[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileEva.write(reinterpret_cast<char *>(eva), sizeof(float)*numLand);
    fileEva.close();
    delete[] eva;

    //  Snow store in landscape elements
    float * swe = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "swe");
    strcat(fileName, timeName);
    fileSwe.open(fileName, ios::out | ios::binary);
    if (!fileSwe.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetSnowStore() > missingData)
        {
            swe[k] = (float)((Dew[k].GetSnowStore() + Dew[k].GetMeltWater()) * 1000.0);
        }
        else
        {
            swe[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    swe[ELEMENT(i,j)] = noData;
    }
    }
    else {
    swe[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileSwe.write(reinterpret_cast<char *>(swe), sizeof(float)*numLand);
    fileSwe.close();
    delete[] swe;

    //  Snowcover fraction in landscape elements
    float * scf = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "scf");
    strcat(fileName, timeName);
    fileScf.open(fileName, ios::out | ios::binary);
    if (!fileScf.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetSnowCoverFraction() > missingData)
        {
            scf[k] = (float)((Dew[k].GetSnowCoverFraction()) * 100.0);
        }
        else
        {
            scf[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    scf[ELEMENT(i,j)] = noData;
    }
    }
    else {
    scf[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileScf.write(reinterpret_cast<char *>(scf), sizeof(float)*numLand);
    fileScf.close();
    delete[] scf;

    //  Snow meltwater in landscape elements
    float * smw = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "smw");
    strcat(fileName, timeName);
    fileSmw.open(fileName, ios::out | ios::binary);
    if (!fileSmw.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetMeltWater() > missingData)
        {
            smw[k] = (float)((Dew[k].GetMeltWater()) * 1000.0);
        }
        else
        {
            smw[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    smw[ELEMENT(i,j)] = noData;
    }
    }
    else {
    smw[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileSmw.write(reinterpret_cast<char *>(smw), sizeof(float)*numLand);
    fileSmw.close();
    delete[] smw;

    //  Runoff from landscape elements
    float * runoff = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "run");
    strcat(fileName, timeName);
    fileRun.open(fileName, ios::out | ios::binary);
    if (!fileRun.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetRunoff() > missingData)
        {
            runoff[k] = (float)((Dew[k].GetRunoff()) * 1000.0);
        }
        else
        {
            runoff[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    runoff[ELEMENT(i,j)] = noData;
    }
    }
    else {
    runoff[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileRun.write(reinterpret_cast<char *>(runoff), sizeof(float)*numLand);
    fileRun.close();
    delete[] runoff;

    //  Discharge from landscape elements
    float * discharge = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "dew");
    strcat(fileName, timeName);
    fileDew.open(fileName, ios::out | ios::binary);
    if (!fileDew.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetDischarge() > missingData)
        {
            discharge[k] = (float)(Dew[k].GetDischarge());
        }
        else
        {
            discharge[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    discharge[ELEMENT(i,j)] = noData;
    }
    }
    else {
    discharge[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileDew.write(reinterpret_cast<char *>(discharge), sizeof(float)*numLand);
    fileDew.close();
    delete[] discharge;

    //  Soil moisture deficit in landscape elements
    float * hsd = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "hsd");
    strcat(fileName, timeName);
    fileHsd.open(fileName, ios::out | ios::binary);
    if (!fileHsd.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetHbvSoilMoistureDeficit() > missingData)
        {
            if (Dew[k].GetHbvSoilMoistureDeficit() >= 0.0)
            {
                hsd[k] = (float)((Dew[k].GetHbvSoilMoistureDeficit()) * 1000.0);
            }
            else
            {
                hsd[k] = 0.0;
            }
        }
        else
        {
            hsd[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    hsd[ELEMENT(i,j)] = noData;
    }
    }
    else {
    hsd[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileHsd.write(reinterpret_cast<char *>(hsd), sizeof(float)*numLand);
    fileHsd.close();
    delete[] hsd;

    //  Soil moisture in landscape elements
    float * hsm = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "hsm");
    strcat(fileName, timeName);
    fileHsm.open(fileName, ios::out | ios::binary);
    if (!fileHsm.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetHbvSoilMoisture() > missingData)
        {
            if (Dew[k].GetHbvSoilMoisture() >= 0.0)
            {
                hsm[k] = (float)((Dew[k].GetHbvSoilMoisture()) * 1000.0);
            }
            else
            {
                hsm[k] = 0.0;
            }
        }
        else
        {
            hsm[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    hsm[ELEMENT(i,j)] = noData;
    }
    }
    else {
    hsm[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileHsm.write(reinterpret_cast<char *>(hsm), sizeof(float)*numLand);
    fileHsm.close();
    delete[] hsm;

    //  Groundwater in landscape elements
    float * hgw = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "hgw");
    strcat(fileName, timeName);
    fileHgw.open(fileName, ios::out | ios::binary);
    if (!fileHgw.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetHbvUpperZone() > missingData && Dew[k].GetHbvLowerZone() > missingData)
        {
            hgw[k] = (float)((Dew[k].GetHbvUpperZone() + Dew[k].GetHbvLowerZone()) * 1000.0);
        }
        else
        {
            hgw[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    hgw[ELEMENT(i,j)] = noData;
    }
    }
    else {
    hgw[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileHgw.write(reinterpret_cast<char *>(hgw), sizeof(float)*numLand);
    fileHgw.close();
    delete[] hgw;

    //  Lower zone groundwater in landscape elements
    float * hlz = new float[numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "hlz");
    strcat(fileName, timeName);
    fileHlz.open(fileName, ios::out | ios::binary);
    if (!fileHlz.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    /*  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {*/
    while (k < numLand)
    {
        //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
        if (Dew[k].GetHbvLowerZone() > missingData)
        {
            hlz[k] = (float)((Dew[k].GetHbvLowerZone()) * 1000.0);
        }
        else
        {
            hlz[k] = (float)noData;
        }
        k++;
    }
    /*        else {
    hlz[ELEMENT(i,j)] = noData;
    }
    }
    else {
    hlz[ELEMENT(i,j)] = noData;
    }
    }
    }*/
    fileHlz.write(reinterpret_cast<char *>(hlz), sizeof(float)*numLand);
    fileHlz.close();
    delete[] hlz;
}


void WriteBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,
                     int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
    int i, j, k;
    char fileName[100];
    char timeName[100];
    unsigned short int noData16 = 65535;
    ofstream fileTemp, filePre, fileEva, fileSwe, fileScf, fileSmw, fileRun, fileDew, fileHsd, fileHsm, fileHuz, fileHlz;

    sprintf(timeName, "_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());

    //  Temperature in landscape elements
    unsigned short int * temp = new unsigned short int[nRows * nCols];
    strcpy(fileName, "tem");
    strcat(fileName, timeName);
    fileTemp.open(fileName, ios::out | ios::binary);
    if (!fileTemp.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedTemperature() > missingData)
                    if (Dew[k].GetTemperature() > missingData)
                    {
                        temp[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetTemperature() * 10.0) + 2731);    //temp10K[ELEMENT(i,j)]-2731)/10.0
                    }
                    else
                    {
                        temp[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    temp[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                temp[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileTemp.write(reinterpret_cast<char *>(temp), sizeof(unsigned short int)*nRows * nCols);
    fileTemp.close();
    delete[] temp;

    //  Precipitation in landscape elements
    unsigned short int * pre = new unsigned short int[nRows * nCols];
    strcpy(fileName, "pre");
    strcat(fileName, timeName);
    filePre.open(fileName, ios::out | ios::binary);
    if (!filePre.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetPrecipitation() > missingData)
                    {
                        pre[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetPrecipitation()) * 1000.0);
                    }
                    else
                    {
                        pre[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    pre[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                pre[ELEMENT(i, j)] = noData16;
            }
        }
    }
    filePre.write(reinterpret_cast<char *>(pre), sizeof(unsigned short int)*nRows * nCols);
    filePre.close();
    delete[] pre;

    //  Evapotranspiration from landscape elements
    unsigned short int * eva = new unsigned short int[nRows * nCols];
    strcpy(fileName, "eva");
    strcat(fileName, timeName);
    fileEva.open(fileName, ios::out | ios::binary);
    if (!fileEva.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
		  //                    if (Dew[k].GetInterceptionLoss() > missingData &&
		  //                            Dew[k].GetTranspSoilEvap() > missingData &&
		  //                            Dew[k].GetLakeEvap() > missingData)
                        //            eva[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetInterceptionLoss() + Dew[k].GetTranspSoilEvap() +
                        //                                                       Dew[k].GetLakeEvap())*1000.0);
		    if (Dew[k].GetEvapotranspiration() > missingData)
                    {
                        eva[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetEvapotranspiration()) * 1000.0);
                    }
                    else
                    {
                        eva[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    eva[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                eva[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileEva.write(reinterpret_cast<char *>(eva), sizeof(unsigned short int)*nRows * nCols);
    fileEva.close();
    delete[] eva;

    //  Snow store in landscape elements
    unsigned short int * swe = new unsigned short int[nRows * nCols];
    strcpy(fileName, "swe");
    strcat(fileName, timeName);
    fileSwe.open(fileName, ios::out | ios::binary);
    if (!fileSwe.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetSnowStore() > missingData)
                    {
                        swe[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetSnowStore() + Dew[k].GetMeltWater()) * 1000.0);
                    }
                    else
                    {
                        swe[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    swe[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                swe[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileSwe.write(reinterpret_cast<char *>(swe), sizeof(unsigned short int)*nRows * nCols);
    fileSwe.close();
    delete[] swe;

    //  Snowcover fraction in landscape elements
    unsigned short int * scf = new unsigned short int[nRows * nCols];
    strcpy(fileName, "scf");
    strcat(fileName, timeName);
    fileScf.open(fileName, ios::out | ios::binary);
    if (!fileScf.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetSnowCoverFraction() > missingData)
                    {
                        scf[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetSnowCoverFraction()) * 100.0);
                    }
                    else
                    {
                        scf[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    scf[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                scf[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileScf.write(reinterpret_cast<char *>(scf), sizeof(unsigned short int)*nRows * nCols);
    fileScf.close();
    delete[] scf;

    //  Snow meltwater in landscape elements
    unsigned short int * smw = new unsigned short int[nRows * nCols];
    strcpy(fileName, "smw");
    strcat(fileName, timeName);
    fileSmw.open(fileName, ios::out | ios::binary);
    if (!fileSmw.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetMeltWater() > missingData)
                    {
                        smw[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetMeltWater()) * 1000.0);
                    }
                    else
                    {
                        smw[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    smw[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                smw[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileSmw.write(reinterpret_cast<char *>(smw), sizeof(unsigned short int)*nRows * nCols);
    fileSmw.close();
    delete[] smw;

    //  Runoff from landscape elements
    unsigned short int * runoff = new unsigned short int[nRows * nCols];
    strcpy(fileName, "run");
    strcat(fileName, timeName);
    fileRun.open(fileName, ios::out | ios::binary);
    if (!fileRun.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetRunoff() > missingData)
                    {
                        runoff[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetRunoff()) * 1000.0);
                    }
                    else
                    {
                        runoff[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    runoff[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                runoff[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileRun.write(reinterpret_cast<char *>(runoff), sizeof(unsigned short int)*nRows * nCols);
    fileRun.close();
    delete[] runoff;

    //  Discharge from landscape elements
    unsigned short int * discharge = new unsigned short int[nRows * nCols];
    strcpy(fileName, "dew");
    strcat(fileName, timeName);
    fileDew.open(fileName, ios::out | ios::binary);
    if (!fileDew.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetDischarge() > missingData)
                    {
                        discharge[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetDischarge()) * 1000.0);
                    }
                    else
                    {
                        discharge[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    discharge[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                discharge[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileDew.write(reinterpret_cast<char *>(discharge), sizeof(unsigned short int)*nRows * nCols);
    fileDew.close();
    delete[] discharge;

    //  Soil moisture deficit in landscape elements
    unsigned short int * hsd = new unsigned short int[nRows * nCols];
    strcpy(fileName, "hsd");
    strcat(fileName, timeName);
    fileHsd.open(fileName, ios::out | ios::binary);
    if (!fileHsd.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetHbvSoilMoistureDeficit() > missingData)
                    {
                        hsd[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetHbvSoilMoistureDeficit()) * 1000.0);
                    }
                    else
                    {
                        hsd[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    hsd[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                hsd[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileHsd.write(reinterpret_cast<char *>(hsd), sizeof(unsigned short int)*nRows * nCols);
    fileHsd.close();
    delete[] hsd;

    //  Soil moisture in landscape elements
    unsigned short int * hsm = new unsigned short int[nRows * nCols];
    strcpy(fileName, "hsm");
    strcat(fileName, timeName);
    fileHsm.open(fileName, ios::out | ios::binary);
    if (!fileHsm.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetHbvSoilMoisture() > missingData)
                    {
                        hsm[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetHbvSoilMoisture()) * 1000.0);
                    }
                    else
                    {
                        hsm[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    hsm[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                hsm[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileHsm.write(reinterpret_cast<char *>(hsm), sizeof(unsigned short int)*nRows * nCols);
    fileHsm.close();
    delete[] hsm;

    //  Upper zone groundwater in landscape elements
    unsigned short int * huz = new unsigned short int[nRows * nCols];
    strcpy(fileName, "huz");
    strcat(fileName, timeName);
    fileHuz.open(fileName, ios::out | ios::binary);
    if (!fileHuz.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetHbvUpperZone() > missingData)
                    {
                        huz[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetHbvUpperZone()) * 1000.0);
                    }
                    else
                    {
                        huz[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    huz[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                huz[ELEMENT(i, j)] = noData16;
            }
        }
        fout << endl;
    }
    fout << endl;
    fileHuz.write(reinterpret_cast<char *>(huz), sizeof(unsigned short int)*nRows * nCols);
    fileHuz.close();
    delete[] huz;

    //  Lower zone groundwater in landscape elements
    unsigned short int * hlz = new unsigned short int[nRows * nCols];
    strcpy(fileName, "hlz");
    strcat(fileName, timeName);
    fileHlz.open(fileName, ios::out | ios::binary);
    if (!fileHlz.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    //          if (Dew[k].GetAccumulatedPrecipitation() > missingData)
                    if (Dew[k].GetHbvLowerZone() > missingData)
                    {
                        hlz[ELEMENT(i, j)] = (unsigned short int) ((Dew[k].GetHbvLowerZone()) * 1000.0);
                    }
                    else
                    {
                        hlz[ELEMENT(i, j)] = noData16;
                    }
                    k++;
                }
                else
                {
                    hlz[ELEMENT(i, j)] = noData16;
                }
            }
            else
            {
                hlz[ELEMENT(i, j)] = noData16;
            }
        }
    }
    fileHlz.write(reinterpret_cast<char *>(hlz), sizeof(unsigned short int)*nRows * nCols);
    fileHlz.close();
    delete[] hlz;
}


void WriteAsciiGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,
                    int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
    int i, j, k;
    char fileName[100];
    char dirName[100];
    char timeName[100];

    //  sprintf(dirName,"./%04d/",datetime.getYear());
    sprintf(dirName, "./");
    sprintf(timeName, "_%04d_%02d_%02d.asc", datetime.getYear(), datetime.getMonth(), datetime.getDay());

    //  fout << "Temperature in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "tem");
    strcat(fileName, timeName);
    ofstream fileTemp(fileName);
    if (!fileTemp.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileTemp.precision(0); 
    fileTemp.setf(ios::fixed); 
    fileTemp << "ncols         " << nCols << endl;
    fileTemp << "nrows         " << nRows << endl;
    fileTemp << "xllcorner     " << xllCorner << endl;
    fileTemp << "yllcorner     " << yllCorner << endl;
    fileTemp << "cellsize      " << cellSize << endl;
    fileTemp << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetTemperature() > missingData)
                    {
                        fileTemp.width(15);
                        fileTemp.precision(5);
                        fileTemp.setf(ios::showpoint);
                        fileTemp.setf(ios::fixed);
                        fileTemp << (Dew[k].GetTemperature()) << endl;
                    }
                    else
                    {
                        fileTemp.width(15);
                        fileTemp << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileTemp.width(15);
                    fileTemp << noData << endl;
                }
            }
            else
            {
                fileTemp.width(15);
                fileTemp << noData << endl;
            }
        }
    }
    fileTemp.close();

    //  fout << "Precipitation in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "pre");
    strcat(fileName, timeName);
    ofstream filePrec(fileName);
    if (!filePrec.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    filePrec.precision(0); 
    filePrec.setf(ios::fixed); 
    filePrec << "ncols         " << nCols << endl;
    filePrec << "nrows         " << nRows << endl;
    filePrec << "xllcorner     " << xllCorner << endl;
    filePrec << "yllcorner     " << yllCorner << endl;
    filePrec << "cellsize      " << cellSize << endl;
    filePrec << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetPrecipitation() > missingData)
                    {
                        filePrec.width(15);
                        filePrec.precision(5);
                        filePrec.setf(ios::showpoint);
                        filePrec.setf(ios::fixed);
                        filePrec << (Dew[k].GetPrecipitation()) * 1000.0 << endl;
                    }
                    else
                    {
                        filePrec.width(15);
                        filePrec << noData << endl;
                    }
                    k++;
                }
                else
                {
                    filePrec.width(15);
                    filePrec << noData << endl;
                }
            }
            else
            {
                filePrec.width(15);
                filePrec << noData << endl;
            }
        }
    }
    filePrec.close();

    //  fout << "Evapotranspiration from landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "eva");
    strcat(fileName, timeName);
    ofstream fileEvap(fileName);
    if (!fileEvap.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileEvap.precision(0);
    fileEvap.setf(ios::fixed); 
    fileEvap << "ncols         " << nCols << endl;
    fileEvap << "nrows         " << nRows << endl;
    fileEvap << "xllcorner     " << xllCorner << endl;
    fileEvap << "yllcorner     " << yllCorner << endl;
    fileEvap << "cellsize      " << cellSize << endl;
    fileEvap << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
		  //                    if (Dew[k].GetInterceptionLoss() > missingData &&
		  //                            Dew[k].GetTranspSoilEvap() > missingData &&
		  //                            Dew[k].GetLakeEvap() > missingData)
		    if (Dew[k].GetEvapotranspiration() > missingData)
                    {
                        fileEvap.width(15);
                        fileEvap.precision(5);
                        fileEvap.setf(ios::showpoint);
                        fileEvap.setf(ios::fixed);
                        //            fileEvap << (Dew[k].GetInterceptionLoss() + Dew[k].GetTranspSoilEvap() +
                        //                         Dew[k].GetLakeEvap())*1000.0 << endl;
                        fileEvap << (Dew[k].GetEvapotranspiration()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileEvap.width(15);
                        fileEvap << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileEvap.width(15);
                    fileEvap << noData << endl;
                }
            }
            else
            {
                fileEvap.width(15);
                fileEvap << noData << endl;
            }
        }
    }
    fileEvap.close();

    //  fout << "Snow store in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "swe");
    strcat(fileName, timeName);
    ofstream fileSwe(fileName);
    if (!fileSwe.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileSwe.precision(0); 
    fileSwe.setf(ios::fixed); 
    fileSwe << "ncols         " << nCols << endl;
    fileSwe << "nrows         " << nRows << endl;
    fileSwe << "xllcorner     " << xllCorner << endl;
    fileSwe << "yllcorner     " << yllCorner << endl;
    fileSwe << "cellsize      " << cellSize << endl;
    fileSwe << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetSnowStore() > missingData)
                    {
                        fileSwe.width(15);
                        fileSwe.precision(5);
                        fileSwe.setf(ios::showpoint);
                        fileSwe.setf(ios::fixed);
                        fileSwe << (Dew[k].GetSnowStore() + Dew[k].GetMeltWater()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileSwe.width(15);
                        fileSwe << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileSwe.width(15);
                    fileSwe << noData << endl;
                }
            }
            else
            {
                fileSwe.width(15);
                fileSwe << noData << endl;
            }
        }
    }
    fileSwe.close();

    //  fout << "Snowcover fraction in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "scf");
    strcat(fileName, timeName);
    ofstream fileScf(fileName);
    if (!fileScf.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileScf.precision(0); 
    fileScf.setf(ios::fixed); 
    fileScf << "ncols         " << nCols << endl;
    fileScf << "nrows         " << nRows << endl;
    fileScf << "xllcorner     " << xllCorner << endl;
    fileScf << "yllcorner     " << yllCorner << endl;
    fileScf << "cellsize      " << cellSize << endl;
    fileScf << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetSnowCoverFraction() > missingData)
                    {
                        fileScf.width(15);
                        fileScf.precision(5);
                        fileScf.setf(ios::showpoint);
                        fileScf.setf(ios::fixed);
                        fileScf << (Dew[k].GetSnowCoverFraction()) * 100.0 << endl;
                    }
                    else
                    {
                        fileScf.width(15);
                        fileScf << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileScf.width(15);
                    fileScf << noData << endl;
                }
            }
            else
            {
                fileScf.width(15);
                fileScf << noData << endl;
            }
        }
    }
    fileScf.close();

    //  fout << "Snow meltwater in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "smw");
    strcat(fileName, timeName);
    ofstream fileSmw(fileName);
    if (!fileSmw.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileSmw.precision(0); 
    fileSmw.setf(ios::fixed); 
    fileSmw << "ncols         " << nCols << endl;
    fileSmw << "nrows         " << nRows << endl;
    fileSmw << "xllcorner     " << xllCorner << endl;
    fileSmw << "yllcorner     " << yllCorner << endl;
    fileSmw << "cellsize      " << cellSize << endl;
    fileSmw << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetMeltWater() > missingData)
                    {
                        fileSmw.width(15);
                        fileSmw.precision(5);
                        fileSmw.setf(ios::showpoint);
                        fileSmw.setf(ios::fixed);
                        fileSmw << (Dew[k].GetMeltWater()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileSmw.width(15);
                        fileSmw << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileSmw.width(15);
                    fileSmw << noData << endl;
                }
            }
            else
            {
                fileSmw.width(15);
                fileSmw << noData << endl;
            }
        }
    }
    fileSmw.close();

    //  fout << "Runoff from landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "run");
    strcat(fileName, timeName);
    ofstream fileRunoff(fileName);
    if (!fileRunoff.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileRunoff.precision(0); 
    fileRunoff.setf(ios::fixed); 
    fileRunoff << "ncols         " << nCols << endl;
    fileRunoff << "nrows         " << nRows << endl;
    fileRunoff << "xllcorner     " << xllCorner << endl;
    fileRunoff << "yllcorner     " << yllCorner << endl;
    fileRunoff << "cellsize      " << cellSize << endl;
    fileRunoff << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetRunoff() > missingData)
                    {
                        fileRunoff.width(15);
                        fileRunoff.precision(5);
                        fileRunoff.setf(ios::showpoint);
                        fileRunoff.setf(ios::fixed);
                        fileRunoff << (Dew[k].GetRunoff()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileRunoff.width(15);
                        fileRunoff << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileRunoff.width(15);
                    fileRunoff << noData << endl;
                }
            }
            else
            {
                fileRunoff.width(15);
                fileRunoff << noData << endl;
            }
        }
    }
    fileRunoff.close();

    //  fout << "Discharge from landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "dew");
    strcat(fileName, timeName);
    ofstream fileDisch(fileName);
    if (!fileDisch.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileDisch.precision(0); 
    fileDisch.setf(ios::fixed); 
    fileDisch << "ncols         " << nCols << endl;
    fileDisch << "nrows         " << nRows << endl;
    fileDisch << "xllcorner     " << xllCorner << endl;
    fileDisch << "yllcorner     " << yllCorner << endl;
    fileDisch << "cellsize      " << cellSize << endl;
    fileDisch << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetDischarge() > missingData)
                    {
                        fileDisch.width(15);
                        fileDisch.precision(5);
                        fileDisch.setf(ios::showpoint);
                        fileDisch.setf(ios::fixed);
                        fileDisch << Dew[k].GetDischarge() << endl;
                    }
                    else
                    {
                        fileDisch.width(15);
                        fileDisch << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileDisch.width(15);
                    fileDisch << noData << endl;
                }
            }
            else
            {
                fileDisch.width(15);
                fileDisch << noData << endl;
            }
        }
    }
    fileDisch.close();

    //  fout << "Soil moisture deficit in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "hsd");
    strcat(fileName, timeName);
    ofstream fileHsd(fileName);
    if (!fileHsd.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileHsd.precision(0); 
    fileHsd.setf(ios::fixed); 
    fileHsd << "ncols         " << nCols << endl;
    fileHsd << "nrows         " << nRows << endl;
    fileHsd << "xllcorner     " << xllCorner << endl;
    fileHsd << "yllcorner     " << yllCorner << endl;
    fileHsd << "cellsize      " << cellSize << endl;
    fileHsd << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetHbvSoilMoisture() > missingData)
                    {
                        fileHsd.width(15);
                        fileHsd.precision(5);
                        fileHsd.setf(ios::showpoint);
                        fileHsd.setf(ios::fixed);
                        fileHsd << (Dew[k].GetHbvSoilMoistureDeficit()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileHsd.width(15);
                        fileHsd << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileHsd.width(15);
                    fileHsd << noData << endl;
                }
            }
            else
            {
                fileHsd.width(15);
                fileHsd << noData << endl;
            }
        }
    }
    fileHsd.close();

    //  fout << "Soil moisture in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "hsm");
    strcat(fileName, timeName);
    ofstream fileHsm(fileName);
    if (!fileHsm.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileHsm.precision(0); 
    fileHsm.setf(ios::fixed); 
    fileHsm << "ncols         " << nCols << endl;
    fileHsm << "nrows         " << nRows << endl;
    fileHsm << "xllcorner     " << xllCorner << endl;
    fileHsm << "yllcorner     " << yllCorner << endl;
    fileHsm << "cellsize      " << cellSize << endl;
    fileHsm << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetHbvSoilMoisture() > missingData)
                    {
                        fileHsm.width(15);
                        fileHsm.precision(5);
                        fileHsm.setf(ios::showpoint);
                        fileHsm.setf(ios::fixed);
                        fileHsm << (Dew[k].GetHbvSoilMoisture()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileHsm.width(15);
                        fileHsm << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileHsm.width(15);
                    fileHsm << noData << endl;
                }
            }
            else
            {
                fileHsm.width(15);
                fileHsm << noData << endl;
            }
        }
    }
    fileHsm.close();

    //  fout << "Upper zone in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "huz");
    strcat(fileName, timeName);
    ofstream fileHuz(fileName);
    if (!fileHuz.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileHuz.precision(0); 
    fileHuz.setf(ios::fixed); 
    fileHuz << "ncols         " << nCols << endl;
    fileHuz << "nrows         " << nRows << endl;
    fileHuz << "xllcorner     " << xllCorner << endl;
    fileHuz << "yllcorner     " << yllCorner << endl;
    fileHuz << "cellsize      " << cellSize << endl;
    fileHuz << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetHbvUpperZone() > missingData)
                    {
                        fileHuz.width(15);
                        fileHuz.precision(5);
                        fileHuz.setf(ios::showpoint);
                        fileHuz.setf(ios::fixed);
                        fileHuz << (Dew[k].GetHbvUpperZone()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileHuz.width(15);
                        fileHuz << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileHuz.width(15);
                    fileHuz << noData << endl;
                }
            }
            else
            {
                fileHuz.width(15);
                fileHuz << noData << endl;
            }
        }
        fileHuz << endl;
    }
    fileHuz << endl;
    fileHuz.close();

    //  fout << "Lower zone in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "hlz");
    strcat(fileName, timeName);
    ofstream fileHlz(fileName);
    if (!fileHlz.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileHlz.precision(0); 
    fileHlz.setf(ios::fixed); 
    fileHlz << "ncols         " << nCols << endl;
    fileHlz << "nrows         " << nRows << endl;
    fileHlz << "xllcorner     " << xllCorner << endl;
    fileHlz << "yllcorner     " << yllCorner << endl;
    fileHlz << "cellsize      " << cellSize << endl;
    fileHlz << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetHbvLowerZone() > missingData)
                    {
                        fileHlz.width(15);
                        fileHlz.precision(5);
                        fileHlz.setf(ios::showpoint);
                        fileHlz.setf(ios::fixed);
                        fileHlz << (Dew[k].GetHbvLowerZone()) * 1000.0 << endl;
                    }
                    else
                    {
                        fileHlz.width(15);
                        fileHlz << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileHlz.width(15);
                    fileHlz << noData << endl;
                }
            }
            else
            {
                fileHlz.width(15);
                fileHlz << noData << endl;
            }
        }
    }
    fileHlz.close();

    //  fout << "Glacier snow store in landscape elements at time step " << timeStep << ":\n";
    /*    strcpy(fileName,dirName);
    strcat(fileName,"gsw");
    strcat(fileName,timeName);
    ofstream fileGlacSnowStore(fileName);
    if (!fileGlacSnowStore.is_open()) 
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit (1);
    }
    fileGlacSnowStore.precision(0); 
    fileGlacSnowStore.setf(ios::fixed); 
    fileGlacSnowStore << "ncols         " << nCols << endl;
    fileGlacSnowStore << "nrows         " << nRows << endl;
    fileGlacSnowStore << "xllcorner     " << xllCorner << endl;
    fileGlacSnowStore << "yllcorner     " << yllCorner << endl;
    fileGlacSnowStore << "cellsize      " << cellSize << endl;
    fileGlacSnowStore << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
        {
            if (k < numLand) 
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i,j)) 
                {
                    if (Dew[k].GetAnnualMassBalance() > missingData) 
		    {
		        fileGlacSnowStore.width(15); 
			fileGlacSnowStore.precision(5); 
			fileGlacSnowStore.setf(ios::showpoint); 
			fileGlacSnowStore.setf(ios::fixed); 
			fileGlacSnowStore << (Dew[k].GetGlacierSnowStore()+Dew[k].GetGlacierSnowMeltWater())*1000.0 << endl;
		    }
		    else {
		        fileGlacSnowStore.width(15); 
			fileGlacSnowStore << noData << endl;
		    }
		    k++;
		}
		else 
		{
		    fileGlacSnowStore.width(15); 
		    fileGlacSnowStore << noData << endl;
		}
	    }
	    else 
	        {
	        fileGlacSnowStore.width(15); 
		fileGlacSnowStore << noData << endl;
	    }
	}
    }
    fileGlacSnowStore.close();*/

    //  fout << "Glacier mass balance in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "gmb");
    strcat(fileName, timeName);
    ofstream fileGlacMassBal(fileName);
    if (!fileGlacMassBal.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileGlacMassBal.precision(0); 
    fileGlacMassBal.setf(ios::fixed); 
    fileGlacMassBal << "ncols         " << nCols << endl;
    fileGlacMassBal << "nrows         " << nRows << endl;
    fileGlacMassBal << "xllcorner     " << xllCorner << endl;
    fileGlacMassBal << "yllcorner     " << yllCorner << endl;
    fileGlacMassBal << "cellsize      " << cellSize << endl;
    fileGlacMassBal << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetAnnualMassBalance() > missingData)
                    {
                        fileGlacMassBal.width(15);
                        fileGlacMassBal.precision(5);
                        fileGlacMassBal.setf(ios::showpoint);
                        fileGlacMassBal.setf(ios::fixed);
                        //            fileGlacMassBal << (Dew[k].GetGlacierMassBalance())*1000.0 << endl;
                        fileGlacMassBal << (Dew[k].GetGlacierMassBalance()) << endl;
                    }
                    else
                    {
                        fileGlacMassBal.width(15);
                        fileGlacMassBal << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileGlacMassBal.width(15);
                    fileGlacMassBal << noData << endl;
                }
            }
            else
            {
                fileGlacMassBal.width(15);
                fileGlacMassBal << noData << endl;
            }
        }
    }
    fileGlacMassBal.close();

    //  fout << "Glacier ice melt in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "gim");
    strcat(fileName, timeName);
    ofstream fileIceMelt(fileName);
    if (!fileIceMelt.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileIceMelt.precision(0); 
    fileIceMelt.setf(ios::fixed); 
    fileIceMelt << "ncols         " << nCols << endl;
    fileIceMelt << "nrows         " << nRows << endl;
    fileIceMelt << "xllcorner     " << xllCorner << endl;
    fileIceMelt << "yllcorner     " << yllCorner << endl;
    fileIceMelt << "cellsize      " << cellSize << endl;
    fileIceMelt << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetAnnualMassBalance() > missingData)
                    {
                        fileIceMelt.width(15);
                        fileIceMelt.precision(5);
                        fileIceMelt.setf(ios::showpoint);
                        fileIceMelt.setf(ios::fixed);
                        fileIceMelt << (Dew[k].GetGlacierIceMelt()) * 1000.0 << endl;
                        //            fileIceMelt << (Dew[k].GetGlacierIceMelt()) << endl;
                    }
                    else
                    {
                        fileIceMelt.width(15);
                        fileIceMelt << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileIceMelt.width(15);
                    fileIceMelt << noData << endl;
                }
            }
            else
            {
                fileIceMelt.width(15);
                fileIceMelt << noData << endl;
            }
        }
    }
    fileIceMelt.close();
}


void WriteAsciiGridWaterBalance(DistributedElement * const Dew, DateTime startSimulationTime, DateTime endSimulationTime, int numLand, int nRows,  
				int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
    int i,j,k;
    double numberDaysPerYear=365.25;
    char fileName[100];
    char dirName[100];
    char timeName[100];
  
    sprintf(dirName,"./");
    sprintf(timeName,"_%04d_%02d_%02d_%04d_%02d_%02d.asc",startSimulationTime.getYear(),startSimulationTime.getMonth(),startSimulationTime.getDay(),
  	  endSimulationTime.getYear(),endSimulationTime.getMonth(),endSimulationTime.getDay());
  
    //  fout << "Mean annual precipitation from landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"meanAnnualPrecipitation");
    strcat(fileName,timeName);
    ofstream filePrec(fileName);
    if (!filePrec.is_open()) 
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit (1);
    }
    filePrec.precision(0); 
    filePrec.setf(ios::fixed); 
    filePrec << "ncols         " << nCols << endl;
    filePrec << "nrows         " << nRows << endl;
    filePrec << "xllcorner     " << xllCorner << endl;
    filePrec << "yllcorner     " << yllCorner << endl;
    filePrec << "cellsize      " << cellSize << endl;
    filePrec << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
            if (k < numLand) 
	    {
	        if (Dew[k].GetGeoIndex() == ELEMENT(i,j)) 
		{
		    if (Dew[k].GetSumPrecipitation() > missingData) 
		    {
		        filePrec.width(15); 
			filePrec.precision(5); 
			filePrec.setf(ios::showpoint); 
			filePrec.setf(ios::fixed); 
			filePrec << (Dew[k].GetSumPrecipitation()/Dew[k].GetNumberSum())*numberDaysPerYear*1000.0 << endl;
		    }
		    else 
		    {
  		        filePrec.width(15); filePrec << noData << endl;
		    }
		    k++;
		}
		else 
		{
		    filePrec.width(15); 
		    filePrec << noData << endl;
		}
	    }
	    else 
	    {
	        filePrec.width(15); 
		filePrec << noData << endl;
	    }
	}
    }
    filePrec.close();
  
    //  fout << "Mean annual evapotranspiration from landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"meanAnnualEvapotranspiration");
    strcat(fileName,timeName);
    ofstream fileEvap(fileName);
    if (!fileEvap.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    fileEvap.precision(0); 
    fileEvap.setf(ios::fixed); 
    fileEvap << "ncols         " << nCols << endl;
    fileEvap << "nrows         " << nRows << endl;
    fileEvap << "xllcorner     " << xllCorner << endl;
    fileEvap << "yllcorner     " << yllCorner << endl;
    fileEvap << "cellsize      " << cellSize << endl;
    fileEvap << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
            if (k < numLand) 
	    {
	        if (Dew[k].GetGeoIndex() == ELEMENT(i,j)) 
		{
		    if (Dew[k].GetSumEvapotranspiration() > missingData) 
		    {
		        fileEvap.width(15); 
			fileEvap.precision(5); 
			fileEvap.setf(ios::showpoint); 
			fileEvap.setf(ios::fixed); 
			fileEvap << (Dew[k].GetSumEvapotranspiration()/Dew[k].GetNumberSum())*numberDaysPerYear*1000.0 << endl;
		    }
		    else 
		    {
		        fileEvap.width(15); 
			fileEvap << noData << endl;
                     }
		    k++;
		}
		else 
		{
		    fileEvap.width(15); 
		    fileEvap << noData << endl;
		}
	    }
	    else 
	    {
	        fileEvap.width(15); 
		fileEvap << noData << endl;
	    }
	}
    }
    fileEvap.close();
  
    //  fout << "Mean annual runoff from landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"meanAnnualRunoff");
    strcat(fileName,timeName);
    ofstream fileRunoff(fileName);
    if (!fileRunoff.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    fileRunoff.precision(0); 
    fileRunoff.setf(ios::fixed); 
    fileRunoff << "ncols         " << nCols << endl;
    fileRunoff << "nrows         " << nRows << endl;
    fileRunoff << "xllcorner     " << xllCorner << endl;
    fileRunoff << "yllcorner     " << yllCorner << endl;
    fileRunoff << "cellsize      " << cellSize << endl;
    fileRunoff << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
	    if (k < numLand) 
	    {
	        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) 
		{
		    if (Dew[k].GetSumRunoff() > missingData) 
		    {
		        fileRunoff.width(15); 
			fileRunoff.precision(5); 
			fileRunoff.setf(ios::showpoint); 
			fileRunoff.setf(ios::fixed); 
			fileRunoff << (Dew[k].GetSumRunoff()/Dew[k].GetNumberSum())*numberDaysPerYear*1000.0 << endl;
		    }
		    else 
		    {
		        fileRunoff.width(15); 
			fileRunoff << noData << endl;
		    }
		    k++;
		}
		else 
		{
		    fileRunoff.width(15); 
		    fileRunoff << noData << endl;
		}
	    }
	    else 
	    {
	        fileRunoff.width(15); 
		fileRunoff << noData << endl;
	    }
	}
    }
    fileRunoff.close();
  
    //  fout << "Change in storage except glacier storage for landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"storageChange");
    strcat(fileName,timeName);
    ofstream fileChange(fileName);
    if (!fileChange.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    fileChange.precision(0); 
    fileChange.setf(ios::fixed); 
    fileChange << "ncols         " << nCols << endl;
    fileChange << "nrows         " << nRows << endl;
    fileChange << "xllcorner     " << xllCorner << endl;
    fileChange << "yllcorner     " << yllCorner << endl;
    fileChange << "cellsize      " << cellSize << endl;
    fileChange << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
            if (k < numLand) 
            {
                if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) 
                {
                    if (Dew[k].GetSumRunoff() > missingData) 
                    {
                        fileChange.width(15); 
                        fileChange.precision(5); 
                        fileChange.setf(ios::showpoint); 
                        fileChange.setf(ios::fixed); 
                        fileChange << (Dew[k].GetFinalStorage()-Dew[k].GetInitialStorage())/Dew[k].GetNumberSum()*numberDaysPerYear*1000.0 << endl;
                    }
                    else 
                    {
                        fileChange.width(15); 
                        fileChange << noData << endl;
                    }
                    k++;
                }
                else 
                {
                    fileChange.width(15); fileChange << noData << endl;
                }
            }
            else 
            {
                fileChange.width(15); fileChange << noData << endl;
            }
	}
    }
    fileChange.close();
}
  
  
void WriteAsciiGridAnnualGlacierValues(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,
                                       int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
    int i, j, k;
    char fileName[100];
    char dirName[100];
    char timeName[100];

    //  sprintf(dirName,"./%04d/",datetime.getYear());
    sprintf(dirName, "./");
    sprintf(timeName, "_%04d_%02d_%02d.asc", datetime.getYear(), datetime.getMonth(), datetime.getDay());

    //  fout << "Glacier ice thickness in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "git");
    strcat(fileName, timeName);
    ofstream fileGlacIceThick(fileName);
    if (!fileGlacIceThick.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileGlacIceThick.precision(0); 
    fileGlacIceThick.setf(ios::fixed); 
    fileGlacIceThick << "ncols         " << nCols << endl;
    fileGlacIceThick << "nrows         " << nRows << endl;
    fileGlacIceThick << "xllcorner     " << xllCorner << endl;
    fileGlacIceThick << "yllcorner     " << yllCorner << endl;
    fileGlacIceThick << "cellsize      " << cellSize << endl;
    fileGlacIceThick << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetGlacierIceThickness() > missingData)
                    {
                        fileGlacIceThick.width(15);
                        fileGlacIceThick.precision(5);
                        fileGlacIceThick.setf(ios::showpoint);
                        fileGlacIceThick.setf(ios::fixed);
                        //            fileGlacIceThick << (Dew[k].GetGlacierIceThickness())*1000.0 << endl;
                        fileGlacIceThick << (Dew[k].GetGlacierIceThickness()) << endl;
                    }
                    else
                    {
                        fileGlacIceThick.width(15);
                        fileGlacIceThick << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileGlacIceThick.width(15);
                    fileGlacIceThick << noData << endl;
                }
            }
            else
            {
                fileGlacIceThick.width(15);
                fileGlacIceThick << noData << endl;
            }
        }
    }
    fileGlacIceThick.close();

    //  fout << "Glacier surface elevation in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "gse");
    strcat(fileName, timeName);
    ofstream fileGlacSurfElev(fileName);
    if (!fileGlacSurfElev.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileGlacSurfElev.precision(0); 
    fileGlacSurfElev.setf(ios::fixed); 
    fileGlacSurfElev << "ncols         " << nCols << endl;
    fileGlacSurfElev << "nrows         " << nRows << endl;
    fileGlacSurfElev << "xllcorner     " << xllCorner << endl;
    fileGlacSurfElev << "yllcorner     " << yllCorner << endl;
    fileGlacSurfElev << "cellsize      " << cellSize << endl;
    fileGlacSurfElev << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetGlacierSurfaceElevation() > missingData)
                    {
                        fileGlacSurfElev.width(15);
                        fileGlacSurfElev.precision(5);
                        fileGlacSurfElev.setf(ios::showpoint);
                        fileGlacSurfElev.setf(ios::fixed);
                        //      fileGlacSurfElev << (Dew[k].GetGlacierSurfaceElevation())*1000.0 << endl;
                        fileGlacSurfElev << (Dew[k].GetGlacierSurfaceElevation()) << endl;
                    }
                    else
                    {
                        fileGlacSurfElev.width(15);
                        fileGlacSurfElev << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileGlacSurfElev.width(15);
                    fileGlacSurfElev << noData << endl;
                }
            }
            else
            {
                fileGlacSurfElev.width(15);
                fileGlacSurfElev << noData << endl;
            }
        }
    }
    fileGlacSurfElev.close();

    //  fout << "Glacier annual mass balance in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "gam");
    strcat(fileName, timeName);
    ofstream fileAnnMassBal(fileName);
    if (!fileAnnMassBal.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileAnnMassBal.precision(0); 
    fileAnnMassBal.setf(ios::fixed); 
    fileAnnMassBal << "ncols         " << nCols << endl;
    fileAnnMassBal << "nrows         " << nRows << endl;
    fileAnnMassBal << "xllcorner     " << xllCorner << endl;
    fileAnnMassBal << "yllcorner     " << yllCorner << endl;
    fileAnnMassBal << "cellsize      " << cellSize << endl;
    fileAnnMassBal << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetAnnualMassBalance() > missingData)
                    {
                        fileAnnMassBal.width(15);
                        fileAnnMassBal.precision(5);
                        fileAnnMassBal.setf(ios::showpoint);
                        fileAnnMassBal.setf(ios::fixed);
                        //            fileAnnMassBal << (Dew[k].GetAnnualMassBalance())*1000.0 << endl;
                        fileAnnMassBal << (Dew[k].GetAnnualMassBalance()) << endl;
                    }
                    else
                    {
                        fileAnnMassBal.width(15);
                        fileAnnMassBal << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileAnnMassBal.width(15);
                    fileAnnMassBal << noData << endl;
                }
            }
            else
            {
                fileAnnMassBal.width(15);
                fileAnnMassBal << noData << endl;
            }
        }
    }
    fileAnnMassBal.close();

    //  fout << "Percentage of grid cells covered by glaciers in landscape elements at time step " << timeStep << ":\n";
    strcpy(fileName, dirName);
    strcat(fileName, "gpe");
    strcat(fileName, timeName);
    ofstream fileGlacPer(fileName);
    if (!fileGlacPer.is_open())
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit(1);
    }
    fileGlacPer.precision(0); 
    fileGlacPer.setf(ios::fixed); 
    fileGlacPer << "ncols         " << nCols << endl;
    fileGlacPer << "nrows         " << nRows << endl;
    fileGlacPer << "xllcorner     " << xllCorner << endl;
    fileGlacPer << "yllcorner     " << yllCorner << endl;
    fileGlacPer << "cellsize      " << cellSize << endl;
    fileGlacPer << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            if (k < numLand)
            {
                if (Dew[k].GetGeoIndex() == ELEMENT(i, j))
                {
                    if (Dew[k].GetGlacier())
                    {
                        fileGlacPer.width(15);
                        fileGlacPer.precision(5);
                        fileGlacPer.setf(ios::showpoint);
                        fileGlacPer.setf(ios::fixed);
                        fileGlacPer << (Dew[k].GetGlacier()->GetGlacierIceAreaFraction()) << endl;
                    }
                    else
                    {
                        fileGlacPer.width(15);
                        fileGlacPer.precision(5);
                        fileGlacPer.setf(ios::showpoint);
                        fileGlacPer.setf(ios::fixed);
                        fileGlacPer << 0.0 << endl;
                        //            fileGlacPer.width(15); fileGlacPer << noData << endl;
                    }
                    k++;
                }
                else
                {
                    fileGlacPer.width(15);
                    fileGlacPer << noData << endl;
                }
            }
            else
            {
                fileGlacPer.width(15);
                fileGlacPer << noData << endl;
            }
        }
    }
    fileGlacPer.close();
}


/*  Objective criteria of model performance  */
void ObjectiveCriteria(int numberTimeSteps, double *obs, double *sim,
                       double *ns, double *rmse, double *bias, double *pears)
{
    int i, sum_count;
    double obs_mean, sim_mean;
    double sum_1;
    double sum_2;
    double sum_3;

    *ns = missingData;            /*  Nash-Sutcliffe efficiency criterion  */
    *rmse = missingData;          /*  Root mean square error criterion  */
    *bias = missingData;          /*  Bias (volume error) criterion  */
    *pears = missingData;         /*  Pearson's product-moment correlation coefficient  */
    /*    for (i = 0; i < numberTimeSteps; i++) sim[i] = obs[i];*/

    /* Mean values */
    obs_mean = 0;
    sim_mean = 0;
    sum_count = 0;
    for (i = 0; i < numberTimeSteps; i++)
    {
        if (sim[i] > missingData && obs[i] > missingData)
        {
            obs_mean = obs_mean + obs[i];
            sim_mean = sim_mean + sim[i];
            sum_count++;
        }
    }
    if (sum_count > 0)
    {
        obs_mean = obs_mean / (double)sum_count;
        sim_mean = sim_mean / (double)sum_count;
        /* Nash-Sutcliffe efficiency criterion */
        sum_1 = 0;
        sum_2 = 0;
        for (i = 0; i < numberTimeSteps; i++)
        {
            if (sim[i] > missingData && obs[i] > missingData)
            {
                sum_1 = sum_1 + pow((obs[i] - sim[i]), 2);
                sum_2 = sum_2 + pow((obs[i] - obs_mean), 2);
            }
        }
        *ns = 1.0 - (sum_1 / sum_2);
        /* Root mean square error criterion */
        sum_1 = 0;
        for (i = 0; i < numberTimeSteps; i++)
        {
            if (sim[i] > missingData && obs[i] > missingData)
            {
                sum_1 = sum_1 + pow((obs[i] - sim[i]), 2);
            }
        }
        *rmse = pow((sum_1 / (double)i), 0.5);
        /* Bias (volume error) criterion */
        sum_1 = 0;
        sum_2 = 0;
        for (i = 0; i < numberTimeSteps; i++)
        {
            if (sim[i] > missingData && obs[i] > missingData)
            {
                sum_1 = sum_1 + (sim[i] - obs[i]);
                sum_2 = sum_2 + obs[i];
            }
        }
        *bias = sum_1 / sum_2;
        /* Pearson's product-moment correlation coefficient criterion */
        sum_1 = 0;
        sum_2 = 0;
        sum_3 = 0;
        for (i = 0; i < numberTimeSteps; i++)
        {
            if (sim[i] > missingData && obs[i] > missingData)
            {
                sum_1 = sum_1 + (obs[i] - obs_mean) * (sim[i] - sim_mean);
                sum_2 = sum_2 + pow((obs[i] - obs_mean), 2);
                sum_3 = sum_3 + pow((sim[i] - sim_mean), 2);
            }
        }
        *pears = pow((sum_1 / sqrt(sum_2 * sum_3)), 2);
    }
    return;
}


void WriteTotalReservoirStorage(TotalReservoirStorage * const TotalReservoirStore, DateTime startSimulationTime,
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                bool modelCalibration, int secondsPerTimeStep)
{
    FILE *fpOut;
    int timeStep;
    DateTime datetime;
    char fileName[100];

    sprintf(fileName, "res_totalvol.var");
    if ((fpOut = fopen(fileName, "w")) == NULL)
    {
        printf("\n File %s not found!\n\n", fileName);
        exit(1);
    }
    if (modelCalibration)
        fprintf(fpOut, "%29.6f\n", TotalReservoirStore->GetTotalReservoirStorage(initialTimeSteps + numberTimeSteps - 1) -
                TotalReservoirStore->GetTotalReservoirStorage(initialTimeSteps));
    timeStep = initialTimeSteps;
    for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
    {
        if (TotalReservoirStore->GetTotalReservoirStorage(timeStep) != missingData)
        {
            fprintf(fpOut, "%4d%02d%02d/%02d%02d%29.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    TotalReservoirStore->GetTotalReservoirStorage(timeStep));
        }
        else
        {
            fprintf(fpOut, "%4d%02d%02d/%02d%02d%29.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
        }
        timeStep++;
    }
    fclose(fpOut);
}


void WriteDistributedElementTimeSeries(DistributedElement * const Dew, int numLand,
                                       DateTime startSimulationTime, DateTime endSimulationTime,
                                       int initialTimeSteps, int numberTimeSteps, int secondsPerTimeStep)
{
    FILE *fpHgw, *fpKgw1, *fpKgw2;
    char fileName[100];
    int i, j, timeStep;
    DateTime datetime;

    for (i = 0; i < numLand; i++)
    {
        //  HBV state variables time series output
        for (j = 0; j < Dew[i].GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++)
        {
            if (Dew[i].GetLandIndex() == Dew[i].GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j))
            {
                sprintf(fileName, "HBV_groundwater_%d.var", Dew[i].GetLandIndex());
                if ((fpHgw = fopen(fileName, "w")) == NULL)
                {
                    printf("\n File %s not found!\n\n", fileName);
                    exit(1);
                }
                timeStep = initialTimeSteps;
                for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
                {
                    if (Dew[i].GetDistributedElementHbvUpperZone(timeStep) != missingData)
                    {
                        fprintf(fpHgw, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                Dew[i].GetSelectedHbvTimeSeriesElements()->GetGroundWaterRef(j) +
                                (Dew[i].GetDistributedElementHbvUpperZone(timeStep) + Dew[i].GetDistributedElementHbvLowerZone(timeStep)) /
                                Dew[i].GetSelectedHbvTimeSeriesElements()->GetEffPor(j));
                    }
                    else
                    {
                        fprintf(fpHgw, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    timeStep++;
                }
                fclose(fpHgw);
                if (timeStep != initialTimeSteps + numberTimeSteps)
                {
                    cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
                    exit(1);
                }
            }
        }
        //  KiWa state variables time series output
        for (j = 0; j < Dew[i].GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++)
        {
            if (Dew[i].GetLandIndex() == Dew[i].GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j))
            {
                sprintf(fileName, "KiWa_groundwater_%d_one.var", Dew[i].GetLandIndex());
                if ((fpKgw1 = fopen(fileName, "w")) == NULL)
                {
                    printf("\n File %s not found!\n\n", fileName);
                    exit(1);
                }
                timeStep = initialTimeSteps;
                for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
                {
                    if (Dew[i].GetDistributedElementKiWaGroundWaterDepthOne(timeStep) != missingData)
                    {
                        fprintf(fpKgw1, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                Dew[i].GetDistributedElementKiWaGroundWaterDepthOne(timeStep));
                    }
                    else
                    {
                        fprintf(fpKgw1, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    timeStep++;
                }
                fclose(fpKgw1);
                if (timeStep != initialTimeSteps + numberTimeSteps)
                {
                    cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
                    exit(1);
                }
                sprintf(fileName, "KiWa_groundwater_%d_two.var", Dew[i].GetLandIndex());
                if ((fpKgw2 = fopen(fileName, "w")) == NULL)
                {
                    printf("\n File %s not found!\n\n", fileName);
                    exit(1);
                }
                timeStep = initialTimeSteps;
                for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
                {
                    if (Dew[i].GetDistributedElementKiWaGroundWaterDepthTwo(timeStep) != missingData)
                    {
                        fprintf(fpKgw2, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                                Dew[i].GetDistributedElementKiWaGroundWaterDepthTwo(timeStep));
                    }
                    else
                    {
                        fprintf(fpKgw2, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
                    }
                    timeStep++;
                }
                fclose(fpKgw2);
                if (timeStep != initialTimeSteps + numberTimeSteps)
                {
                    cout << " timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps + numberTimeSteps << endl << endl;
                    exit(1);
                }
            }
        }
    }
}
