#pragma once
enum MODEL_STRUCTURE { HBV_MODEL, KWA_MODEL };
//enum LANDSURFACE { OPEN, BOG, FOREST, ALPINE, HEATHER, ROCK, GLACIER };
//enum SOIL { OPEN_SOIL, PEAT, FOREST_SOIL, ALPINE_SOIL, HEATHER_SOIL, BEDROCK, GLACIER_BED };
enum LANDSURFACE { SURF0, SURF1, SURF2, SURF3, SURF4, SURF5, SURF6, SURF7, SURF8, SURF9,  
                   SURF10, SURF11, SURF12, SURF13, SURF14, SURF15, SURF16, SURF17, SURF18, SURF19, 
                   GLACIER };
//enum SOIL { SOIL0, SOIL1, SOIL2, SOIL3, SOIL4, SOIL5, SOIL6, SOIL7, SOIL8, SOIL9, 
//            SOIL10, SOIL11, SOIL12, SOIL13, SOIL14, SOIL15, SOIL16, SOIL17, SOIL18, SOIL19, 
//            GLACIER_BED };
enum SOIL { SOIL0, SOIL1, SOIL2, SOIL3, SOIL4, SOIL5, GLACIER_BED }; //NB! Navn i hbv_soil_parameters.dta må stemme overens med disse 
enum GLACIER_TYPE { VALLEY, PLATEAU };
#define _USE_MATH_DEFINES   //pi
#define ELEMENT(a,b) (((a)*nCols)+(b))
const int numberPotentialEvaporationValuesPerYear=12;      // Number of long-term mean potential evaporation values per year
const int numberModelStructures = 2;           // All possible model structures
const int maximumNumberLandClasses = 3;        // Maximum number of land/soil classes in use for each computational element
//const int numberLandSurfaceClasses=7;        // All possible land surface types including glaciers, excluding lakes
const int numberSoilClasses=7;               // All possible soil/subsurface types including glaciers, excluding lakes
const int numberLandSurfaceClasses=21;       // All possible land surface types including glaciers, excluding lakes
//const int numberSoilClasses=21;              // All possible soil/subsurface types including glaciers, excluding lakes
const int numberGlacierClasses = 2;            // All possible glacier types
const int maximumNumberWaterVelocities = 1;
const int maximumNumberManningClasses = 1000;
const int maximumCorrectionCatchments = 1000;
const int numberInputSeries = 7;               // Number of meteorological input series (1 prec + 1 temp = 2)
const int numberSnowClasses = 9;
const int numberCharacteristic = 100;          // Number of characteristic curves in kinematic wave model of saturated subsurface flow
const int minimumTimeStep = 3600;
const double numberSecondsDay=86400.0;
const double missingData = -9999.0;
const double largeMissingData = -1.0e6;
const double probNorm[9] = { 0.01, 0.04, 0.1, 0.2, 0.3, 0.2, 0.1, 0.04, 0.01 };
const double epsilon = 1.0e-4;
//const double minimumGlacierIceThickness=10.0;
//const int numCols1km = 5;
//const int numRows1km = 8;
//const double xllCorner1km = 0.0;
//const double yllCorner1km = 0.0;
//const int cellSize1km = 500;
const int numCols1km = 1195;
const int numRows1km = 1550;
const double xllCorner1km = -75000.0;
const double yllCorner1km = 6450000.0;
const int cellSize1km = 1000;
const double maximumInterceptedWater = 0.2;    // Maximum intercepted water as fraction of leaf area index                          ParametersGeneral
const double vegetationDensity = 0.8;          // Horisontal projection of vegetation density fraction = 1 - sky view fraction      ParametersGeneral
