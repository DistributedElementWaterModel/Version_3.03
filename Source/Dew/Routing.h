#pragma once

#include "Util.h"
#include "stdafx.h"

class DateTime;
class SubCatchment;
class DistributedElement;

using namespace std;


void TraverseWaterCourseTravelTime(SubCatchment * const thisSubCatchment, double travelTime, double * waterVelocity, ofstream &fout);
void TraverseLandScapeTravelTime(DistributedElement * const thisElement, double travelTime, double * waterVelocity, ofstream &fout);
void RouteSourceToSinkDischarge(DistributedElement * Dew, int numLand, double * sourceToSinkDischarge, int initialTimeSteps, int numberTimeSteps, int timeStep, int secondsPerTimeStep);
void TraverseRouteWaterCourse(SubCatchment * const thisSubCatchment, int timeStep, int secondsPerTimeStep, ofstream &fout);
double RoutReferenceValues(double manning, double riverSlope, double width, double Q);
double RoutReferenceDischarge(double manning, double riverSlope, double area, double perimeter);
double RoutReferenceWaterLevel(double manning, double riverSlope, double width, double area, double perimeter);
double RoutOutputValues(double cRef, double epsilon, double secondsPerTimeStep, double spaceDist, double routeInput1, double routeInput2, double routeOutput1);
void WriteSourceToSinkDischarge(SubCatchment * const thisSubCatchment, double * sourceToSinkDischarge, DateTime startSimulationTime,
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps, int secondsPerTimeStep);
double MaxBasWeight(double, double, double);
double MaxBasFunction(double, double);
void TraverseMaxBas(SubCatchment * const thisSubCatchment, int initialTimeSteps, int numberTimeSteps, ofstream &fout);
void TransformMaxBas(SubCatchment * const thisSubCatchment, int initialTimeSteps, int numberTimeSteps, ofstream &fout);

