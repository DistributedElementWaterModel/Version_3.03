/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of river networks,          *
 *  lakes and landscape elements.                                                       *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "classRiverFlow.h"

#define ELEMENT(a,b) (((a)*nCols)+(b))


//class DistributedElement
DistributedElement::DistributedElement()
{
  SetNumUpLand(0);
  for (int i=0; i<7; i++) 
    upLandFlow[i]=0;
}

DistributedElement::~DistributedElement()
{ 
}


//class WaterCourseElement
WaterCourseElement::WaterCourseElement()
{
  SetNumUpStream(0);
  SetLandScapeElement(0); 
  SetLakeElementList(0);
  for (int i=0; i<7; i++) 
    upStreamElement[i]=0;
}

WaterCourseElement::~WaterCourseElement()
{
}
