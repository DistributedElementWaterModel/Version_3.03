#include "GlacierElements.h"
#include "DistributedElement.h"
#include "Glacier.h"
#include "GlacierIce.h"
#include "ParametersGeneral.h"
#include "ParametersGlacierRetreat.h"
#include "SubCatchment.h"
#include "stdafx.h"
#include "DateTime.h"

GlacierElements::GlacierElements(SubCatchment ** OutletList, int numOut) :
    numWatcOut(numOut)
{
    Outlet = OutletList;
    //  SetGeneralPar(0);
    //  SetLandSurfacePar(0);
}

GlacierElements::~GlacierElements()
{
}

void GlacierElements::BuildGlacierLists()
{
    int i;
    for (i = 0; i < numWatcOut; i++)
    {
        TraverseGlacierSubCatchment(Outlet[i], &glacierElementsList);
    }
    //  Test
    /*  list <DistributedElement *>::iterator listIterator;
    for(listIterator=glacierElementsList.begin(); listIterator!=glacierElementsList.end(); listIterator++) {
    if (*listIterator) cout << "BuildGlacierLists (*listIterator) " << (*listIterator) << endl;
    }*/
    //  Test end
    //  SortAllGlacierLists(&glacierElementsList);
}

void GlacierElements::SortGlacierLists()
{
    SortAllGlacierLists(&glacierElementsList);
}

void GlacierElements::RemoveElementsGlacierLists()
{
    RemoveElementsAllGlacierLists(&glacierElementsList);
}


void GlacierElements::TraverseGlacierSubCatchment(SubCatchment * const thisSubCatchment, list <DistributedElement *> *glacierElementsList)
{
    int i;
    bool first = true;
    DistributedElement *thisElement, *lastDew = 0, *distributedElementList = 0;
    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        TraverseGlacierSubCatchment(thisSubCatchment->GetUpStream(i), glacierElementsList);
    }
    //cout << "\n* SubCatchment " << thisSubCatchment->GetSubCatchmentIndex() << "  " << thisSubCatchment->GetIdentifier() << endl;
    thisElement = thisSubCatchment->GetLandScapeElement();
    while (thisElement)
    {
        TraverseGlacierLandScape(thisElement, &distributedElementList, &lastDew);
        thisElement = thisElement->GetNextElement();
    }
    // Sort and store list of glacier elements
    if (distributedElementList)
    {
        //    cout << "3 distributedElementList " << distributedElementList << "  " << distributedElementList->GetGeoIndex() << endl;
        SortOneGlacierList(&distributedElementList);
        (*glacierElementsList).push_back(distributedElementList);
    }
    //  Test
    /*  list <DistributedElement *>::iterator listIterator;
    for(listIterator=(*glacierElementsList).begin(); listIterator!=(*glacierElementsList).end(); listIterator++) {
    if (*listIterator) cout << "TraverseGlacierSubCatchment (*listIterator) " << (*listIterator) << endl;
    }*/
    //  Test end
}


void GlacierElements::TraverseGlacierLandScape(DistributedElement * thisElement, DistributedElement ** distributedElementList, DistributedElement ** lastDew)
{
    int i;
    for (i = 0; i < thisElement->GetNumUpLand(); i++)
    {
        TraverseGlacierLandScape(thisElement->GetUpLandFlow(i), distributedElementList, lastDew);
    }
    // Build list of glacier elements
    if (thisElement->GetGlacier())
    {
        //    cout << "1            thisElement " << thisElement << "  " << thisElement->GetGeoIndex() << endl;
        if (!*distributedElementList)
        {
            *distributedElementList = thisElement;
            *lastDew = thisElement;
            //      cout << "dist 1 " << (*distributedElementList)->GetGeoIndex() << endl;
            //      cout << "last 1 " << (*lastDew)->GetGeoIndex() << endl;
        }
        else
        {
            (*lastDew)->SetNextGlacierElement(thisElement);
            *lastDew = thisElement;
            //      cout << "dist 2 " << (*distributedElementList)->GetGeoIndex() << endl;
            //      cout << "last 2 " << (*lastDew)->GetGeoIndex() << endl;
        }
        //    cout << "2 distributedElementList " << *distributedElementList << "  " << (*distributedElementList)->GetGeoIndex() << endl;
    }
}


void GlacierElements::SortAllGlacierLists(list <DistributedElement *> *glacierElementsList)
{
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = (*glacierElementsList).begin(); listIterator != (*glacierElementsList).end(); listIterator++)
    {
        //  Test
        //  if (*listIterator) cout << "SortAllGlacierLists (*listIterator) " << (*listIterator) << endl;
        //  Test end
        if (*listIterator)
        {
            SortOneGlacierList(&(*listIterator));
        }
    }
}


// Sort list of landscape elements with glaciers
void GlacierElements::SortOneGlacierList(DistributedElement ** distributedElementList)
{
    int i;
    bool first = true;
    DistributedElement *thisDew, *endDew, *tempDew, *thisNext, *endNext, *thisPrevious, *endPrevious;
    //  Test
    /*  thisDew=*distributedElementList;
    cout << "\n";
    while (thisDew) {
    cout << "1 SortOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << endl;
    thisDew=thisDew->GetNextGlacierElement();
    }*/
    //  Test end
    thisDew = *distributedElementList;
    thisPrevious = 0;
    while (thisDew->GetNextGlacierElement())
    {
        endPrevious = thisDew;
        endDew = thisDew->GetNextGlacierElement();
        do
        {
            if (endDew->GetGlacier()->GetGlacierSurfaceElevation() < thisDew->GetGlacier()->GetGlacierSurfaceElevation())
            {
                thisNext = thisDew->GetNextGlacierElement();
                endNext = endDew->GetNextGlacierElement();
                tempDew = thisDew;
                thisDew = endDew;
                endDew = tempDew;
                if (thisDew == thisNext)
                {
                    thisDew->SetNextGlacierElement(endDew);
                }
                else
                {
                    thisDew->SetNextGlacierElement(thisNext);
                    endPrevious->SetNextGlacierElement(endDew);
                }
                endDew->SetNextGlacierElement(endNext);
                if (thisPrevious)
                {
                    thisPrevious->SetNextGlacierElement(thisDew);
                }
            }
            endPrevious = endDew;
            endDew = endDew->GetNextGlacierElement();
        }
        while (endDew);
        if (first)
        {
            *distributedElementList = thisDew;
        }
        thisPrevious = thisDew;
        thisDew = thisDew->GetNextGlacierElement();
        first = false;
    }
    //  Test
    /*  thisDew=*distributedElementList;
    cout << "\n";
    while (thisDew) {
    cout << "2 SortOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << endl;
    thisDew=thisDew->GetNextGlacierElement();
    }
    cout << "\n";*/
    //  Test end
}


void GlacierElements::RemoveElementsAllGlacierLists(list <DistributedElement *> *glacierElementsList)
{
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = (*glacierElementsList).begin(); listIterator != (*glacierElementsList).end(); listIterator++)
    {
        //    if (*listIterator) cout << "\n1 RemoveElementsAllGlacierLists (*listIterator) " << (*listIterator) << endl;
        if (*listIterator)
        {
            RemoveElementsOneGlacierList(&(*listIterator));
        }
    }
    //  Test
    /*  DistributedElement *thisDew;
    list <DistributedElement *>::iterator testIterator;
    for (testIterator=(*glacierElementsList).begin(); testIterator!=(*glacierElementsList).end(); testIterator++) {
    thisDew = (*testIterator);
    if (*testIterator) cout << "2 RemoveElementsAllGlacierLists (*testIterator) " << (*testIterator) << endl;
    while (thisDew) {
    cout << "3 RemoveElementsAllGlacierLists thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << endl;
    thisDew=thisDew->GetNextGlacierElement();
    }
    }*/
    //  Test end
}


// Remove ice free elements from list of landscape elements with glaciers
void GlacierElements::RemoveElementsOneGlacierList(DistributedElement ** distributedElementList)
{
    int i;
    DistributedElement *thisDew, *thisPrevious;
    thisDew = *distributedElementList;
    *distributedElementList = 0;
    thisPrevious = 0;
    while (thisDew)
    {
        //    cout << "1 RemoveElementsOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << endl;
        if (thisDew->GetGlacier()->GetGlacierIceThickness() > 0.0)
        {
            //      cout << "2 RemoveElementsOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << endl;
            if (thisPrevious)
            {
                thisPrevious->SetNextGlacierElement(thisDew);
            }
            else
            {
                *distributedElementList = thisDew;
            }
            thisPrevious = thisDew;
        }
        else
        {
            thisDew->GetGlacier()->SetGlacierIceAreaFraction(0.0);
        }
        thisDew = thisDew->GetNextGlacierElement();
    }
    //  Test
    /*  thisDew=*distributedElementList;
    while (thisDew) {
    cout << "3 RemoveElementsOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << endl;
    thisDew=thisDew->GetNextGlacierElement();
    }*/
    //  Test end
}


void GlacierElements::SetThisYearAnnualGlacierValues(DateTime datetime)
{
    DistributedElement *thisDew;
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = glacierElementsList.begin(); listIterator != glacierElementsList.end(); listIterator++)
    {
        thisDew = (*listIterator);
        while (thisDew)
        {
            //      cout << "1 SetThisYearAnnualGlacierValues thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceAreaFraction() << endl;
            thisDew->SetThisYearAnnualGlacierValues(datetime);
            thisDew = thisDew->GetNextGlacierElement();
        }
    }
    //cout << "\n";
}


/*  Algorithm for redistribution of glacier surface based on Huss et al, Hydrology and Earth System Sciences, 2010  */
void GlacierElements::GlacierSurfaceElevationReDistribution()
{
    int numberGlacierElements;
    double a, b, c, gamma, scaleFactor;
    double elevationMin, elevationMax, elevationNormalized, elevationChangeNormalized;
    double elevationNew, iceThicknessNew;
    double densityIce, massBalanceSum, massBalanceSum2, meanMassBalance, areaSum, areaElevationSum;
    DistributedElement *thisDew;
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = glacierElementsList.begin(); listIterator != glacierElementsList.end(); listIterator++)
    {
        thisDew = (*listIterator);
        if (thisDew)
        {
            densityIce = thisDew->GetGeneralPar()->GetDENSITY_ICE();
            a = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetA();
            b = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetB();
            c = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetC();
            gamma = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetGAMMA();
            elevationMin = (-1) * largeMissingData;
            elevationMax = largeMissingData;
            // Huss et al, eq. 1
            while (thisDew)
            {
                if (thisDew->GetGlacier()->GetGlacierSurfaceElevation() < elevationMin)
                {
                    elevationMin = thisDew->GetGlacier()->GetGlacierSurfaceElevation();
                }
                if (thisDew->GetGlacier()->GetGlacierSurfaceElevation() > elevationMax)
                {
                    elevationMax = thisDew->GetGlacier()->GetGlacierSurfaceElevation();
                }
                thisDew = thisDew->GetNextGlacierElement();
            }
            massBalanceSum = 0.0;
            areaSum = 0.0;
            areaElevationSum = 0.0;
            numberGlacierElements = 0;
            thisDew = (*listIterator);
            while (thisDew)
            {
                numberGlacierElements++;
                if (elevationMax != elevationMin)
                {
                    elevationNormalized = (elevationMax - thisDew->GetGlacier()->GetGlacierSurfaceElevation()) / (elevationMax - elevationMin);
                    elevationChangeNormalized = pow((elevationNormalized + a), gamma) + b * (elevationNormalized + a) + c;
                }
                else
                {
                    elevationChangeNormalized = thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized();
                }
                thisDew->GetGlacier()->SetGlacierSurfaceElevationChangeNormalized(elevationChangeNormalized);
                massBalanceSum = massBalanceSum + thisDew->GetAnnualMassBalance() * thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                areaSum = areaSum + thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                areaElevationSum = areaElevationSum + elevationChangeNormalized * thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                thisDew = thisDew->GetNextGlacierElement();
            }
            if (numberGlacierElements <= 0)
            {
                cout << " Number of glacier elements error:  numberGlacierElements = " << numberGlacierElements << endl;
                exit(1);
            }
            if (massBalanceSum <= 0.0)
            {
                // Huss et al, eq. 2
                if (areaElevationSum != 0.0)
                {
                    scaleFactor = massBalanceSum / (areaElevationSum * densityIce);
                }
                else
                {
                    scaleFactor = 0.0;
                }
                // Huss et al, eq. 3
                massBalanceSum2 = 0.0;
                thisDew = (*listIterator);
                while (thisDew)
                {
                    elevationNew = thisDew->GetGlacier()->GetGlacierSurfaceElevation() + scaleFactor * thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized();
                    iceThicknessNew = thisDew->GetGlacier()->GetGlacierIceThickness() + scaleFactor * thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized();
                    if (iceThicknessNew < 0.0)
                    {
                        elevationNew = elevationNew - iceThicknessNew;
                        iceThicknessNew = 0.0;
                    }
                    thisDew->GetGlacier()->SetGlacierSurfaceElevation(elevationNew);
                    thisDew->GetGlacier()->SetGlacierIceThickness(iceThicknessNew);
                    thisDew->GetGlacier()->SetGlacierIceVolume(iceThicknessNew * thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0);
                    massBalanceSum2 = massBalanceSum2 + thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized() * thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                    thisDew = thisDew->GetNextGlacierElement();
                }
                // Control mass balance sum
                if (elevationMax != elevationMin)
                {
                    massBalanceSum2 = scaleFactor * densityIce * massBalanceSum2;
                    if (massBalanceSum < massBalanceSum2 - epsilon || massBalanceSum > massBalanceSum2 + epsilon)
                    {
                        cout << " Mass balance error:  mass balance sum initial = " << massBalanceSum << " mass balance sum final = " << massBalanceSum2 << endl;
                        exit(1);
                    }
                }
            }
            else
            {
                meanMassBalance = massBalanceSum / areaSum;
                thisDew = (*listIterator);
                while (thisDew)
                {
                    elevationNew = thisDew->GetGlacier()->GetGlacierSurfaceElevation() + meanMassBalance / densityIce;
                    iceThicknessNew = thisDew->GetGlacier()->GetGlacierIceThickness() + meanMassBalance / densityIce;
                    if (iceThicknessNew < 0.0)
                    {
                        elevationNew = elevationNew - iceThicknessNew;
                        iceThicknessNew = 0.0;
                    }
                    thisDew->GetGlacier()->SetGlacierSurfaceElevation(elevationNew);
                    thisDew->GetGlacier()->SetGlacierIceThickness(iceThicknessNew);
                    thisDew->GetGlacier()->SetGlacierIceVolume(iceThicknessNew * thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0);
                    thisDew = thisDew->GetNextGlacierElement();
                }
            }
        }
    }
}


/*  Algorithm for changing elevation of glacier surface  */
void GlacierElements::ChangeGlacierSurfaceElevation()
{
    int numberGlacierElements;
    double elevationNew, iceThicknessNew;
    double densityIce, massBalance;
    DistributedElement *thisDew;
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = glacierElementsList.begin(); listIterator != glacierElementsList.end(); listIterator++)
    {
        thisDew = (*listIterator);
        if (thisDew)
        {
	    densityIce = thisDew->GetGeneralPar()->GetDENSITY_ICE();
            thisDew = (*listIterator);
            while (thisDew)
            {
		elevationNew = thisDew->GetGlacier()->GetGlacierSurfaceElevation() + thisDew->GetAnnualMassBalance() / densityIce;
		iceThicknessNew = thisDew->GetGlacier()->GetGlacierIceThickness() + thisDew->GetAnnualMassBalance() / densityIce;
		if (iceThicknessNew < 0.0)
		{
		    elevationNew = elevationNew - iceThicknessNew;
		    iceThicknessNew = 0.0;
		}
		thisDew->GetGlacier()->SetGlacierSurfaceElevation(elevationNew);
		thisDew->GetGlacier()->SetGlacierIceThickness(iceThicknessNew);
		thisDew->GetGlacier()->SetGlacierIceVolume(iceThicknessNew * thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0);
                thisDew = thisDew->GetNextGlacierElement();
            }
	}
    }
}
