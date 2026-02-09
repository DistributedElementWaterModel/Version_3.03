#pragma once
#include <list>

using namespace std;
class DateTime;
class SubCatchment;
class DistributedElement;

class GlacierElements
{
public:
    GlacierElements(SubCatchment ** OutletList, int numOut);
    ~GlacierElements();
    //  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
    //  ParametersGeneral *GetGeneralPar() const { return commonPar; }
    //  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
    //  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
    void BuildGlacierLists();
    void TraverseGlacierSubCatchment(SubCatchment * const thisSubCatchment, list <DistributedElement *> *glacierElementsList);
    void TraverseGlacierLandScape(DistributedElement * thisElement, DistributedElement ** glacierList, DistributedElement ** lastGlac);
    void SortGlacierLists();
    void SortAllGlacierLists(list <DistributedElement *> *glacierElementsList);
    void SortOneGlacierList(DistributedElement ** glacierList);
    void RemoveElementsGlacierLists();
    void RemoveElementsAllGlacierLists(list <DistributedElement *> *glacierElementsList);
    void RemoveElementsOneGlacierList(DistributedElement ** glacierList);
    void SetThisYearAnnualGlacierValues(DateTime datetime);
    void GlacierSurfaceElevationReDistribution();
    void ChangeGlacierSurfaceElevation();

private:
    int numWatcOut;
    SubCatchment ** Outlet;
    //  ParametersGeneral *commonPar;
    //  ParametersLandSurface *landSurfacePar;
    list <DistributedElement *> glacierElementsList;
};
