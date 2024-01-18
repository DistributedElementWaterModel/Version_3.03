/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
 *  Preprocessing.                                                                      *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "classPreDew.h"

void ReadLandUse(DistributedElement * const Dew, int nRows, int nCols, int
                 noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void ReadLandUseGeneral(DistributedElement * const Dew, int nRows, int nCols, int
                 noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void WriteLandScapeElements(DistributedElement * const Dew, ParametersGeneral * const ParGeneralStore,
                            MeteorologicalStations * const MetStations, int landIndex, 
                            int preSetModelStructure, int nRows, int nCols, 
                            int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void DistanceSort(DistributedElement  * const distElement, ParametersGeneral * const ParGeneralStore, 
                  MeteorologicalStations * const MetStations,
                  int * precStations, int * tempStations, double * precWeights, double * tempWeights);
void GetRowCol(int elementNo, int nCols, int &row, int &col); 
void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout);
void ReadFlowDirection(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
int ControlFlowDirection(int, int);
void FindSubCatchmentIdentifier(DistributedElement * const Dew, SubCatchment * const CatchmentElement, int numWatc, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void WriteSubCatchmentIdentifier(DistributedElement * const Dew, SubCatchment * const CatchmentElement, int numWatc, ofstream &fout);
void FindUpLandFlow(DistributedElement * const Dew, int nRows, int nCols, int noData);
void WriteLandScapeAndSubCatchmentFlow(DistributedElement * const Dew, int nRows, int nCols, int noData, ofstream &fout);

int main(int argc, char *argv[])
{
  std::cout << "\n\n Preprocessing for Distributed Element Water Model \n\n";

  char fileName[80];
  char buffer[1024];
  char hierarchy = 'H';
  char ch;
  int i,j,k;
  int landIndex;
  int nRows, nCols, noData;
  int numWatc, numWatcUp, numWatcOut;
  int subCatchmentId;
  int preSetModelStructure = -1;
  int lakeNumber;
  double value;
  double correction, maxBas;
  double xllCorner, yllCorner, cellSize;

  if (argc != 2) {
    cout << " " << argv[0] << "  <control file name>\n\n";
    exit(1);
  }
  ifstream fileControl(argv[1]);
  if (!fileControl.is_open()) {
    cout << " Error opening file " << argv[1] << endl << endl;
    exit (1);
  }

  /*  while (preSetModelStructure < 0 || preSetModelStructure >= numberModelStructures) {
      cout << " Model structure, HBV (0) or KinematicWave (1): ";
      cin >> preSetModelStructure;
      cout << endl;
      }*/   
  fileControl.ignore(100,':');
  fileControl >> preSetModelStructure;
  fileControl.ignore(1024,'\n');
  if (preSetModelStructure < 0 || preSetModelStructure >= numberModelStructures) {
    cout << "\n Model structure, HBV (0) or KinematicWave (1): \n\n";
    exit (1);
  }

  /*  while (hierarchy != 'N' && hierarchy != 'C' && hierarchy != 'n' && hierarchy != 'c') {
      cout << " Landscape elements hierarchy, flow direction network(N) or nested catchments(C): ";
      cin >> hierarchy;
      cout << endl;
      }*/   
  fileControl.ignore(100,':');
  fileControl >> hierarchy;
  fileControl.ignore(1024,'\n');
  if (hierarchy != 'N' && hierarchy != 'C' && hierarchy != 'n' && hierarchy != 'c') {
    cout << "\n Landscape elements hierarchy, flow direction network(N) or nested catchments(C) \n\n";
    exit (1);
  }

  /*  cout << " Output file: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ofstream fout(fileName);  // Open for writing

  // Object for storing meteorological station information
  MeteorologicalStations * MetStations = new MeteorologicalStations;
  MetStations->SetMeteorologicalStations(fileControl, fout);

  // Read common parameters file and set parameter values
  ParametersGeneral * ParGeneralStore = new ParametersGeneral;
  SetGeneralParameters(ParGeneralStore, fileControl, fout);
  // End read common parameters

  // Read landscape element file and generate landscape element objects  
  /*  cout << " File with geographical analysis area: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGeo(fileName);  // Open for reading
  if (!finGeo.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  } 
  finGeo >> buffer >> nCols;
  finGeo >> buffer >> nRows;
  finGeo >> buffer >> xllCorner;
  finGeo >> buffer >> yllCorner;
  finGeo >> buffer >> cellSize;
  finGeo >> buffer >> noData;
  cout << endl << nCols << endl;
  cout << nRows << endl;
  cout << xllCorner << endl;
  cout << yllCorner << endl;
  cout << cellSize << endl;
  cout << buffer << "    " << noData << endl;
  landIndex = 0;
  DistributedElement * Dew = new DistributedElement [nRows*nCols];
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      Dew[ELEMENT(i,j)].SetGeoIndex(ELEMENT(i,j));
      finGeo >> value;
      if (value!=noData) {
        Dew[ELEMENT(i,j)].SetLandIndex(landIndex++);
      } else {
        Dew[ELEMENT(i,j)].SetLandIndex(noData);
      }
      Dew[ELEMENT(i,j)].SetXCoord(xllCorner + (j+0.5)*cellSize);
      Dew[ELEMENT(i,j)].SetYCoord(yllCorner + (nRows-i-0.5)*cellSize);
    }
  }
  finGeo.close();
  // Read end

  // Read land use for landscape elements
  //  ReadLandUse(Dew, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);
  ReadLandUseGeneral(Dew, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);

  // Open file with sub-catchment hierarchy
  /*  cout << " File with sub-catchment hierarchy: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream fileWCo(fileName);
  if (!fileWCo.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  
  // Read sub-catchment information and generate sub-catchment objects  
  fileWCo.ignore(100,':');
  fileWCo >> numWatc;
  cout << "\n # Number of sub-catchments " << numWatc << endl;
  SubCatchment *CatchmentElement = new SubCatchment [numWatc];
  for (i=0; i<numWatc; i++) {
    fileWCo >> j >> ch >> subCatchmentId >> lakeNumber >> maxBas >> correction;
    if (j != i) {
      cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << j << endl;
      exit (1);
    }
    fileWCo.ignore(1024,'\n');
    CatchmentElement[i].SetSubCatchmentIndex(i);
    CatchmentElement[i].SetIdentifier(subCatchmentId);
    CatchmentElement[i].SetCorrection(correction);
    cout << " Sub-catchment index " << i << "  " << "Sub-catchment identifier " << CatchmentElement[i].GetIdentifier() << endl;
  }

  // Watercourse/sub-catchment outlets
  fileWCo.ignore(100,':');
  fileWCo >> numWatcOut;
  cout << "\n # Number of watercourse/sub-catchment outlets " << numWatcOut << endl;
  SubCatchment ** Outlet = new SubCatchment * [numWatcOut];
  for (i=0; i<numWatcOut; i++) {
    fileWCo >> j;
    Outlet[i] = &CatchmentElement[j];
    cout << " Outlet no. " << i << "\t" << " Sub-catchment no. " << j << "\t" << endl;
    fileWCo.ignore(1024,'\n');
  }

  // Hierarchy of sub-catchments
  fileWCo.getline(buffer, 1024);
  cout << "\n " << buffer << endl;
  while (fileWCo >> i) {
    fileWCo >> numWatcUp;
    CatchmentElement[i].SetNumUpStream(numWatcUp);
    fileWCo.ignore(100,':');
    cout << " Downstream, sub-catchment no.  " << i << "    Identifier  " << CatchmentElement[i].GetIdentifier() << endl; 
    cout << " No. of upstream sub-catchments " << numWatcUp << endl;
    k = 0;
    while (fileWCo.peek() != '\n') {
      fileWCo >> j;
      cout << "\t" << "Upstream, sub-catchment no. " << j ;
      CatchmentElement[i].SetUpStream(k, &CatchmentElement[j]);
      cout  << "\t" << "UpStream[" << k << "]" << "    Identifier  " << CatchmentElement[i].GetUpStream(k)->GetIdentifier() << endl;
      while (fileWCo.peek() == ' ') fileWCo.ignore(1,' ');
      k++;
    }
    fileWCo.ignore(1024,'\n');
    if (numWatcUp!=k) {
      cout << endl << " Error in number of upstream pointers for sub-catchment no. " << i << endl << endl;
      exit (1);
    } 
  }
  fileWCo.close();
  // Read end

  // Read flow direction for landscape elements
  ReadFlowDirection(Dew, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);
    
  // Write landscape element information to input file to dew
  WriteLandScapeElements(Dew, ParGeneralStore, MetStations, landIndex, preSetModelStructure,  
                         nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
  
  // Find watercourse/sub-catchment identifiers for landscape elements and 
  // connect landscape elements to watercourse elements
  FindSubCatchmentIdentifier(Dew, CatchmentElement, numWatc, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout); 
    
  // Write information about watercourse/sub-catchment elements and landscape elements to input file to dew.cpp
  WriteSubCatchmentIdentifier(Dew, CatchmentElement, numWatc, fout);
    
  if (hierarchy == 'N' || hierarchy == 'n') {
    // Find flow direction for landscape elements
    FindUpLandFlow(Dew, nRows, nCols, noData);
    
    // Write flow direction for landscape elements
    WriteLandScapeAndSubCatchmentFlow(Dew, nRows, nCols, noData, fout);
  } 
  
  delete [] Outlet;
  delete [] CatchmentElement;
  
  fout.close();
  fileControl.close();
  delete [] Dew;
  return 0;
}


// Land surface classes based on potential tree level
void ReadLandUse(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[1024];
  int i,j,k;
  int nRo, nCo, noDa;
  double value, notAssigned, previousValue;
  double xllC, yllC, cellS;

  // Area of grid cells 
  for (i=0; i<nRows*nCols; i++) {
    Dew[i].SetArea(cellSize*cellSize);
  }

  // Read elevation of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with grid cell elevations: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finElevation(fileName);  // Open for reading
  if (!finElevation.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finElevation >> buffer >> nCo;
  finElevation >> buffer >> nRo;
  finElevation >> buffer >> xllC;
  finElevation >> buffer >> yllC;
  finElevation >> buffer >> cellS;
  finElevation >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for elevation!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finElevation >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetElevation(value);
    //    previousValue=value;
  }
  finElevation.close();
  // Read end

  // Read slope length of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope lengths: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finSlopeLength(fileName);  // Open for reading
  if (!finSlopeLength.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlopeLength >> buffer >> nCo;
  finSlopeLength >> buffer >> nRo;
  finSlopeLength >> buffer >> xllC;
  finSlopeLength >> buffer >> yllC;
  finSlopeLength >> buffer >> cellS;
  finSlopeLength >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope length!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlopeLength >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetSlopeLength(value);
    //    previousValue=value;
  }
  finSlopeLength.close();
  // Read end

  // Read slope angle of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope angles: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finSlopeAngle(fileName);  // Open for reading
  if (!finSlopeAngle.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlopeAngle >> buffer >> nCo;
  finSlopeAngle >> buffer >> nRo;
  finSlopeAngle >> buffer >> xllC;
  finSlopeAngle >> buffer >> yllC;
  finSlopeAngle >> buffer >> cellS;
  finSlopeAngle >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope angle!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlopeAngle >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetSlopeAngle(value);
    //    previousValue=value;
  }
  finSlopeAngle.close();
  // Read end

  // Read slope aspect of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope aspects: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finAspect(fileName);  // Open for reading
  if (!finAspect.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finAspect >> buffer >> nCo;
  finAspect >> buffer >> nRo;
  finAspect >> buffer >> xllC;
  finAspect >> buffer >> yllC;
  finAspect >> buffer >> cellS;
  finAspect >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for aspect!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finAspect >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetAspect(value);
    //    previousValue=value;
  }
  finAspect.close();
  // Read end

  // Read percentage of grid cells covered by lakes
  //  previousValue=-9999.0;
  /*cout << " File with lake percentages: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finLakePercent(fileName);  // Open for reading
  if (!finLakePercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finLakePercent >> buffer >> nCo;
  finLakePercent >> buffer >> nRo;
  finLakePercent >> buffer >> xllC;
  finLakePercent >> buffer >> yllC;
  finLakePercent >> buffer >> cellS;
  finLakePercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for lake percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finLakePercent >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetLakePercent(value);
    //    previousValue=value;
  }
  finLakePercent.close();
  // Read end

  // Read percentage of grid cells covered by forest
  //  previousValue=-9999.0;
  /*cout << " File with forest percentages: ";
    cin >> fileName;
    cout << endl; */
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finForestPercent(fileName);  // Open for reading
  if (!finForestPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finForestPercent >> buffer >> nCo;
  finForestPercent >> buffer >> nRo;
  finForestPercent >> buffer >> xllC;
  finForestPercent >> buffer >> yllC;
  finForestPercent >> buffer >> cellS;
  finForestPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for forest percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finForestPercent >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetForestPercent(value);
    //    previousValue=value;
  }
  finForestPercent.close();
  // Read end

  // Read percentage of grid cells covered by bogs
  //  previousValue=-9999.0;
  /*cout << " File with bog percentages: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finBogPercent(fileName);  // Open for reading
  if (!finBogPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finBogPercent >> buffer >> nCo;
  finBogPercent >> buffer >> nRo;
  finBogPercent >> buffer >> xllC;
  finBogPercent >> buffer >> yllC;
  finBogPercent >> buffer >> cellS;
  finBogPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for bog percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finBogPercent >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetBogPercent(value);
    //    previousValue=value;
  }
  finBogPercent.close();
  // Read end

  // Read percentage of grid cells covered by glaciers
  //  previousValue=-9999.0;
  /*cout << " File with glacier percentages: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGlacierPercent(fileName);  // Open for reading
  if (!finGlacierPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierPercent >> buffer >> nCo;
  finGlacierPercent >> buffer >> nRo;
  finGlacierPercent >> buffer >> xllC;
  finGlacierPercent >> buffer >> yllC;
  finGlacierPercent >> buffer >> cellS;
  finGlacierPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierPercent >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetGlacierPercent(value);
    //    previousValue=value;
  }
  finGlacierPercent.close();
  // Read end

  // Read elevation of glacier surface
  //  previousValue=-9999.0;
  /*cout << " File with glacier surface elevations: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGlacierSurfaceElevation(fileName);  // Open for reading
  if (!finGlacierSurfaceElevation.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierSurfaceElevation >> buffer >> nCo;
  finGlacierSurfaceElevation >> buffer >> nRo;
  finGlacierSurfaceElevation >> buffer >> xllC;
  finGlacierSurfaceElevation >> buffer >> yllC;
  finGlacierSurfaceElevation >> buffer >> cellS;
  finGlacierSurfaceElevation >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier surface elevation!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierSurfaceElevation >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetGlacierSurfaceElevation(value);
    //    previousValue=value;
  }
  finGlacierSurfaceElevation.close();
  // Read end

  // Read thickness of glacier ice
  //  previousValue=-9999.0;
  /*cout << " File with glacier ice thickness: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGlacierIceThickness(fileName);  // Open for reading
  if (!finGlacierIceThickness.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierIceThickness >> buffer >> nCo;
  finGlacierIceThickness >> buffer >> nRo;
  finGlacierIceThickness >> buffer >> xllC;
  finGlacierIceThickness >> buffer >> yllC;
  finGlacierIceThickness >> buffer >> cellS;
  finGlacierIceThickness >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier ice thickness!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierIceThickness >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetGlacierIceThickness(value);
    //    previousValue=value;
  }
  finGlacierIceThickness.close();
  // Read end

  // Read tree level for grid cells
  //  previousValue=-9999.0;
  /*cout << " File with tree levels: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finTreeLevel(fileName);  // Open for reading
  if (!finTreeLevel.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finTreeLevel >> buffer >> nCo;
  finTreeLevel >> buffer >> nRo;
  finTreeLevel >> buffer >> xllC;
  finTreeLevel >> buffer >> yllC;
  finTreeLevel >> buffer >> cellS;
  finTreeLevel >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for tree levels!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finTreeLevel >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetTreeLevel(value);
    //    previousValue=value;
  }
  finTreeLevel.close();
  // Read end

  // Assign land use classes based on potential tree level
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetLandIndex()!=noData) { 
      notAssigned = 100.0 - Dew[i].GetLakePercent() - Dew[i].GetGlacierPercent() -
          Dew[i].GetForestPercent() - Dew[i].GetBogPercent();
      if (notAssigned < 0.0) notAssigned = 0.0;
      if (Dew[i].GetElevation() > Dew[i].GetTreeLevel()) {
        if (Dew[i].GetForestPercent() > 0.0) {
          if (notAssigned > 0.0) {
            Dew[i].SetAlpinePercent(Dew[i].GetForestPercent() + notAssigned*0.2);
            Dew[i].SetHeatherPercent(notAssigned*0.8);
            Dew[i].SetForestPercent(0.0);
          } 
          else {
            Dew[i].SetAlpinePercent(Dew[i].GetForestPercent());
            Dew[i].SetForestPercent(0.0);
          }
        }
        else if (Dew[i].GetElevation() <= Dew[i].GetTreeLevel()+100.0) {
          Dew[i].SetAlpinePercent(notAssigned*0.2);
          Dew[i].SetHeatherPercent(notAssigned*0.8);
        }
        else if (Dew[i].GetElevation() <= Dew[i].GetTreeLevel()+200.0) {
          Dew[i].SetHeatherPercent(notAssigned*0.5);
          Dew[i].SetBedrockPercent(notAssigned*0.5);
        }
        else { 
          Dew[i].SetBedrockPercent(notAssigned);
        }
      }
      else if (Dew[i].GetElevation() > Dew[i].GetTreeLevel()*0.9) {
        if (notAssigned > 0.0) {
          Dew[i].SetAlpinePercent(Dew[i].GetForestPercent() + notAssigned*0.2);
          Dew[i].SetHeatherPercent(notAssigned*0.8);
          Dew[i].SetForestPercent(0.0);
        }
        else {
          Dew[i].SetAlpinePercent(Dew[i].GetForestPercent()*0.5);
          Dew[i].SetForestPercent(Dew[i].GetForestPercent()*0.5);
        }
      }
      else if (Dew[i].GetElevation() > Dew[i].GetTreeLevel()*0.8) {
        if (notAssigned > 0.0) {
          Dew[i].SetOpenLandPercent(notAssigned);
        }
        else {
          Dew[i].SetAlpinePercent(Dew[i].GetForestPercent()*0.2);
          Dew[i].SetForestPercent(Dew[i].GetForestPercent()*0.8);
        }
      }
      else {
        if (notAssigned > 0.0) {
          Dew[i].SetOpenLandPercent(notAssigned);
        }
      }
    }
  }

  // Control and correct area fractions of land use classes
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetLandIndex()!=noData) { 
      value=Dew[i].GetForestPercent()+Dew[i].GetAlpinePercent()+
        Dew[i].GetHeatherPercent()+Dew[i].GetBedrockPercent()+
        Dew[i].GetLakePercent()+Dew[i].GetGlacierPercent()+Dew[i].GetBogPercent();

      if (value > 100.0) {
        Dew[i].SetForestPercent(Dew[i].GetForestPercent()*100.0/value);
        Dew[i].SetAlpinePercent(Dew[i].GetAlpinePercent()*100.0/value);
        Dew[i].SetHeatherPercent(Dew[i].GetHeatherPercent()*100.0/value);
        Dew[i].SetBedrockPercent(Dew[i].GetBedrockPercent()*100.0/value);
        Dew[i].SetLakePercent(Dew[i].GetLakePercent()*100.0/value);
        Dew[i].SetGlacierPercent(Dew[i].GetGlacierPercent()*100.0/value);
        Dew[i].SetBogPercent(Dew[i].GetBogPercent()*100.0/value);

        value=Dew[i].GetForestPercent()+Dew[i].GetAlpinePercent()+
        Dew[i].GetHeatherPercent()+Dew[i].GetBedrockPercent()+
        Dew[i].GetLakePercent()+Dew[i].GetGlacierPercent()+Dew[i].GetBogPercent();
      }

      if (value > 100.0) {
        GetRowCol(i, nCols, j, k);
        cout << "\n Sum of land use percentages = " << value << "    Element : " << j << "," << k << "\n";
        Dew[i].SetOpenLandPercent(0.0);
      } else {
        Dew[i].SetOpenLandPercent(100.0-value);
      }
    }
  }
}


// Land surface classes general
void ReadLandUseGeneral(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[1024];
  int i,j,k;
  int nRo, nCo, noDa;
  double value, notAssigned, previousValue;
  double xllC, yllC, cellS;

  // Area of grid cells 
  for (i=0; i<nRows*nCols; i++) {
    Dew[i].SetArea(cellSize*cellSize);
  }

  // Read elevation of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with grid cell elevations: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finElevation(fileName);  // Open for reading
  if (!finElevation.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finElevation >> buffer >> nCo;
  finElevation >> buffer >> nRo;
  finElevation >> buffer >> xllC;
  finElevation >> buffer >> yllC;
  finElevation >> buffer >> cellS;
  finElevation >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for elevation!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finElevation >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetElevation(value);
    //    previousValue=value;
  }
  finElevation.close();
  // Read end

  // Read slope length of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope lengths: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finSlopeLength(fileName);  // Open for reading
  if (!finSlopeLength.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlopeLength >> buffer >> nCo;
  finSlopeLength >> buffer >> nRo;
  finSlopeLength >> buffer >> xllC;
  finSlopeLength >> buffer >> yllC;
  finSlopeLength >> buffer >> cellS;
  finSlopeLength >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlopeLength >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetSlopeLength(value);
    //    previousValue=value;
  }
  finSlopeLength.close();
  // Read end

  // Read slope angle of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope angles: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finSlopeAngle(fileName);  // Open for reading
  if (!finSlopeAngle.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlopeAngle >> buffer >> nCo;
  finSlopeAngle >> buffer >> nRo;
  finSlopeAngle >> buffer >> xllC;
  finSlopeAngle >> buffer >> yllC;
  finSlopeAngle >> buffer >> cellS;
  finSlopeAngle >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlopeAngle >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetSlopeAngle(value);
    //    previousValue=value;
  }
  finSlopeAngle.close();
  // Read end

  // Read slope aspect of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope aspects: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finAspect(fileName);  // Open for reading
  if (!finAspect.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finAspect >> buffer >> nCo;
  finAspect >> buffer >> nRo;
  finAspect >> buffer >> xllC;
  finAspect >> buffer >> yllC;
  finAspect >> buffer >> cellS;
  finAspect >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for aspect!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finAspect >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetAspect(value);
    //    previousValue=value;
  }
  finAspect.close();
  // Read end

  // Read percentage of grid cells covered by lakes
  //  previousValue=-9999.0;
  /*cout << " File with lake percentage: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finLakePercent(fileName);  // Open for reading
  if (!finLakePercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finLakePercent >> buffer >> nCo;
  finLakePercent >> buffer >> nRo;
  finLakePercent >> buffer >> xllC;
  finLakePercent >> buffer >> yllC;
  finLakePercent >> buffer >> cellS;
  finLakePercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for lake percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finLakePercent >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetLakePercent(value);
    //    previousValue=value;
  }
  finLakePercent.close();
  // Read end

  // Read percentage of grid cells covered by all land surface classes, number = numberLandSurfaceClasses-1
  for (j=0; j<numberLandSurfaceClasses-1; j++) {
    //    previousValue=-9999.0;
    cout << " File with land surface class percentage: " << j << endl;
    /*    cin >> fileName;
      cout << endl; */
    fileControl.ignore(100,':');
    fileControl >> fileName;
    fileControl.ignore(1024,'\n');
    ifstream finLandSurfacePercent(fileName);  // Open for reading
    if (!finLandSurfacePercent.is_open()) {
      cout << endl << " Error opening file " << fileName << endl << endl;
      exit (1);
    }
    cout << fileName << endl;
    finLandSurfacePercent >> buffer >> nCo;
    finLandSurfacePercent >> buffer >> nRo;
    finLandSurfacePercent >> buffer >> xllC;
    finLandSurfacePercent >> buffer >> yllC;
    finLandSurfacePercent >> buffer >> cellS;
    finLandSurfacePercent >> buffer >> noDa;
    if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
      cout << "\n\n Error in grid header information for forest percent!\n\n";
      cout << " nCols        " << nCo << "\t  " << nCols << endl;
      cout << " nRows        " << nRo << "\t  " << nRows << endl;
      cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
      cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
      cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
      cout << " noData       " << noDa << "\t  " << noData << endl << endl;
      exit (1);
    }
    for (i=0; i<nRows*nCols; i++) {
      finLandSurfacePercent >> value;
      /*      if (j==0) {
	if (i==48921 || i==49341 || i==49760 || i==49761 || i==50180 || i==50181) cout << "land " << i << "  " << value << "\n";
	}*/
      if (value<0.0) value=0.0;
      //    if (value<0.0 && previousValue>=0.0) value=previousValue;
      //      else if (value<0.0) value=0.0;
      Dew[i].SetLandSurfacePercent(j,value);
      //      previousValue=value;
    }
    finLandSurfacePercent.close();
  }
  // Read end

  // Read percentage of grid cells covered by glaciers
  //  previousValue=-9999.0;
  /*cout << " File with glacier percentage: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGlacierPercent(fileName);  // Open for reading
  if (!finGlacierPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierPercent >> buffer >> nCo;
  finGlacierPercent >> buffer >> nRo;
  finGlacierPercent >> buffer >> xllC;
  finGlacierPercent >> buffer >> yllC;
  finGlacierPercent >> buffer >> cellS;
  finGlacierPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierPercent >> value;
    //    if (i==48921 || i==49341 || i==49760 || i==49761 || i==50180 || i==50181) cout << "glac " << i << "  " << value << "\n";
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetGlacierPercent(value);
    //    previousValue=value;
  }
  finGlacierPercent.close();
  // Read end

  // Read elevation of glacier surface
  //  previousValue=-9999.0;
  /*cout << " File with glacier surface elevations: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGlacierSurfaceElevation(fileName);  // Open for reading
  if (!finGlacierSurfaceElevation.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierSurfaceElevation >> buffer >> nCo;
  finGlacierSurfaceElevation >> buffer >> nRo;
  finGlacierSurfaceElevation >> buffer >> xllC;
  finGlacierSurfaceElevation >> buffer >> yllC;
  finGlacierSurfaceElevation >> buffer >> cellS;
  finGlacierSurfaceElevation >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierSurfaceElevation >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetGlacierSurfaceElevation(value);
    //    previousValue=value;
  }
  finGlacierSurfaceElevation.close();
  // Read end

  // Read thickness of glacier ice
  //  previousValue=-9999.0;
  /*cout << " File with glacier ice thickness: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGlacierIceThickness(fileName);  // Open for reading
  if (!finGlacierIceThickness.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierIceThickness >> buffer >> nCo;
  finGlacierIceThickness >> buffer >> nRo;
  finGlacierIceThickness >> buffer >> xllC;
  finGlacierIceThickness >> buffer >> yllC;
  finGlacierIceThickness >> buffer >> cellS;
  finGlacierIceThickness >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierIceThickness >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    Dew[i].SetGlacierIceThickness(value);
    //    previousValue=value;
  }
  finGlacierIceThickness.close();
  // Read end


  // Control and correct area fractions of land use classes
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetLandIndex()!=noData) { 
      //      if (Dew[i].GetLakePercent() > 0) cout << "\n  cell no. " << i << " glac " << Dew[i].GetLakePercent() << "\n";
      //      if (Dew[i].GetGlacierPercent() > 0) cout << "\n cell no. " << i << " glac " << Dew[i].GetGlacierPercent() << "\n";
      value=Dew[i].GetLakePercent()+Dew[i].GetGlacierPercent();
      for (j=0; j<numberLandSurfaceClasses-1; j++) {
	//	if (value > 0 && Dew[i].GetLandSurfacePercent(j) > 0) cout << " cell no. " << i << " valu " << value << " land " << Dew[i].GetLandSurfacePercent(j) << "\n";
	value=value+Dew[i].GetLandSurfacePercent(j);
      }
      //      if (value > 100) cout << " cell no. " << i << "      " << value << "\n";
      if (value < 1.0) {
        Dew[i].SetLandSurfacePercent(0,100.0);
        value = 100.0;
      }
      /*      else if (value < 100.0) {
              Dew[i].SetLandSurfacePercent(2,100.0-value);
              value=value+Dew[i].GetLandSurfacePercent(2);
              }
	      if (value != 100.0) {*/
      else if (value != 100.0) {
        Dew[i].SetLakePercent(Dew[i].GetLakePercent()*100.0/value);
        Dew[i].SetGlacierPercent(Dew[i].GetGlacierPercent()*100.0/value);
        for (j=0; j<numberLandSurfaceClasses-1; j++) Dew[i].SetLandSurfacePercent(j,Dew[i].GetLandSurfacePercent(j)*100.0/value);

        value=Dew[i].GetLakePercent()+Dew[i].GetGlacierPercent();
        for (j=0; j<numberLandSurfaceClasses-1; j++) value=value+Dew[i].GetLandSurfacePercent(j);
      }

      if (value != 100.0) {
        GetRowCol(i, nCols, j, k);
        cout << "\n Sum of land use percentages = " << value << "    Element : " << j << "," << k << "\n";
      }
    }
  }
}


void WriteLandScapeElements(DistributedElement * const Dew, ParametersGeneral * const ParGeneralStore,
                            MeteorologicalStations * const MetStations, int landIndex, 
                            int preSetModelStructure, int nRows, int nCols, int noData, 
                            double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i, k, n, maxIndex;
  double lakePer, glacPer, glacSurfElev, glacIceThick, areaCorr[2], totalArea, tempArea, areaFraction[numberLandSurfaceClasses-1];
  MODEL_STRUCTURE modelStructure[numberModelStructures];
  LANDSURFACE tempLandSurf, landSurfType[numberLandSurfaceClasses-1];
  SOIL soilType[numberSoilClasses-1];
  int * precipitationStation = new int [ParGeneralStore->GetNUM_PREC_SERIES()];
  int * temperatureStation = new int [ParGeneralStore->GetNUM_TEMP_SERIES()];
  double * precipitationWeight = new double [ParGeneralStore->GetNUM_PREC_SERIES()];
  double * temperatureWeight = new double [ParGeneralStore->GetNUM_TEMP_SERIES()];

  for (k=0; k<numberModelStructures; k++) modelStructure[k]=MODEL_STRUCTURE(k);
  // Landscape elements  
  ofstream landScapeOut("dew_landscape.txt");  // Open for writing
  ofstream geoIndexOut("dew_grid_index.txt");  // Open for writing
  landScapeOut.precision(0); landScapeOut.setf(ios::fixed); 
  landScapeOut << "ncols         " << nCols << endl;
  landScapeOut << "nrows         " << nRows << endl;
  landScapeOut << "xllcorner     " << xllCorner << endl;
  landScapeOut << "yllcorner     " << yllCorner << endl;
  landScapeOut << "cellsize      " << cellSize << endl;
  landScapeOut << "NODATA_value  " << noData << endl;
  landScapeOut << "# Number of landscape elements :  " << landIndex << endl;
  geoIndexOut << " Index of grid cells relative to upper left corner of rectangle with no. of rows = " << nRows << " and no. of columns = " << nCols << endl;
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetLandIndex()!=noData) {
      // Sort land surface types based on area
      for (k=0; k<numberLandSurfaceClasses-1; k++) landSurfType[k]=LANDSURFACE(k);
      for (k=0; k<numberLandSurfaceClasses-1; k++) areaFraction[k] = 0.0;
      // Area fraction for land surface classes based on potential tree level
      /*      areaFraction[0]=Dew[i].GetOpenLandPercent();      // Open land (meadows, agriculture)
      areaFraction[1]=Dew[i].GetBogPercent();           // Bogs
      areaFraction[2]=Dew[i].GetForestPercent();        // Forest
      areaFraction[3]=Dew[i].GetAlpinePercent();        // Alpine forest
      areaFraction[4]=Dew[i].GetHeatherPercent();       // Low mountain with vegetation
      areaFraction[5]=Dew[i].GetBedrockPercent();       // Exposed bedrock (high mountain)*/
      // Area fraction for land surface classes general
      for (k=0; k<numberLandSurfaceClasses-1; k++) areaFraction[k] = Dew[i].GetLandSurfacePercent(k);
      // End area fraction for land surface classes
      for (k=0; k<numberLandSurfaceClasses-1; k++) {
        maxIndex=k;
        for (n=k+1; n<numberLandSurfaceClasses-1; n++) {
          if (areaFraction[n] > areaFraction[maxIndex]) {
            tempLandSurf=landSurfType[maxIndex];
            landSurfType[maxIndex]=landSurfType[n];
            landSurfType[n]=tempLandSurf;
            tempArea=areaFraction[maxIndex];
            areaFraction[maxIndex]=areaFraction[n];
            areaFraction[n]=tempArea; 
          }
        }
      }
      for (k=0; k<numberLandSurfaceClasses-1; k++) soilType[k]=SOIL(landSurfType[k]);
      lakePer=Dew[i].GetLakePercent();
      glacPer=Dew[i].GetGlacierPercent();
      glacSurfElev=Dew[i].GetGlacierSurfaceElevation();
      glacIceThick=Dew[i].GetGlacierIceThickness();
      areaCorr[0]=areaFraction[0];
      areaCorr[1]=areaFraction[1];
      //      cout << lakePer << "    " << glacPer << "    "  << areaCorr[0] << "    "
      //           << areaCorr[1] << endl;
      totalArea = lakePer + glacPer + areaCorr[0] + areaCorr[1];
      if (totalArea!=100.0) {
        //      if (areaFraction[0]==0.0 && areaFraction[1]==0.0) {
        //        lakePer = lakePer*100/totalArea;
        //        glacPer = glacPer*100/totalArea;
        //        } else {
        areaCorr[0]=areaFraction[0]+areaFraction[0]*(100.0-totalArea)/(areaFraction[0]+areaFraction[1]);
        areaCorr[1]=areaFraction[1]+areaFraction[1]*(100.0-totalArea)/(areaFraction[0]+areaFraction[1]);
        //       }
        totalArea = lakePer + glacPer + areaCorr[0] + areaCorr[1];
        if (totalArea!=100.0) cout << endl << "totalArea    " << totalArea << endl;
      }

      /* Find meteorological stations */
      DistanceSort(&Dew[i], ParGeneralStore, MetStations, 
                   precipitationStation, temperatureStation, precipitationWeight, temperatureWeight);

      landScapeOut.precision(4); landScapeOut.setf(ios::fixed); landScapeOut.setf(ios::showpoint);
      landScapeOut.width(10); landScapeOut << Dew[i].GetLandIndex() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetGeoIndex() << "  ";
      landScapeOut.width(10); landScapeOut << modelStructure[preSetModelStructure] << "  ";
      landScapeOut.width(15); landScapeOut << Dew[i].GetArea() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetElevation() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetSlopeLength() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetSlopeAngle() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetAspect() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetFlowLandValue() << "  ";
      landScapeOut.width(10); landScapeOut << lakePer << "  ";
      landScapeOut.width(10); landScapeOut << glacPer << "  ";
      landScapeOut.width(10); landScapeOut << glacSurfElev << "  ";
      landScapeOut.width(10); landScapeOut << glacIceThick << "  ";
      for (k=0; k<maximumNumberLandClasses; k++) {
        landScapeOut.width(10); landScapeOut << landSurfType[k] << "  ";
        landScapeOut.width(10); landScapeOut << soilType[k] << "  ";
        landScapeOut.width(10); landScapeOut << areaCorr[k] << "  ";
      }
      for (k=0; k<ParGeneralStore->GetNUM_PREC_SERIES(); k++) {
        landScapeOut.width(10); landScapeOut << precipitationStation[k] << "  ";
        landScapeOut.width(10); landScapeOut << precipitationWeight[k] << "  ";
      }
      for (k=0; k<ParGeneralStore->GetNUM_TEMP_SERIES(); k++) {
        landScapeOut.width(10); landScapeOut << temperatureStation[k] << "  ";
        landScapeOut.width(10); landScapeOut << temperatureWeight[k] << "  ";
      }
      landScapeOut<< endl;
      geoIndexOut.width(10); geoIndexOut << "  " << Dew[i].GetGeoIndex() << "  " << endl;
    }
  }
  landScapeOut << endl;
  landScapeOut.close();
  geoIndexOut.close();

  delete [] precipitationStation;
  delete [] temperatureStation;
  delete [] precipitationWeight;
  delete [] temperatureWeight;
}


void DistanceSort(DistributedElement * const distElement, ParametersGeneral * const ParGeneralStore, 
                  MeteorologicalStations * const MetStations,
                  int * precipitationStation, int * temperatureStation, double * precipitationWeight, double * temperatureWeight)
{
  int i, j;
  int numP, minIndex;
  int temporaryStation;
  double temporaryDistance;
  double sumWeight;
  int * precipitationSort = new int [MetStations->GetNumPrecStations()];
  int * temperatureSort = new int [MetStations->GetNumTempStations()];
  double * distancePrec = new double [MetStations->GetNumPrecStations()];
  double * distanceTemp = new double [MetStations->GetNumTempStations()];

  for (i=0; i<MetStations->GetNumPrecStations(); i++) {
    precipitationSort[i] = i;
    distancePrec[i] = sqrt(pow((distElement->GetXCoord() - MetStations->GetStationCoordX(i)),2.0) 
                           + pow((distElement->GetYCoord() - MetStations->GetStationCoordY(i)),2.0));
    //    cout << "Prec " << precipitationSort[i] << "  " << distancePrec[i] << endl;
  }
  for (i=0; i<MetStations->GetNumPrecStations(); i++) {
    minIndex = i;
    for (j=i+1; j<MetStations->GetNumPrecStations(); j++) {
      if (distancePrec[j] < distancePrec[minIndex]) {
        temporaryDistance = distancePrec[minIndex];
        distancePrec[minIndex] = distancePrec[j];
        distancePrec[j] = temporaryDistance;
        temporaryStation = precipitationSort[minIndex];
        precipitationSort[minIndex] = precipitationSort[j];
        precipitationSort[j] = temporaryStation;
      }
    }
  }
  //  cout << "PREC  " << precipitationSort[0] << "  " << distancePrec[0] << endl;

  numP = MetStations->GetNumPrecStations();
  for (i=0; i<MetStations->GetNumTempStations(); i++) {
    temperatureSort[i] = i;
    distanceTemp[i] = sqrt(pow((distElement->GetXCoord() - MetStations->GetStationCoordX(numP+i)),2.0) 
                    + pow((distElement->GetYCoord() - MetStations->GetStationCoordY(numP+i)),2.0));
    //    cout << "Temp " << temperatureSort[i] << "  " << distanceTemp[i] << endl;
  }
  for (i=0; i<MetStations->GetNumTempStations(); i++) {
    minIndex = i;
    for (j=i+1; j<MetStations->GetNumTempStations(); j++) {
      if (distanceTemp[j] < distanceTemp[minIndex]) {
        temporaryDistance = distanceTemp[minIndex];
        distanceTemp[minIndex] = distanceTemp[j];
        distanceTemp[j] = temporaryDistance;
        temporaryStation = temperatureSort[minIndex];
        temperatureSort[minIndex] = temperatureSort[j];
        temperatureSort[j] = temporaryStation;
      }
    }
  }
  //  cout << "TEMP  "  << temperatureSort[0] << "  " << distanceTemp[0] << endl;

  sumWeight = 0.0;
  for (i=0; i<ParGeneralStore->GetNUM_PREC_SERIES(); i++) {
    precipitationStation[i] = precipitationSort[i];
    precipitationWeight[i] = 1.0/distancePrec[i];
    sumWeight = sumWeight + precipitationWeight[i];
  }
  if (sumWeight != 1.0) {
    for (i=0; i<ParGeneralStore->GetNUM_PREC_SERIES(); i++) {
      precipitationWeight[i] = precipitationWeight[i]/sumWeight;
    }
  }

  sumWeight = 0.0;
  for (i=0; i<ParGeneralStore->GetNUM_TEMP_SERIES(); i++) {
    temperatureStation[i] = temperatureSort[i];
    temperatureWeight[i] = 1.0/distanceTemp[i];
    sumWeight = sumWeight + temperatureWeight[i];
  }
  if (sumWeight != 1.0) {
    for (i=0; i<ParGeneralStore->GetNUM_TEMP_SERIES(); i++) {
      temperatureWeight[i] = temperatureWeight[i]/sumWeight;
    }
  }

  //  cout << distElement->GetLandIndex() << "  " << distElement->GetGeoIndex() << "  "
  //       << distElement->GetXCoord() << "  " << distElement->GetYCoord() << endl;

  delete [] precipitationSort;
  delete [] temperatureSort;
  delete [] distancePrec;
  delete [] distanceTemp;

}


void GetRowCol(int elementNo, int nCols, int &row, int &col) 
{
  row = elementNo/nCols; 
  col = elementNo%nCols; 
}


void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  int numPrec, numTemp;

  /*  cout << " File with common parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGeneralPar(fileName);  // Open for reading
  if (!finGeneralPar.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finGeneralPar.ignore(1024,'\n');
  finGeneralPar.ignore(100,':'); finGeneralPar >> numPrec; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> numTemp; 
  ParGeneralStore->SetNUM_PREC_SERIES(numPrec);
  ParGeneralStore->SetNUM_TEMP_SERIES(numTemp);
  finGeneralPar.close();
  fout << "Common parameters: \n";
  fout << ParGeneralStore->GetNUM_PREC_SERIES() << endl;
  fout << ParGeneralStore->GetNUM_TEMP_SERIES() << endl;
  fout << endl << endl;
}


void ReadFlowDirection(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[1024];
  int i;
  int nRo, nCo, noDa;
  int value;
  double xllC, yllC, cellS;

  // Read flow direction grid for landscape elements
  /*  cout << endl << "\n File with flow direction grid for landscape elements: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finFlowLand(fileName);  // Open for reading
  if (!finFlowLand.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finFlowLand >> buffer >> nCo;
  finFlowLand >> buffer >> nRo;
  finFlowLand >> buffer >> xllC;
  finFlowLand >> buffer >> yllC;
  finFlowLand >> buffer >> cellS;
  finFlowLand >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for flow directions for landscape elements!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finFlowLand >> value;
    if (!ControlFlowDirection(value, noData)) {
      cout << "\n\n Error in flow directions for landscape elements!\n\n";
      cout << " value     " << value << endl;
      cout << " noData    " << noData << endl << endl;
      exit (1);
    }
    Dew[i].SetFlowLandValue(value);
  }
  finFlowLand.close();
  // Read end
}


int ControlFlowDirection(int flowDirection, int noData)
{
  return (flowDirection==1 || flowDirection==2 || flowDirection==4 || flowDirection==8 ||
          flowDirection==16 || flowDirection==32 || flowDirection==64 || flowDirection==128 ||
          flowDirection==noData);
}


void FindSubCatchmentIdentifier(DistributedElement * const Dew, SubCatchment * const CatchmentElement, int numWatc, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[1024];
  int i,j;
  int nRo, nCo, noDa;
  int subCatchmentId;
  double xllC, yllC, cellS;
  bool subCatchmentFound;
  bool noElement=false;
  DistributedElement * lastElement;

  // Read watercourse/sub-catchment identifiers for landscape elements
  /*  cout << "\n File with watercourse/sub-catchment identifiers: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finId(fileName);  // Open for reading
  if (!finId.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finId >> buffer >> nCo;
  finId >> buffer >> nRo;
  finId >> buffer >> xllC;
  finId >> buffer >> yllC;
  finId >> buffer >> cellS;
  finId >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for watercourse/sub-catchment elements!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  // Read end

  // Connect landscape elements to sub-catchment outlets
  for (i=0; i<nRows*nCols; i++) {
    finId >> subCatchmentId;
    Dew[i].SetSubCatchmentValue(subCatchmentId);
    j=0;
    subCatchmentFound = false;
    while (j<numWatc) {
      if (subCatchmentId==CatchmentElement[j].GetIdentifier()) {
        if (!subCatchmentFound) {
          subCatchmentFound=true;
          if (CatchmentElement[j].GetLandScapeElement()) {
            lastElement = CatchmentElement[j].GetLandScapeElement();
            while (lastElement->GetNextElement()) lastElement = lastElement->GetNextElement();
            lastElement->SetNextElement(&Dew[i]);
          } else {
            CatchmentElement[j].SetLandScapeElement(&Dew[i]);
          }
          CatchmentElement[j].SetNumLandScape(CatchmentElement[j].GetNumLandScape()+1);
        } else {
          cout << "\n\n Error in watercourse/sub_catchment identifiers!   ";
          cout << CatchmentElement[j].GetIdentifier() << endl << endl;
          exit (1);
        }
      } 
      j++;
    }
  }

  for (j=0; j<numWatc; j++) {
    if (!(CatchmentElement[j].GetLandScapeElement())) {
      cout << "\n No landscape elements for watercourse/sub-catchment index " << j << "  "  << "Watercourse/sub-catchment identifier " << CatchmentElement[j].GetIdentifier() << endl;
      noElement=true;
    }
  }
  if (noElement) {
    cout << "\n Program is terminated! " << endl;
    exit(1);
  }
  finId.close();
    
}


void WriteSubCatchmentIdentifier(DistributedElement * const Dew, SubCatchment * const CatchmentElement, int numWatc, ofstream &fout)
{
  int i;
  DistributedElement * thisElement;

  // Watercourse/sub-catchment elements
  ofstream subCatchmentOut("dew_waterland.txt");  // Open for writing
  for (i=0; i<numWatc; i++) {
    if (CatchmentElement[i].GetLandScapeElement()) {
      subCatchmentOut << "#  " << CatchmentElement[i].GetIdentifier() << "  #  " << CatchmentElement[i].GetNumLandScape() << endl;
      thisElement = CatchmentElement[i].GetLandScapeElement();
      while (thisElement) {
        subCatchmentOut << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
        thisElement = thisElement->GetNextElement();
      }
    }
  }
  subCatchmentOut << endl;
  subCatchmentOut.close();
}


void FindUpLandFlow(DistributedElement * const Dew, int nRows, int nCols, int noData)
{
  int i,j;
  int nUp;

  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (Dew[ELEMENT(i,j)].GetLandIndex()!=noData) {
        nUp=0;
        if ((j<nCols-1) 
            && (Dew[ELEMENT(i,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i,j+1)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i,j+1)].GetFlowLandValue()==16)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i,j+1)]);
          nUp++;
        }
        if ((i<nRows-1 && j<nCols-1) 
            && (Dew[ELEMENT(i+1,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j+1)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i+1,j+1)].GetFlowLandValue()==32)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j+1)]);
          nUp++;
        }
        if ((i<nRows-1)
            && (Dew[ELEMENT(i+1,j)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i+1,j)].GetFlowLandValue()==64)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j)]);
          nUp++;
        }
        if ((i<nRows-1) && (j>0) 
            && (Dew[ELEMENT(i+1,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j-1)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i+1,j-1)].GetFlowLandValue()==128)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j-1)]);
          nUp++;
        }
        if ((j>0) 
            && (Dew[ELEMENT(i,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i,j-1)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i,j-1)].GetFlowLandValue()==1)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i,j-1)]);
          nUp++;
        }
        if ((i>0 && j>0) 
            && (Dew[ELEMENT(i-1,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j-1)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i-1,j-1)].GetFlowLandValue()==2)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j-1)]);
          nUp++;
        }
        if ((i>0)
            && (Dew[ELEMENT(i-1,j)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i-1,j)].GetFlowLandValue()==4)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j)]);
          nUp++;
        }
        if ((i>0 && j<nCols-1)
            && (Dew[ELEMENT(i-1,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j+1)].GetSubCatchmentValue()==noData)
            && (Dew[ELEMENT(i-1,j+1)].GetFlowLandValue()==8)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j+1)]);
          nUp++;
        }
        if (nUp > 7) {
          cout << "\n\n Error in function FindUpLandFlow!\n\n";
          cout << "\n nUp: " << nUp << endl << endl;
          exit (1);
        }      
        Dew[ELEMENT(i,j)].SetNumUpLand(nUp);
      }
    }
  }
  return;
}


void WriteLandScapeAndSubCatchmentFlow(DistributedElement * const Dew, int nRows, int nCols, int noData, ofstream &fout)
{
  int i,j;

  // Pointers between landscape elements  
  ofstream landUpFlow("dew_landupflow.txt");  // Open for writing
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetNumUpLand()>0) {
      landUpFlow << Dew[i].GetLandIndex() << "  " << Dew[i].GetNumUpLand() << " :  ";
      for (j=0; j<Dew[i].GetNumUpLand(); j++) 
        landUpFlow << Dew[i].GetUpLandFlow(j)->GetLandIndex() << "  ";
      landUpFlow << endl;
    }
  }
  landUpFlow << endl;
  landUpFlow.close();
}
