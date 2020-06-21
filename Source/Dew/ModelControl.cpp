#include "ModelControl.h"
#include "stdafx.h"
#include "DateTime.h"

ModelControl::ModelControl() :
    modelRun('M'),
    hierarchy('H'),
    inputFormat('F'),
    evaporationModelling('T'),
    glacierModelling('G'),
    routingType('R'),
    numberValues(numberPotentialEvaporationValuesPerYear)
{
    int i;
    evaporationArray = new double [numberValues];
    for (i=0; i<numberValues; i++) 
    {
        evaporationArray[i] = missingData;
    }
}

ModelControl::~ModelControl()
{
}

void  ModelControl::SetModelRun(char value)
{
    modelRun = value;
}
char  ModelControl::GetModelRun() const
{
    return modelRun;
}
void  ModelControl::SetHierarchy(char value)
{
    hierarchy = value;
}
char  ModelControl::GetHierarchy() const
{
    return hierarchy;
}
void  ModelControl::SetInputFormat(char value)
{
    inputFormat = value;
}
char  ModelControl::GetInputFormat() const
{
    return inputFormat;
}
void  ModelControl::SetGlacierModelling(char value)
{
    glacierModelling = value;
}
char  ModelControl::GetGlacierModelling() const
{
    return glacierModelling;
}
void  ModelControl::SetRoutingType(char value)
{
    routingType = value;
}
char  ModelControl::GetRoutingType() const
{
    return routingType;
}


void ModelControl::SetModelControl(char mrun, char hier, char form, char evap, char glac, char rout)
{
    SetModelRun(mrun);
    SetHierarchy(hier);
    SetInputFormat(form);
    SetEvaporationModelling(evap);
    SetGlacierModelling(glac);
    SetRoutingType(rout);
}


void ModelControl::SetEvaporationModelling(char value)
{
  char fileNameInput[100];
  char buffer[256];
  int i;
  evaporationModelling = value;
  if (evaporationModelling == 'M' || evaporationModelling == 'm') {
    strcpy(fileNameInput,"monthlyEvaporation.txt");
    cout << " " << fileNameInput << "\n";
    ifstream finInput(fileNameInput);  // Open for reading
    if (!finInput.is_open()) {
      cout << endl << " Error opening file " << fileNameInput << endl << endl;
      exit(1);
    }
    finInput.ignore(256,'\n');
    for (i=0; i<numberValues; i++) {
      finInput.ignore(100,':'); finInput >> evaporationArray[i];
    }
    finInput.close();
    for (i=0; i<numberValues; i++) {
      cout << " " << i << "    " << GetEvaporationArray(i) << "\n";
    }
  }
}


char  ModelControl::GetEvaporationModelling() const
{
    return evaporationModelling;
}
