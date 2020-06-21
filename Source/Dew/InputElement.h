#pragma once

class ParametersGeneral;

class InputElement
{
public:
    InputElement(int value);
    ~InputElement();
    int GetNumberValues() const
    {
        return numberValues;
    }
    void SetInput(int i, double value);
    double GetInput(int i) const;
    double PrecipitationElevationCorrected(ParametersGeneral * const ParGeneralStore, double precipitation, double inputElevation, double outputElevation, double weight) const;
    double TemperatureElevationCorrected(ParametersGeneral * const commonPar, double temperature, double precipitation, double inputElevation, double outputElevation, double weight) const;

private:
    int numberValues;
    double * inputArray;
};
