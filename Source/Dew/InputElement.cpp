#include "InputElement.h"
#include "ParametersGeneral.h"
#include "Dew.h"
#include "Util.h"

#include "stdafx.h"


InputElement::InputElement(int value) :
    numberValues(value)
{
    int i;
    inputArray = new double[numberValues];
    for (i = 0; i < numberValues; i++)
    {
        inputArray[i] = missingData;
    }
}

InputElement::~InputElement()
{
}

void InputElement::SetInput(int i, double value)
{
    inputArray[i] = value;
}
double InputElement::GetInput(int i) const
{
    return inputArray[i];
}

/* Precipitation linear elevation gradient *//*
double InputElement::PrecipitationElevationCorrected(ParametersGeneral * const commonPar, double precipitation, double inputElevation, double outputElevation, double weight) const
{
    double correctedPrecipitation, precGradLow, precGradHigh, gradElevation;
    gradElevation = commonPar->GetGRAD_CHANGE_ALT();
    precGradLow = commonPar->GetPREC_GRAD_LOW() - 1.0;
    precGradHigh = commonPar->GetPREC_GRAD_HIGH() - 1.0;
    if (gradElevation == 0.0 || (outputElevation <= gradElevation && inputElevation <= gradElevation))
    {
        correctedPrecipitation = (precipitation) * (1.0 + precGradLow * (outputElevation - inputElevation) / 100.0) * weight;
    }
    else if (outputElevation > gradElevation && inputElevation < gradElevation)
      correctedPrecipitation = (precipitation) * (1.0 + (precGradLow * (gradElevation - inputElevation) / 100.0) +
						  (precGradHigh * (outputElevation - gradElevation) / 100.0)) * weight;
    else if (outputElevation < gradElevation && inputElevation > gradElevation)
      correctedPrecipitation = (precipitation) * (1.0 + (precGradLow * (outputElevation - gradElevation) / 100.0) +
						  (precGradHigh * (gradElevation - inputElevation) / 100.0)) * weight;
    else
    {
        correctedPrecipitation = (precipitation) * (1.0 + precGradHigh * (outputElevation - inputElevation) / 100.0) * weight;
    }
    if (correctedPrecipitation < 0.0)  correctedPrecipitation = 0.0;
    return correctedPrecipitation;
}*/

/* Precipitation exponential elevation gradient */
double InputElement::PrecipitationElevationCorrected(ParametersGeneral * const commonPar, double precipitation, double inputElevation, double outputElevation, double weight) const
{
    double correctedPrecipitation, precGradLow, precGradHigh, gradElevation;
    gradElevation = commonPar->GetGRAD_CHANGE_ALT();
    precGradLow = commonPar->GetPREC_GRAD_LOW();
    precGradHigh = commonPar->GetPREC_GRAD_HIGH();
    if (gradElevation == 0.0 || (outputElevation <= gradElevation && inputElevation <= gradElevation))
    {
        correctedPrecipitation = (precipitation) * pow(precGradLow, (outputElevation - inputElevation) / 100.0) * weight;
    }
    else if (outputElevation > gradElevation && inputElevation < gradElevation)
        correctedPrecipitation = (precipitation) * pow(precGradLow, (gradElevation - inputElevation) / 100.0) *
                                 pow(precGradHigh, (outputElevation - gradElevation) / 100.0) * weight;
    else if (outputElevation < gradElevation && inputElevation > gradElevation)
        correctedPrecipitation = (precipitation) * pow(precGradLow, (outputElevation - gradElevation) / 100.0) *
                                 pow(precGradHigh, (gradElevation - inputElevation) / 100.0) * weight;
    else
    {
        correctedPrecipitation = (precipitation) * pow(precGradHigh, (outputElevation - inputElevation) / 100.0) * weight;
    }
    if (correctedPrecipitation < 0.0)  correctedPrecipitation = 0.0;
    return correctedPrecipitation;
}

double InputElement::TemperatureElevationCorrected(ParametersGeneral * const commonPar, double temperature, double precipitation, double inputElevation, double outputElevation, double weight) const
{
    double correctedTemperature, lapseRate, lapseRateDry, lapseRateWet;
    lapseRateDry = commonPar->GetLAPSE_DRY();
    lapseRateWet = commonPar->GetLAPSE_WET();
    if (precipitation == 0.0)
    {
        lapseRate = lapseRateDry;
    }
    else
    {
        lapseRate = lapseRateWet;
    }
    correctedTemperature = (temperature + lapseRate * (outputElevation - inputElevation) / 100.0) * weight;
    return correctedTemperature;
}
