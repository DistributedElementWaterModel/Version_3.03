#include "InputTimeSeries.h"
#include "DateTime.h"
#include "Dew.h"

InputTimeSeries::InputTimeSeries(int numberRows, int numberColumns, DateTime firstTime, DateTime lastTime, int secondsPerTimeStep) :
    timeSteps(numberRows),
    numberSeries(numberColumns)
{
    int i;
    SetGeneralPar(0);
    datetime = new DateTime[timeSteps];
    inputArray = new double[timeSteps * numberSeries];
    datetime[0] = firstTime;
    for (i = 1; i < timeSteps; i++) 
    {
        datetime[i] = datetime[i-1] + secondsPerTimeStep; 
    }
    /*    for (i = 0; i < timeSteps; i++)
	  {
	  datetime[i] = firstTime + i * secondsPerTimeStep;
	  }*/
    for (i = 0; i < timeSteps * numberSeries; i++)
    {
        inputArray[i] = missingData;
        //    cout << i << " " << inputArray[i] << "  ";
    }
    if (datetime[0] != firstTime || datetime[timeSteps - 1] != lastTime)
    {
        cout << endl << " DateTime error during initialisation of InputTimeSeries array: " <<
             datetime[0] << "  " << firstTime << "  " << datetime[timeSteps - 1] << "  " << lastTime << endl << endl;
        exit(1);
    }
    //  cout << datetime[0] << "  " << firstTime  << "  " << datetime[timeSteps-1]  << "  " << lastTime << endl << endl;
}

InputTimeSeries::~InputTimeSeries()
{
}

int  InputTimeSeries::GetNumberTimeSteps() const
{
    return timeSteps;
}
int  InputTimeSeries::GetNumberInputSeries() const
{
    return numberSeries;
}
void  InputTimeSeries::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral * InputTimeSeries::GetGeneralPar() const
{
    return commonPar;
}
DateTime  InputTimeSeries::GetDateTime(int i) const
{
    return datetime[i];
}
double  InputTimeSeries::GetInput(int i, int j) const
{
    return inputArray[i * numberSeries + j];
}

void InputTimeSeries::SetInput(ifstream &fin)
{
    char ch;
    char buffer[10240];
    int i, j;
    int date, time, year, mth, day, hour, min;
    double * tempArray = new double[numberSeries];
    fin.getline(buffer, 10240);
    fin.getline(buffer, 10240);
    i = 0;
    while (fin >> date >> ch >> time)
    {
        for (j = 0; j < numberSeries; j++)
        {
            fin >> tempArray[j];
        }
        year = date / 10000;
        mth = (date % 10000) / 100;
        day = date % 100;
        hour = time / 100;
        min = time % 100;
        //    cout << year << " " << mth << " " << day << " " << hour << " " << min << " " << endl;
        while ((datetime[i].getYear() < year ||
                (datetime[i].getYear() == year && datetime[i].getMonth() < mth) ||
                (datetime[i].getYear() == year && datetime[i].getMonth() == mth && datetime[i].getDay() < day) ||
                (datetime[i].getYear() == year && datetime[i].getMonth() == mth && datetime[i].getDay() == day && datetime[i].getHour() < hour) ||
                (datetime[i].getYear() == year && datetime[i].getMonth() == mth && datetime[i].getDay() == day &&
                 datetime[i].getHour() == hour && datetime[i].getMinute() < min)) &&
                i < timeSteps - 1)
        {
            i++;
            //      cout << " i = " << i << endl;
        }
        //    cout << datetime[i].getYear() << " " << datetime[i].getMonth() << " " << datetime[i].getDay() << " " <<
        //      datetime[i].getHour() << " " << datetime[i].getMinute() << endl;
        if (datetime[i].getYear() == year && datetime[i].getMonth() == mth && datetime[i].getDay() == day &&
                datetime[i].getHour() == hour && datetime[i].getMinute() == min)
        {
            //      cout << " data found  " << i << " " << year << " " << mth << " " << day << " " << hour << " " << min << " " << endl;
            for (j = 0; j < numberSeries; j++)
            {
                inputArray[i * numberSeries + j] = tempArray[j];
            }
        }
    }
    delete[] tempArray;
}

void InputTimeSeries::WriteInput()
{
    FILE *fp_out;
    int i, j;
    if ((fp_out = fopen("input_out.txt", "w")) == NULL)
    {
        printf("\n File input_out.txt not found!\n\n");
        exit(1);
    }
    for (i = 0; i < timeSteps; i++)
    {
        fprintf(fp_out, "%04d%02d%02d/%02d%02d", GetDateTime(i).getYear(), GetDateTime(i).getMonth(), GetDateTime(i).getDay(),
                GetDateTime(i).getHour(), GetDateTime(i).getMinute());
        for (j = 0; j < numberSeries; j++)
        {
            fprintf(fp_out, "%15.5f", GetInput(i, j));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
}

