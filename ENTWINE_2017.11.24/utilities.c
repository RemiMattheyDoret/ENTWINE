#include <stdio.h>
#include <math.h>

void sortDouble(double* vec, int length)
{
    if(length > 1)
    {
        int counter = 1;
        while(counter)
        {
            counter = 0;
            for(int i = 1; i < length; i++)
            {
                if(vec[i] < vec[i-1])
                {
                    double temp = vec[i];
                    vec[i] = vec[i-1];
                    vec[i-1] = temp;
                    counter++;
                }
            }
        }
       
    }
    return;
}

unsigned int factorial(int x)
{
    unsigned int retVal= 1;
    for(int i = 2; i <= x; i++) retVal *= i;
    return retVal;
}

double lFactorial(int x)
{
    double retVal = log(1);
    for(int i = 2; i <= x; i++) retVal += log(i);
    return retVal;
}

unsigned int nCK(int n, int k)
{
    unsigned int retVal;
    double lVal = lFactorial(n) - lFactorial(k) - lFactorial(n-k);
    retVal = (int)(rint(exp(lVal)));
    return retVal;
}