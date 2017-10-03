#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int sourceInit(paramList *pList, char *sourceName, double sourceDistance);

double nuFlux(double Enu, paramList *pList, int sourcej, int fluxj);

double fluxIntegral(double ErGeV,  paramList *pList, double Mt, int EnuPow, int sourcej, int fluxj);

double fluxIntegralOsc(double ErGeV,  paramList *pList, double Mt, int EnuPow, int sourcej, int fluxj);
