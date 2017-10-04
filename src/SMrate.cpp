#include <iostream>
#include <cmath>
#ifndef GSL_INTERP_H
    #include <gsl/gsl_interp.h>
#endif
#include "formFactorSI.h"
#include "nuRate.h"
#include "detectorFunctions.h"
#ifndef DETECTORSTRUCT_H
    #include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif    

const double MN = 0.9383;          //mass of nucleon in GeV
const double ME = 0.000510998;     //mass of electron in GeV
const double GeVperKG = 5.6094e26; //conversion of GeV to kg
const double SSW = 0.2387;         //sin^2(theta_w)

const double GVP = 0.02849;//SM vector proton coupling
const double GVN = -0.512213;//SM vector neutron coupling
const double GAP = 0.61689;  //SM axial proton coupling
const double GAN = -0.598426;//SM axial neutron coupling

//returns SM rate per kg/year/keV for the jth flux
double SMrate(double ErKeV, paramList *pL, int detj, int sourcej, int fluxj)                  
{
    
    double rate = 0;
    double targetsPerKG;

    for(int i=0;i<pL->detectors[detj].nIso;i++)
    {
        targetsPerKG = GeVperKG/(MN*pL->detectors[detj].isoA[i]); //how many targets per kg of detector
        
        if(pL->nucScat)
        {
            pL->qA = 4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
            pL->qV = ( GVN * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + GVP * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
            rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
        }
        if(pL->elecScat)
        {
            int Ne=0;
            while(pL->detectors[detj].ionization[i][Ne] > ErKeV && Ne < pL->detectors[detj].isoZ[i]) 
                Ne++;
        
            if(pL->sources[pL->detectors[pL->detj].sourcej].isSolar[fluxj] == 1)
            {
                pL->qA = 0.5;
                pL->qV = 2*SSW+0.5;
                rate += pL->sources[pL->detectors[pL->detj].sourcej].survProb[fluxj] * (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);

                pL->qA = -0.5;
                pL->qV = 2*SSW-0.5;
                rate += (1-pL->sources[pL->detectors[pL->detj].sourcej].survProb[fluxj]) * (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);
            }
            else
            {
                pL->qA = -0.5;
                pL->qV = 0.5+2*SSW;
                rate += (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);
            }
            
        }
        
    }
   
    return rate*detEff(ErKeV,pL->detectors[detj].eff); 
}

//below are functions for total rate of all fluxes
double diffSMrate(double ErkeV, paramList *pL, int detj)              
{   
    double rate=1e-99;
    for(int fluxi=0; fluxi< pL->sources[pL->detectors[detj].sourcej].numFlux; fluxi++)
    {
        if( ErkeV > pL->detectors[detj].ErL && ErkeV < pL->detectors[detj].ErU )//(pL->detectors[detj].ErL + (double)999*(pL->detectors[detj].ErU-pL->detectors[detj].ErL)/900) )
            rate +=  SMrate( ErkeV, pL, detj, pL->detectors[detj].sourcej, fluxi); //gsl_spline_eval(pL->detectors[detj].signalSM[i], ErkeV, pL->detectors[detj].accelSM[i]);
    }
    return rate;
}

double intSMrate(double Er_min, double Er_max, paramList *pL, int detj)                          
{   
    double rate = 1e-99;
    for(int i=0; i < pL->sources[pL->detectors[detj].sourcej].numFlux; i++)
        rate += pL->sources[pL->detectors[detj].sourcej].nuFluxNorm[i] * gsl_spline_eval_integ(pL->detectors[detj].signalSM[i], Er_min, Er_max, pL->detectors[detj].accelSM[i]);
    
    return rate;
}

