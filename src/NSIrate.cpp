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
#include "DEIntegrator.h"
#include "DEIntegrationConstants.h"

const double MN = 0.9383;          //mass of nucleon in GeV
const double ME = 0.000510998;     //mass of electron in GeV
const double GeVperKG = 5.6094e26; //conversion of GeV to kg
const double SSW = 0.2387;         //sin^2(theta_w)

const double GVP = 0.02849;//SM vector proton coupling
const double GVN = -0.512213;//SM vector neutron coupling
const double GAP = 0.61689;  //SM axial proton coupling
const double GAN = -0.598426;//SM axial neutron coupling

//returns SM+NSI rate per kg/year/keV for the jth flux
double NSIrate(double ErKeV, paramList *pL, int detj, int sourcej, int fluxj)          
{
    
    double rate = 0;
    double targetsPerKG;
    int anti = ( 0 < pL->sources[sourcej].nuFluxFlav[fluxj] ) - ( 0 > pL->sources[sourcej].nuFluxFlav[fluxj] );
    
    for(int i=0;i<pL->detectors[detj].nIso;i++)
    {
        targetsPerKG = GeVperKG/(MN*pL->detectors[detj].isoA[i]); //how many targets per kg of detector
        
        if(pL->nucScat)
        {
            switch( abs(pL->sources[sourcej].nuFluxFlav[fluxj]) )
            {
                case 1:
                { 
                    //SM + ee NSI terms
                    pL->qA = 4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (GVN + pL->epEEuV + 2*pL->epEEdV ) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + ( GVP + 2*pL->epEEuV + pL->epEEdV ) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    
                    //em NSI terms
                    pL->qA = 0;//4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (pL->epEMuV + 2*pL->epEMdV) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + (2*pL->epEMuV + pL->epEMdV) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    
                    //et NSI terms
                    pL->qA = 0;//4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (pL->epETuV + 2*pL->epETdV) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + (2*pL->epETuV + pL->epETdV) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    break;
                }
                case 2:
                { 
                    //SM + mm NSI terms
                    pL->qA = 4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (GVN + pL->epMMuV + 2*pL->epMMdV ) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + ( GVP + 2*pL->epMMuV + pL->epMMdV ) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    
                    //em NSI terms
                    pL->qA = 0;//4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (pL->epEMuV + 2*pL->epEMdV) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + (2*pL->epEMuV + pL->epEMdV) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    
                    //mt NSI terms
                    pL->qA = 0;//4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (pL->epMTuV + 2*pL->epMTdV) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + (2*pL->epMTuV + pL->epMTdV) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    break;
                }
                case 3:
                {
                    //SM + tt NSI terms
                    pL->qA = 4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (GVN + pL->epTTuV + 2*pL->epTTdV ) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + ( GVP + 2*pL->epTTuV + pL->epTTdV ) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    
                    //mt NSI terms
                    pL->qA = 0;//4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (pL->epMTuV + 2*pL->epMTdV) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + (2*pL->epMTuV + pL->epMTdV) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    
                    //et NSI terms
                    pL->qA = 0;//4.0/3.0 * (pL->detectors[detj].isoJN[i]+1) / pL->detectors[detj].isoJN[i] * ( pL->detectors[detj].isoSN[i]*GAN + pL->detectors[detj].isoSZ[i]*GAP );     
                    pL->qV = ( (pL->epETuV + 2*pL->epETdV) * (pL->detectors[detj].isoA[i] - pL->detectors[detj].isoZ[i]) + (2*pL->epETuV + pL->epETdV) * pL->detectors[detj].isoZ[i] )* ffactorSI( pL->detectors[detj].isoA[i], ErKeV);     
                    rate += targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, MN*pL->detectors[detj].isoA[i], sourcej, fluxj);
                    break;
                }
            }
         }    
        if(pL->elecScat)
        {
            int Ne=0;
            while(pL->detectors[detj].ionization[i][Ne] > ErKeV && Ne < pL->detectors[detj].isoZ[i]) 
                Ne++;
        
            if(pL->sources[sourcej].isSolar[fluxj] == 1)
            {
                //nu_e component
                pL->qA = 0.5;
                pL->qV = 0.5+2*SSW;
                rate += pL->sources[sourcej].survProb[fluxj] * (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);
                //nu_mu and nu_tau component
                pL->qA = -0.5;
                pL->qV = -0.5+2*SSW;
                rate += (1-pL->sources[sourcej].survProb[fluxj]) * (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);
            }
            else
            {
                switch( abs(pL->sources[sourcej].nuFluxFlav[fluxj]) )
                {
                    case 1:
                    {
                        pL->qA = 0.5*anti+pL->epEEeA;
                        pL->qV = 0.5+2*SSW+pL->epEEeV;
                        rate += (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);
                    }
                    case 2:
                    {
                        pL->qA = -0.5*anti+pL->epMMeA;
                        pL->qV = -0.5+2*SSW+pL->epMMeV;
                        rate += (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);
                    }
                    case 3:
                    {
                        pL->qA = -0.5*anti+pL->epTTeA;
                        pL->qV = -0.5+2*SSW+pL->epTTeV;
                        rate += (pL->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pL->detectors[detj].isoFrac[i] * nuRate( ErKeV, pL, ME, sourcej, fluxj);
                    }
                }
            }
            
        }
        
    }
    
    return rate*detEff(ErKeV,pL->detectors[detj].eff); 
}

//below are functions for total rate of all fluxes
double diffNSIrate(double ErkeV, paramList *pL, int detj)              
{   
    double rate=1e-99;
    for(int fluxi=0; fluxi< pL->sources[pL->detectors[pL->detj].sourcej].numFlux; fluxi++)
    {
        if( ErkeV < (pL->detectors[detj].ErL + (double)999*(pL->detectors[detj].ErU-pL->detectors[detj].ErL)/900) )
           rate += NSIrate( ErkeV, pL, detj, pL->detectors[pL->detj].sourcej, fluxi);
    }
    return rate;
}

class NSIrateIntegral
{
public:
    paramList *pList;
    int detj;
    double operator()(double ErkeV) const
    {
        return diffNSIrate( ErkeV, pList, detj);
    }
};

double intNSIrate(double Er_min, double Er_max, paramList *pList, int detj)                          
{   
    double rate = 0;
    
    NSIrateIntegral NSIint;
    NSIint.detj = detj;
    NSIint.pList = pList;
    
    return DEIntegrator<NSIrateIntegral>::Integrate(NSIint,Er_min,Er_max,1e-6);
}

