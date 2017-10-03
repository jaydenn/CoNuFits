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
double NSIrate(double ErKeV, paramList *pList, int detj, int sourcej, int fluxj)		  
{
    
    double rate = 0;
    paramList pListSM = *pList;
    double targetsPerKG;
    int anti = ( 0 < pList->sources[sourcej].nuFluxFlav[fluxj] ) - ( 0 > pList->sources[sourcej].nuFluxFlav[fluxj] );
    
    for(int i=0;i<pList->detectors[detj].nIso;i++)
    {
        targetsPerKG = GeVperKG/(MN*pList->detectors[detj].isoA[i]); //how many targets per kg of detector
		
    	if(pList->nucScat)
	{

            switch( abs(pList->sources[sourcej].nuFluxFlav[fluxj]) )
            {
	            case 1:
	            {
	                //SM + ee NSI terms
                        pListSM.qA = 4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (GVN + pListSM.epEEuV + 2*pListSM.epEEdV ) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + ( GVP + 2*pListSM.epEEuV + pListSM.epEEdV ) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	                
	                //em NSI terms
	                pListSM.qA = 0;//4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (pListSM.epEMuV + 2*pListSM.epEMdV) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + (2*pListSM.epEMuV + pListSM.epEMdV) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	                
	                //et NSI terms
	                pListSM.qA = 0;//4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (pListSM.epETuV + 2*pListSM.epETdV) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + (2*pListSM.epETuV + pListSM.epETdV) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	            }
	            case 2:
	            {
	                //SM + mm NSI terms
                    pListSM.qA = 4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (GVN + pListSM.epMMuV + 2*pListSM.epMMdV ) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + ( GVP + 2*pListSM.epMMuV + pListSM.epMMdV ) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	                
	                //em NSI terms
	                pListSM.qA = 0;//4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (pListSM.epEMuV + 2*pListSM.epEMdV) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + (2*pListSM.epEMuV + pListSM.epEMdV) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	                
	                //mt NSI terms
	                pListSM.qA = 0;//4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (pListSM.epMTuV + 2*pListSM.epMTdV) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + (2*pListSM.epMTuV + pListSM.epMTdV) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	            }
	            case 3:
	            {
	                //SM + tt NSI terms
                    pListSM.qA = 4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (GVN + pListSM.epTTuV + 2*pListSM.epTTdV ) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + ( GVP + 2*pListSM.epTTuV + pListSM.epTTdV ) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	                
	                //mt NSI terms
	                pListSM.qA = 0;//4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (pListSM.epMTuV + 2*pListSM.epMTdV) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + (2*pListSM.epMTuV + pListSM.epMTdV) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	                
	                //et NSI terms
	                pListSM.qA = 0;//4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
	                pListSM.qV = ( (pListSM.epETuV + 2*pListSM.epETdV) * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + (2*pListSM.epETuV + pListSM.epETdV) * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
	                rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], sourcej, fluxj);
	            }
	        }
	     }	
	    //if(pList->elecScat)
	    {
	        int Ne=0;
	        while(pList->detectors[detj].ionization[i][Ne] > ErKeV && Ne < pList->detectors[detj].isoZ[i]) 
	            Ne++;
	    
	        if(pList->sources[sourcej].isSolar[fluxj] == 1)
	        {
	            //nu_e component
		        pListSM.qA = 0.5;
		        pListSM.qV = 0.5+2*SSW;
		        rate += pList->sources[sourcej].survProb[fluxj] * (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, sourcej, fluxj);
                //nu_mu and nu_tau component
		        pListSM.qA = -0.5;
		        pListSM.qV = -0.5+2*SSW;
		        rate += (1-pList->sources[sourcej].survProb[fluxj]) * (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, sourcej, fluxj);
		    }
		    else
		    {
		        switch( abs(pList->sources[sourcej].nuFluxFlav[fluxj]) )
	            {
    	            case 1:
    	            {
        		        pListSM.qA = 0.5*anti+pListSM.epEEeA;
        		        pListSM.qV = 0.5+2*SSW+pListSM.epEEeV;
        		        rate += (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, sourcej, fluxj);
                    }
                    case 2:
                    {
                        pListSM.qA = -0.5*anti+pListSM.epMMeA;
        		        pListSM.qV = -0.5+2*SSW+pListSM.epMMeV;
        		        rate += (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, sourcej, fluxj);
                    }
                    case 3:
                    {
                        pListSM.qA = -0.5*anti+pListSM.epTTeA;
        		        pListSM.qV = -0.5+2*SSW+pListSM.epTTeV;
        		        rate += (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, sourcej, fluxj);
                    }
                }
		    }
		    
	    }
		
    }
    
	return rate*detEff(ErKeV,pList->detectors[detj].eff); 
}

//below are functions for total rate of all fluxes
double diffNSIrate(double ErkeV, paramList *pList, int detj)			  
{   
    double rate=1e-99;
    for(int fluxi=0; fluxi< pList->sources[pList->detectors[pList->detj].sourcej].numFlux; fluxi++)
    {
        if( ErkeV < (pList->detectors[detj].ErL + (double)999*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/900) )
           rate += NSIrate( ErkeV, pList, detj, pList->detectors[pList->detj].sourcej, fluxi);
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

