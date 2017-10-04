#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include "detectorFunctions.h"
#include "nuFlux.h"
#include "SMrate.h"
#include "BSMrate.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	
#ifndef SOURCESTRUCT_H
	#include "sourceStruct.h"
#endif	

//gets sampling parameters from file
int readConfigFile(paramList *pL, char *filename) 
{

    FILE* input;
    input = fopen(filename,"r");
    if(input==NULL) 
    {
        std::cout << "unable to open parameter file: " << filename << std::endl;
        return -1;	
    }
  
    int mode;
    char *ret;
    char temp[400];
    char root[50];
    ret = fgets(temp,200,input);

    //Mode switch
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&mode);
    
    if ( mode == 2 )
    {
        //Get Multinest sampling parameters
        FILE* input;
        input = fopen("multinest.ini","r");
        if(input==NULL) 
        {
            std::cout << "unable to open multinest.ini\n";
            return -1;	
        }
        ret = fgets(temp,200,input);
        for (int i=0;i<8;i++)
        {
            ret = fgets(temp,200,input);
            sscanf(temp,"%lf",&(pL->sampling[i]));
        }
    }
    
    // root for output files
    ret = fgets(temp,200,input);
    sscanf(temp,"%s %*s",root);
    sprintf(pL->root, "%s", root);		 
  
    //include nuclear scattering?
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->nucScat));
    
    //include electron scattering?
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->elecScat));
    
    if( (pL->nucScat == 0 && pL->elecScat == 0) || (pL->nucScat == 1 && pL->elecScat == 1))
    {
        std::cout << "Must choose one of electron or nuclear scattering" << std::endl;
        return -1;
    }
    
    //which BSM model to consider
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->BSM));
    
    //fiducial coupling
    ret = fgets(temp,200,input);
    sscanf(temp,"%lf",&(pL->C));
    
    switch(pL->BSM)
    {
    
        case 0:
        {
            break;
        }
        case 1:
        {
            pL->gNuS=1;
            switch(pL->nucScat)
            {
                case 1:
                {
                    pL->qNs=pL->qPs=pL->C;
                    break;
                }    
                case 2:
                {
                    pL->qNs=pL->C;
                    break;
                }
                case 3:
                {
                    pL->qPs=pL->C;
                    break;
                }    
            }
            if(pL->elecScat)
            {
                pL->gEs=pL->C;
            }
            break;
        }
        case 2:
        {
            pL->gNuS=1;
            if(pL->elecScat)
            {
                pL->gEp=pL->C;
                break;
            }
            else
            {
                std::cout << "there is no pseudo-scalar scattering in the nuclear case, exiting.\n";
                return -1;
            }
        }
        case 3:
        {
            pL->gNuV=1;
            switch(pL->nucScat)
            {
                case 1:
                {
                    pL->qNv=pL->qPv=pL->C;
                    break;
                }    
                case 2:
                {
                    pL->qNv=pL->C;
                    break;
                }
                case 3:
                {
                    pL->qPv=pL->C;
                    break;
                }    
            }
            if(pL->elecScat)
            {
                pL->gEv=pL->C;
            }
            break;
        }
        case 4:
        {
            pL->gNuV=1;
            switch(pL->nucScat)
            {
                case 1:
                {
                    pL->qNa=pL->qPa=pL->C;
                    break;
                }    
                case 2:
                {
                    pL->qNa=pL->C;
                    break;
                }
                case 3:
                {
                    pL->qPa=pL->C;
                    break;
                }    
            }
            if(pL->elecScat)
            {
                pL->gEa=pL->C;
            }
            break;
        }
        case 6:
        {
           //pL->epMMuV = pL->epMMdV = pL->epEEuV = pL->epEEdV = pL->C;
           break;
        }
        default:
        {
            std::cout << "Must choose a BSM type" << std::endl;
            return -1;
        }
    }
    
    //mediator mass
    ret = fgets(temp,200,input);
    sscanf(temp,"%lf",&(pL->mMed));
    
    //log or linear spaced bins?
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->logBins));
    
    //Detector setup
    double exp;
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);

    char name[30];
    char sourceName[30]="";
    int sourcej;
    double sourceDistance=0;
    
    while(temp[0]=='#')
    {
        sscanf(temp,"%s %lf %s %lf", name, &exp, sourceName, &(sourceDistance));
        sourcej=sourceInit(pL, sourceName, sourceDistance);
        if(sourcej<0) { std::cout << "Could not create source, exiting" << std::endl; return -1; }
        
	    if(newDetector(pL, name, exp, sourcej)) { std::cout << "Could not create detector, exiting" << std::endl; return -1; }
  	    
        std::cout << "Using detector " << name << " with source " << sourceName << std::endl;
        ret = fgets(temp,200,input);
    }    
    
    //asimov or random sim?
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->asimov));
    
    fclose(input); 

    return mode;
}

