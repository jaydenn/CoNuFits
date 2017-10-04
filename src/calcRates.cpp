#include <iostream>
#include <fstream>
#include <iomanip>
#include "SMrate.h"
#include "BSMrate.h"
#include "NSIrate.h"
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int calcRates(paramList *pList)
{
    double ErkeV;
    char filename[90];    
    std::ofstream outfile;
    char BSMname[30];
    int masskeV = (int)(pList->mMed*1e6);    
    
    switch(pList->BSM)
    {

        case 1:
        {
            sprintf(BSMname,"%dkeVscalar",masskeV);
            break;
        }
        case 2:
        {
            sprintf(BSMname,"%dpseudoS",masskeV);
            break;
        }
        case 3:
        {
            sprintf(BSMname,"%dvector",masskeV);
            break;
        }
        case 4:
        {
            sprintf(BSMname,"%daxial",masskeV);
            break;
        }
        case 5:
        {
            sprintf(BSMname,"sterile");
            break;
        }
        case 6:
        {
            sprintf(BSMname,"NSI");
            break;
        }

    }
    
    //format output streams
    std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    outfile   << std::setiosflags(std::ios::scientific) << std::setprecision(5);
    
    //output model
    if( pList->BSM < 5 )    
	    std::cout << "\nBSM rate for " << masskeV << " keV "<< BSMname << " mediator\n";
    char target[5];
    if(pList->nucScat)
        sprintf(target,"%s","N");
    else
        sprintf(target,"%s","E");
        
    for(int detj=0; detj < pList->ndet; detj++)
    {
        //open file
        sprintf(filename,"%s%s_rate%s_%s_%s.dat", pList->root, BSMname, target, pList->detectors[detj].name, pList->sources[pList->detectors[detj].sourcej].name);
        outfile.open(filename,std::ios::out);
        
        if( !outfile )
        {
            std::cout << "output file could not be created (does the directory in specified root path exist?" << std::endl;
            return 1;
        }
        
        std::cout << "------------------------\n";
        std::cout << detj+1 << ". " << pList->detectors[detj].name << std::endl;
        std::cout << "------------------------\n";
        std::cout << "  total rates: \n"; 
        std::cout << "     SM  = " << intSMrate(  pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj)         << " events/kg/day" << std::endl;
        std::cout << "     BG  = " << intBgRate(  pList->detectors[detj], pList->detectors[detj].ErL, pList->detectors[0].ErU) << " events/kg/day" << std::endl;
        if ( pList->BSM < 5 )
            std::cout << "     BSM = " << intBSMrate( pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj,1)       << " events/kg/day" << std::endl;
        else if ( pList->BSM == 6 )
            std::cout << "     NSI = " << intNSIrate( pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj)       << " events/kg/day" << std::endl;
        
        std::cout << "  differential rates: \n"; 
        std::cout << "    Er (keV)        SM dN/dE        BG dN/dE        BSM dN/dE (events/kg/day/keV)" << std::endl;
        outfile   << "    Er (keV)        SM dN/dE        BG dN/dE        BSM dN/dE (events/kg/day/keV)" << std::endl;

        int skip=0;          
        double BSM;  
        for (int i=0; i<501; i+=1)
        {
            if(pList->logBins == 1)
                ErkeV = pow(10, log10(pList->detectors[detj].ErL) + (double)i*(log10(pList->detectors[detj].ErU)-log10(pList->detectors[detj].ErL))/500)+1e-4;
            else
                ErkeV = pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/500;

            if ( pList->BSM < 5 )
                BSM = diffBSMrate( ErkeV, pList, detj, 1);
            else if ( pList->BSM == 6 )
                BSM = diffNSIrate( ErkeV, pList, detj);

            if( skip++ % 5 == 0 ) //only print out every fifth value to terminal    
                std::cout << "    " << ErkeV << "      " << diffSMrate( ErkeV, pList, detj) << "      " << diffBgRate( pList->detectors[detj], ErkeV) << "      " << BSM << std::endl;
            
            outfile   << "    " << ErkeV << "      " << diffSMrate( ErkeV, pList, detj) << "      " << diffBgRate( pList->detectors[detj], ErkeV) << "      " << BSM << std::endl; 
        }
        
        outfile.close();
    }
    
}    

