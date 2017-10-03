//definition of parameter structs
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <string.h>
#ifndef IOSTREAM
	#include <iostream>
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	
#ifndef SOURCESTRUCT_H
    #include "sourceStruct.h"
#endif  
#define PARAMETERSTRUCT_H

struct paramList {

	//root directory for file output
	char root[50];

	//Struct for neutrino sources
	int nSource, sourcej, fluxj;
    sourceStruct sources[10];
	double Er;
    
	//setup for rate integration
	gsl_function F;
	double (*rateFunc)(double, double, paramList *, int);

	double A; //is this used?
	double Qa, qA, qAn, qAp, qAu, qAd;
	double Qv, qV, qVn, qVp, qVu, qVd;

	//initial coupling
	double C, signalNorm;
	double SMinterference1,SMinterference2; //so it can be turned off when needed (in initialization)
	
	//nucleon couplings
	double Qs;  //scalar nuclear
	double qPs; //scalar
	double qNs; //scalar
	double Qvp;  //vector nuclear
    double qPv; //vector
    double qNv; //vector
    double Qap;  //axial nuclear
    double qPa; //axial
    double qNa; //axial
    
    //neutrino couplings
    double gNuV;  //vector neutrino
    double gNuS;  //scalar neutrino
        
    //electron couplings	
	double gEs;   //scalar electron
    double gEp;   //pseudoscalar electron
    double gEv;  //vector electron
    double gEa;  //axial-vector electron
    
    //BSM model parameters
	int BSM;
	int nucScat;
	int elecScat;
    int mediator;
    double mMed;
    double delMsqGeV;
    double ss2Theta14;
    
    //NSI parameters
    double epEEuV, epEEdV, epEEuA, epEEdA, epEEeV, epEEeA;
    double epEMuV, epEMdV, epEMuA, epEMdA, epEMeV, epEMeA;
    double epETuV, epETdV, epETuA, epETdA, epETeV, epETeA;
    double epMMuV, epMMdV, epMMuA, epMMdA, epMMeV, epMMeA;
    double epMTuV, epMTdV, epMTuA, epMTdA, epMTeV, epMTeA;
    double epTTuV, epTTdV, epTTuA, epTTdA, epTTeV, epTTeA;
    
    double maxL; 
    int asimov;
    int logBins;
    int maxBins;
    
	int ndet, detj;
	detector detectors[10];
	
	void printPars()
	{
   		std::cout << "Conudl configuration:" << std::endl;		
    	std::cout << "  root: " << root << std::endl;
    	std::cout << "  BSM: " << BSM << std::endl;
    	std::cout << "  asimov: " << asimov << std::endl;
    	std::cout << "  logBins: " << logBins << std::endl;
    	//std::cout << "  flux: " << nuFlux << " +/- " << nuFluxUn*100 << "%" << std::endl;
		std::cout << ndet << " Detector(s):" << std::endl;
		for(int i=0;i<ndet;i++)
			detectors[i].printDetSpecs();
	    std::cout << ndet << " Source(s):" << std::endl;
	    for(int i=0;i<nSource;i++)
			sources[i].printSourceSpecs();
    	std::cout << "  Mediator mass = " << mMed << std::endl;
        std::cout << "  Couplings: " << std::endl;    
        std::cout << "       electron: " << std::endl; 
        std::cout << "             gEs = " << gEs << std::endl; //Qv=Qa=Qs=Qvp=Qap=
        std::cout << "             gEp = " << gEp << std::endl; 
        std::cout << "             gEv = " << gEv << std::endl; 
        std::cout << "             gEa = " << gEa << std::endl; 
        std::cout << "       neutron:      " << std::endl; 
        std::cout << "             qNs = " << qNs << std::endl; 
        std::cout << "             qNv = " << qNv << std::endl; 
        std::cout << "             qNa = " << qNa << std::endl; 
        std::cout << "       proton:     " << std::endl; 
        std::cout << "             qPs = " << qPs << std::endl; 
        std::cout << "             qPv = " << qPv << std::endl; 
        std::cout << "             qPa = " << qPa << std::endl; 
        std::cout << "       neutrino:     " << std::endl; 
        std::cout << "             gNuS = " << gNuS << std::endl; 
        std::cout << "             gNuV = " << gNuV << std::endl; 
    }
	
	paramList()
	{
	    Qv=Qa=Qs=Qvp=Qap=qNs=qPs=qNv=qPv=qNa=qPa=gNuS=gNuV=gEs=gEp=gEv=gEa=0;
	    epEEuV=epEEdV=epEEuA=epEEdA=epEEeV=epEEeA=0;
	    epEMuV=epEMdV=epEMuA=epEMdA=epEMeV=epEMeA=0;
	    epETuV=epETdV=epETuA=epETdA=epETeV=epETeA=0;
	    
	    epMMuV=epMMdV=epMMuA=epMMdA=epMMeV=epMMeA=0;
	    epMTuV=epMTdV=epMTuA=epMTdA=epMTeV=epMTeA=0;
    	    	    
	    epTTuV=epTTdV=epTTuA=epTTdA=epTTeV=epTTeA=0;
	    
	    C=Er=mMed=maxL=A=0;
	    BSM=elecScat=nucScat=mediator=0;
	    signalNorm=SMinterference1=SMinterference2=1;
		ndet=detj=fluxj=nSource=sourcej=0;
        asimov=logBins=1;
        maxBins = MAXBINS;
	}
	
};
	
