#include <iostream>
#include <cmath>
#ifndef GSL_INTERP_H
	#include <gsl/gsl_interp.h>
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	
#ifdef MPI
    #include "mpi.h"
#endif
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "nuFlux.h"
#include "multinest.h"
#include "likelihood.h"
#include "monteCarlo.h"

void dumper(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void **) {}


void NSIglobalFit(paramList *pL)
{
        
        // set some MultiNest pL->sampling parameters
        int pWrap[12];              // which parameters to have periodic boundary conditions?
        int seed = -1;			   	            // random no. generator seed for MultiNest, if < 0 then take the seed from system clock
        int simSeed = 0;
        double logZero = -DBL_MAX;							  // points with loglike < logZero will be ignored by MultiNest
        int initMPI = 0;								      // initialize MPI routines?, relevant only if compiling with MPI
        int outfile = 1;								      // write output files?
        int updateInt = 1000000;							  // update interval (for calling dumper, which isn't used here)

        int ndims = 4+pL->ndet+pL->nSource;                                // number of parameters in the reconstruction

        int npar  = ndims;                                // global number of parameters

        int nEvents = 0;                                         
        int myrank=0;   
        
        //setup random number generator
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T=gsl_rng_default;
        r=gsl_rng_alloc(T);
        gsl_rng_set(r, simSeed);     
        
        //Initialise MPI
        #ifdef MPI
            int argc=0;
            char **argv;
            MPI_Init(&argc, &argv);
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	        if(myrank==0)
	        {
	            gsl_rng_env_setup();
        		T=gsl_rng_default;
	            r=gsl_rng_alloc(T);
	            gsl_rng_set(r, (int)time(NULL));
	            simSeed = (int) 60000*gsl_rng_uniform(r);		//this is to make sure each MPI thread has the same 'random' data (it would be nice to have more than 60,000 different experiments.. the limit is because it only takes an integer as a seed)
	        }
	        MPI_Bcast(&simSeed,1,MPI_INT,0,MPI_COMM_WORLD);
	        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        void *pointer;                                     //a pointer for passing paramList to Multinest
        pointer = (void *)  pL;
        
        //Generate sim data for each detector - right now only sim SM 
        if(myrank==0)
            std::cout << "Using " << pL->ndet << " detector(s):" << std::endl;
    	for(int detj=0; detj<pL->ndet; detj++)
        {
            
            nEvents = generateBinnedDataSM( pL, detj, 1, simSeed);
            if( strcmp(pL->detectors[detj].name, "CsI") == 0 )
            {
                std::cout << "over-writing CsI observed SM events to 134\n";
                pL->detectors[detj].binnedData[0] -= 40;
            }
       	    if(nEvents<0)
       	    {
       	        std::cout << "error in MC\n";
       	        return;
            }
            else if(myrank==0)
    	        std::cout << "  " << pL->detectors[detj].name << " (" << pL->detectors[detj].exposure << " t.y): " << pL->detectors[detj].nEvents << " events" << std::endl;
        }
        //run multinest sampling
        if(myrank==0) std::cout << "Starting MultiNest sampling..." << std::endl;
        //  nestRun(          mmodal,            ceff,           nlive,              tol,              efr,ndims, nPar,nCdims,  maxModes,    updInt,           nullZ,      root, seed, pWrap,        feedback,          resume, outfile,       initMPI, logZero,                loglike, dumper, context)
        nested::run( pL->sampling[0], pL->sampling[1], pL->sampling[2],  pL->sampling[6],  pL->sampling[5],ndims, npar, ndims,       100, updateInt, pL->sampling[7],  pL->root, seed, pWrap, pL->sampling[3], pL->sampling[4],       1, (bool)initMPI, logZero, logLikelihoodGlobalFit, dumper, pointer);

        #ifdef MPI
            MPI_Finalize();
        #endif
}
