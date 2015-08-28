//
//  PrintInsertSizes.cpp
//  This takes insert sizes and select a subsample of reads with a different insert size
//
//  Created by Pelin Sahlen on 09/05/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cstring>
#include <string>
#include <time.h>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
using namespace std;

//ALGLIB PACKAGE HEADERS
#include "alglibmisc.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "statistics.h"
#include "dataanalysis.h"
#include "specialfunctions.h"
#include "solvers.h"
#include "optimization.h"
#include "diffequations.h"
#include "fasttransforms.h"
#include "integration.h"
#include "interpolation.h"
using namespace alglib_impl;

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

const int N_THREADS = 15;
#include <pthread.h>

pthread_mutex_t poolMutex;
std::vector< BamTools::BamAlignment > pool;
bool die = false;


const int max_coverage = 60;

#include "GetReadStats.h"
#include "computeCoverage.h"
#include "thread_functions.h"

#include "radix.h"
#include "sam.h"


#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>


#include "preseq.cpp"

#include "print_usage.h"

int main(int argc, const char * argv[]) {
//int main(void){

    if (argc < 2) {
        print_usage();
        return -1;
    }
 
    ReadStats stats; //Class that keeps the read stats
    stats.InitialiseClass();
    
	string bamfilename="";
//    bamfilename = "/Users/pelinakan/Documents/WORK/GA_PROJECTS/BAM_QC/BAM_QC/test22.bam";
	bamfilename = argv[1];
    
    BamReader reader;
    CoverageClass cov;
    
    
    cov.Initialise(bamfilename, stats, reader);
    
    // Initialise a mutex thread
    pthread_mutex_init(&poolMutex,NULL);
    //Startup generating threads!
    cov.PickAlignments(bamfilename, stats);
    cov.CalculateCoverage();
    
    pthread_mutex_lock(&poolMutex);
    pool.clear();
    pthread_mutex_unlock(&poolMutex);
    
    cov.GenerateHistogram_Global(bamfilename, "coverage");
//    cov.ComputeCoverage(bamfilename, stats);
    
    cout << "Coverages Generated, Printing read and insert size stats... "<< endl;

    //stats.PrintStats(bamfilename);
    //stats.PrintInsertSizes(bamfilename);
    
    cout << "BAM_QC completed on " << bamfilename << endl;
    
//    lc_extrap("yield.out", bamfilename);
    
    return 0;
	
    //    qaCompute(bamfilename);
 //stats.ExtractStats(bamfilename);

}
