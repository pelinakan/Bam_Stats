//
//  MapQ_stats.h
//  BAM_QC
//
//  Created by Pelin Sahlen on 07/06/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//

#ifndef BAM_QC_MapQ_stats_h
#define BAM_QC_MapQ_stats_h


int main(){
    
    const double bucket_size = 0.05;
    double area=0;
    int number_of_buckets = (int)ceil(1 / bucket_size);
    std::vector<int> histogram(number_of_buckets);
    ifstream theFile("Datafile.txt");
    
    
    while (theFile >> area) {
        int bucket = (int)floor(area / bucket_size);
        histogram[bucket] += 1;
    }
    for(auto loop = histogram.begin(); loop != histogram.end();++loop)
    {
        
        cout<<bucket_size<<"\t"<<number_of_buckets<<endl;
        
    }
    
    return 0;
}


#endif
