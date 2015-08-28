//
//  GetReadStats.h
//  BAM_QC
//
//  Created by Pelin Sahlen on 03/06/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//

#ifndef BAM_QC_GetReadStats_h
#define BAM_QC_GetReadStats_h

class ReadStats{
public:
    long int NumberofPairs;
    long int NumofDuplicates;
    long int NumofMappedReads;
    long int NumofMappedPairs;
    long int NumofUnmappedReads;
    long int NumofUniqMappedPairs;
    long int NumofPrimaryAlignments;
    long int NumofReadsFailedQC;
    long int NumofProperPairs;
    long int NumofPairs_ReadsOverlapping;
    
    vector < double > InsertSizes;
    int max_insert_size;
    
    int mapq_bin_size;
    int max_mapq;
    vector <int> mapq_histogram;
    
    void InitialiseClass();
    void PrintStats(string);
    void PrintInsertSizes(string);
private:
    void ExtractStats(string);
    
};
void ReadStats::InitialiseClass(){

    NumberofPairs = 0;
    NumofDuplicates = 0;
    NumofMappedReads = 0;
    NumofMappedPairs = 0;
    NumofUnmappedReads = 0;
    NumofUniqMappedPairs = 0;
    NumofPrimaryAlignments = 0;
    NumofReadsFailedQC = 0;
    NumofPairs_ReadsOverlapping = 0;
    NumofProperPairs = 0;
    
    max_insert_size = 1000;
    mapq_bin_size = 1;
    max_mapq = 70;
    

}

void ReadStats::PrintStats(string bamfilename){
    
    string bam_out = bamfilename.substr(0,(bamfilename.length() - 3));
    string mapq_out = bamfilename.substr(0,(bamfilename.length() - 3));
    bam_out.append("stats.txt");
    mapq_out.append("mapq.txt");
    
    ofstream ofile(bam_out.c_str());
    ofstream mapqfile(mapq_out.c_str());
    
    ofile << "Number of Pairs" << '\t' << "Number of Uniquely Mapped Pairs (mapq > 0)" << '\t'
          << "Percentage of Uniquely Mapped Pairs" << '\t' << "Number of Duplicates" << '\t' << "Percentage of Duplicates" << endl;
    ofile << NumberofPairs << '\t' << NumofMappedReads << '\t' << double(NumofMappedReads)/double(NumberofPairs)*100.0
          << '\t' << NumofDuplicates << '\t'<< double(NumofDuplicates)/double(NumberofPairs)*100.0 << endl;
    
    cout << "Mapping stats printed" << endl;
    for (int i = 0; i < mapq_histogram.size(); ++i)
        mapqfile << i*mapq_bin_size << '\t' << mapq_histogram[i] << '\t' << endl;
    cout << "Mapping quality histogram printed" << endl;
}

void ReadStats::PrintInsertSizes(string bamfilename){
    
    alglib::real_1d_array AX;
    AX.setcontent(InsertSizes.size(), &(InsertSizes[0]));
    double mean = alglib::samplemean(AX);
    double variance = alglib::samplevariance(AX);
    double skew = alglib::sampleskewness(AX);
    double kurtosis = alglib::samplekurtosis(AX);
    double median = 0;
    samplemedian(AX,median);
    
    int size_of_bins = 1;
    int number_of_bins = (int) ceil(max_insert_size / size_of_bins);
	cout << "Num of bins "  << number_of_bins << "  " << max_insert_size << "  " << size_of_bins << endl; 
    vector <long unsigned int> histogram(number_of_bins);
    for (int i = 0; i < number_of_bins; ++i)
        histogram[i] = 0;

	cout << "here 0"  << endl;     
    
    for(int i = 0; i < InsertSizes.size(); ++i){
        int bucket = (int)floor(InsertSizes[i] / size_of_bins);
        if(bucket >= histogram.size())
            histogram.resize(bucket+1,0);
        histogram[bucket] += 1;
    }
    
    string temp = bamfilename.substr(0,bamfilename.size()-3);
    temp.append("InsertSizeHistogram.txt");
    ofstream outf(temp.c_str());
    
    for (int i = 0; i < histogram.size(); ++i){
        outf << i*size_of_bins << '\t' << histogram[i] << endl;
    }
    
    cout << "Insert size histogram printed" << endl;
    
    string temp2 = bamfilename.substr(0,bamfilename.size()-3);
    temp2.append("InsertSize.stats.txt");
    ofstream outfile(temp2.c_str());
    
    outfile << "Number of pairs processed" << '\t' << "MEAN" << '\t' << "VARIANCE" << '\t' << "SKEW" << '\t' << "KURTOSIS" << '\t' << "MEDIAN" << endl;
    outfile << NumberofPairs << '\t' << mean << '\t' << variance << '\t' << skew << '\t' << kurtosis << '\t' << median << endl;
    
    cout << "Insert size stats are printed" << endl;
    
}
void ReadStats::ExtractStats(string bamfilename){

    BamReader reader;
    //	BamWriter writer;
    
    bamfilename = "/Users/pelinakan/Documents/WORK/GA_PROJECTS/BAM_QC/BAM_QC/test22.bam";
    string bam_out = bamfilename.substr(0,(bamfilename.length() - 3));
    string mapq_out = bamfilename.substr(0,(bamfilename.length() - 3));
    
    cout << bamfilename << endl;
    
    if ( !reader.Open(bamfilename.c_str()) )
        cerr << "Could not open input BAM files." << endl;
    
    cout << bamfilename << "  opened, getting the read stats " << endl;
    const SamHeader sam_header = reader.GetHeader();
    const RefVector bam_ref_data = reader.GetReferenceData();
    
    bam_out.append("stats.txt");
    mapq_out.append("mapq.txt");
    
    ofstream ofile(bam_out.c_str());
    ofstream mapqfile(mapq_out.c_str());
    BamAlignment al;
    // ofile << "Read Name" << '\t' << "AS" << '\t' << "XS" << '\t' << "Mapping Quality" << '\t' << "Is Duplicate" << '\t' <<  "Is Primary Alignment" << endl;
    ofile << "Number of Reads" << '\t' << "Number of Paired Reads" << '\t' << "Number of Proper Pairs" << '\t' << "Number of Mapped Reads" << '\t' << "Number of Unmapped Reads" << '\t' << "Number of Mapped Pairs"  << '\t' << "Number of Uniquely Mapped Reads" << '\t' << "Number of Duplicates" << endl;
    
    cout << "Reading BAM file..." << endl;
    int xs, as;
    string sa; //chimeric alignment tag
    while(reader.GetNextAlignment(al)){
        ++NumberofPairs;
        al.GetTag("XS", xs);
        al.GetTag("AS", as);
        al.GetTag("SA", sa);
        if (al.IsDuplicate())
            ++NumofDuplicates;
        if (!al.IsMapped())
            ofile << al.Name << endl;
        if(al.IsMapped())
            ++NumofMappedReads;
        if (al.IsMapped() && al.IsMateMapped())
            ++NumofMappedPairs;
        if (al.IsPrimaryAlignment())
            ++NumofPrimaryAlignments;
        if (al.IsFailedQC())
            ++NumofReadsFailedQC;
        if (al.IsProperPair())
            ++NumofProperPairs;
        if (al.MapQuality > 0 ) //BWA-MEM reports only one aligment for mapq>0 but some reads align chimerically (breakpoints etc) and chimeric alignments are written in SA tag.
            ++NumofUniqMappedPairs;
        int bucket = (int)floor(al.MapQuality / mapq_bin_size);
        mapq_histogram[bucket] += 1;
        //ofile << al.Name << '\t' << as << '\t' << xs << '\t' << al.MapQuality << '\t' << al.IsDuplicate()  << '\t' << al.IsPrimaryAlignment() << '\t' << endl;
        if(NumberofPairs == 5000000)
            break;
    }
    ofile << NumberofPairs << '\t' << NumofProperPairs << '\t' <<  NumofMappedReads << '\t' <<(NumberofPairs - NumofMappedReads) << '\t' << NumofMappedPairs << '\t' << NumofUniqMappedPairs << '\t' << NumofDuplicates << endl;
    for (int i = 0; i < mapq_histogram.size(); ++i)
        mapqfile << i*mapq_bin_size << '\t' << mapq_histogram[i] << '\t' << endl;
        
        
        
    cout << NumberofPairs << "    Finished" << endl;
    
}


#endif
