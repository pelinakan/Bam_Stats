//
//  computeCoverage.h
//  BAM_QC
//
//  Created by Pelin Sahlen on 02/07/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//
#ifndef BAM_QC_computeCoverage_h
#define BAM_QC_computeCoverage_h

struct ChrStruct{
    int *positionmap; // map->first : base position, map->second : how many times it is covered
    string chrname;
    int chr_length;
};
class CoverageClass{
    friend class ReadStats;
public:
    boost::unordered::unordered_map<int, string> chrmap;
    void ComputeCoverage(string, ReadStats&);
private:
    void AllocateMemorytoCoverageVectors(vector<ChrStruct>&, int, int,string);
    bool fillcoveragevector_span(BamAlignment, vector<ChrStruct>&);  //span, double count overlapping reads, counting_overlapping_bases_once count overlapping reads (base coverage)
    bool fillcoveragevector_countonce_overlappingbases(BamAlignment, vector<ChrStruct>&);  //span, double count overlapping reads, counting_overlapping_bases_once countoverlapping reads (base coverage)
    bool fillcoveragevector_counttwice_overlappingbases(BamAlignment, vector<ChrStruct>&);  //span, double count overlapping reads, counting_overlapping_bases_once count overlapping reads (base coverage)

    void GenerateHistogram_perChr(string, int*, int,string);
    void GenerateHistogram_Global(string, string, vector<ChrStruct>);
};

void CoverageClass::AllocateMemorytoCoverageVectors(vector<ChrStruct> & CoverageperChr, int key, int RefLength,string RefName){
    
    CoverageperChr[key].positionmap = (int*) realloc(CoverageperChr[key].positionmap, (RefLength + 1)*sizeof(int));
    if (CoverageperChr[key].positionmap == NULL)
        cerr << "Could not allocate memory to chr array" << endl;
    
    memset(CoverageperChr[key].positionmap, 0, (RefLength + 1)*sizeof(int));
    CoverageperChr[key].chr_length = RefLength;
    CoverageperChr[key].chrname = RefName;
}

bool CoverageClass::fillcoveragevector_span(BamTools::BamAlignment al, vector<ChrStruct> & CoverageperChr){
    
    if (al.Position < al.MatePosition) { // Reads do not overlap, increase coverage for the bases covered by this alignment
        for (int i = al.Position; i < al.MatePosition + al.Length; ++i)
            CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its coverage by one
    }
    else { // Reads overlap !! Only count the non-overlapping bits
        for (int i = al.Position; i < al.Position + al.Length; ++i)
            CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its coverage by one
        int remainingbases = (al.MatePosition + al.Length) - (al.Position + al.Length);
        for (int i = al.Position + al.Length; i < al.Position + al.Length + remainingbases; ++i)
            CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its
    }
    return 0;
}
bool CoverageClass::fillcoveragevector_countonce_overlappingbases(BamAlignment al, vector<ChrStruct>& CoverageperChr){

    if (al.Position < al.MatePosition) { // Reads do not overlap, increase coverage for the bases covered by this alignment
        for (int i = al.Position; i < al.Position + al.Length; ++i)
            CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its coverage by one
        for (int i = al.MatePosition; i < al.MatePosition + al.Length; ++i)
            CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its coverage by one
    }
    else { // Reads overlap !! Only count the non-overlapping bits
        for (int i = al.Position; i < al.MatePosition + al.Length; ++i)
            CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its coverage by one
    }
    return 0;

}
bool CoverageClass::fillcoveragevector_counttwice_overlappingbases(BamAlignment al, vector<ChrStruct>& CoverageperChr){

    for (int i = al.Position; i < al.Position + al.Length; ++i){
       // cout << al.Position << "  " << al.Length << "  " << al.MatePosition << endl;
        CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its coverage by one
    }
    for (int i = al.MatePosition; i < al.MatePosition + al.Length; ++i)
        CoverageperChr[al.RefID].positionmap[i] += 1; // Find the position and increment its coverage by one
    return 0;
}
void CoverageClass::ComputeCoverage(string bamfilename, ReadStats& stats){
    
    vector< ChrStruct > CoverageperChr_span; // Each element keeps the coverage count of each base per chromosome
    vector< ChrStruct > CoverageperChr_nooverlapcount;
    vector< ChrStruct > CoverageperChr_overlapcount;
    
    BamReader reader;
    cout << bamfilename << endl;
    
    
    int number_of_bins = (int) ceil(stats.max_mapq / stats.mapq_bin_size);
    stats.mapq_histogram.resize(number_of_bins);
    
    if ( !reader.Open(bamfilename.c_str()) )
        cerr << "Could not open input BAM file." << endl;
    
    cout << bamfilename << "  opened" << endl;
    
    const RefVector bam_ref_data = reader.GetReferenceData();

    // Make a map of chr names to RefIDs
    RefVector::const_iterator chrit;
    
    
    
    CoverageperChr_span.resize(bam_ref_data.size()); // Create one for each chromosome
    CoverageperChr_nooverlapcount.resize(bam_ref_data.size()); // Create one for each chromosome
    CoverageperChr_overlapcount.resize(bam_ref_data.size()); // Create one for each chromosome

    for (chrit = bam_ref_data.begin(); chrit != bam_ref_data.end(); ++chrit){
        int key = reader.GetReferenceID(chrit->RefName);
        chrmap[key] = chrit->RefName;
 //       cout << chrit->RefName << '\t' << chrit->RefLength << '\t' << key << '\t' << chrmap[key] << endl;
        AllocateMemorytoCoverageVectors(CoverageperChr_span,key, chrit->RefLength, chrit->RefName);
        AllocateMemorytoCoverageVectors(CoverageperChr_nooverlapcount,key, chrit->RefLength, chrit->RefName);
        AllocateMemorytoCoverageVectors(CoverageperChr_overlapcount,key, chrit->RefLength, chrit->RefName);
    }

    BamAlignment al;
    cout << "Generating coverage" << endl;
    
    while(reader.GetNextAlignmentCore(al)){
        ++stats.NumberofPairs;
        if(!al.IsDuplicate()){ // Do not count if it is a duplicate !!
            if (al.IsMapped() && al.IsMateMapped()) { // if both mapped
                if (al.MapQuality > 0 ) { // Do not count if mapping quality is zero !!
                    ++stats.NumofUniqMappedPairs;
                    fillcoveragevector_span(al,CoverageperChr_span);
                    fillcoveragevector_countonce_overlappingbases(al, CoverageperChr_nooverlapcount);
                    fillcoveragevector_counttwice_overlappingbases(al, CoverageperChr_overlapcount);
                    if (al.RefID == al.MateRefID) { // Do not count if they are on different chromosomes !!
//                        if (al.IsProperPair()) { // Do not count if they are not in valid orientation or invalid insert size, This put hard boundries in the insert size that are too restrictive so it is not used for now.
                        ++stats.NumofProperPairs;
                        stats.InsertSizes.push_back(al.InsertSize);
                        int bucket = (int)floor(al.MapQuality / stats.mapq_bin_size);
                        if (bucket >= stats.mapq_histogram.size())
                            stats.mapq_histogram.resize(bucket+1,0);
                        stats.mapq_histogram[bucket] += 1;
                    }
                }
            }
            else{
                if(!al.IsMapped())
                    ++stats.NumofUnmappedReads;
                if (!al.IsMateMapped())
                    ++stats.NumofUnmappedReads;
            }
        }
        else{
            ++stats.NumofDuplicates;
        }
    }
//Counts for each base are generated, now coverage histograms will be generated

    vector< ChrStruct >::iterator iter;
/*
    for (iter = CoverageperChr.begin(); iter != CoverageperChr.end(); ++iter) {
        GenerateHistogram_perChr(bamfilename, iter->positionmap, iter->chr_length, iter->chrname);
        cout << "Histogram generated and printed for " << iter->chrname << endl;
    }
*/
    //Generate overal coverage histogram
    GenerateHistogram_Global(bamfilename,"Span", CoverageperChr_span);
    GenerateHistogram_Global(bamfilename,"Base_w_ReadOverlap_Count", CoverageperChr_nooverlapcount);
    GenerateHistogram_Global(bamfilename,"Base_wo_ReadOverlap_Count", CoverageperChr_overlapcount);
    
    cout << "Global Coverage Histograms Generated" << endl;
    
}
void CoverageClass::GenerateHistogram_perChr(string bamfilename, int *coverage, int chr_len, string chrname){
    
    int size_of_bins = 1;
    int number_of_bins = (int) ceil(max_coverage / size_of_bins);
    vector <long unsigned int> histogram(number_of_bins);
    for (int i = 0; i < number_of_bins; ++i)
        histogram[i] = 0;
    
    for(int i = 0; i < chr_len; ++i){
        int bucket = (int)floor(coverage[i] / size_of_bins);
        if(bucket >= histogram.size())
            histogram.resize(bucket+1,0);
        histogram[bucket] += 1;
    }

    string temp = bamfilename.substr(0,bamfilename.size()-3);
    temp.append(".");
    temp.append("CoverageHistogram.");
    temp.append(chrname);
    temp.append(".txt");
    ofstream outf(temp.c_str());

    for (int i = 0; i < histogram.size(); ++i){
        unsigned long int coverage = 0;
        for (unsigned long int j = (histogram.size() - 1); j > i; --j)
            coverage += histogram[j];
        coverage += histogram[i];
      outf << i*size_of_bins << '\t' << double(coverage)/(double(chr_len))*100.0 << endl;
    }
 
}
void CoverageClass::GenerateHistogram_Global(string bamfilename, string coveragetype, vector<ChrStruct> CoverageperChr){
    
    
    int size_of_bins = 1;
    int number_of_bins = (int) ceil(max_coverage / size_of_bins);
    vector <long unsigned int> histogram(number_of_bins);
    for (int i = 0; i < number_of_bins; ++i)
        histogram[i] = 0;
    
    
    vector < double > X;
    vector< double > mean;
    vector< double > median;
    
    vector< ChrStruct >::iterator iter;
    for (iter = CoverageperChr.begin(); iter != CoverageperChr.end(); ++iter) {
        for(int i = 0; i < iter->chr_length; ++i){
            int bucket = (int)floor(iter->positionmap[i] / size_of_bins);
            if(bucket >= histogram.size())
                histogram.resize(bucket+1,0);
            histogram[bucket] += 1;
        }
        cout << iter->chrname << "   finished " << endl;
    }
    
    string temp = bamfilename.substr(0,bamfilename.size()-3);
    temp.append(coveragetype);
    temp.append(".CoverageHistogram.txt");
    ofstream outf(temp.c_str());
    
    double sum = 0;
    for (int i = 0; i < CoverageperChr.size(); ++i)
        sum += CoverageperChr[i].chr_length;
    
    for (int i = 0; i < histogram.size(); ++i){
        unsigned long int coverage = 0;
        for (unsigned long int j = (histogram.size() - 1); j > i; --j)
            coverage += histogram[j];
        coverage += histogram[i];
        outf << i*size_of_bins << '\t' << double (coverage/sum)*100 << endl;
    }
}

#endif
