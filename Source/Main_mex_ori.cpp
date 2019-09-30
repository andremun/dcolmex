#include "ComplexityMeasures.h"
#include "InputOptions.h"
#include "ResultsContainer.h"

#include <vector>
#include <set>
#include <algorithm>
#include <omp.h>
#include <map>

using namespace std;

#define NUM_RESULTS 14

int main ( int argc, char** argv ) {

	if(argc != 3) {
		cout << "DCoL arff_in_filename na_option";
		return 1;
	}

	float* data;
	int nrows, ncols;
	vector<Attribute> att_type;
	bool has_na;
	data = Utils::loadDataFromArffFile(argv[1], nrows, ncols, att_type, has_na);
	if(!data) {
		cout << "Fail to load data" << endl;
		return 1;
	}

	cout << "nrows: " << nrows << " ncols: " << ncols << endl;
	for(int i=0; i < att_type.size(); i++)
		cout << att_type[i].type << " ";
	cout << endl << endl;

	int na_option = atoi(argv[2]);
	if(na_option < 0 || na_option > 2) {
		cout << "wrong na_option value" << endl;
		return 1;
	}

	// pre-process data
	float *data2;
	if(has_na && na_option == NA_REMOVE) {
		bool* na_row_mask = new bool[nrows]();
		int na_row_count = 0;
		for(int r=0; r < nrows; r++) {
			for(int c=0; c < ncols; c++) {
				if(data[r*ncols + c] == Dataset::UNKNOWN_VALUE) {
					na_row_mask[r] = true;
					break;
				}
			}
			if(!na_row_mask[r])
				na_row_count++;
		}
		data2 = new float[na_row_count*ncols];
		Utils::removeNA(data, nrows, ncols, na_row_mask, data2);
		nrows = na_row_count;

		delete [] na_row_mask;
		delete [] data;

	}else {
		data2 = data;
	}

	//=========================================================
	vector<int> class_values;
	set<int> class_labels;
	for(int row = 0; row < nrows; row++) {
	  	int c = data2[row*ncols + ncols-1];
	  	class_values.push_back(c);
	  	class_labels.insert(c);
	}
	int num_classes = class_labels.size();
	cout << "[INFO] Number of classes: " << num_classes << endl;
	
	int num_out_cols = 1;
	if(num_classes > 2)
	  	num_out_cols = num_classes;
	cout << "[INFO] Number of loops: " << num_out_cols << endl;

    double dcol_values[NUM_RESULTS*num_out_cols];
    for(int i=0; i < NUM_RESULTS*num_out_cols; i++)
		dcol_values[i] = -1;

	// find att min max
	int att_cols = ncols-1;
	double* att_min = new double[att_cols];
	double* att_max = new double[att_cols];
	double* att_diff = new double[att_cols];
	for(int col = 0; col < att_cols; col++) {
		if(att_type[col].type == NUMERICAL) {
			vector<double> col_values;
			for(int row = 0; row < nrows; row++)
				if(data2[row*ncols + col] != Dataset::UNKNOWN_VALUE)
					col_values.push_back(data2[row*ncols + col]);
			if(col_values.size() == 0) {
				att_min[col] = 0;
				att_max[col] = 0;
				att_diff[col] = 0;
				continue;
			}
			vector<double>::iterator it;
			it = std::min_element(std::begin(col_values), std::end(col_values));
			att_min[col] = *it;
			it = std::max_element(std::begin(col_values), std::end(col_values));
			att_max[col] = *it;
			att_diff[col] = att_max[col] - att_min[col];\

		} else {
			att_min[col] = 0;
			att_max[col] = att_type[col].str_int.size() - 1;
			att_diff[col] = att_max[col] - att_min[col];
		}

		
		cout << "[INFO] min: " << att_min[col] << " max: " << att_max[col] << endl;
	}

	// normalize data
	cout << "[INFO] normalize data" << endl;
	float* n_data = new float[nrows * att_cols]; // row wise matrix !!!
	for(int row = 0; row < nrows; row++) {
		for(int col = 0; col < att_cols; col++) {
		  	double v = data2[row*ncols + col];
		  	if(v == Dataset::UNKNOWN_VALUE)
		  		n_data[row*att_cols+col] = Dataset::UNKNOWN_VALUE;
		  	else {
		  		if(att_diff[col] == 0)
			  		n_data[row*att_cols+col] = 1;
			  	else
			  		n_data[row*att_cols+col] = (v - att_min[col]) /  att_diff[col];
		  	}
		}
	}

	delete []data2;

	// calculate distance matrix to speed up
	float* dmat = NULL;
	if (!has_na || na_option != NA_ESTIMATE) {
		cout << "[INFO] calculate dmat" << endl;
		dmat = new float[nrows * nrows];
		Utils::distanceCPU(n_data, nrows, att_cols, dmat, Dataset::UNKNOWN_VALUE);
	}

	// DCoL
	InputOptions opts;
	int myargc = 10;
	char** myargv = new char*[myargc];
		for(int i=0; i < myargc; i++) {
		myargv[i] = new char[255];
	}
	strcpy(myargv[0], "dcol"); strcpy(myargv[1], "-i"); strcpy(myargv[2], "MATLAB");
	strcpy(myargv[3], "-o"); strcpy(myargv[4], "NONE"); 
	strcpy(myargv[5], "-A"); strcpy(myargv[6], "-d");
	strcpy(myargv[7], "-s"); strcpy(myargv[8], "10");
	strcpy(myargv[9],"-D");
	//strcpy(myargv[10], "-cM"); strcpy(myargv[11], "1");
	//strcpy(myargv[12], "-nM"); strcpy(myargv[13], "3");
	
	opts.parseInput ( myargc, myargv );

	if ( !opts.isAnyOptionSelected () ) {
		printf("[ERROR] Call syntax is incomplete. \n");
		return 1;
	}

	if ( opts.isIncompatibleOptions () ) {
		printf("[ERROR] Incompatible options. \n");
		return 1;
	}

	if ( !opts.isAnyComplexityMeasureSelected () )  {
    	printf("[ERROR] No measure is selected. \n");
    	return 1;
  	}

  	cout << "[INFO]  start calculating..." << endl;

  	bool rmNA = (has_na && na_option == NA_ESTIMATE) ? true : false;

  	map<string, int> mf_map;
  	mf_map["F1"]=0; mf_map["F1v"]=1; mf_map["F2"]=2; mf_map["F3"]=3; mf_map["F4"]=4;
	mf_map["L1"]=5; mf_map["L2"]=6; mf_map["L3"]=7; 
	mf_map["N1"]=8; mf_map["N2"]=9; mf_map["N3"]=10; mf_map["N4"]=11;
	mf_map["T1"]=12; mf_map["T2"]=13;

	int* col_type = new int[ncols];
	for(int col = 0; col < ncols; col++)
		col_type[col] = att_type[col].type;
  	
#ifdef USE_OpenMP
#pragma omp parallel for
#endif
	for (int run = 0; run < num_out_cols; run++ ) {

		int* class_of_sample = new int[nrows];

	    if(num_out_cols > 1) {
	      int c = *std::next(class_labels.begin(), run);
	      for(int row = 0; row < nrows; row++)
            class_of_sample[row] = (class_values[row] == c) ? 0 : 1;
	    }
	    else {
	    	for(int row =0; row < nrows; row++)
	    		class_of_sample[row] = class_values[row];
	    }

	    int res_ind = 0;

	    ComplexityMeasures* dSet;
        printf("[INFO] Run complexity measures... \n");
	    dSet = new ComplexityMeasures ( n_data, class_of_sample, nrows, ncols, col_type, att_min, att_max, NULL,
                                     opts.getTypeOfContinuousDistFunction (), opts.getTypeOfNominalDistFunction (), rmNA, true);

	    //F1
	    int att;
	    res_ind = mf_map["F1"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeFisher ( att );
	    cout << "F1: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    //F1v
	    res_ind = mf_map["F1v"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeFisherVectorized ();
	    cout << "F1v: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    //F2
	    res_ind = mf_map["F2"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeVolumeOverlap();
	    cout << "F2: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    //F3 and F4
	    float* vectorResults;
	    float  measureResult;
	    float  measureResultAux;
	    vectorResults = dSet->computeMaximumEfficiencyOfAttributes ( att, measureResult );
	    measureResultAux = 0;
	    for ( int i = 0; i < dSet->getNumberOfAttributes (); i++ ) {
	      measureResultAux += vectorResults[i] / dSet->getNumberOfExamples ();
	    }
	    delete [] vectorResults;
	    res_ind = mf_map["F3"];
	    dcol_values[run * NUM_RESULTS + res_ind] = measureResult;
	    cout << "F3: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;
	  
	  	res_ind = mf_map["F4"];
	    dcol_values[run * NUM_RESULTS + res_ind] = measureResultAux;
	    cout << "F4: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;
	    
	    //N1
	    int** spanningTree = dSet->computePrim ();
	    res_ind = mf_map["N1"];
        dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeBoundary ( spanningTree );
        cout << "N1: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;
	    for ( int i = 0; i < dSet->getNumberOfExamples () - 1; i++ ) {
	        delete [] spanningTree[i];
	    }
	    delete [] spanningTree;

	    //N2
	    res_ind = mf_map["N2"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeIntraInter ();
	    cout << "N2: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    //N3
	    res_ind = mf_map["N3"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeNonLinearityKNNTrain ();
	    cout << "N3: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;
	   
	    //N4
	    res_ind = mf_map["N4"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeNonLinearityKNNConvexHull ();
	    cout << "N4: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;
	    
	    //T1
	    vectorResults = dSet->computeFractMaxCoveringSpheres();
	    res_ind = mf_map["T1"];
	    dcol_values[run * NUM_RESULTS + res_ind] = vectorResults[0] / dSet->getNumberOfExamples ();
	    cout << "T1: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    delete [] vectorResults;
	    
	    //T2
	    res_ind = mf_map["T2"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->averageNumberOfSamplesPerDimension ();
	    cout << "T2: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    //L1 L2 L3
	    float* w = 0; float b = 0;
	    w = dSet->trainSMO ( b );
	    //L1 
	    res_ind = mf_map["L1"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeNonLinearityLCDistance ( w, b );
	    cout << "L1: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    //L2
	    res_ind = mf_map["L2"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeNonLinearityLCTrain ( w, b );
	    cout << "L2: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    //L3
	    res_ind = mf_map["L3"];
	    dcol_values[run * NUM_RESULTS + res_ind] = dSet->computeNonLinearityLCConvexHull ( w, b );
	    cout << "L3: " << dcol_values[run * NUM_RESULTS + res_ind] << endl;

	    delete [] w;

	    delete dSet;
	    delete [] class_of_sample;
	}

	delete [] n_data;
	delete [] dmat;
	delete [] att_min;
	delete [] att_max;
	delete [] att_diff;

	return 0;
}