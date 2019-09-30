#include <matrix.h>
#include <mex.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <vector>
#include <set>
#include <algorithm>
#include <map>

using namespace std;

#include "../Source/ComplexityMeasures.h"
#include "../Source/InputOptions.h"

#define NUM_FEATURES 14

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //[MFs] = metafeature(double** data, double* col_type)

  mexPrintf("dcol_ori ...\n");
#ifdef PRINT_DEBUG
  mexPrintf("[INFO] Using normalized Eulicdean distance functions for attributes\n");
#endif
  int nrows, ncols;
  double* matlab_data;
  double* matlab_col_type = NULL;
  double* matlab_mask = NULL;

  bool has_na = false;
  int  na_option;
  bool feature_mask[NUM_FEATURES];
  int* col_type = NULL;
  int num_output = 0;
  float*  n_data;    // normalized data, row wise
 
  if(nrhs < 2 || nrhs > 4) {
    mexPrintf("Error: Incorrect input arguments\n");
    mexPrintf("Usage: dcol(data, na_option) or dcol(data, na_option, col_type) or dcol(data, na_option, col_type, mask)\n");
    return;
  }
  const mwSize *dims = mxGetDimensions(prhs[0]);
  nrows = dims[0]; ncols = dims[1];
#ifdef PRINT_DEBUG
  mexPrintf("[INFO] nrows: %d, ncols: %d\n", nrows, ncols);
#endif

  matlab_data = mxGetPr(prhs[0]);
  na_option = (int)mxGetScalar(prhs[1]);
  if(na_option < 0 || na_option > 2) {
    mexPrintf("Error: na_option must be in {0: NA_KEEP, 1: NA_REMOVE, 2: NA_ESTIMATE}\n");
      return;
  }

  if(nrhs == 3 || nrhs == 4) {
    matlab_col_type = mxGetPr(prhs[2]);
    if(mxGetN(prhs[2]) != ncols) {
      mexPrintf("Error: col_type must have %d elements\n", ncols);
      return;
    }
    for(int i=0; i < ncols-1; i++) {
      if(matlab_col_type[i] != 0 && matlab_col_type[i] != 1) {
        mexPrintf("Error: data att type must be 0 (numeric) or 1 (nominal)\n");
        return;
      }
    }
    if(matlab_col_type[ncols-1] != 2) {
      mexPrintf("Error: att type of the last column must be 2 (CLASS)\n");
      return;
    }
  }

  if(nrhs == 4) {
    matlab_mask = mxGetPr(prhs[3]);
    if(mxGetN(prhs[3]) != NUM_FEATURES) {
      mexPrintf("Error: mask must have %d elements\n", NUM_FEATURES);
      return;
    }
  }

  col_type = new int[ncols];
  for(int i=0; i < ncols; i++) {
    if(matlab_col_type)
      col_type[i] = (int)matlab_col_type[i];
    else
      col_type[i] = NUMERICAL;
  }

  for(int i=0; i < NUM_FEATURES; i++) {
    if(matlab_mask)
      feature_mask[i] = matlab_mask[i];
    else
      feature_mask[i] = true;
    if(feature_mask[i])
      num_output++;
  }

#ifdef PRINT_DEBUG
  mexPrintf("[INFO] Number of outputs: %d\n", num_output);
#endif

  // pre-process data (using matlab_data[row+col*nrows])
  float *data2;
  if(na_option == NA_REMOVE) {
    bool* na_row_mask = new bool[nrows]();
    int na_row_count = 0;
    for(int row=0; row < nrows; row++) {
      for(int col=0; col < ncols; col++) {
        if(isnan(matlab_data[row + col*nrows])) {
          na_row_mask[row] = true;
          has_na = true;
          break;
        }
      }
      if(!na_row_mask[row])
        na_row_count++;
    }
  
    data2 = new float[na_row_count*ncols];

    int r_ind = 0;
    for(int row=0; row < nrows; row++) {
        if(!na_row_mask[row]) {
          for(int col=0; col < ncols; col++)
            data2[r_ind*ncols + col] = matlab_data[row + col*nrows];
          r_ind++;
        }
    }
    nrows = na_row_count;
    delete [] na_row_mask;

  } else {
    data2 = new float[nrows*ncols];

    for(int row=0; row < nrows; row++) {
      for(int col=0; col < ncols; col++) {
        if(isnan(matlab_data[row + col*nrows])) {
          has_na = true;
          data2[row*ncols + col] = Dataset::UNKNOWN_VALUE;
        }
        else {
          data2[row*ncols + col] = matlab_data[row + col*nrows];
        }
      }
    }
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
  int num_out_cols = 1;
  if(num_classes > 2)
    num_out_cols = num_classes;
#ifdef PRINT_DEBUG
  mexPrintf("[INFO] Number of classes: %d\n", num_classes);
  mexPrintf("[INFO] Number of loops: %d\n", num_out_cols);
#endif

  plhs[0] = mxCreateDoubleMatrix(num_output, num_out_cols, mxREAL); 
  double *dcol_values = mxGetPr(plhs[0]);
  for(int i=0; i < num_output*num_out_cols; i++) {
    dcol_values[i] = -1;
  }

  // find att min max
  int att_cols = ncols-1;
  double* att_min = new double[att_cols];
  double* att_max = new double[att_cols];
  double* att_diff = new double[att_cols];
  for(int col = 0; col < att_cols; col++) {
    
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
    if((int)col_type[col] == NUMERICAL) {
      it = std::min_element(std::begin(col_values), std::end(col_values));
      att_min[col] = *it;
    }
    else
      att_min[col] = 0;

    it = std::max_element(std::begin(col_values), std::end(col_values));
    if((int)col_type[col] == NUMERICAL)
      att_max[col] = *it;
    else
      att_max[col] = *it - 1;
    
    att_diff[col] = att_max[col] - att_min[col];
    
#ifdef PRINT_DEBUG
    mexPrintf("[INFO] min: %0.2f max: %0.2f\n", att_min[col], att_max[col]);
#endif
  }

  // normalize data
  n_data = new float[nrows * att_cols]; // row wise matrix !!!
  if(!n_data) {
   mexPrintf("[ERROR] Not enough memmory!\n");
   return;
  }
  for(int row = 0; row < nrows; row++) {
    for(int col = 0; col < att_cols; col++) {
        double v = data2[row*ncols + col];
        if(v == Dataset::UNKNOWN_VALUE)
          n_data[row*att_cols+col] = Dataset::UNKNOWN_VALUE;
        else {
          if(col_type[col] == NOMINAL)
            v--;
          if(att_diff[col] == 0)
            n_data[row*att_cols+col] = 1;
          else
            n_data[row*att_cols+col] = (v - att_min[col]) /  att_diff[col];
        }
    }
  }
  delete []data2;

  // DCoL measures
  string l_str[] = {"F1", "F1v", "F2", "F3", "F4",
                    "L1", "L2", "L3", 
                    "N1", "N2", "N3", "N4",
                    "T1", "T2"}; //14
  plhs[1] = mxCreateCellMatrix(num_output, 1);
  int ind = 0;
  for(int i=0; i < NUM_FEATURES; i++) {
    if(feature_mask[i]) {
      mxSetCell(plhs[1], ind, mxCreateString(l_str[i].c_str()));
      ind++;
    }
  }

  if (num_classes < 2 || nrows < 1 || ncols < 2) {
    mexPrintf("Error: number of classes < 2 or nrows < 1 or ncols < 2\n", NUM_FEATURES);
    return;
  }

  // DCoL
  InputOptions opts;
  int argc = 9;
  char** argv = new char*[argc];
  for(int i=0; i < argc; i++) {
    argv[i] = new char[255];
  }
  strcpy(argv[0], "dcol"); strcpy(argv[1], "-i"); strcpy(argv[2], "MATLAB");
  strcpy(argv[3], "-o"); strcpy(argv[4], "NONE"); 
  strcpy(argv[5], "-A"); strcpy(argv[6], "-d");
  strcpy(argv[7], "-s"); strcpy(argv[8], "10");
  //strcpy(argv[9], "-cM"); strcpy(argv[10], "1");
  //strcpy(argv[11], "-nM"); strcpy(argv[12], "3");
  opts.parseInput ( argc, argv );

  if ( !opts.isAnyOptionSelected () ) {
    mexPrintf("Error: Call syntax is incomplete. \n");
    return;
  }

  if ( opts.isIncompatibleOptions () ) {
    mexPrintf("Error: Incompatible options. \n");
    return;
  }

  if ( !opts.isAnyComplexityMeasureSelected () )  {
    mexPrintf("Error: No measure is selected. \n");
    return;
  }

  bool rmNA = (has_na && na_option == NA_ESTIMATE) ? true : false;

  map<string, int> mf_map;
  mf_map["F1"]=0; mf_map["F1v"]=1; mf_map["F2"]=2; mf_map["F3"]=3; mf_map["F4"]=4;
  mf_map["L1"]=5; mf_map["L2"]=6; mf_map["L3"]=7; 
  mf_map["N1"]=8; mf_map["N2"]=9; mf_map["N3"]=10; mf_map["N4"]=11;
  mf_map["T1"]=12; mf_map["T2"]=13;

  unsigned int start_t = Utils::getTime();

#ifdef USE_OpenMP
#pragma omp parallel for
#endif
  for (int run = 0; run < num_out_cols; run++ ) {

    int* class_of_sample = new int[nrows];
    if(num_out_cols > 1) {
#ifdef PRINT_DEBUG
      mexPrintf("[INFO] Processing dataset: %d / %d\n", run, num_out_cols);
#endif
      int c = *std::next(class_labels.begin(), run);
      for(int row = 0; row < nrows; row++)
          class_of_sample[row] = (class_values[row] == c) ? 0 : 1;
    }
    else {
      for(int row =0; row < nrows; row++)
        class_of_sample[row] = class_values[row];
    }

    double mf_values[NUM_FEATURES];

    ComplexityMeasures* dSet;
    //mexPrintf("[INFO] Run complexity measures... \n");
    //mexPrintf("[INFO] con dist func: %d\n", opts.getTypeOfContinuousDistFunction ());
    //mexPrintf("[INFO] nom dist func: %d\n", opts.getTypeOfNominalDistFunction () );
    dSet = new ComplexityMeasures ( n_data, class_of_sample, nrows, ncols, col_type, att_min, att_max, NULL,
                                     opts.getTypeOfContinuousDistFunction (), opts.getTypeOfNominalDistFunction (), rmNA, true);

    //F1
    //mexPrintf("[INFO] F1... \n");
    int att;
    if(feature_mask[mf_map["F1"]] || feature_mask[mf_map["F3"]] || feature_mask[mf_map["F4"]]) {
      float v = dSet->computeFisher ( att );
      if(feature_mask[mf_map["F1"]]) {
        mf_values[mf_map["F1"]] = v;
      }
    } 

    //F1v
    //mexPrintf("[INFO] F1v... \n");
    if(feature_mask[mf_map["F1v"]]) {
      mf_values[mf_map["F1v"]] = dSet->computeFisherVectorized ();
    }

    //F2
    //mexPrintf("[INFO] F2... \n");
    if(feature_mask[mf_map["F2"]]) {
      mf_values[mf_map["F2"]]  = dSet->computeVolumeOverlap();
    }

    //F3 and F4
    //mexPrintf("[INFO] F3 and F4... \n");
    float* vectorResults;
    if(feature_mask[mf_map["F3"]] || feature_mask[mf_map["F4"]]) {
      float  measureResult;
      float  measureResultAux;
      vectorResults = dSet->computeMaximumEfficiencyOfAttributes ( att, measureResult );
      measureResultAux = 0;
      for ( int i = 0; i < dSet->getNumberOfAttributes (); i++ ) {
        measureResultAux += vectorResults[i] / dSet->getNumberOfExamples ();
      }
      delete [] vectorResults;
      if(feature_mask[mf_map["F3"]]) {
        mf_values[mf_map["F3"]] = measureResult;
      }
      if(feature_mask[mf_map["F4"]]) { 
        mf_values[mf_map["F4"]] = measureResultAux;
      } 
    }
    
    //N1
    //mexPrintf("[INFO] N1... \n");
    if(feature_mask[mf_map["N1"]]) {
      int** spanningTree = dSet->computePrim ();
      mf_values[mf_map["N1"]]  = dSet->computeBoundary ( spanningTree );
      for ( int i = 0; i < dSet->getNumberOfExamples () - 1; i++ ) {
          delete [] spanningTree[i];
      }
      delete [] spanningTree;
    }

    //N2
    //mexPrintf("[INFO] N2... \n");
    if(feature_mask[mf_map["N2"]]) {
      mf_values[mf_map["N2"]] = dSet->computeIntraInter ();
    }

    //N3
    //mexPrintf("[INFO] N3... \n");
    if(feature_mask[mf_map["N3"]]) {
      mf_values[mf_map["N3"]] = dSet->computeNonLinearityKNNTrain ();
    }
   
    //N4
    //mexPrintf("[INFO] N4... \n");
    if(feature_mask[mf_map["N4"]]) {
      mf_values[mf_map["N4"]] = dSet->computeNonLinearityKNNConvexHull ();
    }
    
    //T1
    //mexPrintf("[INFO] T1... \n");
    if(feature_mask[mf_map["T1"]]) {
      vectorResults = dSet->computeFractMaxCoveringSpheres();
      mf_values[mf_map["T1"]] = vectorResults[0] / dSet->getNumberOfExamples ();
      delete [] vectorResults;
    }
    
    //T2
    //mexPrintf("[INFO] T2... \n");
    if(feature_mask[mf_map["T2"]]) {
      mf_values[mf_map["T2"]]  = dSet->averageNumberOfSamplesPerDimension ();
    }

    //L1 L2 L3
    if(feature_mask[mf_map["L1"]] || feature_mask[mf_map["L2"]] || feature_mask[mf_map["L3"]]) {
      float* w = 0; float b = 0;
      w = dSet->trainSMO ( b );
      //L1 
      //mexPrintf("[INFO] L1... \n");
      if(feature_mask[mf_map["L1"]]) {
        mf_values[mf_map["L1"]] = dSet->computeNonLinearityLCDistance ( w, b );
      }

      //L2
      //mexPrintf("[INFO] L2... \n");
      if(feature_mask[mf_map["L2"]]) {
        mf_values[mf_map["L2"]] = dSet->computeNonLinearityLCTrain ( w, b );
      }

      //L3
      //mexPrintf("[INFO] L3... \n");
      if(feature_mask[mf_map["L3"]]) {
        mf_values[mf_map["L3"]]  = dSet->computeNonLinearityLCConvexHull ( w, b );
      }

      delete [] w;
    }

    int ind = 0;
    for(int i=0; i < NUM_FEATURES; i++) {
      if(feature_mask[i]) {
        dcol_values[run * num_output + ind] = mf_values[i];
        ind++;
      }
    }

    delete dSet;
    delete [] class_of_sample;
  }

  delete [] col_type;
  delete [] n_data;
  delete [] att_min;
  delete [] att_max;
  delete [] att_diff;
  for(int i=0; i < argc; i++)
    delete [] argv[i];
  delete [] argv;
 
  mexPrintf("dcol_ori - processing time: %d\n", Utils::getTime() - start_t);
}