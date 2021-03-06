type = 2; % 0: normal, 1: legacy to generate identical results compared to C code v1.1, 2 :build all

if (type == 0 || type == 2)

    mex -v -largeArrayDims -I../Source -I../Source/DistanceFunctions -c ...
                                ../Source/DistanceFunctions/EuclideanFunction.cpp ...
                                ../Source/DistanceFunctions/NormalizedEuclideanFunction.cpp ...
                                ../Source/DistanceFunctions/StdWeightedEuclideanFunction.cpp ...
                                ../Source/DistanceFunctions/OverlapFunction.cpp ...
                                ../Source/DistanceFunctions/VDMFunction.cpp ...
                                ../Source/Date.cpp ...
                                ../Source/DateContainer.cpp ...
                                ../Source/StringTokenizer.cpp ...
                                ../Source/Utils.cpp ...
                                ../Source/Matrix.cpp ...
                                ../Source/Dataset.cpp ...
                                ../Source/ExtendedDataset.cpp ...
                                ../Source/ComplexityMeasures.cpp ...
                                ../Source/Statistics.cpp ...
                                ../Source/SMO.cpp ...
                                ../Source/InputOptions.cpp ...
                                dcol.cpp


    mex -v -largeArrayDims -output dcol dcol.obj ...
                EuclideanFunction.obj ...
                NormalizedEuclideanFunction.obj ...
                StdWeightedEuclideanFunction.obj ...
                OverlapFunction.obj ...
                VDMFunction.obj ...
                Date.obj ...
                DateContainer.obj ...
                StringTokenizer.obj ...
                Utils.obj ...
                Matrix.obj ...
                Dataset.obj ...
                ExtendedDataset.obj ...
                ComplexityMeasures.obj ...
                Statistics.obj ...
                SMO.obj ...
                InputOptions.obj
                                                                
end

if( type == 1 || type == 2)
    mex -v -largeArrayDims -I../Source -I../Source/DistanceFunctions -c ...
                                ../Source/DistanceFunctions/EuclideanFunction.cpp ...
                                ../Source/DistanceFunctions/NormalizedEuclideanFunction.cpp ...
                                ../Source/DistanceFunctions/StdWeightedEuclideanFunction.cpp ...
                                ../Source/DistanceFunctions/OverlapFunction.cpp ...
                                ../Source/DistanceFunctions/VDMFunction.cpp ...
                                ../Source/Date.cpp ...
                                ../Source/DateContainer.cpp ...
                                ../Source/StringTokenizer.cpp ...
                                ../Source/Utils.cpp ...
                                ../Source/Matrix.cpp ...
                                ../Source/Dataset.cpp ...
                                ../Source/ExtendedDataset.cpp ...
                                ../Source/ComplexityMeasures.cpp ...
                                ../Source/Statistics.cpp ...
                                ../Source/SMO.cpp ...
                                ../Source/InputOptions.cpp ...
                                dcol_ori.cpp


    mex -v -largeArrayDims -output dcol_ori dcol_ori.obj ...
                EuclideanFunction.obj ...
                NormalizedEuclideanFunction.obj ...
                StdWeightedEuclideanFunction.obj ...
                OverlapFunction.obj ...
                VDMFunction.obj ...
                Date.obj ...
                DateContainer.obj ...
                StringTokenizer.obj ...
                Utils.obj ...
                Matrix.obj ...
                Dataset.obj ...
                ExtendedDataset.obj ...
                ComplexityMeasures.obj ...
                Statistics.obj ...
                SMO.obj ...
                InputOptions.obj
end