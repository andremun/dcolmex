type = 0; % 0: normal, 1: legacy to generate identical results compared to C code v1.1; 2: build both

if (type == 0 || type == 2)
    
    %mex -v -largeArrayDims -I/../Source -I/usr/include -I/usr/local/include -cxx -c ...
    %mex -v -largeArrayDims -DUSE_OpenMP GXX= CXXFLAGS="\$CXXFLAGS -fopenmp"  -I/../Source -I/usr/include -I/usr/local/include -cxx -c ...
    mex -v -largeArrayDims -I/../Source -I/usr/include -I/usr/local/include -cxx -c ...
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


    %mex -v -largeArrayDims -o dcol -cxx dcol.o ...
    %mex -v -largeArrayDims -DUSE_OpenMP CXXFLAGS="\$CXXFLAGS -fopenmp" CXXLIBS="\$CXXLIBS -fopenmp" -o dcol -cxx dcol.o ...
    mex -v -largeArrayDims -o dcol -cxx dcol.o ...
                                EuclideanFunction.o ...
                                NormalizedEuclideanFunction.o ...
                                StdWeightedEuclideanFunction.o ...
                                OverlapFunction.o ...
                                VDMFunction.o ...
                                Date.o ...
                                DateContainer.o ...
                                StringTokenizer.o ...
                                Utils.o ...
                                Matrix.o ...
                                Dataset.o ...
                                ExtendedDataset.o ...
                                ComplexityMeasures.o ...
                                Statistics.o ...
                                SMO.o ...
                                InputOptions.o
                            
end

%
if (type == 1 || type == 2)
    
    fprintf('build ori results version\n');
    %mex -v -largeArrayDims -DUSE_OpenMP GXX= CXXFLAGS="\$CXXFLAGS -fopenmp" -I/../Source -I/usr/include -I/usr/local/include -cxx -c ...
    mex -v -largeArrayDims -I/../Source -I/usr/include -I/usr/local/include -cxx -c ...
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


    %mex -v -largeArrayDims -DUSE_OpenMP CXXFLAGS="\$CXXFLAGS -fopenmp" CXXLIBS="\$CXXLIBS -fopenmp" -o dcol_ori -cxx dcol_ori.o ...
    mex -v -largeArrayDims -o dcol_ori -cxx dcol_ori.o ...
                                EuclideanFunction.o ...
                                NormalizedEuclideanFunction.o ...
                                StdWeightedEuclideanFunction.o ...
                                OverlapFunction.o ...
                                VDMFunction.o ...
                                Date.o ...
                                DateContainer.o ...
                                StringTokenizer.o ...
                                Utils.o ...
                                Matrix.o ...
                                Dataset.o ...
                                ExtendedDataset.o ...
                                ComplexityMeasures.o ...
                                Statistics.o ...
                                SMO.o ...
                                InputOptions.o
end
                                                                
