clear;
close all;

test_case = 1;  % 1: current_1 (958 x 10, 2 classes) 
                % 2: IRIS (150 x 5, 3 classes) 
                % 3: abalone: 4177 x 9, 3 classes, first col is nominal

test_data = load('test_data.mat');

% the last colume of data is CLASS
data = test_data.test{test_case}.data;
col_type = test_data.test{test_case}.col_type;
na_option = 1;

% DCoL features (14): "F1", "F1v", "F2", "F3", "F4", "L1", "L2", "L3", 
%                     "N1", "N2", "N3", "N4","T1", "T2"
dcol_mask = [0,0,0,1,1,0,1,0,1,0,0,1,0,0];

tic
% outputs: row wise matrixes
[dcol_values01, dcol_names01] = dcol(data, na_option);
[dcol_values02, dcol_names02] = dcol(data, na_option, col_type);
[dcol_values03, dcol_names03] = dcol(data, na_option, col_type, dcol_mask);

%[dcol_values11, dcol_names11] = dcol_ori(data, na_option);
%[dcol_values12, dcol_names12] = dcol_ori(data, na_option, col_type);
%[dcol_values13, dcol_names13] = dcol_ori(data, na_option, col_type, dcol_mask);

toc
