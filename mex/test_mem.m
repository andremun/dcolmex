clear;
close all;

N = 10000;
data = randn(N,10);
data(1:N/2, 10) = 0;
data(N/2+1:N, 10) = 1;

for ii=1:10
    disp(ii);
    [dcol_values, dcol_names] = DCoL(data);
end