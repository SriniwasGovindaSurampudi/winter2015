%{
Estimating Functional Connectivity from Structural Connectivity
Dividing randomly into 3 groups of 88 ROI and running for each

Data: Autism TD(Typically Developing i.e, Healthy samples)
Age Group: 4 to 20.
ROI: 264
source: http://umcd.humanconnectomeproject.org
Authors: Abbhinav Venkat, Govinda Sriniwas Surampudi
Creation Date: 17/12/2015
%}

%***************************************************************************************************%
%Reading data

%Adjacency Matrix i.e, Structural Connectivity of size NxN
W = dlmread('UCLA_Autism_TD128B_DTI_connectmat.txt');
D = diag(sum(W, 2));                                    %Diagonal Degree matrix

%Functional Connectivity of size NxN where N = #ROI
FC = dlmread('UCLA_Autism_TD128_rsfMRI_connectmat.txt');    %Ground Truth of FC
FC(isinf(FC)) = 0;                                          %Setting Diag as 0
Fc_n = (FC - min(FC(:)))/(max(FC(:)) - min(FC(:)));      %Normalizing in the range [0 1]

%***************************************************************************************************%
%{
%Dividing randomly into 3 groups of 88 ROI and running for each
r = randperm(size(W, 1));
maxi = 0;

for cntr = 1: 88: 177
    %cntr
    W_s = W(r(cntr:cntr+87), r(cntr:cntr+87));
    Fc_norm = Fc_n(r(cntr:cntr+87), r(cntr:cntr+87));
    D_s = D(r(cntr:cntr+87), r(cntr:cntr+87));
    maxi = maxi + struct_to_func(W_s, Fc_norm, D_s);
    
end

maxi = maxi/3
%}

corr = struct_to_func(W, Fc_n, D);