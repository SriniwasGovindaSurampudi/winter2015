%{
Model Inversion: 
Estimating Structural Connectivity given functional connectivity and degree
of structural graph

Data: Autism TD(Typically Developing i.e, Healthy samples)
Age Group: 4 to 20.
ROI: 264
source: http://umcd.humanconnectomeproject.org
Authors: Abbhinav Venkat, Govinda Sriniwas Surampudi
Creation Date: 12/12/2015
%}

%*************************************************************************%
%SVD

[U, S, V] = svd(H_s2{t});
sigma = sum(S, 2);
[m, n] = find(sigma<sigma(1)*0.1);          %Choosing 90% of max value
sigma(m) = 0;                               %Eliminating values < 0.1*S(1,1)
S = diag(sigma);

%*************************************************************************%
%Re-computing the Functional Connectivity

Cf = U*S*V';
Cf = round(Cf, 7);
symm_check = Cf - Cf';
assert(isempty(find(symm_check~=0, 1)));    %To check for symmetry of Cf

%Cf = H_s2{t};

%We need to ensure that Cf is symmetric to get real eigen values
[nvec, nval] = eig(Cf);
gamma = -log(sum(nval,1))/t;                %!!!!!!!!COMPLEX

%Re-computing L
L = nvec*diag(gamma)*nvec';

%Finding Structural Connectivity
C = (D_s^0.5)*(diag(ones(size(W_s,1),1)) - L)*(D_s^0.5);

%*************************************************************************%
%Correlation between emerical and estimated structure 
%{
mean_calc = mean2(H_s2{i});
mean_ground = mean2(Fc_norm);

H = H_s2{i} - mean_calc;
Fc_norm = Fc_norm - mean_ground;
temp = H.*Fc_norm;
H_sq = H.*H;
Fc_sq = Fc_norm.*Fc_norm;

pear_corr(i) = sum(temp(:))/(sqrt(sum(H_sq(:)))*sqrt(sum(Fc_sq(:))));
fprintf('Per_corr for t=%d', i);
disp(pear_corr(i));
%}