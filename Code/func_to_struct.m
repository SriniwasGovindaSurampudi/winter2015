%{
Estimating Structural Connectivity given functional connectivity and degree
of structural graph
Data: Autism TD(Typically Developing i.e, Healthy samples)
Age Group: 4 to 20.
ROI: 264
source: http://umcd.humanconnectomeproject.org
Authors: Abbhinav Venkat, Govinda Sriniwas Surampudi
Date: 12/12/2015
%}

%SVD
[U, S, V] = svd(H_s2{t});
sigma = sum(S, 2);
[m, n] = find(sigma<sigma(1)*0.1);
sigma(m) = 0;                       %Eliminating values < 0.1*S(1,1)
S = diag(sigma);

%Re-computing the Functional Connectivity
Cf = U*S*V';

[nvec, nval] = eig(Cf);
a = sum(nval,2);

gamma = -log(sum(nval,1))/t;

%Re-computing L
L = nvec*diag(gamma)*nvec';

%Finding Structural Connectivity
C = (D_s^0.5)*(diag(ones(size(W_s,1),1)) - L)*(D_s^0.5);
