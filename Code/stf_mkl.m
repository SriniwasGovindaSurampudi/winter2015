%{
Estimating Functional Connectivity from Structural Connectivity
Using Multiple Kernal Learning

Data: Autism TD(Typically Developing i.e, Healthy samples)
Age Group: 4 to 20.
ROI: 264
source: http://umcd.humanconnectomeproject.org
Authors: Abbhinav Venkat, Govinda Sriniwas Surampudi
Creation Date: 16/12/2015
%}

%***************************************************************************************************%
%Reading data

%Adjacency Matrix i.e, Structural Connectivity of size NxN
W_s = dlmread('UCLA_Autism_TD128B_DTI_connectmat.txt');
D_s = diag(sum(W_s, 2));                                    %Diagonal Degree matrix

%Functional Connectivity of size NxN where N = #ROI
FC = dlmread('UCLA_Autism_TD128_rsfMRI_connectmat.txt');    %Ground Truth of FC
FC(isinf(FC)) = 0;                                          %Setting Diag as 0
Fc_norm = (FC - min(FC(:)))/(max(FC(:)) - min(FC(:)));      %Normalizing in the range [0 1]

%***************************************************************************************************%
%Calculating normalized weight matrix

norm_W = ((D_s^(-0.5)))*W_s*((D_s^(-0.5)));
norm_W = round(norm_W, 10);
symm_check = norm_W - norm_W';
assert(isempty(find(symm_check~=0, 1)));                    %To check for symmetry of norm_W

%Calculating normalized symmetrix Laplacian
L_sym = eye(size(W_s,1)) - norm_W;
symm_check = L_sym - L_sym';
assert(isempty(find(symm_check~=0, 1)));                    %To check for symmetry of L_sym

%***************************************************************************************************%
%Eig decomp. of lap.
[vec, val] = eig(L_sym);

maxi = zeros(size(val, 1), 1);
maxi(:, :) = -1;

scale = zeros(size(val, 1), 1);
Cf = cell(size(val, 1), 1);
Cf_final = zeros(size(val));

%For each ROI, t=1..50
for i = 1:size(val, 1)
    for j = 1:50
       Diag_s_exp = diag(exp(-sum(val, 1)*j));
       H_s = vec * Diag_s_exp * vec';                          %Heat Kernal
        
       %Normalizing in the range [0 1]
       H_s2{j} = (H_s - min(H_s(:)))/(max(H_s(:)) - min(H_s(:))); %Normalized Heat Kernal

       %Assuming each ROI to be active independently
       ei = zeros(size(val, 1), 1);
       ei(i) = 1;
       
       H_s2{j} = H_s2{j}*ei;        %nx1
       
       %Correlation 
       mean_calc = mean2(H_s2{j});
       mean_ground = mean2(Fc_norm(:, i));

       H = H_s2{j} - mean_calc;
       Fc = Fc_norm(:, i) - mean_ground;
       temp = H.*Fc;
       H_sq = H.*H;
       Fc_sq = Fc.*Fc;

       pear_corr = sum(temp(:))/(sqrt(sum(H_sq(:)))*sqrt(sum(Fc_sq(:))));
       
       if pear_corr > maxi(i)
           maxi(i) = pear_corr;
           scale(i) = j;
       end
    end
    Cf{i} = H_s2{scale(i)};
    Cf_final(:, i) = Cf{i}; %Combining with B = 1 (Linear coeff)
end

plot(scale);

%***************************************************************************************************%
%Learning B



%***************************************************************************************************%
%{
Pearson Correlation
Treating 2D matrix as a 1D vector

Fc_norm = round(Fc_norm, 10);
mean_calc = mean2(Cf_final);
mean_ground = mean2(Fc_norm);

C1 = Cf_final - mean_calc;
W1 = Fc_norm - mean_ground;

temp = C1.*W1;
C1_sq = C1.*C1;
W1_sq = W1.*W1;

pear_corr_3 = sum(temp(:))/(sqrt(sum(C1_sq(:)))*sqrt(sum(W1_sq(:))));
%}
Fc_norm = round(Fc_norm, 10);

%Pearson correlation b/w corresponding rows and then taking mean
pear_corr_3 = zeros(size(Cf_final, 1), 1);

for cntr = 1:size(Cf_final, 1)
   obs = Cf_final(cntr, :);
   giv = Fc_norm(cntr, :);
   
   obs = obs - mean(obs);
   giv = giv - mean(giv);
   
   temp = obs.*giv;
   obs_sq = obs.*obs;
   giv_sq = giv.*giv;
   
   pear_corr_3(cntr) = sum(temp(:))/(sqrt(sum(obs_sq(:)))*sqrt(sum(giv_sq(:))));
end

pear_corr_calc = abs(pear_corr_3); %removing negative relationships
final_corr = mean(pear_corr_calc);

%***************************************************************************************************%
