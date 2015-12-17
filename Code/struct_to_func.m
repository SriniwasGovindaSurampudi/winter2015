%{
Estimating Functional Connectivity from Structural Connectivity
Data: Autism TD(Typically Developing i.e, Healthy samples)
Age Group: 4 to 20.
ROI: 264
source: http://umcd.humanconnectomeproject.org
Authors: Abbhinav Venkat, Govinda Sriniwas Surampudi
Creation Date: 11/12/2015
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

%{
symm_check = L_sym - L_sym';
We notice that symm_check(i, j) = -symm_check(i, j)
The values are of the order of 10^-17 which is very small. 
Therefore, L_sym isn't symmetric. 
This is because of a matlab round off error while calculating inverse
So, we round it off to 10 decimal digits.
%}

%***************************************************************************************************%
%Eig decomp. of lap.
[vec, val] = eig(L_sym);

%Initializing
maxi = -1;
t = -1;
H_s2 = cell(50, 1);
pear_corr = size(50, 1);

%{
exp(val*t) : Heat kernal
Emperically found that Corr saturates as t tends to 50. 
Running for different t from 1 to 50.
%}

for i = 1:100
    Diag_s_exp = diag(exp(-sum(val,1)*i));
    H_s = vec * Diag_s_exp * vec';                          %Heat Kernal

    %Normalizing in the range [0 1]
    H_s2{i} = (H_s - min(H_s(:)))/(max(H_s(:)) - min(H_s(:))); %Normalized Heat Kernal

    %{
    Pearson correlation between H_s2 and Fc_norm
    For drawing graph Corr(FC_ground_truth, FC_calculated) v/s Bt
    2D vector as 1D
    
    mean_calc = mean2(H_s2{i});
    mean_ground = mean2(Fc_norm);

    H = H_s2{i} - mean_calc;
    Fc_norm = Fc_norm - mean_ground;
    temp = H.*Fc_norm;
    H_sq = H.*H;
    Fc_sq = Fc_norm.*Fc_norm;

    pear_corr(i) = sum(temp(:))/(sqrt(sum(H_sq(:)))*sqrt(sum(Fc_sq(:))));
    %}

    %Pearson correlation b/w corresponding rows and then taking mean
    Heat = H_s2{i};
    pear_corr_3 = zeros(size(H_s2, 1), 1);

    for cntr = 1:size(H_s2, 1)
       obs = Heat(cntr, :);
       giv = Fc_norm(cntr, :);

       obs = obs - mean(obs);
       giv = giv - mean(giv);

       temp = obs.*giv;
       obs_sq = obs.*obs;
       giv_sq = giv.*giv;

       pear_corr_3(cntr) = sum(temp(:))/(sqrt(sum(obs_sq(:)))*sqrt(sum(giv_sq(:))));
    end

    pear_corr_calc = abs(pear_corr_3); %removing negative relationships

    pear_corr(i) = mean(pear_corr_calc);   
   
    fprintf('Per_corr for t=%d', i);
    disp(pear_corr(i));
    
    %Find Max Pearson correlation coefficient and save best t
    if pear_corr(i)>=maxi
        maxi = pear_corr(i);
        t = i;
    end
    
end

%***************************************************************************************************%
%Results

fprintf('******----------------------------------******\n');
fprintf('Max Pearson correlation = %d', maxi);
fprintf('\nt = %d\n', t);

figure; 
subplot(1, 2, 1); imagesc(Fc_norm); title('Ground Truth of Functional Connectivity');
subplot(1, 2, 2); imagesc(H_s2{t}); title('Obtained Functional Connectivity');
figure; plot(pear_corr);  title('Corr(FC ground truth, FC calculated) v/s Bt'); 
xlabel('Bt'); ylabel('Correlation of FC calculated with empirical FC');

%***************************************************************************************************%