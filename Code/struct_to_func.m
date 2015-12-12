%{
Estimating Functional Connectivity from Structural Connectivity
Data: Autism TD(Typically Developing i.e, Healthy samples)
Age Group: 4 to 20.
ROI: 264
source: http://umcd.humanconnectomeproject.org
Authors: Abbhinav Venkat, Govinda Surampudi
Date: 11/12/2015
%}

%Reading data

%Adjacency Matrix i.e, Structural Connectivity of size NxN
W_s = dlmread('UCLA_Autism_TD128B_DTI_connectmat.txt');
D_s = diag(sum(W_s));                                       %D matrix

%Functional Connectivity of size NxN where N = #ROI
FC = dlmread('UCLA_Autism_TD128_rsfMRI_connectmat.txt');    %Ground Truth of FC
FC(isinf(FC)) = 0;                                          %Setting Diag as 0
Fc_norm = (FC - min(FC(:)))/(max(FC(:)) - min(FC(:)));      %Normalizing in the range [0 1]

%Calculating Laplacian
L_sym = diag(ones(size(W_s,1),1)) - D_s^(-0.5)*W_s*D_s^(-0.5);

%Eig decomp. of lap.
[vec, val] = eig(L_sym);

%Initializing
maxi = -1;
t = -1;
H_s2 = cell(50, 1);
pear_corr = size(50, 1);

%exp(val*t) 
%Running for different t from 1 to 50.
for i = 1:50
    Diag_s_exp = diag(exp(-sum(val,1)*i));
    H_s = vec * Diag_s_exp * vec';                             %Heat Kernal

    %Normalizing in the range [0 1]
    H_s2{i} = (H_s - min(H_s(:)))/(max(H_s(:)) - min(H_s(:))); %Normalized Heat Kernal

    %Pearson correlation between H_s2 and Fc_norm
    %For drawing graph Corr(FC_ground_truth, FC_calculated) v/s Bt
    mean_calc = mean2(H_s2{i});
    mean_ground = mean2(Fc_norm);

    H_s2{i} = H_s2{i} - mean_calc;
    Fc_norm = Fc_norm - mean_ground;
    temp = H_s2{i}.*Fc_norm;
    H_sq = H_s2{i}.*H_s2{i};
    Fc_sq = Fc_norm.*Fc_norm;

    pear_corr(i) = sum(temp(:))/(sqrt(sum(H_sq(:)))*sqrt(sum(Fc_sq(:))));
    fprintf('Per_corr for t=%d', i);
    disp(pear_corr(i));
    
    %Find Max Pearson correlation coefficient and save best t
    if pear_corr(i)>=maxi
        maxi = pear_corr(i);
        t = i;
    end
    
end

%Results
fprintf('******----------------------------------******\n');
fprintf('Max Pearson correlation = %d', maxi);
fprintf('\nt = %d\n', t);

figure; imagesc(Fc_norm); title('Ground Truth of Functional Connectivity');
figure; imagesc(H_s2{t}); title('Obtained Functional Connectivity');
figure; plot(pear_corr);  title('Corr(FC ground truth, FC calculated) v/s Bt'); 
xlabel('Bt'); ylabel('Correlation of FC calculated with empirical FC');