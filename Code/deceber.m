%{
Data: Autism TD(Typically Developing) i.e, Healthy samples
Age Group: 4 to 20.
ROI: 264
source: http://umcd.humanconnectomeproject.org
%}

%Reading data
FC = dlmread('UCLA_Autism_TD128_rsfMRI_connectmat.txt');
FC(isinf(FC)) = 0;

W_s = dlmread('UCLA_Autism_TD128B_DTI_connectmat.txt');
D_s = diag(sum(W_s));

%Calculating Laplacian
L_sym = diag(ones(size(W_s,1),1)) - D_s^(-0.5)*W_s*D_s^(-0.5);
Fc_norm = (FC - min(FC(:)))/(max(FC(:))-min(FC(:)));

%Eig decomp. of lap.
[vec, val] = eig(L_sym);

maxi = -1;
t = -1;

%exp(val*i)
for i = 1:50
    Diag_s_exp = diag(exp(-sum(val,1)*i));
    H_s = vec * Diag_s_exp * vec';

    %Normalizing in the range [0 1]
    H_s2 = (H_s - min(H_s(:)))/(max(H_s(:)) - min(H_s(:))); %Normalized Heat Kernal

    %Pearson correlation between H_s2 and Fc_norm
    mean_calc = mean2(H_s2);
    mean_ground = mean2(Fc_norm);

    H_s2 = H_s2 - mean_calc;
    Fc_norm = Fc_norm - mean_ground;
    temp = H_s2.*Fc_norm;
    H_sq = H_s2.*H_s2;
    Fc_sq = Fc_norm.*Fc_norm;

    pear_corr = sum(temp(:))/(sqrt(sum(H_sq(:)))*sqrt(sum(Fc_sq(:))));
    arr(i) = pear_corr;
    fprintf('Per_corr for t=%d', i);
    disp(pear_corr);
    
    %Find Max Pearson correlation coefficient and save t
    if pear_corr>=maxi
        maxi = pear_corr;
        t = i;
    end
    
end

fprintf('Max Pearson correlation = %d', maxi);
fprintf('\nt = %d\n', t);

%Regenerating for displaying o/p
Diag_s_exp = diag(exp(-sum(val,1)*t));
H_s = vec * Diag_s_exp * vec';

%Normalizing in the range [0 1]
H_s2 = (H_s - min(H_s(:)))/(max(H_s(:)) - min(H_s(:))); %Normalized Heat Kernal

figure; title('True Functional Connectivity'); imagesc(Fc_norm);
figure; imagesc(H_s2);
figure; plot(arr);