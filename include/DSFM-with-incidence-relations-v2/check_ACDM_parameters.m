%% preprocessing/preparing data
%unary potential
a = reshape(unary,[427,640]);
%number of columns
col = 640;
% number of rows
row = 427;
% number of superpixels
region = max(max(sup_img));
% number of pixels
N = length(unary);
% number of decomposed parts of the whole submodular function
R = col-1 + row-1 + region;
% incidence_list{i} is the incidence set of the ith submodular function 
incidence_list = cell(R,1);
% parameter_list{i} is the parameter for the ith submodular function 
parameter_list = cell(R,1);
% submodular_type{i} is the type the ith submodular function 
submodular_type = cell(R,1);
% Submodular functions in columns can be decomposed into parallel edges.
% submodular_type is "para_edge"
for i = 1:col-1
    for j = 1:row
        incidence_list{i} = [incidence_list{i} POS(j,i,row) POS(j,i+1,row)];
        parameter_list{i} = [parameter_list{i} lambda * W2(j,i)];
        submodular_type{i} = 'para_edge';
    end
end
% Submodular functions in rows can be decomposed into parallel edges.
% submodular_type is "para_edge"
for i = 1:row-1
    for j = 1:col
        incidence_list{i+col-1} = [incidence_list{i+col-1} POS(i,j,row) POS(i+1,j,row)];
        parameter_list{i+col-1} = [parameter_list{i+col-1} lambda * W1(i,j)];
        submodular_type{i+col-1} = 'para_edge';
    end
end
% Submodular functions in superpixels can be decomposed into F(S) = |S||e/S|.
% submodular_type is "concave_card", parameter_list{i}[j] = F([j]) - F([j-1])
for i = 1:region
    subset = find(sup_img == i);
    card = length(subset);
    subset_row = rem(subset-1, row) + 1;
    subset_col = (subset-1 - rem(subset-1, row))/row + 1;
    incidence_list{i+col-1+row-1} = reshape(POS(subset_row, subset_col, row),1, card);
    parameter_list{i+col-1+row-1} = lambda2*(card + 1 - 2*(1:card));
    submodular_type{i+col-1+row-1} = 'concave_card';
end
bias_vec = reshape(a,1,N);


%% Grid search for parameters for ACDM
gaptol = -1000.00;
numtimes = 3;
alphas = [0.1 0.05 0.02 0.01 0.005];
cs = [10 25 50 100 200]; % c controls the number of iterations for ADCM for each outloop

aa = length(alphas);
cc = length(cs);
Times = zeros(aa,cc);
Gaps = Times;
fprintf('\n%s\n',dataset);

for i = 1:aa
    alpha = alphas(i);
    for b = 1:cc
        c = cs(b);
        K = round(alpha*R);
        T = 50/alpha;
        per = 200;
        record_dis = T/per; % controls how many iterations to perform between recording 
        
        minigaps = zeros(numtimes,1);
        minitimes = zeros(numtimes,1);
        for k = 1:numtimes
            % ACD for DSFM
            tic
            [y, posx, levelcost, timevec, X, sweepi] = ACDM_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, K, T, record_dis, c, gaptol);
            timer = toc;
            gap = levelcost(end)+posx(end);
            minigaps(k) = gap;
            minitimes(k) = timer;
        end
        gapavg = mean(minigaps);
        timeavg = mean(minitimes);
        Gaps(i,b) = gapavg;
        Times(i,b) = timeavg;
        fprintf('alpha = %d, c = %d, time = %f, gap = %f\n',alpha,c,timeavg,gapavg)
    end
    save(strcat('ACDM_tests/',dataset,'_ACDM_search.mat'),'Times','Gaps','alphas','cs')
end
