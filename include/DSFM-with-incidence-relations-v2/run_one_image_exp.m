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

gaptol = -1000.00;  % run until max iterations, don't terminate early

%% Default parameters
alpha = 0.01;
K = round(alpha*R);     % K projections per iteration. Makes less difference for serial code
T = 100/alpha;          % T is the number iterations of RCD, ACD record_dis 
record_dis = T/400;     % controls how many iterations to perform between recording 
% c = 10;                 % c controls the number of iterations for ADCM for each outloop
Tap = 600;              % iterations for IAP
record_dis_ap = Tap/300; % controls iterations/recording
numtimes = 5;

%% run
for i = 1:numtimes

    % Functions return
    %   y: from which we can get primal variables x
    %   posx: used to calculate discrete gap at each step
    %   levelcost: the best primal objective value F(S_lam), as determined by sweep cut
    %   timevec: vector of runtime at each iteration
    tic
    [y, posx, levelcost, timevec, X, sweepi] = RCDM_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, K, T, record_dis,gaptol);
    timer = toc;
    x = -y-bias_vec;
    save(strcat('ContOutput_',dataset,'/RCDM_',num2str(K),'_',num2str(i),'.mat'),'timer','y','levelcost','timevec','posx','x') 

    % ACD for DSFM
    tic
    [y, posx, levelcost, timevec, X, sweepi] = ACDM_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, K, T, record_dis, c, gaptol);
    timer = toc;
    x = -y-bias_vec;
    save(strcat('ContOutput_',dataset,'/ACDM_',num2str(K),'_',num2str(c),'_',num2str(i),'.mat'),'timer','y','levelcost','timevec','posx','x')

    % IAP for DSFM
    % Tap is the number iterations of parallel AP
    tic
    [y, posx, levelcost, timevec, X, sweepi] = APcompact_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, Tap, record_dis_ap, gaptol);
    timer = toc;
    x = -y-bias_vec;
    save(strcat('ContOutput_',dataset,'/APcompact_T_',num2str(Tap),'_',num2str(i),'.mat'),'timer','y','levelcost','timevec','posx','x')

end
