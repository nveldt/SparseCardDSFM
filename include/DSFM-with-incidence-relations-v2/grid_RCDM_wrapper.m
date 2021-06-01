%% plant 500
clear;
addpath('../../data/data_reflections/')
dataset = 'plant500';
load 020_smallplant_wts
load 020_smallplant_p500
lambda = 0.08;
lambda2 = 0.008;
numtimes = 3;
% check_RCDM_parameters


%% plant 200
clear;
dataset = 'plant200';
addpath('../../data/data_reflections/')
load 020_smallplant_wts
load 020_smallplant_p200
lambda = 0.3;
lambda2 = 0.002;
numtimes = 3;
check_RCDM_parameters


%% oct with 200 superpixels

clear;
dataset = 'oct200';
load 001_oct_wts
load 001_oct_p200
lambda = 0.1;
lambda2 = 0.002;
numtimes = 3;
numtimes = 3;
check_RCDM_parameters


%% oct with 500 superpixels
clear;
dataset = 'oct500';
load 001_oct_wts
load 001_oct_p500
lambda = 0.01;
lambda2 = 0.008;
numtimes = 3;
check_RCDM_parameters
