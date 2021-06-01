addpath('../../data/data_reflections/')

mex RCDM_cversion.cpp;
mex RCDM_cversion_greedypar.cpp;
mex ACDM_cversion.cpp;
mex APcompact_cversion.cpp;
mex APfull_cversion.cpp;

%% smallplant with 500 superpixels
clear;
dataset = 'plant500';
load 020_smallplant_wts
load 020_smallplant_p500
lambda = 0.08;
lambda2 = 0.008;

c = 25
run_one_image_exp

%% oct with 500 superpixels
clear;
dataset = 'oct500';
load 001_oct_wts
load 001_oct_p500
lambda = 0.01;
lambda2 = 0.008;

c = 25
run_one_image_exp

%% smallplant with 200 superpixels

clear;
dataset = 'plant200';
addpath('../../data/data_reflections/')
load 020_smallplant_wts
load 020_smallplant_p200
lambda = 0.3;
lambda2 = 0.002;

c = 25
run_one_image_exp

%% oct with 200 superpixels

clear;
dataset = 'oct200';
load 001_oct_wts
load 001_oct_p200
lambda = 0.1;
lambda2 = 0.002;

c = 25
run_one_image_exp
