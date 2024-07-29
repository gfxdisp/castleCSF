clc; close all; clear all;

if ~contains(path, '../')
    addpath(genpath('../'));
end

csf_model = CSF_castleCSF();

%%%%%%%%%%%%%% Plot opponent colour mechanism %%%%%%%%%%

csf_model.plot_mechanism('col_mech');


%%%%%%%%%%%%%% Plot temporal channels responses %%%%%%%%%%

csf_model.plot_mechanism('sust_trans');

%%%%%%%%%%%%%% Plot luminance peak sensitivity %%%%%%%%%%

csf_model.plot_mechanism('peak_s');

if ~contains(path, '../')
    rmpath(genpath('../'));
end