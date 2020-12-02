% individual parameters all simulations paper. Youngmin multiple model paper. 
% Nicolas Schweighofer 10.30.17

% latest: 3.30.2019
% params_indv.prior = [0.9, 0.1];
% params_indv.retention = [0.98 0.999];
% params_indv.sigma_u = 0.13;
% 
% params_indv.initial_S = [0.25 0.8].^2;
% params_indv.initial_S_one_model = [0.4]^2;
% params_indv.minimum_S = params_indv.initial_S(1);
% 
% params_indv.KF_R = [1.1 2].^2;
% params_indv.KF_Q = [0.04 0.04].^2;


% updated: 4.22.2019
% base_prior = 0.9;
% params_indv.prior = [base_prior, 1-base_prior];
% params_indv.retention = [0.97 0.999];
% params_indv.sigma_u = 0.14;
% 
% params_indv.initial_S = [0.25 0.7].^2;
% params_indv.initial_S_one_model = [0.4]^2;
% params_indv.minimum_S = params_indv.initial_S(1);
% 
% params_indv.KF_R = [0.8 2].^2;
% params_indv.KF_Q = [0.04 0.04].^2;


% updated 04/25/2019
%base_prior = 0.95;
base_prior = 0.9;  %/assume which model is more reliable?
params_indv.prior = [base_prior, 1-base_prior];
% params_indv.retention = [0.96 0.9997];
params_indv.retention = [0.96 0.9];

% params_indv.sigma_u = 0.17;
params_indv.sigma_u = 0.17;

% params_indv.initial_S = [0.2 2].^2;
params_indv.initial_S = [20 4].^2; %/assume proprioception has more noise
params_indv.initia_S_one_model = [0.4]^2;
params_indv.minimum_S(1) = params_indv.initial_S(2);%/// assume minimum S is the same as visual noise
params_indv.minimum_S(2) = params_indv.initial_S(2);
% params_indv.KF_R = [0.8 5].^2;
params_indv.KF_R = [80 20].^2; %/assume proprioception has larger measurement noise
params_indv.sigma_s = 0.05 * 20;
params_indv.KF_Q = [params_indv.sigma_s, params_indv.sigma_s].^2;


