% get grouped data from raw individual data
% 
% close all; clear; clc;
% load environment_experiment
% 
% 
% %% Subject id by group
% 
% % group 1a: 20 perturbation / 0.5 noise
% % group 1b: 20 perturbation / 0.5 noise, triggers
% % group 2a: 20 perturbation / 2.0 noise
% % group 2b: 20 perturbation / 2.0 noise, triggers
% % group 3: 10 perturbation / 0.5 noise
% 
% % new subject numbers from 1 to 54
% ind1a = 1:11;
% ind1b = 12:22;
% ind2a = 23:33;
% ind2b = 34:44;
% ind3 = 45:55;
% ind1b_one = 56:66; % for one-model learner with 20-degree perturbation
% 
% n_sim = 66;
% conditions = zeros(n_sim, 1);
% conditions(ind1a) = 1;
% conditions(ind1b) = 2;
% conditions(ind2a) = 3;
% conditions(ind2b) = 4;
% conditions(ind3) = 5;
% conditions(ind1b_one) = 6;
% 
% %% Blocks
% Blocks = [20 60 40 50 20 30 40 50 120]; % Block lengths: baseline - ... - error-clamp
% indBase = [1:20]; % 20
% indP{1} = [21:80]; % 60 
% indP{2} = [121:170]; % 50
% indP{3} = [191:220]; % 30
% indP{4} = [261:310]; % 50
% indW{1} = [81:120]; % 40
% indW{2} = [171:190]; % 20
% indW{3} = [221:260]; % 40
% indEC = [311:430]; % 120 
%     
% 
% %% Parameters for Data analysis
% fnt = 10; % first and last n trials to calculate mean within each block
% npp = 5; % number of trials before/after trigger to calculate mean response
% maxRatio = 0.8; % upper thresholds to count trials
% minRatio = 0.2; % lower thresholds to count trials
% 
% %% Output variables
% N = n_sim; % number of subjects
% hand = zeros(N, lenTotal);
% firstlastMeans = zeros(N, 14); % 14 = 4 + 4 + 3 + 3; first LBs, last LBs, first WBs, last WBs
% trialsThresholds = zeros(N, 7); % 7 = 4 + 3; trs for LBs, trs for WBs
% prepostClamp = zeros(N, 4); % pre 1st-clamp, post 1st-clamp, pre 2nd-clamp, post 2nd-clamp


%% Main session
clf;% close all; 
clear; clc;

%/do a little simulation; clamp size will affect the relevance, smaller
%clamp, larger relevance;
pSize_array = 5:5:40;
% relevance_by_size = [0.9,0.7,0.6,0.5,0.4,0.4,0.4,0.4];
% pSize_array=2:7;
% relevance_by_size = [1,0.95,0.93,0.91,0.89,0.87];
% pSize_array = 15;
% relevance_by_size = [0.6];
relevance_by_size = linspace(0.9,0.4,length(pSize_array));
% relevance_by_size =1-1./(1+exp(-(pSize_array-40)/10));
% relevance_by_size = 0.999*ones(size(pSize_array));
%----
n_sim = length(pSize_array);
N = n_sim; % number of subjects
conditions(1:n_sim) = 1;
% matPerturbation(1,:) = [zeros(1,20) -1.*ones(1,200)]; %baseline and error clamp
% matCondition(1,:) = [ones(1,20) zeros(1,200)]; %1 for normal trials; 0 for clamp
matPerturbation(1,:) = [zeros(1,20) -1.*ones(1,80) zeros(1,80) -1.*ones(1,80)]; %baseline and error clamp
matCondition(1,:) = [ones(1,20) zeros(1,80) ones(1,80) zeros(1,80)]; %1 for normal trials; 0 for clamp
lenTotal = size(matPerturbation,2);
hand = zeros(N, lenTotal);

figure(1);
set(gcf, 'Position', [0 0 2200 1000], 'PaperPosition', [0 0 22 10]);
%ha = tight_subplot(2,3,0.01,[0.05 0.05],[0.03 0.03]);  % simulate 2 conditions, 3 subjects;
fs = 8;
rng(1)
seeds = randperm(n_sim);

for subject=1:n_sim
    % individual parameters script
    
%     individual_parameters
    individual_parameters_sim
    params_indv.relevance = relevance_by_size(subject);
    pSize = pSize_array(subject); %the size of perturbation in degrees

    
    if conditions(subject)==1
        gr = '1a';
        con = 1;
        rtn = matPerturbation(1, :);
    end
    
    rtn = pSize * rtn; %/make it real angle values
    
    params_exp.environment = matPerturbation(con,:) .* pSize;
    params_exp.condition = matCondition(con,:);
    
    % third parameter controls random seed: =0 for shuffle, >0 for random seed
%     output_variables = fMixExperts_KF(params_indv, params_exp, seeds(subject));
    output_variables = fModelAverage_KF(params_indv, params_exp, seeds(subject));

    hand = output_variables.h; %/make it a real angle unit; original model has data between 0 nad 1.
    h(subject, :) = hand;
    phat = output_variables.phat; %/ the estimated perturbation size 
    y = output_variables.y; %/ the estimated perturbation size 
    w = output_variables.w; %/ the estimated perturbation size 
    sp_err = output_variables.sp_err; %/ the estimated sensory prediction error
%     hand = pSize * output_variables.h; %/make it a real angle unit; original model has data between 0 nad 1.
%     h(subject, :) = hand;
%     phat = pSize * output_variables.phat; %/ the estimated perturbation size 
%     y = pSize * output_variables.y; %/ the estimated perturbation size 
%     w = output_variables.w; %/ the estimated perturbation size 
%     sp_err = pSize * output_variables.sp_err; %/ the estimated sensory prediction error

   
    % individual plot
    subplot(4,2,subject);
    %axes(ha(subject));
    set(gca,'fontsize',fs);
    hold on;
    %fill([310.5 430 430 310.5],[-10 -10 40 40],[0.85 0.85 0.85],'edgecolor','none');
    h1 = plot(rtn,'color','r','linewidth',0.5);
    h1.Color(4) = 1;
    h3 = plot(phat(1,:),'color','g','linewidth',0.5);  % default model
    %h4 = plot(-1.*phat(2,:),'color','b','linewidth',0.5);  % perturbation model
    h4 = plot(phat(2,:),'color','b','linewidth',0.5);  % perturbation model
    h5 = plot(y,'color','y','linewidth',0.5);  % sensory prediction
%     h6 = plot(w(1,:),'--','color','g','linewidth',0.5);  % perturbation model
%     h7 = plot(w(2,:),'--','color','b','linewidth',0.5);  % perturbation model
    h6 = plot(sp_err(1,:),'--','color','g','linewidth',0.5);  % perturbation model
    h7 = plot(sp_err(2,:),'--','color','b','linewidth',0.5);  % perturbation model
  %  text(340, 38, sprintf('subj %d',subject),'fontsize',fs,'fontweight','bold');
  %  text(340, 33, sprintf('%s',gr),'fontsize',fs,'fontweight','bold');
        h2 = plot(hand,'color','k','linewidth',0.5);
    h2.Color(4) = 0.6;

    xlim([1 length(hand)]);
    ylim([-40 40]);
    title(num2str(pSize));
end

% data structure to save
sdata.h = h;

%save mat_Sim sdata

% save mat_groupData Blocks indP indW indEC h ...
%     first_Ntrials_LB first_Ntrials_WB last_Ntrials_LB last_Ntrials_WB ...
%     prepostClamp matPerturbation ind1a ind1b ind2a ind2b ind3

% print('-dtiff','-r300','../Figures/individual_plots');

%sdata