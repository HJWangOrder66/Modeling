% Mixture of Experts model
% Kalman Filter version with two experts (baseline and perturbation models)
% Updated: 4.3.2019

function output_variables = fModelAverage_KF(params_indv, params_exp, random_seed)

% control random seed
if random_seed > 0
    rng(random_seed)
else
    rng('shuffle', 'v5normal')
end


% learning schedule
p = params_exp.environment; % perturbation magnitude (c=1) or error-clamp angle (c=0)
c = params_exp.condition; % 1 for ordinary trial (learning or washout), 0 for error-clamp
nTrial = length(p); % number of total trial


% model parameters
relevance = params_indv.relevance;
prior = params_indv.prior; % prior
initS = params_indv.initial_S; % initial uncertainty (var)
minS1 = params_indv.minimum_S(1); % minimum uncertainty (var) for two models
minS2 = params_indv.minimum_S(2);
a1 = params_indv.retention(1); % retention parameter of baseline model
a2 = params_indv.retention(2); % retention parameter of perturbation

sigma_u = params_indv.sigma_u; % motor noise level (std)
process_noiseQ1 = params_indv.KF_Q(1); % process noise level of baseline model (var)
process_noiseQ2 = params_indv.KF_Q(2); % process noise level of perturbation model (var)
measure_noiseR1 = params_indv.KF_R(1); % measurement noise level of baseline model (var)
measure_noiseR2 = params_indv.KF_R(2); % measurement noise level of perturbation model (var)


% variables
cue = zeros(2, nTrial); % model selection cue
u = zeros(1, nTrial); % motor command
h = zeros(1, nTrial); % hand direction
y = zeros(1,nTrial); % visual cursor position
phat = zeros(2, nTrial); % perturbation estimation
yhat = zeros(2, nTrial); % sensory prediction
sp_err = zeros(2, nTrial); % sensory prediction error
S = zeros(2, nTrial); % uncertainty of model
pdf = zeros(2, nTrial); % probability distribution function (likelihood)
pos = zeros(2, nTrial); % posterior
w = zeros(2, nTrial); % weight


% initial variable assignment
S(:, 1) = [initS(1); initS(2)]; %/uncertianty of both models, variance sigma^2
pos(:,1) = [prior(1); prior(2)]; %/position of both models, initially
w(:, 1) = [prior(1); prior(2)];  %/weight of both models
cue(:, 1) = [1; 0];  %/perturbation model vs baseline model


% simulation loop
for t=2:nTrial
    
    % action selection
    u(t) = 0 - phat(1,t-1);  % for error clamp trials, hand move according to baseline model
    % u(t) = 0 - cue(:,t-1)'*phat(:,t-1); % target at 0; %/perturbation at -1, cue = [1,0], phat starts from [0,0];
    u(t) = u(t) + sigma_u*randn;
    h(t) = u(t);
    
    % prediction and feedback
    %      yhat(1,t) = phat(1,t-1); % for the baseline model, it predict where it moves to (i.e., u(t));
    %      yhat(2,t) = phat(2,t-1); %for visual clamp, its prediction is unchanged from previous trials
    yhat(:,t) = phat(:,t-1) + u(t)*ones(2,1);
    %    y(t) = c(t)*(p(t) + h(t)) + (1 - c(t))*p(t);%/c = 1 for normal trial, c=0 for error clamp; p(t) is [0 -1];
    %    sp_err(:,t) = y(t)*ones(2,1) - yhat(:,t);
    %    y(t) = w(1,t-1).*yhat(1,t) + w(2,t-1).*yhat(2,t); % the actual "feedback" is the weighted average.
    %     y(t) = (1-c(t)) * (relevance*p(t) + (1-relevance)*u(t)) + c(t) * (p(t) + h(t)); % the actual "feedback" is the weighted average.
    y(t) = (1-c(t)) * (relevance * p(t) + (1-relevance) * u(t)) + c(t) * (p(t) + h(t)); % the actual "feedback" is the weighted average.
    
    %     y(t) = w(1,t-1).*p(t) + w(2,t-1).*u(t); % the actual "feedback" is the weighted average.
    sp_err(:,t) = y(t)*ones(2,1) - yhat(:,t);
    
    % model selection
    pdf(1,t) = exp(-sp_err(1,t)^2/(2*S(1,t-1))) / sqrt(2*pi*S(1,t-1));
    pdf(2,t) = exp(-sp_err(2,t)^2/(2*S(2,t-1))) / sqrt(2*pi*S(2,t-1));
    pos(1,t) = prior(1) * pdf(1,t);
    pos(2,t) = prior(2) * pdf(2,t);  %/posterial estimate of probability, Note prior is a constant
    w(:,t) = [pos(1,t); pos(2,t)] ./ (pos(1,t) + pos(2,t));
    %cue(:,t) = [w(1,t)>=w(2,t); w(1,t)<w(2,t)]; % winner-take-all, %/ the two models are given either 0 or 1 weight.
    %instead of using winner-takes-all, we use  weights to get the average
    
    
    % model update
    Ks1 = S(1,t-1)/(S(1,t-1) + measure_noiseR1); %/kalman learning rate; measurement noise is a constant
    Ks2 = S(2,t-1)/(S(2,t-1) + measure_noiseR2);
    
    % there is no winner, simply update both models
    if c(t-1)==1
    phat(1,t) = a1*phat(1,t-1) + Ks1*sp_err(1,t);
    S(1,t) = ((1-Ks1)*S(1,t-1)); %/update estimation variance of the selected model;
    elseif c(t-1)==0
%             phat(1,t) = a1*phat(1,t-1);
%     S(1,t) = ((a1^2)*S(1,t-1) + process_noiseQ1);
    phat(1,t) = a1*phat(1,t-1) + Ks1*sp_err(1,t);
        S(1,t) = S(1,t-1);
    end
%     phat(1,t) = a1*phat(1,t-1);
%     S(1,t) = ((a1^2)*S(1,t-1) + process_noiseQ1);
    phat(2,t) = a2*phat(2,t-1) + Ks2*sp_err(2,t);
    S(2,t) = ((1-Ks2)*S(2,t-1));
    %     if cue(1,t)==1  %/baseline model wins
    %         phat(1,t) = a1*phat(1,t-1) + Ks1*sp_err(1,t);
    %         phat(2,t) = a2*phat(2,t-1);
    %         S(1,t) = ((1-Ks1)*S(1,t-1)); %/update estimation variance of the selected model;
    %         S(2,t) = ((a2^2)*S(2,t-1) + process_noiseQ2);  %/if not selected, variance will continue to changed by forgetting and noise.
    %     else
    %         phat(2,t) = a2*phat(2,t-1) + Ks2*sp_err(2,t);
    %         phat(1,t) = a1*phat(1,t-1);
    %         S(2,t) = ((1-Ks2)*S(2,t-1));
    %         S(1,t) = ((a1^2)*S(1,t-1) + process_noiseQ1);
    %     end
    
    
    if S(1,t) < minS1  %/the estimation variance has a floor
        S(1,t) = minS1;
    end
    if S(2,t) < minS2
        S(2,t) = minS2;
    end
    
    % Distributions
    angle = -4:0.01:4;  %/obtain a distribution as a function of angle.
    likelihood(:,1,t) = exp((-(angle-phat(1,t)).^2)/(2*S(1,t-1))) ./ sqrt(2*pi*S(1,t-1));
    likelihood(:,2,t) = exp((-(angle-phat(2,t)).^2)/(2*S(2,t-1))) ./ sqrt(2*pi*S(2,t-1));
    posterior(:,1,t) = prior(1)*likelihood(:,1,t);
    posterior(:,2,t) = prior(2)*likelihood(:,2,t);
    sumPost(:,1,t) = posterior(:,1,t) + posterior(:,2,t);
    pdfDist(:,1,t) = posterior(:,1,t) ./ sumPost(:,1,t);
    pdfDist(:,2,t) = posterior(:,2,t) ./ sumPost(:,1,t);
    
end

Weight1 = pdfDist(end:-1:1,1,:);
Weight1 = reshape(Weight1,801,nTrial);
Weight2 = pdfDist(end:-1:1,2,:);
Weight2 = reshape(Weight2,801,nTrial);
Like1 = likelihood(end:-1:1,1,:);
Like1 = reshape(Like1,801,nTrial);
Like2 = likelihood(end:-1:1,2,:);
Like2 = reshape(Like2,801,nTrial);
Pos1 = posterior(end:-1:1,1,:);
Pos1 = reshape(Pos1,801,nTrial);
Pos2 = posterior(end:-1:1,2,:);
Pos2 = reshape(Pos2,801,nTrial);

Like = (Like1 + Like2)/2;

output_variables.y=y;
output_variables.h = h;
output_variables.phat = phat;
output_variables.S = S;
output_variables.cue = cue;
output_variables.pdf = pdf;
output_variables.pos = pos;
output_variables.w = w;
output_variables.sp_err= sp_err;
output_variables.Like = Like;
output_variables.Like1 = Like1;
output_variables.Like2 = Like2;
output_variables.Like = Like;
output_variables.Pos1 = Pos1;
output_variables.Pos2 = Pos2;
output_variables.Weight1 = Weight1;
output_variables.Weight2 = Weight2;

end
