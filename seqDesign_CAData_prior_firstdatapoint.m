
function point = seqDesign_CAData_prior_firstdatapoint(nin) 
global M growthParams dose 

load('/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/ABM_cancer_radiotherapy/data/CA_rad1_Cm_22_30_pNR_all.mat')
modelType = 'LOG+DVR'; 
rcv_vec = []; 
RGB = get(groot,"FactoryAxesColorOrder");
dose = 2; % RT dose 
M = 1000; % chain lengths

types = [1, 1;     2, 1;     0, 0;     1, 2;     2, 2]; 

npriortype = types(nin, :); %[2,1]; % 1,1 lognormal delta; 1,2 lognormal delta K; 
                                    % 2,1 unif delta;      2,2 unif delta K  
% result :    27    22     0    27    22

%% prior distribution using Monolix using all 135 data
P = readtable( '../data/Update files 02-28-2025/prostate cancer model - large - general - bounded/populationParameters.txt');

pop_est   = P{1:3,2}; %population level estimates for A, B, gamma
sd_est    = P{4:6,2}; %standard deviation of random effect

newChain = zeros( M, 3 );
newChain(:,1) = pop_est(1); 
newChain(:,2) = lognrnd( log( pop_est(2) ), sd_est(2), [M,1] );

% bounded logit normal for beta 
betascale = 0.0211; 

a = eps; b = 10; %range
logit_gamma  = log((pop_est(3) - a)/(b - (pop_est(3))));
guassian_num = normrnd(logit_gamma,sd_est(3), [1,M]);
% guassian_num = (b.*exp(guassian_num) + a)./(1+exp(guassian_num))/120; 
guassian_num = (b.*exp(guassian_num) + a)./(1+exp(guassian_num))*betascale;  
newChain(:,3) = guassian_num; clear guassian_num

% guassian_num = normrnd(logit_gamma,sd_est(3), [1,1000000]);
% beta_val = mean( (b.*exp(guassian_num) + a)./(1+exp(guassian_num))/120 ); clear guassian_num

beta_val = pop_est(3)*betascale; 
%%%% precomputed beta_val from above is    beta_val = 0.0197; 0.0333; 

if( npriortype(1) == 2 )
% %% uniform prior 
newChain(:,1) = pop_est(1); 

lb_PSI = 0.2; ub_PSI = 1.0; 
newChain(:,2) = rand(M,1)*(ub_PSI-lb_PSI) + lb_PSI; 
pop_est(2) = 0.5*(ub_PSI-lb_PSI) + lb_PSI; 

newChain(:,3) = rand( [M,1] )*b*betascale;
pop_est(3) =  0.5*b*betascale;
beta_val = pop_est(3); 
end 

for npatient = [1:8:135,135] 

fulldata.xdata = data.xdata(1:55);
fulldata.ydata = data.ydata(1:55,npatient); %use only tumor volume during treatment period


%Initial starting data - days 15 were used to get intrinsic growth
%params, then add day 19 to start the RT parameter estimation
currentPts = [ fulldata.xdata(16) fulldata.ydata(16);];
icdata = currentPts;

growthParams = pop_est(1:2)';
params = [growthParams, beta_val ];


%% compute the first data point to collect based on prior distribution
[timeFit, volFit] = tumorModel_postRT_wIC(params,icdata);
lowfi = interp1(timeFit,volFit, fulldata.xdata);

expDesigns = fulldata.xdata((17):end);

% Calculate predicted values at each experimental design with
% each parameter set from newChain
for ii = 1:size(newChain,1)

    [tsol, ysol] = tumorModel_postRT_wIC([newChain(ii,:)],icdata);

    for jj = 1:length(expDesigns)
        lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
    end
end

%Calculate MI for each remaining design
%%%%% with normalized by mean/std data
normChain = (newChain-mean(newChain))./std(newChain);
normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);

miVals = zeros( 1, length(expDesigns)); 
for jj = 1:length(expDesigns)
    % [I1,~] = KraskovMI(normChain,normlowfiOut(:,jj),6); % only beta
    [I1,~] = KraskovMI(normChain,lowfiOut(:,jj),6); % only beta
    
    miVals(jj) = I1;
end
maxMI = max(miVals);
relMI = miVals/maxMI;
lbMI = min(miVals);
ubMI = max(miVals);
if lbMI==ubMI
    unifMI = ones(1,length(miVals));
else
    unifMI = (miVals-lbMI)/(ubMI-lbMI);
end

%Calculate score function by penalizing for points skipped
rcv = abs(lowfi(end)-currentPts(end,2))/(lowfi(end)+currentPts(end,2))*2;

%penalize MI by k*absolute rcg*penalty for skipped points
score = zeros(1, length(relMI)); 
for p = 1:length(relMI)
    %(it is choosing every point because MIs are all nearly
    %identical, so any penalty means choosing the first point)
    score(p) = unifMI(p) - abs(rcv)*sum(unifMI(1:(p-1)))/sum(unifMI);
end

%Choose optimal design
point = expDesigns(1,find(score == max(score))) 

point_save(npatient) = point;

figure(152); hold on; 
plot( expDesigns, score, ':', 'color', RGB(nin,:) ); 


clear relMI miVals lowfiOut % score 

end 


point = mode(point_save([1:8:135,135])); 

end 
