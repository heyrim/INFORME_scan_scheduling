
function finalpoint = seqDesign_ZahidData_prior_firstdatapoint(nin)
global  modelType growthParams M V0 
RGB = get(groot,"FactoryAxesColorOrder");

load('Data_Zahid.mat')
modelType = 'PSI+CCR';

nplot = 0; 

M = 2000; 

types = [1, 1;     2, 1;     0, 0;     1, 2;     2, 2]; 

npriortype = types(nin, :); %[2,1]; % 1,1 lognormal delta; 1,2 lognormal delta K; 
                                    % 2,1 unif delta;      2,2 unif delta K  


growthParams = 0.33; % 


% Logitnormal for PSI 
pop_PSI_est = 0.86;  
sd_PSI_est = 0.65;  

deltascale =    2.2933; 
pop_delta_est = 0.017 *deltascale; 
sd_delta_est = 1.04; 

lb_PSI = 0.1; ub_PSI = 1.0; 
lb_delta = 0; ub_delta = 0.15; 

a = lb_PSI; b = ub_PSI; 

if( npriortype(1) == 1 )
logit_gamma  = log((pop_PSI_est(1) - a)/(b - (pop_PSI_est(1))));
guassian_num = normrnd(logit_gamma,sd_PSI_est(1)^2, [M,1]);
newChain1 = (b.*exp(guassian_num) + a)./(1+exp(guassian_num)); 

elseif( npriortype(1) == 2 )
    lb_PSI = 0.5; 
newChain1 = rand(M,1)*(ub_PSI-lb_PSI) + lb_PSI; 
pop_PSI_est = 0.5*(ub_PSI-lb_PSI) + lb_PSI; 

end


params = [growthParams, pop_PSI_est]; 


for npatient = 1:39 

disp( ['patient #' int2str(npatient)] )

fulldata.xdata = data.xdata{npatient}; 
fulldata.ydata = data.ydata{npatient}; 



%% First fit PSI or K using first two points 

V0 = fulldata.ydata( abs(fulldata.xdata) < 0.1 ); 

fitdata.xdata = fulldata.xdata(1:2); 
fitdata.ydata = fulldata.ydata(1:2); 

lb = [0.1 ];
ub = [1.0 ];
opt = optimset('Display','off');
[param_val,ss01] = fmincon(@(param)ssq_tumorVol(param,fitdata),pop_PSI_est,[],[],[],[],lb,ub, [], opt);
PSI_fit = param_val; 


K = V0 / PSI_fit; 
[time,vol] = ode23(@(t,y)tumorLogistic(t,y,[growthParams, K]), [fitdata.xdata(1),fitdata.xdata(end)], fitdata.ydata(1) );

% if( nplot ); figure(221); subplot( 5, 8, npatient); hold on; plot( time, vol ); plot( fulldata.xdata, fulldata.ydata, 'o' ); end 


%% compute the first data point to collect based on population prior distribution
if( npriortype(1) == 1 )
logit_gamma  = log((pop_delta_est - lb_delta)/(ub_delta - (pop_delta_est)));
guassian_num = normrnd(logit_gamma,sd_delta_est^2, [M,1]);
newChain = (ub_delta.*exp(guassian_num) + lb_delta)./(1+exp(guassian_num)); 

elseif( npriortype(1) == 2 )
newChain = rand(M,1)*ub_delta; 
pop_delta_est = 0.5*ub_delta; 

end

% low fidelity model 
params = [growthParams, PSI_fit, pop_delta_est]; 

[timeFit, volFit] = headnecktumorPSI_RT( params, fulldata ); 
expDesigns = 1:min( 6, round(fulldata.xdata(end)) ); 
lowfi = interp1(timeFit,volFit, expDesigns);


% if( nplot ); figure(221); subplot( 5, 8, npatient ); hold on; plot( timeFit, volFit, '--' ); end 


% Calculate predicted values at each experimental design with
% each parameter set from newChain
for ii = 1:size(newChain,1)

    if( npriortype(2) == 1 )
    %%% only prior in delta 
    [tsol, ysol] = headnecktumorPSI_RT([growthParams, PSI_fit, newChain(ii)],fulldata);

    elseif( npriortype(2) == 2 )
    %%%%  prior in PSI and delta 
    [tsol, ysol] = headnecktumorPSI_RT([growthParams, newChain1(ii), newChain(ii)],fulldata);

    end

    for jj = 1:length(expDesigns)
        lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
    end
end

%Calculate MI for each remaining design
%%%%% with normalized by mean/std data
normChain = (newChain-mean(newChain))./std(newChain);
normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);

if( npriortype(2) == 2 )
    normChain1 = (newChain1-mean(newChain1))./std(newChain1);
    normChain = [normChain1, normChain]; 
end

for jj = 1:length(expDesigns)
    [I1,~] = KraskovMI(normChain,normlowfiOut(:,jj),6);
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
%|currentPt-lowFiEndPrediction|/(currentPt+lowFiEndPrediction)

rcv = abs(lowfi(end)-fitdata.ydata(end))/(lowfi(end)+fitdata.ydata(end)); 

%penalize MI by k*absolute rcg*penalty for skipped points
for p = 1:length(relMI)
    score(p) = unifMI(p) - abs(rcv)*sum(unifMI(1:(p-1)))/sum(unifMI);
end

%Choose optimal design
point = expDesigns(1,min(find(score == max(score)),size(expDesigns,2)))

if length(point)>1
    fprintf('Multiple points have same score function. Last point in tie list chosen.')
    point = point(end);
end


point_save(npatient) = [point]; 


if( nplot ) 
figure(314); subplot( 5, 8, npatient); hold on; plot( unifMI, 'x' ); plot( score, 'o' ); 
ii = find(score == max(score)); 
plot( [ii ii], [0 1], 'r:')
end 
if( nplot ); figure(221); subplot( 5, 8, npatient ); hold on; title( int2str(ii) ); end 

figure(1521); hold on; plot( score, ':', 'color', RGB(nin,:) ); 

end 

finalpoint = mode( point_save ); 




% % use a uniform prior 
% newChain = rand( M, 1 )*0.15; 
% pop_delta_est = 0.15/2; 
% params = [growthParams, PSI_fit, pop_delta_est]; 

end 


function dval = tumorLogistic( t, val, params )

lambda = params(1);
K = params(2);

V = val(1);

dval(1) = lambda * V * ( 1 - V / K );

end


% SSQ function for calibration before RT 
function SSrt = ssq_tumorVol(params,data)
global growthParams V0 

K = V0 / params(1); 

[time,vol] = ode23(@(t,y)tumorLogistic(t,y,[growthParams, K]), [data.xdata(1),data.xdata(end)], data.ydata(1) ); 

tumVol = interp1(time,vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end
