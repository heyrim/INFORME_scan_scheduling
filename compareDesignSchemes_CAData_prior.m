%Compare errors, uncertainties, and fits for design schemes:
% - initial 1, 9, 16 &  one point per week, each Friday 

% - prior with 16 & one point per week, each Friday  

% - prior with 16 & "optimal" six points (algorithm choice)

%%%%% - initial 1, 9, 16 & "optimal" six points (algorithm choice)

%%%%% - prior with 16 & one point per week, each Friday, throughout treatment

%%%% unif prior 
schedule_unif = [ 16, 23, 24, 30]; %, 33, 40]; 
%%%% population prior 
schedule_pop = [ 16, 28, 29, 49 ];  
%%%% group prior 
schedule_group{1} = [ 16, 28, 29 ]; %[ 16, 28, 29, 42 ];
schedule_group{2} = [ 16, 28, 30, 49 ];
schedule_group{3} = [ 16, 28, 32 ]; %[ 16, 28, 32, 34 ];

% scan number 
% mean of unif2.3556
% mean of pop2.1926
% mean of group2.0963 - each group  2.0192e+00,   2.2500e+00,    2.0256e+00



% clear all
% close all

%% load data 
load('/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/ABM_cancer_radiotherapy/data/CA_rad1_Cm_22_30_pNR_all.mat')
% exist npatient
% if( ~ans );  npatient = 86;  end  

fulldata.xdata = data.xdata(1:55)';
fulldata.ydata = data.ydata(1:55,npatient); %use only tumor volume during treatment period


path = 'Figures/LognormalPrior_250602_Comparison/'; 
mkdir(path);

plotind = [1 1 1 1 0 1]; 
plotind = [0 0 0 0 0 1]; 

%%%% First calibrate lambda and K using first three points 1 9 16 
% for nn = 1:length( indnp )
% data.xdata = fulldata.xdata( [1 9 16] ); 
% data.ydata = fulldata.ydata( [1 9 16], indnp(nn) ); 
% 
% lb = [0 0];
% ub = [1 1];
% opt = optimset('Display','off');
% [param_val,ss01] = fmincon(@(param)ssq_log(param,data),[0.1, 0.5],[],[],[],[],lb,ub, [], opt);
% 
% growthparamfit(:,indnp(nn)) = param_val'; 
% 
% % figure(156); subplot( 10, 10, nn ); hold on; plot( data.xdata, data.ydata, 'o' ); 
% % [t,y] = ode23(@(t,y)tumorLogistic(t,y,param_val),[data.xdata(1), data.xdata(end)],data.ydata(1));
% % plot( t, y )
% 
% end 

global growthParams

obsdata = [16 17:55];

firstsix = [17:22];
scheme3pts = [fulldata.xdata(firstsix) fulldata.ydata(firstsix)]; %Days 16-21 (first six)

% scheme2pts = point_list(2:7,:); %Chosen points from algorithm

weeklyIdx = [1 9 16 20 27 34 41 48 55]; %  % (each Friday, once per week)
scheme1pts = [fulldata.xdata(weeklyIdx) fulldata.ydata(weeklyIdx)];

modelType = 'LOG+DVR';
no_smps = 10000;



%% First, calibrate with initial 1, 9, 16 + six scheme_pts - equi-spaced 
if plotind(1) 

scheme_pts = scheme1pts; % weekly 
growthParams = growthparamfit(:,npatient)'; 
IC = [fulldata.xdata(1), fulldata.ydata(1)]; 

params1 = {
    {'beta',0.01,0,1}
    };


for nn = 1:6 
data.xdata = scheme_pts(1:(3+nn),1);
data.ydata = scheme_pts(1:(3+nn),2);


model.ssfun = @(params,data) ssq_compareSchemes(params,data,IC,modelType,growthParams);
options.updatesigma = 1;
options.waitbar = 0; 

options.nsimu = no_smps/5;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);


% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
beta_val(nn) = chain(ind,:); %This is our fitted parameter set
beta_valMean = mean(chain); %This is our mean value
params = [growthParams beta_val(nn)];


% Generate current model trajectory, credible intervals, and plot
[timeFit, volFit] = tumorModel_compareSchemes(modelType,params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata);

ln = fulldata.xdata;
modelfun = @(ln,params) tumorcred(ln,params,modelType,growthParams,IC);

figure(1); subplot( 2, 3, nn ); if( nn==1 ); title( strcat( 'Patient #', int2str(npatient) ) ); end 
pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
pred.obslims = [];
mcmcpredplot(pred)
hold on
plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
plot(fulldata.xdata,lowfi,'-k','Linewidth',1) %This is your optimal model fit
plot(data.xdata,data.ydata(:,1),'ok','MarkerFaceColor','k','Linewidth',1,'MarkerSize',8) %Selected scans
xlabel('Time (days)','FontSize',12);
ylabel('Volume','FontSize',12);
axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
set(gca, 'FontSize',12)
hold off
filename = [path strcat('ModelFit_Scheme1_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)
% close(gcf)


% figure(2); subplot( 2, 3, nn ); hold on; 
% mcmcplot(chain,[],[],'chainpanel')
% filename = [path strcat('Chain_Scheme1_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% % close(gcf)
% 
% figure(3); hold on; 
% mcmcplot(chain,[],[],'denspanel')
% filename = [path strcat('Density_Scheme1_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% % close(gcf)


% Calculate MSE to measure error
mse_scheme1(nn) = sum((lowfi(16:end)-fulldata.ydata(16:end)).^2)/numel(fulldata.ydata(16:end));


% Calculate uncertainty metrics
for n = 1:length(ln)
    widthCI(n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
end
areaCI_scheme1(nn) = sum(widthCI(1:end));

end 

save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], 'mse_scheme1', 'areaCI_scheme1', 'scheme1pts' ) 

close all 

end 



%% Second, calibrate with prior + six scheme_pts - equi-spaced 
if plotind(2) 

P = readtable( './Update files 02-28-2025/prostate cancer model - large - general - bounded/populationParameters.txt');

pop_est   = P{1:3,2}; %population level estimates for A, B, gamma
sd_est    = P{4:6,2}; %standard deviation of random effect

% M = 10000000; 
% a = eps; b = 10; %range
% logit_gamma  = log((pop_est(3) - a)/(b - (pop_est(3))));
% guassian_num = normrnd(logit_gamma,sd_est(3), [1,M]);
% guassian_num = (b.*exp(guassian_num) + a)./(1+exp(guassian_num))/120; 
% fitdist( guassian_num', 'lognormal' )
mu_est(3)  = -4.500; 
sig_est(3) = 1.237; 

%prior lognrnd( log( pop_est(i) ), sd_est(i), [M,1] );

growthParams = pop_est(1:2)'; 

IC0 = [fulldata.xdata(1), fulldata.ydata(1)]; 
IC = [fulldata.xdata(16), fulldata.ydata(16)]; 

% load( strcat( '/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Figures/LognormalPrior_250513/Results_nPatient', int2str(npatient), '.mat') )


% weeklyIdx = [16 20 27 34 41 48 55]; %  % (each Friday, once per week)
for mm = 17:20 

    weeklyIdx = [16, mm:7:55]; 

scheme2pts = [fulldata.xdata(weeklyIdx) fulldata.ydata(weeklyIdx)];

params1 = {
    % {'beta',0.1,0,1}
    {'beta',0.01, eps, 1, mu_est(3), sig_est(3) }
    };


for nn = 1:6  

data.xdata = scheme2pts(1:(1+nn),1);
data.ydata = scheme2pts(1:(1+nn),2);

if( nn == 2 ) % update prior 
    yratio = lowfi(end)/lowfi(1);
    if( yratio <= 0.2 );     nprior = 3;
    elseif( yratio <= 0.85 );     nprior = 4;
    else;     nprior = 5;
    end

    % P = readtable( './Update files 02-27-2025/Monolix and MATLAB - fit and distribution/prostate cancer model - large - Cov/populationParameters.txt');
    mu_est(3:5) = [  -3.3816e+00  -4.0682e+00  -5.4838e+00]; 
    sig_est(3:5) =  [ 1.960e-01   1.960e-01   1.960e-01 ]; 

params1 = {
    {'beta',beta_val(1), eps, 1, mu_est(nprior), sig_est(nprior) }
    };

end 

model.ssfun = @(params,data) ssq_compareSchemes(params,data,IC,modelType,growthParams);

%%%% beta prior 
model.priorfun = @lognpdf;

options.updatesigma = 1;
options.waitbar = 0; 


options.nsimu = no_smps/5;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);



% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
beta_val(nn) = chain(ind,:); %This is our fitted parameter set
beta_valMean = mean(chain); %This is our mean value
params = [growthParams beta_val(nn)];


% Generate current model trajectory, credible intervals, and plot
[timeFit, volFit] = tumorModel_compareSchemes(modelType,params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

ln = fulldata.xdata(16:end);
modelfun = @(ln,params) tumorcred(ln,params,modelType,growthParams,IC);



figure(4); subplot( 2, 3, nn ); if( nn==1 ); title( strcat( 'Patient #', int2str(npatient) ) ); end 
pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
pred.obslims = [];
mcmcpredplot(pred)
hold on
plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',1) %This is your optimal model fit
plot(data.xdata,data.ydata(:,1),'ok','MarkerFaceColor','k','Linewidth',1,'MarkerSize',8) %Selected scans
xlabel('Time (days)','FontSize',12);
ylabel('Volume','FontSize',12);
axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
set(gca, 'FontSize',12)
hold off
filename = [path strcat('ModelFit_Scheme2_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)
% close(gcf)


% figure(2); subplot( 2, 3, nn ); hold on; 
% mcmcplot(chain,[],[],'chainpanel')
% filename = [path strcat('Chain_Scheme2_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% % close(gcf)
% 
% figure(3); hold on; 
% mcmcplot(chain,[],[],'denspanel')
% filename = [path strcat('Density_Scheme2_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% % close(gcf)


% Calculate MSE to measure error
mse_scheme2(nn) = sum((lowfi-fulldata.ydata(16:end)).^2)/numel(lowfi);

mse_scheme2_all(nn,mm) = mse_scheme2(nn); 

% Calculate uncertainty metrics
for n = 1:length(ln)
    widthCI(n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
end
areaCI_scheme2(nn) = sum(widthCI(1:end));


areaCI_scheme2_all(nn,mm) = areaCI_scheme2(nn); 

end 

save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'mse_scheme2', 'areaCI_scheme2', 'scheme2pts', 'mse_scheme2_all', 'areaCI_scheme2_all' ) 

close all 

end 

clear beta_val 

end 

if( sum(plotind(3:6)) )

P = readtable( '../data/Update files 02-28-2025/prostate cancer model - large - general - bounded/populationParameters.txt');

pop_est   = P{1:3,2}; %population level estimates for A, B, gamma
sd_est    = P{4:6,2}; %standard deviation of random effect

% M = 10000000; 
% a = eps; b = 10; %range
% logit_gamma  = log((pop_est(3) - a)/(b - (pop_est(3))));
% guassian_num = normrnd(logit_gamma,sd_est(3), [1,M]);
% guassian_num = (b.*exp(guassian_num) + a)./(1+exp(guassian_num))/120; 
% fitdist( guassian_num', 'lognormal' )
mu_est(3)  = -4.500; 
sig_est(3) = 1.237; 

%prior lognrnd( log( pop_est(i) ), sd_est(i), [M,1] );

growthParams = pop_est(1:2)'; 

IC0 = [fulldata.xdata(1), fulldata.ydata(1)]; 
IC = [fulldata.xdata(16), fulldata.ydata(16)]; 

end 



%% Third, calibrate with prior + chosen six (or less) points
if plotind(3)

load( strcat( '/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Figures/LognormalPrior_250531_betanK/Results_nPatient', int2str(npatient), '.mat') )

scheme3pts = point_list(1:min(size(point_list,1),7),1); 
scheme3pts(:,2) = fulldata.ydata(point_list(1:min(size(point_list,1),7),1)+1);

params1 = {
    % {'beta',0.1,0,1}
    {'beta',0.01, eps, 1, mu_est(3), sig_est(3) }
    };


for nn = 1:(min(size(point_list,1),7)-1) 

data.xdata = scheme3pts(1:(1+nn),1);
data.ydata = scheme3pts(1:(1+nn),2);

if( nn == 2 ) % update prior 
    yratio = lowfi(end)/lowfi(1);
    if( yratio <= 0.2 );     nprior = 3;
    elseif( yratio <= 0.85 );     nprior = 4;
    else;     nprior = 5;
    end

    % P = readtable( './Update files 02-27-2025/Monolix and MATLAB - fit and distribution/prostate cancer model - large - Cov/populationParameters.txt');
    mu_est(3:5) = [  -3.3816e+00  -4.0682e+00  -5.4838e+00]; 
    sig_est(3:5) =  [ 1.960e-01   1.960e-01   1.960e-01 ]; 

params1 = {
    {'beta',beta_val(1), eps, 1, mu_est(nprior), sig_est(nprior) }
    };

end 

model.ssfun = @(params,data) ssq_compareSchemes(params,data,IC,modelType,growthParams);

%%%% beta prior 
model.priorfun = @lognpdf;

options.updatesigma = 1;
options.waitbar = 0; 


options.nsimu = no_smps/5;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);



% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
beta_val(nn) = chain(ind,:); %This is our fitted parameter set
beta_valMean(nn) = mean(chain); %This is our mean value
params = [growthParams beta_val(nn)];


% Generate current model trajectory, credible intervals, and plot
[timeFit, volFit] = tumorModel_compareSchemes(modelType,params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

ln = fulldata.xdata(16:end);
modelfun = @(ln,params) tumorcred(ln,params,modelType,growthParams,IC);



figure(4); subplot( 2, 3, nn ); if( nn==1 ); title( strcat( 'Patient #', int2str(npatient) ) ); end 
pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
pred.obslims = [];
mcmcpredplot(pred)
hold on
plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',1) %This is your optimal model fit
plot(data.xdata,data.ydata(:,1),'ok','MarkerFaceColor','k','Linewidth',1,'MarkerSize',8) %Selected scans
xlabel('Time (days)','FontSize',12);
ylabel('Volume','FontSize',12);
axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
set(gca, 'FontSize',12)
hold off
filename = [path strcat('ModelFit_Scheme3_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)
% close(gcf)


% figure(2); subplot( 2, 3, nn ); hold on; 
% mcmcplot(chain,[],[],'chainpanel'); 
% filename = [path strcat('Chain_Scheme3_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% % close(gcf)

figure(3); hold on; 
yy = mcmcplot(chain,[],[],'denspanel'); 
filename = [path strcat('Density_Scheme3_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)
% close(gcf)

pdf_vec{nn} = yy{1};


% Calculate MSE to measure error
mse_scheme3(nn) = sum((lowfi-fulldata.ydata(16:end)).^2)/numel(lowfi);


% Calculate uncertainty metrics
for n = 1:length(ln)
    widthCI(n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
end
areaCI_scheme3(nn) = sum(widthCI(1:end));

end 

save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'mse_scheme3', 'areaCI_scheme3', 'scheme3pts', 'beta_val', 'beta_valMean', 'pdf_vec', 'nprior'  ) 

% close all 

end 



%% Fourth, calibrate with prior + chosen six (or less) points, until beta conv. 
%%                              if finished early, scan last week Friday 
if plotind(4) 

load( strcat( '/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Figures/LognormalPrior_250531_betanK/Results_nPatient', int2str(npatient), '.mat') )

tmp = abs( beta_val(1:end-1) - beta_val(2:end) ); tmp = [1, tmp]; 
thresh_beta = 0.003; 
indscan = find( tmp < thresh_beta ) ;  
indscan = min( min(indscan), length(scheme3pts) ) ;

mse_scheme4 = mse_scheme3(1:indscan); 
areaCI_scheme4 = areaCI_scheme3(1:indscan); 

if( indscan < 7 ) % if less than 6 scan selected, let us add week 6 Friday scan 

data.xdata = [scheme3pts(1:indscan,1); 54];
data.ydata = [scheme3pts(1:indscan,2); fulldata.ydata(55)];

params = [growthParams beta_val(1)];
[timeFit, volFit] = tumorModel_compareSchemes('LOG+DVR',params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

    % yratio = lowfi(end)/lowfi(1);
    % if( yratio <= 0.2 );     nprior = 3;
    % elseif( yratio <= 0.85 );     nprior = 4;
    % else;     nprior = 5;
    % end
    % 
    % % P = readtable( './Update files 02-27-2025/Monolix and MATLAB - fit and distribution/prostate cancer model - large - Cov/populationParameters.txt');
    mu_est(3:5) = [  -3.3816e+00  -4.0682e+00  -5.4838e+00]; 
    sig_est(3:5) =  [ 1.960e-01   1.960e-01   1.960e-01 ]; 

params1 = {
    {'beta',beta_val(indscan), eps, 1, mu_est(nprior), sig_est(nprior) }
    };


model.ssfun = @(params,data) ssq_compareSchemes(params,data,IC,modelType,growthParams);

%%%% beta prior 
model.priorfun = @lognpdf;

options.updatesigma = 1;
options.waitbar = 0; 


options.nsimu = no_smps/5;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);



% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
beta_val_4 = chain(ind,:); %This is our fitted parameter set
beta_valMean_4 = mean(chain); %This is our mean value
params = [growthParams beta_val_4];


% Generate current model trajectory, credible intervals, and plot
[timeFit, volFit] = tumorModel_compareSchemes(modelType,params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

ln = fulldata.xdata(16:end);
modelfun = @(ln,params) tumorcred(ln,params,modelType,growthParams,IC);


figure(4); subplot( 2, 3, indscan+1 ); 
pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
pred.obslims = [];
mcmcpredplot(pred)
hold on
plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',1) %This is your optimal model fit
plot(data.xdata,data.ydata(:,1),'ok','MarkerFaceColor','k','Linewidth',1,'MarkerSize',8) %Selected scans
xlabel('Time (days)','FontSize',12);
ylabel('Volume','FontSize',12);
axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
set(gca, 'FontSize',12)
hold off
filename = [path strcat('ModelFit_Scheme4_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)


% Calculate MSE to measure error
mse_scheme4(indscan) = sum((lowfi-fulldata.ydata(16:end)).^2)/numel(lowfi);


disp( [mse_scheme3(indscan), mse_scheme4(indscan)] )
disp( mse_scheme3(indscan) > mse_scheme4(indscan) )


% Calculate uncertainty metrics
for n = 1:length(ln)
    widthCI(n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
end
areaCI_scheme4(indscan) = sum(widthCI(1:end));


scheme4pts= [data.xdata , data.ydata];

save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'mse_scheme4', 'areaCI_scheme4', 'scheme4pts', 'beta_val_4', 'indscan' ) 

else 
save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'indscan' ) 
    
end 

close(gcf)

end 



%% Fifth, calibrate with prior + 15 + 22 (assuming we scan on the second week) + chosen five (or less) points, until beta conv. 
if plotind(5) 

load( strcat( '/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Figures/LognormalPrior_250621_start22/Results_nPatient', int2str(npatient), '.mat') )


tmp = abs( beta_val(1:end-1) - beta_val(2:end) ); tmp = [1, tmp]; 
thresh_beta = 0.003; 
indscan = find( tmp < thresh_beta ) ;  
indscan = min( min(indscan), length(scheme3pts) ) ;

mse_scheme4 = mse_scheme3(1:indscan); 
areaCI_scheme4 = areaCI_scheme3(1:indscan); 

if( indscan < 7 ) % if less than 6 scan selected, let us add week 6 Friday scan 

data.xdata = [scheme3pts(1:indscan,1); 54];
data.ydata = [scheme3pts(1:indscan,2); fulldata.ydata(55)];

params = [growthParams beta_val(1)];
[timeFit, volFit] = tumorModel_compareSchemes('LOG+DVR',params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

    % yratio = lowfi(end)/lowfi(1);
    % if( yratio <= 0.2 );     nprior = 3;
    % elseif( yratio <= 0.85 );     nprior = 4;
    % else;     nprior = 5;
    % end
    % 
    % % P = readtable( './Update files 02-27-2025/Monolix and MATLAB - fit and distribution/prostate cancer model - large - Cov/populationParameters.txt');
    mu_est(3:5) = [  -3.3816e+00  -4.0682e+00  -5.4838e+00]; 
    sig_est(3:5) =  [ 1.960e-01   1.960e-01   1.960e-01 ]; 

params1 = {
    {'beta',beta_val(indscan), eps, 1, mu_est(nprior), sig_est(nprior) }
    };


model.ssfun = @(params,data) ssq_compareSchemes(params,data,IC,modelType,growthParams);

%%%% beta prior 
model.priorfun = @lognpdf;

options.updatesigma = 1;
options.waitbar = 0; 


options.nsimu = no_smps/5;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);



% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
beta_val_4 = chain(ind,:); %This is our fitted parameter set
beta_valMean_4 = mean(chain); %This is our mean value
params = [growthParams beta_val_4];


% Generate current model trajectory, credible intervals, and plot
[timeFit, volFit] = tumorModel_compareSchemes(modelType,params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

ln = fulldata.xdata(16:end);
modelfun = @(ln,params) tumorcred(ln,params,modelType,growthParams,IC);


figure(4); subplot( 2, 3, indscan+1 ); 
pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
pred.obslims = [];
mcmcpredplot(pred)
hold on
plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',1) %This is your optimal model fit
plot(data.xdata,data.ydata(:,1),'ok','MarkerFaceColor','k','Linewidth',1,'MarkerSize',8) %Selected scans
xlabel('Time (days)','FontSize',12);
ylabel('Volume','FontSize',12);
axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
set(gca, 'FontSize',12)
hold off
filename = [path strcat('ModelFit_Scheme4_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)
% close(gcf)


% Calculate MSE to measure error
mse_scheme4(indscan) = sum((lowfi-fulldata.ydata(16:end)).^2)/numel(lowfi);


disp( [mse_scheme3(indscan), mse_scheme4(indscan)] )
disp( mse_scheme3(indscan) > mse_scheme4(indscan) )




% Calculate uncertainty metrics
for n = 1:length(ln)
    widthCI(n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
end
areaCI_scheme4(indscan) = sum(widthCI(1:end));


scheme4pts= [data.xdata , data.ydata];

save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'mse_scheme4', 'areaCI_scheme4', 'scheme4pts', 'beta_val_4', 'indscan' ) 

else 
save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'indscan' ) 
    
end 

close(gcf)

end 



%% Sixth, calibrate with prior + 15 + 29 + modes for each group. 
if plotind(6) 

clear mse_scheme5 areaCI_scheme5 

beta_val = 0.0197; 

yratio = fulldata.ydata(end)/fulldata.ydata(16); 
    if( yratio <= 0.2 );     nprior = 3;
        scheme5pts =  [15  29    30 ];
    elseif( yratio <= 0.85 );     nprior = 4;
        scheme5pts =  [15  29    36 ];
    else;     nprior = 5;
        scheme5pts =  [15  29    33    50 ];
    end

    
    % lognormal three class of priors
    P = readtable( '../data/Update files 02-27-2025/Monolix and MATLAB - fit and distribution/prostate cancer model - large - Cov/populationParameters.txt');
    pop_est   = P{1:3,2}; %population level estimates for A, B, gamma
    class_est = P{4:5,2}; %this is used to determine the population gamma for the other groups
    sd_est    = P{6:8,2}; %standard deviation of random effect
    clear P
    
    growthParams = pop_est(1:2)';
    
    %calculate the population estimates of gamma for the other groups
    pop_est(4) = pop_est(3)*exp(class_est(1));
    pop_est(5) = pop_est(3)*exp(class_est(2));
    sd_est(4) = sd_est(3); sd_est(5) = sd_est(3);    
                                    
    mu_est(1:2) = log( pop_est(1:2) ); sig_est(1:2) = sd_est(1:2); 
    betascale = 0.0211; 
    
    mu_est(3:5) = log( pop_est(3:5)' * betascale ); 
    sig_est(3:5) = sd_est(3:5); 

    % mu_est(3:5) = [  -3.3816e+00  -4.0682e+00  -5.4838e+00]; 
    % sig_est(3:5) =  [ 1.960e-01   1.960e-01   1.960e-01 ]; 

    growthParams = pop_est(1:2)';    


% tmp = abs( beta_val(1:end-1) - beta_val(2:end) ); tmp = [1, tmp]; 
% thresh_beta = 0.003; 
% indscan = find( tmp < thresh_beta ) ;  
% indscan = min( min(indscan), length(scheme3pts) ) ;

% mse_scheme5 = mse_scheme3(1:indscan); 
% areaCI_scheme5 = areaCI_scheme3(1:indscan); 

if nprior ~= 5 



for indscan = 2:length(scheme5pts)

data.xdata = [scheme5pts(1:indscan)];
data.ydata = [fulldata.ydata([scheme5pts(1:indscan)]+1)]';

% params = [growthParams beta_val(1)];
% [timeFit, volFit] = tumorModel_compareSchemes('LOG+DVR',params,IC);
% lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

params1 = {
    {'beta',beta_val(indscan-1), eps, 1, mu_est(nprior), sig_est(nprior) }
    };


model.ssfun = @(params,data) ssq_compareSchemes(params,data,IC,modelType,growthParams);

%%%% beta prior 
model.priorfun = @(th,mu,sig)-2*log(prod(lognpdf(th,mu,sig)));

options.updatesigma = 1;
options.waitbar = 0; 


options.nsimu = no_smps/5;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);



% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
beta_val_5 = chain(ind,:); %This is our fitted parameter set
beta_valMean_5 = mean(chain); %This is our mean value
params = [growthParams beta_val_5];
beta_val = [beta_val beta_val_5(end)]; 


% Generate current model trajectory, credible intervals, and plot
[timeFit, volFit] = tumorModel_compareSchemes(modelType,params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

ln = fulldata.xdata(16:end);
modelfun = @(ln,params) tumorcred(ln,params,modelType,growthParams,IC);


figure(4); subplot( 2, 3, indscan-1 ); 
pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
pred.obslims = [];
mcmcpredplot(pred)
hold on
plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',1) %This is your optimal model fit
plot(data.xdata,data.ydata,'ok','MarkerFaceColor','k','Linewidth',1,'MarkerSize',8) %Selected scans
xlabel('Time (days)','FontSize',12);
ylabel('Volume','FontSize',12);
axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
set(gca, 'FontSize',12)
hold off
filename = [path strcat('ModelFit_Scheme5_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)


% Calculate MSE to measure error
mse_scheme5(indscan) = sum((lowfi-fulldata.ydata(16:end)).^2)/numel(fulldata.ydata(16:end));


% disp( [mse_scheme3(indscan), mse_scheme4(indscan)] )
% disp( mse_scheme3(indscan) > mse_scheme4(indscan) )


% Calculate uncertainty metrics
for n = 1:length(ln)
    widthCI(n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
end
areaCI_scheme5(indscan) = sum(widthCI(1:end));

   
end 

else 


growthParams = pop_est(1:2)';

% params = [growthParams beta_val(1)];
% [timeFit, volFit] = tumorModel_compareSchemes('LOG+DVR',params,IC);
% lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));


for indscan = 2:length(scheme5pts)

data.xdata = [scheme5pts(1:indscan)];
data.ydata = [fulldata.ydata([scheme5pts(1:indscan)]+1)]';

params1 = {
    {'K',growthParams(2), eps, 1, mu_est(2), sig_est(2) }    
    {'beta',beta_val(indscan-1), eps, 1, mu_est(nprior), sig_est(nprior) }
    };


model.ssfun = @(params,data) ssq_compareSchemes(params,data,IC,modelType,growthParams(1));
options.updatesigma = 1;

        %%%%%%%%%need to change this to new prior
        % prior
        lognormal = @(th,mu,sig) -2*log(prod( exp( -( (log(th)-mu)./(2*sig) ).^2 ) ./(th.*sig*sqrt(2*pi)) ));
        model.priorfun = lognormal; %@lognpdf;

        options.waitbar = 0;
        options.nsimu = no_smps/5;
        [results,chain,s2chain] = mcmcrun(model,data,params1,options);

        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);


        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        param_fit = chain(ind,:); %This is our fitted parameter set 
        params = [growthParams(1) param_fit];
        beta_val = [beta_val param_fit(end)]; 


% Generate current model trajectory, credible intervals, and plot
[timeFit, volFit] = tumorModel_compareSchemes(modelType,params,IC);
lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));        


ln = fulldata.xdata(16:end);
modelfun = @(ln,params) tumorcred(ln,params,modelType,growthParams(1),IC);        
        % modelfun2 = @(ln,params) tumorfun2(ln,params,icdata);


figure(6); subplot( 2, 3, indscan-1 ); 
pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
pred.obslims = [];
mcmcpredplot(pred)
hold on
plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',1) %This is your optimal model fit
plot(data.xdata,data.ydata,'ok','MarkerFaceColor','k','Linewidth',1,'MarkerSize',8) %Selected scans
xlabel('Time (days)','FontSize',12);
ylabel('Volume','FontSize',12);
axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
set(gca, 'FontSize',12)
hold off
filename = [path strcat('ModelFit_Scheme5_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)


% Calculate MSE to measure error
mse_scheme5(indscan) = sum((lowfi-fulldata.ydata(16:end)).^2)/numel(fulldata.ydata(16:end));


% disp( [mse_scheme3(indscan), mse_scheme4(indscan)] )
% disp( mse_scheme3(indscan) > mse_scheme4(indscan) )


% Calculate uncertainty metrics
for n = 1:length(ln)
    widthCI(n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
end
areaCI_scheme5(indscan) = sum(widthCI(1:end));


end

end

close(gcf)

scheme5pts= [data.xdata' , data.ydata'];

save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'mse_scheme5', 'areaCI_scheme5', 'scheme5pts', 'beta_val_5', 'indscan' ) 

% save( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')], '-append', 'indscan' ) 
 

 
end 



%% plot comparison figures 

% figure(11); subplot( 1, 2, 1); plot( mse_scheme1, 'x' ); hold on; plot( mean(mse_scheme2_all(:,17:20),2), '+' ); plot( mse_scheme3(1:(length(scheme3pts)-1)), 'o' ) 
% xlim( [.5 6.5] ); xlabel( 'Number of scans' ); ylabel( 'error' ) 
% set(gca, 'FontSize',14); grid on; 
% 
% subplot( 1, 2, 2);  plot( scheme1pts(4:end,1), mse_scheme1, 'x' ); hold on; plot( scheme2pts(2:7), mean(mse_scheme2_all(:,17:20),2), '+' ); plot( scheme3pts(2:end,1), mse_scheme3(1:(length(scheme3pts)-1)), 'o' )  
% xlabel( 'Time (days)' ); ylabel( 'error' ) 
% set(gca, 'FontSize',14); grid on; 
% legend( '[1,8,15]+equidistanced', '[15]+equidistanced', '[15]+selected')
% 
% filename = [path strcat('Error_Schemes_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% 
% figure(12); subplot( 1, 2, 1); plot( areaCI_scheme1, 'x' ); hold on; plot( mean(areaCI_scheme2_all(:,17:20),2), '+' ); plot( areaCI_scheme3(1:(length(scheme3pts)-1)), 'o' ) 
% xlim( [.5 6.5] ); xlabel( 'Number of scans' ); ylabel( 'uncertainty' ) 
% set(gca, 'FontSize',14); grid on; 
% 
% subplot( 1, 2, 2);  plot( scheme1pts(4:end,1), areaCI_scheme1, 'x' ); hold on; plot( scheme2pts(2:7), mean(areaCI_scheme2_all(:,17:20),2), '+' ); plot( scheme3pts(2:end,1), areaCI_scheme3(1:(length(scheme3pts)-1)), 'o' )  
% xlabel( 'Time (days)' ); ylabel( 'uncertainty' ) 
% set(gca, 'FontSize',14); grid on; 
% legend( '[1,8,15]+equidistanced', '[15]+equidistanced', '[15]+selected')
% 
% filename = [path strcat('Uncertainty_Schemes_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% 
% close all
% % clear all



% figure(11); subplot( 1, 2, 1); plot( mse_scheme1, 'x' ); hold on; plot( mse_scheme3(1:(length(scheme3pts)-1)), 'o' ); 
% if( indscan < 7 );  plot( mse_scheme4(1:(length(scheme4pts)-1)), '+' ); end 
% plot( mean(mse_scheme2_all(:,17:20),2), '*k' ); 
% xlim( [.5 6.5] ); xlabel( 'Number of scans' ); ylabel( 'error' ) 
% set(gca, 'FontSize',14); grid on; 
% 
% subplot( 1, 2, 2);  plot( scheme1pts(4:end,1), mse_scheme1, 'x' ); hold on; plot( scheme3pts(2:end,1), mse_scheme3(1:(length(scheme3pts)-1)), 'o' );  
% if( indscan < 7 );  plot( scheme4pts(2:end,1), mse_scheme4(1:(length(scheme4pts)-1)), '+' ); end 
% for mmm = 17:20; plot( mmm:7:56, (mse_scheme2_all(:,mmm)), '.:k' ); end
% xlabel( 'Time (days)' ); ylabel( 'error' ) 
% set(gca, 'FontSize',14); grid on; 
% if( indscan < 7 ); legend( '[1,8,15]+equidistanced', '[15]+selected 1', '[15]+selected 2', '[15]+equidistanced'); 
% else; legend( '[1,8,15]+equidistanced', '[15]+selected 1', '[15]+equidistanced'); end 
% 
% filename = [path strcat('Error_Schemes_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)
% 
% figure(12); subplot( 1, 2, 1); plot( areaCI_scheme1, 'x' ); hold on; plot( areaCI_scheme3(1:(length(scheme3pts)-1)), 'o' );  
% if( indscan < 7 );  plot( areaCI_scheme4(1:(length(scheme4pts)-1)), '+' ); end 
% plot( mean(areaCI_scheme2_all(:,17:20),2), '*k' ); 
% xlim( [.5 6.5] ); xlabel( 'Number of scans' ); ylabel( 'uncertainty' ) 
% set(gca, 'FontSize',14); grid on; 
% 
% subplot( 1, 2, 2);  plot( scheme1pts(4:end,1), areaCI_scheme1, 'x' ); hold on; plot( scheme3pts(2:end,1), areaCI_scheme3(1:(length(scheme3pts)-1)), 'o' ); 
% if( indscan < 7 );  plot( scheme4pts(2:end,1), areaCI_scheme4(1:(length(scheme4pts)-1)), '+' );  end 
% plot( scheme2pts(2:7), mean(areaCI_scheme2_all(:,17:20),2), '*k' );
% 
% xlabel( 'Time (days)' ); ylabel( 'uncertainty' ) 
% set(gca, 'FontSize',14); grid on; 
% if( indscan < 7 ); legend( '[1,8,15]+equidistanced', '[15]+selected 1', '[15]+selected 2', '[15]+equidistanced'); 
% else; legend( '[1,8,15]+equidistanced', '[15]+selected 1', '[15]+equidistanced'); end 
% 
% filename = [path strcat('Uncertainty_Schemes_np', int2str(npatient), '.jpg')];
% saveas(gcf,filename)



%% average error 
% npTot = 95; 
% npTot = 135; 
% err1(1:54,1:npTot) = NaN; 
% for mmm = 17:20; scheme2pts_all(:,mmm) = [ mmm:7:56 ]; end
% err2_all(1:54,1:4) = NaN;
% err2(1:54,1:npTot) = NaN; 
% err3(1:54,1:npTot) = NaN; 
% err4(1:54,1:npTot) = NaN; 
% 
% uq1(1:54,1:npTot) = NaN; 
% uq2_all(1:54,1:4) = NaN;
% uq2(1:54,1:npTot) = NaN; 
% uq3(1:54,1:npTot) = NaN; 
% uq4(1:54,1:npTot) = NaN; 
% 
% 
% 
% for np = 1:npTot 
% % npatient = indnp(np); 
% npatient = np; 
% 
% 
% load( [path strcat('Results_Schemes_np', int2str(npatient), '.mat')] )
% 
% for n = 1:6; err1( scheme1pts(n+3):54, np ) = mse_scheme1(n); end 
% 
% for mmm = 17:20; for n = 1:6 
%     err2_all( scheme2pts_all(n,mmm):54, mmm-16) = mse_scheme2_all(n,mmm); 
% end; end
% for n = 1:54 
%     ind = find( ~isnan(err2_all(n,:)) ); 
%     err2( n, np ) = mean( err2_all(n,ind) ); 
% end 
% 
% for n = 1:(length(scheme3pts)-1)
% err3( scheme3pts(n+1):54, np ) = mse_scheme3(n); 
% end 
% 
% for n = 1:(length(scheme4pts)-2) 
% err4( scheme4pts(n+1):54, np ) = mse_scheme4(n); 
% end 
% 
% err1table(:,np) = mse_scheme1; 
% err2table(:,np) = mse_scheme2; 
% err4table(:,np) = NaN; 
% err4table(1:(length(scheme4pts)-2),np) = mse_scheme4(1:(length(scheme4pts)-2)); 
% 
% 
% for n = 1:6; uq1( scheme1pts(n+3):54, np ) = areaCI_scheme1(n); end 
% 
% for mmm = 17:20; for n = 1:6 
%     uq2_all( scheme2pts_all(n,mmm):54, mmm-16) = areaCI_scheme2_all(n,mmm); 
% end; end
% for n = 1:54 
%     ind = find( ~isnan(uq2_all(n,:)) ); 
%     uq2( n, np ) = mean( uq2_all(n,ind) ); 
% end 
% 
% for n = 1:(length(scheme3pts)-1)
% uq3( scheme3pts(n+1):54, np ) = areaCI_scheme3(n); 
% end 
% 
% for n = 1:(length(scheme4pts)-2) 
% uq4( scheme4pts(n+1):54, np ) = areaCI_scheme4(n); 
% end 
% 
% end 
% 
% figure; hold on; 
% plot( mean(err1,2) ); plot( mean(err2,2) ); plot( mean(err3,2) ); plot( mean(err4,2), '--' ); 
% 
% figure; hold on; plot( mean(err1(:,1:30),2) ); plot( mean(err2(:,1:30),2) ); plot( mean(err3(:,1:30),2) ); plot( mean(err4(:,1:30),2), '--' );
% figure; hold on; plot( mean(err1(:,31:60),2) ); plot( mean(err2(:,31:60),2) ); plot( mean(err3(:,31:60),2) ); plot( mean(err4(:,31:60),2), '--' );
% figure; hold on; plot( mean(err1(:,61:95),2) ); plot( mean(err2(:,61:95),2) ); plot( mean(err3(:,61:95),2) ); plot( mean(err4(:,61:95),2), '--' );
% 
% 
% figure; hold on; 
% plot( mean(uq1,2) ); plot( mean(uq2,2) ); plot( mean(uq3,2) ); plot( mean(uq4,2), '--' ); 
% 
% figure; hold on; plot( mean(uq1(:,1:30),2) ); plot( mean(uq2(:,1:30),2) ); plot( mean(uq3(:,1:30),2) ); plot( mean(uq4(:,1:30),2), '--' );
% figure; hold on; plot( mean(uq1(:,31:60),2) ); plot( mean(uq2(:,31:60),2) ); plot( mean(uq3(:,31:60),2) ); plot( mean(uq4(:,31:60),2), '--' );
% figure; hold on; plot( mean(uq1(:,61:95),2) ); plot( mean(uq2(:,61:95),2) ); plot( mean(uq3(:,61:95),2) ); plot( mean(uq4(:,61:95),2), '--' );



% ???? 
% tmp = abs( beta_val(1:end-1) - beta_val(2:end) ); tmp = [1, tmp]; 
% ind = find( tmp < 0.003 ) ;  
% ind = min( min(ind), length(scheme3pts) ) ; 
% if( isempty(ind) ); ind = length(scheme3pts); end 
% 
% plot( scheme3pts(2:(ind),1), mse_scheme3(1:(ind-1)), 'o' ); 
% plot( 20:7:56, mean(mse_scheme2_all(:,17:20),2), '+' );  




close all
% clear all

%% Helper functions


%SSQ function for calibration of RT parameter
function SSrt = ssq_compareSchemes(params, data, IC, modelType, growthParams)

% [time,vol] = tumorModel_compareSchemes(modelType,[growthParams params],IC);
[time,vol] = tumorModel_postRT_wIC([growthParams params],IC);

tumVol = interp1(time,vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end

function SSrt = ssq_tumorVolpostRTnK(params,data)
global currentPts modelType growthParams

[time,vol] = tumorModel_postRT_wIC([growthParams(1) params],[data.xdata(1) data.ydata(1)] );

tumVol = interp1(time,vol, currentPts(:,1));

SSrt = sum((tumVol - currentPts(:,2)).^2);

end


% Function for credible interval plotting
function v=tumorcred(timef,params,modelType,growthParams,IC)

[tsol, ysol] = tumorModel_postRT_wIC([growthParams params'],IC);

v = interp1(tsol,ysol,timef);

end

function v=tumorfun2(timef,params,icdata)

global modelType growthParams


[tsol, ysol] = tumorModel_postRT_wIC([growthParams(1) params(1) params(2)],icdata);

v = interp1(tsol,ysol,timef);

end

function dval = tumorLogistic( t, val, params )

lambda = params(1);
K = params(2);

V = val(1);

dval(1) = lambda * V * ( 1 - V / K );

end 

function SSq = ssq_log(params, data)

[t,y] = ode23(@(t,y)tumorLogistic(t,y,params),[data.xdata(1), data.xdata(end)],data.ydata(1));

tumVol = interp1(t,y, data.xdata);

SSq = sum((tumVol - data.ydata').^2);

end


