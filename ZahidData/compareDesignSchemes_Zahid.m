%Compare errors, uncertainties, and fits for design schemes:

% - prior with -x, 0 & all data  

% - prior with -x, 0 & selected points 

% clear all
% close all

function compareDesignSchemes_Zahid


%% load data 
load('Data_Zahid.mat')
load('patientID.mat')

path = '../Figures/ZahidData_results/'; 
mkdir(path);

RGB = get(groot,"FactoryAxesColorOrder");

nplot = 0; 
plotind = [1 1 1 1]; 

no_smps = 10000; 


global growthParams V0 
growthParams = 0.33;        % new fit all together     %0.13;  % from Zahid's paper 

for npatient = [28, 37, 39] %1:39 

    disp(npatient) 

fulldata.xdata = data.xdata{npatient}; 
fulldata.ydata = data.ydata{npatient}; 

V0 = data.ydata{npatient}( abs(data.xdata{npatient}) < 0.1 ); 


%% fit PSI with two points 
fitdata.xdata = fulldata.xdata(1:2); 
fitdata.ydata = fulldata.ydata(1:2); 

lb = [0.1 ];
ub = [1.0 ];
opt = optimset('Display','off');
[param_val,ss01] = fmincon(@(param)ssq_tumorVol(param,fitdata),0.9,[],[],[],[],lb,ub, [], opt);
PSI_fit = param_val; 


K = V0 / PSI_fit; 
[time,vol] = ode23(@(t,y)tumorLogistic(t,y,[growthParams(1), K]), [fitdata.xdata(1),fitdata.xdata(end)], fitdata.ydata(1) );


lb_delta = 0.0; ub_delta = 0.15; 
lb_PSI = 0.1;   ub_PSI = 1; 


if ( ~ismember( npatient, [2, 3, 4, 5, 6, 7, 10, 24, 35] ) ) 

%% First, calibrate with initial -x, 0 + six scheme_pts - equi-spaced, only delta 

growthParams = [growthParams(1) V0/PSI_fit]; 

for n = 1:min(length(fulldata.ydata)-2, 6)

    fitdata.xdata = fulldata.xdata(1:(2+n)); 
    fitdata.ydata = fulldata.ydata(1:(2+n)); 

    [delta_fit,ss01] = fmincon(@(param)ssq_compareSchemes(param,fitdata),[0.1],[],[],[],[],lb_delta,ub_delta, [], opt);
    params_full = [growthParams delta_fit];

    [timeFit, volFit] = headnecktumor_RT(params_full,fulldata);
    lowfi = interp1(timeFit,volFit, fulldata.xdata);

    % Calculate MSE to measure error
    mse = sum((lowfi-fulldata.ydata).^2)/sum(fulldata.ydata.^2);
    lasterr1(npatient, n) = abs((lowfi(end)-fulldata.ydata(end)))/abs(fulldata.ydata(end));
    err_vec1(npatient, n) = [mse];

    if(nplot)
    figure(npatient*100+1); subplot( 3, 6, n ); hold on; 
    plot(  fulldata.xdata, fulldata.ydata, 'ko' ); 
    plot(  fitdata.xdata, fitdata.ydata, 'ko', 'MarkerFaceColor', 'k' ); 
    plot( timeFit, volFit, 'k-'); 
    xlabel( 'time (weeks)' ); ylabel('Tumor volume'); ylim( [0 max(fulldata.ydata)*1.1])
    end 

end 

else
%% Second, calibrate with initial -x, 0 + six scheme_pts - equi-spaced, PSI and delta 

growthParams = [growthParams(1)]; 

for n = 1:min(length(fulldata.ydata)-2, 6)

    fitdata.xdata = fulldata.xdata(1:(2+n)); 
    fitdata.ydata = fulldata.ydata(1:(2+n)); 


    [param_fit,ss01] = fmincon(@(param)ssq_compareSchemesPSI(param,fitdata),[0.9 0.1],[],[],[],[],[lb_PSI,lb_delta],[ub_PSI,ub_delta], [], opt);

    K = V0 / param_fit(1); 
    params_full = [growthParams K param_fit(2)];


    [timeFit, volFit] = headnecktumor_RT(params_full,fulldata);
    lowfi = interp1(timeFit,volFit, fulldata.xdata);

    % Calculate MSE to measure error
    mse = sum((lowfi-fulldata.ydata).^2)/sum(fulldata.ydata.^2);
    err_vec1(npatient, n) = [mse];
    lasterr1(npatient, n) = abs((lowfi(end)-fulldata.ydata(end)))/abs(fulldata.ydata(end));


    if(nplot)
    figure(npatient*100+1); subplot( 3, 6, n ); hold on; 
    plot(  fulldata.xdata, fulldata.ydata, 'ko' ); 
    plot(  fitdata.xdata, fitdata.ydata, 'ko', 'MarkerFaceColor', 'k' ); 
    plot( timeFit, volFit, 'k-'); 
    xlabel( 'time (weeks)' ); ylabel('Tumor volume'); ylim( [0 max(fulldata.ydata)*1.1])
    end

end 

end

%% Third, selected patient specific 

load( ['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Figures/251221_Zahid_onept/Zahid_Results_nPatient' int2str(npatient) '.mat'])
fulldata_zero.xdata = fulldata.xdata(2:end);
fulldata_zero.ydata = fulldata.ydata(2:end);

% check convergence 
    errcheck = abs( delta_vec(1:end-1) - delta_vec(2:end) ); 
    relerrcheck = abs( delta_vec(1:end-1) - delta_vec(2:end) )./delta_vec(2:end); 
    ind1 = find( errcheck < 0.005, 1 ); 
    ind2 = find( relerrcheck < 0.05, 1 ); 
    if( isempty(ind1) ); ind1=100; end
    if( isempty(ind2) ); ind2=100; end
    if( (ind1+ind2)==200 ); ind1 = length(delta_vec)-1; end
    ind = min( ind1, ind2 );
    nscan(npatient) = ind+1; 

    point_list = point_list(1:(2+ind),1:2); 


for n = 2:size( point_list, 1 )


if( isempty( PSI_vec ) )

    fitdata.xdata = point_list(1:n,1); 
    fitdata.ydata = point_list(1:n,2); 

    pop_PSI_est = 0.86; 
    growthParams = [growthParams(1) V0/pop_PSI_est]; 

    params_full = [growthParams delta_vec(n-1)];

    [timeFit, volFit] = headnecktumor_RT(params_full,fulldata_zero);
    lowfi = interp1(timeFit,volFit, fulldata_zero.xdata);

    % Calculate MSE to measure error
    mse = sum((lowfi-fulldata_zero.ydata).^2)/sum(fulldata_zero.ydata.^2);
    err_vec3(npatient, n) = [mse];
    lasterr3(npatient, n) = abs((lowfi(end)-fulldata_zero.ydata(end)))/abs(fulldata_zero.ydata(end));

    if(nplot)
    figure(npatient*100+1); subplot( 3, 6, 6+n-1 ); hold on; 
    plot(  fulldata.xdata, fulldata.ydata, 'ko' ); 
    plot(  fitdata.xdata, fitdata.ydata, 'ko', 'MarkerFaceColor', 'k' ); 
    plot( timeFit, volFit, 'k-'); 
    xlabel( 'time (weeks)' ); ylabel('Tumor volume'); ylim( [0 max(fulldata.ydata)*1.1])
    end 

else 

    fitdata.xdata = point_list(1:n,1); 
    fitdata.ydata = point_list(1:n,2); 

    
    growthParams = [growthParams(1)]; 
    params_full = [growthParams V0/PSI_vec(n-1) delta_vec(n-1)];

    [timeFit, volFit] = headnecktumor_RT(params_full,fulldata_zero);
    lowfi = interp1(timeFit,volFit, fulldata_zero.xdata);

    % Calculate MSE to measure error
    mse = sum((lowfi-fulldata_zero.ydata).^2)/sum(fulldata_zero.ydata.^2);
    err_vec3(npatient, n) = [mse];
    lasterr3(npatient, n) = abs((lowfi(end)-fulldata_zero.ydata(end)))/abs(fulldata_zero.ydata(end));

    if(nplot)
    figure(npatient*100+1); subplot( 3, 6, 6+n-1 ); hold on; 
    plot(  fulldata.xdata, fulldata.ydata, 'ko' ); 
    plot(  fitdata.xdata, fitdata.ydata, 'ko', 'MarkerFaceColor', 'k' ); 
    plot( timeFit, volFit, 'k-'); 
    xlabel( 'time (weeks)' ); ylabel('Tumor volume'); ylim( [0 max(fulldata.ydata)*1.1])
    end

end

end


% group schedule 
    lowfi = interp1(timeFit,volFit, [0, 6]);
    %%%% ratio to determine
    yratio = lowfi(end)/V0;    
    
    if( yratio <= 1/3 )
        nprior = 1;
        group_schedule = [0, 2, 3, 4]+1; 
    elseif( yratio <= 2/3 )
        nprior = 2;
        group_schedule = [0, 2, 3, 5]+1; 
    else
        nprior = 3;
        group_schedule = [0, 2, 3, 4]+1; 
    end 

    mu_est   = [-2.3212   -3.4585   -4.3293]; 
    sig_est = [ 0.1748    0.3936    0.4564]; 

    if( group_schedule(end) > length(fulldata_zero.xdata) )
        group_schedule = group_schedule(1:end-1); 
    end
    % if( npatient == 22 || npatient == 23 ); group_schedule = [0, 2, 3]+1; end 
 
    ngroup(npatient) = nprior; 
    group_schedule_all(npatient,1:length(group_schedule)) = group_schedule-1; 

fulldata_zero.xdata = fulldata.xdata(2:end);
fulldata_zero.ydata = fulldata.ydata(2:end);

for n = 2:length( group_schedule )

% if( nprior < 3 )
if( isempty( PSI_vec ) )

    fitdata.xdata = fulldata_zero.xdata( group_schedule(1:n) ); 
    fitdata.ydata = fulldata_zero.ydata( group_schedule(1:n) );  

    pop_PSI_est = 0.86; 
    growthParams = [growthParams(1) V0/pop_PSI_est]; 


    [delta_fit,ss01] = fmincon(@(param)ssq_compareSchemes(param,fitdata),delta_vec(min(n-1,length(delta_vec))),[],[],[],[],lb_delta,ub_delta, [], opt);
    params_full = [growthParams delta_fit];

        %%%%%%%%%need to change this to new prior
        params1 = {
            {'delta',delta_fit, eps, ub_delta, mu_est(nprior), sig_est(nprior) }
            };


        model.ssfun = @ssq_compareSchemes;
        % ss0 = ssq_compareSchemes(delta_fit,fitdata); 
        mse = ss01/(length(fitdata.ydata(:))-length(delta_fit));
        model.sigma2 = mse; %Initial guess for error variance

        options.updatesigma = 1;

        model.priorfun = @(th,mu,sig)-2*log(lognpdf(th,mu,sig)); 

        options.waitbar = 0;
        options.nsimu = no_smps/5;
        [results,chain,s2chain] = mcmcrun(model,fitdata,params1,options);

        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,fitdata,params1,options,results);


        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        delta_fit = chain(ind,:); %This is our fitted parameter set\
        params = [growthParams delta_fit];

        %Save off parameter estimates and chains
        delta_vec = [delta_vec; delta_fit];


        % Generate current model trajectory, credible intervals, and plot
        checkdata = fitdata; if( checkdata.xdata(end) < 6 ); checkdata.xdata = [checkdata.xdata, max(6,fulldata.xdata(end))]; end 

        [timeFit, volFit] = headnecktumor_RT(params,checkdata);
        lowfi = interp1(timeFit,volFit, fulldata.xdata);



        ln = fulldata_zero.xdata';
        modelfun = @(ln,params) tumorfun(ln,params,checkdata);

        figure(npatient); subplot( 1, 5, n-1 ); hold on; 
        pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
        % pred.obslims = [];
        mcmcpredplot(pred)
        hold on
        h=gca;
        plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
        plot(fulldata.xdata,lowfi,'-k','Linewidth',2) %This is your optimal model fit
        plot(fitdata.xdata,fitdata.ydata,'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
        xlabel('Time (days)','FontSize',18);
        ylabel('Volume','FontSize',18);
        titleName = ['Patient #' int2str(npatient) ': Iteration ' num2str(n)];
        title(string(titleName),'Interpreter','none','FontSize',14)
        axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
        set(h, 'FontSize',12)
        hold off



    [timeFit, volFit] = headnecktumor_RT(params_full,fulldata_zero);
    lowfi = interp1(timeFit,volFit, fulldata_zero.xdata);

    % Calculate MSE to measure error
    mse = sum((lowfi-fulldata_zero.ydata).^2)/sum(fulldata_zero.ydata.^2);
    err_vec4(npatient, n) = [mse];
    lasterr4(npatient, n) = abs((lowfi(end)-fulldata_zero.ydata(end)))/abs(fulldata_zero.ydata(end));

    if(nplot)
    figure(npatient*100+1); subplot( 3, 6, 12+n-1 ); hold on; 
    plot(  fulldata.xdata, fulldata.ydata, 'ko' ); 
    plot(  fitdata.xdata, fitdata.ydata, 'ko', 'MarkerFaceColor', 'k' ); 
    plot( timeFit, volFit, 'k-'); 
    xlabel( 'time (weeks)' ); ylabel('Tumor volume'); ylim( [0 max(fulldata.ydata)*1.1])
    end 

else 

    fitdata.xdata = fulldata_zero.xdata( group_schedule(1:n) ); 
    fitdata.ydata = fulldata_zero.ydata( group_schedule(1:n) );  
    
    growthParams = [growthParams(1)]; 

    [param_fit,ss01] = fmincon(@(param)ssq_compareSchemesPSI(param,fitdata),[PSI_vec(min(n-1,length(PSI_vec))) delta_vec(min(n-1,length(delta_vec)))],[],[],[],[],[lb_PSI,lb_delta],[ub_PSI,ub_delta], [], opt);

    K = V0 / param_fit(1); 
    params_full = [growthParams K param_fit(2)];

    [timeFit, volFit] = headnecktumor_RT(params_full,fulldata_zero);
    lowfi = interp1(timeFit,volFit, fulldata_zero.xdata);

    % Calculate MSE to measure error
    mse = sum((lowfi-fulldata_zero.ydata).^2)/sum(fulldata_zero.ydata.^2);
    err_vec4(npatient, n) = [mse];
    lasterr4(npatient, n) = abs((lowfi(end)-fulldata_zero.ydata(end)))/abs(fulldata_zero.ydata(end));

    if(nplot)
    figure(npatient*100+1); subplot( 3, 6, 12+n-1 ); hold on; 
    plot(  fulldata.xdata, fulldata.ydata, 'ko' ); 
    plot(  fitdata.xdata, fitdata.ydata, 'ko', 'MarkerFaceColor', 'k' ); 
    plot( timeFit, volFit, 'k-'); 
    xlabel( 'time (weeks)' ); ylabel('Tumor volume'); ylim( [0 max(fulldata.ydata)*1.1])
    end

end

end



if(nplot)
filename = [path strcat('ModelFit_Zahid_Schemeall_np', int2str(npatient), '.jpg')];
saveas(gcf,filename)
filename = [path strcat('ModelFit_Zahid_Schemeall_np', int2str(npatient), '.fig')];
saveas(gcf,filename)


RGB = get(groot,"FactoryAxesColorOrder");

if( length(fulldata.xdata) > 8 ); 
    fulldata.xdata = fulldata.xdata(1:8); 
end 

ii = length( fulldata.xdata(3:end) );
jj = length( point_list(2:end,1) )+1; 
figure(1111); subplot( 5, 8, npatient ); hold on; set(gca,'Yscale', 'log' ); box on; grid on; 
plot( err_vec1(npatient,1:ii), 'x:', 'Color', RGB(1,:))
plot( err_vec3(npatient,2:jj), '+:', 'Color', RGB(2,:))
plot( err_vec4(npatient,2:end), 'o:', 'Color', RGB(4,:))

figure(1212); subplot( 5, 8, npatient ); hold on; set(gca,'Yscale', 'log' ); box on; grid on;  
plot( fulldata.xdata(3:end), err_vec1(npatient,1:ii), 'x:', 'Color', RGB(1,:))
plot( point_list(2:end,1), err_vec3(npatient,2:jj), '+:', 'Color', RGB(2,:))
plot( group_schedule(2:end)-1, err_vec4(npatient,2:length(group_schedule)), 'o:', 'Color', RGB(4,:))

end


end 

figure(1111);
for np = 1:8:39; subplot( 5, 8, np ); ylabel( 'error' ); end
for np = 33:39; subplot( 5, 8, np ); xlabel( 'scan' ); end
for np = 1:39; subplot( 5, 8, np ); title( patientID{np} ); end
figure(1212); 
for np = 1:8:39; subplot( 5, 8, np ); ylabel( 'error' ); end
for np = 33:39; subplot( 5, 8, np ); xlabel( 'time (weeks)' ); end
for np = 1:39; subplot( 5, 8, np ); title( patientID{np} ); end

err_vec4 = err_vec4(:,2:4); 

for nn = 1:3 
ii = find( ngroup == nn ); 
figure(12515); subplot( 1, 3, nn); hold on; plot( mean( err_vec1(ii,:) ), ':x', 'color', RGB(1,:) ); 
figure(12516); subplot( 1, 3, nn); hold on; plot( mean( err_vec1(ii,:) ), ':x', 'color', RGB(1,:) ); 
figure(12515); subplot( 1, 3, nn); hold on; plot( mean( err_vec4(ii,:) ), ':o', 'color', RGB(4,:) );
    if( ngroup == 1 )
        group_schedule = [2, 3, 4]; 
    elseif( ngroup == 2 ) 
        group_schedule = [2, 3, 5]; 
    else
        group_schedule = [2, 3, 4]; 
    end 
figure(12516); subplot( 1, 3, nn); hold on; plot(group_schedule, mean( err_vec4(ii,:) ), ':o', 'color', RGB(4,:) );
end

figure(12515); subplot( 1, 3, 1); ylabel( 'error' );
for nn = 1:3; subplot( 1, 3, nn ); xlabel( 'scan' ); set(gca,'fontsize', 12 ); box on; grid on; xlim( [1 6] ); end 
figure(12516); subplot( 1, 3, 1); ylabel( 'error' );
for nn = 1:3; subplot( 1, 3, nn ); xlabel( 'time (weeks)' ); set(gca,'fontsize', 12 ); box on; grid on; xlim( [1 6] ); end 



for nn = 1:3 
ii = find( ngroup == nn ); 

x =    1:6;  
msemean = mean( err_vec1(ngroup == nn,:) ); 
msestd = std( err_vec1(ngroup == nn,:) ); 

tmp = min(err_vec1(ngroup == nn,:)); tmp(tmp==0) = min( tmp(tmp~=0) ) ; 
figure(12517); subplot( 1, 3, nn ); hold on; 
plot( x, msemean, 'color', RGB(1,:), 'linewidth', 2 ); 
patch([x fliplr(x)], [tmp fliplr(max(err_vec1(ngroup == nn,:)))], RGB(1,:) )
alpha( 0.5 );

% figure(12518); subplot( 1, 3, nn ); hold on; 
% plot( x, msemean, 'color', RGB(1,:), 'linewidth', 2 ); 
% patch([x fliplr(x)], [msemean-msestd fliplr(msemean+msestd)], RGB(1,:) )
% alpha( 0.5 );

    if( nn == 2 ) 
        x = [2, 3, 5]; 
    else 
        x = [2, 3, 4]; 
    end 

msemean = mean( err_vec4(ngroup == nn,:) ); 
msestd = std( err_vec4(ngroup == nn,:) ); 
tmp = min(err_vec4(ngroup == nn,:)); tmp(tmp==0) = min( tmp(tmp~=0) ) ; 

figure(12517); subplot( 1, 3, nn ); hold on; 
plot( x, msemean, 'color', RGB(4,:), 'linewidth', 2 ); 
patch([x fliplr(x)], [tmp fliplr(max(err_vec4(ngroup == nn,:)))], RGB(4,:) )
alpha( 0.5 );
grid on; box on; 
set(gca, 'Yscale', 'log' ); box on; xlabel('time (weeks)'); ylabel('Error'); grid on; %ylim( [0.01 30] ); 
set(gca,'fontsize', 12 ); title( strcat('Group ', int2str(nn)) ); ylim( [0.0006 1] )

% figure(12518); subplot( 1, 3, nn ); hold on; 
% plot( x, msemean, 'color', RGB(1,:), 'linewidth', 2 ); 
% patch([x fliplr(x)], [msemean-msestd fliplr(msemean+msestd)], RGB(4,:) )
% alpha( 0.5 );

end

figure; 
for nn = 1:3 
subplot( 1, 3, nn); hold on; 

pts = ones(sum(ngroup == nn),1) * [1:3];
err = err_vec4(ngroup == nn,:); 
boxplot( err(:), pts(:)-0.2, 'colors', RGB(4,:), 'Positions', unique(pts(:)-0.2),'PlotStyle','compact' ); 
pts = ones(sum(ngroup == nn),1) * [1:6];
err = err_vec1(ngroup == nn,:); 
boxplot( err(:), pts(:), 'colors', RGB(1,:), 'Positions', unique(pts(:)),'PlotStyle','compact' ); 
set(gca, 'Yscale', 'log' ); box on; xlabel('Scan number'); ylabel('Error'); grid on; %ylim( [0.01 30] ); 
set(gca,'fontsize', 12 ); title( strcat('Group ', int2str(nn)) ); ylim( [0.0006 1] )
end 
legend({'Mean Schedule A',  'Mean Schedule C' } )


filename = [path 'ModelFit_Zahid_Schemeall_errors.mat'];

save( filename, 'err_vec1', 'err_vec3', 'err_vec4', 'point_list', 'ngroup', 'group_schedule_all')





end 

%% Helper functions


%SSQ function for calibration of RT parameter
function SSrt = ssq_compareSchemesPSI(params, data)
global growthParams V0 

K = V0 / params(1); 

[time,vol] = headnecktumor_RT( [growthParams(1) K params(end)], data);

tumVol = interp1(time, vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end

%SSQ function for calibration of RT parameter
function SSrt = ssq_compareSchemes(params, data)
global growthParams 
[time,vol] = headnecktumor_RT([growthParams params],data);

tumVol = interp1(time, vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end


%SSQ function for calibration of RT parameter
% function SSrt = ssq_compareSchemes_wIC(params, data)
% global growthParams 
% 
% data.ydata(1) = params(end); 
% 
% [time,vol] = headnecktumor_RT([growthParams params(1:2)],data);
% 
% tumVol = interp1(time, vol, data.xdata);
% 
% SSrt = sum((tumVol - data.ydata(:,1)).^2);
% 
% end

% Function for credible interval plotting
function v=tumorcred(timef,params,modelType,growthParams,IC)

data.xdata(1) = IC(1); data.xdata(2) = timef(end); 
data.ydata(1) = IC(2); 

[tsol, ysol] = headnecktumor_RT([growthParams params], data);

v = interp1(tsol,ysol,timef);

end


% SSQ function for calibration before RT 
function SSrt = ssq_tumorVol(params,data)
global currentPts growthParams V0 

K = V0 / params(1); 

[time,vol] = ode23(@(t,y)tumorLogistic(t,y,[growthParams, K]), [data.xdata(1),data.xdata(end)], data.ydata(1) ); 

tumVol = interp1(time,vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end


function dval = tumorLogistic( t, val, params )

lambda = params(1);
K = params(2);

V = val(1);

dval(1) = lambda * V * ( 1 - V / K );

end 

% Function for credible interval plotting
function v=tumorfun(timef,params,data)

global modelType growthParams


[tsol, ysol] = headnecktumor_RT([growthParams params(end)],data);

v = interp1(tsol,ysol,timef);

end


% function SSq = ssq_log(params, data)
% 
% [t,y] = ode23(@(t,y)tumorLogistic(t,y,params),[data.xdata(1), data.xdata(end)],data.ydata(1));
% 
% tumVol = interp1(t,y, data.xdata);
% 
% SSq = sum((tumVol - data.ydata').^2);
% 
% end
% 

