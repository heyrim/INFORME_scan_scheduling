
function point_list = seqDesign_ZahidData_prior_all_onedata(npatient)

global currentPts modelType growthParams M V0 

if( nargin < 1 ); npatient = 1; end 
disp( strcat( 'Start patient #', int2str(npatient) ) )

modelType = 'PSI+CCR';

nplot = 1; % figures ON or OFF 
npopfit = 1; % population fit ON or OFF 

RGB = get(groot,"FactoryAxesColorOrder");

%% load data
load('Data_Zahid.mat')

path = 'Figures/251202_Zahid/'; 
mkdir( path )


%Settings
no_smps = 10000; %DRAM chain lengths

%Set up vectors for saved quantities;
delta_vec = []; %store parameter estimates at each iteration
K_vec = []; 
err_vec = []; %store MSE at each iteration
widthCI = []; %store width of credible interval at each iteration
areaCI = []; %store area of credible interval at each iteration
rcv_vec = [];


fulldata.xdata = data.xdata{npatient}; 
fulldata.ydata = data.ydata{npatient}; 
fulldata.expDesign = 1:6; 

nclass = data.nclass(npatient);

% clear data

currentPts = [ fulldata.xdata(1:2)', fulldata.ydata(1:2)' ];
icdata = currentPts;


%% Instead of initial calibration to estimate params from pre-treatment data
%% use prior distribution estimated from group data 


% Tin's NLME result 
% growth rate – logitnormal distribution with range , population value: 0.33, SD: 8.28.
% PSI – logitnormal distribution with range , population value: 0.90, SD: 1.77.

growthParams = 0.33; % 

% Logitnormal for PSI 
pop_PSI_est = 0.86; %0.94; %0.98; %0.9; 
sd_PSI_est = 0.65; %1.47; %2.1; %1.77; 

% Logitnormal for delta 
% pop_delta_est = 0.092; 
% sd_delta_est = 1.25; 

% % Group population value: 0.054, SD: 1.71 
% % next group population value: 0.085, SD: 1.73
% 0.017, SD: 1.04
% deltascale = tmp\param_save(2,:)'; % scaling between Monolix and Matlab 
deltascale =    2.2933; 

pop_delta_est = 0.017 *deltascale; 
sd_delta_est = 1.04; 

lb_PSI = 0.1; ub_PSI = 1.0; 
lb_delta = 0; ub_delta = 0.15; 

mu_PSI_est = -0.1615;   sig_PSI_est = 0.0628; 


%% First fit PSI or K using first two points 

V0 = fulldata.ydata( abs(data.xdata{npatient}) < 0.1 ); 

fitdata.xdata = data.xdata{npatient}(1:2); 
fitdata.ydata = data.ydata{npatient}(1:2); 

lb = [0.1 ];
ub = [1.0 ];
opt = optimset('Display','off');
[param_val,ss01] = fmincon(@(param)ssq_tumorVol(param,fitdata),pop_PSI_est,[],[],[],[],lb,ub, [], opt);
PSI_fit = param_val; 

PSI_vec = PSI_fit; 


K = V0 / PSI_fit; 
[time,vol] = ode23(@(t,y)tumorLogistic(t,y,[growthParams, K]), [fitdata.xdata(1),fitdata.xdata(end)], fitdata.ydata(1) );

if( nplot ); figure(221); subplot( 5, 8, npatient); hold on; plot( time, vol ); plot( fulldata.xdata, fulldata.ydata, 'o' ); end 


point_list = 2; % precomputed point 

 

%% using optimal experimental data for the group, calibrate low-fidelity model and classify new group 
%Add point to current list of points for next round of calibration 
currentPts = [currentPts; fulldata.xdata(point_list + 2) fulldata.ydata(point_list + 2)]; 

fitdata.xdata = currentPts(:,1);
fitdata.ydata = currentPts(:,2);


growthParams = [growthParams(1) V0/PSI_fit]; 
lb = [0.0];
ub = [0.15];
opt = optimset('Display','off');
[param_fit,ss01] = fmincon(@(param)ssq_compareSchemes(param,fitdata),[pop_delta_est],[],[],[],[],lb,ub, [], opt);

delta_fit = param_fit; 

checkdata = fitdata; if( checkdata.xdata(end) < 6 ); checkdata.xdata = [checkdata.xdata; max(6,fulldata.xdata(end))]; end 
params = [growthParams(1:2) delta_fit]; 


[time, vol] = headnecktumor_RT( params , checkdata); 
lowfi = interp1(time, vol, [currentPts(:,1);6]);

if( nplot ); figure(221); %subplot( 5, 8, npatient); 
    hold on; plot( time, vol, ':', 'color', RGB(1,:) ); end 

% %%%% Change prior depending on the first data point
% checkdata = fitdata; if( checkdata.xdata(end) < 6 ); checkdata.xdata = [checkdata.xdata; 6]; end 
% [time, vol] = headnecktumor_RT( [growthParams(1) V0/param_fit(1) param_fit(2)], checkdata); 
% lowfi = interp1(time, vol, 6);


%%%% ratio to determine
yratio = lowfi(end)/V0;


if( yratio <= 1/3 )
    nprior = 1;
elseif( yratio <= 2/3 )
    nprior = 2;
else
    nprior = 3;
end

% %%%% prior 
% % % Logitnormal for delta try 1 
% % % range: [0.0,0.1] but increase to 0.1 and over  
% % % three groups 
% % % High responder: population value: 0.1, SD: 0.74	
% % % Medium responder: population value: 0.046, SD: 0.74	
% % % Low responder: population value: 0.019, SD: 0.74
% % % 
% % % Group population value: 0.054, SD: 1.71 
% 
% % pop_delta_est = [0.1, 0.046, 0.019]; 
% % sd_delta_est = 0.74; 
% 
% 
% % try 2 
% % high responder - population value: 0.15, SD: 0.64
% % medium responder - population value: 0.071, SD: 0.64
% % low responder - population value: 0.027, SD: 0.64
% 
% % pop_delta_est = [0.1, 0.046, 0.019]; 
% % sd_delta_est = 0.74; 
% 
% % try 3 
% % high responder - population value: 0.044, SD: 0.47
% % medium responder - population value: 0.014, SD: 0.47
% % low responder - population value: 0.0058, SD: 0.47
% % 0.044 0.014 0.0058
% 
% 
lb_delta = 0.0; ub_delta = 0.15; 

pop_delta_est = [0.044 0.014 0.0058]*deltascale;
sd_delta_est = 0.47*sqrt(deltascale); 

a = lb_delta; b = ub_delta; M = 100000; 
for n = 1:3; 
logit_gamma  = log((pop_delta_est(n) - a)/(b - (pop_delta_est(n))));
guassian_num = normrnd(logit_gamma,sd_delta_est(1)^2, [1,M]);
newChain1 = (b.*exp(guassian_num) + a)./(1+exp(guassian_num)); 
th = lognfit( newChain1 ) ; 
end

mu_est   = [-2.3212   -3.4585   -4.3293]; 
sig_est = [ 0.1748    0.3936    0.4564]; 
    


% fit prior to lognormals 


% Now start the sequential design procedure, estimating beta + other param if needed
% tmp = round( fulldata.xdata ); 
% expDesigns = tmp( tmp > point_list );

expDesigns = (point_list+1):length(fulldata.xdata(3:end));


nIter = 0; %keep track of how many calibration iterations have been performed
ptsLeft = 1; %binary tracker; 1 means there are still points left to choose

%% make the decision to estimate beta only or beta + K
 
fit_thresh = -1; %251202 
%0.05; % 251122 results with thresh  

last_fit_error = max( sum( abs( currentPts(:,2) - lowfi( 1:end-1 ) ) ) / sum( currentPts(:,2) ), ... 
    abs( currentPts(end,2) - lowfi(end-1) )/abs( currentPts(end,2) ) ); 


while ptsLeft == 1 && nIter < 5

    nIter = nIter+1;

    data.xdata = currentPts(:,1);
    data.ydata = currentPts(:,2);

    [delta_fit,ss01] = fmincon(@(param)ssq_compareSchemes(param,data),[delta_fit],[],[],[],[],lb_delta,ub_delta, [], opt);

    [timeFit, volFit] = headnecktumor_RT([growthParams delta_fit],checkdata);        
    lowfi = interp1(timeFit,volFit, currentPts(:,1));

    last_fit_error = max( sum( abs( currentPts(:,2) - lowfi  )./abs(currentPts(:,2)) ) / size(currentPts,1), ... 
        abs( (currentPts(end,2) - lowfi(end))/currentPts(end,2) ) ); 


    % if( last_fit_error < fit_thresh ) %%%% beta only
    if( nprior ~= 3 ) 
        disp( 'delta only' )

        [delta_fit,ss01] = fmincon(@(param)ssq_compareSchemes(param,data),[delta_fit],[],[],[],[],lb_delta,ub_delta, [], opt);
        params = [growthParams delta_fit]; 

        %%%%%%%%%need to change this to new prior
        params1 = {
            % {'lambda',pop_est(1), eps, 1, pop_est(1), sd_est(1) }
            % {'K',pop_est(2), eps, 10, pop_est(2), sd_est(2) }
            {'delta',delta_fit, eps, ub_delta, mu_est(nprior), sig_est(nprior) }
            % {'delta',delta_fit, eps, ub_delta }
            };


        model.ssfun = @ssq_compareSchemes;
        ss0 = ssq_compareSchemes(params,data); 
        mse = ss0/(length(currentPts(:,2))-length(delta_fit));
        model.sigma2 = mse; %Initial guess for error variance

        options.updatesigma = 1;

        %%%%%%%%%need to change this to new prior
        % prior
        % % lognormal = @(th,mu,sig) exp( -sum( ((log(th)-mu+log(120))./sig).^2 ) );
        % lognormal = @(th,mu,sig) exp( -( ((log(th)-mu+log(120))./sig).^2 )/2 ) ./th / sig / sqrt(2*pi); 
        % model.priorfun = lognormal;
        % model.priorfun = @lognpdf;
        model.priorfun = @(th,mu,sig)-2*log(lognpdf(th,mu,sig)); 


        options.waitbar = 0;
        options.nsimu = no_smps/5;
        [results,chain,s2chain] = mcmcrun(model,data,params1,options);

        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);


        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        delta_fit = chain(ind,:); %This is our fitted parameter set\
        params = [growthParams delta_fit];

        %Save off parameter estimates and chains
        delta_vec = [delta_vec; delta_fit];


        % Generate current model trajectory, credible intervals, and plot
        [timeFit, volFit] = headnecktumor_RT(params,checkdata);
        lowfi = interp1(timeFit,volFit, fulldata.xdata);



        ln = fulldata.xdata';
        modelfun = @(ln,params) tumorfun(ln,params,checkdata);

        figure(1)
        pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
        pred.obslims = [];
        mcmcpredplot(pred)
        hold on
        h=gca;
        plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
        plot(fulldata.xdata,lowfi,'-k','Linewidth',2) %This is your optimal model fit
        plot(currentPts(:,1),currentPts(:,2),'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
        xlabel('Time (days)','FontSize',18);
        ylabel('Volume','FontSize',18);
        titleName = ['Patient #' int2str(npatient) ': Iteration ' num2str(nIter)];
        title(string(titleName),'Interpreter','none','FontSize',14)
        axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
        set(h, 'FontSize',18)
        hold off
        filename = [path 'ModelFit_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        saveas(gcf,filename)
        close(gcf)


        % % figure(2)
        % % mcmcplot(chain,[],{'\beta'},'chainpanel')
        % % filename = [path 'Chain_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        % % saveas(gcf,filename)
        % % % filename = [path 'Chain_Iteration' num2str(nIter) '.fig'];
        % % % saveas(gcf,filename)
        % % close(gcf)

        figure(3); hold on; 
        yy = mcmcplot(chain,[],{'\beta'},'denspanel');
        filename = [path 'Density_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        saveas(gcf,filename)
        % % filename = [path 'Density_Iteration' num2str(nIter) '.fig'];
        % % saveas(gcf,filename)
        close(gcf)
        pdf_vec{nIter} = yy{1};


        % Calculate MSE to measure error
        mse = sum((lowfi-fulldata.ydata).^2)/sum(fulldata.ydata.^2);
        err_vec = [err_vec; mse];

        % Calculate uncertainty metrics
        for n = 1:length(ln)
            widthCI(nIter,n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
        end
        areaCI = [areaCI; sum(widthCI(nIter,1:end))];



        % Determine which point should be selected next
        if currentPts(end,1)==fulldata.xdata(end)
            ptsLeft = 0; %this is the final iteration

        else

            knnSet = 10:10:no_smps; %thin out samples for kNN analysis
            newChain = chain(knnSet,:);

            % Calculate predicted values at each experimental design with
            % each parameter set from newChain
            for ii = 1:size(newChain,1)

                [tsol, ysol] = headnecktumor_RT([growthParams newChain(ii,:)],checkdata);

                for jj = 1:length(expDesigns)
                    lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
                end
            end

            %Calculate MI for each remaining design

            %%%%% with normalized data
            %         renewChain = newChain/max(newChain);
            %         relowfiOut = lowfiOut/max(lowfiOut(:));

            %%%%% with normalized by mean/std data
            normChain = (newChain-mean(newChain))./std(newChain);
            normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);

            for jj = 1:length(expDesigns)
                %             [I1,~] = KraskovMI(newChain,lowfiOut(:,jj),6);
                %             [I1,~] = KraskovMI(newChain,relowfiOut(:,jj),6);
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

            %rcVol = (currentPts(end,2)-currentPts(2,2))/currentPts(2,2);
            %delTime = currentPts(end,1)-currentPts(2,1);
            %rcv = rcVol/delTime;
            rcv = abs(lowfi(end)-currentPts(end,2))/(lowfi(end)+currentPts(end,2));
            rcv_vec = [rcv_vec rcv];

            %penalize MI by k*absolute rcg*penalty for skipped points
            for p = 1:length(relMI)
                %score(p) = relMI(p) - abs(rcv)*sum(relMI(1:(p-1)))/sum(relMI);
                %(it is choosing every point because MIs are all nearly
                %identical, so any penalty means choosing the first point)
                score(p) = unifMI(p) - abs(rcv)*sum(unifMI(1:(p-1)))/sum(unifMI);
            end

            disp( strcat( 'abs(rcv):', num2str(abs(rcv)) ) )

            figure(1000)
            plot(expDesigns,score,'ob','MarkerFaceColor','b')
            hold on
            plot(expDesigns,unifMI,'or')
            filename = [path 'ScoreFxn_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
            saveas(gcf,filename)
            close(gcf)

            figure(1002); plot( expDesigns, lowfiOut(:,1:length(expDesigns))' )
            filename = [path 'LowFiOut_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
            saveas(gcf,filename)
            close(gcf)


            %Choose optimal design
            point = expDesigns(1,find(score == max(score)))

            if length(point)>1
                fprintf('Multiple points have same score function. Last point in tie list chosen.')
                point = point(end);
            end


            %Add point to current list of points for next round of calibration
            idx = find(round(fulldata.xdata)==point);
            currentPts = [currentPts; fulldata.xdata(point+2) fulldata.ydata(point+2)]

            %Remove chosen point and all skipped points from experimental design list
            expDesigns(find(expDesigns<=point)) = [];


            clear score relMI miVals lowfiOut
        end

    else %%%% K and beta
        disp( 'K and delta' )

        [param_fit,ss01] = fmincon(@(param)ssq_compareSchemesPSI(param,data),[PSI_fit delta_fit],[],[],[],[],[lb_PSI,lb_delta],[ub_PSI,ub_delta], [], opt);


        %%%%%%%%%need to change this to new prior
        params1 = {  
            {'PSI',param_fit(1), 0.1, 1, mu_PSI_est, sig_PSI_est }
            {'delta',param_fit(2), 0, 0.15, mu_est(nprior), sig_est(nprior) }
            };

        fitdata.xdata = currentPts(:,1);
        fitdata.ydata = currentPts(:,2);
    
        model.ssfun = @ssq_compareSchemesPSI;
        mse = ss01/(length(currentPts(:,2))-2);
        model.sigma2 = mse; %Initial guess for error variance
        options.updatesigma = 1;


        %%%%%%%%%need to change this to new prior
        % prior
        % model.priorfun = @lognpdf;
        lognormal = @(th,mu,sig) prod( exp( -( (log(th)-mu)./(2*sig) ).^2 ) ./(th.*sig*sqrt(2*pi)) );
        model.priorfun = lognormal; %@lognpdf;

        model.priorfun = @(th,mu,sig)-2*log(lognpdf(th(2),mu(2), sig(2))); 
        

        options.waitbar = 0;
        options.nsimu = no_smps/5;
        [results,~,~] = mcmcrun(model,fitdata,params1,options);

        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,fitdata,params1,options,results);


        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        param_fit = chain(ind,:); %This is our fitted parameter set 
        params = [growthParams(1) param_fit];


        %Save off parameter estimates and chains
        delta_val = param_fit(end); 
        delta_vec = [delta_vec; delta_val];
        PSI_vec = [PSI_vec param_fit(1)]; 


        % Generate current model trajectory, credible intervals, and plot
        K = V0 / param_fit(1); 
        params_full = [growthParams(1) K param_fit(2)];
        
        [timeFit, volFit] = headnecktumor_RT(params_full,fulldata);
        lowfi = interp1(timeFit,volFit, fulldata.xdata);

        ln = fulldata.xdata';
        modelfun2 = @(ln,params) tumorfun2(ln,params,fulldata);

        figure(1)
        pred = mcmcpred(results,chain,s2chain,ln,modelfun2,no_smps);
        pred.obslims = [];
        mcmcpredplot(pred)
        hold on
        h=gca;
        plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
        plot(fulldata.xdata,lowfi,'-k','Linewidth',2) %This is your optimal model fit
        plot(currentPts(:,1),currentPts(:,2),'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
        xlabel('Time (days)','FontSize',18);
        ylabel('Volume','FontSize',18);
        titleName = ['Patient #' int2str(npatient) ': Iteration ' num2str(nIter) '*'];
        title(string(titleName),'Interpreter','none','FontSize',14)
        axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
        set(h, 'FontSize',18)
        hold off
        filename = [path 'ModelFit_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        saveas(gcf,filename)
        close(gcf)


        % figure(2)
        % mcmcplot(chain,[],{'K', '\beta'},'chainpanel')
        % filename = [path 'Chain_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        % saveas(gcf,filename)
        % % filename = [path 'Chain_Iteration' num2str(nIter) '.fig'];
        % % saveas(gcf,filename)
        % close(gcf)
        %
        figure(3)
        yy = mcmcplot(chain,[],{'K', '\beta'},'denspanel');
        % filename = [path 'Density_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        % saveas(gcf,filename)
        % % filename = [path 'Density_Iteration' num2str(nIter) '.fig'];
        % % saveas(gcf,filename)
        close(gcf)
        pdf_vec{nIter} = yy{2};



        % Calculate MSE to measure error
        mse = sum((lowfi-fulldata.ydata).^2)/sum(fulldata.ydata.^2);
        err_vec = [err_vec; mse];

        % Calculate uncertainty metrics
        for n = 1:length(ln)
            widthCI(nIter,n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
        end
        areaCI = [areaCI; sum(widthCI(nIter,1:end))];


        % Determine which point should be selected next
        if currentPts(end,1)==fulldata.xdata(end)
            ptsLeft = 0; %this is the final iteration

        else

            knnSet = 10:10:no_smps; %thin out samples for kNN analysis
            newChain = chain(knnSet,:);

            % Calculate predicted values at each experimental design with
            % each parameter set from newChain
            for ii = 1:size(newChain,1)

                %%%%  prior in PSI and delta 

                [tsol, ysol] = headnecktumorPSI_RT([growthParams(1), newChain(ii,1), newChain(ii,2)],fulldata);

                for jj = 1:length(expDesigns)
                    lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
                end
            end

            %Calculate MI for each remaining design

            %%%%% with normalized data
            %         renewChain = newChain/max(newChain);
            %         relowfiOut = lowfiOut/max(lowfiOut(:));

            %%%%% with normalized by mean/std data
            normChain = (newChain-mean(newChain))./std(newChain);
            normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);

            for jj = 1:length(expDesigns)
                %             [I1,~] = KraskovMI(newChain,lowfiOut(:,jj),6);
                %             [I1,~] = KraskovMI(newChain,relowfiOut(:,jj),6);
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

            %rcVol = (currentPts(end,2)-currentPts(2,2))/currentPts(2,2);
            %delTime = currentPts(end,1)-currentPts(2,1);
            %rcv = rcVol/delTime;
            rcv = abs(lowfi(end)-currentPts(end,2))/(lowfi(end)+currentPts(end,2));
            rcv_vec = [rcv_vec rcv];

            %penalize MI by k*absolute rcg*penalty for skipped points
            for p = 1:length(relMI)
                %score(p) = relMI(p) - abs(rcv)*sum(relMI(1:(p-1)))/sum(relMI);
                %(it is choosing every point because MIs are all nearly
                %identical, so any penalty means choosing the first point)
                score(p) = unifMI(p) - abs(rcv)*sum(unifMI(1:(p-1)))/sum(unifMI);
            end

            disp( strcat( 'abs(rcv):', num2str(abs(rcv)) ) )

            % figure(1000)
            % plot(expDesigns,score,'ob','MarkerFaceColor','b')
            % hold on
            % plot(expDesigns,unifMI,'or')
            % filename = [path 'ScoreFxn_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
            % saveas(gcf,filename)
            % close(gcf)
            % 
            % figure(1002); plot( expDesigns, lowfiOut(:,1:length(expDesigns))' )
            % filename = [path 'LowFiOut_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
            % saveas(gcf,filename)
            % close(gcf)


            %Choose optimal design
            point = expDesigns(1,find(score == max(score)))

            if length(point)>1
                fprintf('Multiple points have same score function. Last point in tie list chosen.')
                point = point(end);
            end


            %Add point to current list of points for next round of calibration
            idx = find(round(fulldata.xdata)==point);
            currentPts = [currentPts; fulldata.xdata(point+2) fulldata.ydata(point+2)]

            %Remove chosen point and all skipped points from experimental design list
            expDesigns(find(expDesigns<=point)) = [];


            clear score relMI miVals lowfiOut
        end





    end
    
% Save results for later comparison
file = [path 'Zahid_Results_nPatient' num2str(npatient) '.mat'];
point_list = currentPts; 
save(file, 'point_list','delta_vec','growthParams','PSI_vec','err_vec','rcv_vec','areaCI','pdf_vec')

end


point_list = currentPts; %Final list of points in order of selection
delta_vec = delta_vec(2:end); % remove the first 1. 
PSI_vec = PSI_vec(2:end); %remove the first 10 

% Save results for later comparison
file = [path 'Zahid_Results_nPatient' num2str(npatient) '.mat'];
% save(file, 'point_list','delta_vec','growthParams','err_vec','rcv_vec','areaCI','fitChain','pdf_vec')
save(file, 'point_list','delta_vec','growthParams','PSI_vec','err_vec','rcv_vec','areaCI','pdf_vec')


end



%% Helper functions

function dval = tumorLogistic( t, val, params )

lambda = params(1);
K = params(2);

V = val(1);

dval(1) = lambda * V * ( 1 - V / K );

end


% SSQ function for calibration before RT 
function SSrt = ssq_tumorVol(params,data)
global currentPts growthParams V0 

K = V0 / params(1); 

[time,vol] = ode23(@(t,y)tumorLogistic(t,y,[growthParams, K]), [data.xdata(1),data.xdata(end)], data.ydata(1) ); 

tumVol = interp1(time,vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end

% function [tsol, ysol] = headnecktumorPSI_RT( params, data )
% global V0 
% 
% tsol = data.xdata(1); % convert to days  
% ysol = data.ydata(1); 
% y = ysol; 
% 
% 
% treat = [0:4, 7:11, 14:18, 21:25, 28:32, 35:39, 42:46, 49]/7; 
% treat = [tsol, treat]; 
% 
% K = V0 / params(2); 
% 
% for nt = 1:(length(treat)-1)
% 
%     if( (treat(nt+1) - treat(nt)) ~= 0 )
% 
%             [t,y] = ode23(@(t,y)tumorLogistic(t,y,[params(1), K]),[treat(nt),treat(nt+1)],y(end));
%              K = K*(1-params(3)); 
% 
%     end
%     ysol = [ysol; y(2:end)]; tsol = [tsol; t(2:end)];
% 
% end
% 
% end


%SSQ function for calibration of RT parameter
function SSrt = ssq_compareSchemes(params, data)
global growthParams V0 

[time,vol] = headnecktumor_RT( [growthParams params(end)], data);

tumVol = interp1(time, vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end


%SSQ function for calibration of RT parameter
function SSrt = ssq_compareSchemesPSI(params, data)
global growthParams V0 

K = V0 / params(1); 

[time,vol] = headnecktumor_RT( [growthParams(1) K params(end)], data);

tumVol = interp1(time, vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end


% 
% % SSQ function for calibration of RT parameter
% function SSrt = ssq_tumorVolpostRT(params,data)
% 
% global currentPts modelType growthParams
% 
% [time,vol] = tumorModel_postRT_wIC(modelType,[growthParams params],[data.xdata(1) data.ydata(1)] );
% 
% tumVol = interp1(time,vol, currentPts(:,1));
% 
% SSrt = sum((tumVol - currentPts(:,2)).^2);
% 
% end

% function SSrt = ssq_tumorVolpostRTnK(params,data)
% 
% global currentPts modelType growthParams
% 
% [time,vol] = tumorModel_postRT_wIC(modelType,[growthParams(1) params],[data.xdata(1) data.ydata(1)] );
% 
% tumVol = interp1(time,vol, currentPts(:,1));
% 
% SSrt = sum((tumVol - currentPts(:,2)).^2);
% 
% end


% Function for credible interval plotting
function v=tumorfun(timef,params,data)

global modelType growthParams


[tsol, ysol] = headnecktumor_RT([growthParams params(end)],data);

v = interp1(tsol,ysol,timef);

end

function v=tumorfun2(timef,params,data)

global modelType growthParams V0 

K = V0 / params(1); 

[tsol, ysol] = headnecktumor_RT([growthParams(1) K params(end)],data);

v = interp1(tsol,ysol,timef);

end

function [mu, sig] = convertlogn( mmt )

mu  = log( mmt(:,1).^2 ./ sqrt( mmt(:,1).^2 + mmt(:,2).^2 ) );
sig = sqrt( log( 1 + mmt(:,2).^2 ./ mmt(:,1).^2  ));

end



