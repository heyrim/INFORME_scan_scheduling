
% Heyrim Cho <heyrim.cho@asu.edu>
% edited: 05-31-2024, 02-07-2025 log-normal prior added - start low-fidelity data point from prior
%         03-21-2025 prior of A, B, beta, using single data point

function point_list = seqDesign_CAData_prior_all(npatient,firstpoint)

global currentPts modelType growthParams M dose 

disp( strcat( 'Start patient #', int2str(npatient) ) )

modelType = 'LOG+DVR';
dose = 2; 

ncaseprior = 4; % 1=unif, 2=single lognorm, 3=3 types of lognorm, 4=bdd lognorm
point_list = [];

%% load data

load('/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/ABM_cancer_radiotherapy/data/CA_rad1_Cm_22_30_pNR_all.mat')

fulldata.xdata = data.xdata(1:55);
fulldata.ydata = data.ydata(1:55,npatient); %use only tumor volume during treatment period
clear data

if ncaseprior == 4
    path = ['../Figures/CAdata_RT_betaonly_corrected/']; 
end
mkdir( path )


%Settings
no_smps = 10000; %DRAM chain lengths

%Set up vectors for saved quantities;
beta_vec = []; %store parameter estimates at each iteration
K_vec = []; 
err_vec = []; %store MSE at each iteration
widthCI = []; %store width of credible interval at each iteration
areaCI = []; %store area of credible interval at each iteration
rcv_vec = [];

%Initial starting data - days 15 were used to get intrinsic growth
currentPts = [ fulldata.xdata(16) fulldata.ydata(16);];
icdata = currentPts;


point = firstpoint; % precomputed point 
beta_val = 0.0333; %0.0197; 
beta_vec = [beta_vec; beta_val];

% prior bounded logit normal using all 135 data
P = readtable( '../data/Update files 02-28-2025/prostate cancer model - large - general - bounded/populationParameters.txt');
pop_est   = P{1:3,2}; %population level estimates for A, B, gamma
clear P 
growthParams = pop_est(1:2)';


%Add point to current list of points for next round of calibration
idx = find(fulldata.xdata==point);
currentPts = [currentPts; point fulldata.ydata(idx)];

expDesigns = fulldata.xdata((idx+1):end);


%% Change prior depending on the first data point

data.xdata = currentPts(:,1);
data.ydata = currentPts(:,2);

lb = [0];
ub = [1];
opt = optimset('Display','off');
% opt = optimoptions('fmincon', 'Display', 'iter', 'ObjectiveLimit', 1e-6, 'FunctionTolerance', 1e-8);
[beta_val,ss01] = fmincon( @(param)ssq_tumorVolpostRT_beta(param,data),beta_val,[],[],[],[],lb,ub, [], opt);


[tsol, ysol] = tumorModel_postRT_wIC([growthParams,beta_val],currentPts(1,:));
lowfi = interp1(tsol,ysol, fulldata.xdata(16:end));


yratio = currentPts(2,2)/currentPts(1,2); 
if( yratio < 0.4 ); nprior = 3; elseif( yratio <= 0.85 ); nprior = 4; else; nprior = 5; end


% lognormal three class of priors
P = readtable( '../data/Update files 02-27-2025/Monolix and MATLAB - fit and distribution/prostate cancer model - large - Cov/populationParameters.txt');
pop_est   = P{1:3,2}; %population level estimates for A, B, gamma
class_est = P{4:5,2}; %this is used to determine the population gamma for the other groups
sd_est    = P{6:8,2}; %standard deviation of random effect
clear P

growthParams = pop_est(1:2)';
K_vec = []; 

%calculate the population estimates of gamma for the other groups
pop_est(4) = pop_est(3)*exp(class_est(1));
pop_est(5) = pop_est(3)*exp(class_est(2));
sd_est(4) = sd_est(3); sd_est(5) = sd_est(3);

                                
mu_est(1:2) = log( pop_est(1:2) ); sig_est(1:2) = sd_est(1:2); 
betascale = 0.0211; 

mu_est(3:5) = log( pop_est(3:5)' * betascale ); 
sig_est(3:5) = sd_est(3:5); 

%% Now start the sequential design procedure, estimating beta + other param if needed

nIter = 0; %keep track of how many calibration iterations have been performed
ptsLeft = 1; %binary tracker; 1 means there are still points left to choose


while ptsLeft == 1 && nIter < 5

    nIter = nIter+1;

    data.xdata = currentPts(:,1);
    data.ydata = currentPts(:,2);

    iii = []; for ii = 1:length(currentPts(:,1)); iii = [iii; find( fulldata.xdata(16:end)==currentPts(ii,1) ) ]; end 
    last_fit_error = max( sum( abs( currentPts(:,2) - lowfi( iii )' ) ) / size(currentPts,1), ... 
        abs( currentPts(end,2) - lowfi(iii(end) ) ) ); 

    %% make the decision to estimate beta only or beta + K
    % 250531 results with thresh 0.02 
    % fit_thresh = 0.02;
    fit_thresh = 0.05;

    if( last_fit_error < fit_thresh ) %%%% beta only
        disp( 'beta only' )

        params1 = {
            {'beta',beta_val, eps, 1, mu_est(nprior), sig_est(nprior) }
            };

        model.ssfun = @ssq_tumorVolpostRT_beta;
        
        ss0 = ssq_tumorVolpostRT([growthParams,beta_val],data); 
        mse = ss0/(length(currentPts(:,2))-length(beta_val));
        model.sigma2 = mse; %Initial guess for error variance

        options.updatesigma = 1;

        %%%%%%%%%need to change this to new prior
        % prior
        model.priorfun = @(th,mu,sig)-2*log(prod(lognpdf(th,mu,sig)));
  
        options.waitbar = 0;
        options.nsimu = no_smps/5;
        [results,chain,s2chain] = mcmcrun(model,data,params1,options);

        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);


        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        % beta_val = chain(ind,:); %This is our fitted parameter set\
        % params = [growthParams beta_val];

        params = chain(ind,:); 
        beta_val = params(end); 

        %Save off parameter estimates and chains
        beta_vec = [beta_vec; beta_val];

        % Generate current model trajectory, credible intervals, and plot
        [timeFit, volFit] = tumorModel_postRT_wIC([growthParams,params],icdata);
        
        lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));



        ln = fulldata.xdata(16:end)';
        modelfun = @(ln,params) tumorfun(ln,params,icdata);

        figure(1); subplot( 1, 5, nIter );  hold off; 
        pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
        % pred.obslims = [];
        mcmcpredplot(pred)
        hold on
        h=gca;
        plot(fulldata.xdata(16:end),fulldata.ydata(16:end),'ok','MarkerSize',6) %All possible scans
        plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',2) %This is your optimal model fit
        plot(currentPts(:,1),currentPts(:,2),'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
        xlabel('Time (days)','FontSize',18);
        ylabel('Volume','FontSize',18);
        if( nIter == 1 )
        titleName = ['Patient #' int2str(npatient) ': Iteration ' num2str(nIter)];
        else
        titleName = ['Iteration ' num2str(nIter)];
        end
        title(string(titleName),'Interpreter','none','FontSize',14)
        axis([min(fulldata.xdata(16))-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
        set(h, 'FontSize',18)
        hold off
        filename = [path 'ModelFit_nPatient' num2str(npatient) '.jpg'];
        saveas(gcf,filename)
        % close(gcf)


        % figure(2)
        % mcmcplot(chain,[],{'\beta'},'chainpanel')
        % filename = [path 'Chain_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        % saveas(gcf,filename)
        % % filename = [path 'Chain_Iteration' num2str(nIter) '.fig'];
        % % saveas(gcf,filename)
        % close(gcf)
        %
        % figure(3)
        % yy = mcmcplot(chain,[],{'\beta'},'denspanel');
        % yy = mcmcplot(chain,[],{'\lambda', 'K', '\beta'},'denspanel');
        % filename = [path 'Density_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        % saveas(gcf,filename)
        % % filename = [path 'Density_Iteration' num2str(nIter) '.fig'];
        % % saveas(gcf,filename)
        % close(gcf)
        % pdf_vec{nIter} = yy{1};


        % Calculate MSE to measure error
        mse = sum((lowfi'-fulldata.ydata(16:end)).^2)/numel(fulldata.ydata(16:end));
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

                [tsol, ysol] = tumorModel_postRT_wIC([growthParams newChain(ii,:)],icdata);

                for jj = 1:length(expDesigns)
                    lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
                end
            end

            %Calculate MI for each remaining design
            %%%%% with normalized by mean/std data
            normChain = (newChain-mean(newChain))./std(newChain);
            normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);

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
            idx = find(fulldata.xdata==point);
            currentPts = [currentPts; point fulldata.ydata(idx(1))];

            %Remove chosen point and all skipped points from experimental design list
            expDesigns(find(expDesigns<=point)) = [];


            clear score relMI miVals lowfiOut
        end

    else %%%% K and beta
        disp( 'K and beta' )

        %%%%%%%%%need to change this to new prior
        params1 = {
            % {'lambda',pop_est(1), eps, 1, pop_est(1), sd_est(1) }
            {'K',growthParams(2), eps, 1, mu_est(2), sig_est(2) }
            {'beta',beta_val, eps, 1, mu_est(nprior), sig_est(nprior) }
            };

        model.ssfun = @ssq_tumorVolpostRTnK;
        options.updatesigma = 1;

        %%%%%%%%% need to change this to new prior
        % prior
        % lognormal = @(th,mu,sig) prod( exp( -( (log(th)-mu)./(2*sig) ).^2 ) ./(th.*sig*sqrt(2*pi)) );
        % model.priorfun = @(th,mu,sig)-2*log(lognormal(th,mu,sig) );  
        model.priorfun = @(th,mu,sig)-2*log(prod(lognpdf(th,mu,sig)));

        options.waitbar = 0;
        options.nsimu = no_smps/5;
        [results,chain,s2chain] = mcmcrun(model,data,params1,options);

        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);


        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        param_fit = chain(ind,:); %This is our fitted parameter set 
        params = [growthParams(1) param_fit];


        %Save off parameter estimates and chains
        beta_val = param_fit(end); 
        beta_vec = [beta_vec; beta_val];
        K_vec(nIter) = param_fit(1); 
        growthParams(2) = param_fit(1); 

        % fitChain{nIter}.paramschain = chain;
        % fitChain{nIter}.s2chain = s2chain;
        % fitChain{nIter}.ss2chain = ss2chain;


        % Generate current model trajectory, credible intervals, and plot
        [timeFit, volFit] = tumorModel_postRT_wIC(params,icdata);
        lowfi = interp1(timeFit,volFit, fulldata.xdata(16:end));

        ln = fulldata.xdata(16:end)';
        modelfun2 = @(ln,params) tumorfun2(ln,params,icdata);

        figure(1); subplot( 1, 5, nIter ); hold off; 
        pred = mcmcpred(results,chain,s2chain,ln,modelfun2,no_smps);
        pred.obslims = [];
        mcmcpredplot(pred)
        hold on
        h=gca;
        plot(fulldata.xdata(16:end),fulldata.ydata(16:end),'ok','MarkerSize',6) %All possible scans
        plot(fulldata.xdata(16:end),lowfi,'-k','Linewidth',2) %This is your optimal model fit
        plot(currentPts(:,1),currentPts(:,2),'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
        xlabel('Time (days)','FontSize',18);
        ylabel('Volume','FontSize',18);
        titleName = ['Patient #' npatient ': Iteration ' num2str(nIter) '*'];
        title(string(titleName),'Interpreter','none','FontSize',14)
        axis([min(fulldata.xdata(16))-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
        set(h, 'FontSize',18)
        hold off
        % filename = [path 'ModelFit_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        filename = [path 'ModelFit_nPatient' num2str(npatient) '.jpg'];
        saveas(gcf,filename)
        % close(gcf)


        % figure(2)
        % mcmcplot(chain,[],{'K', '\beta'},'chainpanel')
        % filename = [path 'Chain_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        % saveas(gcf,filename)
        % % filename = [path 'Chain_Iteration' num2str(nIter) '.fig'];
        % % saveas(gcf,filename)
        % close(gcf)
        %
        % figure(3)
        % yy = mcmcplot(chain,[],{'K', '\beta'},'denspanel');
        % filename = [path 'Density_nPatient' num2str(npatient) '_Iteration' num2str(nIter) '.jpg'];
        % saveas(gcf,filename)
        % % filename = [path 'Density_Iteration' num2str(nIter) '.fig'];
        % % saveas(gcf,filename)
        % close(gcf)
        % pdf_vec{nIter} = yy{1};


        % Calculate MSE to measure error
        mse = sum((lowfi'-fulldata.ydata(16:end)).^2)/numel(fulldata.ydata(16:end));
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

                [tsol, ysol] = tumorModel_postRT_wIC([growthParams(1) newChain(ii,:)],icdata);

                for jj = 1:length(expDesigns)
                    lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
                end
            end

            %Calculate MI for each remaining design
            %%%%% with normalized by mean/std data
            normChain = (newChain-mean(newChain))./std(newChain);
            normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);

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
            idx = find(fulldata.xdata==point);
            currentPts = [currentPts; point fulldata.ydata(idx(1))]; 

            %Remove chosen point and all skipped points from experimental design list
            expDesigns(find(expDesigns<=point)) = [];


            clear score relMI miVals lowfiOut
        end



    end

end


point_list = currentPts %Final list of points in order of selection

% Save results for later comparison
file = [path 'Results_nPatient' num2str(npatient) '.mat'];
% save(file, 'point_list','beta_vec','growthParams','err_vec','rcv_vec','areaCI','fitChain','pdf_vec')
save(file, 'point_list','beta_vec','K_vec','err_vec','rcv_vec','areaCI', 'nprior')

figure(1); clf


end



%% Helper functions


% SSQ function for calibration of RT parameter
function SSrt = ssq_tumorVolpostRT(params,data)

global modelType growthParams

% [time,vol] = tumorModel_postRT_wIC([growthParams params],[data.xdata(1) data.ydata(1)] );
[time,vol] = tumorModel_postRT_wIC( params,[data.xdata(1) data.ydata(1)] );

tumVol = interp1(time,vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end

function SSrt = ssq_tumorVolpostRT_beta(params,data)

global modelType growthParams currentPts 

[time,vol] = tumorModel_postRT_wIC([growthParams params],[data.xdata(1) data.ydata(1)] );

% tumVol = interp1(time,vol, data.xdata);
tumVol = interp1(time,vol, currentPts(:,1));

% SSrt = sum((tumVol - data.ydata).^2);
SSrt = sum((tumVol - currentPts(:,2)).^2);

end


function SSrt = ssq_tumorVolpostRTnK(params,data)

global  modelType growthParams currentPts 

[time,vol] = tumorModel_postRT_wIC([growthParams(1) params],[data.xdata(1) data.ydata(1)] );

tumVol = interp1(time,vol, data.xdata);

SSrt = sum((tumVol - data.ydata).^2);

end


% Function for credible interval plotting
function v=tumorfun(timef,params,icdata)

global  growthParams currentPts 


[tsol, ysol] = tumorModel_postRT_wIC([growthParams params],icdata);
% % % % [tsol, ysol] = tumorModel_postRT_wIC(params,icdata);

v = interp1(tsol,ysol,timef);

end

function v=tumorfun2(timef,params,icdata)

global  growthParams


[tsol, ysol] = tumorModel_postRT_wIC([growthParams(1) params(1) params(2)],icdata);

v = interp1(tsol,ysol,timef);

end

function [mu, sig] = convertlogn( mmt )

mu  = log( mmt(:,1).^2 ./ sqrt( mmt(:,1).^2 + mmt(:,2).^2 ) );
sig = sqrt( log( 1 + mmt(:,2).^2 ./ mmt(:,1).^2  ));

end



