


nallcompute = 0; 


%% Compute first data point to collect from group prior 
if( nallcompute )
    firstpoint = seqDesign_CAData_prior_firstdatapoint; 
else
    firstpoint = 27; 

end 

%% Compute adaptive schedule for each patients 
for npatient = 1:135 
    scan_schedule_K_beta{npatient} = seqDesign_CAData_prior_all(npatient,firstpoint); 
end

%% plot results 
load('/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/ABM_cancer_radiotherapy/data/CA_rad1_Cm_22_30_pNR_all.mat')
yratio = data.ydata(27,:)./data.ydata(16,:);
[~,indnp] = sort( yratio );

if( yratio < 0.4 ); nprior = 3; elseif( yratio <= 0.85 ); nprior = 4; else; nprior = 5; end


nprior = zeros( 1, 135 ); 
for n = 1:135 
    if( yratio(n) < 0.4 ); nprior(n) = 3; elseif( yratio(n) <= 0.85 ); nprior(n) = 4; else; nprior(n) = 5; end 
end 

schedulemap = zeros( 55, 135 ); 
for n = 1:135 

    % load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT_betaonly/Results_nPatient' int2str(n) '.mat'])
    % load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT/Results_nPatient' int2str(n) '.mat'])
    % load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT_population/Results_nPatient' int2str(n) '.mat'])
    load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT_unif/Results_nPatient' int2str(n) '.mat'])

    errcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) );
    relerrcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) )./beta_vec(2:end);
    ind1 = find( errcheck < 0.01, 1 );
    ind2 = find( relerrcheck < 0.05, 1 );
    if( isempty(ind1) ); ind1=length(beta_vec)-1;  end
    if( isempty(ind2) ); ind2=length(beta_vec)-1;  end
    ind = max( ind1, ind2 );

    nscan(n) = max(ind,2); 

    %%%%% for just RT 
    % schedulemap( point_list( 1:max(ind+1,3), 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
    %%%%% for _unif and population 
    schedulemap( point_list( 1:max(ind+2,3), 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
    

end

figure; heatmap( schedulemap(:,indnp) )

for nn = 1:3 
ii = find( nprior == nn+2 ); 
jnd{nn} = []; 
for n = 1:length(indnp) 
    if( intersect(ii, indnp(n)) )
        jnd{nn} = [jnd{nn}, indnp(n)]; 
    end
end
schedule_CAdadta{nn} = schedulemap(:,jnd{nn}); 
nscanchoice(nn) = mean( nscan(jnd{nn}) )
end 
% result : nscanchoice =    2.0192    2.2727    2.3846 -> 2, 3, 3 
% % save( 'schedule_CAdadta.mat', 'schedule_CAdadta' ) 

figure; 
for n = 1:3; subplot( 1, 3, n ); 
    heatmap( mean( schedule_CAdadta{n}' )' );
    colorbar off 
    clim( [0 0.6] )
end


nscanchoice = [2, 3, 3]; 
finalschedule_CAdata(1:55,n) = 0; 
for n = 1:3
    meanschedule(:,n) = mean( schedule_CAdadta{n}' )';  
    [ii,jj] = sort( meanschedule(:,n), 'descend' ); 
    finalschedule_CAdata(1:55,n) = 0; 
    for nn = 1:(nscanchoice(n)+1)
        finalschedule_CAdata(jj(nn),n) = 1; 
    end

end 
% final schedule of group 3 is a clustered one 
finalschedule_CAdata(jj(nn),n) = 0; 
finalschedule_CAdata(jj(nn+2),n) = 1; 


% ii = cumsum( [0, 56, 43, 36] ); 
% 
% for nn = 3:5 
%     figure; heatmap( schedulemap(:, indnp( (ii(nn-2)+1):ii(nn-1) ) ) )
% end 

