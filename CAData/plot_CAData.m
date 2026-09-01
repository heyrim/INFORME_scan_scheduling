function plot_CAData (error_thresh, rel_thresh, del_thresh, nfig )

%% data and sorting 
load('/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/ABM_cancer_radiotherapy/data/CA_rad1_Cm_22_30_pNR_all.mat')
yratio = data.ydata(27,:)./data.ydata(16,:);
[~,indnp] = sort( yratio );

nprior = zeros( 1, 135 ); 
for n = 1:135 
    if( yratio(n) < 0.4 ); nprior(n) = 3; elseif( yratio(n) <= 0.85 ); nprior(n) = 4; else; nprior(n) = 5; end 
end 
ngroup = nprior; 


% indnp(1:52), indnp(53:96), indnp(97:135)
ii = [1, 52; 53, 96; 97, 135]; 

schedulemap1_raw = zeros( 55, 135 ); 
for n = 1:135 
    load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT_unif/Results_nPatient' int2str(n) '.mat'])
    point_list( :, 1) = point_list( :, 1) + 1; 
    schedulemap1_raw( point_list( :, 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
end 
figure(1); heatmap( schedulemap1_raw(:,indnp) ); colorbar off; ylabel( 'scan time (days)' ); xlabel( 'patient' )


schedulemap2_raw = zeros( 55, 135 ); 
for n = 1:135 
    load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT_population/Results_nPatient' int2str(n) '.mat'])
    point_list( :, 1) = point_list( :, 1) + 1; 
    schedulemap2_raw( point_list( :, 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
end 
figure(2); heatmap( schedulemap2_raw(:,indnp) ); colorbar off; ylabel( 'scan time (days)' ); xlabel( 'patient' )


schedulemap3_raw = zeros( 55, 135 ); 
for n = 1:135 
    load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT/Results_nPatient' int2str(n) '.mat'])
    point_list( :, 1) = point_list( :, 1) + 1; 
    schedulemap3_raw( point_list( :, 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
end 
figure(3); heatmap( schedulemap3_raw(:,indnp) ); colorbar off; ylabel( 'scan time (days)' ); xlabel( 'patient' )


xlabel( 'scan time (days)' ) 


error_thresh = 0.01; 
rel_thresh = 0.1; 
del_thresh = 0.001; 

schedulemap1 = zeros( 55, 135 ); 
errmap1 = zeros( 55, 135 ); 
for n = 1:135 

    load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT_unif/Results_nPatient' int2str(n) '.mat'])
    point_list( :, 1) = point_list( :, 1) + 1; 

    if( size( point_list, 1 ) - length(beta_vec) == 2 ); beta_vec = [beta_vec;beta_vec(end)]; end 
    if( size( beta_vec, 1) ~= size( err_vec, 1 ) ); err_vec = [err_vec;err_vec(end)]; end 

    errcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) );
    relerrcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) )./beta_vec(2:end);
    delerrcheck = abs( (beta_vec(1:end-1) - beta_vec(2:end)) ./ (point_list(2:end-1,1)-point_list(3:end,1))); 
    
    % figure(152); subplot( 3, 3, ngroup(n)-2 ); hold on; plot( errcheck, 'k:' ); 
    % figure(152); subplot( 3, 3, ngroup(n)-2+3 ); hold on; plot( relerrcheck, 'k:' ); 
    % figure(152); subplot( 3, 3, ngroup(n)-2+6 ); hold on; plot( delerrcheck, 'k:' ); 

    ind1 = find( errcheck < error_thresh, 1 );
    ind2 = find( relerrcheck < rel_thresh, 1 );
    ind3 = find( delerrcheck < del_thresh, 1 );
    
    if( isempty(ind1) ); ind1=length(beta_vec)-1;  end
    if( isempty(ind2) ); ind2=length(beta_vec)-1;  end
    if( isempty(ind3) ); ind3=length(beta_vec)-1;  end
    % ind = max( ind1, ind2 ); %ind = max( ind, ind3 ); 
    ind = min( ind1, ind2 ); ind = min( ind, ind3 ); 

    nscan(n) = ind + 1; 

    %%%%% for just RT 
    % schedulemap( point_list( 1:max(ind+1,3), 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
    %%%%% for _unif and population 
    schedulemap1( point_list( 1:(ind+2), 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 

    errmap1( point_list( 2:(ind+2), 1), n  ) = err_vec( 1:(ind+1) ); 
    errmap1_scan( 1:(ind+1), n ) = err_vec( 1:(ind+1) ); 

    
end
figure(nfig+1); heatmap( schedulemap1(:,indnp) ); colorbar off; ylabel( 'scan time (days)' ); xlabel( 'patient' )

disp( ['mean of unif' num2str( mean(nscan) )] )

schedulemap2 = zeros( 55, 135 ); 
errmap2 = zeros( 55, 135 ); 
for n = 1:135 

    load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT_population/Results_nPatient' int2str(n) '.mat'])
    point_list( :, 1) = point_list( :, 1) + 1; 

    if( size( point_list, 1 ) - length(beta_vec) == 2 ); beta_vec = [beta_vec;beta_vec(end)]; end 
    errcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) );
    relerrcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) )./beta_vec(2:end);
    delerrcheck = abs( (beta_vec(1:end-1) - beta_vec(2:end)) ./ (point_list(2:end-1,1)-point_list(3:end,1))); 

    % figure(153); subplot( 3, 3, ngroup(n)-2 ); hold on; plot( errcheck, 'k:' ); 
    % figure(153); subplot( 3, 3, ngroup(n)-2+3 ); hold on; plot( relerrcheck, 'k:' ); 
    % figure(153); subplot( 3, 3, ngroup(n)-2+6 ); hold on; plot( delerrcheck, 'k:' ); 
    
    ind1 = find( errcheck < error_thresh, 1 );
    ind2 = find( relerrcheck < rel_thresh, 1 );
    ind3 = find( delerrcheck < del_thresh, 1 );
    
    if( isempty(ind1) ); ind1=length(beta_vec)-1;  end
    if( isempty(ind2) ); ind2=length(beta_vec)-1;  end
    if( isempty(ind3) ); ind3=length(beta_vec)-1;  end
    % ind = max( ind1, ind2 ); %ind = max( ind, ind3 ); 
    ind = min( ind1, ind2 ); ind = min( ind, ind3 ); 

    nscan(n) = ind+1; 

    %%%%% for just RT 
    % schedulemap( point_list( 1:max(ind+1,3), 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
    %%%%% for _unif and population 
    schedulemap2( point_list( 1:ind+2, 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
    
    errmap2( point_list( 2:(ind+2), 1), n  ) = err_vec( 1:(ind+1) ); 
    errmap2_scan( 1:(ind+1), n ) = err_vec( 1:(ind+1) ); 
end
figure(nfig+2); heatmap( schedulemap2(:,indnp) ); colorbar off; ylabel( 'scan time (days)' ); xlabel( 'patient' )
disp( ['mean of pop' num2str( mean(nscan) )] )


schedulemap3 = zeros( 55, 135 ); 
errmap3 = zeros( 55, 135 ); 
for n = 1:135 

    load(['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel_public/Figures/CAdata_RT/Results_nPatient' int2str(n) '.mat'])
    point_list( :, 1) = point_list( :, 1) + 1; 

    beta_vec = beta_vec(2:end); 
    if( size( point_list, 1 ) - length(beta_vec) == 2 ); beta_vec = [beta_vec;beta_vec(end)]; 
    elseif( size( point_list, 1 ) - length(beta_vec) == 0 ); beta_vec = beta_vec(1:end-1); 
    end 
    errcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) );
    relerrcheck = abs( beta_vec(1:end-1) - beta_vec(2:end) )./beta_vec(2:end);
    delerrcheck = abs( (beta_vec(1:end-1) - beta_vec(2:end)) ./ (point_list(2:end-1,1)-point_list(3:end,1))); 
    
    % figure(154); subplot( 3, 3, ngroup(n)-2 ); hold on; plot( errcheck, 'k:' ); 
    % figure(154); subplot( 3, 3, ngroup(n)-2+3 ); hold on; plot( relerrcheck, 'k:' ); 
    % figure(154); subplot( 3, 3, ngroup(n)-2+6 ); hold on; plot( delerrcheck, 'k:' ); 


    ind1 = find( errcheck < error_thresh, 1 );
    ind2 = find( relerrcheck < rel_thresh, 1 );
    ind3 = find( delerrcheck < del_thresh, 1 );
    
    if( isempty(ind1) ); ind1=length(beta_vec)-1;  end
    if( isempty(ind2) ); ind2=length(beta_vec)-1;  end
    if( isempty(ind3) ); ind3=length(beta_vec)-1;  end
    % ind = max( ind1, ind2 ); %ind = max( ind, ind3 ); 
    ind = min( ind1, ind2 ); ind = min( ind, ind3 ); 

    nscan(n) = ind+1; 

    %%%%% for just RT 
    % schedulemap( point_list( 1:max(ind+1,3), 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
    %%%%% for _unif and population 
    schedulemap3( point_list( 1:max(ind+2,3), 1), n ) = 1; % later for not growing tumor, result is too similar to prior mean 0.0197, so should include another point 
    
    errmap3( point_list( 2:(ind+2), 1), n  ) = err_vec( 1:(ind+1) ); 
    errmap3_scan( 1:(ind+1), n ) = err_vec( 1:(ind+1) ); 
end
% figure(nfig+3); heatmap( schedulemap3(:,indnp) ); colorbar off; ylabel( 'scan time (days)' ); xlabel( 'patient' )
disp( ['mean of group' num2str( mean(nscan) )] )


figure(nfig+3); subplot( 1, 3, 1 ); ylabel( 'scan time (days)' ); 
for n = 1:3; subplot( 1, 3, n ); 
    heatmap( schedulemap3(:,indnp( ii(n,1):ii(n,2) )), 'CellLabelColor', 'none'); colorbar off;  xlabel( 'patient' )
    ax = gca;
    ax.XData = ii(n,1):ii(n,2);
    ax.YData = 0:54; 
end 

for n = 1:3 
    figure(nfig+n); title( num2str([error_thresh, rel_thresh, del_thresh]) ); 
end 

figure; 
subplot( 1, 5, 1 ); heatmap( mean( schedulemap1, 2 ), 'CellLabelColor', 'none' ); colorbar off; clim( [0 0.1] ); ax = gca; ax.YData = 0:54; 
subplot( 1, 5, 2 ); heatmap( mean( schedulemap2, 2 ), 'CellLabelColor', 'none' ); colorbar off; clim( [0 0.2] ); ax = gca; ax.YData = 0:54;  
for n = 1:3; subplot( 1, 5, n+2 ); heatmap( mean( schedulemap3(:,indnp( ii(n,1):ii(n,2) )), 2 ), 'CellLabelColor', 'none' ); colorbar off; clim( [0 0.2] ); ax = gca; ax.YData = 0:54;  end
ll = {'Uniform', 'Population', 'Group 1', 'Group 2', 'Group 3'};
for n = 1:5; subplot( 1, 5, n ); title( ll{n} ); end 


%%%% error plots 
figure(nfig+11); subplot( 1, 3, 1 ); plot( errmap1_scan, '-xk' ); set(gca,'yscale', 'log' ); ylim( [0.0001, 0.01] )
figure(nfig+11); subplot( 1, 3, 2 ); plot( errmap2_scan, '-xk' ); set(gca,'yscale', 'log' ); ylim( [0.0001, 0.01] )
figure(nfig+11); subplot( 1, 3, 3 ); plot( errmap3_scan, '-xk' ); set(gca,'yscale', 'log' ); ylim( [0.0001, 0.01] )
subplot( 1, 3, 1 ); hold on; for n = 1:6; plot( n, mean( errmap1_scan( n, errmap1_scan( n, :)~= 0 ) ), '+' ); end 
subplot( 1, 3, 2 ); hold on; for n = 1:5; plot( n, mean( errmap2_scan( n, errmap2_scan( n, :)~= 0 ) ), 'x' ); end 
subplot( 1, 3, 3 ); hold on; for n = 1:5; plot( n, mean( errmap3_scan( n, errmap3_scan( n, :)~= 0 ) ), 'o' ); end 


figure(nfig+12); subplot( 1, 5, 1 ); plot( errmap1_scan, '-xk' ); set(gca,'yscale', 'log' );  
subplot( 1, 5, 2 ); plot( errmap2_scan, '-xk' ); set(gca,'yscale', 'log' );  
for m = 1:3; subplot( 1, 5, m+2 ); plot( errmap3_scan(:, indnp( ii(m,1):ii(m,2) ) ), '-xk' ); set(gca,'yscale', 'log' );   end 
ll = {'Uniform', 'Population', 'group 1', 'group 2', 'group 3'}; 
for n = 1:5; subplot( 1, 5, n ); ylim( [0.0001, 0.02] ); grid on; set(gca, 'FontSize',12); title( ll{n} ); xticks(1:3); xlabel( 'Scans'); end
subplot( 1, 5, 1 ); ylabel( 'Error' ); 
subplot( 1, 5, 2 ); xticks(1:4);


for n = 1:3; tmp(n,1) = mean( errmap1_scan( n, errmap1_scan( n, :)~= 0 ) ); end 
for n = 1:4; tmp(n,2) = mean( errmap2_scan( n, errmap2_scan( n, :)~= 0 ) ); end 
for n = 1:3; tmp(n,3) = mean( errmap3_scan( n, errmap3_scan( n, :)~= 0 ) ); end
% for m = 1:3 
%     errtmp = errmap3_scan( :, indnp( ii(m,1):ii(m,2) ) ); 
%     for n = 1:3; tmp(n,m+2) = mean( errtmp( n, errtmp( n, :)~= 0 ) ); end
% end
figure; hold on; 
plot( tmp(1:3,1), '-x', 'linewidth', 2 ); plot( tmp(1:4,2), '--x', 'linewidth', 2 ); 
plot( tmp(1:3,3), '-.x', 'linewidth', 2 ); 
legend( {'Uniform', 'Population', 'Group'})
ylabel( 'Error' ); xlabel( 'Scan number' )
set(gca,'yscale', 'log' ); ylim( [0.0001, 0.01] ); grid on; box on; set(gca, 'FontSize',14)


end 

% ind = max( ind1, ind2 ); with 0.01, 0.1, otherwise max will not work 


% scan number 
% mean of unif 2.3556
% mean of pop 2.1926
% mean of group 2.0963 - each group  2.0192e+00,   2.2500e+00,    2.0256e+00

% error average 
%    unif         pop         group 
   % 3.0029e-03   1.1449e-03   1.1302e-03
   % 1.8181e-03   8.0275e-04   8.5428e-04
   % 1.3492e-03   6.8325e-04   8.4032e-04
   %          0   5.7571e-04            0

% points_group{1} = [ 27, 28 ]; points_group{2} = [ 27, 29, 48 ]; points_group{3} = [ 27, 31 ];
% points_pop{1} = [ 27, 28, 48  ]; points_pop{2} = [ 27, 28, 48  ]; points_pop{3} = [ 27, 28, 48  ];
% points_unif{1} = [ 22, 23, 29  ]; points_unif{2} = [ 22, 23, 29 ]; points_unif{3} = [ 22, 23, 29 ];

% % errmean{1,:}
% errmat{1} = [ 1.5691e-03   1.4545e-03   0 
%               1.5692e-03   1.4544e-03   1.4200e-03
%               1.5402e-03   1.5093e-03   1.4793e-03 ];
% 
% % errmean{2,:}
%  errmat{2} = [  2.3774e-03   1.3205e-03   1.0155e-03
%                2.3776e-03   1.4096e-03   1.2747e-03
%                1.9465e-03   1.4894e-03   1.5601e-03 ];
% 
% % errmean{3,:}
% errmat{3} = [  1.8666e-03   6.7514e-04   0 
%                1.0341e-03   8.8662e-04   7.5924e-04
%                2.7197e-03   1.6594e-03   1.5691e-03 ];

% unif mean [2.0688e-03   1.5527e-03   1.5362e-03] 
% pop mean [1.6603e-03   1.2502e-03   1.1513e-03]

% schpts{1,1} = [1:2]; schpts{2,1} = [1:3]; schpts{3,1} = [1:2]; 
% for nn = 2:3; for m = 1:3; schpts{m,nn} = 1:3; end; end 
% 
% for nn = 1:3  
% subplot( 1, 3, nn); hold on; 
% boxplot( errmean{nn,1}(2:end), schpts{nn,1}(:)-0.25, 'Positions', schpts{nn,1}(:)-0.25, 'PlotStyle','compact' ); 
% boxplot( errmean{nn,3}(2:end), schpts{nn,3}(:)+0.25, 'Positions', schpts{nn,3}(:)+0.25, 'PlotStyle','compact' ); 
% boxplot( errmean{nn,2}(2:end), schpts{nn,2}(:), 'Positions', schpts{nn,2}(:), 'PlotStyle','compact' ); 
% % set(gca, 'Yscale', 'log' ); box on; xlabel('Scan number'); ylabel('error'); grid on; %ylim( [0.01 30] ); 
% % set(gca,'fontsize', 14 ); title( strcat('Group ', int2str(nn)) ); 
% end 

