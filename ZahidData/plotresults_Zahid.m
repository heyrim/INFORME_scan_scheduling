


%% load data 
load('/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Data/Data_Zahid.mat')
% load( 'patientID.mat' )

RGB = get(groot,"FactoryAxesColorOrder");

global growthParams V0 
growthParams = 0.33;        % new fit all together     %0.13;  % from Zahid's paper 

% sort patients 
for npatient = 1:length( data.xdata )
    V0 = data.ydata{npatient}( abs(data.xdata{npatient}) < 0.1 ); 
    reduction(npatient) = data.ydata{npatient}(end)/V0; 
end 
[ii,jj] = sort( reduction ); 


for n = 1:length( data.xdata )

    % np = jj(n);
    np = n; 
    load( [ './Figures/Zahiddata/Zahid_Results_nPatient', int2str(np), '.mat'])

    lastscan = min( length(data.xdata{np}) - 2 , 6 ); 
    fullscan(1:6) = NaN; 
    fullscan(1:lastscan) = 1; 

    adaptscan(1:6,n) = fullscan-1; adaptscan(1:6,n) = NaN; 

    errcheck = abs( delta_vec(1:end-1) - delta_vec(2:end) ); 
    relerrcheck = abs( delta_vec(1:end-1) - delta_vec(2:end) )./delta_vec(2:end); 
    ind1 = find( errcheck < 0.004, 1 ); 
    ind2 = find( relerrcheck < 0.15, 1 ); 
    if( isempty(ind1) ); ind1=100; end
    if( isempty(ind2) ); ind2=100; end
    if( (ind1+ind2)==200 ); ind1 = length(delta_vec)-1; end
    ind = min( ind1, ind2 );
    nscan(n) = ind+1; 

    point_list_loc = point_list(2:(2+ind),1:2); 

    for m = 1:size( point_list_loc )
        ind = find( data.xdata{np}(3:end) == point_list_loc(m,1) ); 
        if( ind > 6 );    ind = 6;    end 
        adaptscan(ind,n) = 1; 
    end 
end 

figure; heatmap( adaptscan(:,jj) ); colorbar off; 
% group 14, 14, 11 
clear ii 
ii{1} = 1:14;   ii{2} = 15:28;   ii{3} = 29:39;  
figure; for n = 1:3; subplot( 1, 3, n ); heatmap( adaptscan(:,jj(ii{n})) ); colorbar off; end 

adaptscantmp = adaptscan; adaptscantmp(isnan(adaptscantmp)) = 0; 
figure; 
for n = 1:3 
    mean( nscan( jj(ii{n})) )
    mean( adaptscantmp( :, jj(ii{n}) )' )
    % subplot( 1, 3, n ); heatmap( mean( adaptscantmp( :, jj(ii{n}) )' ) )
end 

% threshold 0.02 and 0.1 - 2 3 2 
    % 2.2857
    %      0    1.0000    0.5714    0.4286    0.0714    0.2143
    % 2.6429
    %      0    1.0000    0.5000    0.5000    0.4286    0.2143
    % 2.2727
    %      0    1.0000    0.5455    0.5455    0.0909    0.0909

% [0 1 1 0 0 0; 
%  0 1 1 1 0 0; 
%  0 1 1 0 0 0]; 

% threshold 0.01 and 0.1 - 3 3 3 
    % 2.5714
    %      0    1.0000    0.5714    0.4286    0.2857    0.2857
    % 2.8571
    %      0    1.0000    0.5000    0.5000    0.5714    0.2857
    % 2.6364
    %      0    1.0000    0.5455    0.6364    0.3636    0.0909
% [0 1 1 1 0 0; 
%  0 1 1 0 1 0; 
%  0 1 1 1 0 0]; 



% threshold 0.003 and 0.05 - 3 3 3 
    % 2.7857
    %      0    1.0000    0.5714    0.4286    0.2857    0.2857
    % 3.0714
    %      0    1.0000    0.5000    0.5714    0.6429    0.3571
    % 2.9091
    %      0    1.0000    0.5455    0.7273    0.4545    0.1818




    
path = 'Figures/251125_Zahid/';
filename = [path 'ModelFit_Zahid_Schemeall_errors.mat'];
load( filename )

err_vec2( err_vec2==0 ) = NaN;
figure; heatmap( log10(err_vec2') )
clim( [-3 -1] )
set(gca,'fontsize', 12 ); xlabel( 'patient number' ); ylabel( 'week' ); 

err_vec4 = zeros( 39, 6 ); err_vec4(:) = NaN; 
for npatient = 1:39 
    load( ['/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Figures/251221_Zahid_onept/Zahid_Results_nPatient' int2str(npatient) '.mat'])
    
    scanlist = round( point_list(2:end,1) ); 

    %%%% some scans that does not work with round. 
    if( npatient == 16 ) 
        scanlist = [2; 6]; 
    end
    if( npatient == 21 ) 
        scanlist = [2; 5]; 
    end    
    if( npatient == 23 ) 
        scanlist = [2; 3]; 
    end
    if( npatient == 27 ) 
        scanlist = [2; 3]; 
    end

    for nn = 1:length(scanlist)
        err_vec4( npatient, scanlist(nn) ) = err_vec3(npatient, nn+1); 
    end 
    scanlistall( npatient, scanlist ) = 1; 
    nscan(npatient) = length(scanlist); 
end 
figure; heatmap( log10(err_vec4') )
clim( [-3 -1] )
set(gca,'fontsize', 12 ); xlabel( 'patient number' ); ylabel( 'week' ); 

clear ii 
ii{1} = 1:14;   ii{2} = 15:28;   ii{3} = 29:39;  
figure; for n = 1:3; subplot( 1, 3, n ); heatmap( scanlistall(jj(ii{n}),:)' ); colorbar off; end 



% mean(scanlistall(jj(ii{1}),:))
% =         0    1.0000    0.5714    0.5000    0.1429    0.3571
% 
% mean(scanlistall(jj(ii{2}),:))
% =         0    1.0000    0.5714    0.6429    0.5000    0.3571
% 
% mean(scanlistall(jj(ii{3}),:))
% =         0    1.0000    0.5455    0.8182    0.3636    0.1818
% mean( nscan(jj(ii{1}) )) =    2.5714
% mean( nscan(jj(ii{2}) )) =    3.1429
% mean( nscan(jj(ii{3}) )) =    2.9091

priority 
[ 2, 3 ] 
[ 2, 3, 4 ]
[ 2, 4 ]

% err_Zahid' = 
%    1.5414e-01   5.7771e-02   2.9092e-02   1.1267e-02   9.4026e-03   7.2747e-03
%    5.4832e-02   2.1771e-02   1.4700e-02   1.8344e-02   1.1010e-02            0
%    5.5670e-02   2.2042e-02   1.3257e-02            0            0            0

% for n = 1:3 
% ii{n} = find( ngroup == n ); 
% mean( err_vec4( ii{n}, : ) ) 
% end 
%%%% group 
%    3.9259e-02   2.4915e-02   1.7659e-02
%    2.8216e-02   1.5603e-02   8.0721e-03
%    1.1075e-01   1.6186e-02   6.3453e-03





