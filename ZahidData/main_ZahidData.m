function main_ZahidData 

load('/Users/cho/Library/CloudStorage/GoogleDrive-heyrimc@ucr.edu/My Drive/Research/240531_Patient_Treatment_Scheduling/Code/SimpleRadiotherapyModel/Data/Data_Zahid.mat')
load( 'patientID.mat' )


nallcompute = 0; 


%% Compute first data point to collect from group prior 
if( nallcompute )
    firstpoint = seqDesign_ZahidData_prior_firstdatapoint; 
else
    firstpoint = 2; 

end 



%%% Something like 
% %% Compute adaptive schedule for each patients 
% for npatient = 1:length( data.xdata )
%     scan_schedule_K_beta{npatient} = seqDesign_CAData_prior_all(npatient,firstpoint); 
% end





for n = 1:length( data.xdata )

    % np = jj(n);
    np = n; 
    load( [ './Figures/251221_Zahid_onept/Zahid_Results_nPatient', int2str(np), '.mat'])

    lastscan = min( length( data.xdata{np}) - 2 , 6 ); 
    fullscan(1:6) = NaN; 
    fullscan(1:lastscan) = 1; 

    adaptscan(1:6,n) = fullscan-1; 

    errcheck = abs( delta_vec(1:end-1) - delta_vec(2:end) ); 
    relerrcheck = abs( delta_vec(1:end-1) - delta_vec(2:end) )./delta_vec(2:end); 
    % ind = min(find( errcheck < 0.01, 1 ), find( relerrcheck < 0.1, 1 )); if( isempty(ind) ); ind = length(delta_vec)-1; end 
    ind1 = find( errcheck < 0.01, 1 ); 
    ind2 = find( relerrcheck < 0.1, 1 ); 
    if( isempty(ind1) ); ind1=100; end
    if( isempty(ind2) ); ind2=100; end
    if( (ind1+ind2)==200 ); ind1 = length(delta_vec)-1; end
    ind = min( ind1, ind2 );

    point_list_loc = point_list(2:(2+ind),1:2); 

    for m = 1:size( point_list_loc )
        ind = find( data.xdata{np}(3:end) == point_list_loc(m,1) ); 
        if( ind > 6 ); 
            disp( 'what' ); ind = 6; 
        end 
        adaptscan(ind,n) = 1; 
    end 


end 






end 
