function [tsol, ysol] = tumorModel_postRT_wIC(params,IC)
global dose 

%Specify initial condition
tsol = IC(1);
y = IC(2);
ysol = y;

treat = [IC(1) 15:19 22:26 29:33 36:40 43:47 50:54]; 
% treat = unique(treat); 
if( treat(1) == treat(2) )
    treat(2) = treat(2) + 0.1^4; 
end 

beta = params(end); %set alpha, estimate alpha/beta ratio
alpha = 1.5*beta;
d = dose; %dosage level
for nt = 1:(length(treat)-1)
  
    if( (treat(nt+1) - treat(nt)) ~= 0 )
            [t,y] = ode23(@(t,y)tumorLogistic(t,y,params),[treat(nt), treat(nt+1)],y(end));
    end
    ysol = [ysol; y(2:end)]; tsol = [tsol; t(2:end)];

    Radio = (1 - exp( (-alpha * d - beta * d^2 )))  * y(end);
    y(end) = y(end) - Radio;
    
end

end



function dval = tumorLogistic( t, val, params )

lambda = params(1);
K = params(2);

V = val(1);

dval(1) = lambda * V * ( 1 - V / K );

end

