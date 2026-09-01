function [tsol, ysol] = headnecktumor_RT( params, data )


tsol = data.xdata(1); % convert to days  
ysol = data.ydata(1); 
y = ysol; 

RTtime = [0:4, 7:11, 14:18, 21:25, 28:32, 35:39, 42:46]/7; 


if( tsol < 0 )
treat = [tsol RTtime(RTtime<data.xdata(end)) max(data.xdata(end), 6)]; 
else
treat = [RTtime(RTtime<data.xdata(end)) max(data.xdata(end), 6)]; 
end 
treat = unique( treat ); 


for nt = 1:(length(treat)-1)
  
    if( (treat(nt+1) - treat(nt)) ~= 0 )

            [t,y] = ode23(@(t,y)tumorLogistic(t,y,params),[treat(nt),treat(nt+1)],y(end));

             % K = V0 / params(2); 
             params(2) = params(2)*(1-params(3));  % K' = K (1-delta) so PSI' = PSI / (1-delta) 
            
    end
    ysol = [ysol; y(2:end)]; tsol = [tsol; t(2:end)];
    

end

end



function dval = tumorLogistic( t, val, params )

lambda = params(1);
K = params(2);

V = val(1);

dval(1) = lambda * V * ( 1 - V / K );

end

