function [T, Y, output] = multi_tissue_main_VEGF(c,p,m,y0,t_end,runVar)
switch runVar
    
    case 0 % Standard Simulation: from 0 to t_end
        options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-3);
        inits = y0;
        [T,Y] = ode15s(@multi_tissue_eqns_VEGF,0:60:t_end,inits,options,c,p,m); % return every minute
        output = 0;
        
    case 1 % Steady-State Simulation: runs until < 0.01% change in species (steady-state)
        options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-3);
        T=[];
        Y=[];
        
        inits = y0;
        [t,y] = ode15s(@multi_tissue_eqns_VEGF,0:60:3600,inits,options,c,p,m); % return every 60s 
        T = [T;t(2:end,:)];
        Y = [Y;y(2:end,:)];
        frac_change = abs(y(end,:) - y(1,:)) ./ y(1,:);
        
        hour = 1;
        while max(frac_change >= .001) 
            inits = y(end,:);
            [t,y] = ode15s(@multi_tissue_eqns_VEGF,hour*3600:60:(hour+1)*3600,inits,options,c,p,m);
            T = [T;t(2:end,:)];
            Y = [Y;y(2:end,:)];
            frac_change = y(1,:)*0;            
            for i=1:length(frac_change)
                if y(1,i) > 1e-6 % ignore tiny numbers - Matlab issue
                    frac_change(i) = abs(y(end,i) - y(1,i)) ./ y(1,i);
                end
            end
            %disp(perc_change)
            hour = hour + 1
            [maxch,index]= max(frac_change)
        end
        % disp(perc_change) % able to print %change if needed
        output = hour;
       
end
