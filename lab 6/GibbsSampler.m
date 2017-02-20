function [counts] = GibbsSampler(cond_prob,steps)

% Inputs
% cond_prob - conditional probabilities for the random variables, X then Y
% steps - the number of time steps

    states = zeros(2,steps); % hold the state values for the two correlated ion channels

    % pick x0
    x0 = 0; % start X in open state
    % choose y0 from conditional distribution P (Y|X = X0)
    prob = rand;
    if cond_prob(2,1) < cond_prob(2,2)
        if prob < cond_prob(2,1)
            y0 = 0;
        else
            y0 = 1;
        end
    else
        if prob < cond_prob(2,2)
            y0 = 1;
        else
            y0 = 0;
        end
    end
    
    % store the initial states
    x = x0; y = y0;
    states(1,1) = x; states(2,1) = y;
    
    for t = 1:steps-1
        % pick x(t+1) from P(X|Y=y(t))
        prob = rand;
        if cond_prob(1,2*y+1) < cond_prob(1,2*y+2)
            if prob < cond_prob(1,2*y+1)
                x = 0;
            else
                x = 1;
            end
        else
           if prob < cond_prob(1,2*y+2)
               x = 1;
           else
               x = 0;
           end
        end
        
        % pick y(t+1) from P(Y|X=x(t+1))
        prob = rand;
        if cond_prob(2,2*x+1) < cond_prob(2,2*x+2)
            if prob < cond_prob(2,2*x+1)
                y = 0;
            else
                y = 1;
            end
        else
            if prob < cond_prob(2,2*x+2)
                y = 1;
            else
                y = 0;
            end
        end
        
        states(1,t+1) = x; states(2,t+1) = y;
        
    end
    
    % retrieve the total number of open and closed channels at each time step 
    counts = zeros(2,steps); % first row total open, second row total closed
    for t = 1:steps
        counts(1,t) = length(find(states(:,t) == 0));
        counts(2,t) = length(find(states(:,t) == 1));
    end
    
    
end