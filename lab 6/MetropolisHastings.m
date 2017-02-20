function [states] = MetropolisHastings(p,q,x0,steps)

% Inputs
% p - target probabilities of each state
% q - proposal distribution of each state
% x0 - the initial state
% steps - the number of time steps

    states = zeros(1,steps); % states across the time steps
    x = x0; % set the initial state
    states(1) = x;
    
    for t = 2:steps
        
        % generate the candidate state
        % draw random probability from a uniform distribution
        prob = rand;
        if prob < q(x+1)
            candidate = ~x;
        else
            candidate = x;
        end
        
        % decide whether or not to accept the candidate state
        % by drawing another random number from a uniform distribution
        u = rand;
        if u <= min([1 p(candidate+1)*q(candidate+1)/(p(x+1)*q(x+1))])
            x = candidate; % update the state, otherwise stay the same
        end
        
        % record the state at this time step
        states(t) = x;
        
    end

end