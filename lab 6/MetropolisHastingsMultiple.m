function [counts] = MetropolisHastingsMultiple(p,q,x0,steps,M)

% Inputs
% p - target probabilities of each state
% q - proposal distribution of each state
% x0 - the initial state
% steps - the number of time steps
% M - the number of independently opening/closing ion channels

    states = zeros(M,steps); % hold the state values for all M ion channels
    for i = 1:M
        
        % run the Metropolis-Hastings algorithm on the ith ion channel
        x = x0; % set the initial state
        states(i,1) = x;
    
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
            states(i,t) = x;
            
        end
    
    end

    % retrieve the total number of open and closed channels at each time step 
    counts = zeros(2,steps); % first row total open, second row total closed
    for t = 1:steps
        counts(1,t) = length(find(states(:,t) == 0));
        counts(2,t) = length(find(states(:,t) == 1));
    end
    
end