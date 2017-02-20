function [counts] = GibbsSamplerMultiple(cond_prob,steps,states)

% Inputs
% cond_prob - conditional probabilities for being open, 
% give the number of open neigbors
% steps - the number of time steps
% random_states - matrix of randomly assigned states for all ion channels,
% to be modified as the algorithm progresses

    nchannels = size(states,1); % the number of ion channels

    for t = 2:steps % leave initial states as is
        
        for channel = 1:nchannels
            
            % special case for first ion channel
            if channel == 1
                % choose state based on that of following 
                % neighbor from previous time step
                neighbor_sum = 1 - states(channel+1,t-1);
            % special case for last ion channel
            elseif channel == nchannels
                % choose state based on that of previous 
                % neighbor from previous time step
                neighbor_sum = 1 - states(channel-1,t-1);
            % general case
            else
                % choose state based on two surrounding 
                % neighbors from previous time step
                neighbor_sum = 2 - (states(channel-1,t-1) + states(channel+1,t-1));
                
            end
            
            % decide the state of the channel based
            % on a randomly drawn probability
            prob = rand;
            if cond_prob(neighbor_sum+1) <= 0.5
                if prob < cond_prob(neighbor_sum+1)
                    states(channel,t) = 0;
                else
                    states(channel,t) = 1;
                end
            else
                if prob < 1 - cond_prob(neighbor_sum+1)
                    states(channel,t) = 1;
                else
                    states(channel,t) = 0;
                end
            end
                
            
        end
        
        
    end
    
    % retrieve the total number of open and closed channels at each time step 
    counts = zeros(2,steps); % first row total open, second row total closed
    for t = 1:steps
        counts(1,t) = length(find(states(:,t) == 0));
        counts(2,t) = length(find(states(:,t) == 1));
    end
    
end