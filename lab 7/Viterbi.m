function [hidden_states] = Viterbi(states,nucleotides,trans_mat,...
                                obs_mat,init_prob,observed)
           
% a function to calculate the most likely string of hidden states based on
% a given HMM and a string of observations
                            
% Inputs
% states - list of hidden states
% nucleotides - list of observed states
% trans_mat - transition matrix between hidden states
% obs_mat - observation (emission) probability matrix
% init_prob - probabilities for initial hidden state
% observed - string of observed states

    scores = zeros(length(states),length(observed)); % probability scores
    traceback = zeros(length(states),length(observed)); % traceback matrix
    
    % fill in the probability scores and traceback matrix
    idx = find(nucleotides == observed(1));
    scores(:,1) = log(init_prob') + log(obs_mat(:,idx));
    for i = 2:length(observed)
        for j = 1:length(states)
            idx = find(nucleotides == observed(i));
            [max_val,max_idx] = max(scores(:,i-1) + log(trans_mat(j,:))');
            scores(j,i) = log(obs_mat(j,idx)) + max_val;
            traceback(j,i) = max_idx;
        end
    end
    
    % determine the ending most probable path
    % and trace back the path for the hidden states
    z = zeros(1,length(observed)); 
    hidden_states = zeros(1,length(observed));
    [~,final_idx] = max(scores(:,length(observed)));
    z(end) = final_idx;
    hidden_states(end) = states(z(end)); 
    for i = length(hidden_states):-1:2
        z(i-1) = traceback(z(i),i);
        hidden_states(i-1) = states(z(i-1));
    end
    
    % convert result to a character string
    hidden_states = char(hidden_states);
    
end