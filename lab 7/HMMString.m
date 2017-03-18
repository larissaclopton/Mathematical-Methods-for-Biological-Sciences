function [hidden_str,obs_str] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len)

% a function to generate a string of observations from the HMM

% Inputs
% states - list of hidden states
% nucleotides - list of observed states
% trans_mat - transition matrix between hidden states
% obs_mat - observation (emission) probability matrix
% init_prob - initial probabilities for the first state
% str_len - length of string we wish to generate
 
    hidden_str = zeros(1,str_len);
    obs_str = zeros(1,str_len);
    init_hidden = randsample(states, 1, true, init_prob);
    
    % pick an observed state for the initial hidden state
    idx = find(states == init_hidden);
    prob = rand;
    for i = 1:length(nucleotides)
        if prob < sum(obs_mat(idx,1:i))
            init_obs = nucleotides(i);
            break
        end
    end
    
    hidden_str(1) = init_hidden; obs_str(1) = init_obs;
    
    % generate the rest of the hidden and observed states
    for i = 2:str_len
        
        % pick a hidden state given the previous hidden state
        idx = find(states == hidden_str(i-1));
        prob = rand;
        if prob < trans_mat(1,idx)
            hidden_str(i) = 'N';
        else
            hidden_str(i) = 'C';
        end
             
        % pick an observed state for the hidden state
        idx = find(states == hidden_str(i));
        prob = rand;
        for j = 1:length(nucleotides)
            if prob < sum(obs_mat(idx,1:j))
                obs_str(i) = char(nucleotides(j));
                break
            end
        end 
        
    end

    % convert result to a character string
    hidden_str = char(hidden_str); obs_str = char(obs_str);
    
end