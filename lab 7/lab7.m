%% Lab 7: Crouching islands, hidden Markov
%
% Larissa Clopton
%
% A Hidden Markov Model consists of
%
% # list of hidden states
% # list of observed states
% # transition matrix between hidden states
% # observation (emission) probability matrix
%
% The following script defines the set of states (non-CpG 'N' and CpG 'C')
% and observations (nucleotides), the transition probabilities between the
% hidden states and the probabilities of observing each nucleotide for the
% two hidden states. 

clear all; close all;
states = ['N', 'C']; % define the names of the states (Non and CpG regions)
trans_mat = [0.7, 0.1; 0.3, 0.9]; % define a 2 x 2 transition matrix between the hidden states
nucleotides = ['A', 'C', 'G', 'T'];  % define the alphabet of nucleotides
AT_obs_probs = [0.49, 0.01, 0.01, 0.49]; % set the probabilities of observations for the AT-rich state
GC_obs_probs = [0.01, 0.49, 0.49, 0.01]; % set the probabilities of observations for the CG-rich state
obs_mat = [AT_obs_probs; GC_obs_probs]; % create a 2 x 4 matrix of observation probabilities

% Grading: 
% Part 1: 1.1 10 pts (6 pts for code, 4 pts for testing), 1.2 4 pts (for
% testing)
% Part 2: 2.1 16 pts (12 pts for HMM code, 4 pts for testing), 2.2 6 pts (4
% for testing, 2 for answers), 2.3 6 pts (4 for testing, 2 for answers)
% Total: 42 pts
%% Part 1: generating a string of observations from HMM
%
% In this part you will write a function to generate a string of
% observations according to a given HMM. The function requires 6 inputs:
% the 4 listed above plus 5) vector of initial probabilities and 6) the
% length of the string. To generate each character in the string requires two steps:
%
% # pick a hidden state according to the transition probability, given the previous hidden state
% # pick an observed state according the observation probability for the selected hidden state
%
% The hidden state of the first character needs to be selected from an
% initial distribution, since there is no previous state. 
% To randomly choose the initial state from the given
% probability distribution you can use the randsample() function as
% follows: 
% 
% init_prob = [0.2, 0.8];
% states = ['A', 'B'];
% randsample(states, 1, true, init_prob)

% 1.1 Write a function that takes in the 6 inputs listed above, generates a
% string of hidden states and a string of observations and returns them
% both. Test this function by generating a random string of 10 letters with the
% transition and observation matrix given above. Repeat the generation
% once again and report whether the result is different (hint: it most
% likely should be!) 

str_len = 10;
init_prob = [0.2, 0.8];

% generate string of hidden states and observations from the HMM
disp('1.1 - These should be different.');
[hidden_str1,obs_str1] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len);
display(hidden_str1); display(obs_str1);
[hidden_str2,obs_str2] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len);
display(hidden_str2); display(obs_str2);

% The results are in fact different.

% 1.2 Change the transition matrix so that the transition
% probabilites are 1 and 0 (either that the states stay the same or they
% always change) and change the observation matrix in the same fashion and
% run the function again to generate a string of 10 letters. Repeat the
% generation once again and report whether the result is different (hint:
% it should be the same!) 

init_prob = [0, 1];
trans_mat = [1, 0; 0, 1];
AT_obs_probs = [1, 0, 0, 0];
GC_obs_probs = [0, 0, 1, 0];
obs_mat = [AT_obs_probs; GC_obs_probs];

% generate string of hidden states and observations from the HMM
disp('1.2 - These should be the same.');
[hidden_str1,obs_str1] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len);
display(hidden_str1); display(obs_str1);
[hidden_str2,obs_str2] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len);
display(hidden_str2); display(obs_str2);

% The results will be the same, provided the initial states are the same, which
% is why I also defined the initial probability to be 0 for one state and 1
% for the other. The results are in fact the same (i.e. not different).

%% Part 2: implementing and testing the Viterbi algorithm
%
% In this part you will implement the Viterbi algorithm, which calculates
% the most likely string of hidden states based on a given HMM and a string
% of observations. A pseudocode for the algorithm is available here:
% https://en.wikipedia.org/wiki/Viterbi_algorithm#Pseudocode 
%
% The function should have the following inputs: 1) list of hidden states,
% 2) list of observed states, 3) transition matrix between hidden states,
% 4) observation (emission) probability matrix, 5) the vector of initial
% probabilities and 6) string of observed states. The function should
% return the strings of most likely hidden states. 
%
% Comments: 
%
% # the pseudocode uses poor variable names, particularly as it uses T_1 for
% the matrix of probability scores, T_2 for the traceback matrix, and T for
% the length of the sequence. Please choose better variable names!
% # The subscripts in the pseudocode mean indices
% # For better accuracy *add up logs of probabilities instead of multiplying them*.
%
% At each step the algorithm chooses the maximum value of the score out of
% many possibilites (as many as there are states). To do this, I
% recommend that you use the max() function to select the maximum score,
% which also returns the index of the maximum value if you request two
% outputs, e.g.
%
% sample_vector = [10,4,6];
% [max_val, max_ind] = max(sample_vector);
%
% 2.1 test the function by using the same inputs (transition and
% observation matrices) you used in question 1.2.  Generate a
% new string of 30 observations, and then run the Viterbi algorithm on that
% string to reconstruct the string of hidden states. Compare the string of
% hidden states from the Viterbi algorithm with the actual string of hidden
% states. 

clear all; close all;

% redefine the parameters from 1.2
states = ['N', 'C'];
nucleotides = ['A', 'C', 'G', 'T'];
trans_mat = [1, 0; 0, 1];
AT_obs_probs = [1, 0, 0, 0];
GC_obs_probs = [0, 0, 1, 0];
obs_mat = [AT_obs_probs; GC_obs_probs];
str_len = 30;
init_prob = [0, 1];

disp('2.1');
% generate a string of observed states and their hidden states
[actual_hidden,observed] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len);
% obtain the most probable hidden states from the Viterbi algorithm
predicted_hidden = Viterbi(states,nucleotides,trans_mat,obs_mat,init_prob,observed);
% compare the two strings
display(actual_hidden); display(predicted_hidden);

% The Viterbi algorithm reconstructs the string of hidden states perfectly
% for this HMM.

% 2.2 change the transition and observation matrix back to the same
% parameter values provided in the code before part 1 of the
% assignment. Generate a new string of 30 observations, and then run the
% Viterbi algorithm to reconstruct the string of hidden states. Compare the
% string of hidden states from the Viterbi algorithm with the actual string
% of hidden states and report the fraction of correct hidden states.
% Produce another string of 30 observations with the same 
% HMM, and repeat your experiment. Comment on how reliably the Viterbi
% algorithm reconstructs the string of hidden states for this HMM. 

% use parameter values from 1.1
trans_mat = [0.7, 0.1; 0.3, 0.9];
AT_obs_probs = [0.49, 0.01, 0.01, 0.49];
GC_obs_probs = [0.01, 0.49, 0.49, 0.01];
obs_mat = [AT_obs_probs; GC_obs_probs];
init_prob = [0.2, 0.8];

disp('2.2');
for i = 1:2

    disp(['Trial ' num2str(i)]);
    % generate a string of observed states and their hidden states
    [actual_hidden,observed] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len);
    % obtain the most probable hidden states from the Viterbi algorithm
    predicted_hidden = Viterbi(states,nucleotides,trans_mat,obs_mat,init_prob,observed);
    % compare the two strings
    display(actual_hidden); display(predicted_hidden);
    fraction_correct = (1/str_len)*sum(actual_hidden(:) == predicted_hidden(:));
    display(fraction_correct);

end

% For both trials, the string of hidden states returned by the Viterbi
% algorithm very well matches (but not perfectly as before) the actual string
% of hidden states. As such, the Viterbi algorithm fairly reliably reconstructs
% the string of hidden states for this HMM, resulting in a fraction of correct
% hidden states that is close to 1.

% 2.3 Change the matrix of observation probabilities so that for state 'N'
% the emission probabilities of 'A' and 'T' are both 0.4, and the emission
% probabilities of 'G' and 'C' are both 0.1; for state 'C'  the emission
% probabilities of 'G' and 'C' are both 0.4, and the emission probabilities
% of 'A' and 'T' are both 0.1. Generate a new string of 30 observations,
% and then run the Viterbi algorithm to reconstruct the string of hidden
% states and report the fraction of correct hidden states. Produce another
% string of 30 observations with the same HMM, and repeat your experiment.
% Comment on how reliably the Viterbi algorithm reconstructs the string of
% hidden states for this HMM. Explain why the fraction of correct hidden
% states is different in the two cases.

% change the emission probabilities
AT_obs_probs = [0.4, 0.1, 0.1, 0.4];
GC_obs_probs = [0.1, 0.4, 0.4, 0.1];
obs_mat = [AT_obs_probs; GC_obs_probs];

disp('2.3');
for i = 1:2

    disp(['Trial ' num2str(i)]);
    % generate a string of observed states and their hidden states
    [actual_hidden,observed] = HMMString(states,nucleotides,trans_mat,obs_mat,init_prob,str_len);
    % obtain the most probable hidden states from the Viterbi algorithm
    predicted_hidden = Viterbi(states,nucleotides,trans_mat,obs_mat,init_prob,observed);
    % compare the two strings
    display(actual_hidden); display(predicted_hidden);
    fraction_correct = (1/str_len)*sum(actual_hidden(:) == predicted_hidden(:));
    display(fraction_correct);

end
 
% After running the code multiple times, the Viterbi algorithm does not as 
% reliably reconstruct the string of hidden states for this HMM in comparison
% to the one in 2.2. It still gets a relatively high fraction correct, but
% the fraction is lower than that for the previous HMM in 2.2. Granted,
% there is a range of performance depending on the observed string, so how
% close performance is for this HMM with the one in 2.2 can vary across
% runs of the code. The fraction of correct hidden states is different 
% between the HMMs in 2.2 and 2.3 because the emission probabilities of a
% given nucleotide are more balanced across states. For example, the ratio
% of emission probabilities of nucleotides C and G in the CG-rich state to
% C and G in the AT-rich state is lower (or closer to 1) for the HMM in 2.3 
% than that in 2.2. At the same time, the ratio of emission probabilities of
% nucleotides A and T in the AT-rich state to A and T in the CG-rich state 
% is lower (or closer to 1) for the HMM in 2.3 than that in 2.2. This means
% in comparing the HMM in 2.3 to that in 2.2, the algorithm is more likely 
% to guess that a given nucleotide came from the incorrect state, resulting
% in a lower fraction of correct hidden states.