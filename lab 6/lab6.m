%% Lab 6: Markov Chain Monte Carlo simulations
% Larissa Clopton
% Metropolis-Hastings Algorithm
%
% Points breakdown:
%
% 1.1 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 1.2 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 1.3 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 1.4 7 pts (3 pts for code, 1 pt for each plot, 1 pt for answer)
% 2.1 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 2.2 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 2.3 2 pts (1 pt for plot, 1 pt for answer)
% 3.1 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 4.1 6 pts (4 pts for code, 1 pt for plot, 1 pt for answer)
% 4.2 6 pts (4 pts for code, 1 pt for plot, 1 pt for answer)
% Total: 56 points

% # Initialize the first state to X0 and set t = 0
% # Set the target probabilities P(X) (for all states X) and proposal
% probabilities q(i,j) of generating candidate state i from current state j.
% # Generate a candidate state Y using q(i,j) by drawing a random number u 
% from a uniform distribution on (0,1)
% # Draw another random number v; 
%   if v <= min(1,(P(Y)*q(X,Y))/(P(X)*q(Y,X))) then set X = Y 
%   else keep the same state
% # Set t=t+1 and repeat steps 3 through 5
%
%% Part 1: Single 2-state ion channel
% Simulate an ion channel switching between open and closed states, 
% with a number (e.g. 0) denoting the open state and another (e.g. 1) 
% denoting the closed state. Make the target fractions of open
% and closed states a changeable parameter, making sure they add up to 1. 
%
% 1.1 Use the simple proposal distribution: q(i,j) = 0.5 regardless of the states 
% i and j. Start with an open state. Run your simulation for 1000 time 
% steps using the target distribution: P(O) = 0.5; P(C) = 0.5. Plot the 
% states of the channel (0 and 1) over 1000 time steps. 

clear all; close all;

p_o = 0.5; % target open probability
p_c = 0.5; % target closed probability
q_co = 0.5; % proposal distribution, from open to closed
q_oc = 0.5; % proposal distribution, from closed to open
x0 = 0; % initial state, 0 for open and 1 for closed
steps = 1000; % the number of time steps

figure;
states = MetropolisHastings([p_o p_c], [q_co q_oc],x0,steps);
plot((1:1000)-1,states);
xlabel('Time step'); ylabel('State (0 open, 1 closed)');
title('P(O)=P(C)=0.5, q(C,O)=q(O,C)=0.5');

% 1.2 Repeat the simulation with the target distribution of P(O) = 0.1;
% P(C) = 0.9, and plot the states of the channel over 1000 time steps.

% update the target probabilities
p_o = 0.1;
p_c = 0.9;

figure;
states = MetropolisHastings([p_o p_c], [q_co q_oc],x0,steps);
plot((1:1000)-1,states);
xlabel('Time step'); ylabel('State (0 open, 1 closed)');
title('P(O)=0.1, P(C)=0.9, q(C,O)=q(O,C)=0.5');

% 1.3 Change the proposal distribution so that q(O,C) = 0.99 and q(O,O) = 0.99 
% (i.e., the probability of proposing the open state is always 0.99) and 
% report if you see a difference in the simulation for both of the two
% target distributions.

% update the proposal distribution
q_oc = 0.99;
q_co = 1 - 0.99; % as q(O,O) = 0.99

% run with the first target distribution
p_o = 0.5;
p_c = 0.5;

figure;
states = MetropolisHastings([p_o p_c], [q_co q_oc],x0,steps);
plot((1:1000)-1,states);
xlabel('Time step'); ylabel('State (0 open, 1 closed)');
title('P(O)=P(C)=0.5, q(C,O)=0.01, q(O,C)=0.99');

% run with the second target distribution
p_o = 0.1;
p_c = 0.9;

figure;
states = MetropolisHastings([p_o p_c], [q_co q_oc],x0,steps);
plot((1:1000)-1,states);
xlabel('Time step'); ylabel('State (0 open, 1 closed)');
title('P(O)=0.1, P(C)=0.9, q(C,O)=0.01, q(O,C)=0.99');

% There is certainly a difference in the simulation for the two target
% distributions given this proposal distribution. First, comparing the
% first and third figure with the same target distribution but different
% proposal distributions, we see that increasing the probability of proposing
% the open state to 0.99 reduces the amount of switches between open and
% closed states, though the overall proportion of time spent in the open and
% closed states is the same due to the same target distribution. Next,
% comparing the third and fourth figure with the same proposal distribution
% but two target distributions (the topic for this question), we see that
% very few switches are needed to maintain a target distribution of P(O)=0.1
% and P(C)=0.9, whereas more switches are needed to maintain a target
% distribution of P(O)=P(C)=0.5. 

% 1.4 Using P(O)=0.2, q(i,j) = 0.5, repeat the simulation N times (for N 
% at least 100) and plot the histograms of the states (open and closed) at 
% a few of the 1000 time steps (e.g. every 100). Report whether you 
% see convergence to the target distribution and approximately how quickly 
% it happens - this is called the burn-in time of the Markov chain simulation.

% update the target distribution
p_o = 0.2;
p_c = 1 - 0.2; % as P(O) is 0.2

% update the proposal distribution
q_co = 0.5;
q_oc = 0.5;

% run the simulation 200 times
ntrials = 100;
results = zeros(ntrials,steps);
for trial = 1:ntrials
    results(trial,:) = MetropolisHastings([p_o p_c], [q_co q_oc],x0,steps);
end

% plot histograms of the states at every 100 time steps
sample_times = 1:100:steps;
nsamples = length(sample_times);
figure;
for sample = 1:nsamples
    subplot(2,5,sample);
    histogram(results(:,sample_times(sample)),'Normalization','probability');
    axis([-0.5 1.5 0 1]);
    xlabel('State'); ylabel('Frequency');
    title(['t = ' num2str(sample_times(sample)-1)]);
end

% Yes there is convergence to the target distribution, and this happens
% after about 100 time steps. After this, the proportions stay close to the
% target distribution, even if it isn't exact.

%% Part 2: Multiple ion channels
% For this section, always use the parameters: P(O) = 0.2 and q(i,j) = 0.5 
%
% Modify your code to increase the number of ion channels to a changeable 
% parameter M, with each channel opening and closing independently (this 
% only requires adding a loop). Save the total number of open and closed 
% channels at each time step.

clear all; close all;

% define the target and proposal distribution
p_o = 0.2;
p_c = 1 - 0.2;
q_co = 0.5;
q_oc = 0.5;

% define the initial state and number of steps
x0 = 0;
steps = 1000;

% 2.1 Run the simulation at several values of M: 2, 10, 100 for 1000 time
% steps and report the mean number of open channels and the variance for
% each of the three cases. Is this what you expected? 

M = [2 10 100];

for i = 1:length(M)
    
    % retrieve the total number of open and closed channels at each time step
    counts = MetropolisHastingsMultiple([p_o p_c], [q_co q_oc],x0,steps,M(i));
    
    % compute the mean number of open channels
    disp([num2str(M(i)) ' ion channels']);
    mean_open = mean(counts(1,:));
    display(mean_open);
    
    % compute the variance, same value for open and closed
    var_open = var(counts(1,:));
    display(var_open);
    
end

% As expected, the variance increases with more ion channels. This makes
% sense because as more stochastic processes (in this case the opening and
% closing of an ion channel) are chained together the variability in the 
% output will increase. In a similar manner, the mean number of open
% channels tends to stray further from the expected number of open channels
% based on the target distribution, likely due to the greater variability
% of response.

% 2.2 Repeat the multiple channel simulation N times (for N at least 100) and
% plot the histograms of the number of channels at several different time
% points (e.g. plot them at t=1, t=100, t=200, etc.)
% Do they resemble any probability distribution you know?

sample_times = 1:100:steps;
nsamples = length(sample_times);
ntrials = 100;

% save the counts distribution (row 1 open and row 2 closed) 
counts_distribution_select = zeros(2,nsamples,ntrials); % for every 100th time point
counts_distribution_100 = zeros(2,100,ntrials); % for first 100 time points, for 2.3 below
for trial = 1:ntrials
    
    % 100 channel simulation
    counts = MetropolisHastingsMultiple([p_o p_c], [q_co q_oc],x0,steps,100);
    
    % save counts for first 100 time points, 
    % to be used in the following question
    counts_distribution_100(1,:,trial) = counts(1,1:100);
    counts_distribution_100(2,:,trial) = counts(2,1:100);
    
    % save counts for every 100th time point
    for n = 1:nsamples
        counts_distribution_select(1,n,trial) = counts(1,sample_times(n));
        counts_distribution_select(2,n,trial) = counts(2,sample_times(n));
    end
    
end
 
% plot histograms for the number of channels at each time point
for n = 1:nsamples
    figure(1);
    subplot(2,5,n);
    histogram(counts_distribution_select(1,n,:),'Normalization','probability');
    xlabel('# channels'); ylabel('Frequency');
    if n == 1
        axis([99.5 101.5 0 1]);
    else
        axis([0 40 0 0.2]);
    end
    title(['Open, t=' num2str(sample_times(n)-1)]);
    
    figure(2);
    subplot(2,5,n);
    histogram(counts_distribution_select(2,n,:),'Normalization','probability');
    xlabel('# channels'); ylabel('Frequency');
    if n == 1
        axis([-0.5 1.5 0 1]);
    else
        axis([60 100 0 0.2]);
    end
    title(['Closed, t=' num2str(sample_times(n)-1)]);
end

% The histograms for the number of open channels at various time points are
% given in the first figure, and the histograms for the number of closed
% channels at various time points are given in the second figure. One can
% see that the probability distributions are approximately normal.

% 2.3 We would like to find the burn-in time for converging to the
% stationary distribution. To do this, plot the
% *average* number of open/closed channels over the first 100 or so time
% steps. How quickly does this number reach an equilibrium value?

steps = 100; % look at the first 100 time steps
average_open = zeros(1,steps);
average_closed = zeros(1,steps);

for t = 1:steps
    % calculate the average number of open and closed channels
    average_open(t) = mean(counts_distribution_100(1,t,:));
    average_closed(t) = mean(counts_distribution_100(2,t,:));
end

figure(3);
plot(average_open);
xlabel('Time step'); ylabel('Average open channels');
title('Average open channels versus time');
figure(4);
plot(average_closed);
xlabel('Time step'); ylabel('Average closed channels');
title('Average closed channels versus time');

% The equilibrium values of 20 for open and 80 for closed are reached after 
% about 10 time steps.

%% Part 3: 2 correlated ion channels
% Now let us investigate the situation when the probability of opening and 
% closing of a channel depends on the state of its neighbor. Let the 
% target joint probability distribution be:
%   P(O,O) = 0.6; P(C,C) = 0.2; P(O,C) = P(C,O) = 0.1
%
% Calculate the conditional probabilities for states of
% one channel given the state of its neighbor based on the joint distribution
% and use the Gibbs algorithm to compute a random sample. The Gibbs
% sampler does not accept or reject proposed states, but instead uses
% conditional probabilities to generate new states of the two random
% variables X and Y (the two ion channels): 
%
% Gibbs Sampler
%
% # Pick X0 and choose Y0 from the conditional distribution P (Y|X = X0) 
% # Choose X1 from P(X|Y = Y0) and Y1 from P(Y|X = X1)
% # Repeat N times
% 
% 3.1 Run for 1000 steps and plot the histograms for number of open channels 
% for the first 100 time steps and for the second 100 time steps. How
% different are they? Plot the histogram for 201 through 1000 time steps,
% and compare with the first two. Approximately how long does it take for
% the histrograms to converge?

clear all; close all;

% define the target joint probability distribution
p_oo = 0.6;
p_cc = 0.2;
p_oc = 0.1;
p_co = 0.1;

% Pr(X=O|Y=O) = Pr(X=O^Y=O)/Pr(Y=O)
p_o_o = p_oo/(p_oo+p_co);
% Pr(X=C|Y=O) = Pr(X=C^Y=O)/Pr(Y=O)
p_c_o = p_co/(p_oo+p_co);
% Pr(X=O|Y=C) = Pr(X=O^Y=C)/Pr(Y=C)
p_o_c = p_oc/(p_cc+p_oc);
% Pr(X=C|Y=C) = Pr(X=C^Y=C)/Pr(Y=C)
p_c_c = p_cc/(p_cc+p_oc);

% the conditional probabilities for X
cond_x = [p_o_o p_c_o p_o_c p_c_c];

% Pr(Y=O|X=O) = Pr(Y=O^X=O)/Pr(X=O)
p_o_o = p_oo/(p_oo+p_oc);
% Pr(Y=C|X=O) = Pr(Y=C^X=O)/Pr(X=O)
p_c_o = p_oc/(p_oo+p_oc);
% Pr(Y=O|X=C) = Pr(Y=O^X=C)/Pr(X=C)
p_o_c = p_co/(p_cc+p_co);
% Pr(Y=C|X=C) = Pr(Y=C^X=C)/Pr(X=C)
p_c_c = p_cc/(p_cc+p_co);

% the conditional probabilities for Y
cond_y = [p_o_o p_c_o p_o_c p_c_c];

% obtain the counts of open/closed channels at each time point
steps = 1000;
counts = GibbsSampler([cond_x;cond_y],steps);

% plot histograms of open channels for 
% first 100 and second 100 time steps
figure;
histogram(counts(1,1:100),'Normalization','probability');
axis([-0.5 2.5 0 1]);
xlabel('Number of open channels'); ylabel('Frequency');
title('Open channel histogram, t = 1:100');
figure;
histogram(counts(1,101:200),'Normalization','probability');
axis([-0.5 2.5 0 1]);
xlabel('Number of open channels'); ylabel('Frequency');
title('Open channel histogram, t = 101:200');

% plot histogram for the rest of the time steps
figure;
histogram(counts(1,201:end),'Normalization','probability');
axis([-0.5 2.5 0 1]);
xlabel('Number of open channels'); ylabel('Frequency');
title('Open channel histogram, t = 201:1000');

% The probability of 2 ion channels being open should converge to P(O,O)=0.6, 
% the probability of 1 ion channel being open should converge to 
% P(O,C)+P(C,O)=0.1+0.1=0.2, and the probability of 0 ion channels being open
% should converge to P(C,C)=0.2. The 1-100 time step histogram is somewhat around
% these values but not near convergence. On the other hand, the
% 101-200 time step histogram is much closer and appears to have converged
% to these values. As such, the 201-1000 time step histogram more closely
% represents the the 101-200 time step histogram. From this, one can say it
% takes around 200 time steps for these histograms to converge.

%% Part 4: Multiple correlated ion channels
%
% Now we will use the Gibbs Sampler to simulate ion channels switching in
% the case that there are multiple correlated ion channels. The state of
% each ion channel depends on its right-hand and left-hand nearest neighbor: 
%
% P(O|2 open neighbors) = 0.4; P(O|1 open neighbor) = 0.2; P(O|0 open neighbors) = 0.1
% For the two ends of the ion channel array, the state of each ion channel
% can only depend on one nearest neighbor: 
% P (X = O|1 open neighbor) = 0.2 P (X = O|0 open neighbors) = 0.1
% The other probabilities are complementary.
%
% Gibbs Sampler with Multiple Correlated Channels
%
% # Generate an array of 100 ion channels in randomly-assigned states
% (0 or 1). 
% # Choose the state of the first ion channel in the array X(t,1) based on
% the state of its one neighbor X(t?1,2) from the previous time step. 
% # Choose the state of the second ion channel in the array X(t,2) based on 
% its neighbors X(t,1) and X(t?1,3), and do this for channels 2 through 99.
% # Select X(t,100) based on the state of its one neighbor X(t,99). 
% # Repeat N times.
%
% 4.1 Run for 10000 steps and plot the histogram for the number of channels 
% being open over all except for the first 100 time steps. 
% Approximately what number of open channels (total ion current) is
% expected? 

clear all; close all;

% define the conditional probabilities
p_o_0 = 0.1; % P(O|0 open neighbors)
p_o_1 = 0.2; % P(O|1 open neighbor)
p_o_2 = 0.4; % P(O|2 open neighbors)
cond_prob = [p_o_0 p_o_1 p_o_2];

% generate array of 100 ion channels in randomly-assigned states (0 or 1)
steps = 10000;
states = randi(2,100,steps) - 1;

% retrieve the counts at each time step from the Gibbs sampler
counts = GibbsSamplerMultiple(cond_prob,steps,states);

figure;
histogram(counts(1,101:end),'Normalization','probability');
axis([0 70 0 0.14]);
xlabel('Number of open channels'); ylabel('Frequency');
title('Open channel histogram, t=101:end');

expected_open1 = mean(counts(1,:));
display(expected_open1);

% Taking the mean count of the number of open channels across all time
% steps (keeping in mind that this number will vary each run due to 
% randomness), the expected number of open ion channels is approximately in
% the range of 12.5 to 13 ion channels.

% 4.2 Change your conditional probabilities to 
% P(O|2 open neighbors) = 0.2; P(O|1 open neighbor) = 0.5; 
% P(O|0 open neighbors) = 0.7 and re-run the simulation. Compare this 
% histogram to the one you obtained previously and estimate the expected
% current through the membrane.  

% update the conditional probabilities
p_o_0 = 0.7;
p_o_1 = 0.5; 
p_o_2 = 0.2; 
cond_prob = [p_o_0 p_o_1 p_o_2];

% redefine the randomly-assigned states
states = randi(2,100,steps) - 1;

% retrieve the counts at each time step from the Gibbs sampler
counts = GibbsSamplerMultiple(cond_prob,steps,states);

figure;
histogram(counts(1,101:end),'Normalization','probability');
axis([0 70 0 0.14]);
xlabel('Number of open channels'); ylabel('Frequency');
title('Open channel histogram, t=101:end');

expected_open2 = mean(counts(1,:));
display(expected_open2);

% This new histogram not only has a higher mean value but also a larger 
% spread (i.e. more variance and a greater standard deviation). It makes
% sense to have a higher mean value as we have greatly increased some of
% the conditional probabilities. The expected number of open ion channels 
% (total ion current) is around 48.5.