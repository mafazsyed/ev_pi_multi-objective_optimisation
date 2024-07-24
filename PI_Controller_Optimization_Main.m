clear all;
clc;

%% Sampling Plans

%Initializing q for the sampling plans
q = [10, 10];

% Full Factorial Sampling Plan
fullFactorial_sp = fullfactorial(q, 1);
fullFactorial_sp = 10*fullFactorial_sp + eps;

% RLH Sampling Plan
rlh_sp = rlh(100, 2, 5);
rlh_sp = 10*rlh_sp + eps;

% Sobol Sampling Plan
P = sobolset(2);
sobol_sp = net(P,100);
sobol_sp = 10*sobol_sp + eps;

% Scatter plots of sampling plans
figure;
scatter(fullFactorial_sp(:, 1), fullFactorial_sp(:, 2))
xlabel('K_p')
ylabel('K_i')
% title('Full Factorial Sampling')

figure;
scatter(rlh_sp(:, 1), rlh_sp(:, 2))
xlabel('K_p')
ylabel('K_i')
% title('RLH Sampling ')

figure;
scatter(sobol_sp(:, 1), sobol_sp(:, 2))
xlabel('K_p')
ylabel('K_i')
% title('Sobol Sampling')

% Assessing sampling plans using mmphi metric
fullFactorial_metric = mmphi(fullFactorial_sp, 5, 2); % Full Factorial
rlh_metric = mmphi(rlh_sp, 5, 2); % RLH 
sobol_metric = mmphi(sobol_sp, 5, 2); % Sobol

%% Knowledge Discovery
% Evaluating full factorial design
Z = evaluateControlSystem(fullFactorial_sp);
% Y = [fullFactorial_sp, Z];

% Plotting the correlation of performance criteria and design variables using gplotmatrix
figure;
labels = ["Largest Pole CL", "GM", "PM", "RT", "PT", "OS_{max}", "US_{max}", "ST", "SS_{err}","ACE"];
gplotmatrix(Z, [], [], [], [], [], 'on', [], labels(1:10),labels(1:10));

% Plotting the correlation of performance criteria and design variables using parallelplot
figure;
parallelplot(Z);

%% Heat Map
% labels for the axes
labels = ["Largest Pole CL", "GM", "PM", "RT", "PT", "OS_{max}", "US_{max}", "ST", "SS_{err}","ACE"];

% Identify rows that have inf values
hasInf = any(isinf(Z), 2);

% Remove rows that contain any Inf values
Z_inf_removed = Z(~hasInf, :);

% Compute the correlation matrix
corrMatrix = corr(Z_inf_removed);

% Create the heatmap
figure;
h = heatmap(corrMatrix);

% Set the labels for the axes
h.XDisplayLabels = labels;
h.YDisplayLabels = labels;

% Display the correlation values inside the heatmap cells
h.CellLabelFormat = '%.2f';

%% Optimization Process
% Applying post-processing conditions of Gain Margin and Phase Margin before the optimization process
Z_eval = optimizeControlSystem(fullFactorial_sp);

% Initializing an empty vector for hypervolume indicator
HV_values = [];
% Reference points of hypervolume indicator
Reference_values = [max(Z_eval)];

% Loop for NSGA-II optimizer for 250 iterations
parent = fullFactorial_sp;
for i = 1:125
   disp('iteration: ')
   disp(i)
  % Priorities: Hard_constraint =3; High =2; Moderate = 1; Low = 0;
  % Progressive use of preferences over iterations
  if i <= 50
      Goals = [0.8 -inf -inf -inf -inf -inf -inf -inf -inf -inf];
      Priorities = [1 0 0 0 0 0 0 0 0 0];
  elseif i > 50 && i <= 75
      Goals = [0.8 -6 20 -inf -inf -inf -inf -inf -inf 0.67];
      Priorities = [2 1 1 0 0 0 0 0 0 1];
  elseif i > 75 && i <=100
  % Applying all the goals and priorities as per cheif engineer's preferences
      Goals = [0.8 -6 20 2 10 10 8 20 1 0.67];
      Priorities = [3 2 2 1 0 1 0 0 1 2];
     
  else
      Goals = [0.8 -6 20 2 10 10 8 20 1 0.67];
      Priorities = [3 2 2 1 0 1 0 0 1 1];
    
  end
  % Incorporating rank preferability
  [ranking,ClassV] = rank_prf(Z_eval,Goals, Priorities);
 
  % crowd sorting
  crowd = crowding(Z_eval,ranking);
  ranking = max(ranking) - ranking;
  % Binary tournament selection
  binary = btwr([ranking,crowd],100);
  % setting bounds for the cross and mutation
  bounds = [0 0; 1 1];
  % CrossOver and Mutation variation
 
  cross = sbx(parent(binary,:),bounds); % CrossOver
 
  % Post Variation Children
  child = polymut(cross,bounds); % Mutation
  
  % Concatinating Parent with child
  P_C = [parent;child];
 
  % Evaluating the concatenated Parent and child
  Z_eval_PC = optimizeControlSystem(P_C);
  
  % rank preferability for the concatenated Parent and child
  [ranking_PC,Class] = rank_prf(Z_eval_PC,Goals, Priorities);
 
  % crowd sorting for the concatenated Parent and child
  crowd_PC = crowding(Z_eval_PC,ranking_PC);

  % NSGA-II selection for survival operator
  reduceNSGA = reducerNSGA_II(P_C,ranking_PC,crowd_PC,100);

  % reduced designs of two populations
  reduced_design = P_C(reduceNSGA,:);

  % Evaluating the reduced designs of two populations
  Z_eval = optimizeControlSystem(reduced_design);

  % ..... Hypervolume ..... %
  % Indicator
  HV = Hypervolume_MEX(Z_eval,Reference_values);
  % Vector
  HV_values = [HV_values; HV];
  % ....................... %

  % Updating sampling plan for next iteration
  parent = reduced_design;

end

%% Plotting Optimised Results
figure;
scatter(parent(:, 1), parent(:, 2))
xlabel('K_p')
ylabel('K_i')
% title('Full Factorial Sampling After Optimisation')

% Plotting optimised correlation of performance criteria and design variables gplotmatrix)
figure;
labels = ["Largest Pole CL", "GM", "PM", "RT", "PT", "OS_{max}", "US_{max}", "ST", "SS_{err}","ACE"];
gplotmatrix(Z_eval, [], [], [], [], [], 'on', [], labels(1:10),labels(1:10));

% Plotting the optimized correlation of design variables and performance criteria using gplotmatrix
figure;
parallelplot(Z_eval);

%% Plotting Hypervolume Result
figure;
plot(HV_values);
xlabel('Iteration');
ylabel('Hypervolume');

%% Functions
function Z = optimizeControlSystem(x)
  Z = evaluateControlSystem(x);
  % Gain Margin adjustment
  Z(:, 2) = -(20 * log10(Z(:, 2)));
  % Phase Margin adjustment
  Z(:, 3) = abs(Z(:, 3) - 50);
end