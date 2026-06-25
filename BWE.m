%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Bézier walk evolution (BWE) source codes version 1.0
%  
%  Developed in:	MATLAB 23.2.0.2365128 (R2023b)
%  
%  Programmer:		Jinpeng Wang
%                   e-mail: wangjinpengchunuo@163.com &
%                   2211060117@stu.lntu.edu.cn
%  
%  Original paper:	Jinpeng Wang, Xingguo Xu, Yujing Sun, Jiguang Yu,
%                   Kaichen Ouyang, and Yuansheng Gao.
%                   Random Walk on Bézier Curves for
%                   Global Optimization
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Inputs:
%   fun              - Objective function handle
%   nvars            - Number of decision variables (Dimensions)
%   lb               - Lower bound of search space
%   ub               - Upper bound of search space
%   N                - Population size
%   T                - Maximum number of evaluations (iterations)
% Outputs:
%   TargetX          - Best solution found (Global Elite)
%   TargetF          - Best fitness value found
%   ConvergenceCurve - Historical best fitness at each evaluation
% -------------------------------------------------------------------------

function [TargetX, TargetF, ConvergenceCurve] = BWE(fun, nvars, lb, ub, N, T)
%% 1. Initialization Phase
% Uniformly initialize population within the specified bounds
Positions = initialization(N, nvars, ub, lb);
Fitness = zeros(N, 1);
for i = 1:N
    Fitness(i) = fun(Positions(i, :));
end

% Identify the current global best individual (Elite)
[Alpha_Fitness, alphaIdx] = min(Fitness);
Alpha_X = Positions(alphaIdx, :);
ConvergenceCurve = zeros(1, T);

% Core Parameter Settings
S = max(3, round(N * 0.2));          % Size of the random walk sample pool
triu_mask = triu(true(S, S), 1);     % Mask for internal sample distance calculations
[row_idx, col_idx] = find(triu_mask);
t = 1; 
alpha = 0.2;                   % Stochastic fluctuation control parameter
beta = 0.8;                     % Baseline factor for step parameter calculation

%% 2. Main Evolution Loop
while t < T
    % === Adaptive Strategy Probability Calculation ===
    vartheta = t / T;                % Iteration progress factor
    theta = (1 - vartheta)^2;        % Time-dependent decay factor
    gamma_param = beta - theta; % Baseline for movement aggressiveness
    
    % Strategy weights inspired by Bernstein polynomials
    w3 = (1 - vartheta);             % Weight for cubic path (high exploration)
    w2 = 3 * vartheta * (1 - vartheta); % Weight for quadratic path (balanced)
    w1 = vartheta;                   % Weight for linear path (high exploitation)
    
    % Softmax-like normalization for strategy selection probabilities
    strategy_probs = [w3, w2, w1]; 
    strategy_probs = max(strategy_probs, 0); 
    strategy_probs = strategy_probs / sum(strategy_probs);
    
    % === Random Walk Sampling ===
    % Form a sample pool from the current population to represent local structure
    sampleIdxs = randperm(N, S);
    pos_sample = Positions(sampleIdxs, :);
    
    % Map distances to selection probabilities for guidance nodes
    distPopToSample = Distance_Calc(Positions, pos_sample);
    p1_vec = Probability_Mapping(distPopToSample); % Individual-to-sample probability
    
    distMatrix = Distance_Calc(pos_sample, pos_sample);
    p2_vals = Probability_Mapping(distMatrix(triu_mask)); % Inter-sample probability
    
    distSampleToAlpha = Distance_Calc(pos_sample, Alpha_X);
    p3_vec = Probability_Mapping(distSampleToAlpha); % Sample-to-elite probability
    
    newPositions = Positions;
    
    % === Individual Evolution via Bézier Walk ===
    for i = 1:N
        currentPos = Positions(i, :);
        % Select evolution order via Roulette Wheel Selection
        order = RouletteWheelSelection(strategy_probs);
        % Adaptive step parameter with stochastic fluctuation
        tau_i = abs(gamma_param + alpha * rand());
        
        if order == 1 % 3rd-order Bézier Evolution
            % Joint selection of two guidance nodes C1 and C2
            prob_3rd = p1_vec(i, row_idx) .* p2_vals' .* p3_vec(col_idx)';
            pathIdx = RouletteWheelSelection(prob_3rd);
            s1 = row_idx(pathIdx); s2 = col_idx(pathIdx);
            P1 = pos_sample(s1, :); P2 = pos_sample(s2, :);
            
            % Generate new position along a cubic Bézier trajectory
            newX = (1-tau_i)^3 * currentPos + ...
                   3*(1-tau_i)^2 * tau_i * P1 + ...
                   3*(1-tau_i) * tau_i^2 * P2 + ...
                   tau_i^3 * Alpha_X;
            % Approximate path curvature for perturbation scaling
            kappa = norm(P1 - currentPos) + norm(P2 - P1) + norm(Alpha_X - P2);
            
        elseif order == 2 % 2nd-order Bézier Evolution
            % Selection of one guidance node C1
            prob_2nd = p1_vec(i, :) .* p3_vec';
            s1 = RouletteWheelSelection(prob_2nd);
            P1 = pos_sample(s1, :);
            
            % Generate new position along a quadratic Bézier trajectory
            newX = (1-tau_i)^2 * currentPos + ...
                   2*(1-tau_i) * tau_i * P1 + ...
                   tau_i^2 * Alpha_X;
            kappa = norm(P1 - currentPos) + norm(Alpha_X - P1);
            
        else % 1st-order Bézier Evolution
            % Linear interpolation towards the global elite
            newX = (1 - tau_i) * currentPos + tau_i * Alpha_X;
            kappa = norm(Alpha_X - currentPos);
        end
        
        % === Approximate Curvature Perturbation ===
        % Compute normalization factor based on path tortuosity
        rho_i = kappa / (norm(Alpha_X - currentPos) + eps);
        % Apply heavy-tailed perturbation using tangent mapping
        newX = newX + rho_i * theta * tan(randn(1, nvars));
        
        % Boundary handling: Ensure individuals stay within search space
        newPositions(i, :) = max(min(newX, ub), lb);
    end
    
    % === Greedy Selection and Global Best Update ===
    for i = 1:N
        fNew = fun(newPositions(i, :));
        
        % Greedy update: accept only if the new position is better
        if fNew < Fitness(i)
            Positions(i, :) = newPositions(i, :);
            Fitness(i) = fNew;
        end
        
        % Update the global elite solution
        if fNew < Alpha_Fitness
            Alpha_X = newPositions(i, :);
            Alpha_Fitness = fNew;
        end
    end
    ConvergenceCurve(t) = Alpha_Fitness;
    t=t+1;
end

% Final output assignment
TargetX = Alpha_X; 
TargetF = Alpha_Fitness;
end

%% --- Auxiliary Function Modules ---

function idx = RouletteWheelSelection(weights)
    % Standard roulette wheel selection mechanism
    weights = weights / (sum(weights) + eps);
    c = cumsum(weights);
    idx = find(c >= rand(), 1, 'first');
    if isempty(idx), idx = length(weights); end
end

function p = Probability_Mapping(dist)
    % Transform raw distances into selection probabilities via Softmax
    if size(dist, 1) > 1 && size(dist, 2) > 1
        % Matrix processing: row-wise normalization for individual guidance
        d_min = min(dist, [], 2);
        d_max = max(dist, [], 2);
        norm_d = (dist - d_min) ./ (d_max - d_min + eps);
        exp_d = exp(norm_d);
        p = exp_d ./ (sum(exp_d, 2) + eps);
    else
        % Vector processing: normalization for intra-sample relationships
        d_min = min(dist);
        d_max = max(dist);
        norm_d = (dist - d_min) ./ (d_max - d_min + eps);
        exp_d = exp(norm_d);
        p = exp_d ./ (sum(exp_d) + eps);
    end
end

function d = Distance_Calc(X, Y)
    % Compute Euclidean distances between two sets of points
    n = size(X, 1); m = size(Y, 1);
    d = zeros(n, m);
    for j = 1:m
        d(:, j) = sqrt(sum((X - Y(j, :)).^2, 2));
    end
end

function Positions = initialization(N, dim, ub, lb)
    % Initialize population using uniform random distribution
    if length(lb) == 1
        Positions = rand(N, dim) .* (ub - lb) + lb;
    else
        Positions = zeros(N, dim);
        for i = 1:dim
            Positions(:, i) = rand(N, 1) .* (ub(i) - lb(i)) + lb(i);
        end
    end
end