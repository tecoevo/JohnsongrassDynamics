function [P, R, SB, P_dens] = Dynamics_TwoHerbicides(RR1, RW1, ...
    RR2, RW2, herb1, herb2, tillage, seeds, bank, n_years, ...
    cost_seeds1, cost_seeds2, k_cost1, k_cost2, k_herb1, k_herb2)
% Dynamics_TwoHerbicides Gives the deterministic dynamics of a Johnsongrass
% population over a given number of years depending on initial frequencies  
% of resistant types, herbicde application, tillage, seed production,
% seed bank formation, fitness cost on seed production and dominance of 
% resistance and fitness cost.
%
%   Input: 
%   RR(1/2): initial fraction of the RR type in seeds and rhizomes
%   RW(1/2): initial fraction of the RW type in seeds and rhizomes
%   herb1: 1 x n_years vektor of herbicide application. Each entry 
%   corresponds to one season and is a logical value stating whether 
%   ACCase-inhibitor is applied
%   herb2: 1 x n_years vektor of herbicide application. Each entry 
%   corresponds to one season and is a logical value stating whether 
%   ALS-inhibitor is applied
%   tillage: 1 x (n_years+1) vektor of tillage strategy. Each entry 
%   corresponds to one season and is a logical value stating whether
%   the soil is tilled at season start
%   seeds: logical value stating whether seed production is considered
%   bank: logical value stating whether a seed bank is considered
%   n_years: number of years
%   cost_seeds(1/2): fitness cost on seed production associated with 
%   resistance 
%   k_cost(1/2): factor reducing the fitness cost of RS type relative to
%   RR type
%   k_herb(1/2): factor reducing the herbicide efficiency of RS type 
%   relative to SS type
%
%   Output:
%   P: matrix of absolute genotype frequencies in plants
%   R: matrix of absolute genotype frequencies in rhizomes
%   SB: matrix of absolute genotype frequencies in seed bank
%   P_dens: vector of plant densities

%% Simulation:
% Setting parameters:

% Simulation:
% Field size: 
A = 10000;

% Ecological:
% Proportion of selfpollination: 
p_self = 0.95;
% Proportion of seedling germination: 
p_germ = 0.26;

% Fecundity, i.e. number of seeds produced by one plant:
f = 5; 

% Avarage number of secondary rhizomes: 
lambda_sec = 2;
% Avarage number of nodes: 
lambda_nodes = 3;

% Rhizome winter mortality (south, north)(tillage; no tillage): 
p_mort = [0.5, 0.4; 0.25, 0.1];

% Natural seedling mortality:
d_seedlings = 0;
% Natural tiller mortality:
d_tillers = 0;
% Natural seed mortality in the seedbank:
d_seeds = 0.2;

% Highest possible plant density:
dens_max = 250;

% Evolutionary:
% Mutation rate:
mu = 10^(-8);


% Antropogenic:
% Herbicide efficacy herbicide 1 (ACCase-inhibitor): 
% Seedlings/Tillers for early POST: 
E_early1 = 0.95;
% Seedlings/Tillers for late POST: 
E_late1 = 0.9;
% Rhizomes (tillage, no tillage): 
E_rhizomes1 = [0.5, 0.25];
% Herbicide efficacy herbicide 1 (ACCase-inhibitor): 
% Seedlings/Tillers for early POST: 
E_early2 = 0.95;
% Seedlings/Tillers for late POST: 
E_late2 = 0.9;
% Rhizomes (tillage, no tillage): 
E_rhizomes2 = [0.5, 0.25];

% Initialization:
% Initial seedbank density: 
dens_seeds = 10;
% Initial rhizome density: 
dens_rhizomes = 1;


% Parameters to choose:
% Region with 1 corresponding to south and 2 corresponding to north:
region = 1;
% Adaption of praramters according to the region:
% a = a(region);
p_mort = p_mort(:, region);


% Cell array of inheritance matrices for one gene. With the row j column 
% k entry of the cell i matrix giving the fraction of type i seeds 
% produced by a type j plant pollinated by type k pollen. 
% (1 corresponding to genotype SS, 2 corresponding to genotype RS, 
% 3 corresponding to genotype RR) 
MI = {[1 0.5 0; 0.5 0.25 0; 0 0 0], [0 0.5 1; 0.5 0.5 0.5; 1 0.5 0],...
    [0 0 0; 0 0.25 0.5; 0 0.5 1]};
% Add Mutation:
M1{1} = (1 - mu)^2 * MI{1};
M1{2} = 2 * mu * (1 - mu) * MI{1} + (1 - mu) * MI{2};
M1{3} = mu^2 * MI{1} + mu * MI{2} + MI{3};

% Cell array of inheritance matrices. With the row j column k entry of the
% cell i matrix giving the fraction of type i seeds produced by a type j
% plant pollinated by type k pollen. (1 corresponding to genotype S1S1 S2S2, 
% 2 corresponding to genotype R1S1 S2S2, 3 corresponding to genotype R1R1 S2S2,
% 4 corresponding to genotype S1S1 R2S2, 5 corresponding to genotype R1S1 R2S2,
% 6 corresponding to genotype R1R1 R2S2, 7 corresponding to genotype S1S1 R2R2,  
% 8 corresponding to genotype R1S1 R2R2, 9 corresponding to genotype R1R1 R2R2) 
M = cell(1,9);
for i = 1:9
    % Type regarding gene 1:
    j = mod(i-1, 3) + 1;
    % Type regarding gene 2:
    k = ceil(i/3);
    
    % Inheritance matrix for type i seeds:
    M{1, i} = repmat(M1{1, j}, 3, 3) .* repelem(M1{1, k}, 3, 3);
end


% 9 x (n_years+1) array of genotype frequencies in the seed bank. Each  
% column corresponds to one season. Row 1 contains the numbers of 
% S1S1 S2S2 seeds at season start. Row 2 contains the numbers of
% R1S1 S2S2 seeds. Row 3 contains the numbers of R1R1 S2S2 seeds. 
% Row 4 contains the numbers of S1S1 R2S2 seeds. Row 5 contains the 
% numbers of R1S1 R2S2 seeds. Row 6 contains the numbers of R1R1 R2S2 
% seeds. Row 7 contains the numbers of S1S1 R2R2 seeds. Row 8 contains 
% the numbers of R1S1 R2R2 seeds. Row 9 contains the numbers of
% R1R1 R2R2 seeds. 
SB = zeros(9, n_years+1);

% 9 x n_years array of genotype frequencies in the plants. Each column 
% corresponds to one season. Row 1 contains the numbers of 
% S1S1 S2S2 plants at season start. Row 2 contains the numbers of
% R1S1 S2S2 plants. Row 3 contains the numbers of R1R1 S2S2 plants. 
% Row 4 contains the numbers of S1S1 R2S2 plants. Row 5 contains the 
% numbers of R1S1 R2S2 plants. Row 6 contains the numbers of R1R1 R2S2 
% plants. Row 7 contains the numbers of S1S1 R2R2 plants. Row 8 contains 
% the numbers of R1S1 R2R2 plants. Row 9 contains the numbers of
% R1R1 R2R2 plants. 
P = zeros(9, n_years);

% 9 x n_years array of genotype frequencies in the rhizomes. Each column 
% corresponds to one season. Row 1 contains the numbers of 
% S1S1 S2S2 rhizomes at season start. Row 2 contains the numbers of
% R1S1 S2S2 rhizomes. Row 3 contains the numbers of R1R1 S2S2 seeds. 
% Row 4 contains the numbers of S1S1 R2S2 rhizomes. Row 5 contains the 
% numbers of R1S1 R2S2 rhizomes. Row 6 contains the numbers of R1R1 R2S2 
% rhizomes. Row 7 contains the numbers of S1S1 R2R2 rhizomes. Row 8 contains 
% the numbers of R1S1 R2R2 rhizomes. Row 9 contains the numbers of
% R1R1 R2R2 rhizomes.
R = zeros(9, n_years+1);

if seeds
    % Initial seedbank:
    % Absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
    % S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
    % in the initial seed bank:
    SB(:, 1) = dens_seeds * A * repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
        repelem([1-RR2-RW2; RW2; RR2], 3, 1);
end

% Initial rhizomes:
% Absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
% S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
% in the initial rhizomes:
R(:, 1) = dens_rhizomes * A * repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
    repelem([1-RR2-RW2; RW2; RR2], 3, 1);

% Vector of absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
% S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) in the 
% produced seeds: 
S = zeros(9, 1);

% Loop over seasons:
for t = 1:n_years
    % Overall herbicide efficiencies:
    if herb1(t)
        h1 = 1 - (1 - E_early1) * (1 - E_late1);
        h_T1 = 1 - (1 - ((lambda_nodes-1)/lambda_nodes + 1/lambda_nodes * ...
            E_rhizomes1(1 + ~tillage(t))) * E_early1) * (1 - E_late1);
    else 
        h1 = 0;
        h_T1 = 0;
    end
    if herb2(t)
        h2 = 1 - (1 - E_early2) * (1 - E_late2);
        h_T2 = 1 - (1 - ((lambda_nodes-1)/lambda_nodes + 1/lambda_nodes * ...
            E_rhizomes2(1 + ~tillage(t))) * E_early2) * (1 - E_late2);
    else 
        h2 = 0;
        h_T2 = 0;
    end
    
    % Plants emerging from rhizomes and seeds during the season t:
    P(:, t) = diag(repmat([1-h1 (1-k_herb1*h1) 1], 1, 3) .* ...
        repelem([1-h2 (1-k_herb2*h2) 1], 1, 3)) * ...
        (1 - d_seedlings) * p_germ * SB(:, t) + ...
        (1 - d_tillers) * lambda_sec * lambda_nodes * ...
        diag(repmat([1-h_T1 (1-k_herb1*h_T1) 1], 1, 3) .* ...
        repelem([1-h_T2 (1-k_herb2*h_T2) 1], 1, 3)) * R(:, t);
    
    % Self-thinning:
    P(:, t) = P(:, t) * (1 + 1/dens_max * sum(P(:, t)) / A) ^(-1);
    
    if seeds
         % Seeds produced during the season t:
        for i = 1:9
          S(i) = f * P(:, t)' * ...
             diag(repmat([1 1-k_cost1*cost_seeds1 1-cost_seeds1], 1, 3) .* ...
        repelem([1 1-k_cost2*cost_seeds2 1-cost_seeds2], 1, 3)) * ...
             ((1 - p_self) / sum(P(:, t)) * M{1, i} * P(:, t) + ...
             p_self * diag(M{1, i}));
        end
    end
    
    % Rhizomes present at the beginning of next season t+1:
    R(:, t+1) = (1 - p_mort(1 + ~tillage(t+1))) * P(:, t);
    
    % Seed bank present at the beginning of next season t+1:
    if seeds
        SB(:, t+1) = (1 - d_seeds) * (1 - p_germ) * bank * SB(:, t) + S;
    end
end

% 1 x n_years vektor of plant densities. Every entry corresponds to one
% season.
P_dens = sum(P, 1) / A;

end