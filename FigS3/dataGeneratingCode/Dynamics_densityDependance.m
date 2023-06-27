function [P, R, SB, P_dens] = ...
    Dynamics_densityDependance(A, p_self, S0, R0, dens0, ...
    herbs, tillage, seeds, bank, n_years, cost_seeds, k_cost, k_herb)
%Dynamics_selfThinning Gives the deterministic dynamics of a Johnsongrass
% population over 30 years depending on herbicde application, tillage and
% seed production.
%
%   Input: 
%   A: field size
%   p_self: proportion of selfpollination
%   S0: 3 x 1 vector of absolute genotype frequencies in initial seeds
%   R0: 3 x 1 vector of absolute genotype frequencies in initial rhizomes
%   dens0: plant density in preceding season
%   herbs: 1 x n_years vektor of herbicide application. Each entry 
%   corresponds to one season and is a logical value stating whether 
%   ACCase-inhibitor is applied
%   tillage: 1 x (n_years+1) vektor of tillage strategy. Each entry 
%   corresponds to one season and is a logical value stating whether
%   the soil is tilled at season start
%   seeds: logical value stating whether seed production is considered
%   bank: logical value stating whether a seed bank is considered
%   n_years: number of years
%   cost_seeds: fitness cost on seed production associated with resiance 
%   k_cost: factor reducing the fitness cost of RS type relative to RR type
%   k_herb: factor reducing the herbicide efficiency of RS type relative to 
%           SS type
%
%   Output:
%   P: matrix of absolute genotype frequencies in plants
%   R: matrix of absolute genotype frequencies in rhizomes
%   SB: matrix of absolute genotype frequencies in seed bank
%   P_dens: vector of plant densities

%% Simulation:
% Setting parameters:

% Ecological:
% Proportion of seedgermination: 
p_germ = 0.3;
% Proportion of bud sprouting (no tillage, tillage):
p_sprout = [0.2, 0.3];

% Fecundity, i.e. number of seeds produced per plant:
f = 13000; 

% number of rhizome buds produced per plant:
b = 140;

% Rhizome winter mortality (no tillage, tillage): 
d_rhizomes = [0.35, 0.6];

% Loss and natural mortality of fresh seeds over winter:
d_seeds = 0.94;
% Natural yearly seed mortality in the seedbank:
d_bank = 0.48;

% Highest possible plant density:
dens_max = 220;
% Area needed for a plant to produce f seeds and b buds:
a = 0.1;

% Evolutionary:
% Mutation rate:
mu = 10^(-8);


% Antropogenic:
% Herbicide efficacy: 
% Seedlings: 
E_seedlings = 0.998;
% Tillers (no tillage, tillage): 
E_tillers = [0.985, 0.992];


% Cell array of inheritance matrices. With the row j column k entry of the
% cell i matrix giving the fraction of type i seeds produced by a type j
% plant pollinated by type k pollen. (1 corresponding to genotype SS, 2 
% corresponding to genotype RS, 3 corresponding to genotype RR) 
MI = {[1 0.5 0; 0.5 0.25 0; 0 0 0], [0 0.5 1; 0.5 0.5 0.5; 1 0.5 0],...
    [0 0 0; 0 0.25 0.5; 0 0.5 1]};
% Add Mutation:
M{1} = (1 - mu)^2 * MI{1};
M{2} = 2 * mu * (1 - mu) * MI{1} + (1 - mu) * MI{2};
M{3} = mu^2 * MI{1} + mu * MI{2} + MI{3};

% 3 x (n_years+1) array of genotype frequencies in the seed bank. Each  
% column corresponds to one season. Row 1 contains the numbers of SS seeds  
% at season start. Row 1 contains the numbers of RS seeds. Row 1 contains 
% the numbers of RR seeds. 
SB = zeros(3, n_years+1);

% 3 x n_years array of genotype frequencies in the plants. Each column 
% corresponds to one season. Row 1 contains the numbers of SS plants at the 
% end of the season. Row 1 contains the numbers of RS plants. Row 1 
% contains the numbers of RR plants. 
P = zeros(3, n_years);

% 3 x n_years array of genotype frequencies in the rhizomes. Each column 
% corresponds to one season. Row 1 contains the numbers of SS rhizomess at
% the end of the season. Row 1 contains the numbers of RS rhizomes. Row 1 
% contains the numbers of RR rhizomes. 
R = zeros(3, n_years+1);

if seeds
    % Initial seedbank:
    % Absolute genotype frequencies (SS, RS, RR) in the initial seed bank:
    SB(:, 1) = S0;
end

% Initial rhizomes:
% Absolute genotype frequencies (SS, RS, RR) in the initial rhizomes:
R(:, 1) = R0;

% Number of buds on initial rhizomes
b_t = b * (1 + a * dens0) ^(-1);

% Vector of absolute genotype frequencies (SS, RS, RR) in the produced
% seeds: 
S = zeros(3, 1);

% Loop over seasons:
for t = 1:n_years
    % Overall herbicide efficiency:
    if herbs(t)
        h = E_seedlings;
        h_T = E_tillers(1 + tillage(t));
    else 
        h = 0;
        h_T = 0;
    end
    
    % Plants emerging from rhizomes and seeds during the season t:
    P(:, t) = diag([1-h 1-(1-k_herb)*h 1]) * p_germ * SB(:, t) + ...
        diag([1-h_T 1-(1-k_herb)*h_T 1]) * ...
        p_sprout(1 + tillage(t)) * b_t * R(:, t);
    
    % Self-thinning:
    P(:, t) = P(:, t) * (1 + 1/dens_max * sum(P(:, t)) / A) ^(-1);
    
    % Density dependant reproduction:
    f_t = f * (1 + a * sum(P(:, t)) / A) ^(-1);
    b_t = b * (1 + a * sum(P(:, t)) / A) ^(-1);

    if seeds
         % Seeds produced during the season t:
        for i = 1:3
          S(i) = f_t * P(:, t)' * ...
             diag([1 1-k_cost*cost_seeds 1-cost_seeds]) * ...
             ((1 - p_self) / sum(P(:, t)) * M{1, i} * P(:, t) + ...
             p_self * diag(M{1, i}));
        end

        % Seed bank present at the beginning of next season t+1:
        SB(:, t+1) = (1 - d_bank) * (1 - p_germ) * ...
            bank * SB(:, t) + ...
            (1 - d_seeds) * S;
    end
    
    % Rhizomes present at the beginning of next season t+1:
    R(:, t+1) = (1 - d_rhizomes(1 + tillage(t+1))) * P(:, t);
end

% 1 x n_years vektor of plant densities. Every entry corresponds to one
% season.
P_dens = sum(P, 1) / A;

end