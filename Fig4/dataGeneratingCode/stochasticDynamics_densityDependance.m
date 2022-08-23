function [P, R, SB, P_dens] = ...
    stochasticDynamics_densityDependance(A, p_self, S0, R0, dens0, herbs, ...
    tillage, seeds, bank, n_years, cost_seeds, k_cost, k_herb)
% Gives the deterministic dynamics of a Johnsongrass
% population over 30 years depending on herbicde application and tillage.
%
%   Input: 
%   p_self: proportion of selfpollination
%   A: field size
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


% Array of inheritance matrices. With the row j column k entry of the
% matrix i giving the fraction of type i seeds produced by a type j
% plant pollinated by type k pollen. (1 corresponding to genotype SS, 2 
% corresponding to genotype RS, 3 corresponding to genotype RR) 
M = [1 0.5 0; 0.5 0.25 0; 0 0 0];
M(:,:,2) = [0 0.5 1; 0.5 0.5 0.5; 1 0.5 0];
M(:,:,3) = [0 0 0; 0 0.25 0.5; 0 0.5 1];

% 3 x (n_years+1) array of genotype frequencies in the seed bank. Each  
% column corresponds to one season. Row 1 contains the numbers of SS seeds  
% at season start. Row 1 contains the numbers of RS seeds. Row 1 contains 
% the numbers of RR seeds. 
SB = zeros(3, n_years+1);

% 3 x n_years array of genotype frequencies in the tillers. Each column 
% corresponds to one season. Row 1 contains the numbers of SS tillers in 
% the season. Row 1 contains the numbers of RS tillers. Row 1 contains the 
% numbers of RR tillers. 
T = zeros(3, n_years);

% 3 x n_years array of genotype frequencies in the seedlings. Each column 
% corresponds to one season. Row 1 contains the numbers of SS tillers in 
% the season. Row 1 contains the numbers of RS seedlings. Row 1 contains 
% the numbers of RR seedlings. 
L = zeros(3, n_years);

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
    
    % Number of nodes on rhizomes:
    T(:, t) = [poissrnd(b_t * R(1, t), 1); poissrnd(b_t * R(2, t), 1);...
        poissrnd(b_t * R(3, t), 1)];
    % Tillers emerging from rhizome nodes:
    T(:, t) = binomrand(T(:, t), p_sprout(1 + tillage(t)));
    % Mortaility due to herbicide application:
    T(:, t) = [binomrand(T(1, t), 1-h_T);...
        binomrand(T(2, t), 1-(1-k_herb)*h_T); T(3, t)];
    
    % Seedlings emerging from seeds:
    L(:, t) = binomrand(SB(:, t), p_germ);
    % Update of seedbank:
    SB(:, t+1) = SB(:, t) - L(:, t);
    % Mortaility due to herbicide application:
    L(:, t) = [binomrand(L(1, t), 1-h);...
        binomrand(L(2, t), 1-(1-k_herb)*h); L(3, t)];

    % Plants emerging from rhizomes and seeds during the season t:
    P(:, t) = T(:, t) + L(:, t);
    
    % Self-thinning:
    P(:, t) = binomrand(P(:, t), (1 + 1/dens_max * sum(P(:, t)) / A) ^(-1));

    % Density dependant reproduction:
    f_t = f * (1 + a * sum(P(:, t)) / A) ^(-1);
    b_t = b * (1 + a * sum(P(:, t)) / A) ^(-1);
    
    if seeds
        % Number of self-pollinated plants of each type:
        P_self = binomrand(P(:, t), p_self);
        
        % Number of seeds produced by the self-pollinated plants of
        % each type: 
        n_self = [poissrnd(f_t * P_self(1), 1); ...
            poissrnd(f_t * (1-k_cost*cost_seeds) * P_self(2), 1); ...
            poissrnd(f_t * (1-cost_seeds) * P_self(3), 1)];

        % Number of seeds of each type produced by self-pollinated
        % plants:
        S_self = multinomrand(n_self(1), reshape(M(1, 1, :), 3, 1))' + ...
            multinomrand(n_self(2), reshape(M(2, 2, :), 3, 1))' + ...
            multinomrand(n_self(3), reshape(M(3, 3, :), 3, 1))';
      
        % Number of seeds produced by the cross-pollinated plants of
        % each tupe: 
        n_cross = [poissrnd(f * (P(1, t) - P_self(1)), 1); ...
            poissrnd(f * (1-k_cost*cost_seeds) * (P(2, t) - P_self(2)), 1); ...
            poissrnd(f * (1-cost_seeds) * (P(3, t) - P_self(3)), 1)];

        % 3 x 3 array with numbers of seeds procuced by the different
        % mating. With the row j and column k entry giving the number 
        % of seeds produced by plants of type j cross-pollinated with 
        % type k pollen:
        C = multinomrand(n_cross, P(:, t)'/sum(P(:, t)));
        % Seeds produced by the cross-pollinated plants:
        S_cross = zeros(3, 1);
        for j = 1:3
            for k = 1:3
                S_cross = S_cross + multinomrand(C(j, k), ...
                    reshape(M(j, k, :), 3, 1))';
            end
        end
        
        % Seeds produced during the season t:
        S = S_self + S_cross;

        % Mutations in seeds:
        S = multinomrand(S(1),[(1-mu)^2, 2*mu*(1-mu), mu^2])' + ...
            multinomrand(S(2), [0, (1-mu), mu])' + [0; 0; S(3)];

        % Seed bank present at the beginning of next season t+1:
        SB(:, t+1) = bank * binomrand(SB(:, t+1), 1 - d_bank) + ...
            binomrand(S, 1 - d_seeds);
    end
    
    % Rhizomes present at the beginning of next season t+1:
    R(:, t+1) = binomrand(P(:, t), 1 - d_rhizomes(1 + tillage(t+1)));
end

% 1 x n_years vektor of plant densities. Every entry corresponds to one
% season.
P_dens = sum(P, 1) / A;

end