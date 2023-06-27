%% Generation of data for impact of cycling length on resistance
% Generates the time till the R allele reaches fixation (99.5 %) in 
% a Johnsongrass population treated with ACCase-inhibitor
% and ALS-inhibitor with same efficiency cycled depending on the
% cycle length.

%% Simulation:
% Setting parameters:

% Simulation:
% Field size
A = 10^4;
% Number of years
n_years = 1000;

% Initial seedbank density: 
dens_seeds = 10;
% Initial rhizome density: 
dens_rhizomes = 1;

% Logical value stating whether sexual reproduction is considered
seeds = true; 
% Logical value stating whether a seed bank is considered
bank = true;

% Ecological:
% Proportion of selfpollination 
p_self = 0.95;
% Proportion of seedgermination: 
p_germ = 0.3;
% Proportion of bud sprouting (no tillage, tillage):
p_sprout = [0.2, 0.3]; 

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

%   Fitness cost on seed production associated with resistance
%   against herbicide 1
cost_seeds1 = 0.3;
%   Fitness cost on seed production associated with resistance
%   against herbicide 2
cost_seeds2 = 0.3;

%   Factor reducing the fitness cost of RS type relative to RR
%   type regarding gene 1
k_cost1 = 0.5;
%   Factor reducing the fitness cost of RS type relative to RR 
%   type regarding gene 2
k_cost2 = 0.5;

%   Factor reducing the herbicide efficiency of RS type relative  
%   to SS type regarding gene 1
k_herb1 = 0.5;
%   Factor reducing the herbicide efficiency of RS type relative  
%   to SS type regarding gene 2
k_herb2 = 0.5;


% Antropogenic:
% Herbicide efficacy herbicide 1 (ACCase-inhibitor):
% Seedlings: 
E_seedlings1 = 0.998;
% Tillers (no tillage, tillage): 
E_tillers1 = [0.985, 0.992];
% Herbicide efficacy herbicide 2 (ALS-inhibitor): 
% Seedlings: 
E_seedlings2 = 0.998;
% Tillers (no tillage, tillage): 
E_tillers2 = [0.985, 0.992];


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

% Initial population composition:
% Read table with genotype frewuencies at eqilibrium
T = readtable('Table_Equilibrium_tillage0.txt');
%   RR1: initial fraction of the RR type in seeds and rhizomes regarding
%   herbicide 1
RR1 = T.RR(round(T.Cost,4) == cost_seeds1 & round(T.Het,4) == k_cost1);
%   RW1: initial fraction of the RW type in seeds and rhizomes regarding
%   herbicide 1
RW1 = T.RW(round(T.Cost,4) == cost_seeds1 & round(T.Het,4) == k_cost1);
%   RR2: initial fraction of the RR type in seeds and rhizomes regarding
%   herbicide 2
RR2 = T.RR(round(T.Cost,4) == cost_seeds2 & round(T.Het,4) == k_cost2);
%   RW2: initial fraction of the RW type in seeds and rhizomes regarding
%   herbicide 2
RW2 = T.RW(round(T.Cost,4) == cost_seeds2 & round(T.Het,4) == k_cost2);

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
    SB(:, 1) = dens_seeds * A(i1) * repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
        repelem([1-RR2-RW2; RW2; RR2], 3, 1);
end

% Initial rhizomes:
% Absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
% S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
% in the initial rhizomes:
R(:, 1) = dens_rhizomes * A(i1) * repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
    repelem([1-RR2-RW2; RW2; RR2], 3, 1);

% 2 x 10 array with time till R allele Frequncy exceeds 99.5% 
% with (row 2) and without tillage (row 1) for different 
% cycle length (equaling column number) of herbicide cycling,
% with matrix 1 corresponding to ACCase- and matrix 2 to ALS-inhibitor: 
time_R_fixation = zeros(2, 10, 2);

% Two replicates, one with and one without tillage
for j = 1:2
    %   1 x (n_years+1) vektor of tillage strategy. Each entry 
    %   corresponds to one season and is a logical value stating whether
    %   the soil is tilled at season start
    tillage = (j-1) * ones(1, n_years+1);

    % Loop over different cucling lentgh
    for l = 1:10
        temp1 = [ones(1,l), zeros(1,l)];
        temp2 = [zeros(1,l), ones(1,l)];
        %   1 x n_years vektor of herbicide application. Each entry 
        %   corresponds to one season and is a logical value stating whether 
        %   ACCase-inhibitor is applied
        herb1 = [repmat(temp1, 1, floor((n_years+1) / length(temp1))), ...
                temp1(1, 1:mod(n_years+1, length(temp1)))];
        %   1 x n_years vektor of herbicide application. Each entry 
        %   corresponds to one season and is a logical value stating whether 
        %   ALS-inhibitor is applied
        herb2 = [repmat(temp2, 1, floor((n_years+1) / length(temp2))) ...
                temp2(1, 1:mod(n_years+1, length(temp2)))];

        % Number of buds on initial rhizomes
        b_t = b * (1 + a * dens_rhizomes / d_rhizomes(1 + tillage(1))) ^(-1);

        % Vector of absolute genotype frequencies (SS, RS, RR) in the
        % produced seeds: 
        S = zeros(3, 1);
        
        % Loop over seasons:
        for t = 1:n_years
            % Overall herbicide efficiencies:
            if herb1(t)
                h1 = E_seedlings1;
                h_T1 = E_tillers1(1 + tillage(t));
            else 
                h1 = 0;
                h_T1 = 0;
            end
            if herb2(t)
                h2 = E_seedlings2;
                h_T2 = E_tillers2(1 + tillage(t));
            else 
                h2 = 0;
                h_T2 = 0;
            end
            
            % Plants emerging from rhizomes and seeds during the season t:
            P(:, t) = diag(repmat([1-h1 1-(1-k_herb1)*h1 1], 1, 3) .* ...
                repelem([1-h2 1-(1-k_herb2)*h2 1], 1, 3)) * p_germ * SB(:, t) + ...
                diag(repmat([1-h_T1 1-(1-k_herb1)*h_T1 1], 1, 3) .* ...
                repelem([1-h_T2 1-(1-k_herb2)*h_T2 1], 1, 3)) * ...
                p_sprout(1 + tillage(t)) * b_t * R(:, t);
            
            % Self-thinning:
            P(:, t) = P(:, t) * (1 + 1/dens_max * sum(P(:, t)) / A) ^(-1);
        
            % Density dependant reproduction:
            f_t = f * (1 + a * sum(P(:, t)) / A) ^(-1);
            b_t = b * (1 + a * sum(P(:, t)) / A) ^(-1);

            if (((P(2, t)+P(5, t)+P(8, t)+2*(P(3, t)+P(6, t)+P(9, t)))/ ...
                    (2*sum(P(:, t)))) > 0.995) && ~time_R_fixation(j, l, 1)
                % Time till R allele Frequncy exceeds 99%:
                time_R_fixation(j, l, 1) = t;
            end
            if (((P(4, t)+P(5, t)+P(6, t)+2*(P(7, t)+P(8, t)+P(9, t)))/ ...
                    (2*sum(P(:, t)))) > 0.995) && ~time_R_fixation(j, l, 2)
                % Time till R allele Frequncy exceeds 99%:
                time_R_fixation(j, l, 2) = t;
            end
            if logical(time_R_fixation(j, l, 1)) && logical(time_R_fixation(j, l, 2))
                break
            end
                
            if seeds
                 % Seeds produced during the season t:
                for i = 1:9
                  S(i) = f_t * P(:, t)' * ...
                     diag(repmat([1 1-k_cost1*cost_seeds1 1-cost_seeds1], 1, 3) .* ...
                repelem([1 1-k_cost2*cost_seeds2 1-cost_seeds2], 1, 3)) * ...
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
    end
end



% create a table
T = table;
% assign columns to table
T.Tillage = reshape(repmat([0;1], 1, 10), [], 1);
T.CycleLength = reshape(repmat((1:10), 2, 1), [], 1);
T.FixationTimeACCase = reshape(time_R_fixation(:, :, 1), [], 1);
T.FixationTimeALS = reshape(time_R_fixation(:, :, 2), [], 1);
T.FixationTimeAverage = (reshape(time_R_fixation(:, :, 1), [], 1) + ...
    reshape(time_R_fixation(:, :, 2), [], 1)) / 2;
% write table to text file 
writetable(T, 'Table_Deterministic_R_allele_fixation_time_cycle_length');