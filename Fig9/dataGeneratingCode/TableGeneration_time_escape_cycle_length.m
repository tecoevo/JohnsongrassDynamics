%% Generation of data for impact of cycling length on resistance
% Generates the time till the R allele reaches fixation in 
% a Johnsongrass population treated with ACCase-inhibitor
% and ALS-inhibitor with same efficiency cycled depending on the
% cycle length.

%% Simulation:
% Setting parameters:

% Simulation:
% Number of replicates 
n_rep = 10^5;
% Field size
A = 10^4;
% Number of years
n_years = 30;

% Initial seedbank density: 
dens_seeds = 10;
% Initial rhizome density: 
dens_rhizomes = 1;

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

% Array of inheritance matrices for one gene. With the row j column k 
% entry of the matrix i giving the fraction of type i seeds produced 
% by a type j plant pollinated by type k pollen. (1 corresponding 
% to genotype SS, 2  corresponding to genotype RS, 3 corresponding 
% to genotype RR) 
M1 = [1 0.5 0; 0.5 0.25 0; 0 0 0];
M1(:,:,2) = [0 0.5 1; 0.5 0.5 0.5; 1 0.5 0];
M1(:,:,3) = [0 0 0; 0 0.25 0.5; 0 0.5 1];

% Array of inheritance matrices. With the row j column k entry of the
% matrix i giving the fraction of type i seeds produced by a type j
% plant pollinated by type k pollen. (1 corresponding to genotype S1S1 S2S2, 
% 2 corresponding to genotype R1S1 S2S2, 3 corresponding to genotype R1R1 S2S2,
% 4 corresponding to genotype S1S1 R2S2, 5 corresponding to genotype R1S1 R2S2,
% 6 corresponding to genotype R1R1 R2S2, 7 corresponding to genotype S1S1 R2R2,  
% 8 corresponding to genotype R1S1 R2R2, 9 corresponding to genotype R1R1 R2R2) 
M = zeros(9, 9, 9);
for i = 1:9
    % Type regarding gene 1:
    j = mod(i-1, 3) + 1;
    % Type regarding gene 2:
    k = ceil(i/3);
    
    % Inheritance matrix for type i seeds:
    M(:, :, i) = repmat(M1(:, :, j), 3, 3) .* repelem(M1(:, :, k), 3, 3);
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

% 9 x n_years array of genotype frequencies in the tillers. Each column 
% corresponds to one season. Row 1 contains the numbers of 
% S1S1 S2S2 tillers at season start. Row 2 contains the numbers of
% R1S1 S2S2 tillers. Row 3 contains the numbers of R1R1 S2S2 tillers. 
% Row 4 contains the numbers of S1S1 R2S2 tillers. Row 5 contains the 
% numbers of R1S1 R2S2 tillers. Row 6 contains the numbers of R1R1 R2S2 
%tillers. Row 7 contains the numbers of S1S1 R2R2 tillers. Row 8 contains 
% the numbers of R1S1 R2R2 tillers. Row 9 contains the numbers of
% R1R1 R2R2 tillers. 
T = zeros(9, n_years);

% 9 x n_years array of genotype frequencies in the seedlings. Each column 
% corresponds to one season. Row 1 contains the numbers of 
% S1S1 S2S2 seedlings at season start. Row 2 contains the numbers of
% R1S1 S2S2 seedlings. Row 3 contains the numbers of R1R1 S2S2 seedlings. 
% Row 4 contains the numbers of S1S1 R2S2 seedlings. Row 5 contains the 
% numbers of R1S1 R2S2 seedlings. Row 6 contains the numbers of R1R1 R2S2 
% seedlings. Row 7 contains the numbers of S1S1 R2R2 seedlings. Row 8 contains 
% the numbers of R1S1 R2R2 seedlings. Row 9 contains the numbers of
% R1R1 R2R2 seedlings. 
L = zeros(9, n_years);

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

% 2 x 10 array with average time till escape from control, i.e. survival of
% first plant fulli resistant against one of the herbicides, with (row 2)
% and without tillage (row 1) for different cycle length (equaling column 
% number) of herbicide cycling: 
time_escape = zeros(2, 10);
% 2 x 10 array with number of escapes from control, i.e. survival of
% first plant fulli resistant against one of the herbicides, with (row 2)
% and without tillage (row 1) for different cycle length (equaling column 
% number) of herbicide cycling: 
n_escape = zeros(2, 10);

% Two replicates, one with and one without tillage
for til = 1:2
    %   1 x (n_years+1) vektor of tillage strategy. Each entry 
    %   corresponds to one season and is a logical value stating whether
    %   the soil is tilled at season start
    tillage = (til-1) * ones(1, n_years+1);

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

        % Replicates:
        for i = 1:n_rep
            % Initial seedbank:
            % Absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
            % S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
            % in the initial seed bank:
            SB(:, 1) = multinomrand(dens_seeds * A(i1), ...
                       repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
                       repelem([1-RR2-RW2; RW2; RR2], 3, 1));
        
            % Initial rhizomes:
            % Absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
            % S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
            % in the initial rhizomes:
            R(:, 1) = multinomrand(dens_rhizomes * A(i1), ...
                      repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
                      repelem([1-RR2-RW2; RW2; RR2], 3, 1));
    
            % Number of buds on initial rhizomes
            b_t = b * (1 + a * dens_rhizomes / ...
                d_rhizomes(1 + tillage(1))) ^(-1);
    
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
                
                % Number of nodes on rhizomes:
                T(:, t) = [poissrnd(b_t * R(1, t), 1);...
                    poissrnd(b_t * R(2, t), 1);...
                    poissrnd(b_t * R(3, t), 1);...
                    poissrnd(b_t * R(4, t), 1);...
                    poissrnd(b_t * R(5, t), 1);...
                    poissrnd(b_t * R(6, t), 1);...
                    poissrnd(b_t * R(7, t), 1);...
                    poissrnd(b_t * R(8, t), 1);...
                    poissrnd(b_t * R(9, t), 1)];
                % Tillers emerging from rhizome nodes:
                T(:, t) = binomrand(T(:, t), p_sprout(1 + tillage(t)));
                % Mortaility due to herbicide application:
                T(:, t) = [binomrand(T(1, t), (1-h_T1)*(1-h_T2));...
                    binomrand(T(2, t), (1-(1-k_herb1)*h_T1)*(1-h_T2));...
                    binomrand(T(3, t), (1-h_T2));...
                    binomrand(T(4, t), (1-h_T1)*(1-(1-k_herb2)*h_T2));...
                    binomrand(T(5, t), (1-(1-k_herb1)*h_T1)*(1-(1-k_herb2)*h_T2));...
                    binomrand(T(6, t), (1-(1-k_herb2)*h_T2));...
                    binomrand(T(7, t), (1-h_T1));...
                    binomrand(T(8, t), (1-(1-k_herb1)*h_T1));...
                    T(9, t)];
                
                % Seedlings emerging from seeds:
                L(:, t) = binomrand(SB(:, t), p_germ);
                % Update of seedbank:
                SB(:, t+1) = SB(:, t) - L(:, t);
                % Mortaility due to herbicide application:
                L(:, t) = [binomrand(L(1, t), (1-h1)*(1-h2));...
                    binomrand(L(2, t), (1-(1-k_herb1)*h1)*(1-h2));...
                    binomrand(L(3, t), (1-h2));...
                    binomrand(L(4, t), (1-h1)*(1-(1-k_herb2)*h2));...
                    binomrand(L(5, t), (1-(1-k_herb1)*h1)*(1-(1-k_herb2)*h2));...
                    binomrand(L(6, t), (1-(1-k_herb2)*h2));...
                    binomrand(L(7, t), (1-h1));...
                    binomrand(L(8, t), (1-(1-k_herb1)*h1));...
                    T(9, t)];
            
                % Plants emerging from rhizomes and seeds during the season t:
                P(:, t) = T(:, t) + L(:, t);
        
                 % Self-thinning:
                 P(:, t) = binomrand(P(:, t), ...
                     (1 + 1/dens_max * sum(P(:, t)) / A) ^(-1));
        
                % Density dependant reproduction:
                f_t = f * (1 + a * sum(P(:, t)) / A) ^(-1);
                b_t = b * (1 + a * sum(P(:, t)) / A) ^(-1);
    
                if (P(3, t)+P(6, t)+P(9, t))>0 || (P(7, t)+P(8, t)+P(9, t))>0
                    % Time till first plant fully resistant against one of the
                    % herbicides survives
                    time_escape(til, l) = time_escape(til, l) + t;
                    % Number of replicates with a plant fully resistant a
                    % gainst one of the herbicides survives
                    n_escape(til, l) = n_escape(til, l) + 1;
                    break
                end
                    
                % Number of self-pollinated plants of each type:
                P_self = binomrand(P(:, t), p_self);
                
                % Number of seeds produced by the self-pollinated plants of
                % each tupe: 
                n_self = [poissrnd(f_t * P_self(1), 1); ...
                    poissrnd(f_t * (1-k_cost1*cost_seeds1) * P_self(2), 1); ...
                    poissrnd(f_t * (1-cost_seeds1) * P_self(3), 1); ...
                    poissrnd(f_t * (1-k_cost2*cost_seeds2) * P_self(4), 1); ...
                    poissrnd(f_t * (1-k_cost1*cost_seeds1)*(1-k_cost2*cost_seeds2) * P_self(5), 1); ...
                    poissrnd(f_t * (1-cost_seeds1)*(1-k_cost2*cost_seeds2) * P_self(6), 1); ...
                    poissrnd(f_t * (1-cost_seeds2) * P_self(7), 1); ...
                    poissrnd(f_t * (1-k_cost1*cost_seeds1)*(1-cost_seeds2) * P_self(8), 1); ...
                    poissrnd(f_t * (1-cost_seeds1)*(1-cost_seeds2) * P_self(9), 1)]; 
        
                % Number of seeds of each type produced by self-pollinated
                % plants:
                S_self = multinomrand(n_self(1), reshape(M(1, 1, :), 9, 1))' + ...
                    multinomrand(n_self(2), reshape(M(2, 2, :), 9, 1))' + ...
                    multinomrand(n_self(3), reshape(M(3, 3, :), 9, 1))' + ...
                    multinomrand(n_self(4), reshape(M(4, 4, :), 9, 1))' + ...
                    multinomrand(n_self(5), reshape(M(5, 5, :), 9, 1))' + ...
                    multinomrand(n_self(6), reshape(M(6, 6, :), 9, 1))' + ...
                    multinomrand(n_self(7), reshape(M(7, 7, :), 9, 1))' + ...
                    multinomrand(n_self(8), reshape(M(8, 8, :), 9, 1))' + ...
                    multinomrand(n_self(9), reshape(M(9, 9, :), 9, 1))';
              
                % Number of seeds produced by the cross-pollinated plants of
                % each tupe: 
                n_cross = [poissrnd(f * (P(1, t) - P_self(1)), 1); ...
                    poissrnd(f * (1-k_cost1*cost_seeds1) * (P(2, t) - P_self(2)), 1); ...
                    poissrnd(f * (1-cost_seeds1) * (P(3, t) - P_self(3)), 1); ...
                    poissrnd(f * (1-k_cost2*cost_seeds2) * (P(4, t) - P_self(4)), 1); ...
                    poissrnd(f * (1-k_cost1*cost_seeds1)*(1-k_cost2*cost_seeds2) * ...
                    (P(5, t) - P_self(5)), 1); ...
                    poissrnd(f * (1-cost_seeds1)*(1-k_cost2*cost_seeds2) * ...
                    (P(6, t) - P_self(6)), 1); ...
                    poissrnd(f * (1-cost_seeds2) * (P(7, t) - P_self(7)), 1); ...
                    poissrnd(f * (1-k_cost1*cost_seeds1)*(1-cost_seeds2) * ...
                    (P(8, t) - P_self(8)), 1); ...
                    poissrnd(f * (1-cost_seeds1)*(1-cost_seeds2) * (P(9, t) - P_self(9)), 1)]; 
        
                % 9 x 9 array with numbers of seeds procuced by the different
                % mating. With the row j and column k entry giving the number 
                % of seeds produced by plants of type j cross-pollinated with 
                % type k pollen:
                C = multinomrand(n_cross, P(:, t)'/sum(P(:, t)));
                % Seeds produced by the cross-pollinated plants:
                S_cross = zeros(9, 1);
                for j = 1:9
                    for k = 1:9
                        S_cross = S_cross + multinomrand(C(j, k), ...
                            reshape(M(j, k, :), 9, 1))';
                    end
                end
                
                % Seeds produced during the season t:
                S = S_self + S_cross;
        
                % Matrix with mutation probabilities for one gene. The row i column
                % j entry gives the probability by which a mutation from type i to
                % type j occures. 
                mut = [(1-mu)^2, 2*mu*(1-mu), mu^2; 0, (1-mu), mu; 0, 0, 1];
                % Matrix with mutation probabilities for two genes. The row i column
                % j entry gives the probability by which a mutation from type i to
                % type j occures. 
                mut = repmat(mut, 3, 3) .* repelem(mut, 3, 3);
               
                % Mutations in seeds:
                S = multinomrand(S(1), mut(1, :))' + ...
                multinomrand(S(2), mut(2, :))' + ... 
                multinomrand(S(3), mut(3, :))' + ...
                multinomrand(S(4), mut(4, :))' + ...
                multinomrand(S(5), mut(5, :))' + ...
                multinomrand(S(6), mut(6, :))' + ...
                multinomrand(S(7), mut(7, :))' + ...
                multinomrand(S(8), mut(8, :))' + ...
                [zeros(8,1); S(9)];
        
                % Seed bank present at the beginning of next season t+1:
                SB(:, t+1) = bank * binomrand(SB(:, t+1), 1 - d_bank) + ...
                    binomrand(S, 1 - d_seeds);
                  
                % Rhizomes present at the beginning of next season t+1:
                R(:, t+1) = (1 - d_rhizomes(1 + tillage(t+1))) * P(:, t);
            end
        end
    end
end



% create a table
T = table;
% assign columns to table
T.Tillage = reshape(repmat([0;1], 1, 10), [], 1);
T.CycleLength = reshape(repmat((1:10), 2, 1), [], 1);
T.AverageEscapeYear = reshape(time_escape(:, :), [], 1) ./ reshape(n_escape(:, :), [], 1);
T.FractionEscapes = reshape(n_escape(:, :), [], 1) / n_rep;
% write table to text file 
writetable(T, 'Table_Stochastic_EscapeTime_EscapeNumber_cycle_length');