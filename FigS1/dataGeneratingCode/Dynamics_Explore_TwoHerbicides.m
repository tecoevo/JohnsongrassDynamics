% Gives deterministic dynamics of a Johnsongrass population over 30 years
% depending on herbicde application and tillage.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters to choose: 
%   control: 3 x (n_years+1) x array of strategies:
%   The first row gives the tillage strategy. The second and third rows
%   indicates the herbicide apllication. Each column corresponds to one 
%   season and the entries are logical values stating whether the soil is 
%   tilled at season start, whether ACCase-inhibitor is applied or whether 
%   ALS-inhibitor is applied, respectively. 
control = [ones(1, n_years+1); ones(1, n_years+1); ...
                    zeros(1, n_years+1)];
control(:, :, 2) = [ones(1, n_years+1); zeros(1, n_years+1); ...
                    zeros(1, n_years+1)];
control(:, :, 3) = [zeros(1, n_years+1); ones(1, n_years+1); ...
                    zeros(1, n_years+1)];
control(:, :, 4) = [zeros(1, n_years+1); ones(1, n_years+1); ...
                    ones(1, n_years+1)];
control(:, :, 5) = [ones(1, n_years+1); ones(1, n_years+1); ...
                    ones(1, n_years+1)];


%   seeds: logical value stating whether sexual reproduction is considered
seeds = [false, true]; 
%   bank: logical value stating whether a seed bank is considered
bank = [false, true];
%   cost_seeds1: fitness cost on seed production associated with resistance
%   against herbicide 1
cost_seeds1 = [0.001, 0.3];
%   cost_seeds2: fitness cost on seed production associated with resistance
%   against herbicide 2
cost_seeds2 = [0.001, 0.3];
%   k_cost1: factor reducing the fitness cost of RS type relative to RR
%   type regarding gene 1
k_cost1 = [0, 0.5, 1];
%   k_cost2: factor reducing the fitness cost of RS type relative to RR 
%   type regarding gene 2
k_cost2 = [0, 0.5, 1];
%   k_herb1: factor reducing the herbicide efficiency of RS type relative  
%   to SS type regarding gene 1
k_herb1 = [0, 0.5, 1];
%   k_herb2: factor reducing the herbicide efficiency of RS type relative  
%   to SS type regarding gene 2
k_herb2 = [0, 0.5, 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over all parameter sets
%for i2 = 1:length(seeds) 
i2 = 2;
%for i3 = 1:length(bank) 
i3 = 2;
%for i4 = 1:length(cost_seeds1) 
i4 = 2;
%for i5 = 1:length(cost_seeds2) 
i5 = 2;
%for i6 = 1:length(k_cost1)   
i6 = 2;
%for i7 = 1:length(k_cost2) 
i7 = 2;
%for i8 = 1:length(k_herb1) 
i8 = 2;
%for i9 = 1:length(k_herb2) 
i9 = 2;
for i10 = 1:size(control, 3) 

% Initial population composition:
% Read table with genotype frewuencies at eqilibrium
T = readtable('Table_Equilibrium_herbs0_tillage0.txt');
%   RR1: initial fraction of the RR type in seeds and rhizomes regarding
%   herbicide 1
RR1 = T.RR(round(T.Cost,4) == cost_seeds1(i4) & round(T.Het,4) == k_cost1(i6));
%   RW1: initial fraction of the RW type in seeds and rhizomes regarding
%   herbicide 1
RW1 = T.RW(round(T.Cost,4) == cost_seeds1(i4) & round(T.Het,4) == k_cost1(i6));
%   RR2: initial fraction of the RR type in seeds and rhizomes regarding
%   herbicide 2
RR2 = T.RR(round(T.Cost,4) == cost_seeds2(i5) & round(T.Het,4) == k_cost2(i7));
%   RW2: initial fraction of the RW type in seeds and rhizomes regarding
%   herbicide 2
RW2 = T.RW(round(T.Cost,4) == cost_seeds2(i5) & round(T.Het,4) == k_cost2(i7));
  

% gives the dynamics:
%   P: matrix of absolute genotype frequencies in plants
%   R: matrix of absolute genotype frequencies in rhizomes
%   SB: matrix of absolute genotype frequencies in seed bank
%   P_dens: vector of plant densities
[P, R, SB, P_dens] = Dynamics_TwoHerbicides(RR1, RW1, RR2, RW2, ...
    control(2, :, i10), control(3, :, i10), control(1, :, i10), seeds(i2), ...
    bank(i3), n_years, ...
    cost_seeds1(i4), cost_seeds2(i5), k_cost1(i6), k_cost2(i7), ...
    k_herb1(i8), k_herb2(i9));


% create a table
T = table;
% assign columns to table
T.Season = (1:length(P_dens))';
T.PlantDensity = P_dens';
T.SSSSplants = P(1, :)';
T.RSSSplants = P(2, :)';
T.RRSSplants = P(3, :)';
T.SSRSplants = P(4, :)';
T.RSRSplants = P(5, :)';
T.RRRSplants = P(6, :)';
T.SSRRplants = P(7, :)';
T.RSRRplants = P(8, :)';
T.RRRRplants = P(9, :)';
T.R1plants = (T.RSSSplants + T.RSRSplants + T.RSRRplants + ...
    2.*(T.RRSSplants + T.RRRSplants + T.RRRRplants)) ./ ...
    (2.*(T.SSSSplants + T.SSRSplants + T.SSRRplants + ...
    T.RSSSplants + T.RSRSplants + T.RSRRplants +...
    T.RRSSplants + T.RRRSplants + T.RRRRplants));
T.R2plants = (T.SSRSplants + T.RSRSplants + T.RRRSplants + ...
    2.*(T.SSRRplants + T.RSRRplants + T.RRRRplants)) ./ ...
    (2.*(T.SSSSplants + T.SSRSplants + T.SSRRplants + ...
    T.RSSSplants + T.RSRSplants + T.RSRRplants +...
    T.RRSSplants + T.RRRSplants + T.RRRRplants));
T.SSSSseeds = SB(1, 1:length(P_dens))';
T.RSSSseeds = SB(2, 1:length(P_dens))';
T.RRSSseeds = SB(3, 1:length(P_dens))';
T.SSRSseeds = SB(4, 1:length(P_dens))';
T.RSRSseeds = SB(5, 1:length(P_dens))';
T.RRRSseeds = SB(6, 1:length(P_dens))';
T.SSRRseeds = SB(7, 1:length(P_dens))';
T.RSRRseeds = SB(8, 1:length(P_dens))';
T.RRRRseeds = SB(9, 1:length(P_dens))';
T.R1seeds = (T.RSSSseeds + T.RSRSseeds + T.RSRRseeds + ...
    2.*(T.RRSSseeds + T.RRRSseeds + T.RRRRseeds)) ./ ...
    (2.*(T.SSSSseeds + T.SSRSseeds + T.SSRRseeds + ...
    T.RSSSseeds + T.RSRSseeds + T.RSRRseeds +...
    T.RRSSseeds + T.RRRSseeds + T.RRRRseeds));
T.R2seeds = (T.SSRSseeds + T.RSRSseeds + T.RRRSseeds + ...
    2.*(T.SSRRseeds + T.RSRRseeds + T.RRRRseeds)) ./ ...
    (2.*(T.SSSSseeds + T.SSRSseeds + T.SSRRseeds + ...
    T.RSSSseeds + T.RSRSseeds + T.RSRRseeds +...
    T.RRSSseeds + T.RRRSseeds + T.RRRRseeds));
T.SSSSrhizomes = R(1, 1:length(P_dens))';
T.RSSSrhizomes = R(2, 1:length(P_dens))';
T.RRSSrhizomes = R(3, 1:length(P_dens))';
T.SSRSrhizomes = R(4, 1:length(P_dens))';
T.RSRSrhizomes = R(5, 1:length(P_dens))';
T.RRRSrhizomes = R(6, 1:length(P_dens))';
T.SSRRrhizomes = R(7, 1:length(P_dens))';
T.RSRRrhizomes = R(8, 1:length(P_dens))';
T.RRRRrhizomes = R(9, 1:length(P_dens))';
T.R1rhizomes = (T.RSSSrhizomes + T.RSRSrhizomes + T.RSRRrhizomes + ...
    2.*(T.RRSSrhizomes + T.RRRSrhizomes + T.RRRRrhizomes)) ./ ...
    (2.*(T.SSSSrhizomes + T.SSRSrhizomes + T.SSRRrhizomes+ ...
    T.RSSSrhizomes + T.RSRSrhizomes + T.RSRRrhizomes +...
    T.RRSSrhizomes + T.RRRSrhizomes + T.RRRRrhizomes));
T.R2rhizomes = (T.SSRSrhizomes + T.RSRSrhizomes + T.RRRSrhizomes + ...
    2.*(T.SSRRrhizomes + T.RSRRrhizomes + T.RRRRrhizomes)) ./ ...
    (2.*(T.SSSSrhizomes + T.SSRSrhizomes + T.SSRRrhizomes+ ...
    T.RSSSrhizomes + T.RSRSrhizomes + T.RSRRrhizomes +...
    T.RRSSrhizomes + T.RRRSrhizomes + T.RRRRrhizomes));
% write table to text file 
writetable(T, strcat('Table_Deterministic_TwoHerbicides_control', ...
    num2str(i10), '_sexualReproduction', num2str(seeds(i2)), '_seedbank', ...
    num2str(bank(i3)), '_cost1_', num2str(1000*cost_seeds1(i4)), ...
    '_cost2_', num2str(1000*cost_seeds2(i5)), '_kCost1_', ...
    num2str(10*k_cost1(i6)), '_kCost2_', num2str(10*k_cost2(i7)), ...
    '_kHerb1_', num2str(10*k_herb1(i8)), '_kHerb2_', ...
    num2str(10*k_herb2(i9))));

end
%end
%end
%end
%end
%end
%end
%end
%end