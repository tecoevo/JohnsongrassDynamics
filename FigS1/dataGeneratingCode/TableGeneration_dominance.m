%% Generation of data for dominance figure
% Saving deterministic dynamics of a Johnsongrass population over 30 years 
% depending on herbicde application, tillage and seed production, seed bank 
% formation and dominance of resistance allele and associated fitness cost.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters: 
%   A: field size
A = 10^4;
%   n_years: number of years
n_years = 30;

% Initial seedbank density: 
dens_seeds = 10;
% Initial rhizome density: 
dens_rhizomes = 1;

%   control: 2 x (n_years+1) x array of strategies:
%   The first row gives the tillage strategy. The second row indicates the 
%   herbicide apllication. Each column corresponds to one season and the 
%   entries are logical values stating whether the soil is tilled at season 
%   start or whether ACCase-inhibitor is applied, respectively. 
control = zeros(2, n_years+1);
control(:, :, 2) = [ones(1, n_years+1); zeros(1, n_years+1)];
control(:, :, 3) = [zeros(1, n_years+1); ones(1, n_years+1)];
control(:, :, 4) = [ones(1, n_years+1); ones(1, n_years+1)];

%   p_self: proportion of selfpollination 
p_self = [0, 0.95];
%   seeds: logical value stating whether sexual reproduction is considered
seeds = [false, true]; 
%   bank: logical value stating whether a seed bank is considered
bank = [false, true];
%   cost_seeds: fitness cost on seed production associated with resiance
cost_seeds = [0.001, 0.3];
%   k_cost: factor reducing the fitness cost of RS type relative to RR type
k_cost = [0, 0.5, 1];
%   k_herb: factor reducing the herbicide efficiency of RS type relative to 
%   SS type
k_herb = [0, 0.5, 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over all parameter sets
for i4 = 2 %1:length(cost_seeds) 
for i5 = 1:length(k_cost)     

% Initial population composition:
% Read table with genotype frewuencies at eqilibrium
T = readtable('Table_Equilibrium_tillage0.txt');
%   RR: initial fraction of the RR type in seeds and rhizomes
RR = T.RR(round(T.Cost,4) == cost_seeds(i4) & round(T.Het,4) == k_cost (i5));
%   RW: initial fraction of the RW type in seeds and rhizomes
RW = T.RW(round(T.Cost,4) == cost_seeds(i4) & round(T.Het,4) == k_cost (i5));

% Initial seedbank:
% Absolute genotype frequencies (SS, RS, RR) in the initial seed bank:
S0= dens_seeds * A * [1-RR-RW; RW; RR];
% Initial rhizomes:
% Absolute genotype frequencies (SS, RS, RR) in the initial rhizomes:
R0 = dens_rhizomes * A * [1-RR-RW; RW; RR];
% Plant density in presecing season
dens0 = dens_rhizomes / 0.65 ;

for i1 = 2 %1:length(p_self)
for i2 = 2 %1:length(seeds) 
for i3 = 2 %1:length(bank) 
for i6 = 1:length(k_herb) 
for i7 = 3 %1:size(control, 3) 

% gives the dynamics:
%   P: matrix of absolute genotype frequencies in plants
%   R: matrix of absolute genotype frequencies in rhizomes
%   SB: matrix of absolute genotype frequencies in seed bank
%   P_dens: vector of plant densities
[P, R, SB, P_dens] = Dynamics_densityDependance(A, p_self(i1),...
    S0, R0, dens0, control(2, :, i7), control(1, :, i7), seeds(i2), ...
    bank(i3), n_years, cost_seeds(i4), k_cost (i5), k_herb(i6));


% create a table
T = table;
% assign columns to table
T.Season = (1:n_years)';
T.Tillage = control(1, 1:n_years, i7)';
T.ACCase = control(2, 1:n_years, i7)';
T.PlantDensity = P_dens';
T.SSplants = P(1, :)';
T.RSplants = P(2, :)';
T.RRplants = P(3, :)';
T.SSseeds = SB(1, 1:n_years)';
T.RSseeds = SB(2, 1:n_years)';
T.RRseeds = SB(3, 1:n_years)';
T.SSrhizomes = R(1, 1:n_years)';
T.RSrhizomes = R(2, 1:n_years)';
T.RRrhizomes = R(3, 1:n_years)';
% write table to text file 
writetable(T, strcat('Table_Deterministic_control', num2str(i7), ...
    '_selfing', num2str(logical(p_self(i1))), ...
    '_sexualReproduction', num2str(seeds(i2)), '_seedbank', ...
    num2str(bank(i3)), '_cost', num2str(1000*cost_seeds(i4)), ...
    '_kCost', num2str(10*k_cost(i5)), '_kHerb', num2str(10*k_herb(i6))));

end
end
end
end
end
end
end