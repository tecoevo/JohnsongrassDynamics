%% Generation of data for seedbank figure
% Stochastic dynamics of a Johnsongrass population over 30 years 
% depending on herbicde application, tillage and seed production, seed bank 
% formation and dominance of resistance allele and associated fitness cost.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters: 
%   n_rep: number of replicates 
n_rep = 10^3;

%   A: field size
A = [10^4, 10^5];
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
for i5 = 2 %1:length(k_cost)     

% Initial population composition:
% Read table with genotype frewuencies at eqilibrium
T = readtable('Table_Equilibrium_tillage0.txt');
%   RR: initial fraction of the RR type in seeds and rhizomes
RR = T.RR(round(T.Cost,4) == cost_seeds(i4) & round(T.Het,4) == k_cost (i5));
%   RW: initial fraction of the RW type in seeds and rhizomes
RW = T.RW(round(T.Cost,4) == cost_seeds(i4) & round(T.Het,4) == k_cost (i5));


for i1 = 2 %1:length(p_self)
for i2 = 2 %1:length(seeds) 
for i3 = 1:length(bank) 
for i6 = 2 %1:length(k_herb) 
for i7 = 4 %1:size(control, 3) 
for i8 = 1:length(A)

% n_rep x n_years x 10 array with simulation results. 
% Each row corresponds to one replicate. 
% Each column corresponds to one season of the simulation. 
% Pages contain population level variables of plants, rhizomes and 
% seedbank. 
% Page 1 contains the plant densities at the end of the seasons.
% Pages 2-4 contain the relative genotype (SS, RS, RR) frequencies in 
% plants at the end of the seasons.
% Pages 5-7 contain the absolute genotype (SS, RS, RR) frequencies in
% primary rhizomes of the seasons.
% Pages 8-10 contain the absolute genotype (SS, RS, RR) frequencies
% in seedbank of the seasons.
Sim = zeros(n_rep, n_years, 10);

% Replicates:
for i = 1:n_rep

    % Initial seedbank:
    % Absolute genotype frequencies (SS, RS, RR) in the initial seed bank:
    S0= multinomrand(dens_seeds * A(i8), [1-RR-RW; RW; RR])';
    % Initial rhizomes:
    % Absolute genotype frequencies (SS, RS, RR) in the initial rhizomes:
    R0 = multinomrand(dens_rhizomes * A(i8), [1-RR-RW; RW; RR])';
    % Plant density in presecing season
    dens0 = dens_rhizomes / 0.65 ;

    % gives the dynamics:
    %   P: matrix of absolute genotype frequencies in plants
    %   R: matrix of absolute genotype frequencies in rhizomes
    %   SB: matrix of absolute genotype frequencies in seed bank
    %   P_dens: vector of plant densities
    [P, R, SB, P_dens] = stochasticDynamics_densityDependance(A(i8), ...
        p_self(i1), S0, R0, dens0, control(2, :, i7), control(1, :, i7), ...
        seeds(i2), bank(i3), n_years, cost_seeds(i4), k_cost (i5), ...
        k_herb(i6));

    % Save simulation results of the current run:
    Sim(i, :, 1) = P_dens';
    Sim(i, : , 2:4) = P';
    Sim(i, :, 5:7) = R(:, 1:n_years)';
    Sim(i, :, 8:10) = SB(:, 1:n_years)';

end

% n_years x 1 vektor of avarage plant densities. Each entry corresponds to
% one season and gives the average plant density over the replicates. 
avg_dens = sum(Sim(:, 1:n_years, 1), 1)' / n_rep;

% n_years x 3 array of avarage absolute genotype frequencies in plants. 
% Each row corresponds to one season. Column 1 contains the avarage 
% frequency of SS plants. Column 2 contains the avarage frequency of RS 
% plants. Column 3 contains the avarage frequency of RR plants. 
avg_plantsG = reshape(sum(Sim(:, 1:n_years, 2:4), 1) / n_rep, n_years, 3);

% (n_years+1) x 3 array of avarage absolute genotype frequencies in 
% rhizomes. 
% Each row corresponds to one season. Column 1 contains the avarage 
% frequency of SS rhizomes. Column 2 contains the avarage frequency of RS 
% rhizomes. Column 3 contains the avarage frequency of RR rhizomes. 
avg_rhizomesG = reshape(sum(Sim(:, :, 5:7), 1) / n_rep, n_years, 3);

% (n_years+1) x 3 array of avarage absolute genotype frequencies in seeds. 
% Each row corresponds to one season. Column 1 contains the avarage 
% frequency of SS seeds. Column 2 contains the avarage frequency of RS 
% seeds. Column 3 contains the avarage frequency of RR seeds. 
avg_seedsG = reshape(sum(Sim(:, :, 8:10), 1) / n_rep, n_years, 3);

% create a table
T = table;
% assign columns to table
T.Season = [reshape(repmat(1:n_years, n_rep, 1), ...
    n_rep*n_years, 1); (1:n_years)'];
T.Run = [repmat((1:n_rep)', n_years, 1); NaN(n_years, 1)];
T.Tillage = [reshape(repmat(control(1, 1:n_years, i7), n_rep, 1), ...
    n_rep*n_years, 1); (control(1, 1:n_years, i7))'];
T.ACCase = [reshape(repmat(control(2, 1:n_years, i7), n_rep, 1), ...
    n_rep*n_years, 1); (control(2, 1:n_years, i7))'];
T.PlantDensity = [reshape(Sim(:, 1:n_years, 1), n_rep*n_years, 1);...
    avg_dens];
T.SSplants = [reshape(Sim(:, :, 2), n_rep*n_years, 1);...
    avg_plantsG(:, 1)];
T.RSplants = [reshape(Sim(:, :, 3), n_rep*n_years, 1);...
    avg_plantsG(:, 2)];
T.RRplants = [reshape(Sim(:, :, 4), n_rep*n_years, 1);...
    avg_plantsG(:, 3)];
T.SSseeds = [reshape(Sim(:, :, 8), n_rep*n_years, 1);...
    avg_seedsG(:, 1)];
T.RSseeds = [reshape(Sim(:, :, 9), n_rep*n_years, 1);...
    avg_seedsG(:, 2)];
T.RRseeds = [reshape(Sim(:, :, 10), n_rep*n_years, 1);...
    avg_seedsG(:, 3)];
T.SSrhizomes = [reshape(Sim(:, :, 5), n_rep*n_years, 1);...
    avg_rhizomesG(:, 1)];
T.RSrhizomes = [reshape(Sim(:, :, 6), n_rep*n_years, 1);...
    avg_rhizomesG(:, 2)];
T.RRrhizomes = [reshape(Sim(:, :, 7), n_rep*n_years, 1);...
    avg_rhizomesG(:, 3)];
% write table to text file 
writetable(T, strcat('Table_Stochastic_control', num2str(i7), ...
    '_selfing', num2str(logical(p_self(i1))), ...
    '_sexualReproduction', num2str(seeds(i2)), '_seedbank', ...
    num2str(bank(i3)), '_cost', num2str(1000*cost_seeds(i4)), ...
    '_kCost', num2str(10*k_cost(i5)), '_kHerb', num2str(10*k_herb(i6)),...
    '_fieldsize', num2str(A(i8))));

end
end
end
end
end
end
end
end