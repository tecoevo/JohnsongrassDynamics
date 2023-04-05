%% Generation of data for control regimes high ALS Efficiency figure
% Stochastic dynamics of a Johnsongrass population over 30 years depending 
% on herbicde application, tillage, seed production, seed bank formation 
% and dominance of resistance allele and associated fitness cost.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters: 
%   n_rep: number of replicates 
n_rep = 10^3;

%   A: field size
A = [10^4, 10^5];
%   n_years: number of years
n_years = 30;

%   dens_seeds: initial seedbank density: 
dens_seeds = 10;
%   dens_rhizomes: initial rhizome density: 
dens_rhizomes = 1;

%   p_self: proportion of selfpollination
p_self = 0.95;

%   control: 3 x (n_years+1) x array of strategies:
%   The first row gives the tillage strategy. The second and third rows
%   indicates the herbicide apllication. Each column corresponds to one 
%   season and the entries are logical values stating whether the soil is 
%   tilled at season start, whether ACCase-inhibitor is applied or whether 
%   ALS-inhibitor is applied, respectively. 
control = [ones(1, n_years+1); zeros(1, n_years+1); ...
                    zeros(1, n_years+1)];
control(:, :, 2) = [zeros(1, n_years+1); ones(1, n_years+1); ...
                    zeros(1, n_years+1)];
control(:, :, 3) = [ones(1, n_years+1); ones(1, n_years+1); ...
                    zeros(1, n_years+1)];
control(:, :, 4) = [zeros(1, n_years+1); ones(1, n_years+1); ...
                    ones(1, n_years+1)];
control(:, :, 5) = [ones(1, n_years+1); ones(1, n_years+1); ...
                    ones(1, n_years+1)];

temp1 = [0, 1];
temp2 = [1, 0];
temp3 = [zeros(1,2), ones(1,1)];
temp4 = [ones(1,2), zeros(1,1)];
temp5 = [zeros(1,5), ones(1,2)];
temp6 = [ones(1,5), zeros(1,2)];
temp7 = [zeros(1,3), ones(1,1)];
temp8 = [ones(1,3), zeros(1,1)];
temp9 = [zeros(1,2), ones(1,2)];
temp10 = [ones(1,2), zeros(1,2)];
temp11 = [zeros(1,3), ones(1,3)];
temp12 = [ones(1,3), zeros(1,3)];
temp13 = [ones(1,1), zeros(1,2)];
temp14 = [0 1 0];
temp15 = [0 1 0 0];
control(:, :, 6) = [ones(1, n_years+1); ...
    [repmat(temp2, 1, floor((n_years+1) / length(temp2))) ...
    temp2(1, 1:mod(n_years+1, length(temp2)))]; ...
    [repmat(temp1, 1, floor((n_years+1) / length(temp1))), ...
    temp1(1, 1:mod(n_years+1, length(temp1)))]];
control(:, :, 7) = [ones(1, n_years+1); ...
    [repmat(temp10, 1, floor((n_years+1) / length(temp10))) ...
    temp10(1, 1:mod(n_years+1, length(temp10)))]; ...
    [repmat(temp9, 1, floor((n_years+1) / length(temp9))), ...
    temp9(1, 1:mod(n_years+1, length(temp9)))]];
control(:, :, 8) = [ones(1, n_years+1); ...
    [repmat(temp12, 1, floor((n_years+1) / length(temp12))) ...
    temp12(1, 1:mod(n_years+1, length(temp12)))]; ...
    [repmat(temp11, 1, floor((n_years+1) / length(temp11))), ...
    temp11(1, 1:mod(n_years+1, length(temp11)))]];
control(:, :, 9) = [ones(1, n_years+1); ...
    [repmat(temp4, 1, floor((n_years+1) / length(temp4))) ...
    temp4(1, 1:mod(n_years+1, length(temp4)))]; ...
    zeros(1, n_years+1)];
control(:, :, 10) = [ones(1, n_years+1); ...
    [repmat(temp4, 1, floor((n_years+1) / length(temp4))) ...
    temp4(1, 1:mod(n_years+1, length(temp4)))]; ...
    [repmat(temp4, 1, floor((n_years+1) / length(temp4))) ...
    temp4(1, 1:mod(n_years+1, length(temp4)))]];
control(:, :, 11) = [ones(1, n_years+1); ...
    [repmat(temp6, 1, floor((n_years+1) / length(temp6))) ...
    temp6(1, 1:mod(n_years+1, length(temp6)))]; ...
    zeros(1, n_years+1)];
control(:, :, 12) = [ones(1, n_years+1); ...
    [repmat(temp6, 1, floor((n_years+1) / length(temp6))) ...
    temp6(1, 1:mod(n_years+1, length(temp6)))]; ...
    [repmat(temp6, 1, floor((n_years+1) / length(temp6))) ...
    temp6(1, 1:mod(n_years+1, length(temp6)))]];

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
for i4 = 2 %1:length(cost_seeds1) 
for i5 = 2 %1:length(cost_seeds2) 
for i6 = 2 %1:length(k_cost1)  
for i7 = 2 %1:length(k_cost2) 

% Initial population composition:
% Read table with genotype frewuencies at eqilibrium
T = readtable('Table_Equilibrium_tillage0.txt');
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

for i1 = 1 %1:length(A)
for i2 = 2 %1:length(seeds) 
for i3 = 2 %1:length(bank) 
for i8 = 2 %1:length(k_herb1) 
for i9 = 2 %1:length(k_herb2) 
for i10 = 6:8 %6:size(control, 3) 

% n_rep x n_years x 28 array with simulation results. 
% Each row corresponds to one replicate. 
% Each column corresponds to one season of the simulation. 
% Pages contain population level variables of plants, rhizomes and 
% seedbank. 
% Page 1 contains the plant densities at the end of the seasons.
% Pages 2-10 contain the absolute genotype (S1S1 S2S2, R1S1 S2S2, R1R1 
% S2S2, S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
% frequencies in plants at the end of the seasons.
% Pages 11-19 contain the absolute genotype (S1S1 S2S2, R1S1 S2S2, R1R1 
% S2S2, S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2)
% frequencies in primary rhizomes of the seasons.
% Pages 20-28 contain the absolute genotype (S1S1 S2S2, R1S1 S2S2, R1R1 
% S2S2, S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
% frequencies in seedbank of the seasons.
Sim = zeros(n_rep, n_years, 28);

% Replicates:
for i = 1:n_rep

    % Initial seedbank:
    % Absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
    % S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
    % in the initial seed bank:
    S0 = multinomrand(dens_seeds * A(i1), ...
        repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
        repelem([1-RR2-RW2; RW2; RR2], 3, 1));
    % Initial rhizomes:
    % Absolute genotype frequencies (S1S1 S2S2, R1S1 S2S2, R1R1 S2S2,
    % S1S1 R2S2, R1S1 R2S2, R1R1 R2S2, S1S1 R2R2, R1S1 R2R2, R1R1 R2R2) 
    % in the initial rhizomes:
    R0 = multinomrand(dens_rhizomes * A(i1), ...
        repmat([1-RR1-RW1; RW1; RR1], 3, 1) .* ...
        repelem([1-RR2-RW2; RW2; RR2], 3, 1));
    % Plant density in presecing season
    dens0 = dens_rhizomes / 0.65 ;

    % gives the dynamics:
    %   P: matrix of absolute genotype frequencies in plants
    %   R: matrix of absolute genotype frequencies in rhizomes
    %   SB: matrix of absolute genotype frequencies in seed bank
    %   P_dens: vector of plant densities
    [P, R, SB, P_dens] = ...
        stochasticDynamics_TwoHerbicides_densityDependance_highE(A(i1), p_self, ...
        S0, R0, dens0, control(2, :, i10), control(3, :, i10), ...
        control(1, :, i10), seeds(i2), bank(i3), n_years, ...
        cost_seeds1(i4), cost_seeds2(i5), k_cost1(i6), k_cost2(i7), ...
        k_herb1(i8), k_herb2(i9));

    % Save simulation results of the current run:
    Sim(i, :, 1) = P_dens';
    Sim(i, : , 2:10) = P';
    Sim(i, :, 11:19) = R(:, 1:n_years)';
    Sim(i, :, 20:28) = SB(:, 1:n_years)';

end

% n_years x 1 vektor of avarage plant densities. Each entry corresponds to
% one season and gives the average plant density over the replicates. 
avg_dens = sum(Sim(:, 1:n_years, 1), 1)' / n_rep;

% n_years x 9 array of avarage absolute genotype frequencies in plants. 
% Each row corresponds to one season. Column 1 contains the numbers of 
% S1S1 S2S2 plants at season start. Column  2 contains the numbers of
% R1S1 S2S2 plants. Column  3 contains the numbers of R1R1 S2S2 plants. 
% Column 4 contains the numbers of S1S1 R2S2 plants. Column 5 contains the 
% numbers of R1S1 R2S2 plants. Column 6 contains the numbers of R1R1 R2S2 
% plants. Column 7 contains the numbers of S1S1 R2R2 plants. Column 8 
% contains the numbers of R1S1 R2R2 plants. Column 9 contains the numbers 
% of R1R1 R2R2 plants.
avg_plantsG = reshape(sum(Sim(:, 1:n_years, 2:10), 1) / n_rep, n_years, 9);

% (n_years+1) x 9 array of avarage absolute genotype frequencies in 
% rhizomes. 
% Each row corresponds to one season. Column 1 contains the numbers of 
% S1S1 S2S2 rhizomes at season start. Column  2 contains the numbers of
% R1S1 S2S2 rhizomes. Column  3 contains the numbers of R1R1 S2S2 rhizomes. 
% Column 4 contains the numbers of S1S1 R2S2 rhizomes. Column 5 contains 
% the numbers of R1S1 R2S2 rhizomes. Column 6 contains the numbers of 
% R1R1 R2S2 rhizomes. Column 7 contains the numbers of S1S1 R2R2 rhizomes. 
% Column 8 contains the numbers of R1S1 R2R2 rhizomes. Column 9 contains 
% the numbers of R1R1 R2R2 rhizomes. 
avg_rhizomesG = reshape(sum(Sim(:, :, 11:19), 1) / n_rep, n_years, 9);

% (n_years+1) x 9 array of avarage absolute genotype frequencies in seeds. 
% Each row corresponds to one season. Column 1 contains the numbers of 
% S1S1 S2S2 seeds at season start. Column  2 contains the numbers of
% R1S1 S2S2 seeds. Column  3 contains the numbers of R1R1 S2S2 seeds. 
% Column 4 contains the numbers of S1S1 R2S2 seeds. Column 5 contains the 
% numbers of R1S1 R2S2 seeds. Column 6 contains the numbers of R1R1 R2S2 
% seeds. Column 7 contains the numbers of S1S1 R2R2 seeds. Column 8 
% contains the numbers of R1S1 R2R2 seeds. Column 9 contains the numbers 
% of R1R1 R2R2 seeds.
avg_seedsG = reshape(sum(Sim(:, :, 20:28), 1) / n_rep, n_years, 9);


% create a table
T = table;
% assign columns to table
T.Season = [reshape(repmat(1:n_years, n_rep, 1), ...
    n_rep*n_years, 1); (1:n_years)'];
T.Run = [repmat((1:n_rep)', n_years, 1); NaN(n_years, 1)];
T.Herb1 = [reshape(repmat(control(2, 1:n_years, i10), n_rep, 1), ...
    n_rep*n_years, 1); (control(2, 1:n_years, i10))'];
T.Herb2 = [reshape(repmat(control(3, 1:n_years, i10), n_rep, 1), ...
    n_rep*n_years, 1); (control(3, 1:n_years, i10))'];
T.Tillage = [reshape(repmat(control(1, 1:n_years, i10), n_rep, 1), ...
    n_rep*n_years, 1); (control(2, 1:n_years, i10))'];
T.PlantDensity = [reshape(Sim(:, 1:n_years, 1), n_rep*n_years, 1);...
    avg_dens];
T.SSSSplants = [reshape(Sim(:, :, 2), n_rep*n_years, 1);...
    avg_plantsG(:, 1)];
T.RSSSplants = [reshape(Sim(:, :, 3), n_rep*n_years, 1);...
    avg_plantsG(:, 2)];
T.RRSSplants = [reshape(Sim(:, :, 4), n_rep*n_years, 1);...
    avg_plantsG(:, 3)];
T.SSRSplants = [reshape(Sim(:, :, 5), n_rep*n_years, 1);...
    avg_plantsG(:, 4)];
T.RSRSplants = [reshape(Sim(:, :, 6), n_rep*n_years, 1);...
    avg_plantsG(:, 5)];
T.RRRSplants = [reshape(Sim(:, :, 7), n_rep*n_years, 1);...
    avg_plantsG(:, 6)];
T.SSRRplants = [reshape(Sim(:, :, 8), n_rep*n_years, 1);...
    avg_plantsG(:, 7)];
T.RSRRplants = [reshape(Sim(:, :, 9), n_rep*n_years, 1);...
    avg_plantsG(:, 8)];
T.RRRRplants = [reshape(Sim(:, :, 10), n_rep*n_years, 1);...
    avg_plantsG(:, 9)];
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
T.SSSSseeds = [reshape(Sim(:, :, 11), n_rep*n_years, 1);...
    avg_seedsG(:, 1)];
T.RSSSseeds = [reshape(Sim(:, :, 12), n_rep*n_years, 1);...
    avg_seedsG(:, 2)];
T.RRSSseeds = [reshape(Sim(:, :, 13), n_rep*n_years, 1);...
    avg_seedsG(:, 3)];
T.SSRSseeds = [reshape(Sim(:, :, 14), n_rep*n_years, 1);...
    avg_seedsG(:, 4)];
T.RSRSseeds = [reshape(Sim(:, :, 15), n_rep*n_years, 1);...
    avg_seedsG(:, 5)];
T.RRRSseeds = [reshape(Sim(:, :, 16), n_rep*n_years, 1);...
    avg_seedsG(:, 6)];
T.SSRRseeds = [reshape(Sim(:, :, 17), n_rep*n_years, 1);...
    avg_seedsG(:, 7)];
T.RSRRseeds = [reshape(Sim(:, :, 18), n_rep*n_years, 1);...
    avg_seedsG(:, 8)];
T.RRRRseeds = [reshape(Sim(:, :, 19), n_rep*n_years, 1);...
    avg_seedsG(:, 9)];
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
T.SSSSrhizomes = [reshape(Sim(:, :, 20), n_rep*n_years, 1);...
    avg_rhizomesG(:, 1)];
T.RSSSrhizomes = [reshape(Sim(:, :, 21), n_rep*n_years, 1);...
    avg_rhizomesG(:, 2)];
T.RRSSrhizomes = [reshape(Sim(:, :, 22), n_rep*n_years, 1);...
    avg_rhizomesG(:, 3)];
T.SSRSrhizomes = [reshape(Sim(:, :, 23), n_rep*n_years, 1);...
    avg_rhizomesG(:, 4)];
T.RSRSrhizomes = [reshape(Sim(:, :, 24), n_rep*n_years, 1);...
    avg_rhizomesG(:, 5)];
T.RRRSrhizomes = [reshape(Sim(:, :, 25), n_rep*n_years, 1);...
    avg_rhizomesG(:, 6)];
T.SSRRrhizomes = [reshape(Sim(:, :, 26), n_rep*n_years, 1);...
    avg_rhizomesG(:, 7)];
T.RSRRrhizomes = [reshape(Sim(:, :, 27), n_rep*n_years, 1);...
    avg_rhizomesG(:, 8)];
T.RRRRrhizomes = [reshape(Sim(:, :, 28), n_rep*n_years, 1);...
    avg_rhizomesG(:, 9)];
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
writetable(T, strcat('Table_TwoHerbicides_high_Stochastic_control', num2str(i10), ...
    '_sexualReproduction', num2str(seeds(i2)), '_seedbank', ...
    num2str(bank(i3)), '_cost1_', num2str(1000*cost_seeds1(i4)), '_cost2_', ...
    num2str(1000*cost_seeds2(i5)), '_kCost1_', num2str(10*k_cost1(i6)), ...
    '_kCost2_', num2str(10*k_cost2(i7)), '_kHerb1_', ...
    num2str(10*k_herb1(i8)), '_kHerb2_', num2str(10*k_herb2(i9)), ...
    '_fieldsize', num2str(A(i1))));

end
end
end
end
end
end
end
end
end
end