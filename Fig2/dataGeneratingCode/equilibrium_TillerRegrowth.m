% Calculates the long-term population composition in an untreated 
% population using the Perron-Frobenius theorem.

%% Setting parameters:

% Ecological:
% Proportion of seedgermination: 
p_germ = 0.3;
% Proportion of bud sprouting (no tillage, tillage):
p_sprout = [0.2, 0.3];

% Fecundity, i.e. number of seeds produced by one plant:
f = 13000; 

% Number of rhizome buds produced per plant:
b = 140;

% Rhizome winter mortality (no tillage, tillage): 
d_rhizomes = [0.35, 0.6];

% Loss and natural mortality of fresh seeds over winter:
d_seeds = 0.94;
% Natural yearly seed mortality in the seedbank:
d_bank = 0.48;

% Evolutionary:
% Mutation rate:
mu = 10^(-8);

% Fitness cost of resistan on seed production: 
cost_seeds = 0.001:0.001:0.4;
% Factor reducing the fitness cost of RS type relative to RR type
k_cost = [0 0.25 0.5 0.75 1];

% Parameters to choose:
% Control strategies:
% Logical value stating whether the soil is tilled at season start:
tillage = 0;
% Adaption of praramters according to the control strategy:
d_rhizomes = d_rhizomes(1 + tillage);
p_sprout = p_sprout(1 + tillage);


% Inheritence matrix:
% Without backmutation
M = [(1 - mu)^2 2 * mu * (1 - mu) mu^2; ...
    1/4 * (1 - mu)^2 1/2 * (1 - mu^2) 1/4 * (1 + mu)^2; 0 0 1];

% Array with R allele fraction depending on the fitness cost on seed
% production (columns) and the fitness cost reduction in heterozygots:
R_frac = zeros(length(k_cost), length(cost_seeds));

% Arrays with genotype fractions depending on the fitness cost on seed
% production (columns) and the fitness cost reduction in heterozygots:
WW_frac = zeros(length(k_cost), length(cost_seeds));
RW_frac = zeros(length(k_cost), length(cost_seeds));
RR_frac = zeros(length(k_cost), length(cost_seeds));

for j = 1:length(cost_seeds) 
    for i = 1:length(k_cost)
        % Diagonal matrices:
        D = diag([1 1-k_cost(i)*cost_seeds(j) 1-cost_seeds(j)]);
        
        X = [(1 - d_rhizomes) * p_sprout * b * eye(3) +...
            (1 - d_seeds)  * p_germ * f * D * M,...
            (1 - d_seeds)  * f * D * M;...
            (1 - d_bank) * (1 - p_germ)  * p_germ * eye(3),...
            (1 - d_bank) * (1 - p_germ) * eye(3)];
        
        % Eigenvalues and eigenvectors 
        [V,D,W] = eig(X);
        [~, index] = sort(diag(D), 'descend');
        % Right eigenvalue 
        v = V(:, index(1));
        % Left eigenvalue:
        w = W(:, index(1));
        % Normalization:
        w = w / (v' * w);
        v = v * sum(w);
        w = w / sum(w);
        
        % R allele frequency:
        R_frac(i, j) = (w(2)+2*w(3))/(2*sum(w(1:3)));
        % Genotype frequencies:
        WW_frac(i, j) = w(1) / sum(w(1:3));
        RW_frac(i, j) = w(2) / sum(w(1:3));
        RR_frac(i, j) = w(3) / sum(w(1:3));
    end
end

% create a table
T = table;
% assign columns to table
T.Cost = repmat(cost_seeds', length(k_cost), 1);
T.Het = reshape(repmat(k_cost, length(cost_seeds), 1),...
    length(k_cost)*length(cost_seeds) , 1);
T.R = reshape(R_frac', length(k_cost)*length(cost_seeds) , 1);
T.WW = reshape(WW_frac', length(k_cost)*length(cost_seeds) , 1);
T.RW = reshape(RW_frac', length(k_cost)*length(cost_seeds) , 1);
T.RR = reshape(RR_frac', length(k_cost)*length(cost_seeds) , 1);

% write table to text file 
writetable(T, strcat('Table_Equilibrium_tillage', num2str(tillage)));