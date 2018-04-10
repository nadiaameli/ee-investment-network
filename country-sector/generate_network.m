clear all
close all

load main_bipartite_network_ws % Load bipartite network and lists of investors
                               % and projects (with associated info)


SP = unique(projects_subsector); % Uniques list of sectors for projects
CI = unique(investors_country); % Unique list of countries for investors

%%% Construction of network

C = sparse(zeros(length(CI),length(SP)));

for i = 1:size(E,1)
   
    country = investors_country(E(i,1));
    sector = projects_subsector(E(i,2));
    
    ind_1 = strfind(CI,country);
    ind_1 = find(not(cellfun('isempty',ind_1)));
    
    ind_2 = strfind(SP,sector);
    ind_2 = find(not(cellfun('isempty',ind_2)));   
    
    C(ind_1,ind_2) = C(ind_1,ind_2) + 1;
    
end

C = full(C);

%%% Notation as in De Gangi - Lillo paper
%%% (https://arxiv.org/pdf/1509.00607.pdf)

N = size(C,1); % Number of rows
K = size(C,2); % Number of columns

Y = double(C > 0); % Binarized matrix

out_deg = sum(Y')'; % Out-degree
in_deg = sum(Y)'; % In-degree

out_str = sum(C')';
in_str = sum(C)';

V = full([out_str; out_deg; in_str; in_deg]); % Vector of knowns

%%% Fitnesses

load g_ws % Loads workspace from optimization, where the parameters to
          % generate the null models have been computed

phi = g(1:N); lambda = -log(phi);
psi = g(N+1:2*N); rho = -log(psi);
xi = g(2*N+1:2*N+K); eta = -log(xi);
gamma = g(2*N+K+1:end); delta = -log(gamma);

%%% Monte Carlo to generate null models

N_mods = 1000; % Number of null models
N_sweeps = (1000+2)*N_mods; % Number of iterations

% Matrices where the in / out strenghts and degrees of all nodes are stored
% for each null model
out_str_null = []; in_str_null = [];
out_deg_null = []; in_deg_null = [];
top_eig = []; % Largest eigenvalue (needed for nestedness analysis)

M = ones(N,K); 

nm = 0;

for ns = 1:N_sweeps
    
    for moves = 1:N*K
        
        n = randi(N); k = randi(K);
        
        %%% Throwing coin to decide whether to add or subtract unit weight
        if rand > 1/2
            tmp = M(n,k) + 1;
        else
            tmp = M(n,k) - 1;
            if tmp < 0
                tmp = 0;
            end
        end
        
        ind_tmp = 0; ind = 0;
        
        if tmp > 0
            ind_tmp = 1;
        end
        if M(n,k) > 0
            ind = 1;
        end
        
        %%% Change in Hamiltonian
        DeltaH = (tmp - M(n,k))*(lambda(n)+eta(k)) + ...
               + (ind_tmp - ind)*(rho(n)+delta(k));
        
        if rand < min(1,exp(-DeltaH))
           M(n,k) = tmp; 
        end
           
    end
    
    if mod(ns,1000) == 0 & ns > 1000
        
       out_str_null = [out_str_null sum(M')']; 
       in_str_null = [in_str_null sum(M)'];
       
       aux = double(M > 0);
       
       out_deg_null = [out_deg_null sum(aux')'];
       in_deg_null = [in_deg_null sum(aux)'];
       
       big_M = [zeros(size(M,1),size(M,1)) M; M' zeros(size(M,2),size(M,2))];
       
       top_eig = [top_eig; max(eig(big_M))];
       
       nm = nm+1
               
    end
    
end
    
