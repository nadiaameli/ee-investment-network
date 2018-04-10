clear all
close all

load main_bipartite_network_ws % Load bipartite network and lists of investors
                               % and projects (with associated info)


CP = unique(projects_country); % Uniques list of countries for projects
CI = unique(investors_country); % Unique list of countries for investors

%%% Construction of network

C = sparse(zeros(length(CI),length(CP)));

for i = 1:size(E,1)
   
    country_1 = investors_country(E(i,1));
    country_2 = projects_country(E(i,2));
    
    ind_1 = strfind(CI,country_1);
    ind_1 = find(not(cellfun('isempty',ind_1)));
    
    ind_2 = strfind(CP,country_2);
    ind_2 = find(not(cellfun('isempty',ind_2)));   
    
    C(ind_1,ind_2) = C(ind_1,ind_2) + 1;
    
end

%%% Fraction of weight staying within a country %%%%%%%%%%%%%%%%%%%%%%%%%%%

within_country_weight = 0;

for i = 1:length(CI)
   
    country = CI(i);

    ind = strfind(CP,country);
    ind = find(not(cellfun('isempty',ind)));
    
    within_country_weight = within_country_weight + sum(C(i,ind));
    
end

%%% Validation - Hypergeometric test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.01;

L = length(find(C > 0)); % Number of links == number of tests
str = sum(sum(C)); % Total strength

alpha_country = alpha/L; % Bonferroni threshold for over-expressions

p_over_country = [];

for i = 1:size(C,1)
   
    f = find(C(i,:) > 0); % Find links going out of node i
    
    k = length(f); % Out-degree of node i
    s = sum(C(i,f))'; % Strength of node i (sum of all weights)
    s_i = sum(C(:,f))'; % Incoming strength of tags that user i has interacted with

    pvalues = hygecdf(C(i,f)'-ones(length(f),1),str,s,s_i,'upper');
    
    p_over_country = [p_over_country; pvalues i*ones(length(pvalues),1) f']; % +1 adjusts indexing with tags
    
end

%%% FDR
[su,ind] = sort(p_over_country(:,1));
p_over_country_FDR = [su p_over_country(ind,2) p_over_country(ind,3)];
aux = alpha_country*[1:1:length(su)]';
f = find(su < aux);
f = f(end);
FDR_country = p_over_country_FDR(1:f,:);

%%% Plotting network

big_C = [zeros(size(C,1),size(C,1)) C; C' zeros(size(C,2),size(C,2))];

for i = 1:length(CI)
   names{i} = CI{i}; 
end

for i = length(CI)+1:length(CI)+length(CP)
   names{i} = CP{i-length(CI)}; 
end

[i,j] = find(triu(big_C,1) > 0);
w = big_C(find(triu(big_C,1) > 0));

G = graph(i,j,full(w));

LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);

h = plot(G,'LineWidth',LWidths);

for i = 1:size(C,1)+size(C,2)
   
    labelnode(h,i,names{i})
    
end

%%% Make it bipartite
h.XData(1:size(C,1)) = 1;
h.XData((size(C,1)+1):end) = 2;
h.YData(1:size(C,1)) = linspace(0,1,size(C,1));
h.YData((size(C,1)+1):end) = linspace(0,1,size(C,2));

%%% Make it pretty
nl = h.NodeLabel;
h.NodeLabel = '';
xd = get(h,'XData');
yd = get(h,'YData');
text(xd(1:length(CI)),yd(1:length(CI)),nl(1:length(CI)),'FontSize',12,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','middle')
text(xd(length(CI)+1:length(CI)+length(CP)),yd(length(CI)+1:length(CI)+length(CP)),nl(length(CI)+1:length(CI)+length(CP)),'FontSize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle')

