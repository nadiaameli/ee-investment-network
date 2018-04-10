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

%%% False Discovery Rate correction
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

for i = length(CI)+1:length(CI)+length(SP)
   names{i} = SP{i-length(CI)}; 
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
text(xd(length(CI)+1:length(CI)+length(SP)),yd(length(CI)+1:length(CI)+length(SP)),nl(length(CI)+1:length(CI)+length(SP)),'FontSize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle')

%%% Nestedness

aux = double(C > 0);

k_out = sum(aux');
k_in = sum(aux);

[s,ind_out] = sort(k_out,'descend');
[s,ind_in] = sort(k_in,'descend');

figure(2)
colormap(winter)
imagesc(aux(ind_out,ind_in))

