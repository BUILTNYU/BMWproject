% This heuristic code is done by Zhexi (Jesse) Fu on Nov15 2018
% This work is modified by adding capacity constraints at charging stations 
% from the preivous heuristic code, which is named as 
% "EV_relo_heuristic_Nov142018.m"
% This version is for the testing as a function.
% This code is the same as "EV_relo_heuristic_Nov192018.m"

% N is about 304 zones
% Veh # TBD
% CS location TBD


% clear all; close all; clc;
% time_start = cputime; 
function [PathFlow,PathFinDec]=HA(Lamda,H,N, set_Ori,set_J,dist_ij,Cap)
%% Parameters

% load Lamda6; %arrival rates at each node-charge i
% Length of Lamda is 24 (namely N*H)
% Lamda = [3.8 2.3 1.2 2.5 ... 4.5 0.25 5.8 2.3]'
%load Lamda100;

%N = 6; % set of nodes
%N= 100;
%H = 4; % charge levels


%set_Ori = [14,7,3]'; % origin nodes of vehs
%set_Ori = [7 8 16 18 20 27 37 38 45 50]';
F = length(set_Ori); % num of vehs


%set_J = [2 6]'; % set of charging station nodes
%set_J=[15 27 36 49]';
NJ = length(set_J); % num of charging station 

set_Relo = (1:N)'; % set of relocation position nodes


%Cap = [2,2];
%Cap = [5,5,5,5];% capacity at charging station
set_Cap = zeros(1,length(set_J)); % occupied # of vehicles at CS

%dist_ij=5;%time between two grid neighbor nodes
dist_h=25;% time between two charge neighbor nodes 

NN1=N*H; % size_i is total number of nodes-charge in the node-charge graph
theta = 0.2;
bigM = 100000000;
   

%% Generate Demand Values
Lamda_Hor = zeros(N,1);
for i=1:NN1
   index =  rem(i,N);
   if index == 0
       index = N;
   end
   Lamda_Hor(index,1) = Lamda_Hor(index,1)+Lamda(i);
end
% Lamda_Hor = [25.9 14.2 15.8 20.25 28.3 18.7]' check


%% Generate Charging Levels
Level=zeros(NN1,1);
for h=1:H
   h_vec =  (h-1)*N+1: (h-1)*N+N;
   Level(h_vec,1)=h;
end
% Level(1:6)=1, Level(7:12)=2, Level(13:18)=3, Level(19:24)=4 check


%% Build the Distance Matrix between Nodes and CHarging Levels 
% dist_t = zeros(N,N); 
% for i=1:N-1
%    for j=i+1:N
%         dist_t(i,j) = (j-i)*dist_ij;
%         dist_t(j,i) = dist_t(i,j);
%    end
% end
% % horizontal distance between nodes in each row (1~N), size of NxN check 
% 
% dist_tij = zeros(NN1,NN1); 
% for i=1:NN1-1
%     for j=i+1:NN1
%         index_i = rem(i,N);
%         index_j = rem(j,N);
%         if index_i == 0
%            index_i = N; 
%         end
%         if index_j == 0
%            index_j = N; 
%         end
%         dist_tij(i,j) = dist_t(index_i,index_j);
%         dist_tij(j,i) = dist_tij(i,j);
%     end
% end
% horizontal distance between positions (1~NN1), size of NN1xNN1
% cost between charging levels is NOT considered here


%% For each vehicle, find current charging level and the charging cost
%  Veh_Charg is in size of (F,1)
Veh_Charg = zeros(F,1); % charging cost
for f=1:F
    Veh_Charg(f,1) = ( H-Level(set_Ori(f)) )*dist_h;
end
% check


%% For each vehicle, find the cost to all CS 
%  Veh_toCS is in size of (F,NJ)
Veh_toCS = zeros(F,NJ);
for i=1:F
    for j=1:NJ
        Veh_toCS(i,j) = dist_tij(set_Ori(i),set_J(j));
    end
end
% check


%% On the fully charged level, relocation distance from each CS
%  dist_fromCS is in size of (NJ,N)
dist_fromCS = zeros(NJ,N);
for i=1:NJ
    for j=1:N
        dist_fromCS(i,j) = dist_tij(set_J(i),j);
    end
end
% check


%% Recursive Loop
% For each vehicle, find the min cost path to fully charged level

% Check the total capacity enought or not
ite_lim = F;
if sum(Cap)<F
    ite_lim = sum(Cap);
    disp('Total Capacity Not Enough')
end

set_Comp = [];
dist_cus = zeros(F,N);
Veh_CostSum = zeros(F*F,NJ*N);
Fin_Relo = zeros(ite_lim,4); % [veh #, cost, CS, relocation position]
pre_relo_cost = 0;
idx_beg = 1;
idx_end = 0;

% Because of total capacity constraint
% Only allocate the # of "sum(Cap)" vehicles
for ite = 1:ite_lim 

    
    % 1. For each final relocation position, find the total customer cost
    for i=1:N
        dist_cus_idx = zeros(1,N);
        for j=1:N
            dist_cus_idx(j) = min( [dist_tij(j,i),dist_tij(j,set_Comp)] )*Lamda_Hor(j);
        end
        dist_cus(ite,i) = sum(dist_cus_idx);
    end
    % dist_cus check
    
    % 2. For each vehicle via different charging station, for each
    %    final relocation position, sum the customer and operator cost.
    for f=1:length(set_Ori)
        idx_f = idx_end+f;
        for j=1:NJ
              
           for i=1:N  
               index = (j-1)*N+i;

               if ismember(i,set_Relo)==1
                   if set_Cap(j)<Cap(j)
                        Veh_CostSum(idx_f,index) = theta*(Veh_toCS(f,j)+Veh_Charg(f)+dist_fromCS(j,i))+dist_cus(ite,i);
                   elseif set_Cap(j)==Cap(j)
                        Veh_CostSum(idx_f,index) = bigM;
                   else
                       error('Capacity Wrong')
                   end
               else
                   Veh_CostSum(idx_f,index) = bigM;
               end 

           end
                 
        end
    end
    
    % 3. Find the min cost among all vehicles 
    % 4. Record the [veh #, cost, CS, relocation position]
    idx_beg = idx_end+1; % 1; 3+1=4; 5+1=6
    idx_end = idx_end+length(set_Ori); % 3; 3+2=5; 5+1=6
    if ite<F
        [veh_val, veh_idx] = min(Veh_CostSum(idx_beg:idx_end,:), [], 2); % min cost for each vehicle
        [tol_val,tol_idx] = min(veh_val); % min cost among all vehicles
        
        Fin_Relo(ite,1) = set_Ori(tol_idx); % veh #
        Fin_Relo(ite,2) = tol_val+pre_relo_cost; % cost
        Fin_Relo(ite,3) = floor( (veh_idx(tol_idx)-1)/N )+1; % CS
        
        % Update the capacity at charging stations
        set_Cap(Fin_Relo(ite,3)) = set_Cap(Fin_Relo(ite,3))+1; 
        
        Fin_Relo(ite,4) = rem(veh_idx(tol_idx),N); % relocation position
        if Fin_Relo(ite,4)==0
            Fin_Relo(ite,4) = N;
        end
    
    elseif ite==ite_lim % the last iteration
        [tol_val,veh_idx] = min(Veh_CostSum(idx_beg:idx_end,:), [], 2);
        tol_idx = 1;
        
        Fin_Relo(ite,1) = set_Ori(1); % veh #
        Fin_Relo(ite,2) = tol_val+pre_relo_cost; % cost
        Fin_Relo(ite,3) = floor( (veh_idx-1)/N )+1; % CS
        
        % Update the capacity at charging stations
        set_Cap(Fin_Relo(ite,3)) = set_Cap(Fin_Relo(ite,3))+1; 
        
        Fin_Relo(ite,4) = rem(veh_idx,N); % relocation position
        if Fin_Relo(ite,4)==0
            Fin_Relo(ite,4) = N;
        end
        
    else
        error('ite')
    end
    
    
    % 5. Update the relocation cost since only customer cost will change
    pre_relo_cost = pre_relo_cost + theta*(Veh_Charg(tol_idx)...
        +Veh_toCS(tol_idx,Fin_Relo(ite,3))...
        +dist_fromCS(Fin_Relo(ite,3),Fin_Relo(ite,4)));
    
    % 6. Update the fixed relocation position
    set_Comp = [set_Comp,Fin_Relo(ite,4)];
    % test: A=[1,2,3;4,5,6];B=[2,3];A(1,B)=[2,3]
    
    % 7. Eliminate necessary parts
     set_Ori(tol_idx) = [];
     set_Relo( find(set_Relo==Fin_Relo(ite,4)) ) = [];
     Veh_Charg(tol_idx) = [];
     Veh_toCS(tol_idx,:) = [];
    
end


%% Computation Time
time_end = cputime;
time_duration = time_end-time_start


% Convert the sequence # of CS to the node # of CS
Fin_Relo(:,3) = set_J( Fin_Relo(:,3) ); 
% Convert the Relocation Position on the fully charged level to the node 
% # on the graph
graph_node = ((H-1)*N+1):H*N; %19:24;
Fin_Relo(:,4) = graph_node(Fin_Relo(:,4));

disp('Fin_Relo shows [veh #, Total Cost, Charging Station, Relocation Position]')
%Fin_Relo


%% Generate the path flows
PathFlow = cell(F,1); % Empty cell here

% From vehicle to CS
for i=1:ite_lim
    index_veh = Fin_Relo(i,1);
    veh_level = Level(index_veh);
    index_CS = Fin_Relo(i,3)+N*(veh_level-1);
    if index_veh < index_CS
        PathFlow(i,1) = {index_veh :1: index_CS};
    elseif index_veh > index_CS
        PathFlow(i,1) = {index_veh :-1: index_CS};
    elseif index_veh == index_CS
        PathFlow(i,1) = {index_veh};
    else
        error('Path Flow Error from Vehicle to Charging Station')
    end
end

% From current CS to fully charged CS
for i=1:ite_lim
    index_CS = Fin_Relo(i,3);
    index_full = (H-1)*N+index_CS;
    
    index_veh = Fin_Relo(i,1);
    veh_level = Level(index_veh);
    index_vehCS1 = Fin_Relo(i,3)+N*veh_level;
    index_vehCS1 = min(index_vehCS1,index_full);
    
    % PF(i,j) is accessing the cell at i,j
    % PF{i,j} is accessing the value in the cell i,j
    PathFlow(i,1) = {[PathFlow{i,1},index_vehCS1:N:index_full]};
    
end

% From fully charged CS to relocation position
for i=1:ite_lim
    if Fin_Relo(i,4)>PathFlow{i,1}(end)
        PathFlow(i,1) = {[PathFlow{i,1},PathFlow{i,1}(end)+1:1:Fin_Relo(i,4)]};
    elseif Fin_Relo(i,4)<PathFlow{i,1}(end)
        PathFlow(i,1) = {[PathFlow{i,1},PathFlow{i,1}(end)-1:-1:Fin_Relo(i,4)]};
    end
end

if ite_lim<F
    for i=1:length(set_Ori)
            ite_idx = ite+1;
            PathFlow(ite_idx,1) = {set_Ori(i,1)};
    end
end
%celldisp(PathFlow)

PathFinDec = zeros(F,1);
for i=1:F
    PathFinDec(i) = PathFlow{i,1}(1,end);  
end



