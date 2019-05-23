% load Lamda
% 
% [RebalancePath,Out]=MIP(Lamda/20,5,6,AA,[2 6],5,[0 2]);

%---------------------------------------------- MIP PROGRAM ------------------------------------------------%

function [RebalanceDecision,runtime2,fval]=MIP(Lamda, layers, nodes, initial_loc, charg_stations, grid_dist,dist_h, station_capacity,theta)

   temp=cputime;
   N=nodes; %set of nodes
   H=layers; % charge levels
   set_O=initial_loc'; % origin nodes of vehicles
   F=size(set_O,1);
   O=F; % n of origins 
   DestinationMatrix=zeros(N,N);
   DestinationMatrix=[0.1 0.1 0.3 0.2 0.2 0.1;0.1 0.1 0.3 0.2 0.1 0.2;0.2 0.1 0.1 0.3 0.2 0.1;0.3 0.2 0.1 0.1 0.1 0.2;0.1 0.3 0.2 0.1 0.1 0.2; 0.1 0.2 0.1 0.3 0.2 0.1]; 
   VehLocation=zeros(size(set_O,1),2); % location of iddle vehicles
   for i=1:size(set_O,1)
   VehLocation(i,1)=i;
   end
   VehLocation(:,2)=set_O;
   set_J=charg_stations'; % set of charging station nodes
   NJ = size(set_J,1); % num of charging station
   yig =ones(F,1);
   Cj=F; % max num of veh at a node
   rho= generate_rho(0.95,0,F); % rho generated based on b=0;alfa=0.95;
   yig_vec = zeros(N*H,1);
   for i=1:size(set_O,1)
   yig_vec(set_O(i))=yig_vec(set_O(i))+1;%n of idle veh at origins
   end
   dist_ij=grid_dist; %time between two grid neighbor nodes
   activate_queueConstr=1; % 1 means non-myopic model, 0 means myopic model 
   NN1=N*H; % size_i is total number of nodes-charge in the node-charge graph
   B=F;
   big_M = 10000;
   Cap = station_capacity; % capacity of charging station : if hetero, change to vec
   mu=zeros(NN1,1)+144; % service rate : 1 customer per 10 minutes

if O~= size(set_O,1)
    error('O~= size(set_O,1) error');
end
if NJ~= size(set_J,1)
    error('NJ~= set_J(set_O,1) error');
end
if H<2
    error('H<2');
end

%%%%%%%
%generate graph 
%%%%%%%

n_arc=0; % n of arcs in the node-charge graph 
MAX_INOUTARCS=8; % max num of outgoing arcs of a node in G  
r_ij = zeros(NN1*MAX_INOUTARCS,1);
arcset_in  = zeros(NN1,MAX_INOUTARCS);
arcset_out = zeros(NN1,MAX_INOUTARCS);
arcset_charging = zeros(NN1*2,1);
count_tract_in  = zeros(NN1  ,1);
count_tract_out = zeros(NN1  ,1);
n_arc_charging=0;

II=[];JJ=[];VV=[];

G_node1=[];G_node2=[];%for visualization
for i=1:H
    for j=1:N-1
     index_i = (i-1)*N+j;   
     index_j = (i-1)*N+j+1; 
     G_node1= [G_node1;index_i]; G_node2= [G_node2;index_j];

     %set up link (i,j)
     n_arc=n_arc+1;
     r_ij(n_arc)= dist_ij;
     count_tract_in(index_j)=count_tract_in(index_j)+1;
     an_index1= count_tract_in(index_j);
     arcset_in(index_j,an_index1) = n_arc;% n_arc is arc id
     
     count_tract_out(index_i)=count_tract_out(index_i)+1;
     an_index2= count_tract_out(index_i);
     arcset_out(index_i,an_index2)= n_arc;
     II=[II;index_i ];JJ=[JJ;index_j];VV=[VV;n_arc];
     
     %set up link (j,i)
     G_node1= [G_node1;index_j]; G_node2= [G_node2;index_i];
     n_arc=n_arc+1;
     r_ij(n_arc)= dist_ij;
     count_tract_in(index_i)=count_tract_in(index_i)+1;
     an_index1= count_tract_in(index_i);
     arcset_in(index_i,an_index1) = n_arc;% n_arc is arc id
     
     count_tract_out(index_j)=count_tract_out(index_j)+1;
     an_index2= count_tract_out(index_j);
     arcset_out(index_j,an_index2)= n_arc;
     II=[II;index_j];JJ=[JJ;index_i];VV=[VV;n_arc];
    end
end

for i=1:NJ % set up links on charging stations
    for h=1:H-1
        index_i = set_J(i,1)+(h-1)*N;%start node
        index_j = set_J(i,1)+ h*N;   %end node
        G_node1= [G_node1;index_i]; G_node2= [G_node2;index_j];
        %set up link (i,j)
        n_arc=n_arc+1;
        r_ij(n_arc)= dist_h; %go up
        count_tract_in(index_j)=count_tract_in(index_j)+1;
        an_index1= count_tract_in(index_j);
        arcset_in(index_j,an_index1) = n_arc;% n_arc is arc id
        count_tract_out(index_i)=count_tract_out(index_i)+1;
        an_index2= count_tract_out(index_i);
        arcset_out(index_i,an_index2)= n_arc;
        II=[II;index_i];JJ=[JJ;index_j];VV=[VV;n_arc];
        n_arc_charging = n_arc_charging+1;
        arcset_charging(n_arc_charging)= n_arc;
  
    end
end
n_arc; % check
link_dic = sparse(II,JJ,VV,NN1,NN1);% store node-link dictionary

%% checked ok

% generate charging levels 
Level=zeros(NN1,1);
for h=1:H
   h_vec= (h-1)*N+1: (h-1)*N+N;
   Level(h_vec,1)=h;
end

%%%%%%%%%%%%%%%%%%%
%check
%%%%%%%%%%%%%%%%%%%%%%
if size(rho,1)~=Cj
    error(size(rho,1)~=Cj);
end
if F<=0
    error(' F<=0');
end

%compute t_ij % geo distance between nodes 
dist_t=zeros(N,N);
for i=1:N-1
    for j=i+1:N
        dist_t(i,j)=(j-i)*dist_ij;
         dist_t(j,i)=dist_t(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dist_tij is dist between nodes-charge
%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_tij= zeros(NN1,NN1);
for i=1:NN1-1
    for j=i+1:NN1
        index_i = rem(i,N);
        index_j = rem(j,N);
        if index_i==0 
            index_i=N;
        end
        if index_j==0 
            index_j=N;
        end
        dist_tij(i,j)= dist_t(index_i,index_j);
        dist_tij(j,i)= dist_tij(i,j);
    end
end
now=tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the MIP problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_subpath=  H*(H-1)/2;

length_X= NN1^2 + NN1*Cj+n_arc+NJ*N_subpath  ;%xij, yjm, wij
init_rowi=zeros(1, length_X);
% eq (1)  objective function
f= init_rowi;

%objecive function
for i=1:NN1
    for j=1:NN1
        index = (i-1)*NN1 +j;
        f(1,index)=dist_tij(i,j)*Lamda(i,1);
    end
end
f(1, NN1*(NN1+Cj)+1 : NN1*(NN1+Cj)+n_arc) = theta*r_ij(1:n_arc);

%constraints
%Aeq=[];beq=[];A=[];b=[];

N_A = NN1*(Cj-1)+NN1+NN1*NN1+O+2*(NN1-O)+N;%larger than necessary
N_Aeq = 2*NN1+1+NN1+N;

%A=zeros(N_A,length_X);
b=zeros(N_A,1);
%Aeq=zeros(N_Aeg,length_X);
beq=zeros(N_Aeq,1);

Row_Aeq=0;
Row_A=0;
%set up sparse matrix
count_index_A=0;
count_index_Aeq=0;
N_MAX= NN1.*NN1+NN1.*(Cj+5)-O+NN1;
I_A = zeros(N_MAX,1);
J_A = zeros(N_MAX,1);
V_A = zeros(N_MAX,1);
I_Aeq = zeros(N_MAX,1);
J_Aeq = zeros(N_MAX,1);
V_Aeq = zeros(N_MAX,1);

%eq (2)  sum Xij=1 

for i=1:NN1
    %rowi=init_rowi;
    Row_Aeq = Row_Aeq+1;
    for j=1:NN1
        if Level(j)>=Level(i) 
           count_index_Aeq = count_index_Aeq+1;
           index = (i-1)*NN1+j;
         %  rowi(index)=1;
           I_Aeq(count_index_Aeq,1) = Row_Aeq;
           J_Aeq(count_index_Aeq,1) = index;
           V_Aeq(count_index_Aeq,1) = 1;
        end
    end
    %Aeq=[Aeq;rowi];
    %beq=[beq; 1];
    beq(Row_Aeq,1)=1;
end

%eq (3) sum Xij=0

for i=N+1:NN1
  % rowi=init_rowi;
   Row_Aeq = Row_Aeq+1;
   for j=1:NN1
       if Level(j)<Level(i)  
         count_index_Aeq = count_index_Aeq+1;
         index = (i-1)*NN1+j;        
         %   rowi(index)=1;
         I_Aeq(count_index_Aeq,1) = Row_Aeq;
         J_Aeq(count_index_Aeq,1) = index;
         V_Aeq(count_index_Aeq,1) = 1;
        end
    end
    %Aeq=[Aeq;rowi];
  % beq=[beq;0];
    beq(Row_Aeq,1)=0;
end
% test
 
%eq 4  Yjm<Yjm_1 N*m
for j=1:NN1
    for m=2:Cj
        %rowi=init_rowi;
         Row_A = Row_A +1;
         
         count_index_A = count_index_A+1;
         index1= NN1*NN1 + (j-1)*Cj+m-1;
        %rowi(index1)=-1;
         I_A(count_index_A,1) = Row_A;
         J_A(count_index_A,1) = index1;
         V_A(count_index_A,1) = -1;
       
        %rowi(index1+1)=1;
         count_index_A = count_index_A+1;
         I_A(count_index_A,1) = Row_A;
         J_A(count_index_A,1) = index1+1;
         V_A(count_index_A,1) = 1;
        %A=[A;rowi];
        %b=[b; 0];
        b(Row_A,1)=0;
    end
end
%size(A,1) 
%size(b,1)

% eq 5 sum lamda Xij<=mu(Yjm*rho+sum(Yjm*(rho_m-ro_m_1)
if activate_queueConstr==1
   for j=1:NN1
     %rowi=init_rowi;
     Row_A = Row_A +1;
     for i=1:NN1
        index1=(i-1)*NN1+j;
        %rowi(index1)=lamda(i);
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index1;
        V_A(count_index_A,1) = Lamda(i);
    end
    %m loop
    index2= NN1*NN1+(j-1)*Cj+1;
    %rowi(index2)= -mu(j).*rho(1);
     count_index_A = count_index_A+1;
     I_A(count_index_A,1) = Row_A;
     J_A(count_index_A,1) = index2;
     V_A(count_index_A,1) = -mu(j).*rho(1);
    for m=2:Cj 
        index= NN1*NN1+(j-1)*Cj+m;
       % rowi(index)= mu(j).*(-rho(m)+rho(m-1));
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index;
        V_A(count_index_A,1) = mu(j).*(-rho(m)+rho(m-1));
    end
    %A=[A;rowi];
    %b=[b;0];
     % A(Row_A,:)=rowi;
      b(Row_A,1)=0;
   end
end    
    
% eq (6) sum Yjm=B

%rowi=init_rowi;
Row_Aeq = Row_Aeq+1;
for j=1:NN1
    for m=1:Cj
        index_x= NN1*NN1+(j-1)*Cj+m;
        %rowi(index_x)=1;
        count_index_Aeq = count_index_Aeq+1;
        I_Aeq(count_index_Aeq,1) = Row_Aeq;
        J_Aeq(count_index_Aeq,1) = index_x;
        V_Aeq(count_index_Aeq,1) = 1;
    end
end
%Aeq=[Aeq;rowi];
%beq=[beq; B];
beq(Row_Aeq,1)=B;


% eq (7) Xij<=Yjm1
for i=1:NN1
    for j=1:NN1
        %rowi=init_rowi;
         Row_A = Row_A+1;
        index1= (i-1)*NN1+j;
        index2=NN1*NN1+(j-1)*Cj+1;
        %rowi(index1)=1; 
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index1;
        V_A(count_index_A,1) = 1;
        %rowi(index2)=-1;
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index2;
        V_A(count_index_A,1) = -1;
        %A=[A;rowi];
        %b=[b; 0];
         b(Row_A,1)=0;
    end
end
%size(A,1) 
%size(b,1)

% eq (8) sum -inflow + sum outflow <= yjm*bigM
ss=1:NN1;
ss(set_O)=[];
for i=1:size(ss,2)
    %rowi=init_rowi;
    Row_A = Row_A +1;
    nodeid= ss(i);
    %incoming flow
    for j=1:MAX_INOUTARCS
        if arcset_in(nodeid,j) > 0
           linkid= arcset_in(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
           %rowi(index)= 1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = 1;
        else
            break;
        end
    end
     %outgoing flow
     for j=1:MAX_INOUTARCS
        if arcset_out(nodeid,j) > 0
           linkid= arcset_out(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
           %rowi(index)= -1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = -1;
        else
            break;
        end
     end    
    index2= NN1*NN1+ (nodeid-1)*Cj+1;
    %rowi(index2) = - big_M ;%big_M
    count_index_A = count_index_A+1;
    I_A(count_index_A,1) =  Row_A  ;
    J_A(count_index_A,1) =  index2 ;
    V_A(count_index_A,1) = -big_M  ;
   
    %A=[A;rowi];
    %b=[b;  0];
    b(Row_A,1)=0;
end


%eq (9) - (sum inflow - sum outflow) < = yjm*bigM
for i=1:size(ss,2)
   % rowi=init_rowi;
    Row_A = Row_A +1;
    nodeid= ss(i);
    %incoming flow
    for j=1:MAX_INOUTARCS
        if arcset_in(nodeid,j) > 0
           linkid= arcset_in(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
           %rowi(index)= -1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = -1;
        else
            break;
        end
    end
     %outgoing flow
     for j=1:MAX_INOUTARCS
        if arcset_out(nodeid,j) > 0
           linkid= arcset_out(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
          % rowi(index)= 1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = 1;
        else
            break;
        end
     end    
    index2= NN1*NN1+ (nodeid-1)*Cj+1;
    %rowi(index2) = - big_M;%big_M
    count_index_A = count_index_A+1;
    I_A(count_index_A,1) =  Row_A  ;
    J_A(count_index_A,1) =  index2 ;
    V_A(count_index_A,1) = -big_M  ;
    %A=[A;rowi];
    %b=[b;  0];
    b(Row_A,1)=0;
end
%size(A,1) 
%size(b,1)

% eq (10) sum Wij = sum Yjm
for j=1:NN1
    %rowi=init_rowi;
    Row_Aeq = Row_Aeq+1;
    %incoming flow
    for i=1:MAX_INOUTARCS
       if arcset_in(j,i) > 0
           linkid= arcset_in(j,i);
           index = NN1*NN1+NN1*Cj+linkid;
          % rowi(index)= 1;
           count_index_Aeq = count_index_Aeq+1;
           I_Aeq(count_index_Aeq,1) = Row_Aeq;
           J_Aeq(count_index_Aeq,1) = index;
           V_Aeq(count_index_Aeq,1) = 1;
        else
            break;
        end
    end
     %outgoing flow
    for i=1:MAX_INOUTARCS
        if arcset_out(j,i) > 0
           linkid= arcset_out(j,i);
           index = NN1*NN1+NN1*Cj+linkid;
          %rowi(index)= -1;
           count_index_Aeq = count_index_Aeq+1;
           I_Aeq(count_index_Aeq,1) = Row_Aeq;
           J_Aeq(count_index_Aeq,1) = index;
           V_Aeq(count_index_Aeq,1) = -1;
        else
            break;
        end
    end 
        
    for m=1:Cj
        index = NN1*NN1+ (j-1)*Cj+m;
        %rowi(index)= -1;
         count_index_Aeq = count_index_Aeq+1;
         I_Aeq(count_index_Aeq,1) = Row_Aeq;
         J_Aeq(count_index_Aeq,1) = index;
         V_Aeq(count_index_Aeq,1) = -1;
    end
    %Aeq=[Aeq;rowi];
    %beq=[beq; -yig_vec(j) ];
    beq(Row_Aeq,1)= -yig_vec(j);
end
%size(Aeq,1) 
%size(beq,1)

% eq (11) sum prs =wij
id_arc_charge=0;
for i=1:NJ % set up links on charging stations
    for arc_ch=1:H-1 % arc
        Row_Aeq = Row_Aeq+1;
        index_p=0;
        g=arc_ch;
        h=arc_ch+1;
        id_arc_charge=id_arc_charge+1;
        for gg=1:H-1
            for hh=gg+1:H
                index_p = index_p+1;
                if gg<=g && hh>=h
                   index_i = NN1^2 + NN1*Cj + n_arc + (i-1)*N_subpath+index_p; %N_subpath = H(H-1)/2
                   count_index_Aeq = count_index_Aeq+1;
                   I_Aeq(count_index_Aeq,1) = Row_Aeq;
                   J_Aeq(count_index_Aeq,1) = index_i;
                   V_Aeq(count_index_Aeq,1) = 1;
                end
            end
        end
        link_id = arcset_charging(id_arc_charge);
        index_w = NN1^2 + NN1*Cj+ link_id;
        count_index_Aeq = count_index_Aeq+1;
        I_Aeq(count_index_Aeq,1) = Row_Aeq;
        J_Aeq(count_index_Aeq,1) = index_w;
        V_Aeq(count_index_Aeq,1) = -1;
        beq(Row_Aeq,1)=0;
    end
end

% eq (12) % sum pij <=uj 
for i=1:NJ % 
    Row_A = Row_A+1;
    index_p=0;        
    for g=1:H-1
        for h=g+1:H
            index_p =index_p+1;
            index_i = NN1^2 + NN1*Cj + n_arc + (i-1)*N_subpath+index_p; 
            count_index_A = count_index_A+1;
            I_A(count_index_A,1) = Row_A;
            J_A(count_index_A,1) = index_i;
            V_A(count_index_A,1) = 1;
        end
    end
    b(Row_A,1)=Cap(i); %modified
end


intcon =[1:length_X ];
lb = zeros(length_X,1);
ub = [ones(NN1*(NN1+Cj),1); Inf*ones(length_X-NN1*(NN1+Cj),1)];

I_A=I_A(1:count_index_A,1);
J_A=J_A(1:count_index_A,1);
V_A=V_A(1:count_index_A,1);

I_Aeq=I_Aeq(1:count_index_Aeq,1);
J_Aeq=J_Aeq(1:count_index_Aeq,1);
V_Aeq=V_Aeq(1:count_index_Aeq,1);

A   = sparse(I_A,J_A,V_A,      Row_A, length_X);
Aeq = sparse(I_Aeq,J_Aeq,V_Aeq, Row_Aeq, length_X);
b=b(1:Row_A,1);
beq=beq(1:Row_Aeq , 1);

t=cputime;
disp('eq matrix preparation');
diff =t -temp;
temp=t;
length_X; %check

size(A,1)+size(Aeq,1);%check

options = optimoptions(@intlinprog,'MaxNodes',10000,'MaxTime',20000);
[x,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);

t=cputime;
%disp('exec time')
diff =t -temp;
fval;

% get solution
if exitflag  ~= 1 
    disp('no point satisfies the constraints');
else 
Xij=zeros(NN1,NN1);%
Yjm=zeros(NN1,Cj);
Wij=zeros(n_arc,1);%node based, sparse matrix
Pij=zeros(NJ*N_subpath,1);% subpath flow
for i=1:NN1
    for j=1:NN1
        index= (i-1)*NN1+j;
        Xij(i,j)=x(index);
    end
end
Xij;

for j=1:NN1
    for m=1:Cj
        index= NN1*NN1+(j-1)*Cj+m;
        Yjm(j,m) = x(index);
    end
end
Yjm;

for i=1:n_arc
    index = NN1*NN1+NN1*Cj+i;
    Wij(i)=x(index);
end


for i=1:NJ*N_subpath
    index = NN1*NN1+NN1*Cj+n_arc+i;
    Pij(i,1)=x(index);
end
Wij;
z=f*x;
Yjm ;
Pij;
end

tmp=0;
n=1;
%Calculate rebalanced vehicles
RebalanceDecision=zeros(F,1);
for k=1:size(Yjm,2)
    if tmp==1
        break;
    end
for i=1:N*H
    if Yjm(i,k)~=0
        RebalanceDecision(n,1)=i;
        if n==size(set_O,1)
            tmp=1;
        end
        n=n+1;
    end
    if tmp==1
        break
    end
end
end
% ---------decomposition of link flows (Ted 09/29/2018)-----------

tmp=1;
 mark=-1;
for i=1:(N-1)*H*2  
    if mod((i-1),2*(N-1))==0
        mark=mark+1;
    else
        mark=mark+0;
    end
 if Wij(i)>0 
    if mod(i,2)==0
        indx=1;
    else 
        indx=0;
    end
    arc(tmp,1)=i-(i-1+indx)/2+mark+indx;
    if indx~=0
        arc(tmp,2)=arc(tmp,1)-1;
    else
        arc(tmp,2)=arc(tmp,1)+1;
    end
    arc(tmp,3)=Wij(i);
    tmp=tmp+1;
end
end

for j=1:size(set_J,1)
    count=set_J(j);
   for i=(N-1)*H*2+1+(j-1)*(H-1):(N-1)*H*2+(H-1)+(j-1)*(H-1)
    arc(tmp,1)=count;
    arc(tmp,2)=arc(tmp,1)+N;
    count=arc(tmp,2);
    arc(tmp,3)=Wij(i);
    tmp=tmp+1;
   end
end

dest=RebalanceDecision;
q=0;

for e=1:F   
    
% create adjacency matrix 
adj=zeros((N-1)*H*2,(N-1)*H*2);
for i=1:size(arc,1)
    if arc(i,3)>0.1 
    adj(arc(i,1),arc(i,2))=1;
    end
end
D = digraph(adj);
io=2;

for k=1:F 
    if io==1
  break;
    end
    
if io==0
else
    
for j=1:q
  if list(j)==k
     io=0;
  end
end

if io==0
    io=2;
else
    
P=shortestpath(D,set_O(e),dest(k));
if isempty(P)      
else    
    io=1;
     q=q+1;
    list(q)=k;
    M{e}=P;
for i=1:size(P,2)-1
  for j=1:size(arc,1)
      if P(i+1)==arc(j,2)&& P(i)==arc(j,1)
          arc(j,3)=arc(j,3)-1;
      end
  end
end
end
end
end
end
end

RebalancePath=M;
runtime2=toc(now);
%-----------------------------------------------------%

% FUNCTION TO FIND ALL THE POSSIBLE PATHS FROM A SOURCE NODE TO SINK NODE
% PathFinder(B,StartNode,EndNode)
% B is an Nx2 matrix, where N is the number of Edges in the Graph. The data
% is in the form of 'From Node' to 'To Node'.
% StartNode is the source node, and EndNode is the Sink Node.
% Limitation: Works good till N=20. Also as N increses, execution time also increases.

% By- Abhishek Chakraborty
% Dt: 01-May-2010
% For suggestions and queries, please contact the author at: abhishek.piku@gmail.com

function [costs] = mdijkstra(A,C)
%
%A=square matrix (either adjecancy or cost)
%
%if C=1 then A=adjecancy matrix
%    where, element(i,j)=1 when vertex v is directly connected with j
%    else (i,j)=0
%
%if C=2 then A=cost matrix 
%    where, element (i,j) represents positive integer representing cost
%    between vertex i and j
%
% Output: [costs]: calculated cost matrix
% Developed by: Bharat Patel
% Release date: 03/28/2009

n = length(A);
costs=single(A);
costs(costs==0)=Inf;
for k = 1:n
    disp(sprintf('%d/13-%d/%d',C,k,n));
    w_col=single(costs(:,k));
    if C==1
        mx=max(w_col(find(w_col~=Inf)));
        pmx=0;
        while(mx~=pmx)
            cols=find(w_col==mx);
            tmp=min(costs(:,cols),[],2);
            tmp=tmp+mx;
            tmp1=[w_col tmp];
            w_col=min(tmp1,[],2);
            pmx=mx;
            mx=max(w_col(find(w_col~=Inf)));
        end
        costs(:,k)=w_col;
    elseif C==2
        m1=(w_col(find(w_col)));
        m=sort(unique(m1),'ascend');
        for j=1:length(m)
            mx=m(j);
            cols=find(w_col==mx);
            tmp=min(costs(:,cols),[],2);
            tmp=tmp+mx;
            tmp1=[w_col tmp];
            w_col=min(tmp1,[],2);
        end
        costs(:,k)=w_col;
        costs(k,:)=w_col';
    end
end
for k=1:n
    costs(k,k)=0;
end
end

for i=1:size(RebalancePath,2)
    if isempty(RebalancePath{i})
        RebalancePath{i}=RebalanceDecision(i)
    end
end
%-------------------- Rho Generation Function----------------------------%
function rho=generate_rho(alpha,beta,F)
Veh=F;
Cj=Veh;

b=beta; a=alpha;
rho=zeros(Cj,1);
for m=1:Cj 
 [rho(m,1),~]= find_rho(m,b,a);
end
rho;
%save 'rho_6.txt' rho -ascii;
end

function [x fval] = find_rho(m,b,a)
INFINI=1e+20;
x0 = [0.01,INFINI];
f = @(x) find_rho_1(x,m,b,a);
f(x0);
[x,fval] = fzero(f,x0);
end

function y=find_rho_1(x, m, b ,a)
sum_temp=0;
for k=0:m-1
sum_temp = sum_temp + ((m-k)*factorial(m)*m^b)/factorial(k)*(1./(x.^(m+b+1-k)));
end
y = sum_temp - 1/(1-a);
end

%--------------------------------------------------------------------------------------------%
end
%--------------------------------------------------------------------------------------------%
