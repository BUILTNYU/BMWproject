function [adj,dist_tij]=network(layers,nodes,charg_stations,grid_dist,dist_h)
H=layers;
N=nodes;
NN1=N*H;
DestinationMatrix=zeros(N,N);
set_J=charg_stations'; % set of charging station nodes
NJ= size(set_J,1); % num of charging station
dist_ij=grid_dist; % time between two grid neighbor nodes
NN1=N*H; % size_i is total number of nodes-charge in the node-charge graph
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
     
     % set up link (i,j)
     n_arc=n_arc+1;
     r_ij(n_arc)= dist_ij;
     count_tract_in(index_j)=count_tract_in(index_j)+1;
     an_index1= count_tract_in(index_j);
     arcset_in(index_j,an_index1) = n_arc; % n_arc is arc id
     count_tract_out(index_i)=count_tract_out(index_i)+1;
     an_index2= count_tract_out(index_i);
     arcset_out(index_i,an_index2)= n_arc;
     II=[II;index_i];JJ=[JJ;index_j];VV=[VV;n_arc];
     
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
        index_j = set_J(i,1)+ h*N;%end node
        G_node1= [G_node1;index_i]; G_node2= [G_node2;index_j];%set up link (i,j)
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
dist_tij=zeros(NN1,NN1);
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

adj=full(link_dic);
for i=1:NN1
    for j=1:NN1
        if adj(i,j)>0
           adj(i,j)=1;
           if abs(i-j)>1
               adj(j,i)=1;
           end
        end
    end
end
end
