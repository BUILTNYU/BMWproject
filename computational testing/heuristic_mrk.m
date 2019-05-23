function [RebalanceDecision,runtime,objective]=heuristic_mrk(Lamda,N,H,adj,dist_tij,theta,dist_h,station_capacity,charg_stations,initial_loc)
big_M = 9999999999999;
demand=sum(Lamda);
NN1=N*H;
mu=zeros(NN1,1)+144;
set_O=initial_loc';
F=size(set_O,1);
set_J=charg_stations';
rho=generate_rho(0.95,0,F);
Cap=station_capacity;
sum_cap=sum(Cap);
excess=F-sum_cap;
linked_dist=zeros(1,NN1)+big_M;
populate=set_O(F);
used_station=0;
obj_heur=0;
markovian_sum=zeros(NN1,1)+demand^2;
servers=zeros(NN1,1);
temp_link=zeros(NN1,NN1);
dem=zeros(NN1,NN1);
linked=zeros(NN1,1);
served=zeros(NN1,1);
margin=zeros(NN1,2);
placed=zeros(NN1,1);
now=tic();
for i=F:-1:1
    rebalance(i)=populate;
    cost=zeros(NN1,1);
    savings=zeros(NN1,1);
    D=digraph(adj);
    layer=floor(set_O(i)/N)+1;
    if excess>0
       tmp1=(layer-1)*N+1;
       tmp2=(layer-1)*N+N;
       excess=excess-1;
    else
       tmp1=1;
       tmp2=NN1;
    end
    init_path=shortestpath(D,set_O(i),rebalance(i));
    init_cost=0;
    for q=1:size(init_path,2)-1
        if init_path(q+1)-init_path(q)>1
           init_cost=init_cost+dist_h;
        else
           init_cost=init_cost+dist_tij(init_path(q),init_path(q+1));   
        end
    end
    init_cost=init_cost*theta;
    temp=populate;
    for n=tmp1:tmp2
         dem(n,:)=served;
         stop=0;
         temp_list(n,:)=linked_dist;
         
         % Markovian candidate point selection %
         if i~=F   
         tmp_servers(temp)=tmp_servers(temp)-1;
         tmp_servers(n)=tmp_servers(n)+1;
         markovian_sum=markov_rhs(rho,tmp_servers,NN1,mu(1));
         temp=n;
         if sum(markovian_sum)<=demand
            savings(n)=-big_M;
            stop=1;
         end
         end
   if placed(n)>0
      if rho(placed(n))*mu(n)>served(n)
         savings(n)=0;
      else
         gain=(rho(placed(n)+1)-rho(placed(n)))*mu(n);
         if served(n)-rho(placed(n))*mu(n)<gain
            gain=served(n)-rho(placed(n))*mu(n);
         end
         savings(n)=(margin(n,1)*gain)/margin(n,2); 
      end
   end
   if stop==0     
        path=shortestpath(D,set_O(i),n);
        for q=1:size(path,2)-1
            if path(q+1)-path(q)>1
               cost(n)=cost(n)+dist_h;
               used_station=path(q+1)-floor(path(q+1)/N);
            if used_station==0
               used_station=N;
            end
            else
               cost(n)=cost(n)+dist_tij(path(q),path(q+1)); 
            end
         end
         cost(n)=cost(n)*theta;
         savings(n)=savings(n)+init_cost-cost(n);
         %--------------------------------------------%    
         for k=tmp1:tmp2
             layr=zeros(1,2);
             if k<=size(set_J,1)
                if used_station==set_J(k)
                   station_capacity(k)=station_capacity(k)-1;
                   if station_capacity(k)==0
                      for s=1:H-1
                          adj(set_J(k)+N*(s-1),set_J(k)+N*s)=0;
                      end
                   end
                end
             end 
             check=[n k];
             for b=1:2
                 if mod(check(b),N)==0
                    layr(b)=floor(check(b)/N);
                 else
                    layr(b)=floor(check(b)/N)+1;
                 end
             end
             if linked_dist(k)>dist_tij(n,k) && layr(1)>=layr(2)   
                savings(n)=savings(n)+(linked_dist(k)-dist_tij(n,k))*Lamda(k);
                temp_link(n,k)=n;
                if i~=F
                dem(n,n)=dem(n,n)+Lamda(k);
                dem(n,linked(k))=dem(n,linked(k))-Lamda(k);
                end
                temp_list(n,k)=dist_tij(n,k);
             else
                temp_link(n,k)=linked(k);
                temp_list(n,k)=linked_dist(k);
             end
        end
      end
    end
     rebalance(i)=max(find(savings(:)==max(savings(tmp1:tmp2))));
     linked=temp_link(rebalance(i),:);
     served=dem(rebalance(i),:);
     placed(rebalance(i))=placed(rebalance(i))+1; 
     linked_dist=temp_list(rebalance(i),:);
     obj_heur=obj_heur+cost(rebalance(i));
     RebalancePath{i}=shortestpath(D,set_O(i),rebalance(i));
     if i==F
        populate=rebalance(i);
        servers(populate)=F;
        served(populate)=demand;
        rebalance=zeros(F,1)+populate;
        markovian_sum=markov_rhs(rho,servers,NN1,mu(1));
        start_grad=savings;
     else
         if i==F-1
            margin(populate,1)=start_grad(populate)-start_grad(rebalance(i));
            margin(populate,2)=demand;
         end
         if placed(rebalance(i))==1
            margin(rebalance(i),1)=savings(rebalance(i));
            margin(rebalance(i),2)=served(rebalance(i));
         end
         servers(populate)=servers(populate)-1;
         servers(rebalance(i))=servers(rebalance(i))+1;
     end
     tmp_servers=servers;
    end
runtime = toc(now);
RebalanceDecision=rebalance;
if excess>0
   error=1;
end

% Calculate objective value
count=zeros(NN1,1);
reb=sort(unique(rebalance));
for i=1:size(rebalance,1)
    count(rebalance(i))=count(rebalance(i))+1;
end

layr=[];
f=[];
points=size(reb,1);
for i=1:points 
    if mod(reb(i),N)==0
       layr(i)=floor(rebalance(i)/N);
    else
       layr(i)=floor(rebalance(i)/N)+1;
    end
    for j=1:NN1
        f=[f,Lamda(j)*dist_tij(reb(i),j)];
    end
end

Aeq=[];
beq=[];
for i=1:NN1
    for j=1:points
        Aeq(i,(j-1)*NN1+i)=1;   
    end
    beq(i,1)=1;
end
temp=[];
temp2=[];
for i=1:points
    for j=layr(i):H-1
        for k=1:N
            temp(k+j*N+(i-1)*NN1,k+j*N+(i-1)*NN1)=1;
            temp2(k+j*N+(i-1)*NN1,1)=0;
        end
    end
end
Aeq=[Aeq;temp];
beq=[beq;temp2];

A=[];
b=[];
for i=1:points
    for j=1:NN1
        A(i,j+(i-1)*NN1)=Lamda(j);
    end
    b(i,1)=rho(count(reb(i)))*mu(reb(i));
end

intcon =[1:NN1*points];
lb = zeros(NN1*points,1);
ub = ones(NN1*points,1);
[x,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
objective=fval+obj_heur

%-------------------- Markovian RHS function ---------------------%
function markov_sum=markov_rhs(rho,servers,NN1,mu)
for j=1:NN1
    if servers(j)>0
       markov_sum(j)=rho(servers(j))*mu;
    end
end
end
%-----------------------------------------------------------------%

%-------------------- Rho Generation Function --------------------%
function rho=generate_rho(alpha,beta,F)
Veh=F;
Cj=Veh;
b=beta;
a=alpha;
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
% end
%----------------------------------------------------------------%
end