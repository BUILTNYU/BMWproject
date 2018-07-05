% Mixed integer program for EV rebalancing (4 layer example)%

clear all;
clear variables;
clc;

lamdas=[3.8 2.3 1.2 2.5 1.1 0.3; 9.1 0.5 8.3 6.1 9.8 9.4; 10 6 11.4 11.6 1.5 6.7; 3 5.4 4.5 0.3 5.8 2.3]';
theta=0.2;

N=6;
G=4;
x=zeros(N,1);
y=[0,5,10,15,20,25];

idle_loc_layer=[1,2,3];
car_paths=zeros(4,3);
reloc_paths=0;

for g=1:G
for i=1:size(idle_loc_layer,2)  
        car_paths(g,i)=0;
    if ( idle_loc_layer(i)>=g)
reloc_paths=reloc_paths+N;
car_paths(g,i)=car_paths(g,i)+N;
    else
      reloc_paths=reloc_paths+2*N;  
      car_paths(g,i)=car_paths(g,i)+2*N;
    end
    end
end

for i=1:N
     for j=1:N
           distance(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
     end
end

t=1;
for h=1:G
    for j=1:N
        for g=1:G
            for i=1:N
                f(t)=lamdas(i,g)*distance(i,j);
                t=t+1;
            end
        end
    end
end

% r=[10 5 0 5 10 15 35 30 35 40 45 40 60 55 60 65 70 65 85 80 85 90 95 90 0 5 10 15 20 25 0 5 10 15 20 25 35 30 35 40 45 50 60 55 60 65 70 75 5 0 5 10 15 20 5 0 5 10 15 20 5 0 5 10 15 20 30 25 30 35 40 55];

r=[10 5 0 5 10 15 35 65 30 60 35 55 40 50 45 45 50 40 60 90 55 85 60 80 65 75 70 70 75 65 85 115 80 110 85 105 90 100 95 95 100 90 0 5 10 15 20 25 0 5 10 15 20 25 35 75 30 70 35 65 40 60 45 55 50 50 60 100 55 95 60 90 65 85 70 80 75 75  5 0 5 10 15 20 5 0 5 10 15 20 5 0 5 10 15 20 30 70 25 65 30 60 35 55 40 50 45 45]


Z=zeros(72,1);
f=[f,theta*r,Z'];

temp=zeros(1,N*N*G*G+3*G*N+reloc_paths);
Aeq=temp;


for g=1:G
    for i=1:N
        for h=1:G
            for j=1:N
                if  (h>=g)
                    temp(((g-1)*N+i-1)+1+((h-1)*N+j-1)*N*G)=1;
                end
            end
        end
        Aeq=[Aeq;temp];
        temp=zeros(1,N*N*G*G+3*G*N+reloc_paths);
    end
    
end      %x=1

% temp=zeros(1,N*N*G*G+3*G*N+reloc_paths));
% Aeq1=temp;

for g=1:G
    for i=1:N
        for h=1:G
            for j=1:N
                if  (h<g)
                    temp(((g-1)*N+i-1)+1+((h-1)*N+j-1)*N*G)=1;
                end
            end
        end
        Aeq=[Aeq;temp];
        temp=zeros(1,N*N*G*G+3*G*N+reloc_paths);
    end
    
end      %x=0

tempy=zeros(72,24)';
for j=1:24
    for k=1:size(idle_loc_layer,2)
    tempy(j,size(idle_loc_layer,2)*(j-1)+k)=-1;
    end
end    

car_short=zeros(1,3);
for k=1:3
for q=1:4
    car_short(1,k)=car_short(1,k)+car_paths(q,k);
end
end

tempx=zeros(G*N,reloc_paths);
n=zeros(3,1)
 for g=1:G
    for j=1:N
        sum=0;
        for k=1:size(idle_loc_layer,2)          
            if (car_paths(g,k)==N)
        tempx(j+(g-1)*N,j+(g-1)*N+sum)=1;
            else   
        tempx(j+(g-1)*N,(j)+(g-1)*N+sum+n(k))=1;
        tempx(j+(g-1)*N,(j+1)+(g-1)*N+sum+n(k))=1;
             n(k)=n(k)+1;
            end
        sum=sum+car_short(k);
        end
   end
 end

 temp=[zeros(24,576),tempx,tempy];
 Aeq=[Aeq;temp];
% Y=W

temp=zeros(3,N*N*G*G+3*G*N+reloc_paths);
sum=0;
for k=1:3
    tmp=sum;
sum=sum+car_short(k)
    for i=N*N*G*G+1+tmp:N*N*G*G+sum
        temp(k,i)=1;        
    end
   
end
 Aeq=[Aeq;temp];
% W=1

A=zeros(1,N*N*G*G+3*G*N+reloc_paths);
temp=zeros(1,N*N*G*G+3*G*N+reloc_paths);

for i=N*N*G*G+reloc_paths+1:N*N*G*G+reloc_paths+72
    temp(i)=1;
end   %sum y=B

Aeq=[Aeq;temp];


temp=zeros(1,N*N*G*G+3*G*N+reloc_paths);
for i=N*N*G*G+reloc_paths+1:N*N*G*G+reloc_paths+24
    for k=1:2
        temp((i-1-N*N*G*G-reloc_paths)*3+N*N*G*G+reloc_paths+1+k-1)=-1;
        temp((i-1-N*N*G*G-reloc_paths)*3+N*N*G*G+reloc_paths+2+k-1)=1; 
        A=[A;temp];
        temp=zeros(1,N*N*G*G+3*G*N+reloc_paths);
    end  
end   %sum y_m<=y_m-1

for h=1:G
    for j=1:N
        for g=1:G
            for i=1:N
                    temp((h-1)*N*G*N+(j-1)*G*N+(g-1)*N+i)=1;
                    temp(N*N*G*G+reloc_paths+((h-1)*N+j-1)*3+1)=-1;
                    A=[A;temp];
                    temp=zeros(1,N*N*G*G+3*G*N+reloc_paths);
            end
        end
        
    end
end      % sum x<=y

temp=zeros(2,reloc_paths);
s1=zeros(3,1);
s2=zeros(3,1);
 for g=1:G
     sum=0;
        for k=1:size(idle_loc_layer,2)  
             s2(k)=s2(k)+car_paths(g,k);
            if (car_paths(g,k)==2*N)        
                for l=(s1(k)+1+sum):2:s2(k)+sum
        temp(1,l)=1;
        temp(2,l+1)=1;
                end      
            end               
             sum=sum+car_short(k); 
           s1(k)=s1(k)+car_paths(g,k); 
        end
      
 end
 
 temp=[zeros(2,576),temp,zeros(2,72)];
 A=[A;temp];

lb=zeros(1,N*N*G*G+3*G*N+reloc_paths);
ub=ones(1,N*N*G*G+3*G*N+reloc_paths);
intcon=1:N*N*G*G+reloc_paths;

beq=ones(24,1); 
beq=[beq;zeros(24,1)];
beq=[beq;zeros(24,1)];
beq=[0;beq;1;1;1;3];

b=zeros(49+576,1);
b=[b;1;1]

[x,fval,exitflag]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);