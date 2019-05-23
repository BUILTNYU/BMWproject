clear all;
N=6;
H=5;
charg_stations=[2 6];
grid_dist=5;
dist_h=20;
theta=0.02;
load Lamda6;
station_capacity=[2 4];
initial_loc=[26 27 28 29 12 5];
[adj,dist_tij]=network1(H,N,charg_stations,grid_dist,dist_h);
[RebalanceDecision,runtime,objective]=heuristic_mrk(Lamda,N,H,adj,dist_tij,theta,dist_h,station_capacity,charg_stations,sort(initial_loc,'descend'));
% [RebalanceDecision2,runtime2,objective2]=MIP(Lamda, H, N, initial_loc, charg_stations, grid_dist,dist_h, station_capacity,theta);
% Opt_gap=(objective-objective2)*100/objective2;
% faster=(runtime-runtime2)*100/runtime2;