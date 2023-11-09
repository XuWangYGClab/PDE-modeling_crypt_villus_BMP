% Code name: parameter_screening_grid
% Author: Xu Wang
% Last update: 11/09/2023
%

% 
% %-------------- preparation before screening ----------------
% % %%%%%%%%%% load sparse parameter grid for screening
load('para_085.mat')

%%%%%%%%%% Constant model parameters for model
  %%% nodes and time range
n = 61;   % number of nodes to evaluate finite difference equations
tRange = [0 14400];  % time interval to evaluate differential equations for 2 hours
  %%% geometry
parameters.Ltot = 600;    % length of the crypt-villus axis (1D assumption)       microns total length: 600; villus length: 425; crypt length 100; stroma length: 75
parameters.LI = 75;   % length of inhibitor expression region   microns
%-------------- screening operations ----------------
MBMP_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MI_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MBI_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MBR_WT   = zeros(n,1000000); % initialization of a vecotor storing WT model outputs for all parameter combinations
SSE = zeros(1,1000000);  % initialization of a vecotor storing final SSE's for all parameter combinations
j = 1;
for i=1:85
disp(['i = ', int2str(i)]);
%%%%%%%%%% parameters to be screened
parameters.DB= para_grid(1,i);
parameters.DI=para_grid(2,i);
parameters.DBI=para_grid(3,i);
parameters.konI=para_grid(4,i);
parameters.koffI=para_grid(5,i);
parameters.k1=para_grid(6,i);
parameters.k2=para_grid(7,i);
parameters.decB=para_grid(8,i);
parameters.decI=para_grid(9,i);
parameters.decBR=para_grid(10,i);
parameters.j1= para_grid(11,i);
parameters.j2= para_grid(12,i);
parameters.j3= para_grid(13,i);
%%%%%%%%% ODE solver
% solve for WT

[B1, I1, BI1, BR1] = FiniDiffMod_smallintestine(n, tRange, parameters);
[m,~] = size(B1);
B_WT = B1(m,:);                                     
I_WT = I1(m,:);
BI_WT = BI1(m,:);
BR_WT = BR1(m,:);

%%%%%%%%%% model results scaling

MBMP_WT(:,j)  = B_WT(1:n);
MI_WT(:,j)   = I_WT(1:n);
MBI_WT(:,j)   = BI_WT(1:n);
MBR_WT(:,j)   = BR_WT(1:n);

j=j+1;
end

%-------------- save data ----------------

save('modeldata_4hour_1to10w.mat','MBMP_WT','MI_WT','MBI_WT','MBR_WT');
