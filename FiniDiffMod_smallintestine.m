% function name: FiniDiffMod_smallintestine.m
% author: Xu Wang 
%function input:
%    n - number of nodes to evaluate finite difference equations
%    tRange - time interval to evaluate differential equations, eg. [0 5000]
%    parameters - a stucture storing input parameter values like diffusion
%                 rates, and reaction rates, etc. 

%% Extract parameter
function [B, I,BI, BR,time] = FiniDiffMod_smallintestine(n, tRange, parameters)

tic; % start record time

% Create a function that handles faster evaluation of differential equations
fdAdt = @dAdt;

%Capture current time to measure computational speed
%----------------- initialize geometry and vectors --------------------
B0 = zeros(1,n);   % initialize BMP vector
I0 = zeros(1,n);   % initialize inhibitor vector
BI0 = zeros(1,n);   % initialize BMP-inhibitor complex vector
BR0 = zeros(1,n);   % initialize BMP-receptor vector
initial = [B0, I0, BI0, BR0];    % initialize condition vector

options = odeset('RelTol', 1e-9);    % set solving options relative error tolerance
% Set to force parameters
% solve ode using 15s
[time, D] = ode15s(fdAdt, tRange, initial, options, n,parameters);   

%------------------------- store data in vectors --------------------------
B = D(:, 1:n);      
I = D(:, n+1:2*n);
BI = D(:, 2*n+1:3*n);
BR = D(:, 3*n+1:4*n);
%--------------------------- make x data vector --------------------------- 
% start = -parameters.Ltot + parameters.Ltot/n;
% X = start:(parameters.Ltot/n):parameters.Ltot;
end

% This fucntion is the set of differential equations
function dY = dAdt(t, Y, n, parameters) 

%----------------------- obtain parameter values --------------------------
  %%% geometry
Ltot = parameters.Ltot;   % length of the intestine (1D assumption)       microns
LI = parameters.LI;   % length of stromal inhibitor expression region   microns

  %% diffusion rates
DB = parameters.DB;       % diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
DI = parameters.DI;       % diffusion rate of inhibitor        (microns^2*s^-1)*60s*m^-1
DBI = parameters.DBI;     % diffusion rate of [BMP,inhibitor]       (microns^2*s^-1)*60s*m^-1
  %% kinetic parameters
konI = parameters.konI;       % binding rates for BMP ligand and Inhibitor          nM^-1*m^-1
koffI = parameters.koffI;     % unbinding rates for BMP ligand and Inhibitor      m^-1
k1 = parameters.k1;       % binding rates for BMP ligand and Receptor         nM^-1*m^-1
k2 = parameters.k2;     % unbinding rates for BMP ligand and Receptor         m^-1
  %% Positive Feedback

  %% decay/recycle rates 
decB = parameters.decB;   % decay rate of Ligand (BMP)    nM*m^-1  
decI = parameters.decI;   % decay rate of inhibitor            nM*m^-1  
decBR = parameters.decBR;   % decay rate of BMP-receptor complex             nM*m^-1 
  %% production rates 
j1 = parameters.j1;       % production rate of BMP          nM*m^-1  
j2 = parameters.j2;       % production rate of Inhibitor      nM*m^-1
j3 = parameters.j3;       % production rate of BMP receptor      nM*m^-1
  %% Tolloid behavior
%------------- Convert boundaries to positions of nearest node ------------
nI = round(LI*n/Ltot);
%--- Increment size for finite difference method,distance between nodes --- 
dx = Ltot/(n-1);
x_lig = 0:dx:Ltot;% change x from space to node
x_tem=0:75:525;
%----------------- Initialize prepatterns for secretion -------------------
%% inhibitor pattern
etaI = zeros(1,n);
for i = 1:n
    etaI(i) = j2 * double( i>(n-nI) );  % define intestinal bottom region for inhibitor production
end
%% Bmp2 pattern
load('Bmp2_nor.mat')
p1 = polyfit(x_tem,Bmp2,3); % ployfitting of Bmp2 profile
y1 = polyval(p1,x_lig);
y1( y1 <= 0 ) = 0;
%plot(x_lig,y1)
etaB = j1 * y1;      % define intestinal top region for BMP production
%% Bmpr1a pattern
load('Bmpr1a_nor.mat')
p2 = polyfit(x_tem,Bmpr1a,3); % ployfitting of Bmp2 profile
y2 = polyval(p2,x_lig);
y2( y2 <= 0 ) = 0;
%plot(x_lig,y2)
% Bmpr1a uniform distribution
% y2=ones(1,61);
% 
Rtotal = j3 * y2;      % define intestinal top region for Bmpr1a production

%------------ convert Y vector into vectors for each component ------------
Bmp = Y(1:n);
I = Y(n+1:2*n);
BI= Y(2*n+1:3*n);
BR = Y(3*n+1:4*n);
%----------------------- Zero out difference vectors ----------------------
dBmp = zeros(1,n);
dI = zeros(1,n);
dBI = zeros(1,n);
dBR = zeros(1,n);
%--------------------Begin solving for ODEs at each time point--------
%--------------------Index one corresponds to ventral midline---------

dBmp(1) = DB/dx^2*(-2*Bmp(1)+2*Bmp(2)) ...
        - konI*Bmp(1)*I(1) + koffI*BI(1) ...
        - k1*Bmp(1)*(Rtotal(1)-BR(1)) + k2*BR(1) ...
        - decB*Bmp(1) + etaB(1);
    
dI(1) = DI/dx^2*(-2*I(1)+2*I(2)) ...
      - konI*Bmp(1)*I(1) + koffI*BI(1) ...
      - decI*I(1) + etaI(1);
      
dBI(1) = DBI/dx^2*(-2*BI(1)+2*BI(2)) ...
       +konI*Bmp(1)*I(1) - koffI*BI(1);
       
dBR(1) = k1*Bmp(1)*(Rtotal(1)-BR(1)) - k2*BR(1)- decBR*BR(1);
   
%----------------------Internal node points--------------------------    
for i=2:1:n-1
       
dBmp(i) = DB/dx^2*(Bmp(i-1)-2*Bmp(i)+Bmp(i+1)) ...
        - konI*Bmp(i)*I(i) + koffI*BI(i) ...
        - k1*Bmp(i)*(Rtotal(i)-BR(i)) + k2*BR(i) ...
        - decB*Bmp(i) + etaB(i);
    
dI(i) = DI/dx^2*(I(i-1)-2*I(i)+I(i+1)) ...
        - konI*Bmp(i)*I(i) + koffI*BI(i) ...
        - decI*I(i) + etaI(i);

dBI(i) = DBI/dx^2*(BI(i-1)-2*BI(i)+BI(i+1)) ...
         +konI*Bmp(i)*I(i) - koffI*BI(i);
       
dBR(i) = k1*Bmp(i)*(Rtotal(i)-BR(i)) - k2*BR(i)- decBR*BR(i);
       
end
%--------------nth node point corresponds to dorsal midline----------------  
dBmp(n) = DB/dx^2*(2*Bmp(n-1)-2*Bmp(n)) ...
        - konI*Bmp(n)*I(n) + koffI*BI(n) ...
        - k1*Bmp(n)*(Rtotal(n)-BR(n)) + k2*BR(n) ...
        - decB*Bmp(n) + etaB(n);
    
dI(n) = DI/dx^2*(2*I(n-1)-2*I(n)) ...
      - konI*Bmp(n)*I(n) + koffI*BI(n) ...
      - decI*I(n) + etaI(n);
      
dBI(n) = DBI/dx^2*(2*BI(n-1)-2*BI(n)) ...
       +konI*Bmp(n)*I(n) - koffI*BI(n);
   
dBR(n) = k1*Bmp(n)*(Rtotal(n)-BR(n)) - k2*BR(n)- decBR*BR(n);
   
%-----------------------Update solution vector---------------------------
dY = [dBmp dI dBI dBR]';

end