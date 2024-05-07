function [info,compDataEnd,myRunTime]=Hyperbolic_Flow_V3D(gamx,gamy,gamz,eta,kappa,tau0,A1,A2,Cc,rad,THF,GridNum,dt,sig0,RNGmod,len,ttime,samplet)

%% Define parameters
% Version 01
% hyperbolic flow profile. Rate*(-x + y). (compressile in x and stretching
% in y) where rate ~1 (units sec^-1) Flow velocity is uFL
% Flow start time is flowStart, flowSteady is how long it take for the flow to reach a steady state  
%% Define parameters
clear
%% input parameter (to comment out when running using the "run" file)
eta = .1;           % visocity in Pascal seconds (T media is 1.2cP) 0.001 = water
kappa = 1E-15;      % typical 5.0 for spiroplasma. try k=6.8 and 6 for leptospira
tau0 = 1E-15;       % rough value for spiro 6.8 when kappa 5.0
len = 5.0;          % total cell body length
rad = 0.07.*len/5.0;          % 0.08 for spiroplasma, 0.07 for leptospira

GridNum = round(len*100+1,0);   % number of nodes
ds = len./(GridNum-1);

A1 = 0.2;           % Bend modulus in e1 direction
A2 = 0.2;           % Bend modulus in e2 direction
Cc = 0.2*10;   % Twist modulus in e3 direction
sig0 = .8E5.*(A1./0.2); % non-linear spring force along line connecting points
sig0 = sig0.*ds;    % redefine sig0 based on ds

THF = 1;              % Thermal Factor. 1 = on, anything else (like 0) is off
RNGmod = 5;   % makes the random number generator seed random based on time 'shuffle' or 'default' which uses same set of random numbers each run

dt = 1.0*10^-6.*(A1./0.2);       % time step
ttime = .12;      % total time of simulation (used when no motor action.... (else use the final pulse of spiroplasma motor to set final time)
samplet = 2E-4;   % how often to collect the data and plot

gamx = -0.05E3;      % rate of hyperbolic fluid flow (sec^-1). negative is compressile
gamy = -gamx;         
gamz = 0;
% 
plotmystuff = 1;
format shortg
cbegin = clock;
factor = (len.^2./(4.*rad.*rad.*exp(1)));
goldSig = (2*eta*(gamy)*(len.^4))./(pi.^3.*((A1+A2)./2).*log(factor));
alph = rad/len;
factor2 = -log(alph.^2*exp(1));
muHat = 8*pi*eta*(len.^4).*abs(gamx)./(((A1+A2)./2).*factor2);

touch = 0;
%%
flowStart = 0.0005;     % time (seconds) when flow starts
flowSteady = 0.001;   % time (seconds) it takes to establish full steady state flow
%%

cRep = 0.01; % constant of repulsion f = (cRep)/delR^2
ThermStart = 1;
hRad = kappa./(kappa.^2 + tau0.^2);
pitch0 = 2*pi*tau0./(kappa.^2 + tau0.^2);

L = len;
Version = '01-Hyperbolic_3D';
s = linspace(0,L,GridNum);
ds = s(2) - s(1);
closenodes = round(10*rad./ds,0);

xmax = 0; ymax = 0; zmax = 0;

%%
       % time in seconds for first kink to appear
datestring = datestr(now,'yy_mm_dd_HH_MM');
rng(RNGmod) ; % 'default' with each new run of the code sets RNG to SAME set of random numbers. Each j,k loop iteration get different numbers. 'shuffle' sets the rng based on current time  
kT = .004114;      % room temp kT in pN*um
Lpers = (A1+A2)./(2.*kT);

k0 = kappa.*ones(GridNum,1);
a = ones(GridNum,1).*rad;       % radius of the cell cylinder

%% Original Drag and DragRat
zetaPerp = 4*pi*eta;   % Original and matches infinitely long rod in slender body theory. Value at eta = 0.001; ZPerp = 0.0126
zetaPar = 2*pi*eta;        % Makes zetaPar = 0.0063 when eta = 0.001
zetar = (4).*pi.*eta.*rad.^2;     % rotational drag coefficient for the cell cylinder. changed this to 4 (2) instead of 1/2 after 2/12/20 meet

dragrat = zetaPar/zetaPerp;

%% Thermal Force
if THF == 1
    TFe3 = sqrt(2*zetaPerp*dragrat*kT*ds/dt);
    TFe1 = sqrt(2*zetaPerp*kT*ds/dt);
    TFRot = sqrt(2.*zetar*kT*ds/dt);
   
else
    TFe3 = 0;
    TFe1 = 0;
    TFRot = 0;
end

%% create derivative matrices

L1 = FirstDer6(GridNum,ds); % finite element first derivative to 6th order accuracy
L2 = SecondDer6(GridNum,ds);  % finite element first derivative to 6th order accuracy

O1 =    sparse((3:GridNum-2),(2:GridNum-3),1,GridNum,GridNum) ...
      + sparse((3:GridNum-2),(3:GridNum-2),-2,GridNum,GridNum) ...
      + sparse((3:GridNum-2),(4:GridNum-1),1,GridNum,GridNum) ...
      + sparse(1,2,1,GridNum,GridNum) ...
      + sparse(2,(2:3),[-2 1],GridNum,GridNum) ...
      + sparse(GridNum-1,(GridNum-2:GridNum-1),[1 -2],GridNum,GridNum) ...
      + sparse(GridNum,GridNum-1,1,GridNum,GridNum);

O1 = O1./ds;

O2 =    sparse((3:GridNum-1),(2:GridNum-2),-0.5,GridNum,GridNum) ...
      + sparse((2:GridNum-2),(3:GridNum-1),0.5,GridNum,GridNum) ...
      + sparse(1,2,0.5,GridNum,GridNum) ...
      + sparse(GridNum,GridNum-1,-0.5,GridNum,GridNum);

O3 =   sparse((1:GridNum-1),(1:GridNum-1),1,GridNum,GridNum-1) ...
 + sparse((2:GridNum),(1:GridNum-1),-1,GridNum,GridNum-1);


%% Create some matrices

Zero = sparse(GridNum,GridNum);
CV =  sparse([1 GridNum],[1 GridNum],ds./2,GridNum,GridNum) ...
    + sparse((2:GridNum-1),(2:GridNum-1),ds,GridNum,GridNum);

CV2 = [ CV     Zero   Zero;
        Zero   CV     Zero;
        Zero   Zero   CV ];
% Identity matrix with 1/root(2) at ends
TFCM = sparse([1 GridNum],[1 GridNum],1./sqrt(2),GridNum,GridNum) ...
    + sparse((2:GridNum-1),(2:GridNum-1),1,GridNum,GridNum);
unvecT(1,:) = [1 0 0];
%% matrices for Phi equation

tau = tau0.*ones(GridNum,1);


%% Tors, create k01, k02, Om1, and Om2

k01 = k0.*ones(GridNum,1);
k02 = 0.0.*ones(GridNum,1);

Om1 = k01.*ones(GridNum,1);
Om2 = k02.*ones(GridNum,1);

%% Bend - Sets OM3 to zero for all time. Define shape by Om1 and Om2
% use curvatures to define the initial shape of the filament

Om3 = tau0.*ones(GridNum,1);
Om30 = tau0.*ones(GridNum,1);

[r,e1,e2,~]=IntegrateMaterialFrameNew2(Om1,Om2,Om3,L);

%% Rotate cell so that x-axis is aligned with helical axis
%% switching back to larger ds

x = r(:,1);
y = r(:,2);
z = r(:,3);

% rotate shape onto x-axis

qx = polyfit(s',x,1);
qy = polyfit(s',y,1);
qz = polyfit(s',z,1);

qx = qx(1);
qy = qy(1);
qz = qz(1);

c1 = qx./sqrt(qx.^2+qy.^2);
s1 = qy./sqrt(qx.^2+qy.^2);

c2 = sqrt(qx.^2+qy.^2)./sqrt(qx.^2+qy.^2 + qz.^2);
s2 = qz./sqrt(qx.^2+qy.^2 + qz.^2);

xnew = c2.*(c1.*x + s1.*y) + s2.*z;
ynew = -s1.*x + c1.*y;
znew = -s2.*(c1.*x+s1.*y) + c2.*z;

px = c2.*c1.*e1(:,1) + c2.*s1.*e1(:,2) +s2.*e1(:,3);
py = -s1.*e1(:,1) + c1.*e1(:,2);
pz = -s2.*c1.*e1(:,1) - s2.*s1.*e1(:,2) + c2.*e1(:,3);

x = xnew;
y = ynew;
z = znew;

x = x - mean(x);
y = y - mean(y);
z = z - mean(z);

e1(:,1) = px;
e1(:,2) = py;
e1(:,3) = pz;

mag = sqrt(e1(:,1).^2 + e1(:,2).^2 + e1(:,3).^2);

e1(:,1) = e1(:,1)./mag;
e1(:,2) = e1(:,2)./mag;
e1(:,3) = e1(:,3)./mag;

%% Setting up speed of kinks and timing
TotalTime = ttime;
Skip = (samplet)/dt; % will plot 200 times per 1.0 seconds
NSteps = ceil(TotalTime./dt./Skip) + 1; % total number of stored time steps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time(1) = 0.0;
Lee = len;
compData = [Time 0 0 0 0 0 0 Lee 0 0];
%%
j = 1;
while  j <= NSteps && touch == 0
j = j+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x(:,j) = x(:,j-1);
    y(:,j) = y(:,j-1);
    z(:,j) = z(:,j-1);
    Time(j) = Time(j-1);
    Iter = j;
    %%
%%%% Take Full Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:Skip
  %% Stepping forward in time, setting up tau to find kappa 1 and 2
        
Time(j) = Time(j) + dt;
tau(1:GridNum) = tau0;
t0Half = 0.5.*(tau(1:GridNum-1)+tau(2:GridNum));
Om30 = tau;

%% Calculating kappa 1 and 2 and setting up ends

k02 = 0.0.*ones(GridNum,1);

%% define derivatives
xs = L1*x(:,j);
ys = L1*y(:,j);
zs = L1*z(:,j);
x2s = L2*x(:,j);
y2s = L2*y(:,j);
z2s = L2*z(:,j);
         
%% define unit vectors

mag = sqrt( xs.^2 + ys.^2 + zs.^2 );
e3(:,1) = xs./mag;
e3(:,2) = ys./mag;
e3(:,3) = zs./mag;
     
e2(:,1) = e3(:,2).*e1(:,3) - e3(:,3).*e1(:,2);
e2(:,2) = e3(:,3).*e1(:,1) - e3(:,1).*e1(:,3);
e2(:,3) = e3(:,1).*e1(:,2) - e3(:,2).*e1(:,1);

mag2 = sqrt( e2(:,1).^2 + e2(:,2).^2 + e2(:,3).^2 );

e2(:,1) = e2(:,1)./mag2;
e2(:,2) = e2(:,2)./mag2;
e2(:,3) = e2(:,3)./mag2;

%% define curvatures    
Om1 = -( e2(:,1).*x2s + e2(:,2).*y2s + e2(:,3).*z2s ); 
Om2 =  ( e1(:,1).*x2s + e1(:,2).*y2s + e1(:,3).*z2s ); 

%     End point conditions
Om1(1) = k01(1);
Om1(GridNum) = k01(GridNum);
Om2(1) = k02(1);
Om2(GridNum) = k02(GridNum);

% for 2 bend moduli
bnF1 = A1*(k01-Om1) ; bnF2 = A2*(k02-Om2);
MOmD = -(x2s.*(bnF1.*e1(:,1) + bnF2.*e2(:,1)) + y2s.*(bnF1.*e1(:,2) + bnF2.*e2(:,2)) + z2s.*(bnF1.*e1(:,3) + bnF2.*e2(:,3))); 

         %% Aug 14,2020 adding in half step e1 transform. CW ThermalFil3

mTH = zeros(GridNum,1);
if THF == 1
    mTH = TFRot.*randn(GridNum,1);
end
MOmD = MOmD+mTH;

E2xMat = sparse((1:GridNum),(1:GridNum),e2(:,1),GridNum,GridNum);
E2yMat = sparse((1:GridNum),(1:GridNum),e2(:,2),GridNum,GridNum);
E2zMat = sparse((1:GridNum),(1:GridNum),e2(:,3),GridNum,GridNum);        
        
TwX = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,1)/ds,GridNum-1,GridNum) ...
      + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,1)/ds,GridNum-1,GridNum);
TwY = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,2)/ds,GridNum-1,GridNum) ...
      + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,2)/ds,GridNum-1,GridNum);
TwZ = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,3)/ds,GridNum-1,GridNum) ...
      + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,3)/ds,GridNum-1,GridNum);

Tw = [ TwX, TwY, TwZ ];
e1Vec = [ e1(:,1); e1(:,2); e1(:,3) ];

Om3Half = Tw*e1Vec;
Om3(1) = tau(1); % Om3 initial and final equal Om30 (for bend code = 0)
Om3(2:GridNum-1) = 0.5.*(Om3Half(1:GridNum-2) + Om3Half(2:GridNum-1));
Om3(GridNum) = tau(GridNum);

TwMat = [ E2xMat*O3;
          E2yMat*O3;
          E2zMat*O3 ];
      
MatE1 = CV2 - (Cc.*dt./zetar).*(TwMat*Tw);
SourceE1 = CV2*e1Vec - (Cc.*dt./zetar).*(TwMat*t0Half);
SourceE1 = SourceE1 - (dt./zetar).*[ CV*(MOmD.*e2(:,1)); (CV*MOmD.*e2(:,2)); (CV*MOmD.*e2(:,3)) ];

Ans = MatE1\SourceE1;

e1m(:,1) = Ans(1:GridNum);
e1m(:,2) = Ans(GridNum+1:2.*GridNum);
e1m(:,3) = Ans(2.*GridNum+1:3.*GridNum);
       
e1mmag = sqrt(e1m(:,1).^2 + e1m(:,2).^2 + e1m(:,3).^2);
e1m(:,1) = e1m(:,1)./e1mmag;  % Not part of the ThermalFil3 code
e1m(:,2) = e1m(:,2)./e1mmag;
e1m(:,3) = e1m(:,3)./e1mmag;  
    
%% MOmX from e3 x dM/ds
% [e3]x[dM/ds] = d/ds[e3 x M] - [d/ds(e3)]x M = -e3*F + F(vector)
% for 2 bend moduli
MOmX = A1.*(Om1 - k01).*Om1 + A2.*(Om2 - k02).*Om2;       
MOmX = 0.5.*(MOmX(1:GridNum-1) + MOmX(2:GridNum)); 

%% solve for tangential force
      
DelS(1:GridNum-1) = sqrt(   ( x(2:GridNum,j) - x(1:GridNum-1,j) ).^2 ... 
                          + ( y(2:GridNum,j) - y(1:GridNum-1,j) ).^2 ...
                          + ( z(2:GridNum,j) - z(1:GridNum-1,j) ).^2 );

q = sig0.*(1 - ds./DelS');
q(q>0) = q(q>0)./3;
% Linear
% q = sig0.*((DelS')./ds-1);

TanForce = dt.*(   sparse((1:GridNum-1),(1:GridNum-1),(q-MOmX)./ds,GridNum,GridNum) ...
                 - sparse((1:GridNum-1),(2:GridNum),(q-MOmX)./ds,GridNum,GridNum) ...
                 - sparse((2:GridNum),(1:GridNum-1),(q-MOmX)./ds,GridNum,GridNum) ...
                 + sparse((2:GridNum),(2:GridNum),(q-MOmX)./ds,GridNum,GridNum) );
%% Generate random force per ds. 

Ftx = zeros(GridNum,1);
Fty = zeros(GridNum,1);
Ftz = zeros(GridNum,1);

if THF == 1  

    randvar = normrnd(0,1, [GridNum*3,1]);
    Fte1 = randvar(1:GridNum,1).*TFe1; % thermal force causing cylinder to translate in e1
    Fte2 = randvar(GridNum+1:2*GridNum,1).*TFe1;
    Fte3 = randvar(2*GridNum+1:end,1).*TFe3; % thermal force causing cylinder to translate in e3
    Fte1 = TFCM*Fte1; % modifying edges by factor of sqrt(2)
    Fte2 = TFCM*Fte2;
    Fte3 = TFCM*Fte3;

    Ftx = Fte1.*e1(:,1) + Fte2.*e2(:,1) + Fte3.*e3(:,1);
    Fty = Fte1.*e1(:,2) + Fte2.*e2(:,2) + Fte3.*e3(:,2);
    Ftz = Fte1.*e1(:,3) + Fte2.*e2(:,3) + Fte3.*e3(:,3); 
    
end
   
%% Self repulsion force
fRep = zeros(GridNum,3);
rnow = [x(:,j) y(:,j) z(:,j)];
for ii = 1:GridNum-1
    delR = 1E5.*ones(GridNum-ii,1); % set all distances to 10^5
    maxplace = min(closenodes,GridNum-ii); % finding the node 
    for jj = 1:GridNum-ii
        diff = rnow(ii,:) - rnow(jj+ii,:);
        delR(jj,1) = sqrt(diff(1).^2 + diff(2).^2 + diff(3).^2); 
        
    end
    delR(1:maxplace,1) = 1E5; % sets near nodes equal to distance 100
    delR(delR>(2.1*rad)) = 1E5; % sets anything that lies outside of 2.1*rad away from another node as distance = 10^5 
    if min(delR)<10 % checks to see if any nodes outside of exclusion zone are close
        if touch == 0
            touch = touch + 1; % sets the boolean for data extraction
        end
        repWhere = find(delR<2.2*rad); % find close nodes
        for kk = 1:length(repWhere) % cycle through close locations and find the actual distance
            magU = sqrt((rnow(ii+repWhere(kk),1)-rnow(ii,1)).^2+(rnow(ii+repWhere(kk),2)-rnow(ii,2)).^2 +(rnow(ii+repWhere(kk),3)-rnow(ii,3)).^2);
            unR =  (rnow(ii+repWhere(kk),:)-rnow(ii,:))./magU;
            fRep(ii,:) = fRep(ii,:)- unR.*(cRep/((delR(repWhere(kk))).^2)); % determines force of reuplsion
            fRep(ii+repWhere(kk),:) = fRep(ii+repWhere(kk),:) + unR.*(cRep/((delR(repWhere(kk))).^2)); % equal and opposite
        end
    end
        
end

%% integrate position full time step

Mxx = zetaPerp.*(CV*sparse((1:GridNum),(1:GridNum),1-(1-dragrat).*e3(:,1).^2,GridNum,GridNum));
Mxy = zetaPerp.*(CV*sparse((1:GridNum),(1:GridNum),(dragrat-1).*e3(:,1).*e3(:,2),GridNum,GridNum)); 
Mxz = zetaPerp.*(CV*sparse((1:GridNum),(1:GridNum),(dragrat-1).*e3(:,1).*e3(:,3),GridNum,GridNum));

Myy = zetaPerp.*(CV*sparse((1:GridNum),(1:GridNum),1-(1-dragrat).*e3(:,2).^2,GridNum,GridNum));
Myz = zetaPerp.*(CV*sparse((1:GridNum),(1:GridNum),(dragrat-1).*e3(:,2).*e3(:,3),GridNum,GridNum));

Mzz = zetaPerp.*(CV*sparse((1:GridNum),(1:GridNum),1-(1-dragrat).*e3(:,3).^2,GridNum,GridNum));


%% For 2 moduli
if A2 >= A1
    
BigMat = [ Mxx+(A2.*dt).*(O1*L2)  Mxy                       Mxz;
           Mxy                     Myy+(A2.*dt).*(O1*L2)    Myz;
           Mxz                     Myz                       Mzz+(A2.*dt).*(O1*L2) ];

    hrssX = (O1*(Om1.*e2(:,1))); hrssY = (O1*(Om1.*e2(:,2))); hrssZ = (O1*(Om1.*e2(:,3)));
else
BigMat = [ Mxx+(A1.*dt).*(O1*L2)  Mxy                       Mxz;
           Mxy                     Myy+(A1.*dt).*(O1*L2)    Myz;
           Mxz                     Myz                       Mzz+(A1.*dt).*(O1*L2) ];
    hrssX = (O1*(Om2.*e1(:,1))); hrssY = (O1*(Om2.*e1(:,2))); hrssZ = (O1*(Om2.*e1(:,3)));
end

%% Fluid Force

Fxx = zetaPerp.*ds.*gamx.*x(:,j).*(1-(1-dragrat).*e3(:,1).^2);
Fxy = zetaPerp.*ds.*gamx.*x(:,j).*(dragrat-1).*e3(:,1).*e3(:,2); 
Fxz = zetaPerp.*ds.*gamx.*x(:,j).*(dragrat-1).*e3(:,1).*e3(:,3);

Fyx = zetaPerp.*ds.*gamy.*y(:,j).*(dragrat-1).*e3(:,1).*e3(:,2); 
Fyy = zetaPerp.*ds.*gamy.*y(:,j).*(1-(1-dragrat).*e3(:,2).^2);
Fyz = zetaPerp.*ds.*gamy.*y(:,j).*(dragrat-1).*e3(:,2).*e3(:,3);

Fzx = zetaPerp.*ds.*gamz.*z(:,j).*(dragrat-1).*e3(:,1).*e3(:,3); 
Fzy = zetaPerp.*ds.*gamz.*z(:,j).*(dragrat-1).*e3(:,2).*e3(:,3);
Fzz = zetaPerp.*ds.*gamz.*z(:,j).*(1-(dragrat-1)).*e3(:,3).*e3(:,3);

%% Fluid force in x, y, z directions
fluX = Fxx + Fyx + Fzx; % fluid force in the x hat
fluY = Fxy + Fyy + Fzy;
fluZ = Fxz + Fyz + Fzz;

%% changed RHS equations

hx = Cc.*(Om3 - Om30).*(ys.*z2s - zs.*y2s); 
hy = Cc.*(Om3 - Om30).*(zs.*x2s - xs.*z2s); 
hz = Cc.*(Om3 - Om30).*(xs.*y2s - ys.*x2s);

XRHS =   Mxx*x(:,j) + Mxy*y(:,j) + Mxz*z(:,j) - TanForce*(x(:,j)) ...
    + dt.*( -O1*(A1.*k01.*e2(:,1)) + O1*(A2.*k02.*e1(:,1)) + O2*(hx) + Ftx + fluX) ...
    + dt.*(A1 - A2).*hrssX;

YRHS =   Mxy*x(:,j) + Myy*y(:,j) + Myz*z(:,j) - TanForce*(y(:,j)) ...
    + dt.*( -O1*(A1.*k01.*e2(:,2)) + O1*(A2.*k02.*e1(:,2)) + O2*(hy) + Fty + fluY) ...
    + dt.*(A1 - A2).*hrssY;

ZRHS =   Mxz*x(:,j) + Myz*y(:,j) + Mzz*z(:,j) - TanForce*(z(:,j)) ...
    + dt.*( -O1*(A1.*k01.*e2(:,3)) + O1*(A2.*k02.*e1(:,3)) + O2*(hz) + Ftz + fluZ) ...
    + dt.*(A1 - A2).*hrssZ; 

%% Collecting Data      
    
Source =  [ XRHS;
            YRHS;
            ZRHS ];

Answer = BigMat\Source;

x(:,j) = Answer(1:GridNum);
y(:,j) = Answer(GridNum+1:2*GridNum);
z(:,j) = Answer(2*GridNum+1:3*GridNum);


%% new e1b & e1c update
    xs2 = L1*x(:,j);
    ys2 = L1*y(:,j);
    zs2 = L1*z(:,j);
    mag = sqrt( xs2.^2 + ys2.^2 + zs2.^2 );
    
    % find new e3 with latest xyz
    e32(:,1) = xs2./mag;
    e32(:,2) = ys2./mag;
    e32(:,3) = zs2./mag;

    % calculating e1&e2 dotted into the new e3
    e1mdote3 = e1m(:,1).*e32(:,1) + e1m(:,2).*e32(:,2) + e1m(:,3).*e32(:,3);

    costheta = e3(:,1).*e32(:,1) + e3(:,2).*e32(:,2) + e3(:,3).*e32(:,3);
 
    e1b(:,1) =  e1m(:,1) - e1mdote3.*(e3(:,1)+e32(:,1))./(1+costheta);
            
    e1b(:,2) =  e1m(:,2) - e1mdote3.*(e3(:,2)+e32(:,2))./(1+costheta);
            
    e1b(:,3) =  e1m(:,3) - e1mdote3.*(e3(:,3)+e32(:,3))./(1+costheta);
    mag = sqrt(e1b(:,1).^2+e1b(:,2).^2+e1b(:,3).^2);
        
    e1b(:,1) = e1b(:,1)./mag;
    e1b(:,2) = e1b(:,2)./mag;
    e1b(:,3) = e1b(:,3)./mag;
    e1 = e1b;

    end  % end of k loop
    % Collecting data for output excel file. 
    
    unvec1 = polyfit(s',x(:,j),1);
    unvec2 = polyfit(s',y(:,j),1);
    unvec3 = polyfit(s',z(:,j),1);
    unvecT(j,:) = [unvec1(1) unvec2(1) unvec3(1)];
    unvecT(j,:) = unvecT(j,:) ./ norm(unvecT(j,:));
    t2tL(1,1) = len;
    t2tL(j,1) = sqrt((x(GridNum,j)-x(1,j)).^2 + (y(GridNum,j)-y(1,j)).^2 + (z(GridNum,j)-z(1,j)).^2);
 %%   Plotting
if plotmystuff == 1
 if j == 2 || j == 1
    Hfig = figure(1); % nothing on first iter... did pop up on 2nd
    clf(1);
    end
    
    if j>=3
    clf(1);
    % Hfig.WindowState = 'minimized'; % %minimizes plot window so it doesn't interupt typing 
    end
    
    % Om3 = zeros(size(x));
    r(:,1) = x(:,j);
    r(:,2) = y(:,j);
    r(:,3) = z(:,j);
    
    [X,Y,Z] = sphere(60);
    % for linear flow
    % for x-y hyperbolic free flow
        if A1==A2
            tubeplot2(r,Om3,a,50,[0.3 0.7 0.4],[0.3 0.7 0.4],1);
        else
            ribbonplot(r,e1,e2,0.75.*rad,0.2.*rad,[0.8 0.3 0.4],1);
        end
    
    axis equal;
    axis([-L*0.7 L*0.7 -L*0.7 L*0.7 -L*0.7 L*0.7]);
    view([0 90])
    
    hold on
    plot3(0,0,0,'x');
    spacings = 11;
    xP = linspace(-L*.7,L*.7,spacings);
    yP = linspace(-L*.5,L*.5,spacings);
    zP = linspace(-L*.5,L*.5,spacings);
    
    [mgx,mgy,mgz] = meshgrid(xP,yP,zP);
    flPx = gamx.*mgx;
    flPy = gamy.*mgy;
    flPz = gamz.*mgz;
    
    flPx2 = zeros(spacings,spacings,spacings);
    flPy2 = zeros(spacings,spacings,spacings);
    flPz2 = zeros(spacings,spacings,spacings);
    
    flPx2(:,:,ceil(spacings/2)) = flPx(:,:,ceil(spacings/2));
    flPy2(:,:,ceil(spacings/2)) = flPy(:,:,ceil(spacings/2));
    flPz2(:,:,ceil(spacings/2)) = flPz(:,:,ceil(spacings/2));
    
    quiver3(mgx,mgy,mgz,flPx2,flPy2,flPz2)
        if A1==A2
            h1 = patch(surf2patch(a(GridNum).*X+r(GridNum,1),a(GridNum).*Y+r(GridNum,2),a(GridNum).*Z+r(GridNum,3),Z)); % patch last point with sphere
            set(h1,'FaceColor',[0.3 0.7 0.4],'EdgeColor','none'); 
            h2 = patch(surf2patch(a(1).*X+r(1,1),a(1).*Y+r(1,2),a(1).*Z+r(1,3),Z)); % patch 1st point with sphere
            set(h2,'FaceColor',[0.3 0.7 0.4],'EdgeColor','none');box off;axis off; 
        end
    
    lightangle(0,+40); % at -40 for spiroplasma helix stuff but also the view is set at [2 15] normally. y axis is into the board
    hold off
    set(gcf,'Color',[1 1 1])
    
    % set(Hfig, 'Position',  [3500, -300, 800, 500]);     % for 2 screens 
    set(Hfig, 'Position',  [50, 50, 1400, 1000]);         % for 1 screen  (x0,y0,xFin,yFin)
    %% Save figure and display data

mov(:,j-1) = getframe(Hfig);
end
%% compression helicoidal data
ymax = max(max(abs(y(:,j)),ymax));
zmax = max(max(abs(z(:,j)),zmax));
LeeOld = Lee;
Lee = x(GridNum,j) - x(1,j); % x length end to end
Reff = sqrt(mean((y(round(GridNum/3,0):GridNum-round(GridNum/3,0),j)-mean(y(:,j))).^2 + (z(round(GridNum/3,0):GridNum-round(GridNum/3,0),j)-mean(y(:,j))).^2)); % Chakrabarti 2020
LeeDot = abs(Lee-len)./Time(j); % end to end length rate of change
epsDot = abs(gamx); % strain rate (assumed from Chakrabarti 2020 paper)
bFact = ((A1+A2)./(2*eta*epsDot)).^(0.25); % chakrabarti 2020 ... modded to handle two different bend moduli
compData(j,:) = [Time(j) LeeDot ymax zmax Reff bFact epsDot Lee Reff/len muHat  ];

display = ["Iter","Max Del/ds","Length","end2end","maxY","maxZ";...
           strcat(string(Iter),'/',string(NSteps)),max(abs(1-DelS./ds)),sum(DelS./L),Lee,ymax,zmax]
% end of 'j' loop
if mod(Iter,1) == 0 && plotmystuff ==1 
    
    figure(2);
    clf(2)
    r(:,1) = x(:,j);
    r(:,2) = y(:,j);
    r(:,3) = z(:,j);
    
    [X,Y,Z] = sphere(60);
    % for linear flow
    % tubeplot2(r,Om3,a,20,[0.3 0.7 0.4],[0.3 0.7 0.4],1);axis equal;axis([-L*0.5 L*1.2 -hRad*4 hRad*4 -hRad*4 hRad*4]);view([2 15])
    % for x-y hyperbolic free flow
    if A1==A2
        tubeplot2(r,Om3,a,50,[0.3 0.7 0.4],[0.3 0.7 0.4],1);
    else
        ribbonplot(r,e1,e2,0.75.*rad,0.2.*rad,[0.8 0.3 0.4],1);
    end
    ySc = 0.20;
    zSc = 0.20;   
    axis equal;
    axis([-L*0.7 L*0.7 -L*(ySc+.1) L*(ySc+.1) -L*(zSc+.1) L*(zSc+.1)]);
    view([-47 23])
    
    hold on
    plot3(0,0,0,'x');
    
    flPx2 = zeros(spacings,spacings,spacings);
    flPy2 = zeros(spacings,spacings,spacings);
    flPz2 = zeros(spacings,spacings,spacings);
    
    flPx2(:,:,ceil(spacings/2)-2) = flPx(:,:,ceil(spacings/2));
    flPy2(:,:,ceil(spacings/2)-2) = flPy(:,:,ceil(spacings/2));
    flPz2(:,:,ceil(spacings/2)-2) = flPz(:,:,ceil(spacings/2));

    quiver3(mgx,mgy,mgz,flPx2,flPy2,flPz2,'LineWidth',2.0)
    
    if A1==A2
        h1 = patch(surf2patch(a(GridNum).*X+r(GridNum,1),a(GridNum).*Y+r(GridNum,2),a(GridNum).*Z+r(GridNum,3),Z)); % patch last point with sphere
        set(h1,'FaceColor',[0.3 0.7 0.4],'EdgeColor','none'); 
        h2 = patch(surf2patch(a(1).*X+r(1,1),a(1).*Y+r(1,2),a(1).*Z+r(1,3),Z)); % patch 1st point with sphere
        set(h2,'FaceColor',[0.3 0.7 0.4],'EdgeColor','none');box off;axis off; 
    end

    box1p1 = [-L.*0.7 -L.*ySc -L.*zSc];
    box1p2 = [-L.*0.7 +L.*ySc -L.*zSc];
    box1p3 = [-L.*0.7 +L.*ySc +L.*zSc];
    box1p4 = [-L.*0.7 -L.*ySc +L.*zSc];
    box1plot = [box1p1;box1p2;box1p3;box1p4];
    box2p1 = [L.*0.7 -L.*ySc -L.*zSc];
    box2p2 = [L.*0.7 +L.*ySc -L.*zSc];
    box2p3 = [L.*0.7 +L.*ySc +L.*zSc];
    box2p4 = [L.*0.7 -L.*ySc +L.*zSc];
    box2plot = [box2p1;box2p2;box2p3;box2p4];
    
    for oo = 1:4
        plist1 = mod(oo-1,4)+1;
        plist2 = mod(oo,4)+1;
        linenow = [box1plot(plist1,:);box1plot(plist2,:)];
        plot3(linenow(:,1), linenow(:,2), linenow(:,3), 'LineWidth',2.0,'color',[0,0,0]+0.8)
        linenow = [box2plot(plist1,:);box2plot(plist2,:)];
        plot3(linenow(:,1), linenow(:,2), linenow(:,3), 'LineWidth',2.0,'color',[0,0,0]+0.8)
        linenow = [box1plot(oo,:);box2plot(oo,:)];
        plot3(linenow(:,1), linenow(:,2), linenow(:,3), 'LineWidth',2.0,'color',[0,0,0]+0.8)
    
    end  
    
    plot3(x(:,j),y(:,j),zeros(GridNum,1)-L.*zSc, 'LineWidth',3.0,'color',[0,0,0]+0.6);
    plot3(zeros(GridNum,1)+L.*0.7,y(:,j),z(:,j), 'LineWidth',3.0,'color',[0,0,0]+0.6);
    plot3(x(:,j),zeros(GridNum,1)+L.*ySc,z(:,j), 'LineWidth',3.0,'color',[0,0,0]+0.6);
    
    
    lightangle(0,+70); % at -40 for spiroplasma helix stuff but also the view is set at [2 15] normally. y axis is into the board
    hold off
    set(gcf,'Color',[1 1 1])
    
    % set(Hfig, 'Position',  [3500, -300, 800, 500]);     % for 2 screens 
    set(Hfig, 'Position',  [50, 50, 1400, 1000]);         % for 1 screen  (x0,y0,xFin,yFin)
    
    mov2(:,j-1) = getframe(figure(2));

    clf(2)
    figure(1)
    
end

if touch == 1
    "Aborted because filament touched itself"
end
end


%% writing figures to AVI
vidtitle2 = strcat('HyperbolicBox_V',Version,'_',datestring,'_sig',num2str(sig0,'%6.3g'),'_eta', num2str(eta*1000, '%6.3g'),'_Ac1_', num2str(A1, '%6.3g'),'_Ac2_', num2str(A2, '%6.3g'),'_Cc_', num2str(Cc, '%6.3g'),'_gamma', num2str(gamx, '%6.3g'), '_dt', num2str(dt, '%6.3g'),'_GN',num2str(GridNum,'%6.3g'), '_L',num2str(L,'%6.3g'), '_kap_',num2str(kappa,'%6.3g'),'_tau0_',num2str(tau0,'%6.3g'),'_rad',num2str(rad,'%6.3g'), '_RNG',num2str(RNGmod,'%6.3g'));

%% Create tables of data
% Tx = table(x);
% Ty = table(y);
% Tz = table(z);

%%
finL = sum(DelS);
goldY = 1-min(t2tL)./len;
info = [string(eta), num2str(gamx),num2str(gamy), num2str(gamz), num2str(len),num2str(rad),num2str(A1), num2str(A2),num2str(Cc),Version, ...
    num2str(sig0), num2str(GridNum), num2str(dt),num2str(goldSig), num2str(goldY), num2str(finL)];

%% Write tables to excel file

compDataEnd = compData(end,:);

%% write video
% writerObj = VideoWriter(strcat(vidtitle));
% writerObj.FrameRate = 20; %FPS
% %Open Video Writer
% open(writerObj);
% for i=1:length(mov)
%     %convert image to frame of movie
%     frame = mov(i);
%     writeVideo(writerObj,frame);
% end
% %close the video writer
% close(writerObj);
%% write video 2
writerObj = VideoWriter(strcat(vidtitle2));
writerObj.FrameRate = 20; %FPS
%Open Video Writer
open(writerObj);
for i=1:length(mov2)
    %convert image to frame of movie
    frame = mov2(i);
    writeVideo(writerObj,frame);
end
%close the video writer
close(writerObj);

cend = clock;
myRunTime = cend(6)-cbegin(6);
keyboard
end