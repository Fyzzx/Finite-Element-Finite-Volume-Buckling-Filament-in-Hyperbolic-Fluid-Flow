function [r,e_1,e_2,r_s]=IntegrateMaterialFrameNew2(omega1,omega2,omega3,L)

% integrate the material frame using exponentiation of the integral of the
% strain matrix

N=length(omega3);
dx = L./(N-1);

L1 = FirstDer6(N,dx);
L1(1,:) = 0;
L1(1,1) = 1;

% initialize material frame and r

e_1=zeros(N,3);
e_2=zeros(N,3);
r_s = zeros(N,3);

r = zeros(N,3);

r(1,1) = 0;
r(1,2) = 0;
r(1,3) = 0;

r_s(1,1) = 0;
r_s(1,2) = 0;
r_s(1,3) = 1;

e_1(1,1)= 1;
e_1(1,2)= 0;
e_1(1,3)= 0;

e_2(1,1)= 0;
e_2(1,2)= 1;
e_2(1,3)= 0;

%% begin integration

for n = 2:N
        
    if n~=N
    
        Om1Int = (dx./12).*( -omega1(n+1) + 8.*omega1(n) + 5.*omega1(n-1) );
        Om2Int = (dx./12).*( -omega2(n+1) + 8.*omega2(n) + 5.*omega2(n-1) );
        Om3Int = (dx./12).*( -omega3(n+1) + 8.*omega3(n) + 5.*omega3(n-1) );
        
    else
        
        Om1Int = (dx./12).*( -omega1(n-2) + 8.*omega1(n-1) + 5.*omega1(n) );
        Om2Int = (dx./12).*( -omega2(n-2) + 8.*omega2(n-1) + 5.*omega2(n) );
        Om3Int = (dx./12).*( -omega3(n-2) + 8.*omega3(n-1) + 5.*omega3(n) );
        
    end
    
% define matrix

    Iden = [ 1 0 0; 0 1 0; 0 0 1 ];
    
    M = [  0        Om3Int  -Om2Int;
          -Om3Int   0        Om1Int;
          Om2Int   -Om1Int   0 ];
      
   Mexp = Iden + M + 0.5.*M*M + (1./6).*M*M*M + (1./24).*(M*M*M*M);
      
% integrate material frame

    e_1(n,1) = Mexp(1,:)*[ e_1(n-1,1); e_2(n-1,1); r_s(n-1,1) ];
    e_2(n,1) = Mexp(2,:)*[ e_1(n-1,1); e_2(n-1,1); r_s(n-1,1) ];
    r_s(n,1) = Mexp(3,:)*[ e_1(n-1,1); e_2(n-1,1); r_s(n-1,1) ];
    
    e_1(n,2) = Mexp(1,:)*[ e_1(n-1,2); e_2(n-1,2); r_s(n-1,2) ];
    e_2(n,2) = Mexp(2,:)*[ e_1(n-1,2); e_2(n-1,2); r_s(n-1,2) ];
    r_s(n,2) = Mexp(3,:)*[ e_1(n-1,2); e_2(n-1,2); r_s(n-1,2) ];
    
    e_1(n,3) = Mexp(1,:)*[ e_1(n-1,3); e_2(n-1,3); r_s(n-1,3) ];
    e_2(n,3) = Mexp(2,:)*[ e_1(n-1,3); e_2(n-1,3); r_s(n-1,3) ];
    r_s(n,3) = Mexp(3,:)*[ e_1(n-1,3); e_2(n-1,3); r_s(n-1,3) ];
    
    mag = sqrt( e_1(n,1).^2 + e_1(n,2).^2 + e_1(n,3).^2 );
    e_1(n,:) = e_1(n,:)./mag;
    
    mag = sqrt( e_2(n,1).^2 + e_2(n,2).^2 + e_2(n,3).^2 );
    e_2(n,:) = e_2(n,:)./mag;
    
    mag = sqrt( r_s(n,1).^2 + r_s(n,2).^2 + r_s(n,3).^2 );
    r_s(n,:) = r_s(n,:)./mag;
       
end

Source = r_s(:,1);
Source(1) = 0;
r(:,1) = L1\Source;

Source = r_s(:,2);
Source(1) = 0;
r(:,2) = L1\Source;

Source = r_s(:,3);
Source(1) = 0;
r(:,3) = L1\Source;

end