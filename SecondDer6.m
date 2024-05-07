function [Ld]=SecondDer6(LATT,delta_x) 
%
% 
% input --------------------------------------------------
% 
% kappa        -  curvature squared
% f            -  force
% delta_alpha  -  spacing 
%
% output --------------------------------------------------
%
% lam  - solution vector
%----------------------------------------------------------

dx=delta_x;
N=LATT;

% make the big matrix; make it sparse 

ia=zeros(7*N,1);  % i index assigned to zero
ja=zeros(7*N,1);  % j index assigned to zero
za=zeros(7*N,1);  % value index assigned to zero

% compute prefactors for derivative at first point

M = [ 1 0 0  0   0  0   0;
      1 1 1  1   1  1   1;
      1 2 4  8   16 32  64;
      1 3 9  27  81 243 3^6;
      1 4 16 4^3 4^4 4^5 4^6;
      1 5 25 5^3 5^4 5^5 5^6;
      1 6 36 6^3 6^4 6^5 6^6 ];
  
M = 2.*inv(M);

%--------------------------------------------------------------------------

% diagonal	

ia(1)=[1];	
ja(1)=[1];	
za(1)=[ M(3,1)/dx^2 ];

ia(2)=[2];	
ja(2)=[2];	
za(2)=[-147/(180*dx^2)];

ia(3)=[3];
ja(3)=[3];
za(3)=[-21/(9*dx^2)];

ia(4:N-3)=[4:N-3];	
ja(4:N-3)=[4:N-3];	
za(4:N-3)=[-49/dx^2/18];   % supported boundary

ia(N-2)=[N-2];	
ja(N-2)=[N-2];	
za(N-2)=[-21/(9*dx^2)];

ia(N-1)=[N-1];	
ja(N-1)=[N-1];	
za(N-1)=[-147/(180*dx^2)];

ia(N)=[N];	
ja(N)=[N];	
za(N)=[ M(3,1)/dx^2 ];

% 1st lower diagonal

ia(N+1)=[2];	
ja(N+1)=[1];	
za(N+1)=[137/dx^2/180];

ia(N+2)=[3];	
ja(N+2)=[2];	
za(N+2)=[114/dx^2/90];

ia(N+3:2*N-4)=[4:N-3];	
ja(N+3:2*N-4)=[3:N-4];	
za(N+3:2*N-4)=[6/dx^2/4];

ia(2*N-3)=[N-2];	
ja(2*N-3)=[N-3];	
za(2*N-3)=[10/dx^2/9];

ia(2*N-2)=[N-1];	
ja(2*N-2)=[N-2];	
za(2*N-2)=[-255/dx^2/180];

ia(2*N-1)=[N];	
ja(2*N-1)=[N-1];	
za(2*N-1)=[ M(3,2)/dx^2 ];

% 2nd lower diagonal

ia(2*N)=[3];	
ja(2*N)=[1];	
za(2*N)=[-13/dx^2/180];

ia(2*N+1:3*N-6)=[4:N-3];	
ja(2*N+1:3*N-6)=[2:N-5];	
za(2*N+1:3*N-6)=[-3/dx^2/20];

ia(3*N-5)=[N-2];	
ja(3*N-5)=[N-4];	
za(3*N-5)=[15/dx^2/180];

ia(3*N-4)=[N-1];	
ja(3*N-4)=[N-3];	
za(3*N-4)=[47/dx^2/18];

ia(3*N-3)=[N];	
ja(3*N-3)=[N-2];	
za(3*N-3)=[ M(3,3)/dx^2 ];

% 1st upper diagonal

ia(3*N-2)=[1];	
ja(3*N-2)=[2];	
za(3*N-2)=[ M(3,2)/dx^2 ];

ia(3*N-1)=[2];	
ja(3*N-1)=[3];	
za(3*N-1)=[-255/dx^2/180];

ia(3*N)=[3];	
ja(3*N)=[4];	
za(3*N)=[10/dx^2/9];

ia(3*N+1:4*N-6)=[4:N-3];	
ja(3*N+1:4*N-6)=[5:N-2];	
za(3*N+1:4*N-6)=[6/dx^2/4];

ia(4*N-5)=[N-2];	
ja(4*N-5)=[N-1];	
za(4*N-5)=[114/dx^2/90];

ia(4*N-4)=[N-1];	
ja(4*N-4)=[N];	
za(4*N-4)=[137/dx^2/180];

% 2nd upper diagonal

ia(4*N-3)=[1];	
ja(4*N-3)=[3];	
za(4*N-3)=[ M(3,3)/dx^2 ];

ia(4*N-2)=[2];	
ja(4*N-2)=[4];	
za(4*N-2)=[47/dx^2/18];

ia(4*N-1)=[3];	
ja(4*N-1)=[5];	
za(4*N-1)=[15/dx^2/180];

ia(4*N:5*N-7)=[4:N-3];	
ja(4*N:5*N-7)=[6:N-1];	
za(4*N:5*N-7)=[-3/dx^2/20];

ia(5*N-6)=[N-2];	
ja(5*N-6)=[N];	
za(5*N-6)=[-13/dx^2/180];

% 3rd upper diagonal

ia(5*N-5)=[1];	
ja(5*N-5)=[4];	
za(5*N-5)=[ M(3,4)/dx^2 ];

ia(5*N-4)=[2];	
ja(5*N-4)=[5];	
za(5*N-4)=[-285/dx^2/180];

ia(5*N-3)=[3];	
ja(5*N-3)=[6];	
za(5*N-3)=[-6/90/dx^2];

ia(5*N-2:6*N-9)=[4:N-3];	
ja(5*N-2:6*N-9)=[7:N];	
za(5*N-2:6*N-9)=[1/dx^2/90];

% 3rd lower diagonal

ia(6*N-8:7*N-15)=[4:N-3];	
ja(6*N-8:7*N-15)=[1:N-6];	
za(6*N-8:7*N-15)=[1/dx^2/90];

ia(7*N-14)=[N-2];	
ja(7*N-14)=[N-5];	
za(7*N-14)=[-6/dx^2/90];

ia(7*N-13)=[N-1];	
ja(7*N-13)=[N-4];	
za(7*N-13)=[-285/dx^2/180];

ia(7*N-12)=[N];	
ja(7*N-12)=[N-3];	
za(7*N-12)=[ M(3,4)/dx^2 ];

% 4th upper diagonal

ia(7*N-11)=[1];
ja(7*N-11)=[5];
za(7*N-11)=[ M(3,5)/dx^2 ];

ia(7*N-10)=[2];
ja(7*N-10)=[6];
za(7*N-10)=[93/dx^2/180];

ia(7*N-9)=[3];
ja(7*N-9)=[7];
za(7*N-9)=[1/dx^2/90];

% 4th lower diagonal

ia(7*N-8)=[N-2];
ja(7*N-8)=[N-6];
za(7*N-8)=[1/dx^2/90];

ia(7*N-7)=[N-1];
ja(7*N-7)=[N-5];
za(7*N-7)=[93/dx^2/180];

ia(7*N-6)=[N];
ja(7*N-6)=[N-4];
za(7*N-6)=[ M(3,5)/dx^2 ];

% 5th upper diagonal

ia(7*N-5)=[1];
ja(7*N-5)=[6];
za(7*N-5)=[ M(3,6)/dx^2 ];

ia(7*N-4)=[2];
ja(7*N-4)=[7];
za(7*N-4)=[-13/dx^2/180];

% 5th lower diagonal

ia(7*N-3)=[N-1];
ja(7*N-3)=[N-6];
za(7*N-3)=[-13/dx^2/180];

ia(7*N-2)=[N];
ja(7*N-2)=[N-5];
za(7*N-2)=[ M(3,6)/dx^2 ];

% 6th upper diagonal

ia(7*N-1)=[1];
ja(7*N-1)=[7];
za(7*N-1)=[ M(3,7)/dx^2 ];

% 6th lower diagonal

ia(7*N)=[N];
ja(7*N)=[N-6];
za(7*N)=[ M(3,7)/dx^2 ];

% make derivative symmetric

% ia(5*N-5)=[1];
% ja(5*N-5)=[N-1];
% za(5*N-5)=[0];
% 
% ia(5*N-4)=[1];
% ja(5*N-4)=[N];
% za(5*N-4)=[0];
% 
% ia(5*N-3)=[2];
% ja(5*N-3)=[N];
% za(5*N-3)=[0];
% 
% ia(5*N-2)=[N-1];
% ja(5*N-2)=[1];
% za(5*N-2)=[0];
% 
% ia(5*N-1)=[N];
% ja(5*N-1)=[1];
% za(5*N-1)=[0];
% 
% ia(5*N)=[N];
% ja(5*N)=[2];
% za(5*N)=[0];

% now making the sparce matrix
Ld=sparse(ia,ja,za);
