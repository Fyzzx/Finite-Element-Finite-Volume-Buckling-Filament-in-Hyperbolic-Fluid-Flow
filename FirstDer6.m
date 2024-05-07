function [Ld]=FirstDer6(LATT,delta_x) 

% function that computes a sixth order polynomial approximation to the
% first derivative (one dimensional).  LATT is the number of grid points
% and delta_x is the grid spacing.
% 

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
  
M = inv(M);

%----------------------------------------------------------------------------

% diagonal	

ia(1)=[1];	
ja(1)=[1];	
za(1)=[M(2,1)/dx];

ia(2)=[2];	
ja(2)=[2];	
za(2)=[-77/(60*dx)];

ia(3)=[3];
ja(3)=[3];
za(3)=[-7/(12*dx)];

ia(4:N-3)=[4:N-3];	
ja(4:N-3)=[4:N-3];	
za(4:N-3)=[0];   % supported boundary

ia(N-2)=[N-2];	
ja(N-2)=[N-2];	
za(N-2)=[7/(12*dx)];

ia(N-1)=[N-1];	
ja(N-1)=[N-1];	
za(N-1)=[77/(60*dx)];

ia(N)=[N];	
ja(N)=[N];	
za(N)=[-M(2,1)/dx ];

% 1st lower diagonal

ia(N+1)=[2];	
ja(N+1)=[1];	
za(N+1)=[-1/dx/6];

ia(N+2)=[3];	
ja(N+2)=[2];	
za(N+2)=[-2/dx/5];

ia(N+3:2*N-4)=[4:N-3];	
ja(N+3:2*N-4)=[3:N-4];	
za(N+3:2*N-4)=[-45/dx/60];

ia(2*N-3)=[N-2];	
ja(2*N-3)=[N-3];	
za(2*N-3)=[-4/dx/3];

ia(2*N-2)=[N-1];	
ja(2*N-2)=[N-2];	
za(2*N-2)=[-5/dx/2];

ia(2*N-1)=[N];	
ja(2*N-1)=[N-1];	
za(2*N-1)=[ -M(2,2)/dx ];

% 2nd lower diagonal

ia(2*N)=[3];	
ja(2*N)=[1];	
za(2*N)=[1/dx/30];

ia(2*N+1:3*N-6)=[4:N-3];	
ja(2*N+1:3*N-6)=[2:N-5];	
za(2*N+1:3*N-6)=[9/dx/60];

ia(3*N-5)=[N-2];	
ja(3*N-5)=[N-4];	
za(3*N-5)=[1/dx/2];

ia(3*N-4)=[N-1];	
ja(3*N-4)=[N-3];	
za(3*N-4)=[5/dx/3];

ia(3*N-3)=[N];	
ja(3*N-3)=[N-2];	
za(3*N-3)=[ -M(2,3)/dx ];

% 1st upper diagonal

ia(3*N-2)=[1];	
ja(3*N-2)=[2];	
za(3*N-2)=[ M(2,2)/dx ];

ia(3*N-1)=[2];	
ja(3*N-1)=[3];	
za(3*N-1)=[5/dx/2];

ia(3*N)=[3];	
ja(3*N)=[4];	
za(3*N)=[4/dx/3];

ia(3*N+1:4*N-6)=[4:N-3];	
ja(3*N+1:4*N-6)=[5:N-2];	
za(3*N+1:4*N-6)=[45/dx/60];

ia(4*N-5)=[N-2];	
ja(4*N-5)=[N-1];	
za(4*N-5)=[2/dx/5];

ia(4*N-4)=[N-1];	
ja(4*N-4)=[N];	
za(4*N-4)=[1/dx/6];

% 2nd upper diagonal

ia(4*N-3)=[1];	
ja(4*N-3)=[3];	
za(4*N-3)=[ M(2,3)/dx ];

ia(4*N-2)=[2];	
ja(4*N-2)=[4];	
za(4*N-2)=[-5/dx/3];

ia(4*N-1)=[3];	
ja(4*N-1)=[5];	
za(4*N-1)=[-1/dx/2];

ia(4*N:5*N-7)=[4:N-3];	
ja(4*N:5*N-7)=[6:N-1];	
za(4*N:5*N-7)=[-9/dx/60];

ia(5*N-6)=[N-2];	
ja(5*N-6)=[N];	
za(5*N-6)=[-1/dx/30];

% 3rd upper diagonal

ia(5*N-5)=[1];	
ja(5*N-5)=[4];	
za(5*N-5)=[ M(2,4)/dx ];

ia(5*N-4)=[2];	
ja(5*N-4)=[5];	
za(5*N-4)=[5/dx/6];

ia(5*N-3)=[3];	
ja(5*N-3)=[6];	
za(5*N-3)=[2/dx/15];

ia(5*N-2:6*N-9)=[4:N-3];	
ja(5*N-2:6*N-9)=[7:N];	
za(5*N-2:6*N-9)=[1/dx/60];

% 3rd lower diagonal

ia(6*N-8:7*N-15)=[4:N-3];	
ja(6*N-8:7*N-15)=[1:N-6];	
za(6*N-8:7*N-15)=[-1/dx/60];

ia(7*N-14)=[N-2];	
ja(7*N-14)=[N-5];	
za(7*N-14)=[-2/dx/15];

ia(7*N-13)=[N-1];	
ja(7*N-13)=[N-4];	
za(7*N-13)=[-5/dx/6];

ia(7*N-12)=[N];	
ja(7*N-12)=[N-3];	
za(7*N-12)=[ -M(2,4)/dx ];

% 4th upper diagonal

ia(7*N-11)=[1];
ja(7*N-11)=[5];
za(7*N-11)=[ M(2,5)/dx ];

ia(7*N-10)=[2];
ja(7*N-10)=[6];
za(7*N-10)=[-1/dx/4];

ia(7*N-9)=[3];
ja(7*N-9)=[7];
za(7*N-9)=[-1/dx/60];

% 4th lower diagonal

ia(7*N-8)=[N-2];
ja(7*N-8)=[N-6];
za(7*N-8)=[1/dx/60];

ia(7*N-7)=[N-1];
ja(7*N-7)=[N-5];
za(7*N-7)=[1/dx/4];

ia(7*N-6)=[N];
ja(7*N-6)=[N-4];
za(7*N-6)=[ -M(2,5)/dx ];

% 5th upper diagonal

ia(7*N-5)=[1];
ja(7*N-5)=[6];
za(7*N-5)=[ M(2,6)/dx ];

ia(7*N-4)=[2];
ja(7*N-4)=[7];
za(7*N-4)=[1/dx/30];

% 5th lower diagonal

ia(7*N-3)=[N-1];
ja(7*N-3)=[N-6];
za(7*N-3)=[-1/dx/30];

ia(7*N-2)=[N];
ja(7*N-2)=[N-5];
za(7*N-2)=[ -M(2,6)/dx ];

% 6th upper diagonal

ia(7*N-1)=[1];
ja(7*N-1)=[7];
za(7*N-1)=[ M(2,7)/dx ];

% 6th lower diagonal

ia(7*N)=[N];
ja(7*N)=[N-6];
za(7*N)=[ -M(2,7)/dx ];

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
