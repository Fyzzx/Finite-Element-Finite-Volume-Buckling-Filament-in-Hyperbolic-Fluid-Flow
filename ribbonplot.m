function ribbonplot(r,e1,e2,Width,Thick,RibbonColor,transparency)

% plots a 3D ribbon along the curve defined by
% the vector r, assuming that the wide axis of the ribbon
% lies along the orthonormal direction e1.
% e2 is the direction perpendicular to the tangent and e1.  
% The Ribbon width and thickness are 
% Width and Thick, respectively.

N=length(r(:,1));

% define signs for e1 and e2 vertices

c1 = [ -1 1 1 -1];
c2 = [ 1 1 -1 -1];

% define patches

for m=1:N
    for i=1:4
        for j=1:3
        
        vertices(4*(m-1)+i,j) = r(m,j) ...
                  + (c1(i).*Width./2).*e1(m,j) ...
                  + (c2(i).*Thick./2).*e2(m,j);
        end
    end
end

% define connectivity matrix

vec1 = (1:4);
vec2 = circshift(vec1,[0 -1]);

for m=1:4
        
        MT(m,1) = vec1(m);
        MT(m,2) = vec2(m);
        MT(m,4) = vec1(m) + 4;
        MT(m,3) = MT(m,2) + 4;
        
end

for i=2:N-1
    
    for m=1:4
        
        M(m,1) = (i-1)*4+vec1(m);
        M(m,2) = (i-1)*4+vec2(m);
        M(m,4) = M(m,1) + 4;
        M(m,3) = M(m,2) + 4;
        
    end
    
    MT = [ MT ; M ];
end

fsv = zeros(4*N,3);

for k=1:N
    fsv(4*(k-1)+1,:) = RibbonColor;
    fsv(4*(k-1)+2,:) = [ 0 0 0 ];
    fsv(4*(k-1)+3,:) = RibbonColor;
    fsv(4*(k-1)+4,:) = [ 0 0 0 ];  
    
end

view(0,35)
%  lightangle(0,45)
%  lightangle(40,45)
% lightangle(165,30)
set(gcf,'Renderer','zbuffer')

set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',0.1,'DiffuseStrength',.1,...
    'SpecularStrength',0.8,'SpecularExponent',25,...
    'SpecularColorReflectance',0.1,...
    'BackFaceLighting','off')

patch('Vertices',vertices,'Faces',MT,'FaceVertexCData',...
    fsv,'FaceAlpha',transparency,'FaceColor','flat','EdgeColor','none'),axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subroutines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ld]=first_der(LATT,dx) 

% creates a matrix that calculates the second order
% first derivative with unequal spacings dx

N=LATT;

% make the big matrix; make it sparse 

ia=zeros(3*N,1);  % i index assigned to zero
ja=zeros(3*N,1);  % j index assigned to zero
za=zeros(3*N,1);  % value index assigned to zero

%----------------------------------------------------------------------------

% diagonal	

ia(1)=[1];	
ja(1)=[1];	
za(1)=[-(2*dx(1)+dx(2))./dx(1)./(dx(1)+dx(2))];

ia(2:N-1)=[2:N-1];	
ja(2:N-1)=[2:N-1];	
za(2:N-1)=[(dx(2:N-1)-dx(1:N-2))./dx(2:N-1)./dx(1:N-2)];   % supported boundary

ia(N)=[N];	
ja(N)=[N];	
za(N)=[(2*dx(N-1)+dx(N-2))./dx(N-1)./(dx(N-1)+dx(N-2))];

% 1st lower diagonal

ia(N+1:2*N-2)=[2:N-1];	
ja(N+1:2*N-2)=[1:N-2];	
za(N+1:2*N-2)=[-dx(2:N-1)./dx(1:N-2)./(dx(1:N-2)+dx(2:N-1))];

ia(2*N-1)=[N];	
ja(2*N-1)=[N-1];	
za(2*N-1)=[-(dx(N-2)+dx(N-1))./dx(N-1)./dx(N-2)];

% 2nd lower diagonal

ia(2*N)=[N];	
ja(2*N)=[N-2];	
za(2*N)=[dx(N-1)./(dx(N-1)+dx(N-2))./dx(N-2)];

% 1st upper diagonal

ia(2*N+1)=[1];	
ja(2*N+1)=[2];	
za(2*N+1)=[(dx(1)+dx(2))./dx(1)./dx(2)];

ia(2*N+2:3*N-1)=[2:N-1];	
ja(2*N+2:3*N-1)=[3:N];	
za(2*N+2:3*N-1)=[dx(1:N-2)./dx(2:N-1)./(dx(1:N-2)+dx(2:N-1))];

% 2nd upper diagonal

ia(3*N)=[1];	
ja(3*N)=[3];	
za(3*N)=[-dx(1)./(dx(1)+dx(2))./dx(2)];

% now making the sparce matrix
Ld=sparse(ia,ja,za);


function [Ld]=second_der(LATT,dx)

% creates a matrix that calculates the second order second
% derivative with unequal spacings dx

N=LATT;

% make the big matrix; make it sparse 

ia=zeros(3*N,1);  % i index assigned to zero
ja=zeros(3*N,1);  % j index assigned to zero
za=zeros(3*N,1);  % value index assigned to zero

%----------------------------------------------------------------------------


% diagonal	

ia(1)=[1];	
ja(1)=[1];	
za(1)=[2./dx(1)./(dx(1)+dx(2))];

ia(2:N-1)=[2:N-1];	
ja(2:N-1)=[2:N-1];	
za(2:N-1)=[-2./dx(2:N-1)./dx(1:N-2)];   % supported boundary

ia(N)=[N];	
ja(N)=[N];	
za(N)=[2./dx(N-1)./(dx(N-1)+dx(N-2))];

% 1st lower diagonal

ia(N+1:2*N-2)=[2:N-1];	
ja(N+1:2*N-2)=[1:N-2];	
za(N+1:2*N-2)=[2./dx(1:N-2)./(dx(1:N-2)+dx(2:N-1))];

ia(2*N-1)=[N];	
ja(2*N-1)=[N-1];	
za(2*N-1)=[-2./dx(N-1)./dx(N-2)];

% 2nd lower diagonal

ia(2*N)=[N];	
ja(2*N)=[N-2];	
za(2*N)=[2./(dx(N-1)+dx(N-2))./dx(N-2)];

% 1st upper diagonal

ia(2*N+1)=[1];	
ja(2*N+1)=[2];	
za(2*N+1)=[-2./dx(1)./dx(2)];

ia(2*N+2:3*N-1)=[2:N-1];	
ja(2*N+2:3*N-1)=[3:N];	
za(2*N+2:3*N-1)=[2./dx(2:N-1)./(dx(1:N-2)+dx(2:N-1))];

% 2nd upper diagonal

ia(3*N)=[1];	
ja(3*N)=[3];	
za(3*N)=[2./(dx(1)+dx(2))./dx(2)];

% now making the sparce matrix
Ld=sparse(ia,ja,za);

