function tubeplot(r,Omega3,Radius,RadPoints,TubeColor,StripeColor,transparency)

% plots a 3D tube about the curve defined by
% the vector r.  Radius is the user defined radius
% of the tube and RadPoints is the number of
% circumferential patches used to make the tube.

N=length(r(:,1));

% find tangent vector

dr = sqrt(   (r(2:N,1) - r(1:N-1,1)).^2  ...
           + (r(2:N,2) - r(1:N-1,2)).^2  ...
           + (r(2:N,3) - r(1:N-1,3)).^2 );

L1 = first_der(N,dr);       
L2 = second_der(N,dr); 

DrDs(:,1) = L1*r(:,1);
DrDs(:,2) = L1*r(:,2);
DrDs(:,3) = L1*r(:,3);

g = sqrt( DrDs(:,1).^2 + DrDs(:,2).^2 + DrDs(:,3).^2 );

e3(:,1) = DrDs(:,1) ./ g;
e3(:,2) = DrDs(:,2) ./ g;
e3(:,3) = DrDs(:,3) ./ g;

% define e3 and Omega3 at half grids

% e32 = 0.5.*( e3(1:N-1,:) + e3(2:N,:) );
e32(:,1) = ( r(2:N,1) - r(1:N-1,1) )./dr;
e32(:,2) = ( r(2:N,2) - r(1:N-1,2) )./dr;
e32(:,3) = ( r(2:N,3) - r(1:N-1,3) )./dr;

om32 = 0.5.*( Omega3(1:N-1) + Omega3(2:N) );

g2 = sqrt( e32(:,1).^2 + e32(:,2).^2 + e32(:,3).^2 );

e32(:,1) = e32(:,1)./g2;
e32(:,2) = e32(:,2)./g2;
e32(:,3) = e32(:,3)./g2;

% define e1(1) and e2(1)

if e3(1,3)>=1-eps
    
    e1(1,1) = 1;
    e1(1,2) = 0;
    e1(1,3) = 0;
    
else

    e1(1,1) = -e3(1,1).*e3(1,3)./sqrt(1-e3(1,3).^2);
    e1(1,2) = -e3(1,2).*e3(1,3)./sqrt(1-e3(1,3).^2);
    e1(1,3) = sqrt(1-e3(1,3)^2);

end

e2(1,1) = e3(1,2).*e1(1,3) - e3(1,3).*e1(1,2);
e2(1,2) = e3(1,3).*e1(1,1) - e3(1,1).*e1(1,3);
e2(1,3) = e3(1,1).*e1(1,2) - e3(1,2).*e1(1,1);

% find curvatures 

De3Ds(:,1) = L2*r(:,1) ;
De3Ds(:,2) = L2*r(:,2) ;
De3Ds(:,3) = L2*r(:,3) ;

% define half grid curvatures

curv_half(:,1) = 0.5*(De3Ds(1:N-1,1)+De3Ds(2:N,1)) ;
curv_half(:,2) = 0.5*(De3Ds(1:N-1,2)+De3Ds(2:N,2)) ;
curv_half(:,3) = 0.5*(De3Ds(1:N-1,3)+De3Ds(2:N,3)) ;

% integrate e1 and e2

for m=1:N-1
    
    om1(m) = -( De3Ds(m,1)*e2(m,1) + De3Ds(m,2)*e2(m,2) + De3Ds(m,3)*e2(m,3) );
    om2(m) =  ( De3Ds(m,1)*e1(m,1) + De3Ds(m,2)*e1(m,2) + De3Ds(m,3)*e1(m,3) );
    
% integrate e1 and e2 half step

    e1_half(m,:) =  e1(m,:) ...
                  + 0.5*dr(m).*(Omega3(m).*e2(m,:)-om2(m).*e3(m,:));
    e2_half(m,:) =  e2(m,:) ...
                  + 0.5*dr(m).*(om1(m).*e3(m,:)-Omega3(m).*e1(m,:));
 
                  
     m1(m) = sqrt( e1_half(m,1).^2 + e1_half(m,2).^2 + e1_half(m,3).^2);
     m2(m) = sqrt( e2_half(m,1).^2 + e2_half(m,2).^2 + e2_half(m,3).^2);
    
    e1_half(m,1) = e1_half(m,1) ./ m1(m);
    e1_half(m,2) = e1_half(m,2) ./ m1(m);
    e1_half(m,3) = e1_half(m,3) ./ m1(m);

    e2_half(m,1) = e2_half(m,1) ./ m2(m);
    e2_half(m,2) = e2_half(m,2) ./ m2(m);
    e2_half(m,3) = e2_half(m,3) ./ m2(m);
    
% define half grid omegas
 
        om1_half(m) = -( curv_half(m,1).*e2_half(m,1) ...
                       + curv_half(m,2).*e2_half(m,2) ...
                       + curv_half(m,3).*e2_half(m,3) );
                   
        om2_half(m) =  ( curv_half(m,1).*e1_half(m,1) ...
                       + curv_half(m,2).*e1_half(m,2) ...
                       + curv_half(m,3).*e1_half(m,3) );
        
% integrate full step

    e1(m+1,:) = e1(m,:) + dr(m).*(om32(m).*e2_half(m,:)-om2_half(m).*e32(m,:));
    e2(m+1,:) = e2(m,:) + dr(m).*(om1_half(m).*e32(m,:)-om32(m).*e1_half(m,:));

     m1(m+1) = sqrt( e1(m+1,1).^2 + e1(m+1,2).^2 + e1(m+1,3).^2);
     m2(m+1) = sqrt( e2(m+1,1).^2 + e2(m+1,2).^2 + e2(m+1,3).^2);
    
    e1(m+1,1) = e1(m+1,1) ./ m1(m+1);
    e1(m+1,2) = e1(m+1,2) ./ m1(m+1);
    e1(m+1,3) = e1(m+1,3) ./ m1(m+1);

    e2(m+1,1) = e2(m+1,1) ./ m2(m+1);
    e2(m+1,2) = e2(m+1,2) ./ m2(m+1);
    e2(m+1,3) = e2(m+1,3) ./ m2(m+1);
    
end

% define angle about the curve

theta = linspace(0,2*pi,RadPoints);

% % define arcs of fixed theta
% 
% for m = 1:RadPoints
%     arc1(:,:,m) = r + Radius.*(sin(theta(m)).*e1 + cos(theta(m)).*e2);
% 
%     hold on
% 
%     plot3(arc1(:,1,m),arc1(:,2,m),arc1(:,3,m))
% end
% 
% % define circles at fixed arclength
% 
% for m = 1:N
%     for i=1:3
%     arc2(:,i,m) = r(m,i) ...
%                   + Radius.*(sin(theta).*e1(m,i) + cos(theta).*e2(m,i));
%     end
%     hold on
% 
%     plot3(arc2(:,1,m),arc2(:,2,m),arc2(:,3,m))
% end

% define patches

for m=1:N
    for i=1:RadPoints
        for j=1:3
        
        vertices(RadPoints*(m-1)+i,j) = r(m,j) ...
                  + Radius(m).*(cos(theta(i)).*e1(m,j) ...
                  + sin(theta(i)).*e2(m,j));
        end
    end
end

% define connectivity matrix

vec1 = (1:RadPoints);
vec2 = circshift(vec1,[0 -1]);

for m=1:RadPoints
        
        MT(m,1) = vec1(m);
        MT(m,2) = vec2(m);
        MT(m,4) = vec1(m) + RadPoints;
        MT(m,3) = MT(m,2) + RadPoints;
        
end

for i=2:N-1
    
    for m=1:RadPoints
        
        M(m,1) = (i-1)*RadPoints+vec1(m);
        M(m,2) = (i-1)*RadPoints+vec2(m);
        M(m,4) = M(m,1) + RadPoints;
        M(m,3) = M(m,2) + RadPoints;
        
    end
    
    MT = [ MT ; M ];
end

%len = length(MT);

for k=1:N
    for l=1:RadPoints-2
    fsv(RadPoints*(k-1)+l,:) = TubeColor;
    end
    fsv(RadPoints*k-1,:) = StripeColor;
    fsv(RadPoints*k,:) = StripeColor;
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

