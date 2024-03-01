%-----------------------------------------------------------------------
% FD NavierStokes in 2D - - 30 AUG 2019
%-----------------------------------------------------------------------

tic
clear all;
lx = 1;       % width of box
ly = 1;       % height of box
nx = 32;      % number of x-gridpoints
ny = nx;      % number of y-gridpoints
mu=0.01;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
xm = linspace(0.5*hx,lx-0.5*hx,nx); % x - coordinate
ym = linspace(0.5*hy,ly-0.5*hy,ny); % y - coordinate
[X,Y] = meshgrid(x,y);
dt=hx/8;
t_start=0;
t_final=dt*100;
%-----------------------------------------------------------------------
% boundary conditions
U = zeros(nx-1,ny); V = zeros(nx,ny-1);

uN = x*0+1;uS = x*0;
uE = ym*0;uW = ym*0;
vN = xm*0;vS = xm*0;
vE = y*0;vW = y*0;
%------------------------------------------------------------------------
%RHS
Ubc = mu*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hy^2+[uW;zeros(nx-3,ny);uE]/hx^2);
Vbc = mu*([vS' zeros(nx,ny-3) vN']/hy^2+[2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hx^2);
divbc=reshape([-vS' zeros(nx,ny-2) vN']/hy+[-uW;zeros(nx-2,ny);uE]/hx,[],1);

%Operators
Lu = speye((nx-1)*ny)/(dt)+mu*(kron(speye(ny),D(nx-1,hx,2))+kron(D(ny,hy,3),speye(nx-1))); %Laplacian operator for U
Lv = speye((ny-1)*nx)/(dt)+mu*(kron(speye(ny-1),D(nx,hx,3))+kron(D(ny-1,hy,2),speye(nx))); %Laplacian operator for V
Lp = -(kron(speye(ny),D(nx,hx,1))+kron(D(ny,hy,1),speye(nx))); Lp(1,1) = 3/2*Lp(1,1); %Laplacian operator for P
GX=kron(speye(ny),spdiags([-ones(nx-1,1) ones(nx-1,1)], 0:1,nx-1,nx))/hx; %Gradient operator in x
GY=spdiags([-ones(nx*ny,1) ones(nx*ny,1)], 0:nx:nx,nx*(ny-1),nx*ny)/hy; %Gradient operator in y

for t=t_start:dt:t_final
%Advection
Ue = [uW;U;uE]; Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
Ve = [vS' V vN']; Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
Uc = inter(Ue(:,2:end-1)); Ux = diff(Uc)/hx;
Ua = inter(Ue(2:end-1,:)')'; Uy = diff(Ua')'/hy;
Vc = inter(Ve(2:end-1,:)')';Vy = diff(Vc')'/hy;
Va = inter(Ve(:,2:end-1));Vx = diff(Va)/hx;
Vb = inter(Vc);Ub = inter(Uc')';
Uadvect = U.*Ux + Vb.*Uy;
Vadvect = Ub.*Vx + V.*Vy;    
    
%RHS
rhsu = reshape(U/dt+Ubc-Uadvect,[],1);
rhsv = reshape(V/dt+Vbc-Vadvect,[],1);
Uo=U;Vo=V;

%Solvers
U = Lu\rhsu;
V = Lv\rhsv;
P = Lp\(divbc-GX'*U-GY'*V);
U = U - GX*P;
V = V - GY*P;

%Solutions
U = reshape(U,nx-1,ny);
V = reshape(V,nx,ny-1);
P = reshape(P,nx,ny);

subplot(2,2,1)
contourf(x(2:end-1),ym,U');
axis equal, axis([0 lx 0 ly]);
subplot(2,2,2)
contourf( xm,y(2:end-1),V');
axis equal, axis([0 lx 0 ly]);
subplot(2,2,3)
contourf( xm,ym,P');
axis equal, axis([0 lx 0 ly]);
subplot(2,2,4)
Ue = [uS' inter([uW;U;uE]')' uN'];
Ve = [vW;inter([vS' V vN']);vE];
Len = sqrt(Ue.^2+Ve.^2+eps);
hold on; contourf(x,y,Len'); quiver(x,y,(Ue./Len)',(Ve./Len)',0.3,'k');hold off
axis equal, axis([0 lx 0 ly]);
drawnow
end
toc

function L = D(n,h,b)
L = spdiags([-1 b 0;ones(n-2,1)*[-1 2 -1];0 b -1],-1:1,n,n)'/h^2;
end
function Ac = inter(A)
Ac = (A(2:end,:)+A(1:end-1,:))/2;
end


