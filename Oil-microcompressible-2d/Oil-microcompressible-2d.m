%2D single phase with 4 production wells

%Basic Parameters, SI unit
clear all
dx = 5*0.3048;
dy = 5*0.3048;
dz = 5*0.3048;
nx = 20;
ny = 20;
beta_c = 1;
alpha_c = 1;
W = [5,5;5,15;15,5;15,15];      % well position
p0 = 5500*6894.76;      % original pressure 
por0 = 0.1*rand(nx*ny,1) + 0.3;
kx = (0.2*rand(nx*ny,1) + 0.4)*1e-12;       % random permeability in x direction
ky = (0.2*rand(nx*ny,1) + 0.4)*1e-12;       % random permeability in y direction
q = 200*0.3048^3/86400;     % m^3/s

p_mu = [4000;5000;5500;6000;7000;6500;7500;8000]*6894.76;
mu = [0.9100;0.9200;0.9243;0.9372;0.9650;0.9494;0.9812;1.0019]*1e-3;
B0 = 1.1;

c_phi = 1e-7;       % 地层压缩系数
c = 1e-9;           % 流体压缩系数

nt = 60*24*2;
dt = 60;

% number the grids
Grid = zeros();
N = 1;
for i = 1:nx
    for j = 1:ny
        Grid(j,i)= N;
        N = N+1;
    end
end

%Calculation
p = ones(nx*ny,nt)*p0;
T = zeros(nx*ny,nt);
G = zeros(nx*ny,nt);        % Gamma index
X = zeros(nx*ny,nx*ny);     % matrix
Y = zeros(nx*ny,1);        % right vector
u = ones(ny*nx,nt)*0.9E-3;
B = ones(ny*nx,nt)*1.1;
por = ones(ny*nx,nt)*0.3;
p_grid = zeros(ny,nx,nt);
p_grid(:,:,1) = reshape(p(:,1),ny,nx);

for t = 1:nt
    %u(:,t) = (2*10^(-5)*p(:,t)+0.8)*1E-3;      % viscosity
%     u(:,t) = interp1(p_mu,mu,p(:,t),'linear')   ;  % viscosity
    B(:,t) = B0./(1 + c_phi.*(p(:,t)-p0));   % fluid volume index
    por(:,t) = por0.*(1 + c_phi.*(p(:,t) - p0));
    T(:,t) = dz*beta_c*kx./u(:,t)./B(:,t);
    G(:,t) = dx*dy*dz/alpha_c/dt*(por0*c_phi./B(:,t) + por(:,t)*c./B0);
    
    Y = -G(:,t).*p(:,t);
    Y(Grid(W(1,1),W(1,2))) = q + Y(Grid(W(1,1),W(1,2)));
    Y(Grid(W(2,1),W(2,2))) = q + Y(Grid(W(2,1),W(2,2)));
    Y(Grid(W(3,1),W(3,2))) = q + Y(Grid(W(3,1),W(3,2)));
    Y(Grid(W(4,1),W(4,2))) = q + Y(Grid(W(4,1),W(4,2)));
    
    for i = 1:nx
        for j = 1:ny
            if p(Grid(i,j),t) >= 4000*6894.76 && p(Grid(i,j),t) <= 8000*6894.76
                u(Grid(i,j),t) = interp1(p_mu,mu,p(Grid(i,j),t),'linear')   ;  % viscosity
            elseif p(Grid(i,j),t) < 4000*6894.76
                u(Grid(i,j),t) = 0.91*1E-3;
            else
                u(Grid(i,j),t) = 1.0019*1E-3;
            end
            if mod(i,nx) == 1
                if mod(j,ny) == 1
                    X(Grid(j,i),((i-1)*ny+j)) = -2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t)) ...
                        -2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),i*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t));
                    X(Grid(j,i),((i-1)*ny+j+1)) = 2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) ;
                elseif mod(j,ny) == 0
                    X(Grid(j,i),((i-1)*ny+j)) =  -2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t))...
                        -2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),i*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t));
                    X(Grid(j,i),((i-1)*ny+j-1)) = 2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t));
                else
                    X(Grid(j,i),((i-1)*ny+j)) =  -2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t)) -...
                        2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t)) - ...
                        2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),i*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t));
                    X(Grid(j,i),((i-1)*ny+j-1)) = 2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t));
                    X(Grid(j,i),((i-1)*ny+j+1)) = 2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t));
                end
            elseif mod(i,nx) == 0
                if mod(j,ny) == 1
                    X(Grid(j,i),((i-1)*ny+j)) = -2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t)) ...
                        -2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),(i-2)*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t));
                    X(Grid(j,i),((i-1)*ny+j+1)) = 2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) ;
                elseif mod(j,ny) == 0
                    X(Grid(j,i),((i-1)*ny+j)) =  -2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t))...
                        -2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),(i-2)*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t));
                    X(Grid(j,i),((i-1)*ny+j-1)) = 2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t));
                else
                    X(Grid(j,i),((i-1)*ny+j)) =  -2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t)) -...
                        2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t)) - ...
                        2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) - G(Grid(j,i),t);
                   X(Grid(j,i),(i-2)*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t));
                    X(Grid(j,i),((i-1)*ny+j-1)) = 2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t));
                    X(Grid(j,i),((i-1)*ny+j+1)) = 2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t));
                end
            else
                if mod(j,ny) == 1
                    X(Grid(j,i),((i-1)*ny+j)) = -2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t)) - ...
                        2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t)) -...
                        2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),i*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t));
                    X(Grid(j,i),(i-2)*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t));
                    X(Grid(j,i),((i-1)*ny+j+1)) = 2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) ;
                elseif mod(j,ny) == 0
                    X(Grid(j,i),((i-1)*ny+j)) =  -2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t)) - ...
                        2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t)) -...
                        2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),i*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t));
                    X(Grid(j,i),(i-2)*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t));
                    X(Grid(j,i),((i-1)*ny+j-1)) = 2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t));
                else
                    X(Grid(j,i),((i-1)*ny+j)) =  -2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t)) - ...
                        2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t)) -...
                        2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t)) - ...
                        2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t)) - G(Grid(j,i),t);
                    X(Grid(j,i),i*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i+1),t)/(T(Grid(j,i),t)+T(Grid(j,i+1),t));
                    X(Grid(j,i),(i-2)*ny+j) = 2*T(Grid(j,i),t)*T(Grid(j,i-1),t)/(T(Grid(j,i),t)+T(Grid(j,i-1),t));
                    X(Grid(j,i),((i-1)*ny+j-1)) = 2*T(Grid(j,i),t)*T(Grid(j-1,i),t)/(T(Grid(j,i),t)+T(Grid(j-1,i),t));
                    X(Grid(j,i),((i-1)*ny+j+1)) = 2*T(Grid(j,i),t)*T(Grid(j+1,i),t)/(T(Grid(j,i),t)+T(Grid(j+1,i),t));
                end
            end
        end
    end
    
    p(:,t+1) = X\Y;
    p_grid(:,:,t+1) = reshape(p(:,t+1),ny,nx);
end

[m,n] = meshgrid(1:20);
surf(m,n,p_grid(:,:,end))       % final 3D plot
