%% PWEM Simulation

%% INITIALIZE 
clc
close all 
clear 


%% FIGURE 

figure('Color', 'w', 'Units', 'normalized', 'Outerposition', [0 0 1 1]); 

%% DASHBOARD 
L = 1; %% Lattice period
ur = 1; %% permeability 
er1 =9; %% dielectric permittivity
er2 = 1; %% permittivity
r = 0.4*L; %% radius of the unit cell
NG2X = 40; %% Resolution of the Band structures

%% SHAPES TO DETERMINE THEIR RESPECTIVE BANDS 
square = 0; 
triangle =0;
ring_hollow = 0;
ring_resonator = 1; 

%% BUILD A UNIT CELL (DASHBOARD)
N = 100; %% any number 
Nx = 2*N+1; %% Number of grid in x ( 2N + 1); INVERSION SYMMETRY 
Ny = Nx; %%  ---- in y axis
dx = L/Nx; %% stepsize in x-axis
dy = L/Ny; %% stepsize in y-axis
xa = [0:Nx-1]*dx;  xa = xa - mean(xa); 
ya = [0:Nx-1]*dy;  ya = ya - mean(ya);

%% CONDITION TO SELECT ONLY ONE LATTICE
if (square == 1 && triangle == 1 && ring_hollow ==1 && ring_resonator ==1) | ...
(square == 0 && triangle == 0 && ring_hollow ==0 && ring_resonator ==0)
    disp('Select one unit Lattice') 
end 

if square == 1
[Xa, Ya] = meshgrid(ya, xa); 
AA = sqrt(Xa.^2 + Ya.^2)<= r.^2; 
ER = er1 + (er2-er1)*AA; 
end

%ER = er1 + (er2-er1)*AA; 
UR = ur.*ones(Nx, Ny); 

%% TRIANGLE
if triangle == 1 
Sx =1;
Sy = 1; 
ShifterX = floor(Nx/10); 
ShifterY = floor(Ny/10); 

% width and height
w= .8*Sx; 
h= .9*Sy; 

dx = Sx/Nx; 
dy = Sy/Ny; 

ny = round(h/dy); %number of cells for the object

ny1 = ShifterY+ floor((Ny-ny)/2); % centering of object [index]
ny2 = ny1 + ny - ShifterY; % index2

ER = zeros(Nx, Ny); 

for ny= ny1:ny2
    ffactor = (ny-ny1+ShifterY)/(ny2-ny1+ShifterY); % function for different fractions 
    nx = round(ffactor*w/dx); 
    nx1 = ShifterX + floor((Nx-nx)/2); 
    nx2 = nx1 + nx -ShifterX; 
    ER(nx1:nx2, ny) = 1; 
end 
% Refractive index 
ER =  er1 + (er2-er1)*ER'; 
end 

if ring_hollow == 1 
c0 = 3e-8; 
N = [Nx, Ny];    
xrange = [-2 2];  % x boundaries in L0
yrange = 1*[-2 2];  % y boundaries in L0
wvlen_scan = linspace(1,2.6, 20);
wvlen = 1.7;
k0 = 2*pi/wvlen;
omega_p = 0.72*pi*1e15;%3e15; %omega_p was 3 which gave us good results...
gamma = 400e12; %20e12; % (5.5e12 is the default)
omega = 2*pi*c0/wvlen*1e6;
%epsilon_diel = 16;
epsilon_metal =  1 - omega_p^2./(omega^2-1i*gamma*omega); %MINUS SIGN FOR FDFD
epsilon_diel = epsilon_metal;
thickness = 0.15;
fill_factor = 0.1; %half metal, half dielectric

delta_arc = 6*pi/180;
inner_radius = 0.5; outer_radius = 0.7;
eps = ones(N);
eps = curved_stripe(eps, N,xrange, yrange, ...
    inner_radius, outer_radius, delta_arc, epsilon_metal, epsilon_diel);

delta_arc_2 = 3*pi/180;
inner_rad_2 = 1.3; outer_rad_2 = 1.5;
eps = curved_stripe(eps, N,xrange, yrange, ...
    inner_rad_2, outer_rad_2, delta_arc_2, epsilon_metal, epsilon_diel);
ER = er1 + (er2-er1)*eps; 
% figure(2); visabs(eps, xrange, yrange);
% drawnow();
end 

if ring_resonator == 1 
    N = [Nx, Ny]; 
    inner_rad = 60; outer_rad = 90;
    xc = (round(N(1)/2)); yc = (round(N(2)/2));
    xa = -xc+1:xc; ya = -yc+1:yc;
    [X,Y] = meshgrid(xa,ya);
    eps_ring = ((X.^2+Y.^2)<outer_rad^2);
    eps_inner = ((X.^2+Y.^2)<inner_rad^2);
    eps_ring = eps_ring-eps_inner;
    eps_ring(eps_ring == 1) = er1;
    eps_ring(eps_ring == 0) = er2;
    ER = eps_ring; 
end 

%% COMPUTE CONVOLUTION MATRICES
P = 3;
Q = P;
ERC = convmat(ER, P, Q); 
URC = convmat(UR, P, Q); 


%% COMPUTE LIST OF BLOCH VECTORS
t1 = [L; 0]; 
t2 = [0; L]; 

T1 = 2*pi./t1;
T1(2)  = 0; 
T2 = 2*pi./t2; 
T2(1) = 0; 

G = [0; 0];
X = 0.5.*T1;
M = 0.5*T1 + 0.5*T2; 


% GENERATE LIST 
L1 = norm(G-X); % To find the length from Gamma to X to M in the reduced
L3 = norm(M-G); % Brillioun zone of a square lattice. 
L2 = norm(M-X); 

N1 = NG2X; 
N2 = round(N1*L2/L1); 
N3 = round(N1*L3/L1);

%% Bloch wave waves 
BX = [linspace( G(1), X(1), N1), linspace(X(1), M(1), N2), linspace(M(1), G(1), N3)];  
BY = [linspace( G(2), X(2), N1), linspace(X(2), M(2), N2), linspace(M(2), G(2), N3)];

%% REMOVE REDUNDANCY
BX([N1+1, N1+N2+1]) = []; 
BY([N1+1, N1+N2+1]) = [];
BETA = [BX; BY];

%% ITERATE OVER BLOCH WAVE VECTOR

NBETA = length(BETA(1,:)); 


KP =[ 1, N1, N1+N2-1, N1+N2+N3-2]; % keynotes

KL = { '\Gamma','X',  'M', '\Gamma'};
%% PWEM 

% COMPUTE SPATIAL HARMONIC INDICES  %% MILLER INDICES
p= [-floor(P/2): floor(P/2)]; %indices along x
q= [-floor(Q/2): floor(Q/2)]; % indices along y

for nb = 1: NBETA 
    
    bx = BETA (1, nb); 
    by = BETA (2, nb);
    
    KX = bx - p*T1(1)- q*T2(1);
    KY = by - p*T1(2) - q*T2(2);
    
   [KY, KX] = meshgrid(KY, KX); 
    KX = diag(sparse(KX(:))); 
    KY = diag(sparse(KY(:))); 
    
    %% BUILD EIGEN-VALUE PROBLEM AX = B
    
   A = KX/URC*KX+KY/URC*KY; 
   B = ERC; 
   
   A1 = KX/ERC*KX + KY/ERC*KY; 
   B1 = URC; 
   
   %% SOLVE GENERALIZED EIGEN-VALUE PROBLEM
   [V, D] =  eig(full(A), full(B));  
   D = real(sqrt(diag(D))); 
   KO(:, nb) = sort(D); 
   
    %% SOLVE GENERALIZED EIGEN-VALUE PROBLEM
   [V1, D1] =  eig(full(A1), full(B1));  
   D1 = real(sqrt(diag(D1))); 
   KO1(:, nb) = sort(D1); 
   % SHOW BANDS 
    subplot(1,5,1); 
    plot([1:nb], KO(:, 1:nb), 'b--', 'LineWidth', 2.4); 
    title('TE mode', 'Fontsize', 18); 
    xlim([1 NBETA]); 
    
    subplot(1,5,2); 
    plot([1:nb], KO1(:, 1:nb), 'b--', 'LineWidth', 2.4);
    title('TM mode', 'Fontsize', 18); 
    xlim([1 NBETA]); 
    drawnow; 
end 

% CALCULATE NORMALIZED FREQUENCY 
WN = L*KO/(2*pi);
WN1 = L*KO1/(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW BAND DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING OF UNIT CELL
subplot(1,5,3); 
hold on
imagesc(xa, ya, real(ER')); 
%line(BETA(1,:), BETA(2,:))
title('Structure', 'Fontsize', 18); 
colormap('jet'); 
colorbar; 
axis equal tight; 
% CLEAR FIGURE WINDOW 

% PLOT BANDS
subplot(1, 5, 4); 
plot([1:nb], WN, 'LineWidth', 2.5); 
title('TE mode', 'Fontsize', 18); 
% SET AXIS LIMITs
xlim([1 NBETA]); 
%ylim([0, 1]); 

% SET TICK MARKS 
set(gca, 'XTick',KP,'XTickLabel', KL); 

% for n = 1: length(KP)
%     line(KP(n)*[1 1], [0 1], 'Color', 'k', 'LineStyle', '-')
% end 

% LABEL AXES 
xlabel('Bloch Wave VErctor, $\vec{\beta}$', 'Interpreter', 'LaTEx'); 
ylabel('Normalized Frequency, $\omega$', 'Interpreter', 'LaTEx'); 

subplot(1, 5,5);
plot([1:nb], WN1, 'LineWidth', 2.5);
title('TM mode', 'Fontsize', 18); 
% SET AXIS LIMITs
xlim([1 NBETA]); 
%ylim([0, 1]); 

% SET TICK MARKS 
set(gca, 'XTick',KP,'XTickLabel', KL); 

% for n = 1: length(KP)
%     line(KP(n)*[1 1], [0 1], 'Color', 'k', 'LineStyle', '-')
% end 

% LABEL AXES 
xlabel('Bloch Wave VErctor, $\vec{\beta}$', 'Interpreter', 'LaTEx'); 
ylabel('Normalized Frequency, $\omega$', 'Interpreter', 'LaTEx'); 



