% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion

%% Calculate the p-coefficients
% Define the motion: this is for pitching about the leading edge
beta0 = 1; beta1 = beta0/2;
z = @(xVar) beta0/2 + beta1*xVar;
dzdx = @(xVar) beta1 + 0*xVar;

% Define the frequency
k = .5;

% Define the (reciprocal) flow resistivity
Phifun = @(x) 1./(0.05*(x + 1));

% Define the effective density. Here it is a constant
rhoe = @(x) 1.2 + 0*x;

% Set the number of Jacobi polynomials to use
N = 5; 

% Define the grid to plot the values on
xp = cos(flip(linspace(0,pi,1e2))'); xp(1) = []; xp(end) = [];
    
% Now put everything into a structure
struct.rhoe = rhoe;
struct.Phifun = Phifun;
struct.N = N; 
struct.k = k;
struct.z = z;
struct.dzdx = dzdx;

% Now calculate the solutions. This function calculates the solutions for
% the full pressure, non-circulatory pressure and the quasi-static pressure
solStruct = calculateUnsteadyCoefficients(struct); 

% Now compute the pressure from coefficients
[pressure,NCpressure,QSpressure] = presFun(xp,solStruct);
liftValue = lift(solStruct);
circulationValue = circulation(solStruct);

%% Plots
LW = 'LineWidth'; FS = 'FontSize';
INT = 'Interpreter'; LTX = 'Latex';

figure(1)
p1 = plot(xp,real(pressure),'k-',LW,2);
hold on
p2 = plot(xp,real(NCpressure),'r-',LW,2);
p3 = plot(xp,real(QSpressure),'b-',LW,2);
hold off
xlim([-1,1]);
ylim([-5,15])
grid on
xlabel('$x$',INT,LTX)
xlabel('$\Delta p$',INT,LTX)
legend([p1,p2,p3],{'pressure','non-circulatory pressure','quasi-static pressure'},FS,20,INT,LTX)