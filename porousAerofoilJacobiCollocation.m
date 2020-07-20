% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion

delVec = linspace(0,1,7); delVec(1)=[];
dNum = numel(delVec);
    % Set angle of attack and parabolic camber (only use alp=0 for now)
alph = 0.2; bet = 0;

% Set porosities for evaluation

% Defines geometry and derivative
zGeo = @(xVar) -bet/2*(xVar.^2-1);
dzdx = @(xVar) -alph - bet*xVar;


%% Calculate the p-coefficients
numN = 20;
error = zeros(numN-1,dNum);

R = @(x) (1+x).^2;

tic
pCell = cell(numN,dNum);

for k = 1:dNum
    

for N = 1:(numN)

    delta = delVec(k);

    psi = @(x) 2*delta*R(x);

    leSing = 1/pi*acot(psi(-1)); % Leading edge singularity
    teZero = 1/pi*acot(psi( 1)); % Trailing edge zero

    xCol = myJacobiNodes(N+1,teZero,1-leSing); % Generate collocation points

    jP  = myJacobiP(N+1,N-1,teZero,1-leSing,xCol);  
    jQ  = myJacobiQ2(N-1,teZero,1-leSing,xCol);
    jPVar = myJacobiP(N+1,1,teZero,-leSing,xCol);
    jQVar = myJacobiQ2(1,teZero,-leSing,xCol);

    f = weight(xCol,teZero,1-leSing).*(psi(xCol).*jP - 1/pi*jQ);
    fVar0 = weight(xCol,teZero,-leSing).*(psi(xCol).*jPVar(:,1) - 1/pi*jQVar(:,1));
    F = [fVar0,f];
    RHScol = -4*dzdx(xCol);
    pCoefs = F\RHScol;
    
    pCell{N,k} = @(zVar) cpNum(zVar.',pCoefs,teZero,leSing).';
    
    %k2d = 1/pi*acot(2*delta); % Calculates k(2 \delta)
    %analSol = @(zVar) 4/sqrt(1+4*delta.^2).*(alph+bet*(zVar+2*k2d)).*weight(zVar,k2d,-k2d);
    %integrand = @(zVar) abs(p(zVar)-analSol(zVar)).^2;
    %integrand = @(zVar) abs(p(zVar)-analSol(zVar)).^2;
    %error(N,k) = abs(integral(integrand,-1,1));
    %error(N,k) = abs(integrand(-0));
    N
end

end

for k = 1:dNum
   
for N = 1:(numN-1)

    %delta = delVec(k);
    %k2d = 1/pi*acot(2*delta); % Calculates k(2 \delta)
    %analSol = @(zVar) 4/sqrt(1+4*delta.^2).*(alph+bet*(zVar+2*k2d)).*weight(zVar,k2d,-k2d);
    intFun = @(zVar) abs(pCell{N+1,k}(zVar)-pCell{N,k}(zVar)).^2;
    norFun = @(zVar) abs(pCell{N+1,k}(zVar)).^2;
    %integrand = @(zVar) abs(p(zVar)-analSol(zVar)).^2;
    error(N,k) = sqrt(abs(integral(intFun,-1,1))/abs(integral(norFun,-1,1)));
    %error(N,k) = abs(integrand(-0));
    N
end

end


%hold on
%%
semilogy(1:(numN-1),error.')
hold off

return
%% Plot the solution compared to the analytic solution


% Calculate grid with fewer points than integration
xPlot = sin(linspace(-1,1,10)*pi/2).'; xPlot(1) = []; xPlot(end) = []; 

cols = lines; % Get default colours
legStr=[]; % Prepares string

figure(1)
clf
clear pl;


axis([-1,1,-5,5]);
xlabel('non-dimensional distance along chord, $x$','Interpreter','latex')
ylabel('$- p / \beta$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',18)
grid on
legend(pl.',legStr,'Interpreter','latex','Location','southeast');

hold off
return
%% Plot p-coefficients to show rapid decay

figure(2), clf
pl2=plot(abs(pCoefs),'LineWidth',3); 
xlabel('$n$','Interpreter','latex')
ylabel('$|p_{n-1}|$','Interpreter','latex')
title('Jump in pressure along chord','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',18)
title('Absolute value of pressure coefficients','Interpreter','latex')
lgd = legend(pl2.',legStr,'Interpreter','latex','Location','northeast');

%% Plot streamlines

presPlot = p(:,end);
gamPlot = -.5*presPlot;

xmin = -2; xmax = 2; nxp = 50;
ymin = -2; ymax = 2; nyp = 50;

xp = linspace(xmin,xmax,nxp);
yp = linspace(ymin,ymax,nyp);

[X,Y]=meshgrid(xp,yp);
Z = X+1i*Y;
gam3 = permute(gamPlot,[3,2,1]);
Zminx3 = bsxfun(@minus,Z,xCol);
compVel = exp(-1i*alph) + 1/(2i*pi)*trapz(xCol,gam3./Zminx3,3);
figure(3)
clf
pcolor(X,Y,real(compVel)); shading interp; colormap jet; colorbar;
caxis([.9,1.1]);
hold on
plot(xCol,-zGeo(xCol),'LineWidth',3,'Color','k')
streamNum = 20;
starty2=linspace(ymin,ymax,streamNum); starty2(end)=[];
startx2 = xmin*ones(size(starty2));
hlines2=streamline(X,Y,real(compVel),-imag(compVel),startx2,starty2);
set(hlines2,'Color',[1, 1, 1],'LineWidth',2);
title('Streamlines and horizontal velocity','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',18)

hold off
 