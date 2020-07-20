% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion

% Set angle of attack and parabolic camber (only use alp=0 for now)
alph = 0; bet = -.1;

% Set porosities for evaluation
delVec = [0, .1, .2, .3, .4];
dNum = numel(delVec);

% Defines geometry and derivative
zGeo = @(xVar) -bet/2*(xVar.^2-1);
dzdx = @(xVar) -alph - bet*xVar;

% Sets the number of Chebyshev polynomials to use in our series
N = 10;
n =  1:N;
m = (1:N+1).';

% Set up grid for integration
nx = 150; % Number of points on sine grid
xInt = sin(linspace(-1,1,nx+2)*pi/2).'; xInt(1)=[]; xInt(end)=[];
xInt3 = permute(xInt,[3,2,1]); % Permutes grid into third dimension

% Calcualte the 2nd Chebyshev coefficients of the mean line
vCoefs = 2/pi*trapz(xInt,dzdx(xInt3).*sqrt(1-xInt3.^2).*chebU(xInt3,m-1),3);

%% Calculate the p-coefficients
p = zeros(nx,dNum);
pCoefs = zeros(N+1,dNum);

R = @(x) 1;

tic
profile on
for k = 1:dNum

    delta = delVec(k);

    psi = @(x) 2*delta*R(x);

    leSing = 1/pi*acot(psi(-1)); % Leading edge singularity
    teZero = 1/pi*acot(psi( 1)); % Trailing edge zero

    jP  = myJacobiP(nx,N-1,teZero,teZero,xInt);
    jP3 = permute(jP,[3,2,1]);
    
    jQ  = myJacobiQ2(N-1,teZero,teZero,xInt);
    jQ3 = permute(jQ,[3,2,1]);

    f = psi(xInt3).*jP3 - 1/pi*jQ3;

    jPVar = myJacobiP(nx,1,teZero,-leSing,xInt);
    jPVar3 = permute(jPVar,[3,2,1]);
    jQVar = myJacobiQ2(1,teZero,-leSing,xInt);
    jQVar3 = permute(jQVar,[3,2,1]);

    fVar0 = psi(xInt3).*jPVar3(:,1,:) - 1/pi*jQVar3(:,1,:);

    Fker = weight(xInt3,teZero,teZero).*sqrt(1-xInt3.^2).*chebU(xInt3,m-1).*f;
    Fker = padarray(Fker,[0,1,0],0,'pre');
    Fker(:,1,:)= weight(xInt3,leSing,-teZero).*sqrt(1-xInt3.^2).*chebU(xInt3,m-1).*fVar0;

    F = trapz(xInt,Fker,3);

    pCoefs(:,k) = 2*pi*(F\vCoefs);
    p(:,k) = sum(pCoefs(2:end,k).'.*jP.*weight(xInt,teZero,teZero),2) + weight(xInt,teZero,-leSing)*pCoefs(1,k);

end
profile off
disp(['The loop took ',num2str(toc),' seconds.'])  

%% Plot the solution compared to the analytic solution

% Calculate grid with fewer points than integration
xPlot = sin(linspace(-1,1,10)*pi/2).'; xPlot(1) = []; xPlot(end) = []; 

cols = lines; % Get default colours
k2d = 1/pi*acot(2*delVec); % Calculates k(2 \delta)
legStr=[]; % Prepares string

figure(1)
clf
clear pl;
for k = 1:dNum
    pl(k)=plot(xInt,-p(:,k)/bet,'LineWidth',3,'Color',cols(k,:)); %#ok<SAGROW>
    hold on
    scatter(xPlot,4/sqrt(1+4*delVec(k).^2).*(xPlot+2*k2d(k)).*weight(xPlot,k2d(k),-k2d(k)),...
        'Marker','^','MarkerEdgeColor','k','LineWidth',2,'MarkerFaceColor',cols(k,:))
    legStr = [legStr,strcat("$\delta = ", num2str(delVec(k)),"$")]; %#ok<AGROW>
end

axis([-1,1,-5,5]);
xlabel('non-dimensional distance along chord, $x$','Interpreter','latex')
ylabel('$- p / \beta$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',18)
grid on
legend(pl.',legStr,'Interpreter','latex','Location','southeast');

hold off

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
Zminx3 = bsxfun(@minus,Z,xInt3);
compVel = exp(-1i*alph) + 1/(2i*pi)*trapz(xInt,gam3./Zminx3,3);
figure(3)
clf
pcolor(X,Y,real(compVel)); shading interp; colormap jet; colorbar;
caxis([.9,1.1]);
hold on
plot(xInt,-zGeo(xInt),'LineWidth',3,'Color','k')
streamNum = 20;
starty2=linspace(ymin,ymax,streamNum); starty2(end)=[];
startx2 = xmin*ones(size(starty2));
hlines2=streamline(X,Y,real(compVel),-imag(compVel),startx2,starty2);
set(hlines2,'Color',[1, 1, 1],'LineWidth',2);
title('Streamlines and horizontal velocity','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',18)

hold off
 