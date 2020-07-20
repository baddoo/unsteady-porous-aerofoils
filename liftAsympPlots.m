% Calculates the pressure distribution along a porous aerofoil 
% using the Jacobi polynomial expansion
addpath('matlab2tikz/src')
imageFolder = '../unsteady-jacobi-r1/unsteady-jacobi/images/';
LW = 'LineWidth';
INT = 'Interpreter';
FS = 'FontSize';
% Set angle of attack and parabolic camber (only use alp=0 for now)
%beta0 = -1; beta1 = 1;

% z = @(xVar) 1+(xVar-.5);%+ xVar-.5+ 0*(xVar - .5);%beta0/2 + beta1*xVar; 
% struct.z = z;
%dzdx = @(xVar) +0*xVar.^2; struct.dzdx = dzdx;
beta0 = 1; beta1 = 1/2; % pure pitching

z = @(xVar) beta0/2 + beta1*xVar; struct.z = z;
dzdx = @(xVar) (beta1 + 0*xVar); struct.dzdx = dzdx;
%% Calculate the p-coefficients
clf

profile on
nk = 30;
%kVec = [eps,logspace(-2,3,nk)]; kVec(1) = [];
kVec = logspace(-2,2,nk);
kPlot = [];%[eps,.001,.01,.1,1,10];
kVec = sort([kVec,kPlot]); kVec = unique(kVec); %kVec(end) = []; % remove last entry for check
%kVec = 50;
nk = numel(kVec);
locs = [];
for l=1:numel(kPlot); locs = [locs,find(kVec==kPlot(l),1)];end

na = 4;
fullLift = zeros(nk,na,2);
fullCirc = zeros(nk,na,2);
%myncLift = zeros(nk,na);
aVec = flip([-.5,-.25,0,.25,.5]);
psiFunCell = cell(na,1);
ap = zeros(1,na);
thet = zeros(1,na);
for l = 1:2
for m = 1:na
a = aVec(m);
for j = 1:nk
    k = kVec(j);
    struct.N = round(10+12*sqrt(k)); % change 5 to 12.
    struct.k = k;
    
    if l ==1
    PSIfun = @(x) eps + (1+x);%heaviside(real(xVar-a)).*((real(xVar)-a)./(1-min(aVec)));
    Phifun = @(x) 1./PSIfun(x);
    rhoe = @(x) m + 0*x;
    elseif l==2
            if m == 1
                PSIfun = @(x) eps + .5*(1+x);%heaviside(real(xVar-a)).*((real(xVar)-a)./(1-min(aVec)));
            elseif m ==2
                PSIfun = @(x) eps + cos(x*pi/2)/2;
            elseif m ==3
                PSIfun = @(x) eps + .5*max(0,x);
            elseif m == 4
                PSIfun = @(x) eps + exp(-1./(x+1));
            end
    Phifun = @(x) 1./PSIfun(x);
    rhoe = @(x) 1.2 + 0*x;        
    end
    struct.Phifun = Phifun;
    struct.rhoe = rhoe;
    
    solStruct = calculateUnsteadyCoefficients(struct); 

    [fLift,nLift,qLift] = lift(solStruct);
    [fCirc,nCirc,qCirc] = circulation(solStruct);
    fullLift(j,m,l) = fLift;
    fullCirc(j,m,l) = fCirc;
    
    disp([j,m])
end

end
end
profile off
%%
cols = linspecer(na);

for l = 1:2
figure(1)
clf
for m = flip(1:na)
loglog(kVec,abs(fullLift(:,m,l)),'Color',cols(m,:),LW,3)
hold on
end

locsU = find(kVec>1);
locsL = find(kVec<1);

if l==1
    loglog(kVec(locsU),.1*abs(kVec(locsU).^2),'k--',LW,1)
loglog(kVec(locsL),10+0*kVec(locsL),'k--',LW,1)
    text(.1,.1e5,{'same $\Phi$', 'different $\rho_e$'},INT,'Latex',FS,20,'HorizontalAlignment','center');
text(20,5,'$k^2$',INT,'Latex',FS,30,'HorizontalAlignment','center')
text(.1,20,'constant',INT,'Latex',FS,30,'HorizontalAlignment','center')
else
    loglog(kVec(locsU),.1*abs(kVec(locsU).^2),'k--',LW,1)
loglog(kVec(locsL),9+0*kVec(locsL),'k--',LW,1)
text(.1,.1e5,{'different $\Phi$', 'same $\rho_e$'},INT,'Latex',FS,20,'HorizontalAlignment','center');
text(20,2,'$k^2$',INT,'Latex',FS,30,'HorizontalAlignment','center')
text(.1,20,'constant',INT,'Latex',FS,30,'HorizontalAlignment','center')  
end

hold off
xlabel('$\log_{10} k$','Interpreter','latex')
ylabel('$\log_{10}C_L $','Interpreter','latex')
grid on; grid minor; grid minor;


cleanfigure;
matlab2tikz([imageFolder,'asympLift',num2str(l),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,...
             'extratikzpictureoptions','trim axis left, trim axis right');

figure(2)
clf
for m = flip(1:na)
loglog(kVec,abs(fullCirc(:,m,l)),'Color',cols(m,:),LW,3)
hold on
end

locsU = find(kVec>1);
locsL = find(kVec<1);

if l ==1
loglog(kVec(locsU),.5*abs(kVec(locsU)).^.5,'k--',LW,1)
loglog(kVec(locsL),2+0*kVec(locsL),'k--',LW,1)
text(.1,10,{'same $\Phi$', 'different $\rho_e$'},INT,'Latex',FS,20,'HorizontalAlignment','center');
text(10,1,'$\sqrt{k}$',INT,'Latex',FS,30,'HorizontalAlignment','center')
text(.1,2.5,'constant',INT,'Latex',FS,30,'HorizontalAlignment','center')
else
    loglog(kVec(locsU),.5*abs(kVec(locsU)).^.5,'k--',LW,1)
loglog(kVec(locsL),2.5+0*kVec(locsL),'k--',LW,1)
text(.1,7,{'different $\Phi$', 'same $\rho_e$'},INT,'Latex',FS,20,'HorizontalAlignment','center');
text(10,1,'$\sqrt{k}$',INT,'Latex',FS,30,'HorizontalAlignment','center')
text(.1,3,'constant',INT,'Latex',FS,30,'HorizontalAlignment','center')
end
hold off
xlabel('$\log_{10} k$','Interpreter','latex')
ylabel('$\log_{10} \Gamma $','Interpreter','latex')
%ylim([1.5,inf])
grid on; grid minor; grid minor;

cleanfigure;
matlab2tikz([imageFolder,'asympCirc',num2str(l),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
end
