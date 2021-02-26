%
% k1 = linspace(0,1.25);
% k2 = linspace(0,1,50);
% tvec = sort([log(1+k1).*k1,1./k2]);
% plot(C(tvec),'x')
% return

na = 4;

%Cv = @(k) (C(k)-1)./k;
%C = @(k) real(C(k) - 1)./(k);
np = 40;
tp = [linspace(0,1,np/2),linspace(1,5,np/2)]'+1e-5;
wagner = zeros(np,na,2);
wagner2 = zeros(np,na,2);
tp1 = 2*tp;

for n = 1:2
for j = 1:na
%[r,pol,res] = aaa(C(kVec2),kVec2);
%[r,pol,res] = aaa((numTheo(1:end-5,j)-1)./kVec(1:end-5).',kVec(1:end-5));
[r,pol,res] = aaa(real(numTheo(1:end-5,j,n).*abs(qsLift(10,j,n)./kVec(10)))/2/pi,kVec(1:end-5));
%[r,pol,res] = aaa(real(numTheo(1:end,j,2)),kVec(1:end));

%pol = pol(abs(res)>1e-3);
%res = res(abs(res)>1e-3);
%pol = pol(1:end);
%res = res(1:end);
[max(imag(pol)),min(imag(pol))];

rp = @(t) r(Inf) + res.'*(1./(t-pol));
%cosint2 = @(t) cosint(t) + t-eulergamma - log(t);
%wagner(:,j) = + 2/pi*(real(r(Inf))*pi/2 + res.'*((pi+2*cosint(-pol.*tp.').*sin(pol.*tp.')-cos(pol.*tp.').*(pi+2*sinint(pol.*tp.')))./(-2*pol)));
wagner(:,j,n) = 2/pi*(real(r(Inf))*pi/2 + res.'*((pi+2*cosint(-pol.*tp.').*sin(pol.*tp.')-cos(pol.*tp.').*(pi+2*sinint(pol.*tp.')))./(-2*pol)));
%wagner(:,j) = 1+2/pi*real(res.'*(-cosint(-pol.*tp.').*sin(pol.*tp.')+.5*cos(pol.*tp.').*(pi+2*sinint(pol.*tp.'))));
[r,pol,res] = aaa([(wagner(:,j));abs(qsLift(10,j,n)./kVec(10))/2/pi],[tp;1e4]);
wagner2(:,j,n) = r(tp1);
end
%return
end

C = @(k) besselk(1,1i*k)./(besselk(0,1i*k) + besselk(1,1i*k));
%%

[r1,pol1,res1] = aaa(real(C(kVec(1:2:end))),kVec(1:2:end));
wagner1 = 2/pi*(real(r1(Inf))*pi/2 + res1.'*((pi+2*cosint(-pol1.*tp.').*sin(pol1.*tp.')-cos(pol1.*tp.').*(pi+2*sinint(pol1.*tp.')))./(-2*pol1)));
[r1,pol1,res1] = aaa([wagner1,1],[tp;1e4]);
wagner21 = r1(tp1);
%% Plots
for l = 1:2
for n = 1:2
    if l == 1
        fac = 1+0*(1:na);
    elseif l ==2
        fac = abs(qsLift(10,:,n)./kVec(10))/2/pi;
    end
LW = 'LineWidth';
figure(n)
    if n==1
    cols = flip(cmocean('matter',na+1));
    elseif n==2
    cols = flip(cmocean('speed',na+1));
    end
    if n ==1
    plot(tp1,wagner21,'k',LW,2)
    else
    plot(tp,wagner1,'k',LW,2)
    end
hold on
for m = 1:na
%plot(tp,wagner(:,m),'Color',cols(m,:),LW,2);
if n ==1
plot(tp1,wagner2(:,m,1)./fac(m),'Color',cols(m,:),LW,1);
elseif n ==2
plot(tp,wagner(:,m,2)./fac(m),'Color',cols(m,:),LW,1);
end

hold on
end
ylim([0,1])
if n ==1
xlim([0,10])
text(9.95,0.2,'Varying flow resistance, $\Phi$','HorizontalAlignment','right','BackgroundColor','w')
elseif n==2
xlim([0,5])
text(4.95,0.2,'Varying effective density, $\rho_e$','HorizontalAlignment','right','BackgroundColor','w')
end
hold off
grid on
xlabel('convective time, $t$')
if l==1
    ylabel('$L^C$')
elseif l==2
    ylabel('$\phi(t)$')
end

    cleanfigure;
    if l ==1
    matlab2tikz([imageFolder,num2str(n),'dimwagner.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
    elseif l ==2
    matlab2tikz([imageFolder,num2str(n),'wagner.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');    
    end
end
end