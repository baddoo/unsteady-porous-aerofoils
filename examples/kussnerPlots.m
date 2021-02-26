na = 4;

%Cv = @(k) (C(k)-1)./k;
%C = @(k) real(C(k) - 1)./(k);
np = 60;
tp = [linspace(0,1,np/2),linspace(1,10,np/2)]'+1e-5;
kussner = zeros(np,na);
S = @(kVar) -1i./(kVar.*(besselk(0,1i*kVar) + besselk(1,1i*kVar)));
for n = 1:2
for j = 1:na
%[r,pol,res] = aaa(C(kVec2),kVec2);
%[r,pol,res] = aaa((numTheo(1:end-5,j)-1)./kVec(1:end-5).',kVec(1:end-5));
%[r,pol,res] = aaa(real(numTheo(1:2:end-10,j)-1),kVec(1:2:end-10));:
[r,pol,res] = aaa(real(exp(-1i*kVec(5:5:end-15).').*numSearsNorm(5:5:end-15,j,n)),kVec(5:5:end-15));
rp = @(t) r(Inf) + res.'*(1./(t-pol));
kussner(:,j,n) = 2/pi*(real(r(Inf))*pi/2 + res.'*((pi+2*cosint(-pol.*tp.').*sin(pol.*tp.')-cos(pol.*tp.').*(pi+2*sinint(pol.*tp.')))./(-2*pol)));
%J = @(a,k) 1i/2*exp(1i*a*k).*(pi*sign(k) + 2i*cosint(-a.*abs(k)) + 2*sinint(a.*(k)));
%I = @(a,k) 1/2i*(J(a,1+k) - 0*J(a,1-k));
%I0= @(a,k)
%kussner(:,j) = 1 + 2/pi*real(r(Inf)*1i*atanh(tp.') ...
%            + (exp(1i*pol).*res).'*((pi+2*cosint(-pol.*tp.').*sin(pol.*tp.')-cos(pol.*tp.').*(pi+2*sinint(pol.*tp.')))./(-2*pol)));
%kussner(:,j) = 1 + 2/pi*real(res.'*((I(pol,tp.'))));
%wagner(:,j) = 1+2/pi*real(res.'*(-cosint(-pol.*tp.').*sin(pol.*tp.')+.5*cos(pol.*tp.').*(pi+2*sinint(pol.*tp.'))));
end
end
%return

[r,pol,res] = aaa(real(exp(-1i*kVec).*S(kVec)),kVec);
kussner1 = 2/pi*(real(r(Inf))*pi/2 + res.'*((pi+2*cosint(-pol.*tp.').*sin(pol.*tp.')-cos(pol.*tp.').*(pi+2*sinint(pol.*tp.')))./(-2*pol)));


%% Plots
fac = fullLift(1,:,:);
    LW = 'LineWidth';
for l =1:2
for n = 1:2
     if n==1
        cols = flip(cmocean('matter',nP+1));
    elseif n==2
        cols = flip(cmocean('speed',nP+1));
     end  
    
    if l ==1
        fac = 1+ 0*(1:nP);
    else
        fac = abs(fullLift(1,:,n)/4/pi);
    end
figure(n)    
plot(tp,kussner1,'k',LW,2)
grid on
hold on
    for j = 1:nP
        plot(tp,real(fac(j)*kussner(:,j,n)),'Color',cols(j,:),LW,1)
    end
hold off
ylim([0,1])
ylims=get(gca,'ylim');
xlims=get(gca,'xlim');
if n ==1
xlim([0,10])
text(9.95,0.2,'Varying flow resistance, $\Phi$','HorizontalAlignment','right','BackgroundColor','w')
elseif n==2
xlim([0,5])
text(4.95,0.2,'Varying effective density, $\rho_e$','HorizontalAlignment','right','BackgroundColor','w')
end

    if l==1
                ylabel('$\psi(t)$')
    elseif l==2
        ylabel('$L$')
    end
    xlabel('convective time, $t$')
        cleanfigure;
    if l ==1
    matlab2tikz([imageFolder,num2str(n),'kussner.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
    elseif l ==2
    matlab2tikz([imageFolder,num2str(n),'dimkussner.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');    
    end
    drawnow;
    hold off
    
end
end