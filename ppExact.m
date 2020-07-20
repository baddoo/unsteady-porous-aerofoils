function exactSol = ppExact(xVar,a,psi0,alph,bet)

exactSol = -4./sqrt(1+(heaviside(xVar-a).*psi0).^2).*(alph + bet*(1+xVar-(1-a)*atan(psi0)/pi))...
                .*sqrt((1-xVar)./(1+xVar)).*abs((xVar-a)./(1-xVar)).^(atan(psi0)/pi);

end