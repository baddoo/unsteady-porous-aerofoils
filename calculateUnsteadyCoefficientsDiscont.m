function solStruct = calculateUnsteadyCoefficientsDiscont(struct)

solStruct = struct;
z = struct.z;
dzdx = struct.dzdx;

k = struct.k;
na = struct.na;
nf = struct.nf;

a = struct.junction;

wa = @(xVar) 1i*k*z(xVar)+dzdx(xVar);

rhoe = struct.rhoe;
Phifun = struct.Phifun;
psiFun = @(x) eps + 4./(Phifun(x) + 2i*k.*rhoe(x));

fN = @(xVar) -2*wa(xVar);      
fC= @(xVar)  -1i*k*exp(-1i*k*xVar)/pi.*expint(-1i*k*(xVar-1));

leSing = 1/2;
teZero = 1/pi*acot(psiFun(1));
midZero = 1/pi*(acot(psiFun(a-eps))-acot(psiFun(a+eps)));
%if midZero<0 || midZero>.5; warning('The interior zero is not in the correct range. Check that the porosity profile is not negative');end

solStruct.teZero = teZero;
solStruct.leSing = leSing;
solStruct.midZero = midZero;

tcf = real(myJacobiNodes(nf,midZero,1/2)); % Generate forward collocation points
tca = real(myJacobiNodes(na,teZero,midZero)); % Generate aft collocation points

t2fun = @(tVar) -1 + (tVar - 1)*(1+a)/(1-a);
t1fun = @(tVar)  1 + (tVar + 1)*(1-a)/(1+a);

xt1 = @(tVar) -1 + (tVar + 1)*(1+a)/2;
xt2 = @(tVar)  1 + (tVar - 1)*(1-a)/2;

% Define Lambda coefficients
lamCoefs = defineLambdaCoefs(solStruct);
omegCoefs= exp(1i*k).*lamCoefs/(-1i*k*2^-midZero);
% PI coefficients
piCoefs = definePiCoefs(solStruct);

% Set up forward section
jQ0f = weight(tcf,midZero,-1/2).*myJacobiQ2(0,midZero,-1/2,tcf);
jI0f = myJacobiI(0,midZero,-1/2,tcf)*(1+a)/2;
jQAf = weight(tcf,midZero,1/2).*myJacobiQ2(nf-2,midZero,1/2,tcf);
jPAf = weight(tcf,midZero,1/2).*myJacobiP(nf,nf-2,midZero,1/2,tcf);
jIAf = myJacobiI(nf-2,midZero,1/2,tcf)*(1+a)/2;
jQBf = weight(t2fun(tcf),teZero,midZero).*myJacobiQ2(na-1,teZero,midZero,t2fun(tcf));
jQB0f = weight(t2fun(tcf),0,midZero).*myJacobiQ2(0,1e-8,midZero,t2fun(tcf));
jQC0f = weight(tcf,0,1/2).*myJacobiQ2(0,1e-8,1/2,tcf);
jIC0f = myJacobiI(0,0,1/2,tcf)*(1+a)/2;
jQD0f = weight(t2fun(tcf),teZero,0).*myJacobiQ2(0,teZero,1e-8,t2fun(tcf));

Ff = [ psiFun(xt1(tcf)).*(weight(tcf,midZero,-1/2) + 1i*k*jI0f)-jQ0f/pi, psiFun(xt1(tcf)).*(jPAf + 1i*k*jIAf)-jQAf/pi,-jQBf/pi] ...
    +lamCoefs.*-jQB0f/pi...
    -omegCoefs.*fC(xt1(tcf))...
    +piCoefs.*((psiFun(xt1(tcf)).*(weight(tcf,0,1/2) + 1i*k*jIC0f)-jQC0f/pi)/sqrt(2) - jQD0f/2^teZero/pi) ;
    
% Set up aft section
jQ0a = weight(t1fun(tca),midZero,-1/2).*myJacobiQ2(0,midZero,-1/2,t1fun(tca));
jQAa = weight(t1fun(tca),midZero, 1/2).*myJacobiQ2(nf-2,midZero,1/2,t1fun(tca));
jPAa = weight(tca,teZero,midZero).*myJacobiP(na,na-1,teZero,midZero,tca);
jQBa = weight(tca,teZero,midZero).*myJacobiQ2(na-1,teZero,midZero,tca);
jQCa = weight(t1fun(tca),0,1/2).*myJacobiQ2(0,1e-8,1/2,t1fun(tca));
jQDa = weight(tca,teZero,0).*myJacobiQ2(0,teZero,1e-8,tca);
jID  = myJacobiI(0,teZero,0,tca)*(1-a)/2;  

jI  = myJacobiI(na-1,teZero,midZero,tca)*(1-a)/2;  
jIop = myJacobiI(0,0,midZero,tca)*(1-a)/2;
jQop = weight(tca,0,midZero).*myJacobiQ2(0,1e-8,midZero,tca); % address the fact that we can't use zero parameter

Fa = [-jQ0a/pi, -jQAa/pi, psiFun(xt2(tca)).*(jPAa + 1i*k*jI)-jQBa/pi]...
    +lamCoefs.*(psiFun(xt2(tca)).*(weight(tca,0,midZero)+1i*k*jIop) - jQop/pi)... % Put in integral
    -omegCoefs.*fC(xt2(tca))...
    +piCoefs.*((psiFun(xt2(tca)).*(weight(tca,teZero,0) + 1i*k*jID) -jQDa/pi)/2^teZero -jQCa/pi/sqrt(2))...
    -piCoefs.*psiFun(xt2(tca));

% Now invert to get coefficients
F = [Ff;Fa];
rhsFun = fN([xt1(tcf);xt2(tca)]);

solStruct.coefs  = F\rhsFun;

end