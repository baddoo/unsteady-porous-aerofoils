function solStruct = calculateUnsteadyCoefficientsNC(struct)

z = struct.z;
dzdx = struct.dzdx;

k = struct.k;
N = struct.N;

wa = @(xVar) 1i*k*z(xVar)+dzdx(xVar);

psiFun = struct.psiFun;
fN = @(xVar) -2*wa(xVar);    

leSing = 1/pi*acot(psiFun(-1)); % Leading edge singularity
teZero = 1/pi*acot(psiFun( 1)); % Trailing edge zero

xCol = myJacobiNodes(N+1,teZero,1-leSing); % Generate collocation points

jP  = myJacobiP(N+1,N-1,teZero,1-leSing,xCol);  
jI  = myJacobiI(N-1,teZero,1-leSing,xCol);  
jQ  = myJacobiQ2(N-1,teZero,1-leSing,xCol);
jPVar = myJacobiP(N+1,1,teZero,-leSing,xCol);
jIVar = myJacobiI(1,teZero,-leSing,xCol);
jIop = myJacobiI(1,teZero-1,1-leSing,xCol);
jQVar = myJacobiQ2(1,teZero,-leSing,xCol);
jQop =  myJacobiQ2(1,teZero-1,1-leSing,xCol); 

% Make this into a nice function
op =    psiFun(xCol).*(weight(xCol,teZero-1,1-leSing)        +1i*k.*jIop(:,1)) ...
        - weight(xCol,teZero-1,1-leSing)/pi.*jQop(:,1);
f =     psiFun(xCol).*(weight(xCol,teZero,1-leSing).*jP      +1i*k.*jI) ...
        - weight(xCol,teZero,1-leSing)/pi.*jQ;
f0 =    psiFun(xCol).*(weight(xCol,teZero, -leSing).*jPVar(:,1) +1i*k.*jIVar(:,1)) ...
        - weight(xCol,teZero, -leSing)/pi.*jQVar(:,1);

denom = -2^(1+teZero-leSing)*double(beta(sym(teZero),sym(2-leSing)));
omega0 = (2^(1+teZero-leSing)*double(beta(sym(1-leSing),sym(1+teZero))))/denom;
omega1 = (2^(2+teZero-leSing)*double(beta(sym(2-leSing),sym(1+teZero))))/denom;    

F = [f0    + omega0*op,...
    f(:,1) + omega1*op,...
    f(:,2:end)];

% Solve to find the coefficients
coefs  = F\fN(xCol);

% Calculate the coefficient of the wake
Lambda =  0;
Theta = omega0 * coefs(1) + omega1*coefs(2);

% Put the data into a structure
solStruct = struct;
solStruct.Lambda = Lambda;
solStruct.Theta = Theta;
solStruct.teZero = teZero;
solStruct.leSing = leSing;
solStruct.coefs  = coefs;

end