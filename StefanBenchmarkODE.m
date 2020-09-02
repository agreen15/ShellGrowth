function H = StefanBenchmarkODE(Y,nz,K,k,L,Cp,Ht,Tm,rho,Ts)
%Solving conduction equation with custom finite difference algorithms
Y=Y';
zb=Y(end); %Convective layer thickness
T=Y(1:end-1);
Tm=T(end);
z=linspace(0,zb,nz);
[Q,z2]=FluxNumerical(T,z);

Qbot=k(Tm)*Q(end);


%[Qv,zv2]=FluxNumerical(Tv,zv);
%Qt=k(Tc)*Qv(1);
Qb=k(Tm)*Q(end);


[H,z3]=FluxNumerical(Q,z2);
dzmdt=Qb/(rho*L);
dTdz=CentralDifferenceNumerical(T,z);
H=K*H+dzmdt.*[2:(nz-1)].*dTdz/nz; %Conduction eq.
%Movement of freezing front
H=[0,H,0,dzmdt]';