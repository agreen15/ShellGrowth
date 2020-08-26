function Ht=DefineTidalHeat(roh,Tm,E,w,u,q,nm,R,HRad)


HM=u*(E^2)/w; %max dissipation rate;

%Y=Y0;
%T=Y(1:end-1);
%i=[1:nz];
%zm=Y(end); %thickness of lid (km)
%z=linspace(0,zm,nz);%zm.*(i./nz);
%dT=T(end)-T(1);
Yt=(q)/(R*Tm^2); 
A=nm*exp(Yt.*Tm)*w/u;


% X=A*exp(-Yt.*T); %viscosity
% Hm=2*HM./(X+1./X);
X=@(T)A*exp(-Yt.*T); %viscosity
Hm=@(X)(1*(2*HM./(X+1./X)))+HRad; % volumetric dissipation rate
Ht=@(T)Hm(X(T));%/roh; % dissipation rate per mass



% for i2=[1:nz]
%     Ht(i2)=(2.*Nz(i2)*E.^2)/(1+((w.*Nz(i2))/(u)).^2);
% end
%Ht=Hm./roh;
%semilogx(Hm(X(T)),-z)