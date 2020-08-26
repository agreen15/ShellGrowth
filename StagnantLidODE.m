function H = StagnantLidODE(Y,nz,K,k,g,nm,L,Cp,Ht,Tm,rho,C,D,alph,beta,gamma,xi,zeta,Ts)
%Unpack Input Vector
Y=Y';
b=Y(end); %Convective layer thickness
zb=Y(end-1); %Base of conductive lid
zm=b+zb; %Base of shell
Tc=Y(end-2); %Temperature at base of lid
Qs=((k(Ts)*(Tc-Ts))/zb);
Qc=((k(Tc)*(Tm-Ts))/zm);
Nu=Qs/Qc; %Nusselt Number, convects if >=1
% if Nu<1;
%     %Conduction only solver
%     %T-profile for upper shell
%     T=Y(1:end-2);
%     z=linspace(0,zb,nz);
%     %T-profile for lower shell
%     Tv=linspace(Tc,Tm,nz);
%     zv=linspace(zb,zm,nz);
%     DTc=0;
%     %Solving conduction equation with custom finite difference algorithms
%     [Q,z2]=FluxNumerical(T,z);
%    
%     Qbot=k(Tc)*Q(end);
%     
%     
%     [Qv,zv2]=FluxNumerical(Tv,zv);
%     Qt=k(Tc)*Qv(1);
%     Qb=k(Tm)*Qv(end);
% 
%     HBot=(((Qt-Qbot)/((zv2(1))-(2*zb-zv2(1))))/(rho*Cp))+Ht(Tc)/(rho*Cp);
%     if HBot>=0
%         dzbdt=-(k(Tm)*HBot)/Qbot;
%     else
%         dzbdt=-(k(Tm)*HBot)/Qt;
%     end
% 
%     [H,z3]=FluxNumerical(Q,z2);
%     dTdz=CentralDifferenceNumerical(T,z);
%     H=K(273)*H+dzbdt.*[2:(nz-1)].*dTdz/nz+K(273)*Ht(T(2:end-1)); %Conduction eq.
%     dzmdt=Qb/(rho*L); %Movement of freezing front
% else
    %Convective system, defining convection parameters
    DT=Tm-Tc; %temperature drop across the convective layer
    Te=Tc+0.5*DT;
    G=Ht(Te); %heat production in the convective layer
    Ra=(alph*rho*g*DT*b^3)/(K(Te)*nm); % Rayleigh Number
    f=2*beta-3*gamma;
    Ti=Te+(D*DT*Ra^(-gamma)*((G*b^2)/(k(Te)*DT))^beta); %internal temperature
    Qt=((k(Tc)*DT)/b)*C*Ra^(xi)*((Ti-Tc)/DT)^zeta; %Heat flow thru upper TBL
    % Conductive lid, finite differential solution to conduction eq.
    T=Y(1:end-2);
    z=linspace(0,zb,nz);
    [Q,z2]=FluxNumerical(T,z);
    Qbot=k(Tc)*Q(end);
    HBot=(((Qt-Qbot)/(((2*zb-z2(end)))-z2(end)))/(rho*Cp))+Ht(Tc)/(rho*Cp);
    % Movement of interface btwn lid and convective layer
    if HBot>=0
        dzbdt=-(k(Tc)*HBot)/Qbot;
    else
        dzbdt=-(k(Tc)*HBot)/Qt;
    end

    [H,z3]=FluxNumerical(Q,z2);
    dTdz=CentralDifferenceNumerical(T,z);
    H=K(T(2:end-1)).*H+dzbdt.*[2:(nz-1)].*dTdz/nz+K(T(2:end-1)).*Ht(T(2:end-1));

    % coupling
    I=1+(Cp/L)*f*(Ti-Te); %Heat production in shell interior
    dzmdt=(Qt-G*b)/(rho*L*I); %freezing/melting
%end
dbdt=dzmdt-dzbdt; %change of thickness of the convective layer
%Results vector: track change in T for lid, movement of layer boundaries 
H=[0,H,0,dzbdt,dbdt]';


