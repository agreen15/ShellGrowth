function H = StagnantLidODEv2(Y,nz,K,k,g,nm,L,Cp,Ht,Tm0,rho,C,B,D,A,alph,beta,gamma,xi,zeta,Ts)
%Non-dimensional Convection
Y=Y';
%zm=Y(end);zb=Y(end-1);b=zm-zb;
b=Y(end);zb=Y(end-1);zm=b+zb;
Ti=Y(end-2);
p=Y(end-3);
Tm=Y(end-4);
Tc=Y(end-5);
% Ts=Tc;
Qs=((k*(Tc-Ts))/zb);
Qc=((k*(Tm-Ts))/zm);
Nu=Qs/Qc; %Nusselt Number, convects if >=1
if Nu<1;
    %Conduction only solver
    %T-profile for upper shell
    T=Y(1:end-5);
    z=linspace(0,zb,nz);
    %T-profile for lower shell
    Tv=linspace(Tc,Tm,nz);
    zv=linspace(zb,zm,nz);
    Tcn=Tc;
    
    [Q,z2]=FluxNumerical(T,z);
    Qbot=k*Q(end);
    
    
    [Qv,zv2]=FluxNumerical(Tv,zv);
    Qt=k*Qv(1);
    Qb=k*Qv(end);

    HBot=(((Qt-Qbot)/((zv2(1))-(2*zb-zv2(1))))/(rho*Cp))+Ht(Tc)/(rho*Cp);
    if HBot>=0
        dzbdt=-(k*HBot)/Qbot;
    else
        dzbdt=-(k*HBot)/Qt;
    end

    [H,z3]=FluxNumerical(Q,z2);
    dTdz=CentralDifferenceNumerical(T,z);
    H=K*H+dzbdt.*[2:(nz-1)].*dTdz/nz+K*Ht(T(2:end-1));
    dzmdt=Qb/(rho*L);
else
    %Convective system
    DT=Tm-Tc; %temperature drop accros the convective layer
    Te=Tc+0.5*DT;
    G=Ht(Te); %heat production in the convective layer
    Ra=(alph*rho*g*DT*b^3)/(K*nm); % Rayleigh Number
    f=2*beta-3*gamma;
    Ti2=Te+(D*DT*Ra^(-gamma)*((G*b^2)/(k*DT))^beta); %internal temperature
    Qt=((k*DT)/b)*C*Ra^(xi)*((Ti2-Tc)/DT)^zeta;
    % conductive lid
    T=Y(1:end-5);
    z=linspace(0,zb,nz);
    [Q,z2]=FluxNumerical(T,z);
    Qbot=k*Q(end);
    HBot=(((Qt-Qbot)/(((2*zb-z2(end)))-z2(end)))/(rho*Cp))+Ht(Tc)/(rho*Cp);
    if HBot>=0
        dzbdt=-(k*HBot)/Qbot;
    else
        dzbdt=-(k*HBot)/Qt;
    end

    [H,z3]=FluxNumerical(Q,z2);
    dTdz=CentralDifferenceNumerical(T,z);
    H=K*H+dzbdt.*[2:(nz-1)].*dTdz/nz+K*Ht(T(2:end-1));

    % coupling
    I=1+(Cp/L)*f*(Ti2-Te);
    dzmdt=(Qt-G*b)/(rho*L*I); %freezing/melting
end

% dp=(rho*g*dzmdt)*1e-9;
% dTm=-((0.281294*Tm0)/((1-2.53165*p)^(8/9)))*dp;

if Tm>=251.15
    dp=(rho*g*dzmdt)*1e-9;
    dTm=-((0.281294*Tm0)/((1-2.53165*p)^(8/9)))*dp;
else
    dp=0;
    dTm=0;
    dzmdt=0;
end
Tmn=273.2*(1-((p+dp)/0.395))^(1/9);
%dTm=Tmn-Tm;
%dTc=Tcn-Tc;
dbdt=dzmdt-dzbdt; %change of thickness of the convective layer
% dTc=0;
if Nu<1;
    dTi=0;
    dTc=0;
else
    dTi=(1.03/(sqrt((B-0.06*Ts+2.06*Tm)/B)))*dTm;
    dTc=(1-(2*2.24)*(Ti/(A*Tm)))*dTi;
end
H=[0,H,dTc,dTm,dp,dTi,dzbdt,dbdt]';


