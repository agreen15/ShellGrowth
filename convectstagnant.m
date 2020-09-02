clear all;

%% Setup constants
% Physical parameters
rho=917; %density
k=2.6; %thermal conductivity
K=1.47e-6; %thermal diffusivity
Cp=1930; %Specific heat
%Cp=k/(K*rho); %Specific heat
alpha=1.56e-4; %thermal expensivity
L=284; % heat of fusion
nm=5e13; %reference viscosity
E=60*1000; %Activation energy
R=8.3144621; %gas constant
N=1.8; %Stress exponent

%Set the following flags to 1 to produce figures in Green et al. (2020)
%If no flags are activated, only Figure 2 (reference model) will be
%produced.

%Repeat run with increase heat flag
IncHeat=0;
%Pressure-dependent melting flag (for large icy satellites)
pDep=0;
%Produces data and rough figure for Figure 4 in Green et al. (2020)
heatloop=0; 
%Produces raw data for Figure 5 in Green et al. 2020
transsect=0;

k=@(T)2.21-(0.012*(T-273.2));
K=@(T)k(T)./(rho*Cp);
K2=2.21/(rho*Cp);

%Orbital and physical parameters for small and large satellites
if pDep==0;
    g=1.315; %acceleration of gravity: Europa
    e=1.0e-10; %1/s average tidal strain rate: Europa
    w=3.2963e-06; %1/s orbital frequency: Europa
else
    g=1.428; %acceleration of gravity: Ganymede
    e=1e-12; %Avg tidal strain rate: Ganymede (?)
    w=1.6177e-06; %orbital frequency: Ganymede
end

HRad=(4.5e-12)*rho; %Volumetric radiogenic heating rate


u=3.3e9; %Pa shear modulus

Ts=100; %surface temperature
% melting temperature; can be made p-dependent
%p=(rho*g*zm)*1e-9;
Tm=273.2;%*(1-(p/0.395))^(1/9);
DTs=Tm-Ts; %temperature drop across ice shell
Te=(DTs/2)+Ts;
Q=(E/(R*Te^2))*DTs;

%define temperature in well mixed interior according to Deschamps
c1=1.43; %Deschamps Eq16
c2=-0.03;%Deschamps Eq16
B=E/(2*R*c1); %Deschamps Eq18
% C=c2*DTs; %Deschamps Eq18; variable C used later -> explicit substitution
Ti=B*(sqrt(1+(2/B)*(Tm-c2*DTs))-1); %Deschamps Eq18
%define temperature at the base of the conductive layer 
A=E/(R*Tm);
ni=nm*exp(A*((Tm/Ti)-1));
dnidTi=-nm*((A*Tm)/Ti^2)*exp(A*((Tm/Ti)-1));
DTv=-(ni/dnidTi);%*DTs;
DTe=2.24*DTv;
Tc=Ti-DTe;
DT=Tm-Tc; %Temperature drop across the convective layer


% convective parameters
C=0.3446;
D=1;
beta=0.75;
gamma=0.25;
xi=1/3;
zeta=4/3;

% setup heating rate function
Ht=DefineTidalHeat(rho,Tm,e,w,u,E,nm,R,HRad);

% numerical parameters
nz=100; %Discretization of the conductive profile
tmax=((2.0)*1e6)*60*60*24*365;%duration of run, Myrs
timesteps=1000;
dt=tmax/timesteps; %export time step 

%% initial conditions
zm0=400; %initial ice layer thickness
Qp=0;
ConductiveStart=1;
if ConductiveStart==1; %Initial conditions for conduction-only system
   
    ThetaC=(Tc-Ts)/DTs;
    zc=ceil(ThetaC*100);
    T1=linspace(Ts,Tm,nz);
    z1=linspace(0,zm0,nz);
    dTdz=(T1(zc)-T1(zc-1))/(z1(zc)-z1(zc-1));
    zb0=z1(zc-1)+(Tc-T1(zc-1))/dTdz;
else
    zb0=1.5270e+04; %initial conductive layer thickness
end
b0=zm0-zb0; %initial convective layer thickness
TinitStefan=1; %flag for Stefan initial profile
Benchmark=1; %Benchmark comparison between numerical solver and standard analytical stefan problem solution (Turcotte & Schubert, 2002)

if TinitStefan==1;
    if Benchmark==1;
        %Stefan profile
        LHS=(L*sqrt(pi))/(Cp*(Tm-Ts));
        RHS=@(L)(exp(-L.^2))./((L.*erf(L)));
        % Lambda=0
        % RHS=(exp(-Lambda^2))/((Lambda*erf(Lambda)))
        Lambda=fsolve(@(L)RHS(L)-LHS,2);
        %
        % while LHS~=RHS
        % Lambda=Lambda+.00001
        % RHS=(exp(-Lambda^2))/((Lambda*erf(Lambda)))
        % end

        %
        tall=[0:dt:tmax];%linspace(t0,tmax,dt);
        ym=2*sqrt(tall)*Lambda*sqrt(K2);

        eta=@(z,t)z/(2*sqrt(K2*t));
        Ta=@(z,t,zm)Ts+erf(eta(z,t))./erf(Lambda)*(Tm-Ts)./(z<=zm);
        size(tall);
    else
            %Stefan profile
        LHS=(L*sqrt(pi))/(Cp*(Tc-Ts));
        RHS=@(L)(exp(-L.^2))./((L.*erf(L)));
        % Lambda=0
        % RHS=(exp(-Lambda^2))/((Lambda*erf(Lambda)))
        Lambda=fsolve(@(L)RHS(L)-LHS,2);
        %
        % while LHS~=RHS
        % Lambda=Lambda+.00001
        % RHS=(exp(-Lambda^2))/((Lambda*erf(Lambda)))
        % end

        %
        tall=[0:dt:tmax];%linspace(t0,tmax,dt);
        ym=2*sqrt(tall)*Lambda*sqrt(K(273));

        eta=@(z,t)z/(2*sqrt(K(273)*t));
        Ta=@(z,t,zm)Ts+erf(eta(z,t))./erf(Lambda)*(Tc-Ts)./(z<=zm);
        size(tall);
    end

    
    for i=1:1:numel(tall)
        z=linspace(0,ym(i),nz);
        %eta=z./(2*sqrt(K*tall(i)));
        T(i,:)=Ta(z,tall(i),ym(i));%Ts+(erf(eta)./erf(Lambda)).*(Tm-Ts)./(z<=ym(i));
        Z(i,:)=z;
    end
    
    ts=dt;%1e12;%tall(is);
    y0=spline(tall,ym,ts);
    z0=linspace(0,y0,nz);
    ts=interp1(ym,tall,y0);
    %t=0;
    T0=Ta(z0,ts,y0);
    
else
    T0=[linspace(Ts,Tc,100)];
    ts=0;
end

%split initial conditions for pressure-dependent melting
if pDep==0;
    if Benchmark==1;
        Y0=[T0, zm0];
        [t,Y]=ode45(@(t,Y)StefanBenchmarkODE(Y,nz,K2,k,L,Cp,Ht,Tm,rho,Ts), [linspace(ts,tmax,1001)], Y0);
    else
    
    Y0=[T0, zb0, b0];

    [t,Y]=ode45(@(t,Y)StagnantLidODE(Y,nz,K,k,g,nm,L,Cp,Ht,Tm,rho,C,D,alpha,beta,gamma,xi,zeta,Ts),[ts tmax], Y0);%linspace(ts,tmax,timesteps)],Y0);
    end
else
    p0=(rho*g*zm0)*1e-9;
    Tm0=273.2*(1-(p0/0.395))^(1/9);
    Ti0=Ti;
    Y0=[T0, Tm0, p0, Ti0, zb0, b0];
    [t,Y]=ode45(@(t,Y)StagnantLidODEv2(Y,nz,K,k,g,nm,L,Cp,Ht,Tm,rho,C,B,D,A,alpha,beta,gamma,xi,zeta,Ts),[linspace(ts,tmax,1000)],Y0);
end
%% postprocessing

Te=Tc+0.5*DT; %internal temperature w/o internal heating
%Ti=Te+1.236*(((H*b)/k)^0.75)*(((K*nm)/(alpha*g))^0.25);
% Heat generated at Te
H=Ht(Te); 
% Hm=H/rho;
% Hs=(H*b0^2)/(k*DT);

ball=Y(:,end);
zball=Y(:,end-1);
zmall=zball+ball;
Raall=(alpha*rho*g*DT.*ball.^3)./(K(Tm)*nm);
Tiall=Te+(D*DT.*Raall.^(-gamma).*((H*ball.^2)/(k(Te)*DT)).^beta);
Qsall=((k(Ts)*(Tc-Ts))./zball);
Qcall=((k(Tc)*(Tm-Ts))./zmall);
Nuall=Qsall./Qcall;
if pDep==1;
    Tiall=Y(:,end-2);
    pall=Y(:,end-3);
    Tmall=Y(:,end-4);
    Tcall=Y(:,end-5);
end

ty=t/(365.24*24*3600)/1e6; % convert to Millions of years

figure(1)
clf;
%plot(ty,(Y(:,end)/1000)); %b
hold on
plot(ty,(Y(:,end-1)/1000),'r'); %zb
plot(ty,((Y(:,end)+Y(:,end-1))/1000),'g'); %zm
title('Ice Thickness','fontsize',20)
set(gca,'ydir','reverse');
legend('z_b','z_m')
xlabel('Time (Myr)','fontsize',16)
ylabel('Depth (km)','fontsize',16)

if pDep==1;
    figure(2)
    clf;
    hold on
    plot(ty,Tiall,'b')
    plot(ty,Tmall,'r')
    plot(ty,Tcall,'g')
    title('Interior Temperature and Melting Temperature','fontsize',20)
    legend('Ti','Tm','Tc')
    xlabel('Time (Myr)','fontsize',16)
    ylabel('Temperature (K)','fontsize',16)
    
    figure(3)
    clf;
    plot(ty,pall)
    title('Shell Pressure','fontsize',20)
    xlabel('Time (Myr)','fontsize',16)
    ylabel('Pressure (GPa)','fontsize',16)
else
    figure(2)
    clf;
    hold on
    plot(ty,Tiall,'b')
    %plot(ty,Tmall,'r')
    title('Interior Temperature','fontsize',20)
    xlabel('Time (Myr)','fontsize',16)
    ylabel('Temperature (K)','fontsize',16)
end
% %set(gca,'ydir','reverse');
% iti

figure(4)
clf;
plot(ty,Raall)
title('Rayleigh Number','fontsize',20)
xlabel('Time (Myr)','fontsize',16)
ylabel('Ra','fontsize',16)
%%
Shell=[zball, zmall];
figure(5)
clf;
RefModelFig(ty,Tiall,Shell/1000,Raall/1E5)
% subplot(1,3,1);
% hold on
% plot(ty,Tiall,'b')
% %plot(ty,Tmall,'r')
% title('Interior Temperature','fontsize',20)
% xlabel('Time (Myr)','fontsize',16)
% ylabel('Temperature (K)','fontsize',16)
% 
% subplot(1,3,2)
% hold on
% plot(ty,(Y(:,end-1)/1000),'r'); %zb
% plot(ty,((Y(:,end)+Y(:,end-1))/1000),'g'); %zm
% title('Ice Thickness','fontsize',20)
% set(gca,'ydir','reverse');
% legend('z_b','z_m')
% xlabel('Time (Myr)','fontsize',16)
% ylabel('Depth (km)','fontsize',16)
% 
% subplot(1,3,3)
% plot(ty,Raall)
% title('Rayleigh Number','fontsize',20)
% xlabel('Time (Myr)','fontsize',16)
% ylabel('Ra','fontsize',16)

% plot(ty,Nuall)
% title('Nusselt Number','fontsize',20)
% xlabel('Time (Myr)','fontsize',16)
% ylabel('Nu','fontsize',16)

if Benchmark==1;
    tall=tall/(365.24*24*3600)/1e6;
    RelErr=(ym-ball')./ym
    figure(7)
    clf;
    hold on
    plot(ty,ball/1000,'b')
    plot(tall,ym/1000,'r')
    title('Freezing Front Comparison','fontsize',20)
    xlabel('Time (Myr)','fontsize',16)
    ylabel('Layer Thickness (km)','fontsize',16)
    figure(9)
    clf;
    plot(ty,RelErr*100)
    title('Approximation Error','fontsize',20)
    xlabel('Time (Myr)','fontsize',16)
    ylabel('Percent Difference','fontsize',16)
    
end

%% Second Run: Increase in heat
if IncHeat==1;

    if pDep==0;

        e=3.0e-10; %1/s average tidal strain rate: Europa
        Qp=10;
        

    else

        e=6e-10; %Avg tidal strain rate: Ganymede (?)
    end

    Ht=DefineTidalHeat(rho,Tm,e,w,u,E,nm,R,HRad);
    tmax2=((.2)*1e6)*60*60*24*365;

    if pDep==0;

        Y02=Y(end,:);
        %Y02(end)=Y02(end)-4000; %"Lateral flow" thickness perturbation
        [t2,Y2]=ode45(@(t,Y)StagnantLidODE(Y,nz,K,k,g,nm,L,Cp,Ht,Tm,rho,C,D,alpha,beta,gamma,xi,zeta,Ts),[linspace(ts,tmax2,timesteps)],Y02);
    else
        Y02=Y(end,:);
        %Y02(end)=Y02(end)-4000; %"Lateral flow" thickness perturbation
        [t2,Y2]=ode45(@(t,Y)StagnantLidODEv2(Y,nz,K,k,g,nm,L,Cp,Ht,Tm,rho,C,B,D,A,alpha,beta,gamma,xi,zeta,Ts),[linspace(ts,tmax2,1000)],Y02);
        Y02=Y(end,:);
    end


    Te=Tc+0.5*DT; %internal temperature w/o internal heating
    %Ti=Te+1.236*(((H*b)/k)^0.75)*(((K*nm)/(alpha*g))^0.25);
    % Heat generated at Te
    H=Ht(Te); 
    % Hm=H/rho;
    % Hs=(H*b0^2)/(k*DT);

    ball2=Y2(:,end);
    zball2=Y2(:,end-1);
    zmall2=zball2+ball2;
    Raall2=(alpha*rho*g*DT.*ball2.^3)./(K*nm);
    Tiall2=Te+(D*DT.*Raall2.^(-gamma).*((H*ball2.^2)/(k(Te)*DT)).^beta);
    Qsall2=((k(Ts)*(Tc-Ts))./zball2);
    Qcall2=((k(Tc)*(Tm-Ts))./zmall2);
    Nuall2=Qsall2./Qcall2;
    if pDep==1;
        Tiall2=Y(:,end-2);
        pall2=Y(:,end-3);
        Tmall2=Y(:,end-4);
        Tcall2=Y(:,end-5);
    end

    ty2=t2/(365.24*24*3600)/1e6; % convert to Millions of years
    
    ty2k=ty2*1000;
    
    [Lm, Im]=min(zball2);
    tLm=ty2k(Im);
    
    [sm, im]=min(round(zmall2/1000,2));
    tsm=ty2k(im);
    
    Lm, tLm, sm, tsm

    figure(6)
    clf;
    % Create axes
    axes1 = axes('Parent',figure(6));
    hold(axes1,'on');

    % Create multiple lines using matrix input to plot
    hold on
    plot(ty2,(Y2(:,end-1)/1000),'r'); %zb
    plot(ty2,((Y2(:,end)+Y2(:,end-1))/1000),'b'); %zm

    % Create ylabel
    ylabel('Depth (km)','FontSize',16);

    % Create xlabel
    xlabel('Time (Myr)','FontSize',16);

    % Create title
    title('Ice Thickness','FontSize',20);

    axis(axes1,'ij');
    % Set the remaining axes properties
    set(axes1,'XGrid','on','XTick',...
        [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1],...
        'YGrid','on');
    % Create legend
    legend('z_b','z_m');
% %     clf;
    %plot(ty,(Y(:,end)/1000)); %b
%     hold on
%     plot(ty2,(Y2(:,end-1)/1000),'r'); %zb
%     plot(ty2,((Y2(:,end)+Y2(:,end-1))/1000),'g'); %zm
%     title('Ice Thickness','fontsize',20)
%     set(gca,'ydir','reverse');
%     legend('z_b','z_m')
%     xlabel('Time (Myr)','fontsize',16)
%     ylabel('Depth (km)','fontsize',16)

    if pDep==1;
        figure(7)
        clf;
        hold on
        plot(ty2,Tiall2,'b')
        plot(ty2,Tmall2,'r')
        plot(ty2,Tcall2,'g')
        title('Interior Temperature and Melting Temperature','fontsize',20)
        legend('Ti','Tm','Tc')
        xlabel('Time (Myr)','fontsize',16)
        ylabel('Temperature (K)','fontsize',16)

        figure(8)
        clf;
        plot(ty2,pall2)
        title('Shell Pressure','fontsize',20)
        xlabel('Time (Myr)','fontsize',16)
        ylabel('Pressure (GPa)','fontsize',16)
    else
        figure(7)
        clf;
        hold on
        plot(ty2,Tiall2,'b')
        %plot(ty,Tmall,'r')
        title('Interior Temperature','fontsize',20)
        xlabel('Time (Myr)','fontsize',16)
        ylabel('Temperature (K)','fontsize',16)
    end
    % %set(gca,'ydir','reverse');
    % iti

    figure(9)
    clf;
    plot(ty2,Raall2)
    title('Rayleigh Number','fontsize',20)
    xlabel('Time (Myr)','fontsize',16)
    ylabel('Ra','fontsize',16)

    figure(10)
    clf;
    plot(ty2,Nuall2)
    title('Nusselt Number','fontsize',20)
    xlabel('Time (Myr)','fontsize',16)
    ylabel('Nu','fontsize',16)
end

if heatloop==1;

    tmax2=((1.0)*1e6)*60*60*24*365;
    Y02=Y(end,:);
    eall=[1.5e-10, 2.0e-10, 2.5e-10, 3.0e-10];
    YALL=zeros(102,5);
    YALL(:,1)=Y02;
    TALL=zeros(timesteps,4);
    ZBALL=zeros(timesteps,4);
    ZMALL=zeros(timesteps,4);
    for i=1:4;
        e=eall(i);
        Ht=DefineTidalHeat(rho,Tm,e,w,u,E,nm,R,HRad);
        [ti,Yi]=ode45(@(t,Y)StagnantLidODE(Y,nz,K,k,g,nm,L,Cp,Ht,Tm,rho,C,D,alpha,beta,gamma,xi,zeta,Ts),[linspace(ts,tmax2,timesteps)],Y02);
        ti=ti/(365.24*24*3600)/1e6;
        zballi=Yi(:,end-1);
        balli=Yi(:,end);
        zmalli=zballi+balli;
        YALL(:,(i+1))=Yi(end,:);
        TALL(:,i)=ti*1000;
        ZBALL(:,i)=zballi/1000;
        ZMALL(:,i)=zmalli/1000;


    end
 %%   
    e=3e-10;
    Ht=DefineTidalHeat(rho,Tm,e,w,u,E,nm,R,HRad);
    YALL2=zeros(102,4);
    TALL2=zeros(timesteps,4);
    ZBALL2=zeros(timesteps,4);
    ZMALL2=zeros(timesteps,4);
    tmax3=((0.5)*1e6)*60*60*24*365;
    %Comment/uncomment line below to toggle viscosity test
    %YALL(end-1,:)=[YALL(end-1,1),YALL(end-1,1),YALL(end-1,1),YALL(end-1,1),YALL(end-1,1)];
    for i=1:4;
        Y0i=YALL(:,i);
        [ti,Yi]=ode45(@(t,Y)StagnantLidODE(Y,nz,K,k,g,nm,L,Cp,Ht,Tm,rho,C,D,alpha,beta,gamma,xi,zeta,Ts),[linspace(ts,tmax3,timesteps)],Y0i);
        ti=ti/(365.24*24*3600)/1e6;
        zballi=Yi(:,end-1);
        balli=Yi(:,end);
        zmalli=zballi+balli;
        YALL2(:,i)=Yi(end,:);
        TALL2(:,i)=ti*1000;
        ZBALL2(:,i)=zballi/1000;
        ZMALL2(:,i)=zmalli/1000;


    end
    %%
    figure(11)
    clf;
    subplot(2,2,1);
  
    % Create multiple lines using matrix input to plot
    plot1 = plot(TALL,ZBALL,'LineWidth',2);
    set(plot1(1),'LineStyle','-.','DisplayName','');
    set(plot1(2),'LineStyle',':','DisplayName','');
    set(plot1(3),'LineStyle','--','DisplayName','');
    set(plot1(4),'DisplayName','');
    % Create ylabel
    ylabel('Depth (km)','FontSize',14);

    % Create xlabel
    xlabel('Time (Kyrs)','FontSize',14);

    % Create title
    title('Lithospheric Repsonse to Shell Thinning','FontSize',16);

    % Uncomment the following line to preserve the X-limits of the axes
    xlim([0 800]);
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([0 3.500]);
%     % Create legend
%     legend1 = legend(axes1,'show');
%     set(legend1,...
%         'Position',[0.742486832538371 0.727535735888702 0.129156772591475 0.0795819935691319]);
%     
    subplot(2,2,2);
    % Create multiple lines using matrix input to plot
    plot2 = plot(TALL,ZMALL,'LineWidth',2);
    set(plot2(1),'LineStyle','-.','DisplayName','');
    set(plot2(2),'LineStyle',':','DisplayName','');
    set(plot2(3),'LineStyle','--','DisplayName','');
    set(plot2(4),'DisplayName','');
    % Create ylabel
    ylabel('Depth (km)','FontSize',14);

    % Create xlabel
    xlabel('Time (Kyrs)','FontSize',14);

    % Create title
    title('Timescale of Shell Thinning','FontSize',16);

    % Uncomment the following line to preserve the X-limits of the axes
    xlim([0 800]);
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([0 30.000]);
%     % Create legend
%     legend1 = legend('show');
%     set(legend1,...
%         'Position',[0.742486832538371 0.727535735888702 0.129156772591475 0.0795819935691319]);

    subplot(2,2,3);
    % Create multiple lines using matrix input to plot
    plot3 = plot(TALL2,ZBALL2,'LineWidth',2);
    set(plot3(1),'LineStyle','-.','DisplayName','');
    set(plot3(2),'LineStyle',':','DisplayName','');
    set(plot3(3),'LineStyle','--','DisplayName','');
    set(plot3(4),'DisplayName','');
    % Create ylabel
    ylabel('Depth (km)','FontSize',14);

    % Create xlabel
    xlabel('Time (Kyrs)','FontSize',14);

    % Create title
    title('Lithospheric Repsonse to Shell Thinning','FontSize',16);

    % Uncomment the following line to preserve the X-limits of the axes
    xlim([0 400]);
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([0 3.500]);
    
    legend(' ',' ',' ',' ')

%     % Create legend
%     legend1 = legend('show');
%     set(legend1,...
%         'Position',[0.742486832538371 0.727535735888702 0.129156772591475 0.0795819935691319]);
    
    figure14 = subplot(2,2,4);


    % Create multiple lines using matrix input to plot
    plot4 = plot(TALL2,ZMALL2,'LineWidth',2);
    set(plot4(1),'LineStyle','-.','DisplayName','');
    set(plot4(2),'LineStyle',':','DisplayName','');
    set(plot4(3),'LineStyle','--','DisplayName','');
    set(plot4(4),'DisplayName','');
    % Create ylabel
    ylabel('Depth (km)','FontSize',14);

    % Create xlabel
    xlabel('Time (Kyrs)','FontSize',14);

    % Create title
    title('Timescale of Shell Thinning','FontSize',16);

    % Uncomment the following line to preserve the X-limits of the axes
    xlim([0 400]);
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([0 30.000]);

%     % Create legend
%     legend1 = legend('show');
%     set(legend1,...
%         'Position',[0.742486832538371 0.727535735888702 0.129156772591475 0.0795819935691319]);

end
if transsect == 1
    etrans0=[1.1e-10;1.15e-10;1.3e-10;1.62e-10;2.08e-10;2.42e-10;2.59e-10];
    etrans45=[1.59e-10;1.7e-10;1.9e-10;2.14e-10;2.35e-10;2.5e-10;2.59e-10];
    surftemp=[96;95;94;92;90;72;46];
    for i=1:7
        e0=etrans0(i);
        e45=etrans45(i);
        Ts=surftemp(i);

        DTs=Tm-Ts; %temperature drop across ice shell
        Te=(DTs/2)+Ts;
        Q=(E/(R*Te^2))*DTs;

        %define temperature in well mixed interior according to Deschamps
        c1=1.43; %Deschamps Eq16
        c2=-0.03;%Deschamps Eq16
        B=E/(2*R*c1); %Deschamps Eq18
        % C=c2*DTs; %Deschamps Eq18; variable C used later -> explicit substitution
        Ti=B*(sqrt(1+(2/B)*(Tm-c2*DTs))-1); %Deschamps Eq18
        %define temperature at the base of the conductive layer 
        A=E/(R*Tm);
        ni=nm*exp(A*((Tm/Ti)-1));
        dnidTi=-nm*((A*Tm)/Ti^2)*exp(A*((Tm/Ti)-1));
        DTv=-(ni/dnidTi);%*DTs;
        DTe=2.24*DTv;
        Tc=Ti-DTe;
        DT=Tm-Tc; %Temperature drop across the convective layer


        % convective parameters
        C=0.3446;
        D=1;
        beta=0.75;
        gamma=0.25;
        xi=1/3;
        zeta=4/3;

        % setup heating rate function
        Ht0=DefineTidalHeat(rho,Tm,e0,w,u,E,nm,R,HRad);
        Ht45=DefineTidalHeat(rho,Tm,e45,w,u,E,nm,R,HRad);

        % numerical parameters
        nz=100; %Discretization of the conductive profile
        tmax=((1.5)*1e6)*60*60*24*365;%duration of run, Myrs
        timesteps=1000;
        dt=tmax/timesteps; %export time step 

        %initial conditions
        zm0=400; %initial ice layer thickness
        Qp=0;
        ConductiveStart=1;
        if ConductiveStart==1; %Initial conditions for conduction-only system

            ThetaC=(Tc-Ts)/DTs;
            zc=ceil(ThetaC*100);
            T1=linspace(Ts,Tm,nz);
            z1=linspace(0,zm0,nz);
            dTdz=(T1(zc)-T1(zc-1))/(z1(zc)-z1(zc-1));
            zb0=z1(zc-1)+(Tc-T1(zc-1))/dTdz;
        else
            zb0=1.5270e+04; %initial conductive layer thickness
        end
        b0=zm0-zb0; %initial convective layer thickness
        TinitStefan=1; %flag for Stefan initial profile

        if TinitStefan==1;
            %Stefan profile
            LHS=(L*sqrt(pi))/(Cp*(Tc-Ts));
            RHS=@(L)(exp(-L.^2))./((L.*erf(L)));
            % Lambda=0
            % RHS=(exp(-Lambda^2))/((Lambda*erf(Lambda)))
            Lambda=fsolve(@(L)RHS(L)-LHS,2);
            %
            % while LHS~=RHS
            % Lambda=Lambda+.00001
            % RHS=(exp(-Lambda^2))/((Lambda*erf(Lambda)))
            % end

            %
            tall=[0:dt:tmax];%linspace(t0,tmax,dt);
            ym=2*sqrt(tall)*Lambda*sqrt(K(Tm));

            eta=@(z,t)z/(2*sqrt(K(Tm)*t));
            Ta=@(z,t,zm)Ts+erf(eta(z,t))./erf(Lambda)*(Tc-Ts)./(z<=zm);
            size(tall);

            for j=1:1:numel(tall)
                z=linspace(0,ym(j),nz);
                %eta=z./(2*sqrt(K*tall(i)));
                T(j,:)=Ta(z,tall(j),ym(j));%Ts+(erf(eta)./erf(Lambda)).*(Tm-Ts)./(z<=ym(i));
                Z(j,:)=z;
            end

            ts=dt;%1e12;%tall(is);
            y0=spline(tall,ym,ts);
            z0=linspace(0,y0,nz);
            ts=interp1(ym,tall,y0);
            %t=0;
            T0=Ta(z0,ts,y0);

        else
            T0=[linspace(Ts,Tc,100)];
            ts=0;
        end
        Y0=[T0, zb0, b0];

        [to,Yo]=ode45(@(t,Y)StagnantLidODE(Y,nz,K,k,g,nm,L,Cp,Ht0,Tm,rho,C,D,alpha,beta,gamma,xi,zeta,Ts),[ts tmax], Y0);
        [t45,Y45]=ode45(@(t,Y)StagnantLidODE(Y,nz,K,k,g,nm,L,Cp,Ht45,Tm,rho,C,D,alpha,beta,gamma,xi,zeta,Ts),[ts tmax], Y0);

        zball0=Yo(end,end-1);
        ball0=Yo(end,end);
        zmall0=zball0+ball0;
        zball45=Y45(end,end-1);
        ball45=Y45(end,end);
        zmall45=zball45+ball45;
        %YALL0trans(:,i)=Yo(end,:);
        %TALL0trans(:,i)=to*1000;
        ZBALL0trans(i)=zball0/1000;
        ZMALL0trans(i)=zmall0/1000;
        %YALL45trans(:,i)=Y45(end,:);
        %TALL45trans(:,i)=t45*1000;
        ZBALL45trans(i)=zball45/1000;
        ZMALL45trans(i)=zmall45/1000;
    end
    ZBALL0trans
    ZMALL0trans
    ZBALL45trans
    ZMALL45trans
end