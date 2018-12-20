clc
clear all

%THIS SCRIPT IS USED FOR RADIATIVE HEAT EXCHANGER 
% (TUBE SURFACE & ABSORBER PLATE SURFACE ARE BOTH PIANTED)


% I. INPUT PART


Tair=333; % Tair: air temp (K).

h3=35; % h3: external air-side heat transfer coefficient (W/(m^2K)).

r1=0.00665; % r1: internal radius of heat exchanger tube (m).

r2=0.00795; % r2: external radius (m). 

S=0.02385; % S: distance between two tubes centres (m). 

k2=14.2; % k2: thermal conductivity of tube wall(SS310 material propertity)(W/(mK)).

N=100; % N: number of heat exchanger tube.

L=1.93; % L: lenght of tube/plate (m).

H=0.60; % H: height of absorber plate (m). 

e1=0.96; % e1: emissivity at tube surface at 0~0.25um wavelength.

e2=0.96; % e2: emissivity at 0.25um~inf wavelength.

g=9.81; % g: gravity (m/s^2). 

v=0.2; % v: flow velocity of the molten-salt (m/s).

i=1;


% II. FSOLVE PART


for T1=723:10:1373 % T1: independent variable (K).
    
x0=[350,300,250];

x=fsolve(@(x) RadiativeHX_RectangularFins(x,T1,Tair,h3,r1,r2,S,k2,N,L,H,e1,e2,g,v),x0);

% T1: molten-salt temperature.
% T2: tube inner-wall temperature.
% T3: tube outer-wall surface temperature.
% T4: absorber plate temperature.

rx(1)=real(x(1));
rx(2)=real(x(2));
rx(3)=real(x(3));

T2=rx(1);
T3=rx(2);
T4=rx(3);


% 1. Following is molten salt (FLiNaK) physical property correlation


Cp=1884; % Cp: specific heat capacity.

k=0.43+(5e-4)*T1; % k: thermal conductivity.

rho=2730-0.73*T1; % rho: density.

mu=(4e-5)*exp(4170/T1); % mu: dynamic viscosity.

D1=2*r1; % D1: internal diameter.

Re=(rho*v*D1)/mu; % Re: Reynold's number.

Pr=(mu*Cp)/k; % Pr: Prandtl number.

ff=3.03e-12*Re^3-3.67e-8*Re^2+1.46e-4*Re-0.151; % ff: friction factor for fully developed flow in 2300<Re<4000.

if Re<2300 % Laminar flow region.
    Nu=3.65+(0.0668*Re*Pr*(D1/L))/(1+0.04*Re*Pr*(D1/L)); % Nu: Nusselt number.
    
elseif Re>4000 % Turbulent flow region.
    Nu=0.027*(Re^0.8)*(Pr^0.3);   
    
else % Transition region.
    Nu=((ff/8)*(Re-1000)*Pr)/(1+12.7*(ff/8)^0.5*(Pr^(2/3)-1));   
end

h1=(Nu*k)/D1; % h1: heat transfer coefficient of molten salt flow.


% 2. Wall heat conduction


R=log(r2/r1)/(2*pi*k2); % R: heat conduction resistance of tube wall.


% 3. Following is air physical properties correlation (Unforced air).


Tf=(rx(2)+rx(3))/2; % Tf: temperature of air film on tube surface.

mu2=((1.4592e-6)*Tf^(3/2))/(109.10+Tf); % mu2: dynamic viscosity of air film.

Cp2=1030.5-0.19975*Tf+(3.9734e-4)*Tf^2; % Cp2: specific heat capacity.

k3=((2.334e-3)*Tf^(3/2))/(164.54+Tf); % k: thermal conductivity.

rho2=(351.99/Tf)+(344.84/(Tf^2)); %rho2: density of air film.

alpha=k3/(rho2*Cp2); % alpha: thermal diffusivity.

beta=1/Tf; % beta: thermal expansion coefficient of air film.

nu=mu2/rho2; % nu: kinetic viscosity. 

D2=2*r2; % D2: external diameter.

% Wilks’s correlation:

Pr2=nu/alpha;

Gr=(g*beta*(rx(2)-rx(3))*D2^3)/(nu)^2; % Gr: the Grashof number.
 
Ra=Pr2*Gr;

NuD=0.579*(Ra/(1+(0.442/Pr2)^(9/16))^(16/9))^(1/4); % NuD: Nusselt number of unforced air.

h2=(NuD*k3)/D2; % h2: heat transfer coefficient of unforced air.


% 4. View factor of two equal external cylinders


h=S/r2; % h: distance between centres of two tubes (S)/external radius (r2).

F12=((h^2-4)^0.5-h+2*asin(2/h))/(2*pi); % F12: view factor of two equal exteranl diameter tubes.

f=1-2*F12; % f: view factor bewtween heat exchanger tube and radation absorber plate.


% 5. Thermal radiation from tube surface that reaches to radiation absorber plate surface
% Note: these are paticularly used for surface coated with selective
% coating, which has different emissivity values at different range of
% wavelength (0~0.25um; 0.25um~infinity).

fun = @(lambda,T3) 3.742e8./(lambda.^5.*(exp(1.439e4./(lambda.*T3))-1));

q1 = integral(@(lambda)fun(lambda,T3),0,2.5);

q2 = integral(@(lambda)fun(lambda,T3),2.5,inf);

fun = @(lambda,T4) 3.742e8./(lambda.^5.*(exp(1.439e4./(lambda.*T4))-1));

q3 = integral(@(lambda)fun(lambda,T4),0,2.5);

q4 = integral(@(lambda)fun(lambda,T4),2.5,inf);


% 6. Surface areas


% Surface area of radiative heat exchanger tubes:

A1=2*pi*r1*L*N; % A1: internal surface of heat exchanger tubes.

A2=2*pi*r2*L*N; % A2: external surface of heat exchanger tubes.

% Surface area of radiation absorber plates (one assembly consists 25 tubes and two absorber plates):

Ap1=H*L*8; % Ap1: total surface area of tube-side absorber plates (8 plates).

a=100; % a: number of fins on one side of absorber plate.

b=0.012; % b: fin height (m).

Ap2=(H*L+2*a*b*H)*8; % Ap2: total surface area of air-side absorber plates (8 plates).


% 7. Radiation resistances


% At wavelength 0~0.25um:

R1=(1-e1)/(A2*e1); % R1: tube surface resistance to radiation.

R2=(1-e1)/(Ap1*e1); % R2: absorber plate surface resistance to radiation. 

R12=1/(A2*f); % R12: space resistance to radiation.

% At wavelength 0.25um~infinity:

R3=(1-e2)/(A2*e2); % R3: tube surface resistance to radiation.

R4=(1-e2)/(Ap1*e2); % R4: absorber plate surface resistance to radiation. 

R34=1/(A2*f); % R34: space resistance to radiation.


% 8. Heat balances


F(1)=h1*A1*(T1-x(1))-((x(1)-x(2))/R)*L*N; 
% Heat balance (1): heat convection through molten salt fluid = heat
% conduction through tube wall.

F(2)=((x(1)-x(2))/R)*L*N-h3*Ap2*(x(3)-Tair);
% Heat balance (2): heat conduction through tube wall = heat removed by
% external forced-air on absorber plate surface.

F(3)=h2*A2*(x(2)-Tf)+(q1-q3)/(R1+R12+R2)+(q2-q4)/(R3+R34+R4)-h3*Ap2*(x(3)-Tair);
% Heat balance (3): 
% convection heat from tube surface to absorber plate (by unforced air) +
% radiation (at wavelength range 0~0.25um)+ radiation (at wavelength
% 0.25um~inf) = heat removed by external forced-air.


% III. OUTPUT PART


QR(i,1)=h1*A1*(T1-rx(1)) % Overall heat transfer rate (W).

QR(i,2)=h2*A2*(x(2)-Tf) % Heat tranfer by natural heat convection (unforced-air) (W).

QR(i,3)=(q1-q3)/(R1+R12+R2)+(q2-q4)/(R3+R34+R4) % Thermal radiation (W).


i=i+1;


end