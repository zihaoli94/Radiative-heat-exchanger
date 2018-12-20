function F = ConventionalHX(x,T1,h2,Tair,r1,r2,k2,L,N,e,sigma,f,v)

% THIS FUNCTION IS USED FOR CONVENTIONAL HX DESIGN 

% T1: molten-salt temperature.
% T2: tube inner-wall temperature.
% T3: tube outer-wall surface temperature.

T2=x(1);
T3=x(2);


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


% 3. Surface areas


% Surface area of radiative heat exchanger tubes:

A1=2*pi*r1*L*N; % A1: internal surface of heat exchanger tubes.

A2=2*pi*r2*L*N; % A2: external surface of heat exchanger tubes.


% 4. Thermal radiation from tube surface to chimney walls
% Two-surface enclosure radiation assumed.


Eb1=sigma*x(2)^4; % Eb1: black body radiation of tube.

Eb2=sigma*Tair^4; % Eb2: black body radiation of chimney walls.

R1=(1-e)/(A2*e); %R1: surface resistance of tube to radiation.

R2=0; % R2: surface resistance of chimney walls to radiatino (A>>> then R2=0).

R12=1/(A2*f); % R12: space resistance between tubes and chimney walls.

J2=Eb2; % J2: radiosity of chimney (assumed as black body, resistance=0 since large surface area).

% (Eb1-J1)/R1=(J1-J2)/R12, where J2=Eb2, then: 

J1=(Eb1*R12+J2*R1)/(R1+R12); % J1: radiosity of tube (net radiation balance).

% The thermal radiation rate: Qrad = (Eb1-J1)/R1.

% 5. Heat balances


F(1)=h1*(T1-x(1))*A1-((x(1)-x(2))/R)*L*N; 
% Heat balance (1): heat convection through molten salt fluid bulk = heat
% conduction through tube wall.

F(2)=h1*(T1-x(1))*A1-h2*(x(2)-Tair)*A2-(Eb1-J1)/R1; 
% Heat balance (2): heat convection through molten salt fluid bulk =
% convection + radiation.


end
