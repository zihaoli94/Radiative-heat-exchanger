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

% Parameters
[h1,h2,R,A1,A2,Ap2,q1,q2,q3,q4,R1,R2,R12,R3,R4,R34,Tf]=RHX_Parameters(x,T1,r1,r2,S,k2,N,L,H,e1,e2,g,v);


% III. OUTPUT PART


QR(i,1)=h1*A1*(T1-T2); % Overall heat transfer rate (W).

QR(i,2)=h2*A2*(T3-Tf); % Heat tranfer by natural heat convection (unforced-air) (W).

QR(i,3)=(q1-q3)/(R1+R12+R2)+(q2-q4)/(R3+R34+R4); % Thermal radiation (W).


i=i+1;

end