clc
clear all

%THIS SCRIPT IS USED FOR CONVENTIONAL HEAT EXCHANGER 


% I. INPUT PART


Tair=333; % Tair: air temp (K).

h2=35; % h2: external air-side heat transfer coefficient (W/(m^2K)).

r1=0.00665; % r1: internal radius of heat exchanger tube (m).

r2=0.00795; % r2: external radius (m). 

k2=14.2; % k2: thermal conductivity of tube wall(SS310 material propertity)(W/(mK)).

N=100; % N: number of heat exchanger tube.

L=1.93; % L: lenght of tube/plate (m).

e=0.39; % e: emmisivity of tube(ss310).

sigma=5.67e-8; % sigma: constant.

f=0.1; %f: the view factor of heat exchanger tubes to outside (assumed).

v=0.2; % v: flow velocity of the molten-salt (m/s).

i=1;


% II. FSOLVE PART


for T1=723:10:1373 % T1: independent variable (K).
    
x0=[300,290];

x=fsolve(@(x) ConventionalHX(x,T1,h2,Tair,r1,r2,k2,L,N,e,sigma,f,v),x0);


% T1: molten-salt temperature.
% T2: tube inner-wall temperature.
% T3: tube outer-wall surface temperature.


T2=x(1);
T3=x(2);


% Parametes applied

[h1,R,A1,A2,R1,Eb1,J1]=CHX_Parameters(x,T1,Tair,r1,r2,k2,L,N,e,sigma,f,v);



% III. OUTPUT PART


QR(i,1)=h1*(T1-T2)*A1; % Overall heat transfer rate (W).

QR(i,2)=h2*(T3-Tair)*A2; % Heat tranfer by heat convection (forced-air) (W).

QR(i,3)=(Eb1-J1)/R1; % Thermal radiation (W).

i=i+1;

end   
