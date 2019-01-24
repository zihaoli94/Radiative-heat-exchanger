function F=RadiativeHX_RectangularFins(x,T1,Tair,h3,r1,r2,S,k2,N,L,H,e1,e2,g,v)

% THIS FUNCTION IS USED FOR RADIATIVE HX DESIGN 
% (TUBE SURFACE & ABSORBER PLATE SURFACE ARE BOTH PIANTED)


% I. Parameters Applied in the Heat Balance Equations

[h1,h2,R,A1,A2,Ap2,q1,q2,q3,q4,R1,R2,R12,R3,R4,R34,Tf]=RHX_Parameters(x,T1,r1,r2,S,k2,N,L,H,e1,e2,g,v);


% II. Heat balances

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

end