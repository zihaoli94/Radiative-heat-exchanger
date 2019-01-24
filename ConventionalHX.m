function F = ConventionalHX(x,T1,h2,Tair,r1,r2,k2,L,N,e,sigma,f,v)

% THIS FUNCTION IS USED FOR CONVENTIONAL HX DESIGN 


% I. Parameters applied in the heat balance equations

[h1,R,A1,A2,R1,Eb1,J1]=CHX_Parameters(x,T1,Tair,r1,r2,k2,L,N,e,sigma,f,v);


% II. Heat balances

F(1)=h1*(T1-x(1))*A1-((x(1)-x(2))/R)*L*N; 
% Heat balance (1): heat convection through molten salt fluid bulk = heat
% conduction through tube wall.

F(2)=h1*(T1-x(1))*A1-h2*(x(2)-Tair)*A2-(Eb1-J1)/R1; 
% Heat balance (2): heat convection through molten salt fluid bulk =
% convection + radiation.


end
