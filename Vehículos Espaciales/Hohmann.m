%% TRANSFERENCIA DE HOHMANN
%% Datos fijos:
R_T=6378; %km
mu=3.986*10^5; %km^3/s^2
J2=1.0827*10^-3;
%% Ecuaciones útiles:
%Haz copy, paste y sustituye lo que corresponda.
%Velocidad circunferencia: sqrt(mu/r)
%Velocidad elipse: sqrt(mu*( (2/r)-(1/a) ) )
%Periodo: T=2*pi*sqrt(a^3/mu)
%Semieje mayor (si te dan T y tienes que despejar): a=(mu*(T/(2*pi))^2)^(1/3)
%% Datos enunciado:
%Introdúcelos. REVISA MUY BIEN LOS NOMBRES Y TODO, NO LA VAYAMOS A LIAR.
a=;
ri=;
rf=;
a_trans=(ri+rf)/2;
Vi_a=;
Vtrans_a=sqrt(mu*( (2/ri)-(1/a_trans) ) );
Vf_b=;
Vtrans_b=sqrt(mu*( (2/rf)-(1/a_trans) ) );
%% Resultados buscados
%Si te dan alguna de estas cosas, vas a tener que despejar, amiga, casi que
%mejor hazlo en papel...
DeltaV_a=Vtrans_a-Vi_a;
DeltaV_b=Vf_b-Vtrans_b;
DeltaV_total=DeltaV_a+DeltaV_b;
T_trans=2*pi*sqrt(a_trans^3/mu);
t_trans=T_trans/2;