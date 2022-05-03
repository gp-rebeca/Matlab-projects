%% EL TÍPICO DE TEORÍA QUE TE PIDE UNA COMPONENTE DE h O DE OTRO VECTOR RELACIONADO CON SU TRIEDRO
%% Datos que suelen dar:
mu=3.986*10^5; %km^3/s
nu=0;
%Introduce los valores de tu enunciado:
i=;
RAAN=;
w=;
w_deg=w*180/pi;
e=;
rp=;
%% Cálculos
h=sqrt(rp*(1+e)*mu);
hz=h*cos(i);
raiz_denominador=sqrt(h^2-hz^2); %sqrt(hx^2+hy^2)=sqrt(h^2-hz^2)
%¡¡MUCHO CUIDADO!! EL SIGNO QUE TE ESCUPE ESTO PUEDE ESTAR BIEN O MAL. PARA
%SABERLO HAY QUE DIBUJAR EL VECTOR n Y PONER DÓNDE QUEDARÍA EL VECTOR e
%SEPARADO UN ÁNGULO w. EL VECTOR h FORMA TRIEDRO CON ESTOS DOS, ASÍ VES EN
%QUÉ CUADRANTE QUEDA.
hx=sin(RAAN)*raiz_denominador;
hy=cos(RAAN)*raiz_denominador;
h_vector=[hx,hy,hz];
n=[cos(RAAN),sin(RAAN),0];
% En otro de los problemas pedían una componente del radiovector del
% perigeo en sistema de referencia geocéntrico ecuatorial:
% R_inercial=[cos(RAAN)*cos(w)-sin(RAAN)*sin(w)*cos(i), -cos(RAAN)*sin(w)-sin(RAAN)*cos(w)*cos(i), sin(RAAN)*sin(i);...
%             sin(RAAN)*cos(w)+cos(RAAN)*sin(w)*cos(i), -sin(RAAN)*sin(w)+cos(RAAN)*cos(w)*cos(i), -cos(RAAN)*sin(i);...
%             sin(w)*sin(i), cos(w)*sin(i), cos(i)];
% R_orbital=[cos(nu), sin(nu), 0; -sin(nu), cos(nu), 0; 0,0,1];
R1= [cos(RAAN), sin(RAAN), 0;-sin(RAAN), cos(RAAN), 0; 0, 0, 1 ]; 
R2= [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R3=[cos(nu+w), sin(nu+w), 0; -sin(nu+w), cos(nu+w),0; 0, 0, 1];
R=R3*R2*R1;
rp_vector_dato=[rp,0,0];
rp_vector_respuesta=rp_vector_dato*R; %¿Sería al revés?