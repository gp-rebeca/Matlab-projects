%Matlab Probelma estructura continua PA-MEF2019
%introduciendo los nodos manualmente en las líneas de código.
%% DEFINICIÓN MATERIAL. Valores unidad cuando no se conocen.
E=10; A=100; I=5; EA=E*A; EI=E*I; L=1; q0=1000;
%% DEFINICIÓN GEOMÉTRICA DEL MODELO
%grados de libertad por nodo
gdln=3;
%introduce las coordenadas de los nodos-articulaciones de la estructura.
%Cada fila corresponde a un elemento, coordenada en x en primera columna y
%coordenada y en segunda columna
coord_nodos=[0 0; 0 L/2; 0 L; L 1.5*L; 2*L 2*L; 2*L L; 2*L 0];
%coordenadas de los nodos
xx=coord_nodos(:,1); yy=coord_nodos(:,2);
nnodos=size(xx,1); %número de nodos
%alternativa:nnodos==length(coord_nodos(:,1))
nelementos=size(coord_nodos,1)-1;%número de elementos
GDL=nnodos*gdln; %grados de libertad de toda la estructura
%matriz de conectividad entre nodos: cada fila contiene el elemento que los
%une, igual número de filas que elementos...size(conectaelementos,1)
conectaelementos=zeros(nelementos,2);
for c=1:nelementos
    conectae=[c c+1]
    conectaelementos(c,:)=conectae;
end
%% FUERZAS NODALES EQUIVALENTES
q=zeros(GDL,1); %vector inicial vacío 21 filas, 1 columna
q1=q0/2
q2=q0
long=sqrt(5)*L/2
alfa=atan(0.5)
F1=3*q1*long/20
M1=q1*long^2/30
F2=7*q1*long/20
M2=q1*long^2/20
F3=(q1*long/2)+(3*(q2-q1)*long)/20
M3=(q1*long^2/12)+(q2-q1)*long^2/30
F4=(q1+q2)*long/2-F3
M4=((q2-q1)*long^2/6)+(q1*long^2/2)-F3*long+M3 %corrupción
q(7,1)=F1*sin(alfa)
q(8,1)=-F1*cos(alfa)
q(9,1)=-M1
q(10,1)=(F2+F3)*sin(alfa)
q(11,1)=-(F2+F3)*cos(alfa)
q(12,1)=M2-M3
q(13,1)=F4*sin(alfa)
q(14,1)=-F4*cos(alfa)
q(15,1)=M4
%Construir vector de fuerzas nodales equivalentes (de las fuerzas aplicadas
%fuera de los nodos:consultar ERM)
%% Matriz de Rigidez
K=zeros(GDL); %declara-crea matriz de rigidez global, inicio vacío
%Construir-ensamblar la matriz de rigidez
%___________________________
%% Elemento 1
Le1=L/2; alfa1=pi/2;
Ke1=[EA/Le1 0 0 -EA/Le1 0 0;...
    0 (12*EI)/(Le1)^3 (6*EI)/(Le1)^2 0 -(12*EI)/(Le1)^3 (6*EI)/(Le1)^2;...
    0 (6*EI)/(Le1)^2 (4*EI)/Le1 0 -(6*EI)/(Le1)^2 (2*EI)/Le1;...
    -EA/Le1 0 0 EA/Le1 0 0;...
    0 -(12*EI)/(Le1)^3 -(6*EI)/(Le1)^2 0 (12*EI)/(Le1)^3 -(6*EI)/(Le1)^2;...
    0 (6*EI)/(Le1)^2 (2*EI)/Le1 0 -(6*EI)/(Le1)^2 (4*EI)/Le1;];
R1=[cos(alfa1) sin(alfa1) 0;-sin(alfa1) cos(alfa1) 0; 0 0 1]; 
Zeros= [0 0 0; 0 0 0; 0 0 0];
Re1=[R1 Zeros; Zeros R1];
Ke1_= Re1'*Ke1*Re1 %global
L1=zeros(GDL,6);
L1(1,1)=1;
L1(2,2)=1;
L1(3,3)=1;
L1(4,4)=1;
L1(5,5)=1;
L1(6,6)=1;
K_e1_=L1*Ke1_*L1';

%% Elemento 2
Le2=L/2; alfa2=pi/2;
Ke2=[EA/Le2 0 0 -EA/Le2 0 0;...
    0 (12*EI)/(Le2)^3 (6*EI)/(Le2)^2 0 -(12*EI)/(Le2)^3 (6*EI)/(Le2)^2;...
    0 (6*EI)/(Le2)^2 (4*EI)/Le2 0 -(6*EI)/(Le2)^2 (2*EI)/Le2;...
    -EA/Le2 0 0 EA/Le2 0 0;...
    0 -(12*EI)/(Le2)^3 -(6*EI)/(Le2)^2 0 (12*EI)/(Le2)^3 -(6*EI)/(Le2)^2;...
    0 (6*EI)/(Le2)^2 (2*EI)/Le2 0 -(6*EI)/(Le2)^2 (4*EI)/Le2;];
R2=[cos(alfa2) sin(alfa2) 0;-sin(alfa2) cos(alfa2) 0; 0 0 1]; 
Zeros= [0 0 0; 0 0 0; 0 0 0];
Re2=[R2 Zeros; Zeros R2];
Ke2_= Re2'*Ke2*Re2 %global
L2=zeros(GDL,6);
L2(4,1)=1;
L2(5,2)=1;
L2(6,3)=1;
L2(7,4)=1;
L2(8,5)=1;
L2(9,6)=1;
K_e2_=L2*Ke2_*L2';
%% Elemento 3
Le3=L*0.5*sqrt(5); alfa3=atan(0.5);
Ke3=[EA/Le3 0 0 -EA/Le3 0 0;...
    0 (12*EI)/(Le3)^3 (6*EI)/(Le3)^2 0 -(12*EI)/(Le3)^3 (6*EI)/(Le3)^2;...
    0 (6*EI)/(Le3)^2 (4*EI)/Le3 0 -(6*EI)/(Le3)^2 (2*EI)/Le3;...
    -EA/Le3 0 0 EA/Le3 0 0;...
    0 -(12*EI)/(Le3)^3 -(6*EI)/(Le3)^2 0 (12*EI)/(Le3)^3 -(6*EI)/(Le3)^2;...
    0 (6*EI)/(Le3)^2 (2*EI)/Le3 0 -(6*EI)/(Le3)^2 (4*EI)/Le3;];
R3=[cos(alfa3) sin(alfa3) 0;-sin(alfa3) cos(alfa3) 0; 0 0 1]; 
Re3=[R3 Zeros; Zeros R3];
Ke3_= Re3'*Ke3*Re3 %global
L3=zeros(GDL,6);
L3(7,1)=1;
L3(8,2)=1;
L3(9,3)=1;
L3(10,4)=1;
L3(11,5)=1;
L3(12,6)=1;
K_e3_=L3*Ke3_*L3';
%% Elemento 4
Le4=L*0.5*sqrt(5); alfa4=atan(0.5);
Ke4=[EA/Le4 0 0 -EA/Le4 0 0;...
    0 (12*EI)/(Le4)^3 (6*EI)/(Le4)^2 0 -(12*EI)/(Le4)^3 (6*EI)/(Le4)^2;...
    0 (6*EI)/(Le4)^2 (4*EI)/Le4 0 -(6*EI)/(Le4)^2 (2*EI)/Le4;...
    -EA/Le4 0 0 EA/Le4 0 0;...
    0 -(12*EI)/(Le4)^3 -(6*EI)/(Le4)^2 0 (12*EI)/(Le4)^3 -(6*EI)/(Le4)^2;...
    0 (6*EI)/(Le4)^2 (2*EI)/Le4 0 -(6*EI)/(Le4)^2 (4*EI)/Le4;];
R4=[cos(alfa4) sin(alfa4) 0;-sin(alfa4) cos(alfa4) 0; 0 0 1];
Zeros= [0 0 0; 0 0 0; 0 0 0];
Re4=[R4 Zeros; Zeros R4];
Ke4_= Re4'*Ke4*Re4 %global
L4=zeros(GDL,6);
L4(10,1)=1;
L4(11,2)=1;
L4(12,3)=1;
L4(13,4)=1;
L4(14,5)=1;
L4(15,6)=1;
K_e4_=L4*Ke4_*L4';
%% Elemento 5
Le5=L; alfa5=-pi/2;
Ke5=[EA/Le5 0 0 -EA/Le5 0 0;...
    0 (12*EI)/(Le5)^3 (6*EI)/(Le5)^2 0 -(12*EI)/(Le5)^3 (6*EI)/(Le5)^2;...
    0 (6*EI)/(Le5)^2 (4*EI)/Le5 0 -(6*EI)/(Le5)^2 (2*EI)/Le5;...
    -EA/Le5 0 0 EA/Le5 0 0;...
    0 -(12*EI)/(Le5)^3 -(6*EI)/(Le5)^2 0 (12*EI)/(Le5)^3 -(6*EI)/(Le5)^2;...
    0 (6*EI)/(Le5)^2 (2*EI)/Le5 0 -(6*EI)/(Le5)^2 (4*EI)/Le5;];
R5=[cos(alfa5) sin(alfa5) 0;-sin(alfa5) cos(alfa5) 0; 0 0 1]; 
Zeros= [0 0 0; 0 0 0; 0 0 0];
Re5=[R5 Zeros; Zeros R5];
Ke5_= Re5'*Ke5*Re5 %global
L5=zeros(GDL,6);
L5(13,1)=1;
L5(14,2)=1;
L5(15,3)=1;
L5(16,4)=1;
L5(17,5)=1;
L5(18,6)=1;
K_e5_=L5*Ke5_*L5';
%% Elemento 6
Le6=L; alfa6=-pi/2;
Ke6=[EA/Le6 0 0 -EA/Le6 0 0;...
    0 (12*EI)/(Le6)^3 (6*EI)/(Le6)^2 0 -(12*EI)/(Le6)^3 (6*EI)/(Le6)^2;...
    0 (6*EI)/(Le6)^2 (4*EI)/Le6 0 -(6*EI)/(Le6)^2 (2*EI)/Le6;...
    -EA/Le6 0 0 EA/Le6 0 0;...
    0 -(12*EI)/(Le6)^3 -(6*EI)/(Le6)^2 0 (12*EI)/(Le6)^3 -(6*EI)/(Le6)^2;...
    0 (6*EI)/(Le6)^2 (2*EI)/Le6 0 -(6*EI)/(Le6)^2 (4*EI)/Le6;];
R6=[cos(alfa6) sin(alfa6) 0;-sin(alfa6) cos(alfa6) 0; 0 0 1]; 
Zeros= [0 0 0; 0 0 0; 0 0 0];
Re6=[R6 Zeros; Zeros R6];
Ke6_= Re6'*Ke6*Re6 %global
L6=zeros(GDL,6);
L6(16,1)=1;
L6(17,2)=1;
L6(18,3)=1;
L6(19,4)=1;
L6(20,5)=1;
L6(21,6)=1;
K_e6_=L6*Ke6_*L6';
K=K_e1_+K_e2_+K_e3_+K_e4_+K_e5_+K_e6_
%% Condiciones de contorno 
%Crear vector de desplazamiento vacío, una sola columna
U=zeros(GDL,1);
%Condiciones de contorno esenciales: vector de nodos restricciones (1,2,3)
%n_contorno=(1,4);
%grados de libertad restringidos y libres 
u_contorno=[1:gdln, (nnodos*gdln-gdln+1):nnodos*gdln]; %nodos primero y último
%nodos libres de desplazamiento, la función setdiff de matlab los devuelve
u_libres=setdiff((1:length(U))',u_contorno); %elimino los nodo-contorno
%crear vector de fuerzas nodales vacío
f=zeros(GDL,1);
%define vector de fuerzas nodales exteriores vacío
FL=f(u_libres);
%define puntos de carga:por ejemplo en Nodo 3
Px=0; f(3*gdln-2)=Px; %no pongo carga en nodo 3
Py=0; f(3*gdln-1)=Py; %no pongo carga en nodo 3
Pz=0; f(3*gdln)=Pz; %no pongo carga en nodo 3
%% SOLUCIÓN MATRICIAL
%Condensación de matrices
KLL=K(u_libres,u_libres); %matriz reducida por las condiciones de contorno.
KRL=K(u_libres,u_contorno); %matriz de acoplamiento
KRR=K(u_contorno, u_contorno); %matriz de rigidez en gdl restringidos
FL=f(u_libres)+q(u_libres); %vector de fuerzas nodales generalizadas P=f+q
%Cálculo de incógnitas desplazamiento en nodos libres
uL=KLL\FL
U(u_libres)=uL
%Desplazamientos globales "U"
Ux=U(1:3:(nnodos*gdln)); %vector de desplazamientos x
Uy=U(2:3:(nnodos*gdln)); %vector de desplazamientos y
Uz=U(3:3:(nnodos*gdln)); %vector de giros z
%Desplazamientos pedidos
ux=Ux(4);
uy=Uy(4);
uz=Uz(4);
u_quest=[ux;uy;uz]
%% datos pedidos
f=K*U;
fR=f(u_contorno)
%Reacciones pedidas
f_quest=[f(1);f(2);f(3)]


%% Representar gráficamente la estructura original y la deformada
