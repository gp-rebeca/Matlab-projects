%% CONDENSACIÓN ESTÁTICA DE MATRICES

%Diapositiva 86.
%Viga de Timoshenko de longitud L. Se aplican dos cargas: P1, horizontal
%postiva en x=L, y P2, vertical negativa aplicada en x=L/2. Se dispone un
%desplazamiento de valor u=0.1. ¿Cuánto vale la fuerza libre equivalente de
%este resultado?

%                  | P2             |  u
%\|                v                v
%\|o===============o================o---> P1
%\|

%% Datos
E=100; A=100; I=1; L=1;
P1=10; P2=-1; u_carga=-0.1;

%% Resolución 
nelems=2;
nnodos=3;
gdl=3*nnodos; %3gdl por nodo
gdltotales=[1:gdl]; %vector de gdl
%Matrices elementales
Ke1=zeros(gdl);
Ke2=Ke1;
Ke1(1:6,1:6)=Ke; %no tenemos Ke definida
Ke2(4:9,4:9)=Ke;
K=Ke1+Ke2;

%Definimos el vector de desplazamientos y el de fuerzas como vectores
%vacíos de dimensión 9x1
U=zeros(gdl,1);
F=zeros(gdl,1);

%Primera reducción
gdlcontorno=[1 2 3];
gdlcarga=8;
gdllibres1=setdiff(gdltotales,gdlcontorno); 
KLL1=K(gdllibres1,gdllibres1);
KLR1=K(gdllibres1,gdlcontorno);
KRL1=KLR1';
KRR1=K(gdlcontorno,gdlcontorno);
FL1=F(gdllibres);
FL1(4)=P1;
FL1(2)=P2;

%Segunda reducción
U(gdlcarga)=u_carga;
UR=U(gdlcarga);
gdlUrestringidos=(5);
gdlUlibres2=(1:lenght(gdllibres1));
gdlUlibres2=setdiff(gdlUlibres2,gdlUrestringidos);
KLL2=KLL1(gdlUlibres2,gdlUlibres2);
KLR2=KLL1(gdlUlibres2,gdlUrestringidos);
KRL2=KLR2';
KRR2=KLL1(gdlUrestringidos,gdlUrestringidos);
FL2=FL1(gdlUlibres2);
UL2=KLL2\(FL2-KLR2*UR);
FR2=KRL2*UL2+KRR2*UR;

%Recuperación de restricciones
UL1(gdlUlibres2)=UL2;
UL1(gdlUrestringidos2)=UR;
FL1(gdlUrestringidos2); %fuerza equivalente

%Recuperación de contorno
FR1=KRL1*UL1'
F(gdllibres1)=FL1;
U(gdllibres1)=UL1;

%Alternativa
U(8)=u_carga; %en la coordenada y del último nodo
F(7)=P1;
F(5)=P2;
gdlulibres3=[4 5 6 7 9];
gdlurestringidos3=[1 2 3 8];
KLL3=K(gdlulibres3,gdlulibres3);
KLR3=K(gdlulibres3,gdlurestringidos3);
KRL3=KLR3';
KRR3=K(gdlurestringidos3,gdlurestringidos3);
UR3=U(gdlurestringidos3);
FL3=F(gdlulibres3);
UL3=KLL\(FL-KLR*UR);
FR3=KRL3*UL3+KRR3*UR3;


