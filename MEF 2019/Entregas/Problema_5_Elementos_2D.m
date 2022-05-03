%% PROBLEMA 5. ELEMENTOS 2D
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Datos 
E=210*10^9; v=0.25; L=1; h=0.2; P=-1000*10^3;

%% Malla
nelems=10;
nnodos= (nelems/2+1)*3; %18
%Viga de sección constante
nodos=[1 0 -h/2 ;...
       2 L/5 -h/2;...
       3 2*L/5 -h/2;...
       4 3*L/5 -h/2;...
       5 4*L/5 -h/2;...
       6 L -h/2;...
       7 0 0;...
       8 L/5 0;...
       9 2*L/5 0;...
       10 3*L/5 0;...
       11 4*L/5 0;...
       12 L 0;...
       13 0 h/2;...
       14 L/5 h/2;...
       15 2*L/5 h/2;...
       16 3*L/5 h/2;...
       17 4*L/5 h/2;...
       18 L h/2;]; %matriz de coordenadas de nodos. Primera columna,
       %número de nodos; segunda columna, x, y tercera columna, y.
       %Numerados de abajo a arriba y de izquierda a
       %derecha, en sentido antihorario.
elems=[1 1 2 8 7;...
       2 2 3 9 8;...
       3 3 4 10 9;...
       4 4 5 11 10;...
       5 5 6 12 11;...
       6 7 8 14 13;...
       7 8 9 15 14;...
       8 9 10 16 15;...
       9 10 11 17 16;...
       10 11 12 18 17;]; %matriz de elementos con sus nodos
%Puntos de integración
int2=[-1/sqrt(3), 1/sqrt(3)];
w2=[1,1];
%Grados de libertad
ngdl=nnodos*2; %estamos en 2D, 2gdl por nodo 36
%Declaramos matriz vacía de rigidez
K=zeros(ngdl); %36x36
BB=zeros(3,8,10,2,2); %3,8 dimensión de B, 10 número de elementos, 2 pto integración chi, 2 pto integración eta
%Bucles de elementos
for i=1:nelems
    %Primero identificamos en qué elemento estoy y qué nodos lo forman 
    nodos_=elems(i,2:end); %fila i=1 sería (1 1 2 11 12) de la columna 2 a la 5 porque
    %la columna 1 es la que guarda el nombre del elemento
    coords=nodos(nodos_,[2,3]);
    %Matriz de rigidez elemental, 4 nodos, dimensión 8x8
    Ke=zeros(8);
    %Usamos 2 ptos de integración para cada eje, por lo tanto, dos bucles
    for j=1:2
        chi=int2(j);
        wchi=w2(j); %guaaaachiiii d:
        %La K es integral doble, tengo un doble sumatorio
        aux=zeros(8); %matriz auxiliar
        for s=1:2
            eta=int2(s);
            weta=w2(s);
            R=[-0.25*(1-eta),0.25*(1-eta), 0.25*(1+eta),-0.25*(1+eta);...
               -0.25*(1-chi),-0.25*(1+chi),0.25*(1+chi),0.25*(1-chi)]; %matriz que necesito para calcular J'. Está formada por las derivadas de las funciones de forma
           %Matriz jacobiana traspuesta
           JT=R*coords;
           %Jacobiano
           J=det(JT'); %el determinante de J sin trasponer
           %Calculo la derivada de las funciones de forma 
           dN=JT\R; %la inversa de la traspuesta por R (con esto \ ya se entiende que es la inversa)
           %Matriz cinemática
           B=[dN(1,1) 0 dN(1,2) 0 dN(1,3) 0 dN(1,4) 0;...
               0 dN(2,1) 0 dN(2,2) 0 dN(2,3) 0 dN(2,4);...
              dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3) dN(2,4) dN(1,4)];
           BB(:,:,i,j,s)=B; 
           D=E/(1-v^2)*[1 v 0; v 1 0; 0 0 0.5*(1-v)];
          %Con esto ya puedo hacer el PRIMER sumatorio de la matriz K
          %elemental
          aux=aux+weta*B'*D*B*J;
          %Mais tenemos dos sumatorios
        end
        %SEGUNDO sumatorio
        Ke=Ke+wchi*aux;
    end
    %Hasta aquí, tendríamos la matriz K elemental
    %Toca ensamblar, para ello necesitamos los índices correspondientes a
    %los gdl de los nodos. ¿Cómo determinamos esa mierda? Buscamos un
    %algoritmo que pueda relacionar. Un elemento está formado por 4 nodos.
    %En el elemento 1 tenemos los nodos 1, 2, 7 y 8, y los gdls 1,2, 3,4,
    %13,14, 15,16
    gdl1=zeros(8,1); %ocho índices (4 nodos con 2 gdl por nodo)
    for m=1:4 %los cuatro nodos que forman un elemento
        %Vamos rellenando gdl1
        gdl1(2*(m-1)+1)=(nodos_(m)-1)*2 +1;
        gdl1(2*(m-1)+2)=(nodos_(m)-1)*2 +2; %la siguiente componente, los gdl del nodo 2
        %LA RAYADA MÁXIMA
    end
    K(gdl1,gdl1)=K(gdl1,gdl1)+Ke
end
%% Carga aplicada
%Nodo 13 hacia abajo. Posición 26
F=zeros(ngdl,1);
F(36)=P;
%Vamos a buscar los gdl restringidos
gdl2=[1,2,13,14,25,26];
%gdl libres
gdl3=setdiff(1:ngdl,gdl2); %36-6=30
%Reducir la matriz de rigidez por BC
KLL=K(gdl3,gdl3);
FLL=F(gdl3);
%% Desplazamiento
u=zeros(ngdl,1);
u(gdl3)=KLL\FLL;
%% Flecha viga
wmef=u(36);
w=abs(u(36))
%% Reacciones en el nodo 1
Reacciones=zeros(ngdl,1);
Reacciones=K*u;
R1=zeros(2,1);
R1=Reacciones(1:2)

%% Tensión y Deformación
sigma=zeros(3,30);
nodos__=elems(1,2:end);
    for r=1:4
        gdl4(2*(r-1)+1)=(nodos__(r)-1)*2 +1;
        gdl4(2*(r-1)+2)=(nodos__(r)-1)*2 +2;
    end
u_=u(gdl4);
B_=BB(:,:,1,1,1);
eps=B_*u_
sigma=D*eps

%% Flecha obtenida por Resis
b=1; %espesor unitario
I=b*h^3/12;
wr=P*L^3/(3*(E*I));

    
