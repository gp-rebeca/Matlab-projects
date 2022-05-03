%% CLASE 10-10-2019
%% ~~~~~~~~~~~~~~~~~
%% Datos
E=30e10; v=0.3; L=4; P=1000;

%% Malla
nelems=5;
nnodos= (nelems+1)*2;
%Viga de sección variable, tenemos que buscar su ley de áreas
nodos=[1 0 -L/10;...
       2 L/5 -L/5*(1-L/5/(2*L))/2;...
       3 2*L/5 -L/5*(1-2*L/5/(2*L))/2;...
       4 3*L/5 -L/5*(1-3*L/5/(2*L))/2;...
       5 4*L/5 -L/5*(1-4*L/5/(2*L))/2;...
       6 6*L/5 -L/5*(1-5*L/5/(2*L))/2;...
       7 0 -L/10;...
       8 L/5 L/5*(1-L/5/(2*L))/2;...
       9 2*L/5 L/5*(1-2*L/5/(2*L))/2;...
       10 3*L/5 L/5*(1-3*L/5/(2*L))/2;...
       11 4*L/5 L/5*(1-4*L/5/(2*L))/2;...
       12 L L/10/2;]; %matriz de coordenadas de nodos PRIMERA COLUMNA NÚMERO DE NODOS, SEGUNDA COORDENAD X Y TERCERA COORDENADA Y
elems=[1 1 2 8 7;...
       2 3 3 9 8;...
       3 3 4 10 9;...
       4 4 5 11 10;...
       5 5 6 12 11;]; %matriz de elementos con sus nodos
%Puntos de integración
int2=[-1/sqrt(3), 1/sqrt(3)];
w2=[1,1];
%Grados de libertad
ngdl=nnodos*2; %estamos en 2D
%Declaramos matriz vacía de rigidez
K=zeros(ngdl);
BB=zeros(3,8,5,2,2); %3,8 dimensión de B, 5 nelems, 2 pto integración chi, 2 pto integración eta
%Bucles de elementos
for i=1:nelems
    %primero identificamos en qué elemento estoy y qué nodos forman ese
    %elemento
    nodos_=elems(i,2:end); %fila i (1 1 2 8 7) de la columna 2 a la 5 porque la columna 1 es la que guarda el nombre del elemento
    coords=nodos(nodos_,[2,3]); %mirar ese pto y coma
    %Matriz de rigidez elemental, 4 nodos, dimensión 8x8
    Ke=zeros(8);
    %Usamos 2 ptos de integración para cada eje, por lo tanto, dos bucles
    for j=1:2
        chi=int2(j);
        wchi=w2(j);
        %La K es integral doble, tengo un doble sumatorio
        aux=zeros(8); %matriz auxiliar
        for h=1:2
            eta=int2(h);
            weta=w2(h);
            R=[0.25*(eta-1),0.25*(1-eta), 0.25*(1+eta),-0.25*(1+eta);...
               0.25*(chi-1),-0.25*(1+chi),0.25*(1+chi),0.25*(1-chi)]; %matriz que necesito para calcular J'. Está formada por las derivadas de las funciones de forma
           %Matriz jacobiana traspuesta
           JT=R*coords;
           %Jacobiano
           J=det(JT'); %el determinante de J sin trasponer
           %Calculo la derivada de las funciones de forma 
           dN=JT\R; %la inversa de la traspuesta por R (co esto \ ya se entiende que es la inversa)
           %Matriz cinemática
           B=[dN(1,1) 0 dN(1,2) 0 dN(1,3) 0 dN(1,4) 0;...
               0 dN(2,1) 0 dN(2,2) 0 dN(2,3) 0 dN(2,4);...
              dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3) dN(2,4) dN(1,4)];
           BB(:,:,i,i,h)=B;
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
    %algoritmo qye pueda relacionar, con el nodo 1 tenemos 1,2, 7,8, 13,14,
    %15,16
    gdl1=zeros(8,1);%ocho índices
    for g1=1:4 %los cuatro nodos que forman un elemento
        %Vamos rellenando gdl1
        gdl1(2*(g1-1)+1)=(nodos(g1)-1)*2 +1;
        gdl1(2*(g1-1)+2)=(nodos(g1)-1)*2 +2; %la siguiente componente, los gdl del nodo 2
        %LA RAYADA MÁXIMA
    end
    K(gdl1,gdl1)=K(gdl1,gdl1)+Ke;
end
%% Carga aplicada
%Nodo 6 hacia abajo. Posición 12
F=zeros(ngdl,1);
F(12)=-P;
%Vamos a buscar los gdl restringidos
gdl2=[1,2,13,14];
%gdl libres
gdl3=setdiff(1:ngdl,gdl2);
%Reducir la matriz de rigidez por BC
KLL=K(gdl3,gdl3);
FLL=F(gdl3);
%% Desplazamiento
u=zeros(ngdl,1);
u(gdl3)=KLL\FLL
%% Tensión
%Matriz tangente por el tensor de deformación
%% Deformación
%Derivada de desplazamientos, equivalente a la B*û. La u nodal la tenemos,
%nos falta la B para cada
