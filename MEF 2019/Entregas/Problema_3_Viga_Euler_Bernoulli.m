%% Viga Euler-Bernoulli
%% Datos de entrada
L=7; b=0.3; E=30*10^9; 
P1=-20000; P2x=5000; P2y=-5000; q0=-15000; M=10000;
coordA=2; coordM=5;
%% Malla a discretizar
nelems=7;
nnodos=nelems+1; %elementos de dos nodos
coords=0:L/nelems:L;

%% Matrices globales
%Rigidez
K=zeros(3*nnodos); %24x24
%Fuerza nodal equivalente
F=zeros(3*nnodos,1); %24x1

%% Diagrama de momento flector en toda la viga. 

%Matriz para guardar las informaciones de matriz cinemática (B) de cada pto de integración
BB=zeros(1,4,nelems,2);
%Matriz de n dimensiones. En este caso, 4. Las dos primeras para el vector,
%siguiente número de elementos y luego el número de nodos. Vamos a usar la segunda derivada, 
%así que nos da una función de grado LINEAL.

%Matriz para almacenar momento de inercia en los puntos de integración.
In=zeros(nelems,2);
%¿Cuántos valores de In hay que almacenar por elemento? Dos, porque tenemos dos ptos de integración. 
%OJO, no es uno por elemento. El momento de inercia es un escalar.

%Deformación en ptos de integración
%Para integrar la rigidez elemental AXIAL necesitamos un pto de integración.
int1=0.;
w1=2.0;
%Dos ptos de integración para rigidez a FLEXIÓN 
int2=[-1/sqrt(3),1/sqrt(3)];
w2=[1.0,1.0];
%Para la parte de flexión usamos funciones hermíticas.
%Para la q no hermítica, lagrangiana. Como q lineal y N orden 3=> Orden 4, 3 ptos de
%integración.
int3=[-sqrt(3/5),0,sqrt(3/5)];
w3=[5/9;8/9;5/9];

%Bucle
for i=1:nelems
    nodo1=i; nodo2=i+1;
    %Longitud elemental
    dl=coords(nodo2)-coords(nodo1);
    %Jacobiano
    J=dl/2;
    %Matriz de rigidez elemental 6x6 porque 2 nodos y cada uno 3 gdl. En
    %total, 6gdl.
    Ke=zeros(6,6);
    %Calcular separadas la rigidez a flexión y la axial y luego juntarlas en la elemental
    
    %Matriz de rigidez axial
    Ka=zeros(2,2); 
    
    %Matriz de rigidez a flexión
    Kf=zeros(4,4);
    %Calcular rigiez axial
    
    %Necesito la B 
    B1=1/J*[-0.5,0.5];
    
    %Área en nodos. 
    %Buscar la ley de áreas. Nos dan el ancho de sección y tenemos que 
    %ponerlo al cuadrado.
    A1=(0.6-0.3*coords(nodo1)/L)*b;
    A2=(0.6-0.3*coords(nodo2)/L)*b;
    N1=0.5*(1-int1);
    N2=0.5*(1+int1);
    %Área interpolada
    A=[N1,N2]*[A1;A2]; %escalar
    
    Ka=w1*B1'*E*A*B1*J; %2x2
    
    %Calcular rigidez a flexión
    for j=1:2 %necesitamos 2 ptos de integración
        int=int2(j);
        w=w2(j);
        %Vamos a sumar las matrices en un bucle. Mucho más elegante.
        %Calcular el momento de inercia en el punto de integración
        I1=(b*(0.6-0.3*coords(nodo1)/L)^3)/12; %ancho por 4 entre 12 porque es de sección recta cuadrada
        I2=(b*(0.6-0.3*coords(nodo2)/L)^3)/12;
        N1=0.5*(1-int);
        N2=0.5*(1+int);
        I=[N1,N2]*[I1;I2];
        In(i,j)=I; %guardar en matriz In
        
        %Calcular matriz cinemática
        d2N3=0.25*6*int;
        d2N4=0.25*6*int-0.25*2;
        d2N5=-0.25*6*int;
        d2N6=0.25*6*int+0.25*2;
        B2=J^(-2)*[d2N3,J*d2N4,d2N5,J*d2N6]; 
        %Vector cuatro columnas. La J de dentro es por el giro. Lo otro es la flecha. 
        %Segunda derivada de forma hermítica entre x. B flexión diapositiva 59
        %Diapositiva 58 
        
        BB(:,:,i,j)=B2;
        %La guardamos ahí
        
        %Calcular rigiez a flexión 
        Kf=Kf+w*B2'*E*I*B2*J;
    end
    
%Ensamblar en la matriz de rigidez de elementos.
Ke([1,4],[1,4])=Ke([1,4],[1,4])+Ka; %Un elemento tiene 6 gdl 
%El resto de grados de libertad para la rigidez a flexión
Ke([2,3,5,6],[2,3,5,6])=Ke([2,3,5,6],[2,3,5,6])+Kf;
%Fuerza nodal equivalente
Fe=zeros(4,1);

    for h=1:3
        int=int3(h);
        w=w3(h);
        %Función de carga en los nodos y luego las funciones de forma
        %se llaman así porque N1 y N2 son las lagrangianas de antes
        N1=0.5*(1-int);
        N2=0.5*(1+int);
        N3=0.25*int^3-0.75*int+0.5;
        N4=0.25*int^3-0.25*int^2-0.25*int+0.25;
        N5=-0.25*int^3+0.75*int+0.5;
        N6=0.25*int^3+0.25*int^2-0.25*int-0.25;
        if i<=2
            q1=q0*coords(nodo1)/2;
            q2=q0*coords(nodo2)/2;
        elseif i==3
            q1=2*q0-(q0*coords(nodo1))/2;
            q2=2*q0-(q0*coords(nodo2))/2;
        elseif i==4
            q1=2*q0-(q0*coords(nodo1))/2;
            q2=2*q0-(q0*coords(nodo2))/2;
        else
            q1=0;
            q2=0;
        end
        q=[N1,N2]*[q1;q2];
        Fe=Fe+w*[N3,J*N4,N5,J*N6]'*q*J;
    end

%Ensamblaje
%Grados de libertad correspondientes a este elemento
gdl1=[3*(i-1)+1:3*(i+1)]; %vector donde almacenamos los gdl de cada elemento. 
%Los recolocamos de forma global en la K de abajo
K(gdl1,gdl1)=K(gdl1,gdl1)+Ke
gdl2=[3*(i-1)+2,3*(i-1)+3,3*(i-1)+5,3*(i-1)+6]; 
F(gdl2)=F(gdl2)+Fe;
end
F(17)=F(17)+P1; 
F(18)=F(18)+M;
F(22)=F(22)+P2x;
F(23)=F(23)+P2y;
%gdl restringidos
gdl3=[1,2,20];
%gdl libres
gdl4=setdiff(1:3*nnodos,gdl3); %te compara con otro vector o algo así
%reducir matriz
KLL=K(gdl4,gdl4);
FLL=F(gdl4);
%Vector de desplazamiento
u=zeros(3*nnodos,1);

%% Resolver sistema de ecuaciones
uLL=KLL\FLL;
u(gdl4)=uLL; %EN TODOS LOS GDL4 SE AÑADE ULL


%% Reacción apoyo móvil
R=K*u;
R2=R(20)
%% Giro apoyo extremo izquierdo
phi1=u(3)
%% Desplazamiento horizontal del punto A
ua=u(7)
%% Flecha máxima
uflecha=zeros(1,8);
uflecha(1)=u(2);
uflecha(2)=u(5);
uflecha(3)=u(8);
uflecha(4)=u(11);
uflecha(5)=u(14);
uflecha(6)=u(17);
uflecha(7)=u(20);
uflecha(8)=u(23);
umax=max(uflecha)

