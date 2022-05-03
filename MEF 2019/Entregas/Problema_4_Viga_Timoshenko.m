%% PROBLEMA 4: VIGA TIMOSHENKO
%% INTEGRACIÓN SELECTIVA Y REDUCIDA. 
%% Datos de entrada
E=100; v=0.3; G=E/(2*(1+v)); I=6.75e-04;
Ltot=2; L=1; A=0.09; 
P=-5; q1=-5; q2=-10;
M=2.5;
%% Vamos a discretizar 
nelems=8; 
nnodos=nelems+1;
coords=0:Ltot/nelems:Ltot;

%% Integración reducida
%Ptos de integración 
int1=0;
w1=2;
int2=[-1/sqrt(3),1/sqrt(3)];
w2=[1.0,1.0];
int3=[-sqrt(3/5),0,sqrt(3/5)];
w3=[5/9;8/9;5/9];
%Matriz de rigidez global
Kr=zeros(30); %porque en la rótula tenemos 6gdl
Krv1=zeros(15);
Krv2=zeros(15);
%Fuerza
Fr=zeros(30,1); 
%Bucles para mis rigideces
%VANO 1
for i=1:4
    
    %Nodos del elemento i
    nodo1=i; nodo2=i+1;
    %Longitud elemental
    dl=coords(nodo2)-coords(nodo1);
    %Jacobiano
    J=dl/2;
    %Matriz de rigidez elemental
    Kev1=zeros(6,6);
    %Rigidez a cortante
    %Necesita 2 ptos de integración, mi matriz vacía va a ser 4x4
    Kcv1=zeros(4);
    %Funciones de forma
     N1=0.5*(1-int1);
     N2=0.5*(1+int1);
     dN1=-0.5;
     dN2=0.5;
     %Matriz cinemática cortante
     Bc=[J^(-1)*dN1,-N1, J^(-1)*dN2,-N2];
     %Con esto ya podemos integrar Kc
     Kcv1=Kcv1 + w1*Bc'*G*(5/6)*A*Bc*J;
   
    %Rigidez a flexión, 1 pto de integración
    Bf=J^(-1)*[0,-0.5,0,0.5]; %cte
    Kfv1=w1*Bf'*E*I*Bf*J; %lo del momento de inercia de la barra (I=b*h^3/12)
    %Rigidez axial
    Ba=J^(-1)*[-0.5,0.5]; %matriz cte
    Ka=w1*Ba'*E*A*Ba*J; %solo necesito 1 pto de integración
    
    %Ensamblaje rigidez elemental
    Kev1([1,4],[1,4])=Ka;
    Kev1([2,3,5,6],[2,3,5,6])=Kcv1+Kfv1;
    %Ensamblaje rigidez global reducida
    gdl1=3*(i-1)+1:3*(i+1);
    Krv1(gdl1,gdl1)=Krv1(gdl1,gdl1)+Kev1;
    %Fuerza
    Fe=zeros(4,1);
    for h=1:2
        int=int2(h);
        w=w2(h);
        N1=0.5*(1-int);
        N2=0.5*(1+int);
        qa=q1+(q2-q1)*coords(nodo1)/L;
        qb=q1+(q2-q1)*coords(nodo2)/L;
        q=[N1,N2]*[qa;qb];
        Fe=Fe+w*[N1,0,N2,0]'*q*J;
    end
   gdl2=[3*(i-1)+2,3*(i-1)+3,3*(i-1)+5,3*(i-1)+6]; 
   Fr(gdl2)=Fr(gdl2)+Fe;
end

%VANO 2
for i=1:4
    
    %Nodos del elemento i
    nodo1=i; nodo2=i+1;
    %Longitud elemental
    dl=coords(nodo2)-coords(nodo1);
    %Jacobiano
    J=dl/2;
    %Matriz de rigidez elemental
    Kev2=zeros(6,6);
    %Rigidez a cortante
    %Necesita 2 ptos de integración, mi matriz vacía va a ser 4x4
    Kcv2=zeros(4);
    %Funciones de forma
    N1=0.5*(1-int1);
    N2=0.5*(1+int1);
    dN1=-0.5;
    dN2=0.5;
    %Matriz cinemática cortante
    Bc=[J^(-1)*dN1,-N1, J^(-1)*dN2,-N2];
    %Con esto ya podemos integrar Kc
    Kcv2=Kcv2 + w1*Bc'*G*(5/6)*A*Bc*J;
    %Rigidez a flexión, 1 pto de integración
    Bf=J^(-1)*[0,-0.5,0,0.5]; %cte
    Kfv2=w1*Bf'*E*I*Bf*J; %lo del momento de inercia de la barra (I=b*h^3/12)
    %Rigidez axial
    Ba=J^(-1)*[-0.5,0.5]; %matriz cte
    Ka=w1*Ba'*E*A*Ba*J; %solo necesito 1 pto de integración
    
    %Ensamblaje rigidez elemental
    Kev2([1,4],[1,4])=Ka;
    Kev2([2,3,5,6],[2,3,5,6])=Kcv2+Kfv2;
    %Ensamblaje rigidez global reducida
    gdl1=3*(i-1)+1:3*(i+1);
    Krv2(gdl1,gdl1)=Krv2(gdl1,gdl1)+Kev2;
end

%Ensamblaje rigidez global
Kr(1:15,1:15)=Krv1;
Kr(16:30,16:30)=Krv2;
%Fuerza reducida;
Fr(15)=0;
Fr(18)=0;
Fr(20)=Fr(20)+P; 
Fr(27)=Fr(27)+M;
%Vector de desplazamiento reducido
ureduc=zeros(30,1);

%Reducir gdl restringidos
gdl3=[1,2,3,29]; %restringidos
gdl4=setdiff(1:30,gdl3); %libres
%Matrices libres reducidas
KLLr=Kr(gdl4,gdl4);
FLLr=Fr(gdl4); 
ureduc(gdl4)=KLLr\FLLr
%Fuerzas de reacción
%fr=Kr*ureduc %Si el resultado está bien, nuestro vector de fuerzas de reacción 
%tiene que valer lo mismo que P
%KRLr=Kr(gdl2,gdl3);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% MÉTODO DE PENALIZACIÓN

%En vez de quitar columnas y filas, restricción cinemática
C=zeros(6,30);
r=zeros(6,1);
up=zeros(30,1);
kpen=1.e+06; 
k=kpen*eye(6); %eye= matriz identidad de 3x3. 
%Otra forma: k=diag([1000,1000,1000]) 3x3, contiene el parámetro de 
%penalización. Pongo matriz diagonal de valor 100, puede que haya que 
%hacerla más grande o más pequeña, no sabemos ¯\_(?)_/¯.

%Impongo mi restricción
C(1,1)=1;
C(2,2)=1;
C(3,3)=1;
C(4,13)=1;
C(5,14)=1;
C(4,16)=-1;
C(5,17)=-1;
C(6,29)=1;
%Matriz de rigidez penalizada
Kp=Kr+C'*k*C; %¿Aumenta su tamaño? No. Con Lagrange sí, y aumentan 
%las incógnitas, lo veremos más abajo
Fp=Fr+C'*k*r; %Tampoco cambia su tamaño
%Desplazamientos
up=Kp\Fp
%Reacciones
freacc=Kr*up;

%% Resultados
ur=up(14)
phi1=up(15)
phi2=up(18)
Mr=freacc(3)
uc=up(8)
phic=up(9)

