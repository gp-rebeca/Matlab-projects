%% INTEGRACIÓN SELECTIVA Y REDUCIDA. 
%% Datos de entrada
E=100; v=0.3; G=E/(2*(1+v)); 
L=1; b=L/20; h=L/150; A=b*h; %sección rectangular
q0=10;

%% Vamos a discretizar una malla de 5 elementos
nelems=20; %Con 5 elementos no nos va a salir, con 20 sí.
nnodos=nelems+1;
coords=0:L/nelems:L;
%Vamos a usar un elemento Viga de Timoshenko, su diferencia con la viga de
%Euler-Bernoulli es que tiene en cuenta la cortante. Tendremos rigidez a
%flexión y rigidez a cortante.
%Rigidez a flexión se calcula de forma parecida a la Euler-Bernoulli. 
%Usamos la flexión de forma lineal de un elemento de dos nodos, tanto para 
%flexión como para giros, elemento isoparamétrico.
%Las funciones de forma cortante de un elemento viga de Timoshenko son
%gamma=derivada de la flecha menos el giro (gamma=dw/dx-theta).

%% Integración selectiva: la flexión se integra con 1 pto y la cortante con 2
%Punto de integración para rigidez a flexión
int1=0;
w1=2;
%Punto de integración para rigidez a cortante y fuerza
int2=[-1/sqrt(3), 1/sqrt(3)];
w2=[1,1];
%Matriz de rigidez global
K=zeros(3*nnodos);
%Fuerza
F=zeros(3*nnodos,1);

%Bucles para mis rigideces
for i=1:nelems
    
    %Nodos del elemento i
    nodo1=i; nodo2=i+1;
    %Longitud elemental
    dl=coords(nodo2)-coords(nodo1);
    %Jacobiano
    J=dl/2;
    %Matriz de rigidez elemental
    Ke=zeros(6,6);
    
    %Rigidez a cortante
    %Necesita 2 ptos de integración, mi matriz vacía va a ser 4x4
    Kc=zeros(4);
    
    for j=1:2 %bucle para la integración con 2 ptos
        
        int=int2(j);
        w=w2(j);
        %Funciones de forma
        N1=0.5*(1-int);
        N2=0.5*(1+int);
        dN1=-0.5;
        dN2=0.5;
        %Matriz cinemática cortante
        Bc=[J^(-1)*dN1,-N1, J^(-1)*dN2,-N2];
        %Con esto ya podemos integrar Kc
        Kc=Kc + w*Bc'*G*(5/6)*A*Bc*J;
        
    end
    
    %Rigidez a flexión, 1 pto de integración
    Bf=J^(-1)*[0,-0.5,0,0.5];
    Kf=w1*Bf'*E*b*h^3/12*Bf*J; %lo del momento de inercia de la barra (I=b*h^3/12)
    %Rigidez axial
    Ba=J^(-1)*[-0.5,0.5]; %matriz cte
    Ka=w1*Ba'*E*A*Ba*J; %solo necesito 1 pto de integración
    
    %Ensamblaje rigidez elemental
    Ke([1,4],[1,4])=Ka;
    Ke([2,3,5,6],[2,3,5,6])=Kc+Kf;
    %Ensamblaje rigidez global
    gdl1=3*(i-1)+1:3*(i+1);
    K(gdl1,gdl1)=K(gdl1,gdl1)+Ke;
    
end

%Vamos a aplicar una carga puntual antes de afrontar la complejidad de la
%distribuida juas juas
F(end-1)=-10; % P vertical hacia abajo en el último nodo

%Vector de desplazamiento
u=zeros(3*nnodos,1); 

%Reducir gdl restringidos
gdl2=[1,2,3]; %restringidos
gdl3=setdiff(1:3*nnodos,gdl2); %libres
KLL=K(gdl3,gdl3);
FLL=F(gdl3);
u(gdl3)=KLL\FLL %uL
plot(coords,u(2:3:end)) %solo quiero dibujar la parte vertical que guarda u

%Flecha
flecha= (-10*L^3)/(3*E*(b*h^3/12)) %w=PL^3/(3EI)
u(end-1) %Con 5 elementos me da -1,39*10^5, más pequeño de lo que debería. 
%La viga de Timoshenko tendría que darnos valores más grandes. 
%Si aumentamos el número de elementos, el resultado sigue dando 100 veces 
%menos de lo que debería (u=flecha=-2,7*10^7). Usar integración reducida 
%(sin que llegue a ser un mecanismo).
%¡¡¡PREGUNTA DE EXAMEN!!! Si con elemento hermítico usamos integración
%reducida, la matriz K arroja ceros en la diagonal, se convierte en
%mecanismo y no funciona. Con la viga de Timoshenko no pasa, sale bien.


%% Integración reducida
%Usamos 1 pto de integración (int1, ya definido)
%Matriz de rigidez global
Kr=zeros(3*nnodos);
%Fuerza
Fr=zeros(3*nnodos,1);

%Bucles para mis rigideces
for i=1:nelems
    
    %Nodos del elemento i
    nodo1=i; nodo2=i+1;
    %Longitud elemental
    dl=coords(nodo2)-coords(nodo1);
    %Jacobiano
    J=dl/2;
    %Matriz de rigidez elemental
    Ke=zeros(6,6);
    %Rigidez a cortante
    Kc=zeros(4);
    N1=0.5*(1-int1);
    N2=0.5*(1+int1);
    dN1=-0.5;
    dN2=0.5;
    %Matriz cinemática cortante
    Bc=[J^(-1)*dN1,-N1, J^(-1)*dN2,-N2];
    %Rigidez cortante
    Kc= w1*Bc'*G*(5/6)*A*Bc*J;
    %Rigidez a flexión, 1 pto de integración
    Bf=J^(-1)*[0,-0.5,0,0.5]; %cte
    Kf=w1*Bf'*E*b*h^3/12*Bf*J; %lo del momento de inercia de la barra (I=b*h^3/12)
    %Rigidez axial
    Ba=J^(-1)*[-0.5,0.5]; %matriz cte
    Ka=w1*Ba'*E*A*Ba*J; %solo necesito 1 pto de integración
    
    %Ensamblaje rigidez elemental
    Ke([1,4],[1,4])=Ka;
    Ke([2,3,5,6],[2,3,5,6])=Kc+Kf;
    %Ensamblaje rigidez global reducida
    gdl1=3*(i-1)+1:3*(i+1);
    Kr(gdl1,gdl1)=Kr(gdl1,gdl1)+Ke;
    
end

%Fuerza reducida 
Fr(end-1)=-10; %end es momento, end-1 es coordenada y

%Vector de desplazamiento reducido
ur=zeros(3*nnodos,1)

%Reducir gdl restringidos
gdl2=[1,2,3]; %restringidos
gdl3=setdiff(1:3*nnodos,gdl2); %libres
%Matrices libres reducidas
KLLr=Kr(gdl3,gdl3);
FLLr=Fr(gdl3); 
ur(gdl3)=KLLr\FLLr
hold on %retener una gráfica cuando queremos que nos escupa otra
plot(coords,ur(2:3:end),'r') %la r es de rojo :)
ur(end-1) %nos tiene que dar la calculada en papel
%El número de elementos para el valor exacto tiene que estar entre 5 y 20.

%Fuerzas de reacción
fr=Kr*ur %Si el resultado está bien, nuestro vector de fuerzas de reacción 
%tiene que valer lo mismo que P
%KRLr=Kr(gdl2,gdl3);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% MÉTODO DE PENALIZACIÓN

%En vez de quitar columnas y filas, restricción cinemática
C=zeros(3,3*nnodos);
r=zeros(3,1);
kpen=1000; 
k=kpen*eye(3); %eye= matriz identidad de 3x3. 
%Otra forma: k=diag([1000,1000,1000]) 3x3, contiene el parámetro de 
%penalización. Pongo matriz diagonal de valor 100, puede que haya que 
%hacerla más grande o más pequeña, no sabemos ¯\_(?)_/¯.

%Impongo mi restricción
C(1,1)=1;
C(2,2)=1;
C(3,3)=1;
%Matriz de rigidez penalizada
Kp=Kr+C'*k*C; %¿Aumenta su tamaño? No. Con Lagrange sí, y aumentan 
%las incógnitas, lo veremos más abajo
Fp=Fr+C'*k*r; %Tampoco cambia su tamaño
%Desplazamientos
up=Kp\Fp
plot(coords,up(2:3:end),'g--')
%Aumentando kpen sigue funcionando, incluso poniendo kpen=1000000. Pero es
%porque tenemos pocos gdl, en otros problemas nos vamos a cagar.
%TRUQUITO: 1000*valor máximo de los términos de la diagonal ppal
%(kpen=1000*max(max(Kr)), genera un número que va a estar dentro de lo
%correcto. Si con 1000 no fufa, probamos con 10000, 100000, etc.

up-ur %no son iguales. Un método y otro no dan el mismo resultado, pero es 
%cercano. Cambiando la penalización, la diferencia cambia. 

%% MULTIPLICADORES DE LAGRANGE
%Misma C y misma r. Nueva matriz de rigidez que incluye la original, la de
%coeficientes, C', C y matriz de ceros de tamaño 3x3
%El inconveniente de este método es que sí aumenta el tamaño de las
%matrices. MUY CARO COMPUTACIONALMENTE. Peeeeero es exactooo. SÍ FUNCIONA.
%KML=zeros(66,66);
%KML(1:63,1:63)=K;
%KML(1:18,64:66)=C';
%KML(64:66,1:18)=C;
%KML(64:66,64:66)=zeros(3);
KML=[Kr C'; C zeros(3)];
FML=[Fr;r];
uML=KML\FML %contiene los multiplicadores de Lagrange
uM=uML(1:end-3) %tres son los multiplicadores de Lagrange/restricciones que hemos puesto
plot(coords,uM(2:3:end),'m') %la m es de maaaalvaaaa
uML(end-2:end) %NUESTRAS REACCIONES
