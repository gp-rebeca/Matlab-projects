%% Viga
%% Datos de entrada
L=1; a0=(L/20); E=100; q0=10; P=1;
%Hola, la I cambia con xi
%% Malla a discretizar
nelems=5;
nnodos=nelems+1; %elementos de dos nodos
coords=0:L/nelems:L;

%% Matrices globales
%Rigidez
K=zeros(3*nnodos); %18x18
%Fuerza nodal equivalente
F=zeros(3*nnodos,1); %18x1

%% Diagrama de momento flector en toooda la viga. 

%�C�mo determinamos el momento? Se saca de la ec de equilibrio fuerte, 
%con la ecuaci�n de la El�stica. M=EI*d^2w/dx^2.
%Guardar informaci�n en B. �El momento lo obtenemos en nodos o en ptos de 
%integraci�n? En xi, pto de integraci�n. Si queremos diagrama de momentos 
%flectores, B cambia, la almacenamos para utilizarla seg�n la necesitamos 
%para cada elemento. Como I de xi tambi�n cambia con cada elemento, tambi�n
%hay que almacenar la que calculamos de cada elemento. Se almacenan en los 
%ptos de integraci�n. La tensi�n normal de una viga se calcula con el 
%momento flector. La tensi�n normal de una viga en los ptos de integraci�n 
%se calcula con el momento flector en los ptos de integraci�n.

%Matriz para guardar las informaciones de matriz cinem�tica (B) de cada pto de integraci�n
BB=zeros(1,4,nelems,2);
%Matriz de n dimensiones. En este caso, 4. Las dos primeras para el vector,
%siguiente n�mero de elementos y luego el n�mero de nodos. Vamos a usar la segunda derivada, 
%as� que nos da una funci�n de grado LINEAL.

%Matriz para almacenar momento de inercia en los puntos de integraci�n.
In=zeros(nelems,2);
%�Cu�ntos valores de In hay que almacenar por elemento? Dos, porque tenemos dos ptos de integraci�n. 
%OJO, no es uno por elemento. El momento de inercia es un escalar.

%Deformaci�n en ptos de integraci�n
%Para integrar la rigidez elemental AXIAL necesitamos un pto de integraci�n.
int1=0.;
w1=2.0;
%Dos ptos de integraci�n para rigidez a FLEXI�N 
int2=[-1/sqrt(3),1/sqrt(3)];
w2=[1.0,1.0];
%Para la parte de flexi�n usamos funciones herm�ticas.
%Para la q no herm�tica, lagrangiana. Como q lineal y N orden 3=> Orden 4, 3 ptos de
%integraci�n.
int3=[-sqrt(3/5),0,sqrt(3/5)];
w3=[-5/9;8/9;5/9];

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
    %Calcular separadas la rigidez a flexi�n y la axial y luego juntarlas en la elemental
    
    %Matriz de rigidez axial
    Ka=zeros(2,2); 
    
    %Matriz de rigidez a flexi�n
    Kf=zeros(4,4);
    %Calcular rigiez axial
    
    %Necesito la B 
    B1=1/J*[-0.5,0.5];
    
    %���CALCULAR EL �REA ES PREGUNTA DE EXAMEN!!! �rea en nodos. 
    %Buscar la ley de �reas. Nos dan el ancho de secci�n y tenemos que 
    %ponerlo al cuadrado.
    A1=(a0-2*a0*coords(nodo1)/(3*L))^2;
    A2=(a0-2*a0*coords(nodo2)/(3*L))^2;
    N1=0.5*(1-int1);
    N2=0.5*(1+int1);
    %�rea interpolada
    A=[N1,N2]*[A1;A2]; %escalar
    
    Ka=w1*B1'*E*A*B1*J; %2x2
    
    %Calcular rigidez a flexi�n
    for j=1:2 %necesitamos 2 ptos de integraci�n
        int=int2(j);
        w=w2(j);
        %Vamos a sumar las matrices en un bucle. Mucho m�s elegante.
        %Calcular el momento de inercia en el punto de integraci�n
        I1=(a0-2*a0*coords(nodo1)/(3*L))^4/12; %ancho por 4 entre 12 porque es de secci�n recta cuadrada
        I2=(a0-2*a0*coords(nodo2)/(3*L))^4/12;
        N1=0.5*(1-int);
        N2=0.5*(1+int);
        I=[N1,N2]*[I1;I2];
        In(i,j)=I; %guardar en matriz In
        
        %Calcular matriz cinem�tica
        B2=J^(-2)*[0.25*6*int,J*(0.25*6*int-0.25*2),...
            -0.25*6*int,J*(0.25*6*int+0.25*2)]; 
        %Vector cuatro columnas. La J de dentro es por el giro. Lo otro es la flecha. 
        %Segunda derivada de forma herm�tica entre x. B flexi�n diapositiva 59
        %Diapositiva 58 ���EL A�O PASADO HICERON UNA PREGUNTA SOBRE EL
        %RECUADRO DE LOS POLINOMIOS HERM�TICOS!!!
        
        BB(:,:,i,j)=B2;
        %La guardamos ah�
        
        %Calcular rigiez a flexi�n 
        Kf=Kf+w*B2'*E*I*B2*J;
    end
    
%Ensamblar en la matriz de rigidez de elementos.
Ke([1,4],[1,4])=Ke([1,4],[1,4])+Ka %Un elemento tiene 6 gdl 
%El resto de grados de libertad para la rigidez a flexi�n
Ke([2,3,5,6],[2,3,5,6])=Ke([2,3,5,6],[2,3,5,6])+Kf;
%Fuerza nodal equivalente
Fe=zeros(4,1);

for h=1:3
    int=int3(h);
    w=w3(h);
    %Calcular q interpolando
    q1=q0*coords(nodo1)/L;
    q2=q0*coords(nodo2)/L;
    %Funci�n de carga en los nodos y luego las funciones de forma
    %se llaman as� porque N1 y N2 son las lagrangianas de antes
    N3=0.25*int^3-0.75*int+0.5;
    N4=0.25*int^3-0.25*int^2-0.25*int+0.25;
    N5=-0.25*int^3+0.75*int+0.5;
    N6=0.25*int^3+0.25*int^2-0.25*int-0.25;
    q=[N1,N2]*[q1;q2]; %lagrangianas sin trasponer, las herm�ticas en traspuesta
    Fe=Fe+w*[N3,J*N4,N5,J*N6]'*q*J;
end

%Comprobar en casa que, de este bucle, el vector de fuerzas nodales equivalentes resultante
%es el que he copiado en mis apuntes ��

%Ensamblaje
%Grados de libertad correspondientes a este elemento
gdl1=[3*(i-1)+1:3*(i+1)]; %vector donde almacenamos los gdl de cada elemento. 
%Los recolocamos de forma global en la K de abajo
K(gdl1,gdl1)=K(gdl1,gdl1)+Ke;
gdl2=[3*(i-1)+2,3*(i-1)+3,3*(i-1)+5,3*(i-1)+6]; 
F(gdl2)=F(gdl2)+Fe;
end

%gdl restringidos
gdl3=[1,2,3];
%gdl libres
gdl4=setdiff(1:3*nnodos,gdl3); %te compara con otro vector o algo as�
%reducir matriz
KLL=K(gdl4,gdl4);
FLL=F(gdl4);
%Vector de desplazamiento
u=zeros(3*nnodos,1);

%% Resolver sistema de ecuaciones
uLL=KLL\FLL;
u(gdl4)=uLL; %EN TODOS LOS GDL4 SE A�ADE ULL

%% Calcular momento flector
M=zeros(nelems*2,1);
k=1;
for i=1:nelems
    ww=u([3*(i-1)+2,3*(i-1)+3,3*(i-1)+5,3*(i-1)+6]);
    for j=1:2
        M(k)=E*In(i,j)*BB(:,:,i,j)*ww;
        k=k+1;
    end
end

%% Representar diagrama de momento
%plot(M) %Lo pongo as� porque no quiero que me salga la maldita gr�fica x.x
        
