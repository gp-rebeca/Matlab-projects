%% Datos
L=1; %longitud total
%% Definir primero la geometría
%Coordenadas nodos
nelems=4; %número elementos
nnodes=nelems+1; %número nodos
%tenemos una barra horizontal origen 0,0, longitud L. 
coords=[0:L/nnodes:L;] %Coordenadas nodos, incremento longitud del elemneto
%%  Propiedades barra
E=100; %elástico
A=1; %área sección cte
%¿Cuántos ptos de integraCIÓN NECESITAMOS PARA LA MATRIZ DE RIGIDEZ?
%Un elemento tiene 2 nodos. Función de forma lineal. B matriz de
%compatibilidad. ESTO ENTRA EN EL EXAMEN. Un pto de integración para la K.
%Ptos de integración
int1=0.;
w1=2.0;
%Cuadratura de orden 2
int2=[-1/sqrt(3), 1/sqrt(3)]; %almacena la función de integración
%¿Cuántos grados de libertad tenemos si tenemos 5 nodos? Cada nodo tiene 1 gdl. 
%(Tiene dos pero solo uno es de deformación, el axil. El otro es de sólido
%rígido).
ngdl=nnodes*1;
%Matriz de rigidez total
K=zeros(ngdl); %5x5
%La fuerza nodal total
F=zeros(ngdl,1);
%Declarar fuerzas sobre estructura
q0=10; P=2;
%Bucles para generar y ensamblar la matriz de rigidez
for i=1:nelems %si i es 1, tenemos 1 y 2
    node1=i; node2=i+1;
    %longitud de ese elemento i
    dl=coords(node2)-coords(node1);
    %Jacobiano-para transformar de x a chi
    J=dl/2; %longitud elemento entre 2
    %Determinar la matriz de rigidez elemental. Necesitamos B, E y A.
    B=1/J*[-1/2, 1/2];
    Ke=w1*B'*E*A*B*J;%matriz elemental
    %Fuerzas nodales equivalentes
    Fe=zeros(2,1);
    for j=1:2 %porque son dos ptos de integración
        int=int2(j);
        w=w2(j); %el pto de integración que estoy evaluando ahora
        N1=0.5*(1-int);%FUNCIÓN DE FORMA PARA ESE PTO DE INTEGRACIÓN
        N2=0.5*(1+int);
        N=[N1,N2];
        %Determinar las fuerzas distribuidas en los nodos
        q1=q0*coords(node1)/L;
        q2=q0*coords(node2)/L;
        q=[q1;q2]; %vector columna
        %Vector fuerza equivalente
        Fe=Fe+w*N'*(N*q)*J;
    end
    %Ensamblaje de matriz de rigidez y fuerza
    %Elemento i, gdl es i+1. ???
    gdl1=[i,i+1]; %grados de libertad correspondientes al elemento i
    K(gdl1,gdl1)=K(gdl1,gdl1)+Ke; %el equivalente a usar la matriz L booleana QUE NO SE USA NUNCA
    %la inicial más el paso anterior, el paso actual, el nodo común
    %Fuerza
    F(gdl1)=F(gdl1)+Fe; %subíndices variables
end
%Reducir las matrices de rigidez y fuerza
%Grados de libertad libres (del nodo 2 al 5 libre, no restringido)
gdl2=2:nnodes;
%Rigidez
KLL=K(gdl2,gdl2);
%Fuerza
F(end)=F(end)+P; %nodo final
FLL=F(gdl2);
%Crear vector de desplazamiento 
u=zeros(ngdl,1);
%Calcular
u(gdl2)=KLL\FLL; %solo cojo los gdl libres uL



