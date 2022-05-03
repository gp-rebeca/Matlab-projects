%% Datos
L=1; %longitud total
a0=L/20; 
E=100; 
P=1; f0=25; 
%% Desplazamiento en el extremo libre. Solución fuerte 
x=L %en el extremo libre la x=L
a=a0*exp(-x/(2*L));
A=a^2;
f=f0*(x/L)^3;
ua=(exp(x/L)*(L*P + f0*(-(23*L^2)/4 + 6*L*x - 3*x^2 + x^3/L - x^4/(4*L^2))))/(100*a0^2)-(f0*(-23*L^2/4)+L*P)/(100*a0^2)
%% Coordenadas nodos
nelems=3; %número elementos
nnodos=7; %número nodos
%tenemos una barra horizontal origen (0,0), longitud L. 
coords=[0:(L/6):L]; %Coordenadas nodos. De 0 a L, incremento L/número de nodos
%Ptos de integración para la K
%Cuadratura 
int3=[-0.774596697, 0, 0.774596697]; %almacena la función de integración
w3=[0.5555555556,0.8888888889,0.5555555556];
%¿Cuántos grados de libertad tenemos si tenemos 7 nodos? Cada nodo tiene 1 gdl. 
%(Tiene dos pero solo uno es de deformación, el axil. El otro es de sólido
%rígido).
ngdl=nnodos*1;
%Matriz de rigidez total
K_3=zeros(ngdl); %7x7
%La fuerza nodal total
F_3=zeros(ngdl,1); %7x1

%Bucles para generar y ensamblar la matriz de rigidez
for i=1:nelems %si i es 1, tenemos 1 y 2
    nodo1=i*2-1; nodo2=i*2; nodo3=i*2+1;
    %longitud de ese elemento i
    dl=coords(nodo3)-coords(nodo1);
    %Jacobiano-para transformar de x a xi
    J=dl/2; %longitud elemento entre 2
    Ke=zeros(3);
    
     for h=1:3
        B=zeros(1,3);
        int=int3(h)
        B(1)=1/J*[int-0.5];
        B(2)=1/J*[-2*int];
        B(3)=1/J*[int+0.5];
        a1=a0*exp(-coords(nodo1)/(2*L));
        a2=a0*exp(-coords(nodo2)/(2*L));
        a3=a0*exp(-coords(nodo3)/(2*L));
        A_=[a1^2,a2^2,a3^2];
        
        N1=0.5*int*(int-1);%FUNCIÓN DE FORMA PARA ESE PTO DE INTEGRACIÓN 
        N2=1-int^2;
        N3=0.5*int*(int+1);
        N=[N1,N2,N3];
        A=A_*N'
        w=w3(h);
        Ke=w*B'*E*A*B*J+Ke %matriz elemental 3x3
        
     end
       Fe=zeros(3,1);
     for j=1:3
        xi=int3(j);
        w_=w3(j); %el pto de integración que estoy evaluando ahora
    %Fuerzas nodales equivalentes
        N1_=0.5*xi*(xi-1);%FUNCIÓN DE FORMA PARA ESE PTO DE INTEGRACIÓN 
        N2_=1-xi^2;
        N3_=0.5*xi*(xi+1);
        N_=[N1_,N2_,N3_];
        %Determinar las fuerzas distribuidas en los nodos
        f1=f0*(coords(nodo1)/L)^3;
        f2=f0*(coords(nodo2)/L)^3;
        f3=f0*(coords(nodo3)/L)^3;
        ftot=[f1;f2;f3]; %vector columna
        %Vector fuerza equivalente
        Fe=Fe+w_*N_'*(N_*ftot)*J; %es un sumatorio, vas guardando la suma anterior en Fe
    end
    %Ensamblaje de matriz de rigidez y fuerza
    gdl1=[i*2-1,i*2,i*2+1]; %grados de libertad correspondientes al elemento i
    K_3(gdl1,gdl1)=K_3(gdl1,gdl1)+Ke %el equivalente a usar la matriz L booleana QUE NO SE USA NUNCA
    %la inicial más el paso anterior, el paso actual, el nodo común
    %Fuerza
    F_3(gdl1)=F_3(gdl1)+Fe %subíndices variables
end
%% Reducir las matrices de rigidez y fuerza
%Grados de libertad libres (del nodo 2 al 7 libre, no restringido)
gdl2=2:nnodos;
%Rigidez

KLL=K_3(gdl2,gdl2);
%Fuerza
F_=zeros(7,1);
F_(end)=F_3(end)+P; %nodo final
F_(1:6)=F_3(1:6);
FLL=F_(gdl2);
%Crear vector de desplazamiento 
u=zeros(ngdl,1);
%Calcular
u(gdl2)=KLL\FLL %solo cojo los gdl libres uL
%% Sigma
Bxi1=zeros(1,3);
Bxi1(1)=1/J*[int3(1)-0.5];
Bxi1(2)=1/J*[-2*int3(1)];
Bxi1(3)=1/J*[int3(1)+0.5];
Bxi2=zeros(1,3);
Bxi2(1)=1/J*[int3(2)-0.5];
Bxi2(2)=1/J*[-2*int3(2)];
Bxi2(3)=1/J*[int3(2)+0.5];
Bxi3=zeros(1,3);
Bxi3(1)=1/J*[int3(3)-0.5];
Bxi3(2)=1/J*[-2*int3(3)];
Bxi3(3)=1/J*[int3(3)+0.5];    
eps11=Bxi1*u(1:3);
eps12=Bxi1*u(3:5);
eps13=Bxi1*u(5:7);
eps21=Bxi2*u(1:3);
eps22=Bxi2*u(3:5);
eps23=Bxi2*u(5:7);
eps31=Bxi3*u(1:3);
eps32=Bxi3*u(3:5);
eps33=Bxi3*u(5:7);
eps=[eps11,eps12,eps13,eps21,eps22,eps23,eps31,eps32,eps33];
epsmax=max(eps)
sigma_max=E*epsmax
%% Número de elementos necesario

