%Elemento viga
E=1; A=1; I=1; L=1; x1=0; x2=L; J=L/2; %(IGUAL SIEMPRE, SEAN LOS NODOS QUE SEAN)
xi1=-1/sqrt(3); w1=1;
B1=zeros(1,4); %matriz diapositiva 58
B1(1)=(6*xi1)/(x1-x2)^2;
B1(2)=-(3*xi1-1)/(x1-x2);
B1(3)=-(6*xi1)/(x1-x2)^2;
B1(4)=-(3*xi1+1)/(x1-x2);
%Necesitamos 2 ptos de integración. El primero nos da K1
K1=w1*B1'*E*I*B1*J
%Segundo pto de integración
xi2=1/sqrt(3); w2=1;
B2=zeros(1,4);
B2(1)=(6*xi2)/(x1-x2)^2;
B2(2)=-(3*xi2-1)/(x1-x2);
B2(3)=-(6*xi2)/(x1-x2)^2;
B2(4)=-(3*xi2+1)/(x1-x2);
K2=w2*B2'*E*I*B2*J
%PREGUNTA DE EXAMEN: Calcula el término 1,1 de la matriz de rigidez de un elemento viga. (B1*B1).
K=K1+K2