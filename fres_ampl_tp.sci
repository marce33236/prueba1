function [r,C,A]=fres_ampl_tp(a,q,p,c)
//Función que calcula los residuos, aplicable a ácidos tripróticos.
//a vector , parametros iniciales.
//q matriz de absorbancias; filas, cant. de med de pH; col., cant. de long. de onda medidas.
//p vector fila, pH.
// c, concentración inicial
m=length(p);
for i=1:m //cant. de med de pH.
     h(i)=10^(-p(i))
     C(i,1)=(c*h(i)^3)/(h(i)^3+a(1)*h(i)^2+a(2)*h(i)+a(3));  //[FH3]
     C(i,2)=(c*h(i)^2*a(1))/(h(i)^3+a(1)*h(i)^2+a(2)*h(i)+a(3)); //[FH2-]
     C(i,3)=(c*h(i)*a(1)*a(2))/(h(i)^3+a(1)*h(i)^2+a(2)*h(i)+a(3));   //[FH-2]
     C(i,4)=(c*a(3))/(h(i)^3+a(1)*h(i)^2+a(2)*h(i)+a(3));   //[F-3]
end
A=C\q; // matriz de absortividades molares; filas, especie qca.;col., long. de onda
R=q-C*A
r=R(:)  //vectorización de R
endfunction