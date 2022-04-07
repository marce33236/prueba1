function[a,ssq,k,A,C,curv]=gn_lm2(a,c0,q,p,maxetapas,tol,strg)
//a, parametros iniciales vector; ka1, ka1*ka2 y ka1*ka2*ka3.
//c0, concentración inicial, escalar.
//q, absorbancias, matriz
//p, pH´s, vector
//a, parametros optimizados;ssq,suma de cuadrados;k,número de iteraciones
exec(strg+'\fres_ampl_tp.sci', -1) //Llamada a la función que calcula ssq
m=length(p)
n=length(a)
//tol=1e-5  podria seleccionarse entre dos criterios de convergncia
mp=0   //parametro de marquardt
delta=1e-6;  //tamaño del paso para difernciación numérica
k=0
ssq_old=1e50
while k<maxetapas
    [r0,C,A]=fres_ampl_tp(a,q,p,c0)  //cálculo d residuos
    ssq=sum(r0.*r0);
    crit_conv=(ssq_old-ssq)/ssq_old;
//    fprintf(1,'it=%d, ssq=%f, mp=%f, crit_conv=%f\n',k,ssq,mp,crit_conv)
    if abs(crit_conv)<=tol then
        if mp==0 then
            break
        else
            mp=0
            r0_old=r0
        end
    elseif crit_conv>tol    //convergencia
        mp=mp/3;
        ssq_old=ssq;
        r0_old=r0
        for i=1:n
            a(i)=(1+delta)*a(i)
            for uk=1:length(r0)
                r=fres_ampl_tp(a,q,p,c0)
                J(uk,i)=(r(uk)-r0(uk))/(delta*a(i));
            end
            a(i)=a(i)/(1+delta);
        end
    elseif crit_conv <-tol
        if mp==0 then
            mp=1;
        else
            mp=mp*5;
        end
        a=a-delta_a
    end
    w=mp*eye(n,n)
    J_mp=[J;w];                   //matriz jacobiana aumentada
    [sa,sb]=size(a)
    r0_mp=[r0_old;zeros(sa,sb)];  //vector de residuos aumentado
    delta_a=-J_mp\r0_mp;              // cálculo del desplazamiento de los parámetros
    a=a+delta_a
    k=k+1
end
curv=J'*J;
endfunction