//Introducción de datos
//etiq_1=['Introduca el directorio de trabajo:']
etiq_1=['Introduca el directorio de trabajo:']
tit_1=['Directorio de trabajo:';' Es el directorio en el que estan contenidos los archivos necesrios para el cálculo']
direct=x_mdialog(tit_1,etiq_1,['D:\Documentos\Química\Mi carpeta de trabajo'])
messagebox("El siguiente paso requiere el nombre de los archivos de espectros realizados. Estos deben nombrarse previamente en serie. Ej: Mis archivos son: esp_nar1,esp_nar2,...,esp_nar77,.... Por lo tanto, en Nombre de archivo coloco esp_nar. No coloco la numeración!.', "modal", "info", ["Ok"])
etiq3=['Introduzca el nombre de archivo:','Introduzca la extension de los archivos:']
tit3=['Archivos:';'Datos experimentales de espectros.']
arch=x_mdialog(tit3,etiq3,['esp_nar','.asc'])
etiq2=['Número de espectros realizados:','Longitud de onda inicial:','Longitud de onda final:','Número de puntos o longitudes de onda']
tit2='Introduzca los siguientes datos'
vector=x_mdialog(tit2,etiq2,['0','0','0','0'])
ner=eval(vector(1))
np=eval(vector(4))
long=linspace(eval(vector(2)),eval(vector(3)),np)
B=mopen(direct+'\long_de_onda.txt','wt');
nph=x_matrix('Introduzca los valores de pH:',zeros(ner,1))
C=mopen(direct+'\pH.txt','wt');
absor=[]
for i=1:ner
    dirt(i)=direct+"\"+arch(1)+string(i)+arch(2)
    ur=mopen(dirt(i),'rt')
    u=mgetl(ur,-1)
    c=part(u(13),[11:13])
    for j=1:np
        absor(i,j)=eval(u(j+15))
        mfprintf(B,'%f\n',long(j))
    end
    mfprintf(C,'%f\n',nph(i))
    //Estaba intentando armar un archivo que muestre los datos que se utilizan.
   //mfprintf(A,'###########\n')
    //mfprintf(A,'## Datos ##\n')
    //mfprintf(A,'###########\n')
    //mfprintf(A,'Abs.(UA)   long. de onda(nm)\n')
    mclose(ur)
end
mclose(B)
mclose(C)
fprintfMat(direct+'\absorbancias.txt',absor,'%f')

//parametros iniciales
k=x_matrix('Introduzca los valores iniciales de los parámetros:',[1e-5;1e-10;1e-15])

//Definición de datos necesarios para el cálculo
etiq4=['Concentración inicial:','Máximo de etapas:','Criterio de convergencia:']
tit4='Introduzca los siguientes datos para el cálculo.'
vector_2=x_mdialog(tit2,etiq2,['3.5E-5','50','1E-5'])

//
//Cálculo de las constantes de equilibrio.
conc=eval(vector_2(1))
maxiter=eval(vector_2(2))
crit=eval(vector_2(3))
PH=fscanfMat(direct+'\pH.txt');//pH experimentales
ABS=fscanfMat(direct+'\absorbancias.txt');//matriz de absorbancias
exec(direct+'\gn_lm2.sci', -1)
[k,ssq,it,A,C,curv]=gn_lm2(k,conc,ABS,PH,maxiter,crit,direct)
[abs1,abs2]=size(ABS)
sig_r=sqrt(ssq/(abs1*abs2-length(k)))
sig_k=sig_r*sqrt(diag(inv(curv)))
//Tengo q corregir esta parte para que muestre un arcchivo tipo log de salida
//for i=1:length(k)
//    fprintf(1,'k(%d): %.0f +- %f\n',i,k(i),sig_k(i))
//end
//fprintf(1,'sig_r: %f\n',sig_r)


//Resulados
ka1=k(1)
pka1=-log10(ka1)
ka2=k(2)/k(1)
pka2=-log10(ka2)
ka3=k(3)/k(2)
pka3=-log10(ka3)
disp(pka3,'pka3',pka2,'pka2',pka1,'pka1','Constantes de acidez',it,'Número de iteraciones:')
