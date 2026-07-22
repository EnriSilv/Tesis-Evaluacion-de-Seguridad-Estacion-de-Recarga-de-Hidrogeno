
%-Montecarlo
pi=3*(10^6); %-Presión mínima de operación del hidrógeno (Pa).
po=70*(10^6); %-Presión máxima de operación del hidrógeno (Pa).
patm=101325; %-Presión atmosférica (Pa).
xi=0.04; %-Fracción molar mínima de hidrógeno consumido.
xs=1; %-Fracción molar máxima de hidrógeno consumido.
ne=12000; %-repeticiones.
st=zeros(50,50); %-matriz donde guardar los valores de velocidad de quema turbulenta (m/s).
xh=zeros(50,50);%-matriz donde guardar los valores de fracción molar del hidrógeno.
ph=zeros(50,50);%-matriz donde guardar los valores de presión del hidrógeno (MPa).
u1=zeros(50,50);%-matriz donde guardar los valores de intensidad de la turbulencia (m/s).
u2=zeros(50,50);%-matriz donde guardar los valores de turbulencia generada por la llama (m/s).
su=zeros(50,50);%-matriz donde guardar los valores de velocidad de quema laminar (m/s).
x1=0.04; %-punto mínimo de inflamación. 
x2=0.75; %-punto máximo de inflamación.
xm1=ones(50,50); %-matriz donde generar el punto mínimo de inflamación.
xm2=ones(50,50); %-matriz donde generar el punto máximo de inflamación.
ym1=zeros(50,50); %-matriz donde guardar los valores de presión de hidrógeno del punto mínimo de inflamación.
ym2=zeros(50,50); %-matriz donde guardar los valores de presión de hidrógeno del punto mínimo de inflamación.
zm1=zeros(50,50); %-matriz donde guardar los valores de velocidad de combustión del punto mínimo de inflamación.
zm2=zeros(50,50); %-matriz donde guardar los valores de velocidad de combustión del punto mínimo de inflamación.
i=0; j=0; %-filas (i) y columnas (j).
tu1=0; tu2=0; tsu=0; tvrms=0; tph=0; damb=0; dcomb=0; tz1=0; tz2=0;
mH=0.002016; %-masa molar del hidrógeno (kg/mol).
mO=0.0032; %-masa molar del oxígeno (kg/mol).
mN=0.02801; %-masa molar del nitrógeno (kg/mol).
mvapor=0.01802; %-masa molar del vapor de agua (kg/mol).
Rg=8.31446; %-constante ideal de los gases (J/mol*K).
fc=0.1; %-factor de compresibilidad gas real (Z).
cont=0; %-contador total.
cont1=0; %-contador datos menores al punto mínimo de inflamación.
cont2=0; %-contador datos mayores al punto máximo de inflamación.
Tamb=298.15; %-temperatura ambiental (K).
Tllama=2400; %-temperatura que alcanza la llama (K).
for i=1:50 
    for j=1:50 
        for n=1:ne
            txh=normrnd(xi,xs); %-generación aleatoria no ordenada de valores de fracción molar de hidrógeno.
            tph=normrnd(pi,po); %-generación aleatoria no ordenada de valores de presión de hidrógeno.
            damb=(tph*(txh*mH+(1-txh)*0.21*mO+(1-txh)*0.79*mN))/(fc*Rg*Tamb); %-densidad del hidrógeno a temperatura ambiente.
            tvrms=sqrt((3*tph)/(fc*damb)); %-velocidad rms del hidrógeno.
            dcomb=(tph*(txh*mH+(79/21)*txh*mN+(1/2)*txh*mO+txh*mvapor))/(fc*Rg*Tllama);
            vard=damb/dcomb; %-factor de expansión.
            tsu=(15*tu1/2); %-velocidad laminar de combustión.
            tu1=tvrms/sqrt(2); %-intensidad de turbulencia.
            tu2=tsu*(1-exp(-tu1/tsu))*((vard-1)/sqrt(3)); %-turbulencia generada por la llama
            tst=tsu+tu1+tu2; %-velocidad turbulenta de combustión.
            if (tph>patm)
               cont=cont+1; 
               ph(i,j)=tph/(10^6); %-matriz presión de hidrógeno.
               if ((txh>=xi)&&(txh<=xs)) 
                   xh(i,j)=txh; %-matriz fraccion molar.
                   if (tu1>0)
                       u1(i,j)=tu1; %-matriz intensidad turbulencia.
                       if (tsu>0)
                           su(i,j)=tsu; %-matriz velocidad laminar de quema
                           if (tu2>0)
                               u2(i,j)=tu2; %-matriz turbulencia generada por llama
                               st(i,j)=tst; %-matriz velocidad turbulenta de quema
                           end
                       end
                   end
               end
           end
           if (txh<=x1)
               ym1(i,j)=tph/(10^6); %-matriz valores de relacion molar del punto mínimo de inflamación.
           end
           if (txh>=x2)
               ym2(i,j)=tph/(10^6); %-matriz valores de relacion molar del punto máximo de inflamación.
           end
        end
        if ((txh<=x1)&&(tph>patm))
            cont1=(cont1)+1;
            tz1=normrnd(0,100000);
            zm1(i,j)=tz1; %-matriz valores de velocidad de combustión del punto mínimo de inflamación
        end
        if ((txh>=x2)&&(tph>patm))
            cont2=(cont2)+1;
            tz2=normrnd(0,100000);
            zm2(i,j)=tz2; %-matriz valores de velocidad de combustión del punto máximo de inflamación
        end
    end
end
prob=((cont-(cont1+cont2))*100)/cont; %-probabilidad de combustión
xm1=(ones(50,50))*x1; %-matriz punto mínimo de inflamación
xm2=(ones(50,50))*x2; %-matriz punto máximo de inflamación
scatter3(xh,ph,st,'filled'); %-gráfica principal
hold on
surf(xm1,ym1,zm1) %-Hiperplano punto mínimo de inflamación
surf(xm2,ym2,zm2) %-Hiperplano punto máximo de inflamación


