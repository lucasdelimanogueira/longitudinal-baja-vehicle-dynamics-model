%Rotina para simulacao da dinâmica longitudinal do veículo. 

%parametros constantes do veiculo e ambiente
m = 220;%massa do veiculo + piloto
zfix = 7.85; %relacao de transmissão do redutor fixo
rp = 0.2921; %raio do pneu 
ef = 0.865; %eficiencia da cvt nao contabilizada pelo coastdown
h = 0.525; %altura do CG em relacao ao solo 
lt= 0.583; %distancia do CG ao eixo traseiro 
ld = 0.806; %distancia do CG ao eixo dianteiro
tetagraus = 10;%angulo da rampa (graus)
teta = tetagraus*(pi/180);
isaida = 1.9353; %inercia rotacional saida do redutor fixo 
                 % 4 rodas/discos de freio/semieixos/
isaida = isaida/2; %dividindo pois o total é referente às 4 rodas

imotora = 0.0025; %inercia rotacional polia motora
ientrada = 0.0055; %inercia rotacional polia movida/entrada do redutor fixo

%configuracoes da rotina
rotrefparainterp = [1620:3.375:4100];  
rotref = [1620:0.2:4100]; %vetor de referencia das possiveis rotacoes
t1 = 0; %tempo inicial 
tn = 15; %tempo final
dt = 0.01; %incremento de tempo
tempo = t1:dt:tn; %vetor tempo

%relacao de transmissão da CVT
zcvtparainterp = cvt';
zcvt = interp1(rotrefparainterp,zcvtparainterp2,rotref,'spline');

rotrodaref = rotref./(zcvt.*zfix);


%inicializacao de matrizes dos plots
vveiculoplot = zeros(length(tempo),1);
sveiculoplot = zeros(length(tempo),1);
aveiculoplot = zeros(length(tempo),1);
rotmotorplot = zeros(length(tempo),1);
rotrodaplot = zeros(length(tempo),1);
slipplot = zeros(length(tempo),1);


Tmotor = interp1(rotrefparainterp,Torque,rotref,'spline'); %torque do motor


%parametros iniciais do veiculo
aveiculo = 0; %aceleracao do veiculo (inicial)
vveiculoini = 0; %velocidade inicial do veiculo
sveiculoini = 0; %espaço inicial
rotmotor = 2200; %rotacao do motor (RPM)

%forças normais nas rodas Dianteira e Traseira
FzD = (-m*aveiculo*h-m*9.81*h*sin(teta)+m*9.81*lt*cos(teta))/(lt+ld);
FzT = (m*aveiculo*h+m*9.81*h*sin(teta)+m*9.81*ld*cos(teta))/(lt+ld);

%força trativa inicial
Fx = 0;

%rotacao  inicial da roda traseira
rotrodaini = (rotmotor/(zfix*interp1(rotref,zcvt,rotmotor,'nearest')))*((2*pi)/60);%rad/s

%inicio do laço de repetição
for t = 1:length(tempo)
 
  %Motor & transmissão
    if rotmotor>4100 %caso rotacao do motor >4100 torque é zero
        
        %torque na roda traseira
        Troda = 0;
        %inercia rotacional equivalente no eixo da motora 
        Imoteq = imotora*(zcvt(length(zcvt)).*zfix)^2;
        
    elseif rotmotor<1620 %caso rotacao do motor <1620 torque é zero
        
        %torque na roda traseira
        Troda = 0;
        %inercia rotacional equivalente no eixo da motora 
        Imoteq = imotora*(zcvt(1)*zfix)^2;
        
    else %caso rotacao dentro da faixa normal, calcula torque na roda
        
        %torque na roda traseira
        Troda = interp1(rotref,Tmotor,rotmotor,'nearest')*zfix*interp1(rotref,zcvt,rotmotor,'nearest'); %Nm
        %inercia rotacional equivalente no eixo da motora 
        Imoteq = imotora*(interp1(rotref,zcvt,rotmotor,'nearest')*zfix)^2;
    end
    
    %inercia rotacional no eixo de entrada do redutor fixo em relação à ao
    Ientradaeq = ientrada*((zfix)^2); 
    Isaidaeq = isaida; %rodas traseiras+semieixos+eixo de saída do redutor+disco de freio etc
    Ieq = Imoteq+Ientradaeq+Isaidaeq; %inercia rotacional total equivalente

    
    
  %sistema roda
    %torque resistivo (resistência a rolagem + fricção) na roda traseira
    Tres = rp*(m*0.0029.*vveiculoini.*3.6 + ((FzT*0.41)/9.81)); %Nm 
    
     if Tres < 0
        Tres = 0;
     end
    
    %calculo da rotacao da roda
    rotrodafin = (1/Ieq)*(Troda*ef-Fx*rp-Tres)*dt+rotrodaini; %rad/s
    rotrodaini = rotrodafin;   
    
    
  %slip ratio
    if vveiculoini < 0 %se velocidade < 0 zera slip ratio
        slip = 0;
    else
        slip = (rotrodafin*rp)/vveiculoini - 1; %slip ratio
    end
  
    
  %força trativa  
    %coeficiente de adesão
    
    %formula mágica de pacejka
    %x = slip+shx
    %yx = (D*sin(C*atan(B*x-E*(B*x-atan(B*x)))));
    %tracf = yx + svx
    %mi = tracf/Fz
    
    %dados obtidos pelo artigo
    D = 1.2377*(0.640353+0.261665*exp(-0.080955*slip*3.6*(vveiculoini+0.1)))*1177.2; %asfalto seco
    %D = 1.2377*(0.590189-0.185632*exp(-0.192696*slip*3.6*(vveiculoini+0.1)))*1177.2; %terra batida seca
    x = slip-0.0037;
    mi = ((D*sin(1.3971*atan(20.0573*(x)-(-1.3078)*(20.0573*(x)-atan(20.0573*(x)))))) + 92.5868)/1177.2;
   
    if abs(Troda-Tres) < mi*FzT*rp
        Fx = abs((Troda-Tres)/rp); %força trativa
    else
        Fx = mi*FzT; %força trativa
    end
       
  %veiculo
    if vveiculoini < 0
        Faero = m*(2.38*10^(-4).*(3.6.*vveiculoini).^2); %força de arrasto
        Ffricloss = m*0.0029.*vveiculoini.*3.6; %perdas na transmissão por fricção
        Frollres = 0.41*((FzD+FzT)/9.81); %resistencia a rolagem
        Fres = Faero+Ffricloss+Frollres; %força resistiva
        Fres = -Fres; %muda sentido pois velocidade < 0
    else
        Faero = m*(2.38*10^(-4).*(3.6.*vveiculoini).^2); %força de arrasto
        Ffricloss = m*0.0029.*vveiculoini.*3.6; %perdas na transmissão por fricção
        Frollres = 0.41*((FzD+FzT)/9.81); %resistencia a rolagem
        Fres = Faero+Ffricloss+Frollres; %força resistiva 
    end
        aveiculo = (Fx-Fres-m*9.81*sin(teta))/m; %aceleração do veículo
        
        %velocidade
        vveiculo = vveiculoini + aveiculo*dt;
        vveiculoini = vveiculo;
        
        %espaço percorrido
        sveiculo = sveiculoini + vveiculo*dt;
        sveiculoini = sveiculo;
        
        
  %transferencia de peso
    %Traseira 
    FzT = (m*aveiculo*h+m*9.81*h*sin(teta)+m*9.81*ld*cos(teta))/(lt+ld);
    %Dianteira
    FzD = (-m*aveiculo*h-m*9.81*h*sin(teta)+m*9.81*lt*cos(teta))/(lt+ld);
    
    
  %atualizacao da rotacao do motor
    rotrodafinRPM = rotrodafin*(60/(2*pi)); %rotacao da roda em RPM

    if rotrodafinRPM>rotmotor/(zcvt(length(zcvt))*zfix) %limitação da relação CVT máxima
        %rotação do motor
        rotmotor = rotrodafinRPM*(zcvt(length(zcvt))*zfix);
        
    elseif rotrodaTfinRPM<rotmotor/(zcvt(1)*zfix) %limitação da relação CVT 
        rotmotor = rotrodafinRPM*(zcvt(1)*zfix);
        
    else
    %rotação do motor
    rotmotor = interp1(rotrodaref,rotref,rotrodafinRPM,'nearest'); 
    
    end
 
    
  %armazenar dados para plot
    vveiculoplot(t,:) = vveiculo;
    sveiculoplot(t,:) = sveiculo;
    aveiculoplot(t,:) = aveiculo;
    rotmotorplot(t,:) = rotmotor;
    slipplot(t,:) = slip;
    
    
end

%calculo do Tempo em 30 m e Velocidade em 100 m
tempo30mplot = max(tempo((sveiculoplot<=30)))
velocidade100m = 3.6*min(vveiculoplot(sveiculoplot>=100))


%plot dos gráficos
  subplot(2,2,1);plot(temporeal,vveiculoreal,tempo,3.6*vveiculoplot);title('Velocidade do veículo');ylabel('Velocidade (km/h)');xlabel('Tempo (s)');legend({'Real','Teórico'},'Location','southeast');xlim([0 10])
  subplot(2,2,2);plot(temporeal,sveiculoreal,tempo,sveiculoplot);title('Distância percorrida pelo veículo');ylabel('Distância (m)');xlabel('Tempo (s)');legend({'Real','Teórico'},'Location','southeast');xlim([0 10])
  % subplot(1,2,3);plot(tempo,slipplot);;title('Deslizamento ao longo do tempo');ylabel('Slip ratio');xlabel('Tempo (s)');
 
%  subplot(3,2,4);plot(tempo,aveiculoplot);title('Aceleração do veículo');ylabel('Aceleraçao (m/s²)');xlabel('Tempo (s)');
%  
% 
%     
