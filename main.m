clc; clear; close all;
warning('off', 'all');

n_arm =350;
fmax=1.1468;% frecuencia maxima en la señal
f_N_max=n_arm*fmax;
T_N_max=1/f_N_max;

stepsize=T_N_max/20

tsim=8;
%% Pre procesamiento de la señal ECG
Fs=500; G=200; load('data/rec_1.mat') %sana --------------------------- 0
%Fs = 360;  G = 200; load('data_arritmia/ar1.mat') %arritmia ------------ 1

ECGsignal= val/G;
ecg=(ECGsignal-mean(ECGsignal))/std(ECGsignal);

Ts=1/Fs;
t=(1:1:length(ECGsignal))*Ts;
ecg =[0 ecg]; 
t= [0 t]; 

plot(t,ecg)
title('Señal ECG')
xlabel('Tiempo(seg)') 
ylabel('Amplitud(mV)') 

z=60/200;
umbral_y=6*mean(abs(ecg)); %%sana
 %umbral_y=4*mean(abs(ecg)); %%arritmia


umbral_x= z*Fs;
[pks,locs]=findpeaks(ecg,'MinPeakHeight',umbral_y,'MinPeakDistance',umbral_x);
hold on;
scatter(t(locs),pks);
num_ciclos = length(locs) - 1;
max_length_cycle = max(diff(locs));
max_length_spectrum = floor(max_length_cycle / 2) + 1; 
v_periodos = zeros(1, num_ciclos);
v_ecg_ciclos = cell(1, num_ciclos);
v_t = cell(1, num_ciclos);

%% Identificar cada inicio de ciclo y fin para posteriormente guardarlos 
 for j = 1:num_ciclos
    HRV = diff(t(locs(j:j+1))); 
    periodo = HRV(1);
    v_periodos(j)=periodo;
    period = locs(j+1) - locs(j);
    inicio_ciclo = t(locs(j));
    fin_ciclo = t(locs(j)) + periodo;
    ciclo = ecg(locs(j):locs(j) + period - 1);
    v_ecg_ciclos{j} = ciclo;
    v_t{j}=t(1:length(ciclo));
%  figure();
%  plot(tn,ciclo)
end
%figure();
%plot(v_ecg_ciclos)

%% Esta parte unicamente es para concatenar los ciclos y comprobar que se tiene la señal original 
ciclo_concatenado=[];
tn_concatenado=[];
for j = 1:num_ciclos
    ciclo = v_ecg_ciclos{j};
   % figure(); plot(ciclo) 
    tn = t(locs(j):locs(j) + length(ciclo) - 1 );
    ciclo_concatenado = [ciclo_concatenado, ciclo];
    tn_concatenado = [tn_concatenado, tn];
end
figure;
plot(tn_concatenado, ciclo_concatenado);

%% Procesamiento para mandar la información a Simulink
f_new=1./v_periodos; %% Frecuencia de cada ciclo
for j = 1:num_ciclos
%for j = 1:1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ecg_ciclo =v_ecg_ciclos{j}; 
   % tnn=  v_t{j}; %% tiempo
    ecg_nuevo=[ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo ecg_ciclo];
    tnuevo=0:length(ecg_nuevo)-1;
    tnuevo=tnuevo/Fs; %%%%%%%%%%%%%%%%%%%
    %tnuevo=[0 tnuevo]; % Se concatena un 0 debido a que al mandar los datos a Simulink se elimina el primer dato
    %ecg_nuevo=[0 ecg_nuevo]; % Descomentar estas dos lineas para cuando sean nuevas señales
    tnuevo=[tnuevo]; % Aqui no se concatena el 0 porque ya se hizo anteriormente cuando se guardó cada ciclo
    ecg_nuevo=[ecg_nuevo]; % 
    %figure();
    %plot(tnuevo(1:length(ecg_nuevo)),ecg_nuevo)
    %xlabel('Tiempo(seg)') 
    %ylabel('Amplitud(mV)') 

    data.time = tnuevo;
    data.signals.values = ecg_nuevo';
    data.signals.dimensions = 1;
    assignin('base', 'data', data);

    %% FOURIER
    if j == 1
    L = length(ecg_ciclo);
    Y = fft(ecg_ciclo);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    P1=P1(1:68);
    HarmonicF1=P1(1);
    HarmonicF2=P1(2);
    HarmonicF3=P1(3);
    HarmonicF4=P1(4);
    f=f(1:68);
    f_graficar=f;
    figure()
    bar(f,P1,'g');
    title('Espectro de Fourier con comando');
    end
    %% Definicion y calculos%%
    %n_arm =320;
    f=f_new(j);
    w=2*pi*f;
    matriz = cell(1, n_arm);

    for i = 1:n_arm
        matriz{i} = i * [0, w; -w, 0];
    end
    
    A = blkdiag(matriz{:});
    C = zeros(1,n_arm*2);

    for i = 1:2:n_arm*2
        eval(['k' num2str(i) ' = 0.5;']);
        eval(['k' num2str(i+1) ' = 0.01;']);
        C(i) = 0.5;
         C(i+1) = 0.01;

    end

    Q =eye(2*n_arm)*10;
    R = 0.001;
    [P, L, G] = care(A', C', Q, R);
    L=G;
    x0=zeros(n_arm*2,1);

   % %{
     sim("armonicosv2"); %% Aqui se llama al identificador de armonicos de Simulink
     Aedc = Arm_e_DC(end,:);
     Aedc = abs( Aedc );
     Arm_e= Arm_e(end,:);
     Arm_e_tot = [0 Aedc Arm_e];
    %Arm_e_tot= [0 P1 0];
    
    %% Clasificador
    sim("clasificador_armonicos");
    prediccion = prediccion(1);
    
    if prediccion == 0
        fprintf('Resultado %d: Sano\n',j);
    elseif prediccion == 1
       fprintf('Resultado %d: Arritmia\n',j);
    else
       fprintf('Resultado %d: Error',j);
    end
 %%}
    
    
end