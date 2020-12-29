% QAM y FDMA

M = 256;
k = log2(M);

data = randi([0 1],1000*k,1); % [range], (matrix) N,M

data; % binarios

txSig = qammod(data,M,'InputType','bit','UnitAveragePower',true);

rxSig = awgn(txSig,30); % Añadir ruido gaussiano con SNR 30

cd = comm.ConstellationDiagram('ShowReferenceConstellation',false);
step(cd,rxSig); %Graficar coordenadas de señal recibida

A=1;
fc = 1000; %Frecuencia central o de portadora
f_muestreo = fc*2*1.2; %Nyquist
fase = angle(txSig); % QAM
%txSig(1)

x = []; %Inicializamos vector

for n=1:(1000)
    
    t = linspace((n-1)*8/(f_muestreo),(8*n-1)/(f_muestreo),8);
    A = norm(txSig(n)); %Magnitud
    phi = fase(n); % Desfase
    if isempty(x)
      x = A*cos(2*pi*fc*t + phi); % Asignamos primer valor y dimensiones del vector
    else
      x = [x, A*cos(2*pi*fc*t + phi)]; % En caso de no estar vacío de agregan los siguientes valores
    end
    %data(8*(n-1)+1:8*n)
end

t = (0:1/(f_muestreo):(8000-1)/(f_muestreo)); % Tiempo

%size(t)
%size(x)
%disp("hola")

figure("name", "QAM", 'Position', [200 80 900 600]); % Gráfica señal modulada
subplot(1,1,1);
plot(t, x, 'g');
title('Señal modulada');
xlabel('Tiempo[s]');
ylabel('Magnitud');
legend('x(t)');
grid;


