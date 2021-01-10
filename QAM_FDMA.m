% QAM y FDMA

clear all;

disp("############################################################")
disp("##                   Bienvenidos                          ##")
disp("##            Proyecto ELO241 - QAM FDMA                  ##")
disp("############################################################")
disp(" ")

while 1

    M = input('En qué tipo de QAM desea modular?(Deben ser potencias de 2): ');
    %M = 256;
    k = log2(M);
    if mod(k,1) ~= 0 && M < 4
        fprintf("\n Ingrese un valor válido entre 8 y 256 que sea potencia de 2.\n")
        continue
    end

    n_sim = input('¿Cuantos datos(bytes) desea enviar?(Generados aleatoriamente)[1-inf]: ');
    if n_sim <= 0
        fprintf("\nEl valor no puede ser cero ni negativo.\n")
        continue
    end

    break

end


 
data = randi([0 1],n_sim*k,1); % [range], (matrix) N,M % Binarios

txSig= qammod(data,M, 'bin','InputType','bit','UnitAveragePower', true);

txSig = awgn(txSig,50); % Añadir ruido gaussiano con SNR 50

Q=1;
T = 0.005; %Periodo moduladoras
fc = 1000; %Frecuencia central o de portadora
f_muestreo = fc*2*1.2; %Nyquist
fase = angle(txSig); % QAM

x = []; %Inicializamos vector
Q = [];
I = [];

for n=1:(n_sim)
    
    t = ((n-1)*T:1/(f_muestreo):(n*T)-1/(f_muestreo));
    %t = (n-1)/(f_muestreo);
    
    if isempty(x)
      I = real(txSig(n));
      Q = imag(txSig(n));
      %A = norm(txSig(n));
      %I = fase(n); % Desfase
      x = I*cos(2*pi*fc*t) + Q*sin(2*pi*fc*t);
    else
      I = [I, real(txSig(n))]; % Desfase
      Q = [Q, imag(txSig(n))]; %Magnitud
      x = [x, I(n)*cos(2*pi*fc*t) + Q(n)*sin(2*pi*fc*t)];
    end
    %data(8*(n-1)+1:8*n)
end

t = (0:1/(f_muestreo):(n_sim)*T - 1/(f_muestreo));
%t = (0:1/(f_muestreo):(n_sim-1)/(f_muestreo)); % Tiempo

figure("name", "QAM", 'Position', [200 80 900 600]); % Gráfica señal modulada
subplot(1,1,1);
plot(t, x, 'g');
title('Señal modulada');
xlabel('Tiempo[s]');
ylabel('Magnitud');
legend('x(t)');
grid;

% se multiplican por 2 porque el resultado es I/2 y Q/2.
I_demod = 2.*x.*cos(2*pi*fc*t);
Q_demod = 2.*x.*sin(2*pi*fc*t);

I_demod = lowpass(I_demod,fc,f_muestreo);
Q_demod = lowpass(Q_demod,fc,f_muestreo);

n_muestras = T*f_muestreo;

new_I = [];
new_Q = [];

for n = 1: n_sim
    if n==1
        new_I = sum(I_demod(n_muestras*(n-1)+1:n_muestras*n))/n_muestras;
        new_Q = sum(Q_demod(n_muestras*(n-1)+1:n_muestras*n))/n_muestras;
    else
        new_I = [new_I, sum(I_demod(n_muestras*(n-1)+1:n_muestras*n))/n_muestras];
        new_Q = [new_Q, sum(Q_demod(n_muestras*(n-1)+1:n_muestras*n))/n_muestras];
    end
end

rxSig = transpose(new_I + new_Q*j);

data_demod = qamdemod(rxSig,M,'bin','OutputType','bit','UnitAveragePower', true); 

tf = isequal(data,data_demod); %,'UnitAveragePower', true

if tf == 1
    disp('La data recibida es exactamente igual a la enviada')
else
    disp('Hay errores en la recepción de los datos')
end


cd = comm.ConstellationDiagram('ShowReferenceConstellation',false);
step(cd,rxSig); %Graficar coordenadas de señal recibida



figure("name", "QAM", 'Position', [200 80 900 600]); % Gráfica señal modulada
subplot(2,1,1);
plot(t, I_demod, 'g');
title('Señal demodulada');
xlabel('Tiempo[s]');
ylabel('Magnitud');
legend('I(t)');
grid;
subplot(2,1,2);
plot(t, Q_demod, 'g');
title('Señal demodulada');
xlabel('Tiempo[s]');
ylabel('Magnitud');
legend('Q(t)');
grid;

