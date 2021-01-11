% QAM y FDMA

clear all;

disp("############################################################")
disp("##                   Bienvenidos                          ##")
disp("##            Proyecto ELO241 - QAM FDMA                  ##")
disp("############################################################")
disp(" ")

while 1

    M = input('En qué tipo de QAM desea modular?(Deben ser potencias de 2): ');
    
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
    
    n_user = input('¿Cuantos usuarios estan conectados?[1-254]: ');
    if n_user <= 0 && n_user >254
        fprintf("\nEl valor no puede ser cero, negativo ni mayor a 254.\n")
        continue
    end
    
    snr = input('Indique la SNR con la cual transmitir[0-inf]: ');
    if n_user <= 0
        fprintf("\nEl valor no puede ser cero, negativo ni mayor a 254.\n")
        continue
    end

    break

end

data = randi([0 1],n_sim*k,n_user); % [range], (matrix) N,M % Binarios

[rows,columns] = size(data);

txSig = [];%Inicializamos vector
fc = [];

for i=1:columns
    txSig= [txSig qammod(data(:,i),M, 'bin','InputType','bit','UnitAveragePower', true)]; % Añadir columnas de los datos de cada usuario
    fc = [fc 10000 + (i-1)*1000]; %Frecuencia central o de portadora
end

txSig = awgn(txSig,snr); % Añadir ruido gaussiano con SNR X

T = 0.005; %Periodo moduladoras

[r, c] = size(fc);

f_muestreo = fc(c)*2*1.2; %Nyquist

x = []; %Inicializamos vector
Q = [];
I = [];

for user=1:(n_user)
    x_mod = []; %Inicializamos vector
    q = [];
    i = [];
    for n=1:(n_sim)

        t = ((n-1)*T:1/(f_muestreo):(n*T)-1/(f_muestreo));

        if isempty(x_mod)
          i = real(txSig(n,user)); % Abscisa
          q = imag(txSig(n,user)); % Ordenada
          x_mod = i*cos(2*pi*fc(user)*t) + q*sin(2*pi*fc(user)*t);
        else
          i = [i, real(txSig(n,user))]; % Abscisa
          q = [q, imag(txSig(n,user))]; % Ordenada
          x_mod = [x_mod, i(n)*cos(2*pi*fc(user)*t) + q(n)*sin(2*pi*fc(user)*t)];
        end
    end
    x = [x transpose(x_mod)]; % Se añade data de un usuario como columna en la matriz
    Q = [Q transpose(q)];
    I = [I transpose(i)];
end
    
t = (0:1/(f_muestreo):(n_sim)*T - 1/(f_muestreo));

x = sum(x,2); % Se juntan todas las señales en una

fft_x = abs(fftshift(fft(x)/length(x))); % FFT de la señal
hz_x = linspace(-f_muestreo/2, f_muestreo/2, length(x));

figure("name", "QAM", 'Position', [200 80 900 600]); % Gráfica señal modulada
subplot(2,1,1);
plot(t, x, 'g');
title('Señal modulada');
xlabel('Tiempo[s]');
ylabel('Magnitud');
legend('x(t)');
grid;
subplot(2,1,2);
plot(hz_x, fft_x, 'r');
title('Espectro señal modulada');
xlabel('Frecuencia[hz]');
ylabel('Magnitud');
legend('x(f)');
grid;

% Demodulación
for user=1:n_user
    % se multiplican por 2 porque el resultado es I/2 y Q/2.
    I_demod = 2.*transpose(x).*cos(2*pi*fc(user)*t);
    Q_demod = 2.*transpose(x).*sin(2*pi*fc(user)*t);

    I_demod = lowpass(I_demod,fc(user),f_muestreo);
    Q_demod = lowpass(Q_demod,fc(user),f_muestreo);

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

    rxSig = transpose(new_I + new_Q*1i);

    data_demod = qamdemod(rxSig,M,'bin','OutputType','bit','UnitAveragePower', true); 

    tf = isequal(data(:,user),data_demod); %,'UnitAveragePower', true

    if tf == 1
        disp('No se han encontrado bits errados en la recepción de los datos')
    else
        disp('Hay bits errados en la recepción de los datos')
    end


    cd{user} = comm.ConstellationDiagram('ShowReferenceConstellation',false);
    cd{user}(rxSig); %Graficar coordenadas de señal recibida

    figure("name", "QAM", 'Position', [200 80 900 600]); % Gráfica señal modulada
    subplot(2,1,1);
    plot(t, I_demod, 'g');
    title('Componente I demodulada');
    xlabel('Tiempo[s]');
    ylabel('Magnitud');
    legend('I(t)');
    grid;
    subplot(2,1,2);
    plot(t, Q_demod, 'b');
    title('Componenete Q demodulada');
    xlabel('Tiempo[s]');
    ylabel('Magnitud');
    legend('Q(t)');
    grid;
end
