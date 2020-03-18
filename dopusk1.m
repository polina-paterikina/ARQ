clear all;
close all;
m = 8; %длина информационной последовательности
crc = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1]; %crc16 m=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
% crc = [1 0 1 0 1 0 1 1 1]; %crc8 m=4 x8+x7+x6+x4+x2+1
% crc = [1 1 1 0 1]; % x4+x2+x+1 m=3
% crc = [1 1 0 1]; % x3+x+1 m=4
r = length(crc)-1;
n = 250000; %количество экспериментов
snr_arr_pract = -30:5:10; %сигнал/шум в дБ для экспериментальных вычислений
snr_arr_theor = -30:10; %сигнал/шум в дБ для теоретических вычислений


codewords = zeros(2^m, m+r); %массив кодовых слов
mod_codewords = zeros(2^m, m+r); %массив промодулированных кодовых слов
d_arr = zeros((2^m), 1); %массив весов всех кодовых слов
for word = 0:2^m-1 % цикл по всем информационным последовательностям
    mxr = bitshift(word, r); % m(x)*x^r
    [~, c] = gfdeconv(de2bi(mxr), crc); % контрольная сумма c(x) = m(x)*x^r mod(g(x))
    a = de2bi(bitxor(mxr, bi2de(c)), m+r); % кодовое слово a(x) = m(x)*x^r + c(x)
    codewords(word+1, :) = a;
    d_arr(word+1) = sum(a); % подсчёт веса кодового слова
    a_m = a.*2-1; % модуляция кодового слова 1->+1, 0->-1
    mod_codewords(word+1, :) = a_m;
end

Pe_decode_Theor = ones(1,length(snr_arr_theor)).*(1/2^r); % верхняя граница ошибки декодирования CRC
Pb_theor = qfunc(sqrt(2*(10.^(snr_arr_theor./10)))); % теоретическая вероятность ошибки на бит для BPSK

d_min = min(d_arr(2:end)); % нахождение минимального расстояния кода
P_ed = zeros(1, length(Pb_theor)); % точная вероятность ошибки декодирования CRC
index = 1;
for p = Pb_theor % цикл по теоретической вероятности ошибки на бит для BPSK
    tmp = 0;
    for i = d_min:(m+r) % цикл подсчёта суммы по весам кодовых слов от минимального расстояния кода d до длины кодового слова n
        Ai = sum(d_arr == i); % количество кодовых слов весом i
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i)); % подсчёт суммы согласно формуле
    end
    P_ed(index) = tmp;
    index = index + 1;
end
index = 1;

Pb = zeros(1, length(snr_arr_pract));
Pe_decode = zeros(1, length(snr_arr_pract));
for SNRdB = snr_arr_pract % цикл по массиву сигнал/шум в дБ для экспериментальных вычислений
    N_b = 0; % количество ошибочных бит
    Ne_decode = 0; % количество ошибок декодирования CRC
    SNR = 10^(SNRdB/10); % отношение сигнал/шум
    sigma = sqrt(1/(2*SNR));
    parfor i = 1:n % цикл делающий заданное количество экспериментов
        rand_mess_ind = randi(2^m, 1); % генерация случайного индекса сообщения
        a = codewords(rand_mess_ind, :); % кодовое слово
        a_m = mod_codewords(rand_mess_ind, :); % промодулированное кодовое слово
        b_m = a_m +sigma*randn(1, length(a_m)); % добавление шума в канале
        b = b_m > 0; % демодуляция принятого сообщения if >0 -> 1 Else ->0
        e = xor(a,b); % вычисление вектора ошибок
        N_b = N_b + sum(e); % прибавление количества ошибок в данном эксперименте
        [~, s] = gfdeconv(double(b), crc); % вычисление синдрома s(x) = b(x) mod(g(x))
        if (bi2de(s) == 0) & (sum(e) ~= 0) % проверка на ошибку декодирования CRC
            Ne_decode = Ne_decode + 1;     % если синдром равен нулю, а вектор ошибок
        end                                % не равен нулю, значит декодер не смог обнаружить ошибку 
    end
    Pb(index) = N_b/(n*(m+r));
    Pe_decode(index) = Ne_decode/n;
    index = index + 1;
end

figure(1);
title('Вероятность ошибки на бит для BPSK');
semilogy(snr_arr_theor, Pb_theor, 'r', snr_arr_pract, Pb, 'bo');
legend('Теоретическая вероятность', 'Экспериментальная оценка');
figure(2);
title('Вероятность ошибки декодирования CRC');
semilogy(snr_arr_theor, Pe_decode_Theor, 'g-', snr_arr_pract, Pe_decode, 'bo', snr_arr_theor, P_ed, 'r-');
legend('Верхняя граница', 'Экспериментальная оценка', 'Точная оценка');
