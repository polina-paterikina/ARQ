clear all;
close all;
m = 8; %����� �������������� ������������������
crc = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1]; %crc16 m=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
% crc = [1 0 1 0 1 0 1 1 1]; %crc8 m=4 x8+x7+x6+x4+x2+1
% crc = [1 1 1 0 1]; % x4+x2+x+1 m=3
% crc = [1 1 0 1]; % x3+x+1 m=4
r = length(crc)-1;
n = 250000; %���������� �������������
snr_arr_pract = -30:5:10; %������/��� � �� ��� ����������������� ����������
snr_arr_theor = -30:10; %������/��� � �� ��� ������������� ����������


codewords = zeros(2^m, m+r); %������ ������� ����
mod_codewords = zeros(2^m, m+r); %������ ����������������� ������� ����
d_arr = zeros((2^m), 1); %������ ����� ���� ������� ����
for word = 0:2^m-1 % ���� �� ���� �������������� �������������������
    mxr = bitshift(word, r); % m(x)*x^r
    [~, c] = gfdeconv(de2bi(mxr), crc); % ����������� ����� c(x) = m(x)*x^r mod(g(x))
    a = de2bi(bitxor(mxr, bi2de(c)), m+r); % ������� ����� a(x) = m(x)*x^r + c(x)
    codewords(word+1, :) = a;
    d_arr(word+1) = sum(a); % ������� ���� �������� �����
    a_m = a.*2-1; % ��������� �������� ����� 1->+1, 0->-1
    mod_codewords(word+1, :) = a_m;
end

Pe_decode_Theor = ones(1,length(snr_arr_theor)).*(1/2^r); % ������� ������� ������ ������������� CRC
Pb_theor = qfunc(sqrt(2*(10.^(snr_arr_theor./10)))); % ������������� ����������� ������ �� ��� ��� BPSK

d_min = min(d_arr(2:end)); % ���������� ������������ ���������� ����
P_ed = zeros(1, length(Pb_theor)); % ������ ����������� ������ ������������� CRC
index = 1;
for p = Pb_theor % ���� �� ������������� ����������� ������ �� ��� ��� BPSK
    tmp = 0;
    for i = d_min:(m+r) % ���� �������� ����� �� ����� ������� ���� �� ������������ ���������� ���� d �� ����� �������� ����� n
        Ai = sum(d_arr == i); % ���������� ������� ���� ����� i
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i)); % ������� ����� �������� �������
    end
    P_ed(index) = tmp;
    index = index + 1;
end
index = 1;

Pb = zeros(1, length(snr_arr_pract));
Pe_decode = zeros(1, length(snr_arr_pract));
for SNRdB = snr_arr_pract % ���� �� ������� ������/��� � �� ��� ����������������� ����������
    N_b = 0; % ���������� ��������� ���
    Ne_decode = 0; % ���������� ������ ������������� CRC
    SNR = 10^(SNRdB/10); % ��������� ������/���
    sigma = sqrt(1/(2*SNR));
    parfor i = 1:n % ���� �������� �������� ���������� �������������
        rand_mess_ind = randi(2^m, 1); % ��������� ���������� ������� ���������
        a = codewords(rand_mess_ind, :); % ������� �����
        a_m = mod_codewords(rand_mess_ind, :); % ����������������� ������� �����
        b_m = a_m +sigma*randn(1, length(a_m)); % ���������� ���� � ������
        b = b_m > 0; % ����������� ��������� ��������� if >0 -> 1 Else ->0
        e = xor(a,b); % ���������� ������� ������
        N_b = N_b + sum(e); % ����������� ���������� ������ � ������ ������������
        [~, s] = gfdeconv(double(b), crc); % ���������� �������� s(x) = b(x) mod(g(x))
        if (bi2de(s) == 0) & (sum(e) ~= 0) % �������� �� ������ ������������� CRC
            Ne_decode = Ne_decode + 1;     % ���� ������� ����� ����, � ������ ������
        end                                % �� ����� ����, ������ ������� �� ���� ���������� ������ 
    end
    Pb(index) = N_b/(n*(m+r));
    Pe_decode(index) = Ne_decode/n;
    index = index + 1;
end

figure(1);
title('����������� ������ �� ��� ��� BPSK');
semilogy(snr_arr_theor, Pb_theor, 'r', snr_arr_pract, Pb, 'bo');
legend('������������� �����������', '����������������� ������');
figure(2);
title('����������� ������ ������������� CRC');
semilogy(snr_arr_theor, Pe_decode_Theor, 'g-', snr_arr_pract, Pe_decode, 'bo', snr_arr_theor, P_ed, 'r-');
legend('������� �������', '����������������� ������', '������ ������');
