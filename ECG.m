clear, clc, close all
% Arduino setup and data acquisition
a = arduino();
% Reduce the number of samples taken
num_samples = 1000; % Number of samples
ecg = zeros(1, num_samples);
% Start acquiring ECG data from Arduino
figure;
for k = 1:num_samples
   b = readVoltage(a, 'A0');
   ecg(k) = b;
   subplot(2, 1, 1);
   plot(1:k, ecg(1:k)); % Plot every time new data is added
   grid on;
   title('Sinyal ECG Asli');
   xlabel('Sampel');
   ylabel('Amplitudo');
   drawnow;
end
% Calculation of sample interval and frequency
TS = 1; % because original data is sampled every 1 sample
Freq = 1 / TS;
%%%%% Kompresi sinyal ECG menggunakan transformasi DCT %%%%
dct_ecg = dct(ecg);
% Membuang koefisien DCT yang memiliki magnitudo kecil
threshold_dct = 0.001 * max(abs(dct_ecg));  % Threshold pada magnitudo DCT
dct_ecg_compressed = dct_ecg;
dct_ecg_compressed(abs(dct_ecg) < threshold_dct) = 0;
% Dekompresi sinyal ECG
ecg_decompressed_dct = idct(dct_ecg_compressed);
% Menampilkan hasil kompresi dan dekompresi
figure;
% Plot hasil DCT
subplot(4, 1, 1);
plot(1:num_samples, ecg);
title('Sinyal ECG Asli');
xlabel('Sampel');
ylabel('Amplitudo');
subplot(4, 1, 2);
stem(dct_ecg_compressed), axis([0 num_samples -2 2]);
title('Koefisien DCT Setelah Kompresi');
xlabel('Koefisien');
ylabel('Amplitudo');
subplot(4, 1, 3);
plot(1:num_samples, ecg_decompressed_dct);
title('Sinyal ECG Setelah Dekompresi DCT');
xlabel('Sampel');
ylabel('Amplitudo');
% Perhitungan Compression Ratio dan PRD untuk DCT
% Perhitungan counterB_DCT
counterB_DCT = sum(dct_ecg ~= 0);
% Menghitung counterA_DCT
counterA_DCT = sum(dct_ecg_compressed ~= 0);
% Menghitung Compression Ratio DCT
CompRatio_DCT = ((counterB_DCT - counterA_DCT) / counterB_DCT) * 100;
% Menghitung PRD DCT
error_DCT = ecg - ecg_decompressed_dct;
PRD_DCT = sqrt(sum(error_DCT.^2) / sum(ecg.^2)) * 100;
% Tampilkan Compression Ratio dan PRD DCT
fprintf('Threshold: %.4f\n', threshold_dct);
fprintf('Compression Ratio DCT: %.2f%%\n', CompRatio_DCT);
fprintf('PRD DCT: %.2f%%\n', PRD_DCT);
%%%%% Kompresi sinyal ECG menggunakan transformasi FFT %%%%%
fft_ecg = fft(ecg);
% Membuang komponen frekuensi yang lebih rendah
threshold_fft = 0.001 * max(abs(fft_ecg)); % Threshold pada magnitudo FFT
fft_ecg_compressed = fft_ecg;
fft_ecg_compressed(abs(fft_ecg) < threshold_fft) = 0;
% Dekompresi sinyal ECG
ecg_decompressed_fft = ifft(fft_ecg_compressed);
% Plot hasil FFT
figure;
subplot(4, 1, 1);
plot(1:num_samples, ecg);
title('Sinyal ECG Asli');
xlabel('Sampel');
ylabel('Amplitudo');
subplot(4, 1, 2);
stem(fft_ecg_compressed), axis([0 num_samples -2 2]);
title('Komponen Koefisien FFT Setelah Kompresi');
xlabel('Frekuensi');
ylabel('Amplitudo');
subplot(4, 1, 3);
plot(1:num_samples, ecg_decompressed_fft);
title('Sinyal ECG Setelah Dekompresi FFT');
xlabel('Sampel');
ylabel('Amplitudo');
% Perhitungan Compression Ratio dan PRD untuk FFT
% Perhitungan counterB_FFT
counterB_FFT = sum(abs(fft_ecg) ~= 0);
% Perhitungan counterA_FFT
counterA_FFT = sum(abs(fft_ecg_compressed) ~= 0);
% Menghitung Compression Ratio FFT
CompRatio_FFT = ((counterB_FFT - counterA_FFT) / counterB_FFT) * 100;
% Menghitung PRD FFT
error_FFT = ecg - ecg_decompressed_fft;
PRD_FFT = sqrt(sum(error_FFT.^2) / sum(ecg.^2)) * 100;
% Tampilkan Compression Ratio dan PRD FFT
fprintf('Threshold_FFT: %.4f\n', threshold_fft);
fprintf('Compression Ratio FFT: %.2f%%\n', CompRatio_FFT);
fprintf('PRD FFT: %.2f%%\n', PRD_FFT);
