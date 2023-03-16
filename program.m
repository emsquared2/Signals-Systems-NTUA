% EMMANOUILIDIS EMMANOUIL    
% 03119435

% Exercise 1.1: Echo - Reverb

% 1.1 a)

a = 1;
b2 = [0.55 0 0.45];
b5 = [0.55 0 0 0 0 0.45];


% 1.1 b)


figure (1)
freqz(b2,a)
title("Frequency Response of Echo Effect Filter for P = 2")
figure (2)
freqz(b5,a)
title("Frequency Response of Echo Effect Filter for P = 5")


%1.1 c) 


figure (3)
zplane(b2, a)
title("Zero-Pole Diagram of Echo Effect Filter for P = 2")
[z_2, p_2, k2] = tf2zpk(b2, a); % H tf2zp de douleuei kala me arnhtikous ekthetes
figure (4) 
zplane(b5, a)
title("Zero-Pole Diagram of Echo Effect Filter for P = 5")
[z5, p_5, k5] = tf2zpk(b5, a);


%1.1 d)


figure (5)
impz(b2, a)
title("Impulse Response of Echo Effect Filter for P = 2")
figure (6)
impz(b5, a)
title("Impulse Response of Echo Effect Filter for P = 5")


% 1.1. e)


[h, t_1e] = impz(b2, a);
h_sec = conv(h , h);
h_reverb = conv(h_sec , h);

Thn = toeplitz(h_reverb);

ai = Thn\h_reverb;

bi = h_reverb;
% epeidh ai = 0 gia oles tis times ektos ths prwths 
% prokuptei oti bi = h(i)  => bi = h_reverb

figure (7)
freqz(bi, ai)
title("Frequency Response of Reverb Effect Filter for P = 2")

[zi, pi, ki] = tf2zpk(bi, ai);

figure (8)
zplane(zi, pi)
title("Zero-Pole Diagram of Reverb Effect Filter for P = 2")


figure (9)
impz(bi, ai)
title("Impulse Response of Reverb Effect Filter for P = 2")

% 1.1 st)

H_reverb = fft(h_reverb, 512);
H_dereverb = 1 ./ H_reverb;
h_dereverb = ifft(H_dereverb, 512);

Thn_dereverb = toeplitz(h_dereverb);
ai_dereverb = Thn_dereverb\h_dereverb; 

bi_dereverb = h_dereverb;
% ai_dereverb = 1 0 ... 0 => bi_dereverb = h_dereverb

% x[n] = u[n] - u[n-5]
n = -5:1:10;
x = heaviside(n) - heaviside(n-5);
figure (10)
stem(n, x);
title("input signal x[n]")
xlabel("Samples")
ylabel("Amplitude")

x_reverb = conv(x, h_reverb);
x_dereverb = conv(x_reverb, h_dereverb);

n1 = -5:1:527;
figure (11)
stem(n1, x_dereverb);
title("output signal")
xlabel("Samples")
ylabel("Amplitude")
xlim([-5, 10]);
ylim([0, 1]);


% 1.2 Bandpass Filters 


%1.2 a)

p1 = [(0.65 + 0.65i) (0.65 - 0.65i) (0.65 + 0.65i) (0.65 - 0.65i)]';
z1 = [0.8 0.8 0.8i -0.8i]';
figure (12)
zplane(z1, p1)
title("Zero-Pole Diagram of Bandpass Filter")
k = 1;  %thetoume kerdos k = 1
[b ,a1] = zp2tf(z1, p1, k);  


%1.2 b)


figure (13)
freqz(b, a1);
title("Frequency Response of Bandpass Filter")


%1.2 c)


figure (14)
impz(b, a1)
title("Impulse Response of Bandpass Filter")
figure (15)
stepz(b, a1)
title("Step Response of Bandpass Filter")


%1.2 d)

p2 = [(0.7 + 0.7i) (0.7 - 0.7i) (0.7 + 0.7i) (0.7 - 0.7i)]';
[b_2, a2] = zp2tf(z1, p2, k); 
figure (16)
zplane(z1,p2)
figure (17)
freqz(b_2, a2);
title("Frequency Response of Bandpass Filter (p2)")
figure (18)
impz(b_2, a2);
title("Impulse Response of Bandpass Filter (p2)")


p3 = [(0.707 + 0.707i) (0.707 - 0.707i) (0.707 + 0.707i) (0.707 - 0.707i)]';
figure (19)
zplane(z1,p3)
[b_3, a3] = zp2tf(z1, p3, k);
figure (20)
freqz(b_3, a3);
title("Frequency Response of Bnadpass Filter (p3)")
figure (21)
impz(b_3, a3);
title("Impulse Response of Bandpass Filter (p3)")


p4 = [(0.75 + 0.75i) (0.75 - 0.75i) (0.75 + 0.75i) (0.75 - 0.75i)]';
figure (22)
zplane(z1,p4)
[b_4, a4] = zp2tf(z1, p4, k);
figure (23)
freqz(b_4, a4);
title("Frequency Response of Bandpass Filter (p4)")
figure (24)
impz(b_4, a4);
title("Impulse Response of Bandpass Filter (p4)")


%1.2 e) 


p5 = [(0.4 + 0.7i) (0.4 - 0.7i) (0.4 + 0.7i) (0.4 - 0.7i)]';
figure (25)
zplane(z1, p5);
title("Zero-Pole Diagram of Banpass Filter (p5)")

[b_5, a5] = zp2tf(z1, p5, k);
figure (26)
freqz(b_5, a5)
title("Frequency Response of Bandpass Filter (p5)")


%1.2 st)


z2 = [(0.77 + 0.2i) (0.77 - 0.2i) (0.2 + 0.77i) (0.2 - 0.77i)]';
figure (27)
zplane(z2, p1) 
title("Zero-Pole Diagram of BandPass Filter (z2)")

[b_6, a6] = zp2tf(z2, p1, k);
figure (28)
freqz(b_6, a6)
title("Frequency Response of Bandpass Filter (z2)")


z3 = [(0.4 + 0.7i) (0.4 - 0.7i) (0.7 + 0.4i) (0.7 - 0.4i)]';
figure (29)
zplane(z3, p1)
title("Zero-Pole Diagram of Bandpass Filter (z3)")

[b_7, a7] = zp2tf(z3, p1, k);
figure (30)
freqz(b_7, a7) 
title("Frequency Response of Bandpass Filter (z3)")



% ------------------------------------------------------------


% Exercise 2
% 2.1 Analysis of Musical Sounds


% 2.1 a)


[y_viola, Fs_viola] = audioread("viola_series.wav");

sound(y_viola, Fs_viola);
dt = 1/Fs_viola;
t = 0:dt:(length(y_viola)*dt)-dt;
plot(t,y_viola);
xlabel("Seconds");
ylabel("Amplitude");
title("viola Series")


%2.1 b)


y_viola_norm = normalize(y_viola, 'range', [-1 1]);
figure (31)
plot(t,y_viola_norm);
xlabel("Seconds");
ylabel("Amplitude");
title("viola Series");


w = hamming(1000);
E = conv(y_viola_norm.^2, w);
t1 = 0:dt:(length(E)*dt)-dt;
y2_viola_norm = normalize(y_viola, "range", [-60 60]);
figure (32)
plot(t1, E, 'magenta');
hold on
plot(t, y2_viola_norm, 'blue');
hold off
xlabel("Seconds");
ylabel("Amplitude");
title("Normalized Viola Series & Energy");


% 2.1 c)


Y_viola = fft(y_viola);
N = length(Y_viola);
fy_viola = 0:Fs_viola/N:Fs_viola/2-Fs_viola/N;
figure (33)
plot(fy_viola, abs(Y_viola(1:N/2)));
xlim([0 6000]);
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Spectrum of Viola Series");


% 2.1 d)


figure (34)
plot(t, y_viola);
xlim([0 2]);
xlabel("Seconds")
ylabel("Amplitude")
title("viola Series - Note")


figure (334) 
plot(t, y_viola);
xlim([0.1 0.12]);
xlabel("Seconds")
ylabel("Amplitude")
title("viola Series - Note (Zoomed)")
% Ypologizoume To = 3.88 ms (peripou)


% 2.1 e)

timeForNote = (t > 0) & (t < 2);
Y_note = fft(y_viola(timeForNote), Fs_viola);
N_note = length(Y_note);
fy_note = 0:Fs_viola/N_note:Fs_viola/2-Fs_viola/N_note;
figure (35)
plot(fy_note, abs(Y_note(1:N_note/2)));
xlabel("Frequency (Hz)");
ylabel("Amplitude")
title("Spectrum of Viola Series - Note")
xlim([0 6000]);

% H megisth timh tou fasmatos einai sta 259 Hz, 
% ara Fo = 259 Hz => To = 3.86 ms. (polu konta sthn arxikh timh 3.88)


% 2.1 st)


[y_viola_note, Fs_viola_note] = audioread("viola_note.wav");
dt_viola_note = 1/Fs_viola_note;

t_viola_note = 0:dt_viola_note:(length(y_viola_note)*dt_viola_note)-dt_viola_note;
figure (36)
plot(t_viola_note, y_viola_note);
xlabel("Seconds");
ylabel("Amplitude");
title("Viola Note");

Y_viola_note = fft(y_viola_note);
N_viola_note = length(Y_viola_note)-1; %-1 wste to N/2 na einai akeraios
fy_viola_note = 0:Fs_viola_note/N_viola_note:Fs_viola_note/2-Fs_viola_note/N_viola_note;
figure (37)
plot(fy_viola_note, abs(Y_viola_note(1:N_viola_note/2)));
xlim([0 6000]);
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Spectrum of Viola Note");


% 2nd Harmonic


p_viola = [(0.993 + 0.09i) (0.993 - 0.09i) (0.993 + 0.09i ) (0.993 - 0.09i)]';
z_viola = [0.8 -0.8 0.85i -0.85i]';
figure (38)
zplane(z_viola, p_viola);
title("Zero-Pole Diagram of Bandpass Filter (2nd Harmonic)");

[b_viola, a_viola] = zp2tf(z_viola, p_viola, 0.00001);
figure (39)
freqz(b_viola, a_viola);
title("Frequency Response of Bandpass Filter");

h_viola = impz(b_viola, a_viola);
y_filter_1 = conv(y_viola_note, h_viola);
figure (40)
plot(y_filter_1);
xlabel("Time");
ylabel("Amplitude");
title("Filtered Viola Note (2nd Harmonic)");


Y_filter_1 = fft(y_filter_1);
N_filter_1 = length(Y_filter_1)-1;
fy_filter_1 = 0:Fs_viola_note/N_filter_1:Fs_viola_note/2-Fs_viola_note/N_filter_1;
figure (41)
plot(fy_filter_1, abs(Y_filter_1(1:N_filter_1/2)));
xlim([0 6000]);
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Spectrum of 2nd Harmonic");



% 4th Harmonic


p_viola_2 = [(0.98 + 0.185i) (0.98 - 0.185i) (0.98 + 0.185i ) (0.98 - 0.185i )]';
z_viola_2 = [0.8 -0.8 0.85i -0.85i]';
figure (42)
zplane(z_viola_2, p_viola_2);
title("Zero-Pole Diagram of Bandpass Filter (4th Harmonic)");

[b_viola_2, a_viola_2] = zp2tf(z_viola_2, p_viola_2, 0.0001);
figure (43)
freqz(b_viola_2, a_viola_2);
title("Frequency Response of Bandpass Filter");

h_viola_2 = impz(b_viola_2, a_viola_2);
y_filter_2 = conv(y_viola_note, h_viola_2);
figure (44)
plot(y_filter_2);
xlabel("Time");
ylabel("Amplitude");
title("Filtered Viola Note (4th Harmonic)");


Y_filter_2 = fft(y_filter_2);
N_filter_2 = length(Y_filter_2);
fy_filter_2 = 0:Fs_viola_note/N_filter_2:Fs_viola_note/2-Fs_viola_note/N_filter_2;
figure (45)
plot(fy_filter_2, abs(Y_filter_2(1:N_filter_2/2)));
xlim([0 6000]);
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Spectrum of 4th Harmonic");



%2.2 Reverb Echo Effect on Musical Sounds


%2.2 a)


[y_piano, Fs_piano] = audioread("piano_note.wav");

sound(y_piano, Fs_piano);
dt_piano = 1/Fs_piano;
t_piano = 0:dt_piano:(length(y_piano)*dt_piano)-dt_piano;
plot(t_piano,y_piano);
xlabel("Seconds");
ylabel("Amplitude");
title("Piano Note")


%2.2 b)


% time delay = 0.15s
% (1/Fs_piano) * P = 0.15 => P = 6615

b_echo = zeros(1, 6615);
b_echo(1) = 0.6;
b_echo(6615) = 0.4;

a_echo = zeros(1, 6615);
a_echo(1) = 1;

[h_echo, t_echo] = impz(b_echo, a_echo);

echoed_piano_note = filter(b_echo, a_echo, y_piano);
echoed_t_piano = 0:dt_piano:(length(echoed_piano_note)*dt_piano) - dt_piano;

sound(echoed_piano_note, Fs_piano);

figure (46)
plot(echoed_t_piano, echoed_piano_note);
xlabel("Seconds");
ylabel("Amplitude");
title("Echoed Piano Note");
ylim([-0.03, 0.03])


second_h = conv(h_echo, h_echo);
reverb_h = conv(h_echo, second_h);  % 3 filtra se seira


%Thn = toeplitz(reverb_h);
%ai_piano_reverb = Thn\reverb_h;

%apo thn parapanw entolh prokuptei pws ai_piano_reverb = [1 0 0 ... 0]' ,
%opote gia na mh xreiazetai na upologizetai synexeia o pinakas toeplitz,
%(kati pou pairnei xrono) orizw ai_piano_reverb = 1 apeutheias:
ai_piano_reverb = 1;

% epeidh ai_piano_reverb = 0 gia oles tis times ektos ths prwths
% prokuptei oti bi_piano_reverb = h(i)  => bi_piano_reverb = reverb_h
bi_piano_reverb = reverb_h;

reverbed_piano_note = filter(bi_piano_reverb, ai_piano_reverb, y_piano);
reverbed_t_piano = 0:dt_piano:(length(reverbed_piano_note)*dt_piano) - dt_piano;

sound(reverbed_piano_note, Fs_piano);

figure (47)
plot(reverbed_t_piano, reverbed_piano_note);
xlabel("Seconds");
ylabel("Amplitude");
title("Reverbed Piano Note");


%2.2 c)


Y_piano = fft(y_piano);
figure (48)
plot(20*log10(abs(Y_piano(2:end))));
xlabel("Frequency (Hz)");
ylabel("Amplitude (db)");
title("Piano Note Spectrum");


echoed_Y_piano = fft(echoed_piano_note);
figure (49)
plot(20*log10(abs(echoed_Y_piano(2:end))));
xlabel("Frequency (Hz)");
ylabel("Amplitude (db)");
title("Echoed Piano Note Spectrum");


reverbed_Y_piano = fft(reverbed_piano_note);
figure (50)
plot(20*log10(abs(reverbed_Y_piano(2:end))));
xlabel("Frequency (Hz)");
ylabel("Amplitude (db)");
title("Reverbed Piano Note Spectrum");


%2.2 d)


b2_echo = zeros(1, 2021);
b2_echo(1) = 0.6;
b2_echo(2021) = 0.4;

a2_echo = zeros(1, 2021);
a2_echo(1) = 1;

[h2_echo, t2_echo] = impz(b2_echo, a2_echo);

echoed_piano_note_2 = conv(y_piano, h2_echo);
echoed_t_piano_2 = 0:dt_piano:(length(echoed_piano_note_2)*dt_piano) - dt_piano;

sound(echoed_piano_note_2, Fs_piano);

figure (51)
plot(echoed_t_piano_2, echoed_piano_note_2);
xlabel("Seconds");
ylabel("Amplitude");
title("Echoed Piano Note 2");
ylim([-0.03, 0.03])


%2.2 e)

audiowrite("echoed_piano_note.wav", echoed_piano_note, Fs_piano);
audiowrite("reverbed_piano_note.wav", reverbed_piano_note, Fs_piano);


%2.2 st)


%de xrhsimopoihthke idia diadikasia me to erwthma 1.1 st), opou
%xrhsimopoihsame ton pinaka toeplitz, dioti to megethos htan arketa megalo
%kai o upologismos tou pinaka toeplitz htan xronovora


ai_piano_dereverb = bi_piano_reverb;
bi_piano_dereverb = ai_piano_reverb;

dereverbed_piano_note = filter(bi_piano_dereverb, ai_piano_dereverb, reverbed_piano_note);

dereverbed_t_piano = 0:dt_piano:(length(dereverbed_piano_note)*dt_piano) - dt_piano;

figure (52)
plot(dereverbed_t_piano, dereverbed_piano_note);
hold on
plot(t_piano, y_piano);
xlabel("Seconds");
ylabel("Amplitude");
title("Piano Note & Dereverbed Piano Note");
hold off


%2.2 z)


% 5 samples

b_echo_5 = zeros(1, 6620);
b_echo_5(1) = 0.6;
b_echo_5(6620) = 0.4;

a_echo_5 = zeros(1, 6620);
a_echo_5(1) = 1;

[h_echo_5, t_echo_5] = impz(b_echo_5, a_echo_5);

h_second_5 = conv(h_echo_5, h_echo_5);
h_piano_reverb_5 = conv(h_echo_5, h_second_5);

ai_piano_reverb_5 = 1;
bi_piano_reverb_5 = h_piano_reverb_5;

ai_piano_dereverb_5 = bi_piano_reverb_5;
bi_piano_dereverb_5 = ai_piano_reverb_5;

dereverbed_piano_note_5 = filter(bi_piano_dereverb_5, ai_piano_dereverb_5, reverbed_piano_note);

dereverbed_t_piano_5 = 0:dt_piano:(length(dereverbed_piano_note_5)*dt_piano)-dt_piano;

figure (53)
plot(dereverbed_t_piano_5, dereverbed_piano_note_5);
xlabel("Seconds");
ylabel("Amplitude");
title("Dereverbed Piano Note - 5 Samples Off");

ai_final_5 = ai_piano_dereverb_5 .* ai_piano_reverb;
bi_final_5 = bi_piano_dereverb_5 .* bi_piano_reverb;

figure (54)
freqz(bi_final_5, ai_final_5);
title("Frequency response for 5 Samples Off")


%10 samples

b_echo_10 = zeros(1, 6625);
b_echo_10(1) = 0.6;
b_echo_10(6625) = 0.4;

a_echo_10 = zeros(1, 6625);
a_echo_10(1) = 1;

[h_echo_10, t_echo_10] = impz(b_echo_10, a_echo_10);

h_second_10 = conv(h_echo_10, h_echo_10);
h_piano_reverb_10 = conv(h_echo_10, h_second_10);

ai_piano_reverb_10 = 1;
bi_piano_reverb_10 = h_piano_reverb_10;

ai_piano_dereverb_10 = bi_piano_reverb_10;
bi_piano_dereverb_10 = ai_piano_reverb_10;

dereverbed_piano_note_10 = filter(bi_piano_dereverb_10, ai_piano_dereverb_10, reverbed_piano_note);

dereverbed_t_piano_10 = 0:dt_piano:(length(dereverbed_piano_note_10)*dt_piano)-dt_piano;

figure (55)
plot(dereverbed_t_piano_5, dereverbed_piano_note_10);
xlabel("Seconds");
ylabel("Amplitude");
title("Dereverbed Piano Note - 10 Samples Off");

ai_final_10 = ai_piano_dereverb_10 .* ai_piano_reverb;
bi_final_10 = bi_piano_dereverb_10 .* bi_piano_reverb;

figure (56)
freqz(bi_final_10, ai_final_10);
title("Frequency response for 10 Samples Off");


%50 samples

b_echo_50 = zeros(1, 6665);
b_echo_50(1) = 0.6;
b_echo_50(6665) = 0.4;

a_echo_50 = zeros(1, 6665);
a_echo_50(1) = 1;

[h_echo_50, t_echo_50] = impz(b_echo_50, a_echo_50);

h_second_50 = conv(h_echo_50, h_echo_50);
h_piano_reverb_50 = conv(h_echo_50, h_second_50);

ai_piano_reverb_50 = 1;
bi_piano_reverb_50 = h_piano_reverb_50;

ai_piano_dereverb_50 = bi_piano_reverb_50;
bi_piano_dereverb_50 = ai_piano_reverb_50;

dereverbed_piano_note_50 = filter(bi_piano_dereverb_50, ai_piano_dereverb_50, reverbed_piano_note);

dereverbed_t_piano_50 = 0:dt_piano:(length(dereverbed_piano_note_50)*dt_piano)-dt_piano;

figure (57)
plot(dereverbed_t_piano_50, dereverbed_piano_note_50);
xlabel("Seconds");
ylabel("Amplitude");
title("Dereverbed Piano Note - 50 Samples Off");

ai_final_50 = ai_piano_dereverb_50 .* ai_piano_reverb;
bi_final_50 = bi_piano_dereverb_50 .* bi_piano_reverb;

figure (58)
freqz(bi_final_50, ai_final_50);
title("Frequency response for 50 Samples Off")

%---------------------------------------------------

ai_final = ai_piano_reverb .* ai_piano_dereverb;
bi_final = bi_piano_reverb .* bi_piano_dereverb;

figure(59)
freqz(bi_final, ai_final);
title("Frequency Response of Initial System")

