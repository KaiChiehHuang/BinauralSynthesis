%% Audio Application of FFT Project - Binaural Synthesis %%
[input,fs] = audioread('audio/JustTonight.wav');
input_Left = input(:,1);
input_Right = input(:,2);

%% Upmixing Stage %%
center = (input_Left + input_Right) * 0.7071;
surround = (input_Left - input_Right) * 0.7071;
fc = 4000;
[b,a] = butter(10,fc/(fs/2),'low');
C = filter(b,a,center);
fc = 200;
[b,a] =butter(8,fc/(fs/2),'low');
LFE = filter(b,a,center);
fc = 7000;
[b,a] = butter(10,fc/(fs/2),'low');
delaySamples = round(0.012 * 44100);
surround = [zeros(delaySamples,1);surround];
surround = filter(b,a,surround);
RL = imag(hilbert(surround));
RL = imag(hilbert(RL)); % Make RL 180 degree out of phase with RR
RR = surround;
FR = input_Right;
FL = input_Left;

%% Make all the channels same length %%
C = [C;zeros(delaySamples,1)];
FR = [FR;zeros(delaySamples,1)];
FL = [FL;zeros(delaySamples,1)];
LFE = [LFE;zeros(delaySamples,1)];

%% Binaural Rendering Stage %%

hrir90_L = readraw('full/elev0/L0e090a.dat'); 
hrir90_R = readraw('full/elev0/R0e090a.dat'); 
hrir60_L = readraw('full/elev0/L0e060a.dat'); 
hrir60_R = readraw('full/elev0/R0e060a.dat'); 
hrir40_L = readraw('full/elev0/L0e040a.dat'); 
hrir40_R = readraw('full/elev0/R0e040a.dat'); 
hrir90_L = 0.5*[zeros(1,round(0.008*fs)) hrir90_L zeros(1,round(0.015*fs))];
hrir90_R = 0.5*[zeros(1,round(0.008*fs)) hrir90_R zeros(1,round(0.015*fs))];
hrir60_L = 0.3*[zeros(1,round(0.015*fs)) hrir60_L zeros(1,round(0.008*fs))];
hrir60_R = 0.3*[zeros(1,round(0.015*fs)) hrir60_R zeros(1,round(0.008*fs))];
hrir40_L = 0.2*[zeros(1,round(0.023*fs)+1) hrir40_L ];
hrir40_R = 0.2*[zeros(1,round(0.023*fs)+1) hrir40_R ];
hrir90_L = filter(b,a,hrir90_L);
hrir90_R = filter(b,a,hrir90_R);
hrir60_L = filter(b,a,filter(b,a,hrir60_L));
hrir60_R = filter(b,a,filter(b,a,hrir60_R));
hrir40_L = filter(b,a,filter(b,a,filter(b,a,hrir40_L)));
hrir40_R = filter(b,a,filter(b,a,filter(b,a,hrir40_R)));
% Late reverb generation %
wNoise = wgn(round(0.2*44100),1,1);
S = spectrogram(wNoise,512,256,1024,fs);
t60s = logspace(0,1,length(S(:,1)));
t60s = flip(t60s)/max(t60s) * 0.2;
for i = 1:length(S(:,1))
    S(i,:) = S(i,:) .* exp(-6.91*(256+512*(0:length(S(1,:))-1))/t60s(i)/44100);
end
[x,t] = istft(S,256,1024,fs);
lateReverb = 0.8*max(hrir40_L)*[zeros(1,round(0.024*fs)+1) x];
% Front Right
hrir30_L = readraw('full/elev0/L0e030a.dat'); 
hrir30_R = readraw('full/elev0/R0e030a.dat'); 
hrir30_L = [hrir30_L zeros(1,round(0.023*fs)+1)] + hrir90_L + hrir60_L + hrir40_L;
hrir30_R = [hrir30_R zeros(1,round(0.023*fs)+1)] + hrir90_R + hrir60_R + hrir40_R;
hrir30_L = [hrir30_L zeros(1,length(lateReverb)-length(hrir30_L))];
hrir30_R = [hrir30_R zeros(1,length(lateReverb)-length(hrir30_R))];
hrir30_L = hrir30_L + lateReverb;
hrir30_R = hrir30_R + lateReverb;
FR_left = filter(hrir30_L,1,FR);
FR_right = filter(hrir30_R,1,FR);
out_FR = [FR_left FR_right];
% Front Left
hrir330_L = readraw('full/elev0/L0e330a.dat'); 
hrir330_R = readraw('full/elev0/R0e330a.dat'); 
hrir330_L = [hrir330_L zeros(1,round(0.023*fs)+1)] + hrir90_L + hrir60_L + hrir40_L;
hrir330_R = [hrir330_R zeros(1,round(0.023*fs)+1)] + hrir90_R + hrir60_R + hrir40_R;
hrir330_L = [hrir330_L zeros(1,length(lateReverb)-length(hrir330_L))];
hrir330_R = [hrir330_R zeros(1,length(lateReverb)-length(hrir330_R))];
hrir330_L = hrir330_L + lateReverb;
hrir330_R = hrir330_R + lateReverb;
FL_left = filter(hrir330_L,1,FL);
FL_right = filter(hrir330_R,1,FL);
out_FL = [FL_left FL_right];
% Center
hrir000_L = readraw('full/elev0/L0e000a.dat'); 
hrir000_R = readraw('full/elev0/R0e000a.dat'); 
hrir000_L = [hrir000_L zeros(1,round(0.023*fs)+1)] + hrir90_L + hrir60_L + hrir40_L;
hrir000_R = [hrir000_R zeros(1,round(0.023*fs)+1)] + hrir90_R + hrir60_R + hrir40_R;
hrir000_L = [hrir000_L zeros(1,length(lateReverb)-length(hrir000_L))];
hrir000_R = [hrir000_R zeros(1,length(lateReverb)-length(hrir000_R))];
hrir000_L = hrir000_L + lateReverb;
hrir000_R = hrir000_R + lateReverb;
C_left = filter(hrir000_L,1,C);
C_right = filter(hrir000_R,1,C);
out_C = [C_left C_right];
% Rear Left
hrir250_L = readraw('full/elev0/L0e250a.dat'); 
hrir250_R = readraw('full/elev0/R0e250a.dat');
hrir250_L = [hrir250_L zeros(1,round(0.023*fs)+1)] + hrir90_L + hrir60_L + hrir40_L;
hrir250_R = [hrir250_R zeros(1,round(0.023*fs)+1)] + hrir90_R + hrir60_R + hrir40_R;
hrir250_L = [hrir250_L zeros(1,length(lateReverb)-length(hrir250_L))];
hrir250_R = [hrir250_R zeros(1,length(lateReverb)-length(hrir250_R))];
hrir250_L = hrir250_L + lateReverb;
hrir250_R = hrir250_R + lateReverb;
RL_left = filter(hrir250_L,1,RL);
RL_right = filter(hrir250_R,1,RL);
out_RL = [RL_left RL_right];
% Rear Right
hrir110_L = readraw('full/elev0/L0e110a.dat'); 
hrir110_R = readraw('full/elev0/R0e110a.dat'); 
hrir110_L = [hrir110_L zeros(1,round(0.023*fs)+1)] + hrir90_L + hrir60_L + hrir40_L;
hrir110_R = [hrir110_R zeros(1,round(0.023*fs)+1)] + hrir90_R + hrir60_R + hrir40_R;
hrir110_L = [hrir110_L zeros(1,length(lateReverb)-length(hrir110_L))];
hrir110_R = [hrir110_R zeros(1,length(lateReverb)-length(hrir110_R))];
hrir110_L = hrir110_L + lateReverb;
hrir110_R = hrir110_R + lateReverb;
RR_left = filter(hrir110_L,1,RR);
RR_right = filter(hrir110_R,1,RR);
out_RR = [RR_left RR_right];
% LFE 
LFE = [LFE LFE];
%% Play %%
out = out_FR + out_FL + out_C + out_RL + out_RR + LFE;
%soundsc(out,fs);
%audiowrite('audio/JustTonight_Surr.wav',out/max(max(abs(out))),fs);

%%  Plot brir and its frequency response
% figure;

% subplot(211);
% plot((0:length(hrir30_L)-1)*1/fs,hrir30_L,'-b');
% hold on;
% plot((0:length(hrir30_R)-1)*1/fs,hrir30_R,'-r');
% xlim([0 0.06]);
% ylabel('Magnitude (linear)');
% xlabel('Time (seconds)');
% legend('Left','Right');
% subplot(212);
% [h,w] = freqz(hrir30_L,1);
% plot(w/pi*fs/2,20*log10(abs(h)),'-b');
% hold on;
% [h,w] = freqz(hrir30_R,1);
% plot(w/pi*fs/2,20*log10(abs(h)),'-r');
% xlim([0 fs/2]);
% ylabel('Magnitude (dB)');
% xlabel('Frequency (Hz)');
% legend('Left','Right');

%% Make the loudness of the unprocessed match the processed
% fileNames = {'HowBlue','Yellow','YellowBrick','JesusEtc','RockTribe','WindBlow','YouShook','Angelina','BoldAsLove'};
% 
% for i = 1:9
%     fileName = fileNames{i};
%     [input,fs] = audioread(strcat(strcat('audio/',fileName),'.wav'));
%     audiowrite(strcat(strcat('audio/',fileName),'_Low.wav'),0.6*input(1:60*fs,:),fs);
% end