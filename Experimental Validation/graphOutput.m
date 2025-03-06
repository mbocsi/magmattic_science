function graphOutput(fileName)
% Graphs the voltage signal and fourier transform derived from the
% experimental values
%
% fileName = the json filename
%  

str = fileread(fileName); 
struct0 = jsondecode(str);
N = 600;%Turns

spinFreq = 120;%Hz
omega = 2 * pi * spinFreq;
x = str2double(extractBefore(fileName, "cm")) / 100; %m

l = 0.0889; %m
A = l^2;

Br = 1 / 13200; %T Remanent magnetic flux density 
mu0 = (4 * pi) * 10^-7;
R  = 0.009525; %m
T = 0.009525; %m
volume = (pi * R^2) * T;

M = Br/mu0;
mom = M*volume;

if x == 0
    dist = x + R/2;
else
    dist = x;
end

B = (mu0 * mom) / (2 * pi * dist^3);

VcoilSim = B * omega * A * N;

Gain = 100;
VcoilSimAmp = Gain * VcoilSim;

%Voltage Signals
volt0 = struct0.value.voltageData;
t0 = zeros(length(volt0),1);
voltageVals0 = zeros(length(volt0),1);

for i = 1:length(volt0)
    t0(i) = volt0(i).value(1) - volt0(1).value(1);
    voltageVals0(i) = volt0(i).value(2);
end

%FFT Signal
fft0 = struct0.value.fftData;
f0 = zeros(length(fft0),1);
fftVals0 = zeros(length(fft0),1);

for i = 1:length(fft0)
    f0(i) = fft0(i).value(1) - fft0(1).value(1);
    fftVals0(i) = fft0(i).value(2);
end

%Normalizing the magnitude
fftVals0Max = max(fftVals0);
fftVals0 = fftVals0 / max(fftVals0);

%Matlab fft
fftGraphMat = fft(voltageVals0);
%Normalizing
fftMatMax = max(fftGraphMat);
fftGraphMat = fftGraphMat / fftMatMax;

figure 

subplot(3,1,1)
plot(t0, voltageVals0)
xlabel('Time (ms)')
ylabel('Voltage (V)')
title("Voltage Signal")

subplot(3,1,2)
plot(f0, abs(fftVals0))
xlabel('Frequency (Hz)')
ylabel('Magnitude Normalized')
xlim([0 500])
title("FFT Graph")

subplot(3,1,3)
plot(abs(fftGraphMat))
xlabel('Frequency (Hz)')
ylabel('Magnitude Normalized')
xlim([0 500])
title("FFT Graph Derived from Signal")

graphsTitle = extractBefore(fileName, ".");
sgtitle(graphsTitle)

pk2pk = peak2peak(voltageVals0);
fprintf("Delta V (" + graphsTitle +") = " + string(pk2pk) + " V\n")
fprintf("Delta V Sim (" + graphsTitle +") = " + string(VcoilSimAmp) + " V\n")

ENBW = 1/T;
simNoise = 12.66e-6; %V/sqrt(Hz)
noiseFloor = simNoise * sqrt(ENBW);
fprintf("Noise Floor (" + graphsTitle +") = " + string(noiseFloor) + " V\n")



end