%===========================================================================
% Noise analysis on phase recovery based on pilot-tone. The phase noise is
% first estimated from a lowpass filtered pilot-tone, then the desired
% signal is compensated with estimated phase.
%===========================================================================
clear

fc = 1e6;
fs = 20e6;
nsample = 10^5;
t = (0 : (1/fs) : (nsample-1)/fs)';
pn = phase_noise(nsample, 1e-3, 0);
an = gaussian_noise(nsample, 1, .3, 'linear', 'complex');
x = exp(1i * pn(:)) + an;
H = frequency_response(nsample, fs, 0.01, 6e6, 'rc');
xf = ifft(fft(x) .* H);
xc = x .* conj(xf) ./ abs(xf);
nc = an .* conj(xf) ./ abs(xf);

% observe there is a peak in nc, and xc has lower WGN PSD
spectrumAnalyzer([x, xc], [], fs);
spectrumAnalyzer([an(:), nc(:)], [], fs);
