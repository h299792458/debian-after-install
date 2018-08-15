%============================================================================
% This script is dedicated to section 11.7 of reference [1]. The signal
% model is an unfiltered AR(1) signal with WGN. Various techniques are
% demonstrated to estimate the AR(1) signal.
% 
% [1] Steven M. Kay, "Fundamentals of statistical signal processing:
% estimation theory." (1993).
%============================================================================
clear
PLOT = true;

%----------------------------------------------------------------------------
% Scalar Wiener estimator with scalar unknown and scalar observation.
% Multiple realizations are demonstrated as an iid sequence. The data model
% is given by X = S + W, where S and X are random variables. The estimation
% of S tends to be X when SNR is high and zero when SNR is low.
%----------------------------------------------------------------------------
% N = 50;
% snr = idbw(5);
% var_s = 1;
% var_w = var_s./snr;
% s = gaussian_noise(N, 1, var_s, 'linear', 'real');
% w = gaussian_noise(N, 1, var_w, 'linear', 'real');
% x = s + w;
% s_hat = snr / (snr + 1) * x;
% 
% if PLOT, 
%     figure;
%     plot(1:N, s, '--', 'linewidth', 1); hold on;
%     plot(1:N, x, 'linewidth', 1); 
%     plot(1:N, s_hat, 'linewidth', 2); grid on; 
%     xlabel('Sample number, n')
%     title('Scalar Wiener Filter Example');
% end

%----------------------------------------------------------------------------
% Finite Wiener smoother with vector unknowns and vector observations. The
% data model is X = S + W, where X and S are vector random variables or
% blocks of WSS random processes. The finite Wiener smoother performs
% estimation of each unknown in the block by using all the data samples in
% the block. Therefore, it should give the best results when compared with
% Wiener filter of lengh N in which case only N samples are used to
% estimate the current data
% 
% The Wiener smoother is implemented in time domain as a smoothing matrix,
% which is constructed based on the autocorrelation matrix of X and S. In
% this example, those AC matrices are given explicitly. Otherwise,
% estimated autocorrelation matrix could be used. Note that the smoothing
% matrix grows substantially with the block size and also matrix inversion
% is needed.
% 
% This method should be compared to the finite wiener smoother implemented
% in frequency domain via DFT, i.e., using the infinite wiener smoother
% given by PSD but implementing only discrete DFT points. The two methods
% have virtually the same performance with exceptional points at two ends
% of the data block. They might have theoretical equivalence yet to be
% proved...
% 
% Another method should also be compared to, e.g., the infinite wiener
% smoother implemented in time domain via linear convolution, but with
% finite impulse response. Multiple design methods could be applied to
% obtain the FIR filter in time domain. 
%----------------------------------------------------------------------------
L = 2500;
snr = idbw(0);

% Generate AR(1) process with WGN, note that the power of AR(1) is r[0]
a = -0.95;
var_u = 1;
s = arma(L, 1, [1, a], [1], var_u);
var_w = acfar1(var_u, a, 0) / snr;
w = gaussian_noise(L, 1, var_w, 'linear', 'real');
x = s + w;
Css = zeros(L);
for ii = 1:L
    for jj = 1:L
        Css(ii, jj) = acfar1(var_u, a, ii - jj);
    end
end
t_1 = Css / (Css + var_w * eye(L)) * x;
if PLOT, figure; mesh(1:L, 1:L, Css); xlabel('Sample number, n'); end

%----------------------------------------------------------------------------
% Wiener filter corresponding to the case of scalar unknown and vector
% observations. Could achieve the performance of finite Wiener smoother
% with a much shorter filter length.
%----------------------------------------------------------------------------
N = 11;
Nf = (N - 1) / 2;
Css = zeros(N);
for ii = 1:N
    for jj = 1:N
        Css(ii, jj) = acfar1(var_u, a, ii - jj);
    end
end
k = -Nf : Nf;
rss = acfar1(var_u, a, k); rss = rss.';
h_2 = (Css + var_w * eye(N)) \ rss;
t_2 = conv(x, h_2, 'same');

%----------------------------------------------------------------------------
% Finite Wiener smoother implemented in frequency domain
%----------------------------------------------------------------------------
omega = get_fft_grid(L, 2*pi);
Pss = var_u ./ abs(1 + a * exp(-1i * omega)).^2;
Hw = Pss ./ (Pss + var_w); Hw = transpose(Hw);
t_3 = real(ifft(fft(x) .* Hw));

%----------------------------------------------------------------------------
% Infinite Wiener filter implemented in time domain via digital filter
% design. Note that "Hw" is a sampled version of true Wiener smoother
% frequency response which is derived from theory in this example. However,
% taking directly the inverse DFT of "Hw" does not generate the sampled
% version of true continuous-time impluse response due to time aliasing
% (negligible in this case).
%----------------------------------------------------------------------------
M = 1024;
isodd = mod(M, 2);
htmp = ifft(transpose(Hw));
h_4 = fftshift([htmp(1 : Nf + isodd), htmp(end - Nf + isodd : end)]);
t_4 = conv(x, h_4, 'same');

if PLOT,
    figure; semilogy(fftshift(omega/pi), fftshift(abs(Hw)), 'o-'); grid on;
    figure; plot(1:L, s, 1:L, t_1, 1:L, t_2, 1:L, t_3, 1:L, t_4); grid on;
    xlabel('Sample number, n');
end
fprintf('Estimation MSE: %f, %f, %f, %f\n', mean(abs(s-t_1).^2), mean(abs(s-t_2).^2), mean(abs(s-t_3).^2), mean(abs(s-t_4).^2));
