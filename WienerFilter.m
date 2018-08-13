%============================================================================
% This script is dedicated to section 11.7 of the following reference:
% 
% Steven M. Kay, "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%============================================================================
clear
% close all
PLOT = true;

%----------------------------------------------------------------------------
% Scalar wiener estimator with scalar unknown and scalar observation.
% Multiple realizations are demonstrated as a sequence. The data model is
% given by x = s + w, where s and x are random variables. The estimation of
% s tends to be x when SNR is high and zero when SNR is low.
%----------------------------------------------------------------------------
N = 31;
snr = idbw(5);
var_s = 1;
var_w = var_s./snr;

% Random variable with normal distribution, generate a sequence of iid
% samples to represent multiple realizations, generate WGN as well
s = gaussian_noise(N, 1, var_s, 'linear', 'real');
w = gaussian_noise(N, 1, var_w, 'linear', 'real');
x = s + w;

% Scalar wiener filter: s_hat approaches x when SNR is high, becomes 0 when
% SNR is low
s_hat = snr / (snr + 1) * x;

if PLOT, 
    figure; plot(1:N, s, '--', 'linewidth', 1); hold on;
    plot(1:N, x, 'linewidth', 1); 
    plot(1:N, s_hat, 'linewidth', 2); grid on; 
    xlabel('Sample number, n')
    title('Scalar Wiener Filter Example');
end

%----------------------------------------------------------------------------
% Finite Wiener smoother with vector unknowns and vector observations. The
% data model is x = s + w, where x and s are vector random variables or
% blocks of random processes. In this particular example, the unknown is
% assumed to be a section of AR(1) process, the autocorrelation matrix of
% which is given explicitly. Otherwise, estimated autocorrelation matrix
% could be used for finite wiener smoother. Note that the matrix grows
% substantially with the block size and also matrix inversion is needed.
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
N = 31;
snr = idbw(5);

% Generate AR(1) process with WGN, note that the power of AR(1) is r[0]
a = -0.95;
var_u = 1;
s = arma(N, 1, [1, a], [1], var_u);
var_w = var_u / (1 - a^2) / snr;
w = gaussian_noise(N, 1, var_w, 'linear', 'real');
x = s + w;

% Finite Wiener smoother implemented in time domain
Css = zeros(N);
for ii = 1:N
    for jj = 1:N
        Css(ii, jj) = var_u / (1 - a^2) * (-a)^abs(ii - jj);
    end
end
if PLOT, figure; mesh(1:N, 1:N, Css); xlabel('Sample number, n'); end
s_hat_td = Css / (Css + var_w * eye(N)) * x;

if PLOT, 
    figure; plot(1:N, s, '--', 'linewidth', 1); hold on;
    plot(1:N, x, 'linewidth', 1); 
    plot(1:N, s_hat_td, 'linewidth', 2); grid on; 
    xlabel('Sample number, n');
    title('Vector Wiener Filter Example');
end

% Finite Wiener smoother implemented in frequency domain
omega = get_fft_grid(N, 2*pi);
Pss = var_u ./ abs(1 + a * exp(-1i*omega)).^2;
Hw = Pss ./ (Pss + var_w); Hw = transpose(Hw);
s_hat_fd = real(ifft(fft(x) .* Hw));

if PLOT,
    figure; semilogy(fftshift(omega/pi), fftshift(Hw), 'o-'); grid on;
    figure; plot(1:N, s, 1:N, s_hat_td, 1:N, s_hat_fd); grid on;
end
fprintf('Estimation MSE: TD %f, FD %f\n', mean(abs(s-s_hat_td).^2), mean(abs(s-s_hat_fd).^2));

% Infinite wiener filter implemented in time domain via digital filter
% design
N = 31;
Nf = (N - 1) / 2;
% Note that "Hw" is a sampled version of true Wiener smoother frequency
% response which is derived from theory in this example. However, taking
% directly the inverse DFT of "Hw" does not generate the sampled version of
% true continuous-time impluse response due to time aliasing (negligible in
% this case).
N = 1024;
isodd = mod(N, 2);
omega = get_fft_grid(N, 2*pi);
Pss = var_u ./ abs(1 + a * exp(-1i*omega)).^2;
Hw = Pss ./ (Pss + var_w); 
htmp = ifft(Hw);
h = fftshift([htmp(1 : Nf + isodd), htmp(end - Nf + isodd : end)]);
h = h ./ sum(abs(h));
if PLOT, figure; stem(h); grid on; end
s_hat_if = conv(x, h, 'same');

fprintf('Estimation MSE: TD %f, IF %f\n', mean(abs(s-s_hat_td).^2), mean(abs(s-s_hat_if).^2));
