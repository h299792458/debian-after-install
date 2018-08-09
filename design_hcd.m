function hcd = design_hcd(K, N)
% Design the time-domain impulse response of the chromatic dispersion
% compensation filter

if nargin < 2, N = 111; end
if nargin < 1, K = 2; end

PLOT = false;

N = N + 1 - mod(N, 2);
Nf = (N - 1) / 2;
M = 8193;
isodd = mod(M, 2);
a = 0.22;
L = 1730;
epsilon = 1e-8;
eta = 0.01;

%----------------------------------------------------------------------------
% Filter #1 designed by the time sampling method (TSM), i.e., taking the
% inverse Fourier transform of the all-pass chromatic dispersion frequency
% response with angular frequency from -inf to inf. This gives an impulse
% response with constant magnitude and a corresponding DTFT with
% non-constant magnitude
%----------------------------------------------------------------------------
hcd_A1 = ones(1, N);
for n = -Nf : Nf
    hcd_A1(n + Nf + 1) = sqrt(1i / (4 * K * pi)) * exp(-1i * n^2 / (4 * K));
end

%----------------------------------------------------------------------------
% Filter #2 is designed by the frequency sampling method (FSM), i.e.,
% taking the inverse DFT of the conjugate of discrete chromatic dispersion
% frequency response H(w) with digital frequencies from -pi to pi. This
% will produce an impulse response with non-constant magnitude and a
% corresponding DTFT with N points specified. This corresponds to the
% zero-forcing solution of the linear model Y(w) = H(w)X(w) + W(w), where
% X(w) is the pulseshaped original signal spectrum with two-fold
% oversampling. Note that the operation of fftshift introduces a linear
% phase to the designed frequency response.
%----------------------------------------------------------------------------
W = [0 : Nf, -Nf : -1] * 2 * pi / N;
H = exp(1i * K * W.^2);
hcd_A2 = fftshift(ifft(H));

%----------------------------------------------------------------------------
% Filter #3 is designed by FSM with oversampling. Note that this is
% equivalent to solving the lease squares problem of D*h[n] = H[k], where D
% is an MxN DFT matrix, h has N points, and H has M points. The solution is
% given by h = inv(D'*D)*D'*H. When M is large and the desired frequency
% response is sampled from -pi to pi, the least squares criteria approaches
% the mean squared error criteria over the full band, and according the
% convergence of Fourier series, the solution of minimizing the MSE is the
% ideal Fourier series truncated with an N-point rectangular windown. On
% the other hand, oversampling the ideal frequency response heavily reduces
% time aliasing and therefore the inverse DFT produces an impulse response
% close enough to the ideal impulse response.
%----------------------------------------------------------------------------
W = get_fft_grid(M, 2 * pi);
H = exp(1i * K * W.^2);
hcd_A3 = ifft(H);
hcd_A3 = fftshift([hcd_A3(1 : Nf + isodd), hcd_A3(end - Nf + isodd : end)]);

% This can be formulated as a overdetermined system problem, i.e., design
% an N points time domain impulse response h[n] such that the DTFT of h[n]
% matches the desired DTFT at M (M > N) points in the least squares sense

% D = dftmtx(M);
% D = [D(:, 1 : Nf + isodd), D(:, end - Nf + isodd : end)];
% hcd_A3 = (D' * D + 1e-6 * eye(N)) \ D' * H;
% hcd_A3 = fftshift(transpose(hcd_A3));

%----------------------------------------------------------------------------
% Filter #4 is designed as a frequency weighted version of filter 3, i.e.,
% instead of minimizing the MSE of DTFT within the full band, only the
% passband of the pulseshaper is considered in the optimization
%----------------------------------------------------------------------------
% W = get_fft_grid(M, 2 * pi);
% H = exp(1i * K * W.^2);
% H = [H(1 : M/2 - L), H(end - M/2 + L + 1 : end)];
% D = dftmtx(M);
% D = [D(1 : M/2 - L, :); D(end - M/2 + L + 1 : end, :)];
% D = [D(:, 1 : Nf + isodd), D(:, end - Nf + isodd : end)];
% hcd_A4 = (D' * D + epsilon * eye(N)) \ D' * H(:);
% hcd_A4 = fftshift(transpose(hcd_A4));

W = get_fft_grid(M, 2 * pi);
H = exp(1i * K * W.^2);
H(M/2 - L + 1 : end - M/2 + L) = 0;
D = dftmtx(M);
D = [D(:, 1 : Nf + isodd), D(:, end - Nf + isodd : end)];
hcd_A4 = (D' * D + epsilon * eye(N)) \ D' * H(:);
hcd_A4 = fftshift(transpose(hcd_A4));

%----------------------------------------------------------------------------
% filter #4 and #5 designed by windowing filter #1 by raised-cosine with
% various roll-off
%----------------------------------------------------------------------------
% arc = 0.15;
% hrc = rcshape(Nf, arc);
% hcd_A4 = hcd_A1 .* hrc;
% 
% arc = 0.85;
% hrc = rcshape(Nf, arc);
% hcd_A5 = hcd_A1 .* hrc;

%----------------------------------------------------------------------------
% Filter #5 designed by applying FSM to the Wiener deconvolution filter in
% frequency domain. It has the same group delay as the direct FSM but with
% magnitude optimized for noise suppression. First, the two-fold
% oversampled Wiener deconvolution filter is obtained as the LMMSE solution
% of estimating X(w), the pulseshaped original signal spectrum with
% two-fold oversampling. The linear model is written as Y(w) = H(w)X(w) +
% W(w). Note that some system models consider a WGN channel noise W with a
% scaled identity covariance matrix, however, some models consider a WGN
% channel noise W filtered by the receiver filter and hence with a general
% diagonal covariance matrix.
%----------------------------------------------------------------------------
G = frequency_response(N, 2, a, 1, 'rc'); G = transpose(G);
W = [0 : Nf, -Nf : -1] * 2 * pi / N;
% V1: assuming WGN
% H = exp(1i * K * W.^2) .* G.^2 ./ (G.^2 + eta);
% V2: assuming RRC filtered Gaussian noise
H = exp(1i * K * W.^2) .* G ./ (G + eta);

% Obtain the N points time domain impulse response, the DTFT of which
% matches the desired DTFT at N points exactly. Due to time aliasing, an
% overall MSE between the DTFTs cannot be claimed

hcd_A5 = fftshift(ifft(H));

% This is can be formulated as a linear system problem...

% D = dftmtx(N);
% hcd_A5 = fftshift(transpose(D \ H(:)));

%----------------------------------------------------------------------------
% Filter #6 is designed by applying the oversampling FSM to the Wiener
% deconvolution filter in frequency domain. First, the two-fold oversampled
% Wiener deconvolution filter is obtained as the LMMSE solution of
% estimating X(w), the pulseshaped original signal spectrum with two-fold
% oversampling. The linear model is written as Y(w) = H(w)X(w) + W(w). Note
% that some system models consider a WGN channel noise W with a scaled
% identity covariance matrix, however, some models consider a WGN channel
% noise W filtered by the receiver filter and hence with a general diagonal
% covariance matrix.
%----------------------------------------------------------------------------
G = frequency_response(M, 2, a, 1, 'rc'); G = transpose(G);
W = get_fft_grid(M, 2 * pi);
% V1: assuming WGN
% H = exp(1i * K * W.^2) .* G.^2 ./ (G.^2 + eta);
% V2: assuming RRC filtered Gaussian noise
H = exp(1i * K * W.^2) .* G ./ (G + eta);

% Obtain the M (M > N) points time domain impulse response followed by
% an N points truncation in the center. The large M ensures little time
% aliasing in time domain and therefore a rectangular truncation gives the
% minimum MSE in the resulting DTFT, according to the convergence theorem
% of Fourier series

hcd_A6 = ifft(H);
hcd_A6 = fftshift([hcd_A6(1 : Nf + isodd), hcd_A6(end - Nf + isodd : end)]);

% This can be formulated as a overdetermined system problem, i.e., design
% an N points time domain impulse response h[n] such that the DTFT of h[n]
% matches the desired DTFT at M (M > N) points in the least squares sense

% D = dftmtx(M);
% D = [D(:, 1 : Nf + isodd), D(:, end - Nf + isodd : end)];
% hcd_A6 = (D' * D + epsilon * eye(N)) \ D' * H(:);
% hcd_A6 = fftshift(transpose(hcd_A6));

%----------------------------------------------------------------------------
% Eghbali, Amir, et al. "Optimal least-squares FIR digital filters for
% compensation of chromatic dispersion in digital coherent optical
% receivers." Journal of Lightwave Technology 32.8 (2014): 1449-1456.
%----------------------------------------------------------------------------
% erfcmp = @(x) double(erf(sym(x)));
% hcd_A6 = ones(1, N);
% for n = -Nf : Nf
%     hcd_A6(n + Nf + 1) = exp(-1i * (n^2 / (4*K) + 3*pi/4)) / (4*sqrt(pi*K)) ...
%         * (erfcmp(exp(1i*3*pi/4) * (2*K*pi-n) / (2*sqrt(K))) ...
%         + erfcmp(exp(1i*3*pi/4) * (2*K*pi+n) / (2*sqrt(K))));
% end

%----------------------------------------------------------------------------
% Group them together
%----------------------------------------------------------------------------
hcd = [hcd_A1; hcd_A2; hcd_A3; hcd_A4; hcd_A5; hcd_A6];

%----------------------------------------------------------------------------
% Plot the impulse response
%----------------------------------------------------------------------------
if PLOT
    figure;
    stem(-Nf : Nf, transpose(abs(hcd))); grid on;
    xlabel('Time index'); ylabel('Amplitude');
end

%----------------------------------------------------------------------------
% Plot the frequency response
%----------------------------------------------------------------------------
if PLOT
    figure('color', 'w');
    for ii = 1 : size(hcd, 1)
        hh = fftshift(hcd(ii, :));
        hh = [hh(1 : Nf), zeros(1, M - Nf - Nf), hh(Nf + 2 : end)];
        xx = fftshift(get_fft_grid(M, 2 * pi) / pi);
        yy = fftshift(abs(fft(hh)));
        semilogy(xx, yy, 'linewidth', 2); hold on;
    end
    xlim([0, 1]); grid on;
    xlabel('Normalized frequency (\times \pi)'); ylabel('Magnitude');
end

%----------------------------------------------------------------------------
% Raised cosine function
%----------------------------------------------------------------------------
function h = rcshape(Nf, a)
h = zeros(1, Nf + Nf + 1);
for n = -Nf : Nf
    if abs(n) < (1-a) / 2 * Nf
        h(n + Nf + 1) = 1;
    elseif abs(n) > (1+a) / 2 * Nf
        h(n + Nf + 1) = 0;
    else
        h(n + Nf + 1) = 0.5 * (1 + cos(pi / (a * Nf) * (abs(n) - (1-a) / 2 * Nf)));
    end
end