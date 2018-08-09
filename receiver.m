function ber = receiver(sps, G, hcd, x, rawbits)

x_inph = real(ifft(fft(real(x)).* G));
x_quad = real(ifft(fft(imag(x)).* G));

% using two-fold oversampling
L = 2;
x = x_inph(1:sps/L:end) + 1i * x_quad(1:sps/L:end);

ALPHABET_SIZE = 16;

N = 2 * size(hcd, 2);

x_fltd = zeros(size(x));
x_hat  = zeros(size(x(1:2:end)));
decbits = cell(1, size(hcd, 1));
ber = zeros(1, size(hcd, 1));
for ii = 1 : size(hcd, 1)
    hcd_A = hcd(ii, :);
    x_fltd(:, ii) = conv(x, hcd_A, 'same');
    x_hat(:, ii) = x_fltd(1:2:end, ii);
    decbits{ii} = qamdemod(x_hat(:, ii), ALPHABET_SIZE);
    ber(ii) = nnz(decbits{ii}(:, N+1 : end-N) - rawbits(:, N+1 : end-N)) / numel(rawbits(:, N+1 : end-N));
end
