clear

transmitter;

SNR = 10;
sigma2 = idbw(dbw(calcrms(x).^2) - SNR);
x_full = x + gaussian_noise(nsap, 1, sigma2, 'linear', 'complex');

ts = 1 / (2*baudrate);
K = fiber.DL * wavelength * wavelength / (4 * pi * 299792458 * ts * ts);
N = 101 : 2 : 201;

for ii = 1 : length(N)
    hcd = design_hcd(K, N(ii));
    ber(ii,:) = receiver(sps, G, hcd, x_full, rawbits)
end
BER = snr2ber(SNR + dbw(sps), bitpersym, 'db');
figure; semilogy(N, BER*ones(1,ii), 'k', N, ber); grid on; xlim([min(N), max(N)]);

keyboard;