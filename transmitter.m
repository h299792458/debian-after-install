% QAM transmitter

baudrate = 10.7e9;
bitpersym = 4;
sps = 2;
fs = sps * baudrate;
nsym = 65536 * 4;
nsap = sps * nsym;
ALPHABET_SIZE = 2^bitpersym;
wavelength = 1553e-9;
fiber.DL = 17e-6 * 4000e3;
fiber.SL = 0;

G = frequency_response(nsap, sps, 0.22, 1, 'rrc');

rawbits = randi([0 1], bitpersym, nsym);
rawsyms = mqammod(rawbits) / (sqrt(ALPHABET_SIZE) - 1);
x_inph = real(ifft(fft(upsample(real(rawsyms), sps)) .* G));
x_quad = real(ifft(fft(upsample(imag(rawsyms), sps)) .* G));

x = x_inph + 1i * x_quad;

ts = 1 / fs;

K = - fiber.DL * wavelength * wavelength / (4 * pi * 299792458 * ts * ts);
W = [0 : 1 : nsap/2-1, -nsap/2 : 1 : -1] * 2 * pi / nsap;
H = exp(1i * K * W.^2);

x = ifft(fft(x) .* transpose(H));