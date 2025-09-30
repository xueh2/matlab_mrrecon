function psnr_value = compute_psnr(a, b)

diff = a-b;
res = diff .* conj(diff);
mse = mean(res(:));
psnr_value = 10 * log(2048^2/mse^2);

