clear all; close all; clc;

NSNR = 12;
Ncw = 10000;
EbN0dB = linspace(0,10,NSNR);
EbN0lin = 10.^(EbN0dB/10);
Eb = 1; 

err_unc = zeros(1,NSNR); 
err_rep = zeros(1,NSNR);
err_cro = zeros(1,NSNR);

% Define the generator matrix G for the (9,4) cross-parity-check code
G = [1 0 0 0 1 0 1 0 1;
     0 1 0 0 1 1 0 1 0;
     0 0 1 0 0 1 1 1 0;
     0 0 0 1 0 0 1 1 1];

for SNRidx = 1:NSNR % loop over SNR values

    fprintf("-------------------------------------------\n");
    fprintf("SNR iteration %d / %d, Eb/N0 = %5.2f dB \n", SNRidx, NSNR, EbN0dB(SNRidx));
    fprintf("-------------------------------------------\n");
    EbN0curr = EbN0lin(SNRidx); % pick current SNR value from the vector
    N0 = Eb/EbN0curr;
    sigma_n2 = N0/2;

    %theoretical BER for the uncoded case
    BER_the(SNRidx) = 1/2 * erfc(sqrt(EbN0curr));
    %BER_the_cro(SNRidx) = 4/4 * Q(sqrt(2*4*4/9*EbN0curr)); % approximation

    for CWidx = 1:Ncw

      % 1. uncoded transmission
      u_unc = round(rand(1,1));
      Ec_unc = Eb;
      a_unc = ((u_unc*2)-1) * sqrt(Ec_unc);
      n_unc = sqrt(sigma_n2)*randn(1,1);
      y_unc = a_unc + n_unc;
      a_hat_unc = sign(y_unc);
      u_hat_unc = (a_hat_unc+1)/2;

      if u_unc ~= u_hat_unc % then an error occured
         err_unc(SNRidx) = err_unc(SNRidx) + 1;
      end

      if CWidx == 1
         fprintf("Uncoded case: \n");
         fprintf("   Info word:             u     = [%d] <----\n",u_unc);
         fprintf("   Transmitted code word: a     = [%4.1f] \n", a_unc);
         fprintf("   Noise:                 n     = [%4.1f] \n", n_unc);
         fprintf("   Received word:         y     = [%4.1f] \n", y_unc);
         fprintf("   Estimated code word:   a_hat = [%d] \n", a_hat_unc);
         fprintf("   Estimate info word:    u_hat = [%d] <----\n", u_hat_unc);
      end

      % 2. (3,1) repetition code
      u_rep = round(rand(1,1));
      Ec_rep = 1/3 * Eb;
      a_rep = ((u_rep*2)-1) * sqrt(Ec_rep) .* [1 1 1];
      n_rep = sqrt(sigma_n2)*randn(1,3);
      y_rep = a_rep + n_rep;
      d1 = sum((y_rep - (sqrt(Ec_rep) .* [1 1 1])).^2);
      d2 = sum((y_rep - (sqrt(Ec_rep) .* [-1 -1 -1])).^2);
      if d1<d2
        a_hat_rep = sqrt(Ec_rep)*[1 1 1]; u_hat_rep = 1;
      else
        a_hat_rep = sqrt(Ec_rep)*[-1 -1 -1]; u_hat_rep = 0;
      end
      if u_rep ~= u_hat_rep % then an error occured
         err_rep(SNRidx) = err_rep(SNRidx) + 1;
      end

      if CWidx == 1
         fprintf("\n(3,1) Repetition code: \n");
         fprintf("   Info word:             u     = [%d] <----\n",u_rep);
         fprintf("   Transmitted code word: a     = [%4.1f %4.1f %4.1f] \n", a_rep(1), a_rep(2), a_rep(3));
         fprintf("   Noise:                 n     = [%4.1f %4.1f %4.1f] \n", n_rep(1), n_rep(2), n_rep(3));
         fprintf("   Received word:         y     = [%4.1f %4.1f %4.1f] \n", y_rep(1), y_rep(2), y_rep(3));
         fprintf("   Euklid. dist to a_1:   d1    = %5.2f \n",d1);
         fprintf("   Euklid. dist to a_2:   d2    = %5.2f \n",d2);
         fprintf("   Estimated code word:   a_hat = [%4.1f %4.1f %4.1f] \n", a_hat_rep(1), a_hat_rep(2), a_hat_rep(3));
         fprintf("   Estimated info word:   u_hat = [%d] <----\n",u_hat_rep);
      end

      % 3. (9,4) cross-parity check code
      u_cro = round(rand(1,4));
      a_cro = mod(u_cro * G, 2);
      Ec_cro = 4/9 * Eb;
      a_cro = ((a_cro*2)-1) * sqrt(Ec_cro);
      n_cro = sqrt(sigma_n2)*randn(1,9);
      y_cro = a_cro + n_cro;
      distances = zeros(1,16);
      for i = 1:16
          u_temp = de2bi(i-1, 4, 'left-msb');
          a_temp = mod(u_temp * G, 2);
          a_temp = ((a_temp*2)-1) * sqrt(Ec_cro);
          distances(i) = sum((y_cro - a_temp).^2);
      end
      % Find the codeword with the minimum distance
      [~, min_idx] = min(distances);
      u_hat_cro = de2bi(min_idx-1, 4, 'left-msb');
      % Count errors
      if any(u_cro ~= u_hat_cro)
          err_cro(SNRidx) = err_cro(SNRidx) + 1;
      end

      if CWidx == 1
         fprintf("\n(9,4) Cross-parity check code: \n");
         fprintf("   Info word:             u     = [%d %d %d %d] <----\n",u_cro(1), u_cro(2), u_cro(3), u_cro(4));
         fprintf("   Transmitted code word: a     = [%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f] \n", a_cro(1), a_cro(2), a_cro(3), a_cro(4), a_cro(5), a_cro(6), a_cro(7), a_cro(8), a_cro(9));
         fprintf("   Noise:                 n     = [%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f] \n", n_cro(1), n_cro(2), n_cro(3), n_cro(4), n_cro(5), n_cro(6), n_cro(7), n_cro(8), n_cro(9));
         fprintf("   Received word:         y     = [%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f] \n", y_cro(1), y_cro(2), y_cro(3), y_cro(4), y_cro(5), y_cro(6), y_cro(7), y_cro(8), y_cro(9));
         fprintf("   Estimated info word:   u_hat = [%d %d %d %d] <----\n",u_hat_cro(1), u_hat_cro(2), u_hat_cro(3), u_hat_cro(4));
      end

    end
end

BER_unc = err_unc / Ncw; % get bit error ratio from error count
BER_rep = err_rep / Ncw;
BER_cro = err_cro / (Ncw * 4); % Normalize by the number of info bits per codeword

% Plotting results
figure;
semilogy(EbN0dB, BER_the, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0dB, BER_unc, 'r--', 'LineWidth', 2);
semilogy(EbN0dB, BER_rep, 'g-.', 'LineWidth', 2);
semilogy(EbN0dB, BER_cro, 'm:', 'LineWidth', 2);
grid on;
set(gca, 'FontSize', 18);
xlabel('E_b / N_0 in dB');
ylabel('Bit Error Rate (BER)');
title('BER Comparison');
legend('Uncoded, Theory', 'Uncoded, Simulation', '(3,1) Rep.-Code, Soft-Dec.', '(9,4) Cross-Parity, Soft-Dec.', 'Location', 'southwest');