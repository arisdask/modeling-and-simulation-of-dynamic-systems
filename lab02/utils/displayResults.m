%% Function to display parameter estimation results
function displayResults(m_true, b_true, k_true, m_est, b_est, k_est, case_name)
    fprintf('\n%s:\n', case_name);
    fprintf('True parameters: m = %.3f, b = %.3f, k = %.3f\n', m_true, b_true, k_true);
    fprintf('Final estimates: m = %.3f, b = %.3f, k = %.3f\n', m_est, b_est, k_est);
    fprintf('Relative errors: m = %.2f%%, b = %.2f%%, k = %.2f%%\n', ...
        100*abs(m_est-m_true)/m_true, ...
        100*abs(b_est-b_true)/b_true, ...
        100*abs(k_est-k_true)/k_true);
end
