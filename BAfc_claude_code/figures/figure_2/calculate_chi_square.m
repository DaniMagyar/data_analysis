function [chi2_stat, p_value] = calculate_chi_square(contingency_table)
    % Calculate chi-square statistic from contingency table
    row_totals = sum(contingency_table, 2);
    col_totals = sum(contingency_table, 1);
    grand_total = sum(contingency_table(:));

    % Expected frequencies
    expected = (row_totals * col_totals) / grand_total;

    % Chi-square statistic (only for cells with expected > 0)
    valid_mask = expected > 0;
    chi2_stat = sum(((contingency_table(valid_mask) - expected(valid_mask)).^2) ./ expected(valid_mask));

    % Degrees of freedom
    df = (size(contingency_table, 1) - 1) * (size(contingency_table, 2) - 1);

    % Parametric p-value
    p_value = 1 - chi2cdf(chi2_stat, df);
end
