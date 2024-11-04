% ez barmosag, mindig 1.96 zscore a p=0.05 es 2.58 a p=0.01



% Assuming your z-scores are stored in a vector called zscores
zscores = response_zscore1(1:3,:);
% Assuming your z-scores are stored in a vector called zscores

% Plot the distribution of zscores
figure;
histogram(zscores, 3000, 'Normalization', 'pdf'); % 30 bins, normalized as a PDF
title('Z-score Distribution');
xlabel('Z-score');
ylabel('Probability Density');

% Find the z-score values for p = 0.05 on both sides
p_value_lower = 0.05;  % Lower tail
p_value_upper = 1 - p_value_lower;  % Upper tail

z_critical_lower = norminv(p_value_lower);  % Lower critical z-score
z_critical_upper = norminv(p_value_upper);  % Upper critical z-score

% Overlay the critical z-score lines on the plot
hold on;
yLimits = ylim;  % Get the current y-axis limits for plotting the lines

% Plot the lower critical z-score (left tail)
plot([z_critical_lower z_critical_lower], yLimits, 'r--', 'LineWidth', 2);
text(z_critical_lower, yLimits(2)*0.9, sprintf('p = 0.05, z = %.2f', z_critical_lower), 'Color', 'red', 'FontSize', 12);

% Plot the upper critical z-score (right tail)
plot([z_critical_upper z_critical_upper], yLimits, 'r--', 'LineWidth', 2);
text(z_critical_upper, yLimits(2)*0.9, sprintf('p = 0.95, z = %.2f', z_critical_upper), 'Color', 'red', 'FontSize', 12);

hold off;

% Output the critical z-scores
fprintf('The z-score value where p = 0.05 (lower tail) is: %.4f\n', z_critical_lower);
fprintf('The z-score value where p = 0.95 (upper tail) is: %.4f\n', z_critical_upper);