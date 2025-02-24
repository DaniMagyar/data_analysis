
% Function to update x-axis limits based on slider value
function updateXLimits(sliderValue, plotWindow)
    newXLimit = [sliderValue, sliderValue + (plotWindow(2) - plotWindow(1))];
    set(findall(gcf, 'Type', 'axes'), 'XLim', newXLimit);
end