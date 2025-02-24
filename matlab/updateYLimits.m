% Function to update y-axis limits based on slider value
function updateYLimits(sliderValue)
    newYLim = [sliderValue, sliderValue + (initialYLim(2) - initialYLim(1))];
    set(findall(gcf, 'Type', 'axes'), 'YLim', newYLim);
end