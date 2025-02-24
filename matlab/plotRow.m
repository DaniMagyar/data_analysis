% Nested anonymous function for plotting a row
    function plotRow(j, viewIndex)
        rowIdx = (viewIndex - 1) * rowsPerView + j;
        axes(ax(j));
        if rowIdx <= numRows
            plot(timeVector, dataMatrix1(rowIdx, :), 'b');  % Plot data from dataMatrix1
            hold on;
            plot(timeVector, dataMatrix2(rowIdx, :), 'r');  % Plot data from dataMatrix2
            hold off;
            ylabel(['Row ' num2str(rowIdx)]);
        else
            cla(ax(j));  % Clear axis if there's no corresponding data
        end
    end