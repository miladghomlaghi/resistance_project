function normalizedData = quantileNormalization(data)
    % Replace NaN with a large negative value to ignore in sorting
    nan_mask = isnan(data);
    data(nan_mask) = -Inf;

    % Sort each column and compute the mean rank
    [sortedData, sortIdx] = sort(data, 1, 'ComparisonMethod', 'auto');
    rowMeans = nanmean(sortedData, 2);

    % Reconstruct normalized data by replacing sorted values with the mean rank
    normalizedData = NaN(size(data));
    for col = 1:size(data, 2)
        % Map the mean rank back to the original positions
        normalizedData(sortIdx(:, col), col) = rowMeans;
    end

    % Restore NaN to the original missing positions
    normalizedData(nan_mask) = NaN;
end
