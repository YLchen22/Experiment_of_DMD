function matcompare(mats, varargin)
    p = inputParser;
    addRequired(p, 'mats');
    addParameter(p, 'fixmat', 0);
    parse(p, mats, varargin{:});
    fixmat = p.Results.fixmat;

    if fixmat
        M = mats{fixmat};
    else
        M = cell2mat(mats);
    end
    n = length(mats);
    minvalue = min(M, [], "all");
    maxvalue = max(M, [], "all");
    figure()
    for i = 1:n
        subplot(1, n, i)
        mat = (mats{i} - minvalue) / (maxvalue - minvalue);
        imshow(mat, 'Colormap', cool, 'InitialMagnification', 'fit')
        title([num2str(i) '-th matrix'])
    end
end