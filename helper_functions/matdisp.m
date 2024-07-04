function matdisp(mats)
    n = length(mats);
    figure()
    for i = 1:n
        subplot(1, n, i)
        mat = mats{i};
        imagesc(mat)
        colorbar()
        title([num2str(i) '-th matrix'])
    end
end