function [num_field, res_map] = field_residual(vals, featx, featy)
%FIELD_RESIDUAL: compute a residual map for visualization
%   vals: solved lambdas for reference
%   featx, featy: data projected to the feature space
%   RETURN:
%   num_field: complex number range to plot (same for real and imag, use it to meshgrid)
%   res_map: residual map corresponding to the field to plot

%   use output:
%   [re_map, im_map] = meshgrid(num_field);
%   pcolor(re_map, im_map, res_map)

%   weight matrix W is not included (considered as identity)

    min_num = min([real(vals); imag(vals); -1.1], [], 'all');
    max_num = max([real(vals); imag(vals); 1.1], [], 'all');

    reso = 41;
    num_field = linspace(min_num, max_num, reso);

    [re_map, im_map] = meshgrid(num_field);
    res_map = zeros(reso, reso);

    for x = 1: reso
        for y = 1: reso
            lambda = re_map(x, y) + 1i*im_map(x, y);
            [res, ~] = mini_residual(lambda, featx, featy);
            res_map(x, y) = res;
        end
        disp(['Solving residual map: ', num2str(x), ' / ', num2str(reso)])
    end

end

