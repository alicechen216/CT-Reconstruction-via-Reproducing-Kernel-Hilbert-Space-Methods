% generate_phantoms.m
function img = generate_phantoms(name, img_size)
    % Generates phantom images as described in Fig. 1 [cite: 485]
    switch lower(name)
        case 'shepp-logan'
            img = phantom(img_size);
        case 'bullseye'
            [X, Y] = meshgrid(linspace(-1, 1, img_size));
            R = sqrt(X.^2 + Y.^2);
            img = sin(10 * pi * R);
            img(R > 1) = 0;
        case 'crescent'
             [X, Y] = meshgrid(linspace(-1, 1, img_size));
             c1 = (X-0.2).^2 + Y.^2 < 0.7^2;
             c2 = X.^2 + Y.^2 < 0.5^2;
             img = double(c1 & ~c2);
        case 'bubbles'
            img = zeros(img_size);
            [X, Y] = meshgrid(linspace(-1, 1, img_size));
            centers = (rand(10, 2) - 0.5) * 1.5;
            radii = rand(10, 1) * 0.15 + 0.05;
            for i = 1:10
                img = img + double((X-centers(i,1)).^2 + (Y-centers(i,2)).^2 < radii(i)^2);
            end
        otherwise
            error('Unknown phantom name.');
    end
    img = rescale(img, 0, 255);
end