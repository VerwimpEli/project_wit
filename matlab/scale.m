function S = scale(X)
    m = sqrt(size(X, 1));
    S = reshape(X, m, m);
    S(:, 1)   = 1/2 * S(:, 1);
    S(:, end) = 1/2 * S(:, end);
    S(1, :)   = 1/2 * S(1, :);
    S(end, :) = 1/2 * S(end, :);
    S = reshape(S, size(X));
end