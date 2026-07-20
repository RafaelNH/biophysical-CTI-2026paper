function op_map = order_parameter(kani, fa, background_value)
% Computes the Order Parameter (OP) from K_aniso and FA.
%
%     Parameters
%     ----------
%     kani : array (Nx, Ny, Nz)
%         Anisotropic kurtosis map.
%     fa : array (Nx, Ny, Nz)
%         Fractional anisotropy map.
%     background_value : float, optional
%         The value assigned to OP where K_aniso is <= 0. Default = 0.0
%
%     Returns
%     -------
%     op_map : array (Nx, Ny, Nz)
%         OP map strictly bounded between [0, 1].

if nargin < 3
    background_value = 0.0;
end

% Clean FA and Kani maps
fa_clean = max(min(fa, 1.0), 0.0);
kani_clean = max(kani, 0.0);

% Compute numerator and denominator of OP equation
numerator = 2.4 * (fa_clean .^ 2);
denominator = kani_clean .* (3.0 - 2.0 * (fa_clean .^ 2));

% Create a mask to avoid dividing by zero where K_aniso is near zero
valid_mask = kani_clean > 1e-8;

% Initialize OP with the background value
op_map = background_value * ones(size(kani));

% Safely compute the division
op_sq = zeros(size(kani));
op_sq(valid_mask) = numerator(valid_mask) ./ denominator(valid_mask);

% Apply the square root strictly where the mask is valid
op_map(valid_mask) = sqrt(op_sq(valid_mask));

% Final safety clip to ensure OP does not exceed 1.0 
% (Happens if noise causes macroscopic FA to artificially exceed microscopic FA)
op_map = max(min(op_map, 1.0), 0.0);

end