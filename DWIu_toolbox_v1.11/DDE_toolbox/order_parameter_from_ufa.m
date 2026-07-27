function op_map = order_parameter_from_ufa(ufa, fa, background_value)
% Computes the Order Parameter (OP) directly from uFA and FA.
%
%     Parameters
%     ----------
%     ufa : array (Nx, Ny, Nz)
%         Microscopic fractional anisotropy map.
%     fa : array (Nx, Ny, Nz)
%         Macroscopic fractional anisotropy map.
%     background_value : float, optional
%         The value assigned to OP where uFA is near zero. Default = 0.0
%
%     Returns
%     -------
%     op_map : array (Nx, Ny, Nz)
%         OP map strictly bounded between [0, 1].

if nargin < 3
    background_value = 0.0;
end

% Clean FA and uFA maps to physically valid ranges [0, 1]
fa_clean = max(min(fa, 1.0), 0.0);
ufa_clean = max(min(ufa, 1.0), 0.0);

% Compute numerator and denominator of the rearranged OP equation
numerator = (3.0 - 2.0 * (ufa_clean .^ 2)) .* (fa_clean .^ 2);
denominator = (ufa_clean .^ 2) .* (3.0 - 2.0 * (fa_clean .^ 2));

% Create a mask to avoid dividing by zero where uFA is near zero 
% (It is physically impossible to calculate dispersion of isotropic spheres)
valid_mask = ufa_clean > 1e-6;

% Initialize OP with the background value
op_map = background_value * ones(size(ufa));

% Safely compute the division
op_sq = zeros(size(ufa));
op_sq(valid_mask) = numerator(valid_mask) ./ denominator(valid_mask);

% Apply the square root strictly where the mask is valid
op_map(valid_mask) = sqrt(op_sq(valid_mask));

% Final safety clip to ensure OP does not exceed 1.0 
% (Happens if noise causes macroscopic FA to artificially exceed microscopic FA)
op_map = max(min(op_map, 1.0), 0.0);

end