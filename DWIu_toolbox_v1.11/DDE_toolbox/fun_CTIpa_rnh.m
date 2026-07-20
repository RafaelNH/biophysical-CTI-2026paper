function [Dtotal, Ktotal, Kaniso, Kiso, Kintra, S0, sum_res2] = fun_CTIpa_rnh(data, gtab, ...
    mask, norm_data, rep, nls)
% Mean Signal correlation tensor imaging
%
%     Parameters
%     ----------
%     data: array (x, y, z, ndwis), (x, y, ndwis), (x, ndwis) or (ndwis)
%         Array containing the diffusion-weighted signals, note that last
%         dimension corresponds to different diffusion-weigthed experiment
%         which may have been acquired for different gradient direction
%         or b-value
%
%     gtab: GradientTable
%        A GradientTable with all the gradient information
%        (see read/read_diffusion)
%
%     mask : array (..., )
%         Array containing true values for voxels to be processed
%
%     norm_data : bool
%         If true, data is assumed to be normalized by S0 and S0 is not
%         computed
%
%     rep : int
%         Number of iterations to compute the weights LLS fit
%         if 0 computes un-weighted LLS fit. Defaut 0
%
%     nls : bool
%         if true non-linear least square fit is computed after weighted
%         linear least square solution
%         Defaut: false
%
%     Returns
%     -------
%     Dtotal : array (..., )
%         Diffusivity of the mean signal
%     Ktotal : array (..., )
%         Total Kurtosis of the mean signal
%     Kaniso : array (..., )
%         Anisotropic Kurtosis Source
%     Kiso : array (..., )
%         Isotropic Kurtosis Source
%     Kintra : array (..., )
%         Intracompartmental Kurtosis Source
%     S0 : array (..., )
%         Estimate of signal at b-value=0
%
%     References
%     ----------
%    Henriques, R., Jespersen, S., Shemesh, N. (2021). Evidence for
%    microscopic kurtosis in neural tissue revealed by Correlation Tensor
%    MRI . MRM

% prepare
if ~exist('norm_data', 'var')
    norm_data = false;
end
if ~exist('rep', 'var')
    rep = 0;
end
if ~exist('nls', 'var')
    nls = false;
end

b1 = gtab.bval1;
b2 = gtab.bval2;

cos2 = gtab.cos2;

SIZ = size(data);
Nvol = SIZ(end);
Nvox = prod(SIZ(1:end-1));

data = reshape(data, Nvox, Nvol);

Z0 = zeros(Nvox, 1);
Dtotal = Z0;
Ktotal = Z0;
Kaniso = Z0;
Kiso = Z0;
Kintra = Z0;
S0 = Z0;
sum_res2 = Z0;


% design matrix
if norm_data
    A = [-(b1+b2), ((b1.^2+b2.^2))/6, (b1.*b2.*cos2)/2, (b1.*b2)/6];
    lowc = zeros(4, 1);
    uppc = [3; 3; 3; 3];
else
    A = [-(b1+b2), ((b1.^2+b2.^2))/6, (b1.*b2.*cos2)/2, (b1.*b2)/6, ones(Nvol, 1)];
    lowc = zeros(5, 1);
    uppc = [3; 3; 3; 3; inf];
end

for v=1:size(data,1)
    if(mask(v)==1)
        %B
        S = double(squeeze(data(v, :)));
        Sc = S;
        Sc(Sc<0.0001) = 0.0001;
        B=log(Sc(:));
        
        %ULS
        piA=pinv(A); %piA pseudoinverse of A
        X=piA*B;
        
        for r=1:rep
            Spred = exp(diag(A * X));
            X = pinv(A' * Spred^2 * A) * A' * Spred^2 * B;
        end
        
        D = X(1);
        K = X(2)/(D^2);
        Ka = X(3)/(D^2);
        KiKa = X(4)/(D^2);
        Ki = (KiKa + Ka) / 2;
        
        if nls
            Xnls = X;
            Xnls(2) = K;
            Xnls(3) = Ka;
            Xnls(4) = Ki;
            
            Xnls(Xnls<lowc) = lowc(Xnls<lowc);
            Xnls(Xnls>uppc) = uppc(Xnls>uppc);
            
            ydata = S(:);
            optionsF = optimset('Display', 'off');
            Xnls =  lsqcurvefit(@fun_pacti, Xnls, A, ydata, lowc, uppc,...
                optionsF);
            
            X = Xnls;
            D = Xnls(1);
            K = Xnls(2);
            Ka = Xnls(3);
            Ki = Xnls(4);
            
            X(4) = 2*X(4) - X(3);
            X(2:4) = X(2:4) * (D^2);
        end
        
        Dtotal(v) = D;
        Ktotal(v) = K;
        Kaniso(v) = Ka;
        Kiso(v) = Ki;
        Kintra(v) = K - Ka - Ki;
        
        if norm_data
            S0(v) = 1;
        else
            S0(v) = exp(X(5));
        end
        
        Spred = exp(A * X);
        residual = sum((Spred(:) - S(:)).^2);
        stotal = sum((S(:) - mean(S(:))).^2);
        sum_res2(v) = 1 - residual/stotal;
        
    end
    %disp(v)
end

Dtotal = reshape(Dtotal, SIZ(1:end-1));
Ktotal = reshape(Ktotal, SIZ(1:end-1));
Kaniso  = reshape(Kaniso, SIZ(1:end-1));
Kiso  = reshape(Kiso, SIZ(1:end-1));
Kintra  = reshape(Kintra, SIZ(1:end-1));
S0  = reshape(S0, SIZ(1:end-1));
sum_res2  = reshape(sum_res2, SIZ(1:end-1));
end

function S = fun_pacti(Xnls, A)
% convert Xls to X
X = Xnls;
D = X(1);
X(4) = 2*X(4) - X(3);
X(2:4) = X(2:4) * (D^2);

S = exp(A*X);
end
