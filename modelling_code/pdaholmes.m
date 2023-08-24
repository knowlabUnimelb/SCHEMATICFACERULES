function lik = pdaholmes(varargin)
% Generic function to simulate from a given logical rule model
% data is assumed to be [resp, rt] or [resp, rt, item]

%% Variable input parameters
optargs = {@simlba, 1e5, [], [], '', 10};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[model, n, data, parms, stoppingrule, maxrtCutoff] = optargs{:}; % Place optional args in memorable variable names

%%
N_grid=2^10; % Initialize the number of grid points to bin the data into. Note this is a VERY fine grid

%% Simulate from model
[rt, resp] = model(n, parms, stoppingrule);
resp(rt > maxrtCutoff) = [];
rt(rt > maxrtCutoff) = [];

rts1 = rt(resp == 1); % target
rts2 = rt(resp == 2); % contrast

%%
if numel(rts1) > 1
    % Initialize a histogram grid and construct a filter
    % Comute the min and max of the input data for binning purposes. Pad the
    % data to the left and right to satisfy the periodic wrap around condition
    % needed to apply a FFT.
    h1 = getbandwidth(rts1); % Get bandwidths for each alternative (Turner2014, Eq. 14)
    m1 = min(data(:,2))-3 * h1;
    M1 = max(data(:,2))+3 * h1;
    
    % Setup the histogram bins
    grid_centers1 = linspace(m1,M1,N_grid);
    dt1=grid_centers1(2)-grid_centers1(1);
    bin_edges1=[grid_centers1-dt1/2 , grid_centers1(end)+dt1/2];
    
    % Construct the Gaussian filter that will be used to filter the data.
    freq1=2*pi*N_grid/(M1-m1)*0.5*linspace(0,1,N_grid/2+1);
    filter1=exp(-0.5*freq1.^2*h1^2);
    filter1=[filter1 , fliplr(filter1(2:end-1))];
    filter1=filter1';
    
    % Transform, smooth, and transform back to compute the PDF
    % Create a discrete histogram representation of the PDF
%     PDF_hist1 = 1*histc(rts1,bin_edges1)/(dt1*numel(rts1));
    PDF_hist1 = 1*histc(rts1,bin_edges1)/(dt1*numel(rt));
    
    % Transform into the spectral domain.
    PDF_fft1=fft(PDF_hist1(1:end-1));
    
    % Apply the smoothing filter. This is effectively the convolution step in
    % the Fourier domain.
    PDF_fft_filt1=filter1.*PDF_fft1;
    
    % Transform back into data space. This produces s smoothed PDF on the very
    % fine grid.
    PDF_smoothed1=ifft(PDF_fft_filt1);
    
    % Linearly interpolate the smoothed PDF to the data values to obtain an
    % approximation of the likelihood of each observation.
    PDF1=real(interp1(grid_centers1,PDF_smoothed1,data(data(:,1) == 1, 2)));
    
    % Replace any very small elements with a minimum value of the PDF. Usually
    % taken to be ~ 1/Nsample. This is required for taking the subsequent log
    % in computing the log likelihood.
    PDF1= max(PDF1,10^(-10));
    
    lik(data(:,1) == 1,1) = PDF1;
    
else
    lik(data(:,1) == 1,1) = 10^-10;
end

%%
if numel(rts2) > 1
    h2 = getbandwidth(rts2); % Get bandwidths for each alternative (Turner2014, Eq. 14)
    m2 = min(data(:,2))-3 * h2;
    M2 = max(data(:,2))+3 * h2;
    
    % Setup the histogram bins
    grid_centers2 = linspace(m2,M2,N_grid);
    dt2=grid_centers2(2)-grid_centers2(1);
    bin_edges2=[grid_centers2-dt2/2 , grid_centers2(end)+dt2/2];
    
    % Construct the Gaussian filter that will be used to filter the data.
    freq2=2*pi*N_grid/(M2-m2)*0.5*linspace(0,1,N_grid/2+1);
    filter2=exp(-0.5*freq2.^2*h2^2);
    filter2=[filter2 , fliplr(filter2(2:end-1))];
    filter2=filter2';
    
    % Create a discrete histogram representation of the PDF
%     PDF_hist2 = 1*histc(rts2,bin_edges2)/(dt2*numel(rts2));
    PDF_hist2 = 1*histc(rts2,bin_edges2)/(dt2*numel(rt));
    
    % Transform into the spectral domain.
    PDF_fft2=fft(PDF_hist2(1:end-1));
    
    % Apply the smoothing filter. This is effectively the convolution step in
    % the Fourier domain.
    PDF_fft_filt2=filter2.*PDF_fft2;
    
    % Transform back into data space. This produces s smoothed PDF on the very
    % fine grid.
    PDF_smoothed2=ifft(PDF_fft_filt2);
    
    % Linearly interpolate the smoothed PDF to the data values to obtain an
    % approximation of the likelihood of each observation.
    PDF2=real(interp1(grid_centers2,PDF_smoothed2,data(data(:,1) == 2, 2)));
    
    % Replace any very small elements with a minimum value of the PDF. Usually
    % taken to be ~ 1/Nsample. This is required for taking the subsequent log
    % in computing the log likelihood.
    PDF2= max(PDF2,10^(-10));
    
    lik(data(:,1) == 2,1) = PDF2;
else
    lik(data(:,1) == 2, 1) = 10^-10;
end

% lik = sum(log(PDF1))+sum(log(PDF2)); % Correct rt + error rt