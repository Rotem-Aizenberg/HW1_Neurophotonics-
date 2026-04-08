function [dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, ...
    extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
%% CalcNIRS - Compute hemoglobin concentration changes from fNIRS data
%
% Applies the Modified Beer-Lambert Law (MBLL) to dual-wavelength intensity
% measurements, solving for Delta[HbO] and Delta[HbR] at each time point
% across all measurement channels.
%
% Input:
%   dataFile - .mat file containing the fNIRS recording structure:
%              SD.Lambda : two measurement wavelengths [nm]
%              t         : time vector [s]
%              d         : intensity matrix (n x 40), where columns 1-20
%                          are wavelength 1 and columns 21-40 are wavelength 2
%   SDS      - Source-Detector Separation [cm]
%   tissueType - tissue label matching a row in DPFperTissueFile
%                (e.g. 'adult_forearm', 'baby_head', 'adult_head', 'adult_leg')
%   plotChannelIdx - indices of channels to plot, in [1,20]. Empty = no plot.
%                    (default = [])
%   extinctionCoefficientsFile - CSV with columns: wavelength, Water, HbO2,
%                HHb, FatSoybean. Values are specific extinction coefficients,
%                i.e. epsilon/ln(10). (default = '.\ExtinctionCoefficientsData.csv')
%   DPFperTissueFile - TXT with columns: Tissue, DPF (at 807nm reference).
%                (default = '.\DPFperTissue.txt')
%   relDPFfile - CSV with columns: wavelength, relDPFcoeff. Wavelength-
%                dependent scaling of the base DPF.
%                (default = '.\RelativeDPFCoefficients.csv')
%
% Output:
%   dHbR - deoxyhemoglobin concentration change (n x 20)
%   dHbO - oxyhemoglobin concentration change (n x 20)
%   fig  - figure handle, or [] if no channels were plotted

%% Default arguments
if nargin < 4 || isempty(plotChannelIdx),  plotChannelIdx = [];  end
if nargin < 5 || isempty(extinctionCoefficientsFile)
    extinctionCoefficientsFile = 'ExtinctionCoefficientsData.csv';
end
if nargin < 6 || isempty(DPFperTissueFile)
    DPFperTissueFile = 'DPFperTissue.txt';
end
if nargin < 7 || isempty(relDPFfile)
    relDPFfile = 'RelativeDPFCoefficients.csv';
end

%% Input validation
if ~ischar(dataFile) && ~isstring(dataFile)
    error('dataFile must be a string.');
end
if ~endsWith(dataFile, '.mat')
    error('dataFile must be a .mat file.');
end
if ~isfile(dataFile)
    error('dataFile not found: %s', dataFile);
end
if ~isnumeric(SDS) || ~isscalar(SDS) || SDS <= 0
    error('SDS must be a positive scalar.');
end
if (~ischar(tissueType) && ~isstring(tissueType)) || strlength(tissueType) == 0
    error('tissueType must be a non-empty string.');
end
if ~isnumeric(plotChannelIdx)
    error('plotChannelIdx must be a numeric vector.');
end
if any(plotChannelIdx < 1 | plotChannelIdx > 20 | mod(plotChannelIdx,1) ~= 0)
    error('plotChannelIdx values must be integers in [1, 20].');
end
for fPath = {extinctionCoefficientsFile, DPFperTissueFile, relDPFfile}
    if ~isfile(fPath{1})
        error('File not found: %s', fPath{1});
    end
end

%% Load and validate .mat structure
matData = load(dataFile);

for ii = 1:numel({'SD','t','d'})
    fieldList = {'SD','t','d'};
    if ~isfield(matData, fieldList{ii})
        error('Required field "%s" missing from .mat file.', fieldList{ii});
    end
end
if ~isfield(matData.SD, 'Lambda')
    error('SD.Lambda field missing from .mat file.');
end

lambdas = double(matData.SD.Lambda(:)');
if numel(lambdas) ~= 2
    error('Expected exactly 2 wavelengths in SD.Lambda, got %d.', numel(lambdas));
end

timeVec = double(matData.t(:));
intensityData = double(matData.d);
nTimePoints = length(timeVec);
nChannels = 20;

if size(intensityData,1) ~= nTimePoints || size(intensityData,2) ~= 2*nChannels
    error('Data matrix dimensions [%d,%d] inconsistent with expected [%d,%d].', ...
        size(intensityData,1), size(intensityData,2), nTimePoints, 2*nChannels);
end

%% Load specific extinction coefficients
% These are epsilon/ln(10), tabulated per wavelength. Since the OD is
% defined with log10, the ln(10) factor cancels and they can be used directly.
extTable = readtable(extinctionCoefficientsFile, 'CommentStyle', '%');
extinctionCoeffs = zeros(2, 2);   % [wavelength, chromophore(HbO2|HHb)]

for ii = 1:2
    matchRow = extTable(extTable.wavelength == lambdas(ii), :);
    if isempty(matchRow)
        error('Wavelength %d nm not found in extinction coefficients file.', lambdas(ii));
    end
    extinctionCoeffs(ii, 1) = matchRow.HbO2;
    extinctionCoeffs(ii, 2) = matchRow.HHb;
end

%% Load tissue DPF at 807 nm reference wavelength
fid = fopen(DPFperTissueFile, 'r');
baseDPF = [];
while ~feof(fid)
    currentLine = fgetl(fid);
    if isempty(currentLine) || startsWith(strtrim(currentLine), '%')
        continue;
    end
    tokens = strsplit(strtrim(currentLine));
    if strcmp(tokens{1}, tissueType)
        baseDPF = str2double(tokens{2});
        break;
    end
end
fclose(fid);
if isempty(baseDPF)
    error('Tissue type "%s" not found in %s.', tissueType, DPFperTissueFile);
end

%% Load wavelength-dependent relative DPF
% The base DPF is measured at 807 nm; relDPF(lambda) corrects it for the
% wavelength dependence of scattering in tissue.
relDPFTable = readtable(relDPFfile, 'CommentStyle', '%');
relDPF = zeros(1, 2);
for ii = 1:2
    matchRow = relDPFTable(relDPFTable.wavelength == lambdas(ii), :);
    if isempty(matchRow)
        error('Wavelength %d nm not found in relative DPF file.', lambdas(ii));
    end
    relDPF(ii) = matchRow.relDPFcoeff;
end

%% Effective path length
% Due to multiple scattering, photons travel a path much longer than the
% geometric SDS. The DPF captures this ratio: Leff = SDS * DPF(lambda).
effectivePathLength = SDS * baseDPF * relDPF;   % [cm], one per wavelength

%% Optical density
% I use log10 rather than ln because the provided extinction coefficients
% already incorporate the 1/ln(10) factor (specific extinction convention).
% Baseline is the first time sample.
intensityWL1 = intensityData(:, 1:nChannels);
intensityWL2 = intensityData(:, nChannels+1 : end);

baselineWL1 = intensityWL1(1, :);
baselineWL2 = intensityWL2(1, :);

odWL1 = log10(baselineWL1 ./ intensityWL1);   % (nTimePoints x nChannels)
odWL2 = log10(baselineWL2 ./ intensityWL2);

%% MBLL inversion
% For each wavelength:  OD(lam)/Leff(lam) = e_HbO2(lam)*dHbO + e_HHb(lam)*dHbR
% Stacking the two wavelengths gives a 2x2 system [E]*[dHbO; dHbR] = [OD/Leff].
% I invert E once and apply it to every time-channel pair (vectorized).
extinctionMatrix = [extinctionCoeffs(1,1), extinctionCoeffs(1,2);
                    extinctionCoeffs(2,1), extinctionCoeffs(2,2)];

if abs(det(extinctionMatrix)) < 1e-12
    error('Extinction coefficient matrix is near-singular; check wavelength separation.');
end

invE = inv(extinctionMatrix);

odNormWL1 = odWL1 / effectivePathLength(1);
odNormWL2 = odWL2 / effectivePathLength(2);

dHbO = invE(1,1)*odNormWL1 + invE(1,2)*odNormWL2;   % (nTimePoints x nChannels)
dHbR = invE(2,1)*odNormWL1 + invE(2,2)*odNormWL2;

%% Plot
fig = [];
if ~isempty(plotChannelIdx)
    fig = figure;
    numSubplots = length(plotChannelIdx);
    for ii = 1:numSubplots
        chIdx = plotChannelIdx(ii);
        subplot(numSubplots, 1, ii);
        plot(timeVec, dHbR(:, chIdx), 'b', 'DisplayName', '\DeltaHbR');
        hold on;
        plot(timeVec, dHbO(:, chIdx), 'r', 'DisplayName', '\DeltaHbO');
        xlabel('Time [s]');
        ylabel('Concentration change');
        title(sprintf('Channel %d', chIdx));
        legend('show');
        grid on;
    end
    sgtitle(dataFile, 'Interpreter', 'none');
end

end
