# HW1 - Neurophotonics

MATLAB implementation of the Modified Beer-Lambert Law (MBLL) for computing
ΔHbR and ΔHbO from dual-wavelength fNIRS intensity data.

## Files

- `CalcNIRS.m` — main function. Loads a `.mat` recording, computes ΔHbR / ΔHbO
  for all 20 channels, and optionally plots selected channels.
- `run_HW1.m` — execution script. Runs `CalcNIRS` on both recordings, saves the
  channel plots, and performs the FFT / SNR analysis on channel 1 of file 1.

## Required data files (place in the same folder)

- `FN_031_V2_Postdose2_Nback.mat`
- `FN_032_V1_Postdose1_Nback.mat`
- `ExtinctionCoefficientsData.csv`
- `DPFperTissue.txt`
- `RelativeDPFCoefficients.csv`

## How to run

In MATLAB, from this folder:

```matlab
run_HW1
```

Or call the function directly:

```matlab
[dHbR, dHbO, fig] = CalcNIRS('FN_031_V2_Postdose2_Nback.mat', 3, 'adult_head', [1 2]);
```

### Function signature

```
[dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, ...
                             extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
```

- `SDS` — source-detector separation in cm (3)
- `tissueType` — key in `DPFperTissue.txt` (e.g. `adult_head`)
- `plotChannelIdx` — vector of channel indices to plot (use `[]` for no plot)

The script resolves all paths relative to its own location, so it works from
any directory as long as the data files are next to it.

## Outputs

- Channel plots: `FN_031_V2_Postdose2_Nback_HbR_HbO.png`,
  `FN_032_V1_Postdose1_Nback_HbR_HbO.png`
- FFT plot: `FFT_Channel1_File1.png`
- Console: detected heartbeat frequency (Hz / BPM) and SNR
