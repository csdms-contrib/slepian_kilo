# Changelog

All changes to Slepian_kilo repo will be documented in this file

## [0.0.2] - 03-27-2026
To run demo7 of mleros.m, I need hes2cov function, couldn't find it on your website. I need demo1 to try the other demos. I put MaxIter at 1 to check convergence. 
- **[MINOR]**  'src/mleros.m": Line 271 had an error because there were "too many output arguments" in osopen file function. I truncated the last 4 arguments into fmti (this was to run demo1, which still hasn't converged).
- **[MINOR]** 'src/mleros.m": Throuhgout the script I changed mleros0 and simulros0 and removed the 0. This was for demo 6 but it still hasn't converged.
- 

## [0.0.1] - 03-21-2026
There was an issue in the plotit function of simulros.m file. The lines I have commented out (lines 261-265) made limc = [0, 0] because the values of limc were already very small (on the order of 10^(-4). Because limc was the zero array, the colorbar wasn't working because the min and max were both 0. 
- **[MINOR]** 'src/simulros.m': commented out limc rounding lines


## [0.0.0] - 03-20-2026
Had to change ''FontS'' to "FontSize"

### Changed
- **[MINOR]** 'src/admittance.m': FontS was changed to FontSize.
- **[MINOR]** 'src/xaxis1d.m": FontS was changed to FontSize.

