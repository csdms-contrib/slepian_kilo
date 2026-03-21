# Changelog

All changes to Slepian_kilo repo will be documented in this file

## [0.0.1] - 03-21-2026
There was an issue in the plotit function of simulros.m file. The lines I have commented out (lines 261-265) made limc = [0, 0] because the values of limc already to small (on the order of 10^(-4). Because limc was the zero array, the colorbar wasn't working because the min and max were both 0. 
- **[MINOR}** 'src/simulros.m': commented out limc rounding lines


## [0.0.0] - 03-20-2026
Had to change ''FontS'' to "FontSize"

### Changed
- **[MINOR]** 'src/admittance.m': FontS was changed to FontSize.
- **[MINOR]** 'src/xaxis1d.m": FontS was changed to FontSize.

