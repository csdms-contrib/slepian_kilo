# Changelog

All changes to Slepian_kilo repo will be documented in this file

## [0.0.3] - 04-04-2026
I started debugging mleros.m (simulros.m seems to be working 100%). I used th0=[7e22 0.4 -0.75 0.0025 2 2e4], N = 4. I had to replace ProbDistUnivParam on line 619 to makedist('Gamma', 'a', df/2, 'b', 2) (MATLAB uses upper and lower bounds as opposed to a vector on newer versions). 
- **[MINOR]** 'src/mleros.m': Line 216 was changed to the makedist function to be compatible with newer MATLAB versions
- ** [MINOR]** 'src/mleros.m': Line 603 was changed from hist() to histogram() to support newer MATLAB version
- ** [MINOR]** 'src/mleros.m': Lines 372 & 373, 392 & 393, 410 & 411, I deleted because np was already defined in line 268.CORRECTION: I think this is actually not a good thing to do because of parameter scope
- ** [MINOR]** 'src/mleros.m': Removed a keyboard from what is now line 750.
- ** [MINOR]** 'src/mleros.m': Removed a repetitive definition of nperturb in what is now line 660. 
- ** [MINOR]** 'src/mleros.m': Removed a keyboard from what is now like 490.
 
>> mleros('demo1',N,th0)
## [0.0.2] - 03-27-2026
To run demo7 of mleros.m, I need hes2cov function, couldn't find it on your website. I need demo1 to try the other demos. I put MaxIter at 1 to check convergence. 
- **[MINOR]**  'src/mleros.m': Line 271 had an error because there were "too many output arguments" in osopen file function. I truncated the last 4 arguments into fmti (this was to run demo1, which still hasn't converged).
- **[MINOR]** 'src/mleros.m': Throuhgout the script I changed mleros0 and simulros0 and removed the 0. This was for demo 6 but it still hasn't converged.
- 

## [0.0.1] - 03-21-2026
There was an issue in the plotit function of simulros.m file. The lines I have commented out (lines 261-265) made limc = [0, 0] because the values of limc were already very small (on the order of 10^(-4). Because limc was the zero array, the colorbar wasn't working because the min and max were both 0. 
- **[MINOR]** 'src/simulros.m': commented out limc rounding lines


## [0.0.0] - 03-20-2026
Had to change ''FontS'' to "FontSize"

### Changed
- **[MINOR]** 'src/admittance.m': FontS was changed to FontSize.
- **[MINOR]** 'src/xaxis1d.m": FontS was changed to FontSize.

