# kinkLUT
A lookup table (LUT) for the analysis of coronal loop kink oscillations, as described in:

_Coronal loop seismology using standing kink oscillations with a lookup table_<br/>
David James Pascoe, Alan William Hood, Tom Van Doorsselaere<br/>
Frontiers in Astronomy and Space Sciences, 2019, 6, 22

Version 1.0 of the LUT is stored in the IDL save file `kinklut_v1_0.sav` and its use is demonstrated in the routine `kinklut_example`.

The strong damping of kink oscillations is attributed to resonant absorption, which occurs in loops which have a transition layer between the high density core and lower density background.

Analytical descriptions for damping by resonant absorption demonstrate that the damping profile depends on the density contrast ratio and layer width, but presently only describe the behaviour for thin inhomogeneous layers.

This project aims to provide a simple method of estimating the damping profile of kink oscillations without using the thin boundary (TB) approximation and without prior assumption of the form of the damping profile.

The results of hundreds of numerical simulations for various combinations of density contrast ratios and layer widths have been summarised in the LUT.
