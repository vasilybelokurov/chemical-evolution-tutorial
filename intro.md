# Chemical evolution of dwarf galaxies

The abundances of a galaxy's stars record its history. Core-collapse supernovae enrich the gas quickly, mostly with
α-elements. Type Ia supernovae arrive later, mostly with iron. So [α/Fe] against [Fe/H] traces how fast a galaxy
formed its stars and how much gas it gained and lost. This tutorial builds that idea into a working model and
compares the model with five dwarf galaxies observed by APOGEE: the Large and Small Magellanic Clouds, Sagittarius,
Fornax and Gaia-Sausage/Enceladus.

| Chapter | What you do | Needs VICE |
|---|---|---|
| 0. Five dwarf galaxies in APOGEE | load the members of Hasselquist et al. (2021), check them, draw the [Mg/Fe] and [Si/Fe] tracks | no |
| 1. A one-zone model from scratch | write the model in numpy; check it against exact solutions | no |
| 2. Yields, delays and the [α/Fe] knee | add iron from type Ia supernovae; find what sets the plateau and the knee | no |
| 3. VICE | make VICE reproduce Chapter 2, then relax the approximations; bursts; [Fe/H] is not a clock | yes |
| 4. Fitting five dwarf galaxies | test the fit on mocks; fit 3403 precomputed models; compare with Hasselquist et al. | no |
| 5. What the track cannot see | mass, metallicity distribution and gas: what abundances miss | no |

**Prerequisites.** Python with numpy and matplotlib, and what [X/Fe] means. Each chapter runs on its own in Google
Colab (use the badge at its top) and downloads the small data files it needs.

**Reference.** Hasselquist, S., et al. 2021, ApJ, 923, 172, https://doi.org/10.3847/1538-4357/ac25f9.
