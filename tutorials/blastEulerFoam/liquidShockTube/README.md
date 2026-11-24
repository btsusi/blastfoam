# Validation of polydispserse granular

## Notes

This case uses six granular phases within a continuous air phase. This a dilute simulation where packing of particles does no occur. The continuous gas phase is represented by air. Five of the particles phases are used to represent a continuous size distribution with a density 2700 kg/m^3 and a coefficient of restitution of 0.9. The last granular phase is that of coal with a density of 1300 kg/m^3, coefficient of restitution of 0.9, and a diameter of 10e-6 m.

## STUB: SchnerrSauer cavitation phase (not validated)

An eighth phase, "liquid", and a `(gas in liquid)` mass-transfer pair using
the `SchnerrSauer` cavitation model have been added to `constant/phaseProperties`
purely to exercise that part of the code (see `libblastInterfacialModels`'s
`massTransferModels/SchnerrSauer`). "gas" plays the role of the vapour
(dispersed) phase and "liquid" the continuous phase.

Everything associated with this addition is a placeholder:
- The liquid's Murnaghan equation-of-state coefficients (`rho0`, `pRef`, `K0`,
  `n`, `Gamma`) are rough, untuned, water-like order-of-magnitude guesses.
- The cavitation coefficients (`pSat`, `n`, `dNuc`, `Cc`, `Cv`) are untuned.
- `gas`'s diameter was bumped from `0` to `1e-4` m so it has a sensible value
  when acting as the dispersed phase of the new pair; this does not affect
  the original dust-in-gas physics, which uses the dispersed solid's own
  diameter.
- The `(gas in liquid)` drag/heat-transfer models reuse the same
  `GidaspowErgunWenYu`/`RanzMarshall` types as the existing dust pairs only
  because `phaseSystem` requires *some* entry for every pair; they have not
  been chosen or validated for a gas/liquid cavitation interface.
- `alpha.liquid` initializes to 0, so `liquid`'s volume fraction is derived
  as the complement of the other seven phases at start-up (no cavitation
  will occur until real liquid is actually present in the domain via
  `setFieldsDict` or new initial conditions).

None of these numbers should be trusted for a real cavitation simulation -
they exist only so the dictionary/solver wiring for `SchnerrSauer` can be
compiled and run. Tuning the EOS, initial fields (e.g. `setFieldsDict`) and
cavitation coefficients to a real test case (e.g. a liquid rarefaction/shock
tube) is still required.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

The calculation took approx. 3 min to run on a single core laptop at the reference resolution (e.g. 500 cells). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.

## Reference

```
Lai, S., Houim, R.W., Oran, E.S., 2018. Effects of particle size and density on dust dispersion behind a moving shock. Phys. Rev. Fluids 3, 064306. https://doi.org/10.1103/PhysRevFluids.3.064306


```
