# Monte Carlo simulations of diffusion through 1d permeable membranes (CUDA C++)

The code re-implements 1d Monte Carlo simulations originally developed in [Fieremans, et al., NMR Biomed, 2010](https://doi.org/10.1002/nbm.1577) and [Novikov, et al., Nature Physics, 2011](https://doi.org/10.1038/nphys1936), demonstrating the power spectrum of bead distribution along axons (Figure 8) and diffusivity and kurtosis time-dependence (Figure 9) in [Lee and Papaioannou, et al., NeuroImage 2019](https://doi.org/).

* **Demo 1, bead statistics:** Concatenate the bead distributions along 33 axons coming from [Hellwig, et al., Biological Cybernetics, 1994](https://doi.org/10.1007/BF00198906), and calculate its power spectrum (Figure 8a) and the number of beads within a sliding window (Figure 8c).
* **Demo 2, packing generation:** Generation of randomly distributed membranes based on the bead distance along axons in [Hellwig, et al., Biological Cybernetics, 1994](https://doi.org/10.1007/BF00198906).
* **Demo 3, simulations:** Perform Monte Carlos simulations of diffusion through 1d permeable membranes. The code is implemented in CUDA C++, and you need an Nvidia GPU to run the code.
* **Demo 4, analysis:** Calculate diffusivity and kurtosis time-dependence based on displacement cumulants and diffusion signals.

These are good exercises if you just start your own MC simulation codes.
Some results can suprise you, even if you are well experienced!!

## References
* **Monte Carlo simulation**
  - [Fieremans, et al., NMR Biomed, 2010](https://doi.org/10.1002/nbm.1577)
  - [Novikov, et al., Nature Physics, 2011](https://doi.org/10.1038/nphys1936)
  - [Fieremans and Lee, NeuroImage 2018](https://doi.org/10.1016/j.neuroimage.2018.06.046)
  - [Lee and Papaioannou, et al., NeuroImage 2019](https://doi.org/)

* **Bead distribution along axons**
  - [Hellwig, et al., Biological Cybernetics, 1994](https://doi.org/10.1007/BF00198906)

## Authors
* [Hong-Hsi Lee](http://www.diffusion-mri.com/people/hong-hsi-lee)
* [Antonios Papaioannou](http://www.diffusion-mri.com/people/antonios-papaioannou)
* [Dmitry S Novikov](http://www.diffusion-mri.com/people/dmitry-novikov)
* [Els Fieremans](http://www.diffusion-mri.com/people/els-fieremans)

## License
This project is licensed under the [LICENSE](https://github.com/NYU-DiffusionMRI/monte-carlo-simulation-recipes/blob/example1/LICENSE).
 
