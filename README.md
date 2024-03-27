## Installation
### Install Julia 
To install the latest version of Julia, follow the `juliaup` installation guide here:
https://github.com/JuliaLang/juliaup/blob/61c09c460ca3a150e542522be8b8b8a5041e718a/README.md


### Clone this repository
git@github.com:jonathanfischer97/InteractiveOscillator.git


### Download solution data CSV
Use this Google Drive link, and move the CSV file into the cloned repository.


### Instantiate Julia environment
In order to reproduce the necessary package environment, first:
1. `cd` into the cloned repository.
2. Run `julia --project -e 'using Pkg; Pkg.instantiate()'`


## Running interactive visualization
To open the interactive visualization window, run the main script with 'julia --project visualizer.jl`.

