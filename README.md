# Genetic diversity in human tissues

Axel Rosendahl Huber, Ferran Muiños, Joan Enric Ramis-Zaldivar, Maria Andrianova, Abel Gonzalez Perez, Núria Lopez-Bigas

## Installation

Clone the repository to a folder in your PC

Open `mutrisk.manuscript.Rrpoj` to start

``` r
# if not installed, install devtools
if (!"pak" %in% rownames(installed.packages())) {
  install.packages("pak")
}

pak::pak("AxelRosendahlHuber/wintr")
pak::pak("AxelRosendahlHuber/mutrisk")
pak::pak("gersteinlab/siglasso")

# install all other packages required: 
install.packages("renv")
renv::restore()
```
