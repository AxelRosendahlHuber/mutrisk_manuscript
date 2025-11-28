# Genetic diversity in human tissues

Axel Rosendahl Huber, Ferran Muiños, Joan Enric Ramis-Zaldivar, Maria Andrianova, Abel Gonzalez Perez, Núria Lopez-Bigas

## Installation

Clone the repository to a folder in your PC

Open `mutrisk.manuscript.Rrpoj` to start

``` r
# if not installed, install devtools
if (!"devtools" %in% rownames(installed.packages())) {
  install.packages("devtools")
}

devtools::install_github("AxelRosendahlHuber/wintr")
devtools::install_github("AxelRosendahlHuber/mutrisk")

# install all other packages required: 
install.packages("renv")
renv::restore()
```
