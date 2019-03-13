# OmaDB: R wrapper for the OMA REST API

R package providing access to the data in the OMA browser, using the REST API. As such, it requires a stable internet connection to operate. We also provide a similar _wrapper for python (https://github.com/DessimozLab/pyOmaDB)_.

An notebook containing examples of how to use the package is available <a href="https://github.com/DessimozLab/omadb/blob/master/manuscript_examples.ipynb">here</a>.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/DessimozLab/omadb/master?filepath=manuscript_examples.ipynb)

## Citation
If you use our package in your work, please consider citing:

_Kaleb K, Warwick Vesztrocy A, Altenhoff A and Dessimoz C. Expanding the Orthologous Matrix (OMA) programmatic interfaces: REST API and the OmaDB packages for R and Python. F1000Research 2019, 8:42
(https://doi.org/10.12688/f1000research.17548.1)_

## Installation

The package is available via Bioconductor (https://bioconductor.org) and can be installed as follows:

```
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install('OmaDB')

#load the package
library(OmaDB)

```

Note that the v2.0 of the package requires R version >= 3.6 and Bioconductor version >=3.9. At the time of writting this, this is the Bioconductor development version and so the package can be installed as follows:

```
BiocManager::install('OmaDB', version = 'devel')
```
Alternatively, the latest version of the OmaDB package can also be directly installed from the github repository as below:

```
install.packages('devtools')
library(devtools)
install_github('dessimozlab/omadb')
```


## License

OmaDB is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OmaDB is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License
along with OmaDB.  If not, see <http://www.gnu.org/licenses/>.

