<h1 align="center">
  <img alt="Logo" src="https://raw.githubusercontent.com/adosar/moxel/master/docs/source/images/moxel_logo.svg"/>
</h1>

<h4 align="center">
  
[![Requires Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python&logoColor=yellow&label=Python&labelColor=black&color=blue)](https://www.python.org/downloads/)
[![Licensed under GPL-3.0-only](https://img.shields.io/badge/GPL--3.0--only-gold?label=License&labelColor=black)](https://spdx.org/licenses/GPL-3.0-only.html)
[![Read the Docs](https://img.shields.io/badge/stable-green?logo=readthedocs&logoColor=blue&label=Read%20the%20Docs&labelColor=black)](https://moxel.readthedocs.io)
[![pip install pymoxel](https://img.shields.io/badge/install-blue?logo=pypi&logoColor=yellow&label=PyPI&labelColor=black)](https://pypi.org/project/pymoxel/)
[![Documentation Status](https://readthedocs.org/projects/moxel/badge/?version=stable)](https://moxel.readthedocs.io/en/stable/?badge=stable)
[![PyPI version](https://badge.fury.io/py/pymoxel.svg)](https://badge.fury.io/py/pymoxel)

</h4>

MOXελ is a Python package for **parallel calculation of energy voxels**, with
emphasis on reticular chemistry.

The majority of time in a ML workflow goes into constructing the inputs and
making sure they are clean, rather than focusing on the ML part itself.

MOXελ aims to provide a **simple and fast interface to generate energy voxels in
a ML-ready format**, minimizing as much as possible the time spent on these
preprocessing steps.

<p align="center">
  <img alt="Voxels" src="https://raw.githubusercontent.com/adosar/moxel/master/docs/source/images/voxels.gif" width="25%"/>
</p>

## ⚙️  Installation
It is strongly recommended to **perform the installation inside a virtual environment**.

Check the [installation steps](https://moxel.readthedocs.io/en/stable/installation.html).

Assuming an activated virtual environment:
```sh
pip install pymoxel
```

## 🚀 Usage

```
moxel path/to/CIFs path/to/voxels_data/ --grid_size=5
```

You can also use a configuration file:

```
moxel --config=path/to/config.yaml
```

> [!NOTE]
> For more information, please refer to the [📚 Documentation](https://moxel.readthedocs.io/en/stable/index.html).

## 📰 Citing MOXελ
If you use ΜΟΧελ in your research, please consider citing the following work:

    @article{Sarikas2024,
    title = {Gas adsorption meets deep learning: voxelizing the potential energy surface of metal-organic frameworks},
    volume = {14},
    ISSN = {2045-2322},
    url = {http://dx.doi.org/10.1038/s41598-023-50309-8},
    DOI = {10.1038/s41598-023-50309-8},
    number = {1},
    journal = {Scientific Reports},
    publisher = {Springer Science and Business Media LLC},
    author = {Sarikas,  Antonios P. and Gkagkas,  Konstantinos and Froudakis,  George E.},
    year = {2024},
    month = jan
    }

## 📇 TODO
1. Improve performance
2. Improve voxelization scheme
3. Improve modeling of interactions

## 📑 License
MOXελ is released under the [GNU General Public License v3.0 only](https://spdx.org/licenses/GPL-3.0-only.html).
