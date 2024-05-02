# Documentation

Some notes about the Sphinx and Read-The-Docs (RTD) documentation builders:

- The documentation has two sections, the user documentation section with files located in the 'wiki' directory and developer documentation in the `dev_docs` directory. When a PDF is generated, only the wiki version is included. 
- Requirements for the docs can be installed via `pip install -r doc/rtd_requirements.txt`. 
    - For MacOS pdf builds you may need to install `pango` as follows `brew install pango`. See `sphinx_simplepdf` docs at the link below. 
- The build can be configured in `doc/conf.py`.
     - To build pdf locally run `sphinx-build -M html doc` or in the `doc` folder `make html`. The html will be built in `_build/html/index.html`
- To use markdown files, i.e `.md` files, the `myst_parser` extension is used. [Documentation for  `myst_parser`](https://myst-parser.readthedocs.io/en/latest/).
- To convert the docs to PDF, `sphinx_simplepdf` extension is used.
    - The RTD build is configured in `.readthedocs.yaml` under the `commands` list. 
    - If `sphinx_simplepdf` ever breaks, you can remove the custom build instructions there. 
    - To build pdf locally run `sphinx-build -M simplepdf doc`.
    - [Documentation for  `sphinx_simplepdf`](https://sphinx-simplepdf.readthedocs.io/en/latest/index.html).
- As of this version, the color scheme for the docs are the following (can be changed in `doc/conf.py`):
    - `primary_color = '#3061f3'`
    - `secondary_color = '#1A43A4'`
    - `text_color = '#ffffff'`
