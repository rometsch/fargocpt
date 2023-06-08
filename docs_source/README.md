# Documentation

## Build the html documentation

There are two options to install the requirements.

1) run `scripts/setup_sphinx_venv.sh` to create a new venv and install sphinx in it.
2) install the sphinx requirements with `pip install -r requirements.txt`

Then run
```
make html
```
to compile the html version of the documentation.

## View the html documentation

Either open the file `build/html/index.html` with your webbrowser or run
```
./scripts/serve.sh
```
to start a local webserver. This has the added benefit of autoreloading when you save a change to the source files.
