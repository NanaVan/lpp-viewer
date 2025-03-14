# lpp-viewer
An online tool to display spectrum simulation according to .lpp file from LISE++
* Constant values: CODATA 2018
* Nuclear data: NUBASE2020
* Python library: `SQLite3`, `Numpy 1.22`, `Scipy`, `Bokeh 3.2.2, Panel 1.3.8 (online-pyscript)`, `Bokeh 3.6.2, Panel 1.5.4`

## Web version: deploy from .\pyscript-lpp-viewer
webside address: [China/Overseas](https://lpp-viewer.pages.dev/)
* additional Python library: `PyScript`
* note: this version is extremely slow for first loading and data updating because of WASM service

## Stand alone version: deploy from .\stand-alone
* local requirement: Bokeh(Panel) server 
1. open `powershell`
* Usage: 
```shell
> bokeh serve .\stand-alone
```
2. open browser and visit [http://localhost:5006/stand-alone](https://localhost:5006/stand-alone)

## Additional tool: ion data generator from .\database-maker
generate ion data `ionic_data.db` from NUBASE2020 and NIST binding energy file
* Usage:
```shell
> cd .\database-maker
> python gen_nubase.py
``` 
