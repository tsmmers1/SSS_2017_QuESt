# QuESt: Quantum Energy and Stuff 
A hybrid QM/MM project built on the principles of the MolSSI 2017 software
summer school.

## Installation
To install you first need the MolSSI 2017 Software Summer School stack.
Directions can be found
[here](https://molssi-sss.github.io/Logistics_SSS_2017/Setup.html)/

To build the C++ side of this project, please run:
```
python setup.py cmake
```

To clean the CMake files,
```
python setup.py clean
```

To do a local install of the Python directory,
```
pip install -e .
```


## Testing
Tests can be run using `py.test-v`.

