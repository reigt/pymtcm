# Modified Tension Chord Model

## What is the Modified Tension Chord Model (MTCM)
The MTCM provides the tensile response of a reinforced concrete (RC) tie and is a modified version of the original Tension Chord Model (TCM). The main difference is that the second order differential equation for the slip is solved using the nonlinear bond slip law in *fib* Model Code 2010 instead of the stepped rigid perfect plastic bond slip law in the case of the TCM. The solution of the differential equation for the slip provides thus two solution regimes 

- Comparatively Lightly Loaded Members (CLLM): Strains become compatible at an abcissa *x* from the loaded end 
- Comparatively Heavily Loaded Members (CHLM): Incompatible strains all over the member length

This package provides the solutions to the MTCM in an object oriented programming language in Python. See application examples below. 

## Installation

MTCM is installed directly from source by typing:

```pwsh
py -m pip install git+https://github.com/reigt/pymtcm@main
```

## Application example

### Import
Import the package.

```python
import pymtcm
```

### MTCM class
An object of the MTCM class is instantiated by providing the necessary input for the geometry of a RC tie as 

```python
tie = pymtcm.mtcm(phi_s,n_phi_s,hc_ef,wc_ef,Es,Ecm,fctm)
```

Note that the following parameters are necessary arguments in order to instantiate the object 

- `phi_s`: Rebar diameter [mm]
- `n_phi_s`: Number of rebars with diameter `phi_s`
- `hc_ef`: Effective height [mm]
- `wc_ef`: Effective width [mm]
- `Es`: Young's modulus steel [MPa]
- `Ecm`: Young's modulus concrete [MPa]
- `fctm`: Mean tensile strength concrete [MPa]

The response of the tie caused by a reinforcement stress of e.g. `200 MPa` at the crack is found by 

```python
tie.stress(200)
```

The response of the tie caused by an imposed mean strain of e.g. `0.3 ‰` 

```python
tie.strain(0.3e-3)
```

### PlotMTCM class
An object of the PlotMTCM class is instantiated as 

```python
plot = pymtcm.PlotMTCM(tie)
```

The distribution of slip, bond stresses and strains can be visualized by applying the following method

```python
plot.chord_distribution()
```

The tie response, i.e. stress vs. strain curves, crack widths and transfer lengths, is found by applying the following method

```python
plot.tie_response()
```
### MISO class
MISO is a child class of the MTCM class. It discretizes the MTCM to a multilinear isotropic (MISO) hardening curve and generates an input file that can be used for ANSYS analyses. An object of the MISO class is thus instantiated in a similar manner as the MTCM class

```python
miso = pymtcm.miso(phi_s,n_phi_s,hc_ef,wc_ef,Es,Ecm,fctm)
```

The ANSYS input file is generated by applying the following method

```python
miso.ansys_input()
```

The discretized hardening curve can be visulized by instantiating and applying the following method in the PlotMTCM class

```python
plot = pymtcm.PlotMTCM(miso)
plot.plot_miso()
```
