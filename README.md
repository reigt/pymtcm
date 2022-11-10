# Modified Tension Chord Model

## What is Modified Tension Chord Model (MTCM)
The MTCM provides the tensile response of a reinforced concrete (RC) tie and is a modified version of the Tension Chord Model (TCM). The main difference is that the second order differential equation for the slip is solved using the nonlinear bond slip law in *fib* Model Code 2010 instead of the stepped rigid perfect plastic bond slip law in the case of the TCM. The solution of the differential equation for the slip provides thus two solution regimes 

- Comparatively Lightly Loaded Members (CLLM): Strains become compatible at an abcissa *x* from the loaded end 
- Comparatively Heavily Loaded Members (CHLM): Incompatible strains all over the member length

This package provides the solutions to the MTCM in an object oriented programming format in Python. See application examples below. 

## Installation

MTCM is installed directly from source by typing:

```pwsh
py -m pip install git+https://github.com/reigt/mtcm@main
```

## Application example

### Import
Import the package.

```python
import mtcm
```

### MTCM class
An object of the MTCM class is instantiated as 

```python
tie = mtcm.mtcm(As_tot,n_phi_s,phi_s,rho_s,Es,Ecm,fctm)
```

Note that in order to instantiate the object the necessary parameters such as 

- `As_tot`: Total reinforcement area in RC tie
- `n_phi_s`: Number of rebars with diameter `phi_s`
- `phi_s`: 
