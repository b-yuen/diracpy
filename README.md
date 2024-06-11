# diracpy

Diracpy is a python package for building quantum models quickly and effectively using Dirac notation. Diracpy provides a natural syntax for quantum models, enabling developers to focus on the physical model rather than the task of translating quantum physics to python code. Diracpy contains a diverse toolbox for solving the quantum dynamics of the model once it is built. It is therefore well suited to physicists with a good grasp of quantum theory and Dirac notation who want to translate this into a model that can be rapidly built and efficiently solved numerically.

Diracpy's functionality is split between four modules, `states_operators`, `subspaces`, `quantum_systems`, and `quantum_dynamics`. The `states_operators` module provides the basic definitions for `bra`, `ket` and `qop` (quantum operator) objects, together with all the neccesary binary composition rules and relations needed to define the Hilbert space of `ket` vectors, the dual space of `bra` vectors which maps `ket` vectors to scalars, and the `qop` operators that act on these. 

The `subspaces` module contains several pre-defined state spaces, consiting of a rule for how `ket` and `bra` vectors are indexed, and the rules for how the action of the relevent `qop` operatos modify these indexes and any complex coeffecient when mapping onto the output `bra` or `ket`.

The `quantum_systems` module builds a full quantum system from the user defined Hamiltonian `qop`, and a list of initial states. It uses the Hamiltonian to find states connected to the initial states via the interactions of the Hamiltonian. The contruction of the state space is automatically truncated after a user specified number of interactions. When the number of interactions specified is sufficiently large, the entire closed state space is constructed. The `quantum_system` objects allow states and operators to be automatically vectorized for subsequent numerical simulations of the quantum dynamics.

The `quantum_dynamics` module solves the quantum dynamics of the quantum system once it is built. A wide variety of methods are available, allowing one to choose the most efficient method for the particular problem. This includes most prevalent dynamical equations and methods typically found in quantum electrodynamics. This module includes dynamical methods based on numerical intergration, methods based on eigen-decomposition, and stochastic methods. It supports both open and closed quantum systems.

See also the full API documentation <https://diracpy.readthedocs.io/en/latest/index.html>.

## Installation

Source code and binary files can be downloaded from <https://github.com/b-yuen/diracpy>. To install from the source, navigate to the desired local directory that contains the project. Either download the project into this direcotry, or use

```bash
git clone https://github.com/b-yuen/diracpy.git
```

For more details see <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>.

To install the package binary:
1. Download `diracpy-1.0.0-py3-none-any.whl` from <https://github.com/b-yuen/diracpy>.
2. Create a virtual environment using `venv`. In a terminal window navigate to the directory where you want to keep the environment configuration (it does not matter where this). 
    ```bash
    python3 -m venv path/to/venv
    source path/to/venv/bin/activate
    python3 -m pip install path/to/downloaded/package/diracpy-1.0.0-py3-none-any.whl
    ```
    Here, `path/to/venv` can be any path, but is typically `.venv` which is hidden and often included in `.gitignore`.
3. Install the depedencies `numpy` and `scipy`:
    ```
    python3 -m pip install numpy scipy
    ```
    
## Usage

Simple operations with bra's and ket's:

```python
import diracpy as dp
import numpy as np
# define a ket
ket0 = dp.ket(['g'])

# define a superposition ket
ket1 = ( dp.ket(['e']) + dp.ket(['g']) )/np.sqrt(2)

# define a bra
bra0 = dp.bra(['g'])

# Compare vectors
bra0 == ket0.conj()

# returns comparison
True

# find overlap
bra0 * ket1

# returns overlap
0.7071067811865475

# Define opeator using outer product
g_projection = ket0 * ket0.conj()
print(g_projectoin * ket1)

# returns projection
0.7071067811865475 * ket['g']

```

Build tensor product spaces and Hamiltonian operators from pre-defined state spaces in the `subspaces` module.

```python
# Define state space
atom = dp.two_level_subspace(index=0)
cavity = dp.fock_subspace(index=1)

# define Hamiltonian
Delta, g = 0.1, np.pi
H0 = Delta * atom.sigma_z
V = g * (cavity.adag * atom.sigma_minus + cavity.a * atom.sigma_plus)
H = H0 + V

# Action of V on a the ket |g,1>.
psi0 = dp.ket(['g',1])
print(V * psi0)

#returns
3.141592653589793 * ket['e', 0]
```

Build the quantum system:

```python
# Build system from Hamiltonian and initial state
system = dp.qsys(H, initialstates=[psi0], n_int=2)

# List the basis states found by dp.qsys
system.print_basis()

# returns the basis for the closed subspace
1  *  ket ('g', 1)
1  *  ket ('e', 0)
```

Solve dynamics:

```Python
# Define inital density matrix and list of times to solve for
rho0 = np.zeros([2,2],complex)
rho0[0,0] = 1
times = np.linspace(0,1,5)

# Define unitary evolution solver object for example
usolver = dp.unitaryevolution(rho0, times, system)
usolver.solve

# Solution
usolver.soln

# returns density matrices for 5 specified times
array([[[ 1.00000000e+00+0.00000000e+00j,
          0.00000000e+00+0.00000000e+00j],
        [ 0.00000000e+00+0.00000000e+00j,
          0.00000000e+00+0.00000000e+00j]],

       [[ 5.00027179e-01+0.00000000e+00j,
         -7.95731459e-03+4.99936676e-01j],
        [-7.95731459e-03-4.99936676e-01j,
          4.99972821e-01+0.00000000e+00j]],

       [[ 2.53278377e-04-1.73472348e-18j,
         -1.59114633e-02-1.98905887e-04j],
        [-1.59114633e-02+1.98905887e-04j,
          9.99746722e-01-2.98008989e-19j]],

       [[ 5.00424940e-01-1.73472348e-18j,
         -7.95098402e-03-4.99936597e-01j],
        [-7.95098402e-03+4.99936597e-01j,
          4.99575060e-01-7.54442641e-19j]],

       [[ 9.99999842e-01-3.46944695e-18j,
         -2.51869390e-09+3.97811742e-04j],
        [-2.51869390e-09-3.97811742e-04j,
          1.58254207e-07+1.19961535e-19j]]])

```

For full API documentation, examples and tutorial see <https://github.com/b-yuen/diracpy>.

## Support

For support see first documentation at <https://diracpy.readthedocs.io/en/latest/index.html>. Issues can be raised here also. If you are interested in using this package and need additional support you can contact the developers via the research group <https://www.birmingham.ac.uk/research/activity/physics/quantum/metamaterials>.

## Roadmap

The next minor version will be released in June 2024 and includes full API documentation and complete integration of the quantum_dynamics module into the front-end. A future major release will address optimisation using C.

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

All changes can be tested using `tests/unit_tests.py` which should result in `All Tests Passed`.


## Authors and acknowledgment

The author of diracpy is Dr Ben Yuen <https://github.com/b-yuen>. We acknowledge contributions from <https://github.com/L-Hands> and <https://github.com/ishitajena> to create test script `tests/unit_tests.py` and <https://github.com/AngusCrookes> who developed the `quantum_dynamics.liouville` class.

## License

[MIT](https://choosealicense.com/licenses/mit/)
