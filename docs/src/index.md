# Introduction

`SelfAvoidingBilliards` is a Julia package to simulating the a mathematical billiard with an added memory.
This package is largely copied from `DynamicalBilliards.jl`, please also refer to that [documentation](https://juliadynamics.github.io/DynamicalBilliards.jl/dev/) for more information.

This documentation will cover the important differences from the original code, as well as some practical examples on how to use the code.

## Features

* Automatically adding walls when bouncing a particle
* Reducing the billiard to remove unnessecary walls

## Installation

The package is not (yet) featured on the registry, so for now it can be installed by typing `] add https://github.com/stijndelnoij/SelfAvoidingBilliards` in shell mode, or `Pkg.add("https://github.com/stijndelnoij/SelfAvoidingBilliards")` using `Pkg`.
