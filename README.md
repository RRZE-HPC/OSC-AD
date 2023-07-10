# OSC-AD
"Physical Oscillator Model for Supercomputing"

A parallel program together with the parallel hardware it is running on is 
not only a vehicle to solve numerical problems, it is also a
complex system with interesting dynamical behavior: resynchronization and
desynchronization of parallel processes, propagating phases of idleness, and
the peculiar effects of noise and system topology are just a few
examples. We propose a physical oscillator model (POM) to describe 
aspects of the dynamics of interacting parallel processes. Motivated
by the well-known Kuramoto Model, a process with its
regular compute-communicate cycles is modeled as an oscillator which is
coupled to other oscillators (processes) via an interaction potential. 
Instead of a simple all-to-all connectivity as in the standard Kuramoto
model, we employ a sparse topology matrix mapping the communication structure 
and thus the inter-process dependencies of the parallel program onto the 
oscillator setup.  The particular form of the interaction 
potential is decisive for the dynamics. We propose two potentials
that are suitable for different scenarios in parallel computing: 
one for resource-scalable and one for resource-bottlenecked 
applications. The former are not limited by a resource bottleneck such
as memory bandwidth or network contention, while the latter are. 
As opposed to the original Kuramoto model, which has a periodic
sinusoidal potential that is attractive for small angles, a parallel
program does not allow for arbitrary 
phase shifts between pairs of communicating processes; hence, 
the characteristic potentials are always attractive for large 
angles and only differ in the short-distance behavior. 
We show that the POM with appropriate potentials can mimic
the propagation of delays and the synchronizing and desynchronizing 
behavior of scalable and bottlenecked parallel programs, respectively.
