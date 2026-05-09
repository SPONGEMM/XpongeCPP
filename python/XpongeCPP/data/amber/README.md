# Amber Force-Field Data

This directory vendors the Amber force-field files shipped with the reference
Xponge repository under `Xponge/forcefield/amber`.

The files are copied into the Python package so the C++ core can load residue
templates and parameter data through stable package resources instead of relying
on the original Xponge source tree at runtime.
