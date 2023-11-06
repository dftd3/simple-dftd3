#!/bin/bash


meson compile -C _build

_build/app/s-dftd3 $1 $2 $3 $4 $5 $6 $7 $8 $9
