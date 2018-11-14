#!/bin/bash

make
test_dirichlet
python3 test_data/plot.py
test_periodic
python3 test_data/plot.py
