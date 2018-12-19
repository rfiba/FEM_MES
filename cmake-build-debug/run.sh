#!/bin/bash
echo "Settings: " $1
echo "Results: " $2
./FEM <$1 >$2
