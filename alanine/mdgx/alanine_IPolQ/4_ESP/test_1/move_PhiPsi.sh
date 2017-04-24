#!/bin/bash

cp -r ../../3_QM/output/PhiPsi_1 .

rm PhiPsi_1/normal_termination.py
rm PhiPsi_1/problem.dat

rm -rf PhiPsi_1/120/30
rm -rf PhiPsi_1/30/120

#now remove te files we don't need
rm PhiPsi_1/*/*/Conf*
rm PhiPsi_1/*/*/fort.7
rm PhiPsi_1/*/*/For*
rm PhiPsi_1/*/*/IPolQinp.solv
rm PhiPsi_1/*/*/IPolQinp.vacu
rm PhiPsi_1/*/*/ipolq.out
rm PhiPsi_1/*/*/IPolQout.*
rm PhiPsi_1/*/*/mdgxQM.in
rm PhiPsi_1/*/*/output-*
