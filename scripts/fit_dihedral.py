#!/usr/bin/python
# $Id: fit_dihedral.py,v 1.6 2008/08/26 18:37:58 oguvench Exp $
# Monte Calo fit of dihedral parameters
# includes variation in dihedral phase angle
# Olgun Guvench, 2006-2008
# This software is released under the GNU General Public License
# Please cite:
#   Guvench, O.; MacKerell, A. D., Jr. Automated conformational energy
#   fitting for force-field development. J Mol Model (2008) 14:667-679

import math
#import os
import random
import string
import sys

#datapath    = sys.argv[1]
#paramscheme = int( sys.argv[2] )
#phases      = int( sys.argv[3] )
phases      = 0
#doweight    = int( sys.argv[4] )
#kmax        = float( sys.argv[5] )
run         = int( sys.argv[1] )
run         = "run%i" % ( run )

# version log
# 0.07  activated phase angle sampling support
# 0.06  added constant-temperature support
print "Monte Carlo dihedral fitting v 0.07"
print "Olgun Guvench"
print "MacKerell Lab"
print "University of Maryland, Baltimore"
print "2006-2008"
print "This software is released under the GNU General Public License"
print "Please cite:"
print "   Guvench, O.; MacKerell, A. D., Jr. Automated conformational energy"
print "   fitting for force-field development. J Mol Model in press."
print ""
print "Input the maximum value for k (recommended: 3.0)"
kmax = float( string.strip( sys.stdin.readline() ) )
print "  kmax = %5.2f" % ( kmax)
i = 1
while i:
  print "Would you like to do exponential cooling? [y/n]"
  answer = string.strip( sys.stdin.readline() )
  print "  %s" % ( answer )
  if answer == "y" or answer == "n":
    if answer == "y":
     use_exp_cooling = 1
     print "What is the desired initial temperature?"
    elif answer == "n":
     use_exp_cooling = 0
     print "What is the desired temperature?"
    i = 0
answer = string.strip( sys.stdin.readline() )
print "  %s" % ( answer )
tempr0 = float( answer )
print "Input the maximum number of Monte Carlo steps (recommended: 5000)"
nstep = int( string.strip( sys.stdin.readline() ) )
print "  nstep = %i" % ( nstep)
i = 1
while i:
  print "Would you like to allow the phases to vary? [y/n] (recommended: n)"
  answer = string.strip( sys.stdin.readline() )
  print "  %s" % ( answer )
  if answer == "y" or answer == "n":
    if answer == "y":
     phases = 1
    elif answer == "n":
     phases = 0
    i = 0
if phases:
  kmin   =    0.0
  phimin = -180.0
  phimax =  180.0
else:
  kmin   =  -kmax
  phimin =    0.0
  phimax =    0.0

# read in qm energies
print "What is the file path for the qm energy file?"
path = string.strip( sys.stdin.readline() )
print "What is the file name of the qm energy file?"
file = string.strip( sys.stdin.readline() )
print "  qm energy file = %s/%s" % ( path, file )
qme = []
nconfqm = 0
in_file = open( "%s/%s" % ( path, file ), "r" )
for line in in_file.readlines():
  field = string.split( string.strip( line ) )
  #here he appends the QM energies
  qme.append( float( field[0] ) )
  nconfqm += 1
print "Read %i qm energies." % ( nconfqm )

# read in mm energies
print "What is the file path for the mm energy file?"
path = string.strip( sys.stdin.readline() )
print "What is the file name of the mm energy file?"
file = string.strip( sys.stdin.readline() )
print "  mm energy file = %s/%s" % ( path, file )
nconf = 0
mme = []
in_file = open( "%s/%s" % ( path, file ), "r" )
for line in in_file.readlines():
  field = string.split( string.strip( line ) )
  #here the mm
  mme.append( float( field[0] ) )
  nconf += 1
if nconf != nconfqm:
  raise "Number of mm energies does not equal number of qm energies"

print "How many dihedral angles are to be fit?"
nfile = int( string.strip( sys.stdin.readline() ) )
print " number of angles to fit = %i" % ( nfile )


ifile = 0
dihedral = []
d = {}
dihedral_charmm = []
while ifile < nfile:
  print "What is the file path for dihedral file %i?" % ( ifile + 1 )
  path = string.strip( sys.stdin.readline() )
  print "What is the file name of dihedral file %i?" % ( ifile + 1 )
  file = string.strip( sys.stdin.readline() )
  dihedral.append( file )
  print " file %i = %s/%s" % ( ifile + 1, path, file )
  in_file = open( "%s/%s" % ( path, file ), "r" )
  field = string.split( string.strip( in_file.readline() ) )
  if len( field ) != 4:
    raise "The first line in %s must contain four fields corresponding to atom types." % ( file )
  if field[0] > field[3]:
    field.reverse()
  elif field[0] == field[3]:
    if field[1] > field[2]:
      field.reverse()
  dihedral_charmm.append( field )
  d[ file ] = []
  iline = 0
  for line in in_file.readlines():
    field = string.split( string.strip( line ) )
    d[ file ].append( float( field[0] ) )
    iline += 1
  if iline != nconfqm:
    raise "The number of dihedral values in %s (%i)\ndoes not equal the number of qm energies (%)." % ( file, iline, nconfqm )
  ifile = ifile + 1

#og1 print dihedral
#og1 print dihedral_charmm

# equivalence dihedrals with same atom types
i = 0
iparamscheme = []
while i < len( dihedral ):
  iparamscheme.append( -1 )
  i += 1
i = 0
eq_tag = 0
#og1 print len( dihedral )
while i < ( len( dihedral ) ) :
  if iparamscheme[ i ] == -1:
    iparamscheme[ i ] = eq_tag
    j = i + 1
    while j < len( dihedral ) :
      if dihedral_charmm[j] == dihedral_charmm[i]:
        iparamscheme[ j ] = eq_tag
      j += 1
    eq_tag += 1
  i += 1
#og1  print i-1
#og1  print iparamscheme[i-1]

#og1 print iparamscheme
print "The following equivalences are in effect based on atom types:"
i = 0
while i < eq_tag:
  print "  equivalence-group #%i" % ( i + 1 )
  j = 0
  while j < len( dihedral ):
    if iparamscheme[ j ] == i:
      print "    %s" % dihedral[ j ]
    j +=1
  i +=1


i = 1
while i:
  print "Would you like to equivalence any other dihedral angles? [y/n]"
  answer = string.strip( sys.stdin.readline() )
  print "  %s" % ( answer )
  if answer == "y" or answer == "n":
    i = 0
while answer == "y":
  print "What are the names of the two files containing data for the dihedral angles you wish to equivalence? (separated by a space)"
  answer = string.split( string.strip( sys.stdin.readline() ) )
  print "  equivalencing %s and %s" % ( answer[0], answer[1] )
  ifile0 = -1
  ifile1 = -1
  j = 0
  for filename in dihedral:
    if answer[0] == filename:
      ifile0 = j
    if answer[1] == filename:
      ifile1 = j
    j += 1

# print ifile0
# print ifile1
  if ifile0 == -1:
    raise "%s is not a valid filename." % ( answer[0] )
  elif ifile1 == -1:
    raise "%s is not a valid filename." % ( answer[1] )
  else:
    eq_tag0 = iparamscheme[ifile0]
    eq_tag1 = iparamscheme[ifile1]
    if eq_tag0 < eq_tag1:
      eq_tagnew = eq_tag0
      eq_tagold = eq_tag1
    elif eq_tag0 > eq_tag1:
      eq_tagnew = eq_tag1
      eq_tagold = eq_tag0
    else:
      raise "%s is already equivalenced to %s." % ( answer[0], answer[1] )
    j = 0
    while j < len(iparamscheme):
      if iparamscheme[j] == eq_tagold:
        iparamscheme[j] = eq_tagnew
      j += 1
#   print iparamscheme
    print "The following equivalences are in effect:"
    j = 0
    eq_tagmax = 0
    while j < len( iparamscheme ):
      k = 0
      while k < len( dihedral ):
        if iparamscheme[ k ] == j:
          eq_tagmax += 1
          print "  equivalence-group #%i" % eq_tagmax
          k = len( dihedral ) + 1
        else:
          k += 1
      k = 0
      while k < len( iparamscheme ):
        if iparamscheme[ k ] == j:
          print "    %s" % dihedral[ k ]
        k +=1
      j +=1
  j = 1
  while j:
    print "Would you like to equivalence any other dihedral angles? [y/n]"
    answer = string.strip( sys.stdin.readline() )
    print "  %s" % ( answer )
    if answer == "y" or answer == "n":
      j = 0

#og1 print iparamscheme

# fill in the qme and mme arrays
qme_min = 9999999999999.9
mme_min = 9999999999999.9
i = 0
while i < nconf:
  if qme[ i ] < qme_min:
    qme_min = qme[ i ] + 0.0
  if mme[ i ] < mme_min:
    mme_min = mme[ i ] + 0.0
  i = i + 1
i = 0
while i < nconf:
  qme[ i ] = qme[ i ] - qme_min
  mme[ i ] = mme[ i ] - mme_min
  i = i + 1

###### initialize the parameters
k = []
k_old = []
k_best = []
use_multiplicity = []

phi = []
phi_old = []
phi_best = []

npar = 3
i = 0

#iparamscheme is the number of the terms to fit?
while i < len( iparamscheme ):
  j = 0
  k.append( [] )
  k_old.append( [] )
  k_best.append( [] )
  phi.append( [] )
  phi_old.append( [] )
  phi_best.append( [] )
  use_multiplicity.append( [] )
  while j < npar:
      #so for example we have 3 multiplicites we have to fit K1 K2 and K3
    k[i].append( random.uniform( kmin, kmax ) )
    k_old[i].append( random.uniform( kmin, kmax ) )
    k_best[i].append( 0.0 )
    use_multiplicity[i].append( 0.0 )
    if phases == 1:
      phi[i].append( random.uniform( phimin, phimax ) )
      phi_old[i].append( random.uniform( phimin, phimax ) )
      phi_best[i].append( 0.0 )
    else:
      phi[i].append( 0.0 )
      phi_old[i].append( 0.0 )
      phi_best[i].append( 0.0 )
    j = j + 1
  i = i + 1

k_to_dihe = {}
dihe_to_k = {}
i = 0
while i < len( dihedral ):
  k_to_dihe[ iparamscheme[i] ] = dihedral[i]
  dihe_to_k[ dihedral[i] ] = iparamscheme[i]
  i = i + 1

i = 0
eq_tagmax = 0
print "Equivalence-group multiplicities"
while i < len( iparamscheme ):
  j = 0
  while j < len( dihedral ):
    if iparamscheme[ j ] == i:
      eq_tagmax += 1
      print "  equivalence-group #%i" % eq_tagmax
      j = len( dihedral ) + 1
    else:
      j += 1
  if j > len( dihedral ):
    j = 0
    while j < len( iparamscheme ):
      if iparamscheme[ j ] == i:
        print "    %s" % dihedral[ j ]
      j +=1
    j = 0
    while j < npar:
      ask = 1
      while ask:
        print "    Use multiplicity %i? [y/n]" % (j+1)
        answer = string.strip( sys.stdin.readline() )
        print "      %s" % ( answer )
        if answer == "y":
          use_multiplicity[i][j] = 1.0
          ask = 0
        elif answer == "n":
          ask = 0
        else:
          ask = 1
      j += 1
  i += 1

###### get weights for conformations
i = 1
while i:
  print "Would you like to apply a weight factor to the energies? [y/n]"
  answer = string.strip( sys.stdin.readline() )
  print "  %s" % ( answer )
  if answer == "y":
    doweight = 1
    print "What is the file path for the weight file?"
    path = string.strip( sys.stdin.readline() )
    print "What is the file name of the weight file?"
    file = string.strip( sys.stdin.readline() )
    print "  weight file = %s/%s" % ( path, file )
    weight_file = open( "%s/%s" % ( path, file ) , "r" )
    i = 0
  elif answer == "n":
    doweight = 0
    i = 0
weight = []
tempr = 298.0
tempr = 2000.0
i = 0
weight_sum = 0.0
while i < nconf:
  delta_qm = qme[i] + 0.0
  if delta_qm > 8.0:
    delta_qm = 8.0
  weight.append( math.exp( -1.0 * delta_qm / ( 0.001987 * tempr ) ) )
  if doweight == 0:
    weight[i] = 1.0
  elif doweight == 1:
    weight[i] = float( string.split( string.strip( weight_file.readline() ) )[0] )
#g0
# print weight[i]
  weight_sum = weight_sum + weight[i]
#g1
  i = i + 1


###### calculations toward c
die_best = []
qme_average = 0.0
mme_average = 0.0
c_best = 0.0
i = 0
while i < nconf:
  qme_average = qme_average + weight[i]*qme[ i ] / weight_sum
  mme_average = mme_average + weight[i]*mme[ i ] / weight_sum
  die_best.append( 0.0 )
  i = i + 1
del_average = mme_average - qme_average
#print del_average


###### do the monte carlo simulated annealing
#tempr0 = 1000.0
print "%s" % ( "Starting Monte Carlo fitting" )
print "%11s%11s%11s%11s%11s%11s" % ( "MC step", "tempr", "p", "accepted?", "RMSE", "RMSE_best" )
step = 0
#if kmax == 0.0:
#  nstep = 10
#else:
#  nstep = 5000
#g0
#nstep = 10
#g1
rmsd_old = 999.0
rmsd_best = 999.0
while step < nstep:
  if use_exp_cooling:
    tempr = tempr0 * math.exp(-1.0* ( float( step ) / ( float( nstep ) / 4.0 ) ))
  else:
    tempr = float( tempr0 )
# print tempr
  i = 0
  while i < len( iparamscheme ):
    j = 0
    while j < npar:
      k[i][j] = k_old[i][j] + random.uniform(-0.5, 0.5)
      if k[i][j] > kmax:
        k[i][j] = kmax
      elif k[i][j] < kmin:
        k[i][j] = kmin
#g0
#     if i != paramscheme-1 or j < npar-1:
#     k[i][j] = 0.0
#g1
      if phases == 1:
        phi[i][j] = phi_old[i][j] + random.uniform(-18, 18)
        if phi[i][j] > phimax:
          phi[i][j] = phimax
        elif phi[i][j] < phimin:
          phi[i][j] = phimin
      j = j + 1
    i = i + 1

# dihedral energy to be fit
  die = []
  die_average = 0.0
  iconf = 0
  while iconf < nconf:
    die.append( 0.0 )
    i = 0
    while i < len( dihedral ):
      j = 0
      while j < npar:
#       print die
#       print die[0]
#       print iconf
#       print dihe_to_k[dihedral[i]]
#       print k[dihe_to_k[dihedral[i]]][j]
#       print d[dihedral[i]][iconf]
#       print phi[dihe_to_k[dihedral[i]]][j]
        die[ iconf ] = die[ iconf ] \
          + use_multiplicity[dihe_to_k[dihedral[i]]][j] \
          * k[dihe_to_k[dihedral[i]]][j] \
          * ( 1.0 + math.cos( ( float( j + 1 ) \
          * d[dihedral[i]][iconf] - phi[dihe_to_k[dihedral[i]]][j] ) \
          * math.pi / 180.0 ) )
        j = j + 1
      i = i + 1
    die_average = die_average + weight[iconf]*die[ iconf ] / weight_sum
#   print "%10i%20.2f" % ( iconf, die[iconf] )
    iconf = iconf + 1

#   calculate the value of c such that sum_i( ( eqm_i - emm_i + c )**2 )
#   is minimized:
#     d sum_i( w_i( eqm_i - emm_i + c )**2 ) / dc = 0
#       sum_i( 2(w_i*eqm_i - w_i*emm_i + w_i*c ) )         = 0
#       sum_i(  (w_i*eqm_i - w_i*emm_i + w_i*c ) )         = 0
#       sum_i(  (w_i*eqm_i - w_i*emm_i     ) )             = -sum_i( w_i*c )
#       sum_i(  (w_i*emm_i - w_i*eqm_i     ) ) / w_sum     = c
#                  <emm_i> -   <eqm_i>                     = c
#   where <emm_i> and <eqm_i> are weighted averages and w_sum is sum_i( w_i )

  c = del_average + die_average
#g0
# print c
# c = 0.0
# i = 0
# while i < nconf:
#   print mme[i]
#   print qme[i]
#   c = c + weight[i]*(mme[i]+die[i]) - weight[i]*qme[i]
#   i = i + 1
# c = c / float(nconf)
# print c
# print
#g1
  rmsd = 0.0
  i = 0
  while i < nconf:
    rmsd = rmsd  \
     + ( weight[i]*( qme[ i ] - ( mme[ i ] + die[ i ] )  + c )**2  )
    i = i + 1
  rmsd = math.sqrt( rmsd / weight_sum )

  drmsd = rmsd - rmsd_old
  if step == 0:
    p  = 0.0  # required to prevent overflow
  else:
    drmsd = rmsd - rmsd_old
    boltz = -1.0 * drmsd / ( 0.001987 * tempr )
    p  = math.exp( boltz )
  p0 = random.uniform(0.0,1.0)
# print "%10.5f%10.5f%10.5f" % (p, p0, drmsd)
  if drmsd < 0.0:
    accepted = 1
    p = 1.0
  elif p0 < p:
    accepted = 1
  else:
    accepted = 0
# print accepted
  if accepted:
    rmsd_old = rmsd
    i = 0
    while i < len( iparamscheme ):
      j = 0
      while j < npar:
        k_old[i][j]   = k[i][j] + 0.0
        phi_old[i][j] = phi[i][j] + 0.0
        j = j + 1
      i = i + 1
    if rmsd < rmsd_best:
      rmsd_best = rmsd
      c_best = c
      i = 0
      while i < len( iparamscheme ):
        j = 0
        while j < npar:
          k_best[i][j] = k[i][j] + 0.0
          phi_best[i][j] = phi[i][j] + 0.0
          j = j + 1
        i = i + 1
      iconf = 0
      while iconf < nconf:
        die_best[ iconf ] = die[ iconf ] + 0.0
        iconf = iconf + 1
  if ( step % 100 ) == 0:
#og  if ( step % 10 ) == 0:
    print "%11i%11.1f%11.4f%11i%11.2f%11.2f" % ( step, tempr, p, accepted, rmsd, rmsd_best )
  step = step + 1
print "%11i%11.1f%11.4f%11i%11.2f%11.2f" % ( step, tempr, p, accepted, rmsd, rmsd_best )

out_file = open( "fit_dihedral.%s.str" % (run), "w")
out_file.write( "* foo\n*\nread param card append\n* %10.5f%20.5f\n*\nDIHE\n" % ( rmsd_best, c_best ) )
i = 0
while i < len( dihedral ):
  skip = 0
  string = "%9s%9s%9s%9s" % ( dihedral_charmm[ i ][0],  dihedral_charmm[ i ][1],  dihedral_charmm[ i ][2],  dihedral_charmm[ i ][3] )
  if i > 0:
    j = 0
    while j < i:
      string2 = "%9s%9s%9s%9s" % ( dihedral_charmm[ j ][0],  dihedral_charmm[ j ][1],  dihedral_charmm[ j ][2],  dihedral_charmm[ j ][3] )
      if string == string2:
        skip = 1
      j += 1
  if not skip:
    j = 0
    while j < npar:
      if use_multiplicity[dihe_to_k[dihedral[i]]][j] > 0.1:
        par = k_best[dihe_to_k[dihedral[i]]][j]
        parphi = phi_best[dihe_to_k[dihedral[i]]][j]
        if par < 0.0:
          par = -par
          parphi = 180.0
        out_file.write( "%s %8.2f %2i %8.1f\n" % ( string, par, j+1, parphi ) )
      j += 1
  i += 1
out_file.write( "END\nRETURN" )
out_file.close()


out_file = open( "fit_dihedral.%s.ene" % (run), "w" )
mme_tot = []
i = 0
while i < nconf:
  mme_tot.append( mme[i] + die_best[ i ] )
  if kmax == 0.0:
    out_file.write( "%20.4f%20.4f%20.4f%20.4f\n" % ( qme[i]+qme_min, mme[i]+mme_min, mme_tot[ i ]+mme_min, qme[i] - mme_tot[ i ] + c_best ) )
  else:
    out_file.write( "%20.4f%20.4f%20.4f%20.4f\n" % ( qme[i], mme[i], mme_tot[ i ], qme[i] - mme_tot[ i ] + c_best ) )
  i = i + 1
out_file.close()
