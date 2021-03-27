#!/usr/bin/env python

import numpy
import argparse

import itertools

import sys

import copy
import glob

# script worked ok
error=0

parser = argparse.ArgumentParser(description="Calculate the IPR of given VASP output")
parser.add_argument("-IBZ", dest="ibzkpt", default="", type=str)
parser.add_argument("-OUT", dest="outcar", default="", type=str)
parser.add_argument("-occup", dest="occupation", default="occupation.ground", type=str)
parser.add_argument("-spin1", dest="spin1", default="gamma.1.ground", type=str)
parser.add_argument("-spin2", dest="spin2", default="gamma.2.ground", type=str)
parser.add_argument("-CONT", dest="contcar", default="", type=str)
parser.add_argument("-ipr", dest="ipr", default="ipr.npy", type=str)
parser.add_argument("-iprrange", dest="ipr_range", default=15, type=int)
function_group = parser.add_argument_group("Task arguments", description="Arguments for which tasks to perform (None = all tasks)")
tasks = ["trans", "charge", "hyperfine", "D", "parchg", "spin"]
for task in tasks:
    function_group.add_argument("--" + task, default=False, action="store_true")
input = parser.parse_args()

# if no task is given, assume standard (do all)
do_all = not any(map(lambda x: getattr(input, x), tasks))
if do_all:
    for task in tasks:
        setattr(input, task, True)

## read input file defaults
if input.ibzkpt == "":
    input.ibzkpt = glob.glob("*IBZKPT*")[0]
if input.outcar == "":
    input.outcar = glob.glob("*OUTCAR*")[0]
if input.contcar == "":
    input.contcar = glob.glob("*CONTCAR*")[0]

# read files
ibzkpt=open(input.ibzkpt,'r')
outcar=open(input.outcar,'r')
occupation=open(input.occupation,'r')
spin_1=open(input.spin1,'r')
spin_2=open(input.spin2,'r')
contcar=open(input.contcar,'r')
ipr_range = input.ipr_range

# load IPR
#ipr = wav.inverse_participation_ratio()
ipr=numpy.load(input.ipr)

# find number of kpoints
ibzkpt.readline()
kpoints=int(ibzkpt.readline())

# find number of bands
occupation_lines = occupation.readlines()
spin_1_lines = spin_1.readlines()
spin_2_lines = spin_2.readlines()

bands = int(spin_2_lines[-1].split()[0])

# find number of electrons
outcar_lines = outcar.readlines()
electrons=0

for line in outcar_lines:
	if 'NELECT' in line:
		electrons=int(line.split()[2][:-5])
		break

# find highest occupied band (HOB)
HOB_1 = 0
HOB_2 = 0

for line in spin_1_lines:
	if line.split()[2]=='0.00000':
		HOB_1=int(line.split()[0])-1
		break

for line in spin_2_lines:
	if line.split()[2]=='0.00000':
		HOB_2=int(line.split()[0])-1
		break

# find which atoms are in cell
contcar_lines = contcar.readlines()
atoms = contcar_lines[5].split()

# close files
ibzkpt.close()
outcar.close()
occupation.close()
spin_1.close()
spin_2.close()
contcar.close()

# calculate ipr average over k-points
spin_1_average_ipr = numpy.zeros((bands))
spin_2_average_ipr = numpy.zeros((bands))

for band in range(bands):
	spin_1_average_ipr[band] = numpy.average(ipr[0,:,band,2])
	spin_2_average_ipr[band] = numpy.average(ipr[1,:,band,2])

# find local state close to HOB
spin_1_localstates = []
spin_2_localstates = []

IPR_limit=0.000015

for band in range(HOB_1-ipr_range,HOB_1+ipr_range):
	if spin_1_average_ipr[band]>IPR_limit:
		spin_1_localstates.append(band+1)
		if spin_1_average_ipr[band]<IPR_limit+0.000005:
			print 'WARNING, localstate could be wrong in spin1'

for band in range(HOB_2-ipr_range,HOB_2+ipr_range):
	if spin_2_average_ipr[band]>IPR_limit:
		spin_2_localstates.append(band+1)
		if spin_2_average_ipr[band]<IPR_limit+0.000005:
			print 'WARNING, localstate could be wrong in spin2'

def remove_discontinuous_elements(lst):
	to_remove=[]
	# keeps the elements in the middle that is continuous
	# which should correspond to local state in the band gap
	# if more advanced check is needed, compare with the other spin channel
	# there should be some states overlapping
	for i in range(len(lst)/2):
		if lst[i]+1 != lst[i+1]:
			to_remove.append(lst[i])
		if lst[-(1+i)]-1 != lst[-(2+i)]:
			to_remove.append(lst[-(1+i)])
	for value in to_remove:
		lst.remove(value)
	if to_remove == []:
		return 0
	else:
		remove_discontinuous_elements(lst)

if (len(spin_1_localstates) > 1):
	if spin_1_localstates[0]+len(spin_1_localstates)-1 != spin_1_localstates[-1]:
		print 'WARNING, local states are not continuous in spin 1, attempting to fix'
		old_spin_1_localstates=copy.deepcopy(spin_1_localstates)
		remove_discontinuous_elements(spin_1_localstates)
		if spin_1_localstates[0]+len(spin_1_localstates)-1 != spin_1_localstates[-1]:
			print 'ERROR, local states still not continuous in spin 1'
			error=1
		elif len(spin_1_localstates) == 0:
			print 'ERROR, removed too many items in list'
			error=1
		else:
			for value in spin_1_localstates:
				old_spin_1_localstates.remove(value)
			state=open('lone_states.txt','a+')
			state.write('Spin 1: '+str(old_spin_1_localstates)+'\n')
			state.close()
			print '   it worked, removed states are in lone_states.txt'

if (len(spin_2_localstates) > 1):
	if spin_2_localstates[0]+len(spin_2_localstates)-1 != spin_2_localstates[-1]:
		print 'WARNING, local states are not continuous in spin 2, attempting to fix'
		old_spin_2_localstates=copy.deepcopy(spin_2_localstates)
		remove_discontinuous_elements(spin_2_localstates)
		if len(spin_2_localstates) == 0:
			print 'WARNING, removed too many items in list, comparing with spin 1'
			for value in old_spin_2_localstates:
				if value in spin_1_localstates:
					spin_2_localstates.append(value)
		elif spin_2_localstates[0]+len(spin_2_localstates)-1 != spin_2_localstates[-1]:
			print 'ERROR, local states still not continuous in spin 2'
			error=1
		else:
			for value in spin_2_localstates:
				old_spin_2_localstates.remove(value)
			state=open('lone_states.txt','a+')
			state.write('Spin 2: '+str(old_spin_2_localstates)+'\n')
			state.close()
			print '   it worked, removed states are in lone_states.txt'

# find occupation in local states
spin_1_occupation = []
spin_2_occupation = []
spin_1_energy=numpy.zeros((bands+1))
spin_2_energy=numpy.zeros((bands+1))

for i in range(len(spin_1_localstates)):
	for index, line in enumerate(spin_1_lines,1):
		spin_1_energy[index]=float(line.split()[1])
		if int(line.split()[0])==spin_1_localstates[i]:
			spin_1_occupation.append(int(line.split()[2][0]))

for i in range(len(spin_2_localstates)):
	for index, line in enumerate(spin_2_lines,1):
		spin_2_energy[index]=float(line.split()[1])
		if int(line.split()[0])==spin_2_localstates[i]:
			spin_2_occupation.append(int(line.split()[2][0]))

# Excations
# find local excitations
spin_1_excitations=[]
spin_2_excitations=[]

spin_1_excitations = list(set(itertools.permutations(spin_1_occupation)))
spin_1_excitations.remove(tuple(spin_1_occupation))

for i in reversed(range(len(spin_1_excitations))):
	limit=numpy.dot(spin_1_occupation,spin_1_occupation)-1
	if (numpy.dot(spin_1_excitations[i],tuple(spin_1_occupation)) != limit):
		spin_1_excitations.remove(spin_1_excitations[i])

spin_2_excitations = list(set(itertools.permutations(spin_2_occupation)))
spin_2_excitations.remove(tuple(spin_2_occupation))

for i in reversed(range(len(spin_2_excitations))):
	limit=numpy.dot(spin_2_occupation,spin_2_occupation)-1
	if (numpy.dot(spin_2_excitations[i],tuple(spin_2_occupation)) != limit):
		spin_2_excitations.remove(spin_2_excitations[i])

# number of full and empty states
spin_1_full_local_state=numpy.count_nonzero(spin_1_occupation)
spin_2_full_local_state=numpy.count_nonzero(spin_2_occupation)
spin_1_empty_local_state=len(spin_1_occupation)-spin_1_full_local_state
spin_2_empty_local_state=len(spin_2_occupation)-spin_2_full_local_state

# find CB excitations
if input.trans:
	spin_1_CB_excitations=[]
	spin_2_CB_excitations=[]

	for i in range(len(spin_1_occupation)):
		if (spin_1_occupation[i] == 1):
			temp=copy.deepcopy(spin_1_occupation)
			temp[i]=0
			spin_1_CB_excitations.append(tuple(temp))

	for i in range(len(spin_2_occupation)):
		if (spin_2_occupation[i] == 1):
			temp=copy.deepcopy(spin_2_occupation)
			temp[i]=0
			spin_2_CB_excitations.append(tuple(temp))

# find VB excitations
if input.trans:
	spin_1_VB_excitations=[]
	spin_2_VB_excitations=[]

	for i in range(len(spin_1_occupation)):
		if (spin_1_occupation[i] == 0):
			temp=copy.deepcopy(spin_1_occupation)
			temp[i]=1
			spin_1_VB_excitations.append(tuple(temp))

	for i in range(len(spin_2_occupation)):
		if (spin_2_occupation[i] == 0):
			temp=copy.deepcopy(spin_2_occupation)
			temp[i]=1
			spin_2_VB_excitations.append(tuple(temp))

# calcuate the excitations input for vasp
if input.trans:
	incar=open('analyse.txt','w+')

	incar.write(str(len(spin_1_excitations)+len(spin_2_excitations))+'\t'+str(spin_1_full_local_state+spin_2_full_local_state)+'\t'+str(spin_1_empty_local_state+spin_2_empty_local_state)+'\t'+'NBANDS='+str(bands)+'\n')

	# local excitations
	for i in range(len(spin_1_excitations)):
		# make excitation occupation
		string1=str(spin_1_localstates[0]-1)+'*1.0 '
		for j in range(len(spin_1_excitations[i])):
			if spin_1_excitations[i][j]==1:
				string1+='1*1.0 '
			elif spin_1_excitations[i][j]==0:
				string1+='1*0.0 '
			else:
				print 'Error in analyse.py'
				error=1
		string1+=str(bands-spin_1_localstates[-1])+'*0.0 '
		incar.write('FERWE=')
		for k in range(kpoints):
			incar.write(string1)
		incar.write('\n')
		# keep default occupation
		string2=str(HOB_2)+'*1.0 '+str(bands-HOB_2)+'*0.0 '
		incar.write('FERDO=')
		for k in range(kpoints):
			incar.write(string2)
		incar.write('\n')

	for i in range(len(spin_2_excitations)):
		# keep default occupation
		string1=str(HOB_1)+'*1.0 '+str(bands-HOB_1)+'*0.0 '
		incar.write('FERWE=')
		for k in range(kpoints):
			incar.write(string1)
		incar.write('\n')
		# make excitation occupation
		string2=str(spin_2_localstates[0]-1)+'*1.0 '
		for j in range(len(spin_2_excitations[i])):
			if spin_2_excitations[i][j]==1:
				string2+='1*1.0 '
			elif spin_2_excitations[i][j]==0:
				string2+='1*0.0 '
			else:
				print 'Error in setting up the excitation'
				error=1
		string2+=str(bands-spin_2_localstates[-1])+'*0.0 '
		incar.write('FERDO=')
		for k in range(kpoints):
			incar.write(string2)
		incar.write('\n')

	# CB excitations
	for i in range(len(spin_1_CB_excitations)):
		# make excitation occupation
		string1=str(spin_1_localstates[0]-1)+'*1.0 '
		for j in range(len(spin_1_CB_excitations[i])):
			if spin_1_CB_excitations[i][j]==1:
				string1+='1*1.0 '
			elif spin_1_CB_excitations[i][j]==0:
				string1+='1*0.0 '
			else:
				print 'Error in analyse.py'
				error=1
		string1+='1*1.0 '+str(bands-spin_1_localstates[-1]-1)+'*0.0 '
		incar.write('FERWE=')
		for k in range(kpoints):
			incar.write(string1)
		incar.write('\n')
		# keep default occupation
		string2=str(HOB_2)+'*1.0 '+str(bands-HOB_2)+'*0.0 '
		incar.write('FERDO=')
		for k in range(kpoints):
			incar.write(string2)
		incar.write('\n')

	for i in range(len(spin_2_CB_excitations)):
		# keep default occupation
		string1=str(HOB_1)+'*1.0 '+str(bands-HOB_1)+'*0.0 '
		incar.write('FERWE=')
		for k in range(kpoints):
			incar.write(string1)
		incar.write('\n')
		# make excitation occupation
		string2=str(spin_2_localstates[0]-1)+'*1.0 '
		for j in range(len(spin_2_CB_excitations[i])):
			if spin_2_CB_excitations[i][j]==1:
				string2+='1*1.0 '
			elif spin_2_CB_excitations[i][j]==0:
				string2+='1*0.0 '
			else:
				print 'Error in setting up the excitation'
				error=1
		string2+='1*1.0 '+str(bands-spin_2_localstates[-1]-1)+'*0.0 '
		incar.write('FERDO=')
		for k in range(kpoints):
			incar.write(string2)
		incar.write('\n')

	# VB excitations
	for i in range(len(spin_1_VB_excitations)):
		# make excitation occupation
		string1=str(spin_1_localstates[0]-2)+'*1.0 '+'1*0.0 '
		for j in range(len(spin_1_VB_excitations[i])):
			if spin_1_VB_excitations[i][j]==1:
				string1+='1*1.0 '
			elif spin_1_VB_excitations[i][j]==0:
				string1+='1*0.0 '
			else:
				print 'Error in analyse.py'
				error=1
		string1+=str(bands-spin_1_localstates[-1])+'*0.0 '
		incar.write('FERWE=')
		for k in range(kpoints):
			incar.write(string1)
		incar.write('\n')
		# keep default occupation
		string2=str(HOB_2)+'*1.0 '+str(bands-HOB_2)+'*0.0 '
		incar.write('FERDO=')
		for k in range(kpoints):
			incar.write(string2)
		incar.write('\n')

	for i in range(len(spin_2_VB_excitations)):
		# keep default occupation
		string1=str(HOB_1)+'*1.0 '+str(bands-HOB_1)+'*0.0 '
		incar.write('FERWE=')
		for k in range(kpoints):
			incar.write(string1)
		incar.write('\n')
		# make excitation occupation
		string2=str(spin_2_localstates[0]-2)+'*1.0 '+'1*0.0 '
		for j in range(len(spin_2_VB_excitations[i])):
			if spin_2_VB_excitations[i][j]==1:
				string2+='1*1.0 '
			elif spin_2_VB_excitations[i][j]==0:
				string2+='1*0.0 '
			else:
				print 'Error in setting up the excitation'
				error=1
		string2+=str(bands-spin_2_localstates[-1])+'*0.0 '
		incar.write('FERDO=')
		for k in range(kpoints):
			incar.write(string2)
		incar.write('\n')

	# Save some extra info
	incar.write('\n')
	incar.write('=== Analyse.py_info ===\n')
	incar.write('spin 1 HOB: '+str(HOB_1)+'\n')
	incar.write('spin 2 HOB: '+str(HOB_2)+'\n')
	incar.write('\n')
	incar.write('spin 1 local states: '+str(spin_1_localstates)+'\n')
	incar.write('spin 2 local states: '+str(spin_2_localstates)+'\n')
	incar.write('\n')
	incar.write('spin 1 occupation: '+str(spin_1_occupation)+'\n')
	incar.write('spin 2 occupation: '+str(spin_2_occupation)+'\n')
	incar.write('\n')
	incar.write('spin 1 excitations: '+str(spin_1_excitations)+'\n')
	incar.write('spin 2 excitations: '+str(spin_2_excitations)+'\n')
	incar.write('\n')
	incar.write('spin 1 CB excitations: '+str(spin_1_CB_excitations)+'\n')
	incar.write('spin 2 CB excitations: '+str(spin_2_CB_excitations)+'\n')
	incar.write('\n')
	incar.write('spin 1 VB excitations: '+str(spin_1_VB_excitations)+'\n')
	incar.write('spin 2 VB excitations: '+str(spin_2_VB_excitations)+'\n')
	incar.write('\n')

	# Save which transition: name, spin, from band, to band
	incar.write('transitions\n')
	for i in range(len(spin_1_excitations)):
		from_band=0
		to_band=0
		for j in range(len(spin_1_occupation)):
			difference=spin_1_occupation[j]-spin_1_excitations[i][j]
			if difference == 1:
				from_band=spin_1_localstates[j]
			elif difference == -1:
				to_band=spin_1_localstates[j]
		incar.write('l'+str(i+1)+'\t'+str(1)+'\t'+str(from_band)+'\t'+str(to_band)+'\n')

	for i in range(len(spin_2_excitations)):
		from_band=0
		to_band=0
		for j in range(len(spin_2_occupation)):
			difference=spin_2_occupation[j]-spin_2_excitations[i][j]
			if difference == 1:
				from_band=spin_2_localstates[j]
			elif difference == -1:
				to_band=spin_2_localstates[j]
		incar.write('l'+str(i+1+len(spin_1_excitations))+'\t'+str(2)+'\t'+str(from_band)+'\t'+str(to_band)+'\n')

	for i in range(len(spin_1_CB_excitations)):
		from_band=0
		to_band=spin_1_localstates[-1]+1
		for j in range(len(spin_1_occupation)):
			difference=spin_1_occupation[j]-spin_1_CB_excitations[i][j]
			if difference == 1:
				from_band=spin_1_localstates[j]
		incar.write('c'+str(i+1)+'\t'+str(1)+'\t'+str(from_band)+'\t'+str(to_band)+'\n')

	for i in range(len(spin_2_CB_excitations)):
		from_band=0
		to_band=spin_2_localstates[-1]+1
		for j in range(len(spin_2_occupation)):
			difference=spin_2_occupation[j]-spin_2_CB_excitations[i][j]
			if difference == 1:
				from_band=spin_2_localstates[j]
		incar.write('c'+str(i+1+len(spin_1_CB_excitations))+'\t'+str(2)+'\t'+str(from_band)+'\t'+str(to_band)+'\n')

	for i in range(len(spin_1_VB_excitations)):
		from_band=spin_1_localstates[0]-1
		to_band=0
		for j in range(len(spin_1_occupation)):
			difference=spin_1_occupation[j]-spin_1_VB_excitations[i][j]
			if difference == -1:
				to_band=spin_1_localstates[j]
		incar.write('v'+str(i+1)+'\t'+str(1)+'\t'+str(from_band)+'\t'+str(to_band)+'\n')

	for i in range(len(spin_2_VB_excitations)):
		from_band=spin_2_localstates[0]-1
		to_band=0
		for j in range(len(spin_2_occupation)):
			difference=spin_2_occupation[j]-spin_2_VB_excitations[i][j]
			if difference == -1:
				to_band=spin_2_localstates[j]
		incar.write('v'+str(i+1+len(spin_1_VB_excitations))+'\t'+str(2)+'\t'+str(from_band)+'\t'+str(to_band)+'\n')

	incar.close()

	print 'analyse.txt done'

# Check which different charges are available
if input.charge:
	charge=open('charge.txt','w+')

	charge_string=""
	charges=0

	# -- charge
	if (spin_1_empty_local_state+spin_2_empty_local_state >= 2):
		charge_string+="n2"+'\t'+"NELECT="+str(electrons+2)+'\n'
		charges+=1

	# - charge
	if (spin_1_empty_local_state+spin_2_empty_local_state >= 1):
		charge_string+="n1"+'\t'+"NELECT="+str(electrons+1)+'\n'
		charges+=1

	# + charge
	if (spin_1_full_local_state+spin_2_full_local_state >= 1):
		charge_string+="p1"+'\t'+"NELECT="+str(electrons-1)+'\n'
		charges+=1

	# ++ charge
	if (spin_1_full_local_state+spin_2_full_local_state >= 2):
		charge_string+="p2"+'\t'+"NELECT="+str(electrons-2)+'\n'
		charges+=1

	charge.write(str(charges)+'\t'+'NBANDS='+str(bands)+'\n')
	charge.write(charge_string)

	charge.close()

	print 'charge.txt done'

# hyperfine
if input.hyperfine:
	hyperfine=open('hyperfine.txt','w+')
	gyromagnetic_ratio={}
	number_of_ratios=[]
	try:
		easyspin=open('../easyspin_isotopedata.txt','r')
	except:
		easyspin=open('../../../easyspin_isotopedata.txt','r')
	easyspin_lines=easyspin.readlines()

	# if spin is 0 no hyperfine
	if (HOB_1 != HOB_2):
		number_of_hyperfines=1
		for atom in atoms:
			atom_gyromagnetic_ratio=[]
			# find gyromatic ratio for atom from easyspin_data
			for line in easyspin_lines:
				if (line[0] != '%'):
					split_line=line.split()
					if (split_line[3] == atom and split_line[6] != '0.0'):
						# split_line[6] is the g factor
						# convert g factor to Mhz T^-1 (multiply with 7.622593285, the nulcear magnetron)
						atom_gyromagnetic_ratio.append(7.622593285*float(split_line[6]))
			gyromagnetic_ratio[atom]=atom_gyromagnetic_ratio
			number_of_ratios.append(len(gyromagnetic_ratio[atom]))
			number_of_hyperfines=number_of_hyperfines*len(gyromagnetic_ratio[atom])

		# make all possible hyperfines
		temp=[]
		for i in range(len(number_of_ratios)):
			for j in range(number_of_ratios[i]):
				temp.append(j)
		
		# all combinations
		output=list(set(itertools.combinations(temp,len(number_of_ratios))))

		# find extra that should not be included
		extra=[]
		for i in range(len(output)):
			for j in range(len(number_of_ratios)):
				if output[i][j] >= number_of_ratios[j]:
					extra.append(output[i])
					break

		# remove extra
		for e in extra:
			output.remove(e)
		
		# write data
		hyperfine.write(str(number_of_hyperfines)+'\t'+'NBANDS='+str(bands)+'\t'+'NELECT='+str(electrons)+'\n')
		
		for i in range(len(output)):
			hyperfine_string=""
			for j in range(len(number_of_ratios)):
				hyperfine_string+=str(gyromagnetic_ratio[atoms[j]][output[i][j]])
				hyperfine_string+=" "
			hyperfine.write("NGYROMAG= "+hyperfine_string+'\n')

	hyperfine.close()

	print 'hyperfine.txt done'

# zero-field splitting
if input.D:
	d_tensor=open('D-tensor.txt','w+')

	# need unpair electrons to produce a zero field splitting
	if (HOB_1 != HOB_2):
		d_tensor.write(str(HOB_1-HOB_2))

	d_tensor.close()

	print 'D-tensor.txt done'

# which bands should be visualized
if input.parchg:
	parchg=open('parchg.txt','w+')

	max_visband=[]
	min_visband=[]

	# add CB and VB to be visulaized
	if (spin_1_localstates != [] or spin_1_localstates != []):
		max_visband=[max(spin_1_localstates+spin_2_localstates)+1]
		min_visband=[min(spin_1_localstates+spin_2_localstates)-1]

	# visulaize all local states
	vis_bands=list(set().union(min_visband,spin_1_localstates,spin_2_localstates,max_visband))
	vis_bands.sort()

	if (vis_bands != []):
		parchg.write('IBAND = ')

	for b in vis_bands:
		parchg.write(str(b)+' ')

	if (vis_bands != []):
		parchg.write('\n')

	parchg.close()

	print 'parchg.txt done'

# Check which different spin should be calcuated
if input.spin:
	spin=open('spin.txt','w+')

	if (len(spin_1_localstates) == 0 or len(spin_2_localstates) == 0):
		spin.write("0")
		spin.close()
		print 'spin.txt done'
		sys.exit(error)

	if (len(spin_1_localstates) == len(spin_2_localstates) and (spin_1_occupation.count(1)+spin_2_occupation.count(1)) == 1 ):
		# symmetric, only one electron
		spin.write("0")
		spin.close()
		print 'spin.txt done'
		sys.exit(error)

	if (len(spin_1_localstates) == len(spin_2_localstates)) and (HOB_1-HOB_2-2 == HOB_2-HOB_1):
		# symmetric, no need to do move spin
		spin.write("0")
		spin.close()
		print 'spin.txt done'
		sys.exit(error)

	move_electron_from_spin_1_to_spin_2=True
	move_electron_from_spin_2_to_spin_1=True

	if (HOB_1 not in spin_1_localstates):
		move_electron_from_spin_1_to_spin_2=False

	if (HOB_2 not in spin_2_localstates):
		move_electron_from_spin_2_to_spin_1=False

	if (HOB_1+1 not in spin_1_localstates):
		move_electron_from_spin_2_to_spin_1=False

	if (HOB_2+1 not in spin_2_localstates):
		move_electron_from_spin_1_to_spin_2=False

	if (len(spin_1_localstates) == len(spin_2_localstates)) and (HOB_1 == HOB_2):
		# symmetric, spin zero move only one spin
		move_electron_from_spin_1_to_spin_2=False

	if (not move_electron_from_spin_1_to_spin_2 and not move_electron_from_spin_2_to_spin_1):
		spin.write("0")
		spin.close()
		print 'spin.txt done'
		sys.exit(error)

	if (move_electron_from_spin_1_to_spin_2 and move_electron_from_spin_2_to_spin_1):
		spin_1_energy_increase=spin_1_energy[HOB_1+1]-spin_1_energy[HOB_1]
		spin_1_energy_decrease=spin_1_energy[HOB_1-1]-spin_1_energy[HOB_1]
		spin_2_energy_increase=spin_2_energy[HOB_2+1]-spin_2_energy[HOB_2]
		spin_2_energy_decrease=spin_2_energy[HOB_2-1]-spin_2_energy[HOB_2]

		spin.write("2"+'\t'+'NBANDS='+str(bands)+'\n')
		
		if (spin_1_energy_increase+spin_2_energy_decrease < spin_2_energy_increase+spin_1_energy_decrease):
			# print lowest energy frist
			string1=str(HOB_1+1)+"*1.0 "+str(bands-(HOB_1+1))+"*0.0 "
			spin.write('FERWE=')
			for k in range(kpoints):
				spin.write(string1)
			spin.write('\n')
			string2=str(HOB_2-1)+"*1.0 "+str(bands-(HOB_2-1))+"*0.0 "
			spin.write('FERDO=')
			for k in range(kpoints):
				spin.write(string2)
			spin.write('\n')
			# print next
			string1=str(HOB_1-1)+"*1.0 "+str(bands-(HOB_1-1))+"*0.0 "
			spin.write('FERWE=')
			for k in range(kpoints):
				spin.write(string1)
			spin.write('\n')
			string2=str(HOB_2+1)+"*1.0 "+str(bands-(HOB_2+1))+"*0.0 "
			spin.write('FERDO=')
			for k in range(kpoints):
				spin.write(string2)
			spin.write('\n')
			spin.close()
			print 'spin.txt done'
			sys.exit(error)
		else:
			# print lowest energy frist
			string1=str(HOB_1-1)+"*1.0 "+str(bands-(HOB_1-1))+"*0.0 "
			spin.write('FERWE=')
			for k in range(kpoints):
				spin.write(string1)
			spin.write('\n')
			string2=str(HOB_2+1)+"*1.0 "+str(bands-(HOB_2+1))+"*0.0 "
			spin.write('FERDO=')
			for k in range(kpoints):
				spin.write(string2)
			spin.write('\n')
			# print next
			string1=str(HOB_1+1)+"*1.0 "+str(bands-(HOB_1+1))+"*0.0 "
			spin.write('FERWE=')
			for k in range(kpoints):
				spin.write(string1)
			spin.write('\n')
			string2=str(HOB_2-1)+"*1.0 "+str(bands-(HOB_2-1))+"*0.0 "
			spin.write('FERDO=')
			for k in range(kpoints):
				spin.write(string2)
			spin.write('\n')
			spin.close()
			print 'spin.txt done'
			sys.exit(error)

	if (move_electron_from_spin_1_to_spin_2):
		spin.write("1"+'\t'+'NBANDS='+str(bands)+'\n')
		string1=str(HOB_1-1)+"*1.0 "+str(bands-(HOB_1-1))+"*0.0 "
		spin.write('FERWE=')
		for k in range(kpoints):
			spin.write(string1)
		spin.write('\n')
		string2=str(HOB_2+1)+"*1.0 "+str(bands-(HOB_2+1))+"*0.0 "
		spin.write('FERDO=')
		for k in range(kpoints):
			spin.write(string2)
		spin.write('\n')
		spin.close()
		print 'spin.txt done'
		sys.exit(error)

	spin.close()
	print 'spin.txt done'

# exit script
sys.exit(error)
