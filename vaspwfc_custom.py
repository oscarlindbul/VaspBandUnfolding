#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from math import sqrt
from .vasp_constant import *
from multiprocessing import cpu_count
from scipy.fftpack import fftfreq, fftn, ifftn
from .vaspwfc import vaspwfc as oldwfc


############################################################


def save2vesta(phi=None, poscar='POSCAR', prefix='wfc',
			   lgam=False, lreal=False, ncol=10):
	'''
	Save the real space pseudo-wavefunction as vesta format.
	'''
	nx, ny, nz = phi.shape
	try:
		pos = open(poscar, 'r')
		head = ''
		for line in pos:
			if line.strip():
				head += line
			else:
				break
		head += '\n%5d%5d%5d\n' % (nx, ny, nz)
	except:
		raise IOError('Failed to open %s' % poscar)

	# Faster IO
	nrow = phi.size // ncol
	nrem = phi.size % ncol
	fmt = "%16.8E"

	psi = phi.copy()
	psi = psi.flatten(order='F')
	psi_h = psi[:nrow * ncol].reshape((nrow, ncol))
	psi_r = psi[nrow * ncol:]

	with open(prefix + '_r.vasp', 'w') as out:
		out.write(head)
		out.write(
			'\n'.join([''.join([fmt % xx for xx in row])
					   for row in psi_h.real])
		)
		out.write("\n" + ''.join([fmt % xx for xx in psi_r.real]))
	if not (lgam or lreal):
		with open(prefix + '_i.vasp', 'w') as out:
			out.write(head)
			out.write(
				'\n'.join([''.join([fmt % xx for xx in row])
						   for row in psi_h.imag])
			)
			out.write("\n" + ''.join([fmt % xx for xx in psi_r.imag]))

############################################################


class vaspwfc(oldwfc):
	'''
	Class for processing VASP Pseudowavefunction stored in separated WAV files.

	The format is the same as the standard WAVECAR, but is split up into
	files based on either spin, ion relaxation step, k-points or bands.
	The user may have specified that a smaller number of bands should be printed
	as well.

	'''

	def __init__(self, fnm='WAVEHEAD', E_files=["EIG_S1","EIG_S2"], lsorbit=False, lgamma=None, gamma_half='x', file_check=True, omp_num_threads=1):
		'''
		Initialization.
		'''

		self._fname = fnm
		self._efiles = E_files
		# the directory containing the input file
		self._dname = os.path.dirname(fnm)
		self._efiles = [self._dname + "/" + s for s in E_files]
		if self._dname == '':
			self._dname = '.'

		self._lsoc = lsorbit
		self._lgam = lgamma
		self._gam_half = gamma_half.lower()

		# It seems that some modules in scipy uses OPENMP, it is therefore
		# desirable to set the OMP_NUM_THREADS to tune the parallization.
		os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)

		# read the basic information
		self.readWFHeader()
		# read the band information
		self.readWFBand()

		if self._lsoc:
			assert self._nspin == 1, "NSPIN = 1 for noncollinear version WAVECAR!"
		assert not (lsorbit and lgamma), 'The two settings conflict!'
		assert self._gam_half == 'x' or self._gam_half == 'z', \
			'Gamma_half must be "x" or "z"'

		# calculate record positions of each index combination
		self.calc_rec_pos()
		self._missing_configs = set()
		if file_check:
			self.check_format()
	
	def get_file_format(self, inverse=False):
		if not inverse:
			file_format = "SKB"
			if self._wav_format != "SINGLE":
				for c in self._wav_format:
					file_format = file_format.replace(c, "")
			return file_format
		else:
			return self._wav_format.replace("SINGLE", "")

	#def get_config_ID(self, ion=-1, ispin=-1, ikpt=-1, iband=-1):
	#	file_format = self.get_file_format()
	#	ID = []
	#	ID.append(I+1)
	#	if "S" in file_format:
	#		ID.append(ispin)

	def get_config_name(self, ion, ispin, ikpt, iband, form=None):
		if form is None:
			form = self.get_file_format()
		s = "WAV"
		inds = [ispin,ikpt,iband]
		header ="SKB"
		if self._nion > 1:
			s += "_I{}".format(ion)
		for ind,a in zip(inds, header):
			if a in form:
				s += "_{}{}".format(a,ind)
		return s

	def check_format(self):
		file_format = self.get_file_format()
		
		I_range = [self._nion] 
		S_range = [self._nspin] if "S" in file_format else [-1, 0]
		K_range = [self._nkpts] if "K" in file_format else [-1, 0]
		for I in range(*I_range):
			for S in range(*S_range):
				if self._nbands[S] <= 0:
					continue
				for K in range(*K_range):
					B_range = [self._nbands[S]] if "B" in file_format else [-1,0]
					for B in range(*B_range):
						ID_s = self.get_config_name(I+1,S+1,K+1,B)
						ID = []
						ID.append(I+1)
						if S >= 0:
							ID.append(S+1)
						if K >= 0:
							ID.append(K+1)
						band = 0
						if B >= 0:
							if self._bandnum[S] is None:
								band = B+1
							else:
								band = self._bandnum[S][B]
							ID.append(band)
						ID_s = self.get_config_name(I+1,S+1,K+1,band,form=file_format)
						ID_s = self._dname + "/" + ID_s
						if not os.path.exists(ID_s):
							self._missing_configs.add(tuple(ID))
		
		if len(self._missing_configs) != 0:
			print("WARNING: the number of wavefunction files does not match given header info.")
			print("Missing configurations:")
			print("\t".join(file_format))

			for ID in self._missing_configs:
				ind = 0
				s = []
				for i,a in enumerate(file_format):
					s.append(str(ID[i]))
				s = "\t".join(s)
				print(s)
			print("Configurations will be limited, proceed with caution.")

	def calc_rec_pos(self):
		file_format = self.get_file_format(inverse=True)
		self._recpos = np.zeros((self._nspin, self._nkpts),dtype=int)
		pos = 0
		for S in range(self._nspin):
			if self._nbands[S] <= 0:
				continue
			if "S" not in file_format:
				pos = 0
			if "K" in file_format:
				for K in range(self._nkpts):
					self._recpos[S,K] = pos
					if "B" in file_format:
						pos += self._nbands[S]
					else:
						pos += 1
			else:
				if "B" in file_format:
					self._recpos[S,:] = np.ones(self._nkpts)*pos
					pos += self._nbands[S]
				else:
					self._recpos[S,:] = np.ones(self._nkpts)*pos
					pos += 1
	
	def readWFHeader(self):
		'''
		Read the system information from WAVHEAD header file
		'''

		self._Acell = np.zeros((3,3))
		self._kvecs = None
		self._bandnum = None #[None for i in range(self._nspin)]
		with open(self._fname, "r") as reader:
			for l in reader.readlines():
				l = l.strip()
				if "RECORD LENGTH" in l:
					self._recl = int(l.split(":")[-1])
				elif "NSPIN" in l:
					self._nspin = int(l.split(":")[-1])
					self._bandnum = [None for i in range(self._nspin)]
				elif "RTAG" in l:
					self._rtag = int(l.split(":")[-1])
					self._WFPrec = self.setWFPrec()
				elif "NKPOINTS" in l:
					self._nkpts = int(l.split(":")[-1])
				elif "BANDS" == l.split(":")[0]:
					self._nbands = list(map(int,l.split()[1:]))
				elif "ENMAX" in l:
					self._encut = float(l.split(":")[-1])
				elif "LATTVEC" in l:
					vals = l.split()
					num = int(vals[0][-2]) - 1
					vec = np.array(list(map(float, vals[1:])))
					self._Acell[num,:] = vec
				elif "BINARY FORMAT" in l:
					binformat = l.split(":")[-1].strip()
					gam_val = None
					if binformat == "standard":
						gam_val = False
					elif binformat == "gamma":
						gam_val = True
					else:
						raise Exception("Binary format specified in header is in incorrect format, must be either 'standard' or 'gamma'")
					if self._lgam is None:
						self._lgam = gam_val
					if gam_val != self._lgam:
						print("Provided gamma setting is inconsistent with value found in header, continuing with user setting")
				elif "NPLANEWAVE" in l:
					self._nplws = list(map(int,l.split(":")[1:]))
				elif "KVEC" in l:
					if self._kvecs is None:
						self._kvecs = np.zeros((self._nkpts, 3))
					vals = l.split()
					num = int(vals[0][-2]) - 1
					vec = np.array(list(map(float, vals[1:])))
					self._kvecs[num,:] = vec
				elif "UPBANDS" in l or "DOWNBANDS" in l:
					spin = 0 if "UP" in l else 1
					bands = list(map(int,l.split()[1:]))
					self._bandnum[spin] = np.array(bands)
				#	self._bands[spin] = np.zeros((self._nkpts, self._nbands[spin]), dtype=float)
				#	self._occs[spin] = np.zeros((self.nkpts, self._nbands[spin]), dtype=float)
				elif "WAV FORMAT" in l:
					self._wav_format = l.split(":")[-1].strip()
				
		self._Omega = np.linalg.det(self._Acell)
		self._Bcell = np.linalg.inv(self._Acell).T
		
		Anom = np.linalg.norm(self._Acell, axis=1)
		CUTOF = np.ceil(
				np.sqrt(self._encut / RYTOEV) / (TPI / (Anom /AUTOA))
				)
		self._ngrid = np.array(2*CUTOF + 1, dtype=int)


		# goto the start of the file and read the first record
		#self._wfc.seek(0)
		#self._recl, self._nspin, self._rtag = np.array(
		#	np.fromfile(self._wfc, dtype=np.float, count=3),
		#	dtype=np.int64
		#)
		#self._WFPrec = self.setWFPrec()
		## the second record
		#self._wfc.seek(self._recl)
		#dump = np.fromfile(self._wfc, dtype=np.float, count=12)

		#self._nkpts = int(dump[0])					   # No. of k-points
		#self._nbands = int(dump[1])						# No. of bands
		#self._encut = dump[2]						   # Energy cutoff
		## real space supercell basis
		#self._Acell = dump[3:].reshape((3, 3))
		## real space supercell volume
		#self._Omega = np.linalg.det(self._Acell)
		## reciprocal space supercell volume
		#self._Bcell = np.linalg.inv(self._Acell).T

		## Minimum FFT grid size
		#Anorm = np.linalg.norm(self._Acell, axis=1)
		#CUTOF = np.ceil(
		#	sqrt(self._encut / RYTOEV) / (TPI / (Anorm / AUTOA))
		#)
		#self._ngrid = np.array(2 * CUTOF + 1, dtype=int)

	def readWFBand(self):
		'''
		Extract KS energies and Fermi occupations from WAVECAR.
		'''
		
		self._kvecs = np.zeros((self._nkpts, 3), dtype=float)
		self._bands = [None for i in range(self._nspin)]
		self._occs = [None for i in range(self._nspin)]

		#self._nplws = np.zeros(self._nkpts, dtype=int)
		#self._kvecs = np.zeros((self._nkpts, 3), dtype=float)
		#self._bands = np.zeros(
		#	(self._nspin, self._nkpts, self._nbands), dtype=float)
		#self._occs = np.zeros(
		#	(self._nspin, self._nkpts, self._nbands), dtype=float)

		for spin in range(self._nspin):
			if self._nbands[spin] <= 0:
				continue
			e_file = self._efiles[spin]
			with open(e_file, "r") as reader:
				lines = reader.readlines()
			line_ind = 0
			N_ion = 0
			band_E = []
			band_occ = []
			while line_ind < len(lines):
				if "I = " in lines[line_ind]:
					line_ind += 1
					E_by_bands = []
					occ_by_bands = []
					for i in range(self._nbands[spin]):
						E_occ_vals = list(map(float,lines[line_ind].split()))
						E_vals = E_occ_vals[::2]
						occ_vals = E_occ_vals[1::2]
						E_by_bands.append(E_vals)
						occ_by_bands.append(occ_vals)
						line_ind += 1
					N_ion += 1
					band_E.append(np.array(E_by_bands).T)
					band_occ.append(np.array(occ_by_bands).T)
				line_ind += 1
			self._bands[spin] = np.zeros((N_ion, self._nkpts, self._nbands[spin]), dtype=float)
			self._occs[spin] = np.zeros((N_ion, self._nkpts, self._nbands[spin]), dtype=float)
			self._nion = N_ion
			for n in range(N_ion):
				self._bands[spin][n,:,:] = band_E[n]
				self._occs[spin][n,:,:] = band_occ[n]

		#for ii in range(self._nspin):
		#	for jj in range(self._nkpts):
		#		rec = self.whereRec(ii+1, jj+1, 1) - 1
		#		self._wfc.seek(rec * self._recl)
		#		dump = np.fromfile(self._wfc, dtype=np.float,
		#						   count=4+3*self._nbands)
		#		if ii == 0:
		#			self._nplws[jj] = int(dump[0])
		#			self._kvecs[jj] = dump[1:4]
		#		dump = dump[4:].reshape((-1, 3))
		#		self._bands[ii, jj, :] = dump[:, 0]
		#		self._occs[ii, jj, :] = dump[:, 2]

		if self._nkpts > 1:
			tmp = np.linalg.norm(
				np.dot(np.diff(self._kvecs, axis=0), self._Bcell), axis=1)
			self._kpath = np.concatenate(([0, ], np.cumsum(tmp)))
		else:
			self._kpath = None
		return self._kpath, self._bands

	def wfc_r(self, ion=None, ispin=1, ikpt=1, iband=1,
			  gvec=None, Cg=None, ngrid=None,
			  rescale=None,
			  norm=True):
		'''
		Obtain the pseudo-wavefunction of the specified KS states in real space
		by performing FT transform on the reciprocal space planewave
		coefficients.  The 3D FT grid size is determined by ngrid, which
		defaults to self._ngrid if not given.  Gvectors of the KS states is used
		to put 1D planewave coefficients back to 3D grid.

		Inputs:
			ispin : spin index of the desired KS states, starting from 1
			ikpt  : k-point index of the desired KS states, starting from 1
			iband : band index of the desired KS states, starting from 1
			gvec  : the G-vectors correspond to the plane-wave coefficients
			Cg	  : the plane-wave coefficients. If None, read from WAVECAR
			ngrid : the FFT grid size
			norm  : normalized Cg?

		The return wavefunctions are normalized in a way that

						\sum_{ijk} | \phi_{ijk} | ^ 2 = 1

		'''
		self.checkIndex(ion, ispin, ikpt, iband)

		if ngrid is None:
			ngrid = self._ngrid.copy() * 2
		else:
			ngrid = np.array(ngrid, dtype=int)
			assert ngrid.shape == (3,)
			assert np.alltrue(ngrid >= self._ngrid), \
				"Minium FT grid size: (%d, %d, %d)" % \
				(self._ngrid[0], self._ngrid[1], self._ngrid[2])

		# The default normalization of np.fft.fftn has the direct transforms
		# unscaled and the inverse transforms are scaled by 1/n. It is possible
		# to obtain unitary transforms by setting the keyword argument norm to
		# "ortho" (default is None) so that both direct and inverse transforms
		# will be scaled by 1/\sqrt{n}.

		# default normalization factor so that
		# \sum_{ijk} | \phi_{ijk} | ^ 2 = 1
		normFac = rescale if rescale is not None else np.sqrt(np.prod(ngrid))

		if gvec is None:
			gvec = self.gvectors(ikpt)

		if self._lgam:
			if self._gam_half == 'z':
				phi_k = np.zeros(
					(ngrid[0], ngrid[1], ngrid[2]//2 + 1), dtype=np.complex128)
			else:
				phi_k = np.zeros(
					(ngrid[0]//2 + 1, ngrid[1], ngrid[2]), dtype=np.complex128)
		else:
			phi_k = np.zeros(ngrid, dtype=np.complex128)

		gvec %= ngrid[np.newaxis, :]

		if self._lsoc:
			wfc_spinor = []
			if Cg:
				dump = Cg
			else:
				dump = self.readBandCoeff(ion, ispin, ikpt, iband, norm, check=False)
			nplw = dump.shape[0] // 2

			# spinor up
			phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = dump[:nplw]
			wfc_spinor.append(ifftn(phi_k) * normFac)
			# spinor down
			phi_k[:, :, :] = 0.0j
			phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = dump[nplw:]
			wfc_spinor.append(ifftn(phi_k) * normFac)

			del dump
			return wfc_spinor

		else:
			if Cg is not None:
				phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = Cg
			else:
				phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = self.readBandCoeff(ion, ispin, ikpt, iband, norm, check=False)

			if self._lgam:
				# add some components that are excluded and perform c2r FFT
				if self._gam_half == 'z':
					for ii in range(ngrid[0]):
						for jj in range(ngrid[1]):
							fx = ii if ii < ngrid[0] // 2 + \
								1 else ii - ngrid[0]
							fy = jj if jj < ngrid[1] // 2 + \
								1 else jj - ngrid[1]
							if (fy > 0) or (fy == 0 and fx >= 0):
								continue
							phi_k[ii, jj, 0] = phi_k[-ii, -jj, 0].conjugate()

					phi_k /= np.sqrt(2.)
					phi_k[0, 0, 0] *= np.sqrt(2.)
					return np.fft.irfftn(phi_k, s=ngrid) * normFac
				elif self._gam_half == 'x':
					for jj in range(ngrid[1]):
						for kk in range(ngrid[2]):
							fy = jj if jj < ngrid[1] // 2 + \
								1 else jj - ngrid[1]
							fz = kk if kk < ngrid[2] // 2 + \
								1 else kk - ngrid[2]
							if (fy > 0) or (fy == 0 and fz >= 0):
								continue
							phi_k[0, jj, kk] = phi_k[0, -jj, -kk].conjugate()

					phi_k /= np.sqrt(2.)
					phi_k[0, 0, 0] *= np.sqrt(2.)
					phi_k = np.swapaxes(phi_k, 0, 2)
					tmp = np.fft.irfftn(
						phi_k, s=(ngrid[2], ngrid[1], ngrid[0])) * normFac
					return np.swapaxes(tmp, 0, 2)
			else:
				# perform complex2complex FFT
				return ifftn(phi_k * normFac)

	def write_gamma_wavecar(self, gamma_half="x", out="WAVECAR_gamma"):
		print("gamma conversion of custom format currently prohibited")
		return
		assert not self._lgam

		# get gvectors for gamma representation
		gamma_gvec = self.gvectors(force_Gamma=True, check_consistency=False)
		gvec = self.gvectors()

		# open wavecar file
		with open(out, "w") as new_wc:
			# write header details
			new_nplws = (self._nplws[0] - 1) // 2 + 1
			plws_rec_size = np.max(new_nplws)*np.dtype(self.setWFPrec()).itemsize
			band_rec_size = np.dtype(np.float64).itemsize*(self._nbands*3+1)
			# record needs to be large enough to contain both plane waves and bands
			new_rec_size = max(plws_rec_size, band_rec_size)
			nfloat = new_rec_size // 8 # number of float64s per record
			# header line
			rec = np.zeros(nfloat, dtype=np.float64)
			rec[0:3] = new_rec_size, self._nspin, self._rtag
			rec.tofile(new_wc)
			# header line 2 (nkpts, nbands, encut)
			rec[0:3] = 1, self._nbands, self._encut
			rec[3:3+9] = self._Acell.reshape((1,-1))
			rec.tofile(new_wc)
			wave_rec = np.zeros(new_nplws, dtype=self.setWFPrec())
			for spin in range(self._nspin):
				rec[0] = new_nplws
				rec[1:1+3] = self._kvecs[0]
				rec[4: 4+3*self._nbands : 3] = self._bands[spin, 0, :]
				rec[4+1 : 4+3*self._nbands : 3] = 0.0 # so far energies always real?
				rec[4+2 : 4+3*self._nbands : 3] = self._occs[spin, 0, :]
				rec.tofile(new_wc)
				
				ngrid = self._ngrid.copy()*2
				if gamma_half == "z":
					ordered_grid = ngrid
				elif gamma_half == "x":
					ordered_grid = ngrid[[2,1,0]]
				else:
					raise ValueError("Gamma reduction direction must be z or x")
				
				for band in range(self._nbands):
					# convert band function to real-space
					phi_r = self.wfc_r(spin+1, 1, band+1, norm=False)
					# if gamma-symmetry is to hold, function should be real and have same (total) density -> make real-space function fully real with the same single-state density
					phi_r = np.sqrt(phi_r.real**2 + phi_r.imag**2)*np.sign(phi_r.real)

					# fourier transform back to reciprocal space
					if gamma_half == "z":
						phi_k = np.fft.rfftn(phi_r, s=ordered_grid)
					else: # if x
						phi_r = np.swapaxes(phi_r, 0, 2)
						tmp = np.fft.rfftn(phi_r, s=ordered_grid)
						phi_k = np.swapaxes(tmp, 0, 2)

					# format components in correct format
					phi_k[0,0,0] /= np.sqrt(2.)
					k_coeffs = phi_k[gamma_gvec[:,0], gamma_gvec[:,1], gamma_gvec[:,2]]	
					k_coeffs *= np.sqrt(2.)
					
					# write to wavecar
					wave_rec[:new_nplws] = k_coeffs

					wave_rec.tofile(new_wc)

	def poisson(self, rho=None, iband=1, ikpt=1, ispin=1, ngrid=None, norm=False):
		"""
		Given a charge density "rho", solve the Poisson equation with periodic
		boundary condition to find out the corresponding electric potential and
		field.

		When "rho" is None, construct the charge density from a chosen Kohn-Sham
		state, i.e. rho(r) = phi_n(r).conj() * phi_n(r).

		In SI units, the real space Poisson equation:

					\nabla^2 V = - \rho / \varepsilon_0
							 E = - \nabla V

		the reciprocal space Poisson equation:

					G**2 * V_q = - rho_q / \varepsilon_0
						   E_q = -1j * G * V_q

		Note that the G=(0,0,0) entry is set to 1.0 instead of 0 to avoid
		divergence.
		"""

		if rho is not None:
			rho = np.asarray(rho)
			ngrid = np.array(rho.shape, dtype=int)
			assert ngrid.shape == (3,)
		else:
			ngrid = self._ngrid * 2
			# normalization factor so that
			# \sum_{ijk} | \phi_{ijk} | ^ 2 * volume / Ngrid = 1
			normFac = np.prod(ngrid) / self._Omega
			if self._lsoc:
				rho = np.zeros(ngrid, dtype=float)
				phi_spinor = self.wfc_r(iband=iband, ikpt=ikpt, ispin=ispin,
										ngrid=ngrid, norm=norm)
				# negative charges, hence the minus sign
				for phi in phi_spinor:
					rho += -(phi.conj() * phi).real * normFac
			else:
				phi = self.wfc_r(iband=iband, ikpt=ikpt, ispin=ispin,
								 ngrid=ngrid, norm=norm)
				# negative charges, hence the minus sign
				rho = -(phi.conj() * phi).real * normFac

		fx = [ii if ii < ngrid[0] // 2 + 1 else ii - ngrid[0]
			  for ii in range(ngrid[0])]
		fy = [jj if jj < ngrid[1] // 2 + 1 else jj - ngrid[1]
			  for jj in range(ngrid[1])]
		fz = [kk if kk < ngrid[2] // 2 + 1 else kk - ngrid[2]
			  for kk in range(ngrid[2])]

		# plane-waves: Reciprocal coordinate
		# indexing = 'ij' so that outputs are of shape (ngrid[0], ngrid[1], ngrid[2])
		Dx, Dy, Dz = np.meshgrid(fx, fy, fz, indexing='ij')
		# plane-waves: Cartesian coordinate
		Gx, Gy, Gz = np.tensordot(
			self._Bcell * np.pi * 2, [Dx, Dy, Dz], axes=(0, 0))
		# the norm squared of the G-vectors
		G2 = Gx**2 + Gy**2 + Gz**2
		# Note that the G=(0,0,0) entry is set to 1.0 instead of 0.
		G2[0, 0, 0] = 1.0

		# permittivity of vacuum [F / m]
		_eps0 = 8.85418781762039E-12
		# charge of one electron, in unit of Coulomb [1F * 1V]
		_e = 1.6021766208E-19

		# charge density in reciprocal space, rho in unit of [Coulomb / Angstrom**3]
		rho_q = np.fft.fftn(1E10 * _e * rho / _eps0, norm='ortho')
		# the electric potential in reciprocal space
		# V_q = -rho_q / (-G2)
		V_q = rho_q / G2
		# the electric potential in real space in unit of 'Volt'
		V_r = np.fft.ifftn(V_q, norm='ortho').real
		# the electric field in x/y/z in real space in unit of 'Volt / Angstrom'
		E_x = np.fft.ifftn(-1j * Gx * V_q, norm='ortho').real
		E_y = np.fft.ifftn(-1j * Gy * V_q, norm='ortho').real
		E_z = np.fft.ifftn(-1j * Gz * V_q, norm='ortho').real

		return rho, V_r, E_x, E_y, E_z

	def readBandCoeff(self, ion=None, ispin=1, ikpt=1, iband=1, norm=False, check=True):
		'''
		Read the planewave coefficients of specified KS states.
		'''
		if check:
			self.checkIndex(ion, ispin, ikpt, iband)

		rec = self.whereRec(ion, ispin, ikpt, iband, False)
		path = self._dname + "/" + self.get_config_name(ion,ispin,ikpt,iband)
		with open(path, "rb") as reader:
			reader.seek(rec * self._recl)

			nplw = self._nplws[ikpt - 1]
			dump = np.fromfile(reader, dtype=self._WFPrec, count=nplw)

		cg = np.asarray(dump, dtype=np.complex128)
		if norm:
			cg /= np.linalg.norm(cg)
		return cg

	def whereRec(self, ion=None, ispin=1, ikpt=1, iband=1, check=True):
		'''
		Return the rec position for specified KS state.
		'''
		if check:
			self.checkIndex(ion, ispin, ikpt, iband)

		rec = self._recpos[ispin-1,ikpt-1]
		if "B" in self.get_file_format(inverse=True):
			if self._bandnum[ispin-1] is None:
				band_loc = iband
			else:
				band_loc = np.where(self._bandnum[ispin-1] == iband)[0][0]
			rec += band_loc
		
		return rec

	def checkIndex(self, ion, ispin, ikpt, iband):
		'''
		Check if the index is valid!
		'''
		if ion is None:
			ion = self._nion
		assert 1 <= ion <= self._nion,	'Invalid ion step!'
		assert 1 <= ispin <= self._nspin and self._nbands[ispin-1] > 0,  'Invalid spin index!'
		assert 1 <= ikpt <= self._nkpts,  'Invalid kpoint index!'
		if self._bandnum[ispin-1] is None:
			assert 1 <= iband <= self._nbands[ispin-1]
		else:
			assert iband in self._bandnum[ispin-1], 'Invalid band index!'

		file_format = self.get_file_format()
		ID = []
		if self._nion > 1:
			ID.append(ion)
		if "S" in file_format:
			ID.append(ispin)
		if "K" in file_format:
			ID.append(ikpt)
		if "B" in file_format:
			ID.append(iband)
		assert tuple(ID) not in self._missing_configs, 'Wavefunction file for index configuration not present!'

	def TransitionDipoleMoment(self, ks_i, ks_j, norm=True,
							   realspace=False):
		'''
		calculate Transition Dipole Moment (TDM) between two KS states.

		If "realspace = False", the TDM will be evaluated in momentum space
		according to the formula in:
		https://en.wikipedia.org/wiki/Transition_dipole_moment

									i⋅h
		<psi_a | r | psi_b> =  -------------- ⋅ <psi_a | p | psi_b>
								m⋅(Eb - Ea)

									  2		   ____
									 h			╲
							=  ------------- ⋅	 ╲	 Cai⋅Cbi⋅Gi
								m⋅(Eb - Ea)		 ╱
												╱
											   ‾‾‾‾
												 i

		Otherwise, the TDM will be evaluated in real space.

		Note: |psi_a> and |psi_b> should be bloch function with
			  the same k vector.

		The KS states ks_i (ks_j) is specified by list of index (ispin, ikpt, iband).
		'''

		ks_i = list(ks_i)
		ks_j = list(ks_j)
		assert len(ks_i) == len(ks_j) == 3, 'Must be three indexes!'
		assert ks_i[1] == ks_j[1], 'k-point of the two states differ!'
		self.checkIndex(*ks_i)
		self.checkIndex(*ks_j)

		# energy differences between the two states
		E1 = self._bands[ks_i[0]-1, ks_i[1]-1, ks_i[2]-1]
		E2 = self._bands[ks_j[0]-1, ks_j[1]-1, ks_j[2]-1]
		dE = E2 - E1

		if realspace:
			fx = np.linspace(0, 1, self._ngrid[0], endpoint=False)
			fy = np.linspace(0, 1, self._ngrid[1], endpoint=False)
			fz = np.linspace(0, 1, self._ngrid[2], endpoint=False)

			Dx, Dy, Dz = np.meshgrid(fx, fy, fz, indexing='ij')
			Rx, Ry, Rz = np.tensordot(self._Acell, [Dx, Dy, Dz], axes=[0, 0])

			fac = np.sqrt(np.prod(self._ngrid) / self._Omega)
			phi_i = self.wfc_r(*ks_i, norm=True, ngrid=self._ngrid)
			phi_j = self.wfc_r(*ks_j, norm=True, ngrid=self._ngrid)

			pij = phi_i.conjugate() * phi_j
			tdm = np.array([
				np.sum(pij * Rx),
				np.sum(pij * Ry),
				np.sum(pij * Rz)
			])
			ovlap = pij.sum()
		else:
			# according to the above equation, G = 0 does NOT contribute to TDM.
			gvec = np.dot(self.gvectors(ikpt=ks_i[1]), self._Bcell*TPI)
			# planewave coefficients of the two states
			phi_i = self.readBandCoeff(*ks_i, norm=norm)
			phi_j = self.readBandCoeff(*ks_j, norm=norm)

			tmp1 = phi_i.conjugate() * phi_j
			ovlap = np.sum(tmp1)
			if self._lgam:
				tmp2 = phi_i * phi_j.conjugate()
				# according to the above equation, G = 0 does NOT contribute to TDM.
				tdm = (np.sum(tmp1[:, np.newaxis] * gvec, axis=0) -
					   np.sum(tmp2[:, np.newaxis] * gvec, axis=0)) / 2.
			else:
				tdm = np.sum(tmp1[:, np.newaxis] * gvec, axis=0)

			tdm = 1j / (dE / (2*RYTOEV)) * tdm * AUTOA * AUTDEBYE

		return E1, E2, dE, ovlap, tdm

	def TransitionDipoleMomentBetweenDifferentWAVECAR(self, other, ks_i, ks_j, norm=False):
		'''
		calculate Transition Dipole Moment between two KS states.
		TDM in momentum representation
											 ___			  
								  i⋅h		 ╲				  
		<psi_a| r | psi_b> =	--------- ⋅   ╲   Cai⋅Cbi⋅Gi
								 Eb - Ea	  ╱				  
											 ╱				  
											 ‾‾‾			  
											  i		  
		Note: |psi_a> and |psi_b> should be bloch function with 
			  the same k vector.

		The KS states ks_i (ks_j) is specified by list of index (ispin, ikpt, iband).
		'''
		
		ks_i = list(ks_i); ks_j = list(ks_j)
		assert len(ks_i) == len(ks_j) == 3, 'Must be there indexes!'
		assert ks_i[1] == ks_j[1], 'k-point of the two states differ!'
		self.checkIndex(*ks_i)
		self.checkIndex(*ks_j)

		# according to the above equation, G = 0 does NOT contribute to TDM.
		gvec = np.dot(self.gvectors(ikpt=ks_i[1]), self._Bcell*TPI)
		# planewave coefficients of the two states
		phi_i = self.readBandCoeff(*ks_i, norm=norm)
		phi_j = other.readBandCoeff(*ks_j, norm=norm)
		# energy differences between the two states
		dE = other._bands[ks_j[0]-1, ks_j[1]-1, ks_j[2]-1] - \
			 self._bands[ks_i[0]-1, ks_i[1]-1, ks_i[2]-1]

		tmp1 = phi_i.conjugate() * phi_j
		ovlap = np.sum(tmp1)
		if self._lgam:
			tmp2 = phi_i * phi_j.conjugate()
			tdm = (np.sum(tmp1[:,np.newaxis] * gvec, axis=0) -
				   np.sum(tmp2[:,np.newaxis] * gvec, axis=0)) / 2.
		else:
			tdm = np.sum(tmp1[:,np.newaxis] * gvec, axis=0)

		tdm = 1j / (dE / (2*RYTOEV)) * tdm * AUTOA * AUTDEBYE

		return dE, ovlap, tdm
	
	def inverse_participation_ratio(self, norm=True, bands=None):
		'''
		Calculate Inverse Paticipation Ratio (IPR) from the wavefunction. IPR is
		a measure of the localization of Kohn-Sham states. For a particular KS
		state \phi_j, it is defined as

							\sum_n |\phi_j(n)|^4
			IPR(\phi_j) = -------------------------
						  |\sum_n |\phi_j(n)|^2||^2

		where n iters over the number of grid points.

		bands = bands to calculate the IPR for (starting from index 1)
		'''

		self.ipr = np.zeros((self._nspin, self._nkpts, self._nbands, 3))

		if bands is None:
			bands = [i for i in range(self._nbands)]
		else:
			bands = [b - 1 for b in bands]
			assert len(bands) != 0, "Given list of bands must not be empty"
			assert min(bands) >= 0, "Invalid band indices given (index starting from 1)"

		for ispin in range(self._nspin):
			for ikpt in range(self._nkpts):
				for iband in bands:
					phi_j = self.wfc_r(ispin+1, ikpt+1, iband+1,
									   norm=norm)
					phi_j_abs = np.abs(phi_j)

					print('Calculating IPR of #spin %4d, #kpt %4d, #band %4d' %
						  (ispin+1, ikpt+1, iband+1))
					self.ipr[ispin, ikpt, iband,
							 0] = self._kpath[ikpt] if self._kpath is not None else 0
					self.ipr[ispin, ikpt, iband,
							 1] = self._bands[ispin, ikpt, iband]
					self.ipr[ispin, ikpt, iband, 2] = np.sum(
						phi_j_abs**4) / np.sum(phi_j_abs**2)**2

		np.save('ipr.npy', self.ipr)
		return self.ipr

	def elf(self, kptw, ngrid=None, warn=True):
		'''
		Calculate the electron localization function (ELF) from WAVECAR.

		The following formula was extracted from VASP ELF.F:
					 _
					 h^2	*	 2		T.........kinetic energy
		  T    =  -2 --- Psi grad Psi	T+TCORR...pos.definite kinetic energy
				   ^ 2 m				TBOS......T of an ideal Bose-gas
				   ^
				   I am not sure if we need to times 2 here, use 1 in this
				   script.

				   _								(=infimum of T+TCORR)
				 1 h^2		2			DH........T of hom.non-interact.e- - gas
		  TCORR= - ---	grad rho					(acc.to Fermi)
				 2 2 m					ELF.......electron-localization-function
				   _			 2
				 1 h^2 |grad rho|
		  TBOS = - --- ----------		D = T + TCORR - TBOS
				 4 2 m	  rho
				   _								\				 1
				 3 h^2		  2/3  5/3			=====>	  ELF = ------------
		  DH   = - --- (3 Pi^2)  rho				/					D	2
				 5 2 m											 1 + ( ---- )
																		DH

		REF:
			1. Nature, 371, 683-686 (1994)
			2. Becke and Edgecombe, J. Chem. Phys., 92, 5397(1990)
			3. M. Kohout and A. Savin, Int. J. Quantum Chem., 60, 875-882(1996)
			4. http://www2.cpfs.mpg.de/ELF/index.php?content=06interpr.txt
		'''

		if warn:
			warning = """
			###################################################################
			If you are using VESTA to view the resulting ELF, please rename the
			output file as ELFCAR, otherwise there will be some error in the
			isosurface plot!

			When CHG*/PARCHG/*.vasp are read in to visualize isosurfaces and
			sections, data values are divided by volume in the unit of bohr^3.
			The unit of charge densities input by VESTA is, therefore, bohr^−3.

			For LOCPOT/ELFCAR files, volume data are kept intact.

			You can turn off this warning by setting "warn=False" in the "elf"
			method.
			###################################################################
			"""
			print(warning)

		# the k-point weights
		kptw = np.array(kptw, dtype=float)
		assert kptw.shape == (self._nkpts,), "K-point weights must be provided \
											  to calculate charge density!"
		# normalization
		kptw /= kptw.sum()

		if ngrid is None:
			ngrid = self._ngrid * 2
		else:
			ngrid = np.array(ngrid, dtype=int)
			assert ngrid.shape == (3,)
			assert np.alltrue(ngrid >= self._ngrid), \
				"Minium FT grid size: (%d, %d, %d)" % \
				(self._ngrid[0], self._ngrid[1], self._ngrid[2])

		fx = [ii if ii < ngrid[0] // 2 + 1 else ii - ngrid[0]
			  for ii in range(ngrid[0])]
		fy = [jj if jj < ngrid[1] // 2 + 1 else jj - ngrid[1]
			  for jj in range(ngrid[1])]
		fz = [kk if kk < ngrid[2] // 2 + 1 else kk - ngrid[2]
			  for kk in range(ngrid[2])]

		# plane-waves: Reciprocal coordinate
		# indexing = 'ij' so that outputs are of shape (ngrid[0], ngrid[1], ngrid[2])
		Dx, Dy, Dz = np.meshgrid(fx, fy, fz, indexing='ij')
		# plane-waves: Cartesian coordinate
		Gx, Gy, Gz = np.tensordot(
			self._Bcell * np.pi * 2, [Dx, Dy, Dz], axes=(0, 0))
		# the norm squared of the G-vectors
		G2 = Gx**2 + Gy**2 + Gz**2
		# k-points vectors in Cartesian coordinate
		vkpts = np.dot(self._kvecs, self._Bcell * 2 * np.pi)

		# normalization factor so that
		# \sum_{ijk} | \phi_{ijk} | ^ 2 * volume / Ngrid = 1
		normFac = np.sqrt(np.prod(ngrid) / self._Omega)

		# electron localization function
		ElectronLocalizationFunction = []
		# Charge density
		rho = np.zeros(ngrid, dtype=complex)
		# Kinetic energy density
		tau = np.zeros(ngrid, dtype=complex)

		for ispin in range(self._nspin):
			# initialization
			rho[...] = 0.0
			tau[...] = 0.0

			for ikpt in range(self._nkpts):

				# plane-wave G-vectors
				igvec = self.gvectors(ikpt+1)
				# for gamma-only version, complete the missing -G vectors
				if self._lgam:
					tmp = np.array([-k for k in igvec[1:]], dtype=int)
					igvec = np.vstack([igvec, tmp])
				# plane-wave G-vectors in Cartesian coordinate
				rgvec = np.dot(igvec, self._Bcell * 2 * np.pi)

				k = vkpts[ikpt]						  # k
				gk = rgvec + k[np.newaxis, :]			# G + k
				gk2 = np.linalg.norm(gk, axis=1)**2		# | G + k |^2

				for iband in range(self._nbands):
					# omit the empty bands
					if self._occs[ispin, ikpt, iband] == 0.0:
						continue

					rspin = 2.0 if self._nspin == 1 else 1.0
					weight = rspin * kptw[ikpt] * \
						self._occs[ispin, ikpt, iband]

					# if self._lgam:
					#	  ########################################
					#	  # slower
					#	  ########################################
					#	  # wavefunction in real space
					#	  # VASP does NOT do normalization in elf.F
					#	  phi_r  = self.wfc_r(ispin=ispin+1, ikpt=ikpt+1,
					#						  iband=iband+1,
					#						  ngrid=ngrid,
					#						  norm=False) * normFac
					#	  # wavefunction in reciprocal space
					#	  phi_q  = np.fft.fftn(phi_r, norm='ortho')
					#	  # grad^2 \phi in reciprocal space
					#	  lap_phi_q = -gk2 * phi_q
					#	  # grad^2 \phi in real space
					#	  lap_phi_r = np.fft.ifftn(lap_phi_q, norm='ortho')
					# else:

					########################################
					# faster
					########################################
					# wavefunction in reciprocal space
					# VASP does NOT do normalization in elf.F
					phi_q = self.readBandCoeff(ispin=ispin+1, ikpt=ikpt+1,
											   iband=iband+1,
											   norm=False)
					# pad the missing planewave coefficients for -G vectors
					if self._lgam:
						tmp = [x.conj() for x in phi_q[1:]]
						phi_q = np.concatenate([phi_q, tmp])
						# Gamma only, divide a factor of sqrt(2.0) except for
						# G=0
						phi_q /= np.sqrt(2.0)
						phi_q[0] *= np.sqrt(2.0)
					# wavefunction in real space
					phi_r = self.wfc_r(ispin=ispin+1, ikpt=ikpt+1,
									   iband=iband+1,
									   ngrid=ngrid,
									   gvec=igvec,
									   Cg=phi_q) * normFac
					# grad^2 \phi in reciprocal space
					lap_phi_q = -gk2 * phi_q
					# grad^2 \phi in real space
					lap_phi_r = self.wfc_r(ispin=ispin+1, ikpt=ikpt+1,
										   iband=iband+1,
										   ngrid=ngrid,
										   gvec=igvec,
										   Cg=lap_phi_q) * normFac

					# \phi* grad^2 \phi in real space --> kinetic energy density
					tau += -phi_r * lap_phi_r.conj() * weight

					# charge density in real space
					rho += phi_r.conj() * phi_r * weight

			# charge density in reciprocal space
			rho_q = np.fft.fftn(rho, norm='ortho')

			# grad^2 rho: laplacian of charge density
			lap_rho_q = -G2 * rho_q
			lap_rho_r = np.fft.ifftn(lap_rho_q, norm='ortho')

			# charge density gradient: grad rho
			########################################
			# wrong method for gradient using FFT
			########################################
			# grad_rho_x = np.fft.ifft(1j * Gx * np.fft.fft(rho, axis=0), axis=0)
			# grad_rho_y = np.fft.ifft(1j * Gy * np.fft.fft(rho, axis=1), axis=1)
			# grad_rho_z = np.fft.ifft(1j * Gz * np.fft.fft(rho, axis=2), axis=2)

			########################################
			# correct method for gradient using FFT
			########################################
			grad_rho_x = np.fft.ifftn(1j * Gx * rho_q, norm='ortho')
			grad_rho_y = np.fft.ifftn(1j * Gy * rho_q, norm='ortho')
			grad_rho_z = np.fft.ifftn(1j * Gz * rho_q, norm='ortho')

			grad_rho_sq = np.abs(grad_rho_x)**2 \
				+ np.abs(grad_rho_y)**2 \
				+ np.abs(grad_rho_z)**2

			rho = rho.real
			tau = tau.real
			lap_rho_r = lap_rho_r.real

			Cf = 3./5 * (3.0 * np.pi**2)**(2./3)
			Dh = np.where(rho > 0.0,
						  Cf * rho**(5./3),
						  0.0)
			eps = 1E-8 / HSQDTM
			Dh[Dh < eps] = eps
			# D0 = T + TCORR - TBOS
			D0 = tau + 0.5 * lap_rho_r - 0.25 * grad_rho_sq / rho

			ElectronLocalizationFunction.append(1. / (1. + (D0 / Dh)**2))

		return ElectronLocalizationFunction

	def write_std_wavecar(self, out="WAVECAR_std"):
		assert self._lgam

		with open(out, "w") as new_wc:
			new_nplws = 2*(self._nplws[0] - 1) + 1
			plws_rec_size = np.max(new_nplws)*np.dtype(self.setWFPrec()).itemsize
			band_rec_size = np.dtype(np.float64).itemsize*(self._nbands*3+1)
			# record needs to be large enough to contain both plane waves and bands
			new_rec_size = max(plws_rec_size, band_rec_size)
			nfloat = new_rec_size//8 # number of float64s per record
			# header line
			rec = np.zeros(nfloat, dtype=np.float64)
			rec[0:3] = new_rec_size,self._nspin,self._rtag
			rec.tofile(new_wc)
			# header line 2
			rec[0:3] = self._nkpts,self._nbands,self._encut
			rec[3:3+9] = self._Acell.reshape((1, -1))
			rec.tofile(new_wc)
			wave_rec = np.zeros(new_nplws, dtype=self.setWFPrec())
			for spin in range(self._nspin):
				rec[0] = new_nplws
				rec[1:1+3] = self._kvecs[0]
				rec[4 : 4+3*self._nbands : 3] = self._bands[spin, 0, :]
				rec[4+1 : 4+3*self._nbands : 3] = 0.0 # so far energies always real?
				rec[4+2 : 4+3*self._nbands : 3] = self._occs[spin, 0, :]
				rec.tofile(new_wc)

				### Expand plane wave coefficients
				ngrid = self._ngrid.copy() * 2
				if self._gam_half == "z":
					ordered_grid = ngrid[[2,1,0]]
				elif self._gam_half == "x":
					ordered_grid = ngrid
				else:
					raise ValueError("Gamma reduction direction must be z or x")

				Gvec = self.gvectors(ikpt=1, check_consistency=True)
				full_Gvec = self.gvectors(ikpt=1, force_Gamma=False, check_consistency=False)
				phi_k = np.zeros(ngrid, dtype=wave_rec.dtype)
				for band in range(self._nbands):
					phi_k[Gvec[:, 0], Gvec[:,1], Gvec[:,2]] = self.readBandCoeff(spin+1,1,band+1)
					## Upper sphere
					for ii in range(1, ordered_grid[0] // 2 + 1):
						for jj in range(ordered_grid[1]):
							for kk in range(ordered_grid[2]):
								phi_k[-ii,-jj,-kk] = phi_k[ii,jj,kk].conjugate()

					## Upper part of x-y-plane
					for jj in range(1, ordered_grid[1] // 2 + 1):
						for kk in range(ordered_grid[2]):
							phi_k[0,-jj,-kk] = phi_k[0,jj,kk].conjugate()

					## Upper part of x-axis
					for kk in range(1, ordered_grid[2] // 2 + 1):
						phi_k[0,0,-kk] = phi_k[0,0,kk].conjugate()

					phi_k /= np.sqrt(2.)
					phi_k[0,0,0] *= np.sqrt(2.)

					wave_rec[:new_nplws] = phi_k[full_Gvec[:, 0], full_Gvec[:, 1], full_Gvec[:, 2]]
					wave_rec.tofile(new_wc)

############################################################


if __name__ == '__main__':
	# xx = vaspwfc('wavecar')
	# phi = xx.wfc_r(1, 30, 17, ngrid=(28, 28, 252))
	# xx.save2vesta(phi, poscar='POSCAR')

	# xx = vaspwfc('./gamma/WAVECAR')
	# phi = xx.wfc_r(1, 1, 317, ngrid=(60, 108, 160),
	#				 gamma=True)
	# xx.save2vesta(phi, poscar='./gamma/POSCAR',gamma=True)

	# xx = vaspwfc('WAVECAR')
	# dE, ovlap, tdm = xx.TransitionDipoleMoment([1,30,17], [1,30,18], norm=True)
	# print dE, ovlap.real, np.abs(tdm)**2

	# print xx._recl, xx._nspin, xx._rtag
	# print xx._nkpts, xx._nbands, xx._encut
	# print xx._Acell, xx._Bcell
	# # print np.linalg.norm(xx._Acell, axis=1)
	# print xx._ngrid
	# print xx._bands[0,0,:]
	# print xx._kvecs
	# print xx._kpath
	# b = xx.readBandCoeff(1,1,1)
	# xx = np.savetxt('kaka.dat', xx.gvectors(2), fmt='%5d')
	# gvec = xx.gvectors(1)
	# gvec %= xx._ngrid[np.newaxis, :]
	# print gvec

	# ngrid=(28, 28, 252)
	# phi = xx.wfc_r(1, 30, 17, ngrid=(28, 28, 252))
	# header = open('POSCAR').read()
	# with open('wave_real.vasp', 'w') as out:
	#	  out.write(header)
	#	  out.write('%5d%5d%5d\n' % (ngrid[0], ngrid[1], ngrid[2]))
	#	  nwrite=0
	#	  for kk in range(ngrid[2]):
	#		  for jj in range(ngrid[1]):
	#			  for ii in range(ngrid[0]):
	#				  nwrite += 1
	#				  out.write('%22.16f ' % phi.real[ii,jj,kk])
	#				  if nwrite % 10 == 0:
	#					  out.write('\n')
	# with open('wave_imag.vasp', 'w') as out:
	#	  out.write(header)
	#	  out.write('%5d%5d%5d\n' % (ngrid[0], ngrid[1], ngrid[2]))
	#	  nwrite=0
	#	  for kk in range(ngrid[2]):
	#		  for jj in range(ngrid[1]):
	#			  for ii in range(ngrid[0]):
	#				  nwrite += 1
	#				  out.write('%22.16f ' % phi.imag[ii,jj,kk])
	#				  if nwrite % 10 == 0:
	#					  out.write('\n')

	# xx = vaspwfc('wave_tyz')
	# ipr = xx.inverse_participation_ratio()
	# print xx._nbands, xx._nkpts
	#
	# import matplotlib as mpl
	# import matplotlib.pyplot as plt
	#
	# fig = plt.figure()
	# ax = plt.subplot()
	#
	# ax.scatter(ipr[...,0], ipr[..., 1], s=ipr[..., 2] / ipr[..., 2].max() * 10, c=ipr[..., 2],
	#			 cmap='jet_r')
	#
	# plt.show()

	wfc = vaspwfc('WAVECAR', lgamma=True, gamma_half='x')
	# ngrid = [80, 140, 210]
	phi = wfc.wfc_r(iband=190)

	rho = np.abs(phi)**2
	# rho2 = VaspChargeDensity('PARCHG.0158.ALLK').chg[0]
	# rho /= rho.sum()
	# rho2 /= rho2.sum()
	# rho3 = rho - rho2

	wfc.save2vesta(rho, lreal=True)

	pass
