# cosmo relared inputs

from headers_constants import *


class cosmo_var(object):

    def __init__(self, mass, z):  # , ell):
        self.mass = mass
        self.z = z
        # self.ell = ell

        self.deltah = 200.
        # ########## reading in the matter power spectra #############
        redshifts = np.loadtxt('data_files/redshifts.txt')
        nr = len(redshifts)
        self.zpk = redshifts

        if min(self.z) < min(redshifts) or max(self.z) > max(redshifts):
            print ("If the redshift range is outside of [%s to %s], then " +
                   "values of the matter power spectrum " +
                   "is extrapolated and might be incorrect.") % (min(redshifts), max(redshifts))
        ll = [str(x) for x in range(1, 211)]
        addr = 'data_files/matter_power_spectra'
        pkarray = np.loadtxt('%s/test_highk_lin_matterpower_210.dat' % (addr))
        self.k = pkarray[:, 0]*cosmo.h
        self.Pk = np.zeros((len(self.k), len(redshifts)))
        for i in range(len(redshifts)):
            pkarray = np.loadtxt("%s/test_highk_lin_matterpower_%s.dat" % (addr, ll[209-i]))
            self.Pk[:, i] = pkarray[:, 1]/cosmo.h**3

        self.pkinterpz = interp1d(redshifts, self.Pk, kind='linear', bounds_error=False, fill_value="extrapolate")
        # self.pkinterpk = interp1d(k, Pk.T, kind='linear', bounds_error=False, fill_value="extrapolate")

    # ######## hmf, bias, nfw ###########
    def hmf(self):
        print ("Calculating the halo mass function " +
               "for given mass and redshift for CIB mean calculations.")
        nm = len(self.mass)
        nz = len(self.z)
        hmfmz = np.zeros((nm, nz))
        delta_h = self.deltah

        for r in range(nz):
            pkz = self.pkinterpz(self.z[r])
            instance = hmf_unfw_bias.h_u_b(self.k, pkz, self.z[r],
                                           cosmo, delta_h, self.mass)
            hmfmz[:, r] = instance.dn_dm()
        """
        for r in range(nr):
            instance = hmf_unfw_bias.h_u_b(k, Pk[:, r], redshifts[r],
                                           cosmo, delta_h, self.mass)
            hmfr[:, r] = instance.dn_dlogm()
        self.hmf_r = interp1d(redshifts, hmfr, kind='linear',
                              bounds_error=False, fill_value=0.)
        """
        return hmfmz

    def nfw_u(self):
        """
        nfwaddr = addr+'/cib_hod/data/'
        file_nfw = 'u_nfw_precalculated.fits'
        hdulist = fits.open("%s" % (nfwaddr+file_nfw))
        u_nfw = hdulist[0].data
        hdulist.close()
        """
        nm = len(self.mass)
        nz = len(self.z)
        nfwu = np.zeros((nm, len(self.k), nz))
        delta_h = self.deltah
        for r in range(nz):
            pkz = self.pkinterpz(self.z[r])
            instance = hmf_unfw_bias.h_u_b(self.k, pkz,
                                           self.z[r],
                                           cosmo, delta_h, self.mass)
            nfwu[:, :, r] = instance.nfwfourier_u()
        return nfwu

    def bias_m_z(self):
        nm = len(self.mass)
        nz = len(self.z)
        biasmz = np.zeros((nm, nz))
        delta_h = self.deltah
        for r in range(nz):
            pkz = self.pkinterpz(self.z[r])
            instance = hmf_unfw_bias.h_u_b(self.k, pkz,
                                           self.z[r],
                                           cosmo, delta_h, self.mass)
            biasmz[:, r] = instance.b_nu()
        return biasmz

    def dchi_dz(self, z):
        a = c_light/(cosmo.H0*np.sqrt(cosmo.Om0*(1.+z)**3 + cosmo.Ode0))
        return a.value

    def chi(self, z):
        return cosmo.comoving_distance(z).value

    def karray(self, ell, z):
        nl = len(ell)
        nz = len(z)
        k_array = np.zeros((nl, nz))

        for i in range(nl):
            k_array[i, :] = ell[i]/self.chi(z)

        return k_array

    def interp_bias(self, z):
        nm, nz = len(self.mass), len(z)
        bias = np.zeros((nm, nz))
        for m in range(nm):
            bias[m, :] = np.interp(z, self.zpk, self.bias_m_z[m, :])
        return bias

    def interp_nfw(self, ell, z):
        nm, nl, nz = len(self.mass), len(ell), len(z)
        nfw_ureds = np.zeros((nm, len(self.k), nz))
        for i in range(len(self.k)):
            for m in range(nm):
                nfw_ureds[m, i, :] = np.interp(z, self.zpk, self.nfw_u[m, i, :])

        u_nfw = np.zeros((nm, nl, nz))
        k_array = self.karray(ell, z)
        for m in range(nm):
            for j in range(nz):
                u_nfw[m, :, j] = np.interp(k_array[:, j], self.k,
                                           nfw_ureds[m, :, j])
        return u_nfw

    def Pk_array(self, ell, z):
        nl = len(ell)
        nz = len(z)
        nreds = len(self.zpk)
        pk1 = np.zeros((nl, nreds))
        Pk_int = np.zeros((nl, nz))

        for i in range(nreds):
            ell_chi = ell/self.chi(self.zpk[i])
            pk1[:, i] = np.interp(ell_chi, self.k, self.Pk[:, i])

        for i in range(nl):
            Pk_int[i, :] = np.interp(z, self.zpk, pk1[i, :])

        return Pk_int
