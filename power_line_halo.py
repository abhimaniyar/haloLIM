from headers_constants import *


class power_halo:

    def __init__(self, data_var_iv, cosmo_var_iv):
        self.dv = data_var_iv
        self.uni = cosmo_var_iv
        self.z = self.uni.z  # self.dv.z
        self.z_c = self.dv.z_c
        self.mh = self.uni.mass

        self.line = self.dv.line
        self.nu0 = self.dv.nu0

        # self.deltah = deltah
        self.Meffmax = self.dv.Meffmax
        self.etamax = self.dv.etamax
        self.sigmaMh = self.dv.sigmaMh
        self.tau = self.dv.tau
        self.dndm = self.uni.hmf()  # self.dv.hmf
        # self.biasmz = self.uni.bias_m_z
        # self.biasmz = self.uni.interp_bias(self.z)
        self.biasmz = self.uni.bias_m_z()
        self.sig_z = np.array([max(self.z_c - r, 0.) for r in self.z])
        self.sigpow = self.sigmaMh - self.tau*self.sig_z
        # self.unfw = self.uni.interp_nfw(self.ell, self.z)

    def sfr_mhdot(self, mhalo):
        """ SFR/Mhdot lognormal distribution wrt halomass """
        if hasattr(mhalo, "__len__"):
            a = np.zeros((len(mhalo), len(self.z)))
            for i in range(len(mhalo)):
                if mhalo[i] < self.Meffmax:
                    a[i, :] = self.etamax * np.exp(-(np.log(mhalo[i]) - np.log(self.Meffmax))**2 / (2 * self.sigmaMh**2))
                else:
                    a[i, :] = self.etamax * np.exp(-(np.log(mhalo[i]) - np.log(self.Meffmax))**2 / (2 * self.sigpow**2))
        else:
            if mhalo < self.Meffmax:
                a = self.etamax * np.exp(-(log(mhalo) - log(self.Meffmax))**2 / (2 * self.sigmaMh**2))
            else:
                a = self.etamax * np.exp(-(log(mhalo) - log(self.Meffmax))**2 / (2 * self.sigpow**2))
        return a

    def Mdot(self, mhalo):
        use_mean = True
        if use_mean:
            a = 46.1*(1 + 1.11*self.z) * \
                np.sqrt(cosmo.Om0 * (1 + self.z)**3 + cosmo.Ode0)
            b = (mhalo / 1.0e12)**1.1
            return np.outer(b, a)
        else:
            a = 25.3*(1 + 1.65*self.z) * \
                np.sqrt(cosmo.Om0*(1 + self.z)**3 + cosmo.Ode0)
            b = (mhalo / 1.0e12)**1.1
            return np.outer(b, a)

    def sfr(self, mhalo):
        sfrmhdot = self.sfr_mhdot(mhalo)
        mhdot = self.Mdot(mhalo)
        f_b = cosmo.Ob(self.z)/cosmo.Om(self.z)
        return mhdot * f_b * sfrmhdot

    def subhmf(self, mhalo, ms):
        # subhalo mass function from (https://arxiv.org/pdf/0909.1325.pdf)
        return 0.13*(ms/mhalo)**(-0.7)*np.exp(-9.9*(ms/mhalo)**2.5)*np.log(10)
    # np.log(10) added in the end as in the integration we are integrating with
    # respect to dlogm to the base 10.

    def msub(self, mhalo):
        """
        for a given halo mass mh, the subhalo masses would range from
        m_min to mh. For now, m_min has been taken as 10^5 solar masses
        """
        log10msub_min = 5
        if np.log10(mhalo) <= log10msub_min:
            raise ValueError, "halo mass %d should be greater than subhalo mass \
%d." % (np.log10(mhalo), log10msub_min)
        else:
            logmh = np.log10(mhalo)
            logmsub = np.arange(log10msub_min, logmh, 0.1)
            return 10**logmsub

    def L_line(self, model):
        if self.line == '[CII]':
            if model == 'Silva_15':
                # https://iopscience.iop.org/article/10.1088/0004-637X/806/2/209/pdf
                # m2 parameterization from the paper
                a_LCII = 1.0000
                b_LCII = 6.9647
                res = a_LCII*np.log10(self.sfr(self.mh)) + b_LCII
                result = 10**res
                return result
            elif model == 'Chung_20':
                # https://arxiv.org/pdf/1812.08135.pdf
                alpha = 1.40 - 0.07*self.z
                beta = 7.1 - 0.07*self.z
                logL = alpha*np.log10(self.sfr(self.mh)) + beta
                # result = np.exp(logL)
                result = 10**logL
                return result
            elif model == 'Schaerer_20':
                # a and b are average values taken from Tab A.1 of
                # https://arxiv.org/pdf/2002.00979.pdf
                a = 6.90
                b = 1.02
                logL = a+b*np.log10(self.sfr(self.mh))
                result = 10**logL
                return result
        elif self.line == 'CO10':
            if model == 'Li_16':
                KC = 1.0e-10
                """
                Kennicutt constant which will be multiplied by
                # d_MF factor below to get the constant corresponding to
                # Chabrier
                # or
                # some other IMF
                # dMF, alpha and beta are average values taken from values
                # mentioned in blue color in Fig. 7 of
                # arXiv:1503.08833v2
                # we are going to calculate the LCO different for central and
                # satellite halos based on their respective SFRs and hmfs
                """
                log_dMF = 0.0
                d_MF = 10**log_dMF
                alpha = 1.24
                beta = -0.48
                LIR_c = self.sfr(self.mh)/KC/d_MF
                logLco_c = (np.log10(LIR_c)-beta)/alpha
                Lco_c_prime = 10**logLco_c
                """
                # CO luminosities obtained above are in units of K Km/sec pc^2
                # Need to multiply by following factors to get it in L_sun
                """
                # nu_co_rest = 115.27  # GHz
                nu_co_rest = self.nu0  # GHz
                result = 4.9e-5*(nu_co_rest/115.27)**3*Lco_c_prime
                return result
            elif model == 'Lidz_11':
                # https://arxiv.org/pdf/1104.4800.pdf
                result = 3.2e4*(self.sfr(self.mh))**(3./5)
                return result
            elif model == 'Righi_08':
                # https://arxiv.org/pdf/1104.4800.pdf
                result = 3.7e3*self.sfr(self.mh)
                return result

    def rho_L(self, mod):
        # note here that hmf is dn_dm and not dn_dlogm as we normally use
        L_c = self.L_line(mod)
        # print (self.dndm)
        integrand = self.dndm*L_c
        rhoL_c = intg.simps(integrand, x=self.mh, axis=0, even='avg')
        return rhoL_c

    def I_line(self, mod):
        # final units here are L_sol/Mpc/Mpc/GHz, if we need units of Jy/sr
        # we need this convfac => 1 L_sun/Mpc**2/GHz = 4.0204e-2 Jy/sr
        convfac = 4.0204e-2  # Jy/sr per Lsol/Mpc/Mpc/GHz
        H = cosmo.H(self.z).value
        a = 4*np.pi*self.nu0*H
        b = c_light/a
        rhoL_c = self.rho_L(mod)
        I_c = rhoL_c*b*convfac
        return I_c

    def beff(self, mod):
        hmf = self.dndm
        bmz = self.biasmz
        L = self.L_line(mod)
        bmzhmf = bmz*hmf
        num_int = L*bmzhmf
        den_int = L*hmf
        numerator = intg.simps(num_int, x=self.mh, axis=0, even='avg')
        # den_int = Lcii*hmf
        denominator = intg.simps(den_int, x=self.mh, axis=0, even='avg')
        return numerator/denominator

    def Pclust(self, mod):
        """
        units for I => Jy/sr
        Thus units for Pclust = (Jy/sr)^2 MPc^3
        """
        Plin = self.uni.pkinterpz(self.z)
        I2 = (self.I_line(mod))**2
        # print (np.shape(I2))
        beff2 = (self.beff(mod))**2
        # print (np.shape(beff2))
        res = I2*beff2*Plin
        """
        res = res[:, 0]
        result = interp1d(self.uni.k, res, kind='linear', bounds_error=False,
                          fill_value=0.)
        """
        result = res
        return result

    def Pshot(self, mod):
        convfac = 4.0204e-2  # Jy/sr per Lsol/Mpc/Mpc/GHz
        """
        final units are Lsol**2/Mpc/GHz/GHz. We want (Jy/sr)^2 Mpc^3
        Lsol**2/Mpc/GHz/GHz = (Lsol/Mpc/Mpc/GHz)^2 Mpc^3
        Therefore we multiply the final answer by convfac^2 to get the units
        of (Jy/sr)^2 Mpc^3
        """
        H = cosmo.H(self.z).value
        a = 4*np.pi*self.nu0*H
        b = c_light/a
        # geo = (a/b)**2
        L = self.L_line(mod)
        hmf = self.dndm
        integrand = hmf*L**2
        result1 = intg.simps(integrand, x=self.mh, axis=0, even='avg')
        res = result1*b**2*convfac**2
        result = res
        # result = np.repeat(res, len(self.k))
        return result

    """
    def P_1h(self):
        I_c, I_sub = self.I_line()
        fc = I_c
        fsat = I_sub[:, None]*unfw
        integrand = 2*fc[:, None]*fsat
        integrand += fsat**2
        integrand /= self.cosmology.dn_dm[:, None]
        result = integrate.simps(integrand, x=self.Mh, axis=0, even='avg')
        return result

    def D_2h(self, unfw):
        I_c, I_sub = self.I_line()
        fc = I_c
        fsat = I_sub[:, None]*unfw
        bmz = self.cosmology.b_nu()
        integrand = (fc[:, None] + fsat)*bmz[:, None]
        result = integrate.simps(integrand, x=self.Mh, axis=0, even='avg')
        return result

    def P_2h(self, unfw):
        D2h = self.D_2h(unfw)
        Plin = self.cosmology.Pk
        a = D2h**2
        result = a*Plin
        return result

    def Pshot(self):
        convfac = 4.0204e-2  # Jy/sr per Lsol/Mpc/Mpc/GHz
        # final units are Lsol**2/Mpc/GHz/GHz. We want (Jy/sr)^2 Mpc^3
        # Lsol**2/Mpc/GHz/GHz = (Lsol/Mpc/Mpc/GHz)^2 Mpc^3
        # Therefore we multiply the final answer by convfac^2 to get the units
        # of (Jy/sr)^2 Mpc^3
        H = self.cosmology.cosmo.H(self.z).value
        a = 4*np.pi*self.nu_rest*H
        b = c_light/a
        # geo = (a/b)**2
        L_c, L_sub = self.L.L_line(self.mod)
        # print np.shape(Lcii_c), np.shape(Lcii_sub)
        Ltot = L_c + L_sub
        integrand = self.cosmology.dn_dm*Ltot**2
        result1 = integrate.simps(integrand, x=self.Mh, axis=0, even='avg')
        result = result1*convfac**2*b**2
        return result
    """
