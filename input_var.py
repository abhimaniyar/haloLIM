from headers_constants import *


class data_var_iv(object):

    def __init__(self, linedir):  # , mass):  # , z):  # , ell):
        # ############### cib data #########################
        self.line = linedir['name']

        self.z_c = 1.5

        if self.line == '[CII]':
            self.nu0 = 1900.  # GHz rest freq
        elif self.line == 'CO10':
            self.nu0 = 115.271

        # ######### CIB halo model parameters ###################
        cibparresaddr = linedir['cibpar_resfile']
        self.Meffmax, self.etamax, self.sigmaMh, self.tau = np.loadtxt(cibparresaddr)[:4, 0]
        self.shotpl = np.loadtxt(cibparresaddr)[4:8, 0]  # 217, 353, 545, 857
        # self.Meffmax, self.etamax, self.sigmaMh, self.tau = 8753289339381.791, 0.4028353504978569, 1.807080723258688, 1.2040244128818796

        # if name == 'Planck':
            # self.fc[-4:] = np.loadtxt(cibparresaddr)[-4:, 0]
