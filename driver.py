from headers_constants import *
import power_line_halo
import cosmo_related
import input_var

logmass = np.arange(6, 15.005, 0.1)
mass = 10**logmass1
# logmass1 = np.arange(10, 14.705, 0.1)  # SAM
# mass1 = 10**logmass1
z = np.array([3.0, 4.0])

line_name = '[CII]'  # [CII] or CO10

model = 'Schaerer_20'
"""
Choose model name for a given line from following list
# cii_modelnames: Silva_15, Chung_20, Schaerer_20
# co_modelname: Li_16, Lidz_11, Righi_08
"""

# location for the best fit value of the lognormal parameterization for SFR
strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nolens_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"
# strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nospire_fcpl_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"

cibres = "data_files/one_halo_bestfit_"+strfig
linedir = {'name': line_name,
           'cibpar_resfile': cibres}


def plot_Pk(line, mod):
    datavar = input_var.data_var_iv(line)
    cosmovar = cosmo_related.cosmo_var(mass, z)
    powerhalo = power_line_halo.power_halo(datavar, cosmovar)
    # P1h = line_powerhalo.P_1h(u_nfw_int)
    k_plot = np.linspace(1e-2, 3, 100)
    # print (powerhalo.Pshot(mod))
    # print (np.shape(powerhalo.Pshot(mod)))
    Pclust = interp1d(cosmovar.k, powerhalo.Pclust(mod)[:, 0], kind='linear',
                      bounds_error=False, fill_value=0.)
    Pshot = powerhalo.Pshot(mod)[0]

    # P_clust = powerhalo.Pclust(mod)(k_plot)
    # P_shot = powerhalo.Pshot(mod)(k_plot)
    P_clust = Pclust(k_plot)
    P_shot = np.repeat(Pshot, len(k_plot))
    P_tot = P_clust + P_shot

    fig = plt.figure(figsize=(11.5, 7))
    ax = fig.add_subplot(111)

    ax.plot(k_plot, k_plot**3*P_clust/(2*np.pi**2), label='Clustering', ls='--')
    ax.plot(k_plot, k_plot**3*P_shot/(2*np.pi**2), label='Shot', ls='-.')
    ax.plot(k_plot, k_plot**3*P_tot/(2*np.pi**2), label='Total', ls='-')

    ax.legend(fontsize='18')  # , bbox_to_anchor=(0.38, 0.45))  # , labelspacing=0.1)
    ax.set_xscale('log')
    ax.set_yscale('log', nonposy='mask')
    ax.set_xlabel(r'$k [\frac{1}{\rm Mpc}]$', fontsize=24)
    ax.set_ylabel(r'$\Delta_k^2 [{\rm Jy}^2/{\rm sr}]$', fontsize=26)
    # ax.set_ylim((1e-10))  # , 4.e-6))
    # ax.set_xlim((5., 2.e3))
    ax.tick_params(axis='both', labelsize=20)
    plt.show()


plot_Pk(linedir, model)
