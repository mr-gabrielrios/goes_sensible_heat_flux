from matplotlib import rcParams
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime as dt

def agg_plot(data, loc, **kwargs):

    font_sz = 10
    rcParams['font.family'] = 'FreeSans'
    rcParams['font.size'] = font_sz
    
    ax_num = 6
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(nrows=6, ncols=1, figsize=(6, ax_num*2))    
    axs = (ax1, ax2, ax3, ax4, ax5, ax6)
    
    for ax in axs:
        ax.grid(axis='x', linestyle=':', which='both')
        ax.set_xlim([data['date'].iloc[0], data['date'].iloc[-1]])
    
    data['date'] = [t.to_pydatetime() for t in data['date']]
    
    title_string = loc + " Data, " + str(data['date'].iloc[0]) + " to " + str(data['date'].iloc[-1])
    fig.subplots_adjust(top=0.8)
    fig.tight_layout(pad=0.75)
    fig.suptitle(title_string, fontsize=font_sz*1.5, y=1.05)
    
    ax1.set_title('Sensible heat flux')
    ax1.set_ylabel('Q_H [W/m^2]')
    ax1.plot(data['date'], data['q_h'])
    ax1.plot(data['date'], data['q_h_obs'])
    ax1.legend(["Model", "Observed"], loc='best')
    
    ax2.set_title('Model/observation error')
    ax2.plot(data['date'], data['q_h_error'])
    ax2.set_ylabel('%')
    ax2.set_ylim([-100, 100])
    
    ax3.set_title(r'Atmospheric stability ($\zeta$)')
    ax3.plot(data['date'], data['zeta'])
    ax3.set_ylabel('z/L')
    
    ax4.set_title('Temperatures')
    ax4.set_ylabel('T [K]')
    ax4.plot(data['date'], data['T_lst'])
    ax4.plot(data['date'], data['T_air'])
    ax4.plot(data['date'], data['T_air_obs'])
    ax4.legend(['T_lst', 'T_air', 'T_air_obs'], loc='best')
    
    ax5.set_title('T_lst - T_air')    
    ax5.set_ylabel('T [K]')
    ax5.plot(data['date'], data['T_lst']-data['T_air'])
    ax5.plot(data['date'], data['T_lst']-data['T_air_obs'])
    ax5.legend(['T_lst-T_air', 'T_lst-T_air_obs'], loc='best')
    
    ax6.set_title('Wind velocities')
    # ax6.plot(data['date'], data['u_r'])
    ax6.plot(data['date'], data['u_star'])
    ax6.plot(data['date'], data['u_star_obs'])
    ax6.set_ylabel('u* [m/s]')
    ax7 = ax6.twinx()  # instantiate a second axes that shares the same x-axis
    ax7.plot(data['date'], data['u_r'], color='k', linestyle=':')
    ax7.set_ylabel('u_r [m/s]')
    ax6.legend(['Model u*', 'Observed u*'], loc='0')
    
    plt.setp(ax1, xticklabels=[])
    plt.setp(ax2, xticklabels=[])
    plt.setp(ax3, xticklabels=[])
    plt.setp(ax4, xticklabels=[])
    plt.setp(ax5, xticklabels=[])
    plt.setp(ax6.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.show()
    
    # Save figures with the following format
    # LOCATION_sYYYYmmddHHMM_eYYYYmmddHHMM_rYYYYmmddHHMM_T-AIR-SOURCE
    # 's' indicates model scope start time, 'e' indicates model scope end time, 'r' indicates model run time
    if kwargs['savefig']:
        fn = 'plots/' + loc + '_s' + str(data['date'].iloc[0].strftime('%Y%m%d%H%M')) + '_e' +  str(data['date'].iloc[-1].strftime('%Y%m%d%H%M')) + '_r' + str(dt.now().strftime('%Y%m%d%H%M'))
        if kwargs['t_air_src']:
            fn = fn + '_t-air-' + kwargs['t_air_src']
        fig.savefig(fn, bbox_inches='tight')