#%% Cell: Define the imports and needed functions

import numpy as np
import matplotlib.pyplot as plt
#import colorcet as cc

from func_fill_time_gaps_with_Nan_profiles import func_fill_time_gaps_with_Nan_profiles

#%% Cell: Start the function

def func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                 ylabels, plot_xlabel_flag, xlabels, input_top_data,
                 input_mid_data, input_bot_data, input_time,
                 top_color_values, mid_color_values, bot_color_values,
                 top_color_ticks, mid_color_ticks, bot_color_ticks,
                 title_top_str, title_mid_str, title_bot_str,
                 plot_filename):

    #%% Cell: Define the initial values

    # Define some plot attributes
    max_range = int(np.ceil(np.amax(plot_range/1000)))

        
    f = (input_time >= start_hour) & (input_time <= end_hour)

    if(np.sum(f) > 3):
        plot_top = input_top_data[f, :]
        plot_mid = input_mid_data[f, :]
        plot_bot = input_bot_data[f, :]
        plot_time = input_time[f]

        plot_time, plot_top, plot_mid, plot_bot = func_fill_time_gaps_with_Nan_profiles(
                plot_time, plot_top, plot_mid, plot_bot)

        ## Top Panel ##
        #cmap = plt.get_cmap('viridis')
        #  'viridis', 'plasma', 'inferno', 'magma', 'cividis'
        #cmap = plt.get_cmap('inferno')
        cmap = plt.get_cmap('viridis')
        
        fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)

        zdb = ax0.pcolormesh(plot_time, (plot_range - plot_drange/2) /
                             1000, np.transpose(plot_top[:-1, :-1]), shading='flat', 
                             vmin=top_color_values[0], vmax=top_color_values[2], cmap=cmap)

        cbar = fig.colorbar(zdb, ax=ax0, ticks=np.arange(top_color_values[0], 
                            top_color_values[2]+top_color_values[1], top_color_values[1]), aspect=4)
        cbar.ax.set_yticklabels(top_color_ticks)
        ax0.axis([start_hour, end_hour, 0, max_range])
        if(plot_xlabel_flag):
            x = np.arange(len(xlabels))
            ax0.set_xticks(x)
            ax0.set_xticklabels(xlabels, fontsize=8)
        # end if(plot_xlabel_flag):
        y = np.arange(0, max_range+1, 1)
        ax0.set_yticks(y)
        ax0.set_yticklabels(ylabels[0:max_range+1], fontsize=8)
        ax0.set_title(title_top_str, fontsize=10)
            # f'SGP {site_id}, RWP, Short Pulse (%im), SNR [dB], {date_plot_str}' %plot_drange, fontsize=10)
        ax0.set_axisbelow(True)
        ax0.grid(True)
        ax0.set_ylabel('Height [km]', fontsize=8)
        ax0.text(.05, .8, '(a)', fontsize=8, horizontalalignment='center',
                 verticalalignment='center', transform=ax0.transAxes, 
                 bbox=dict(facecolor='white', edgecolor='black'))

        ## Middle Panel ##
        Vmean = ax1.pcolormesh(plot_time, (plot_range - plot_drange/2) /
                           1000, np.transpose(plot_mid[:-1, :-1]), shading='flat', 
                           vmin=mid_color_values[0], vmax=mid_color_values[2], cmap=cmap)

        cbar = fig.colorbar(Vmean, ax=ax1, ticks=np.arange(mid_color_values[0], 
                            mid_color_values[2]+mid_color_values[1], mid_color_values[1]), aspect=4)
        cbar.ax.set_yticklabels(mid_color_ticks)
        ax1.axis([start_hour, end_hour, 0, max_range])
        if(plot_xlabel_flag):
            ax1.set_xticks(x)
            ax1.set_xticklabels(xlabels, fontsize=8)
        # end if(plot_xlabel_flag):
        ax1.set_yticks(y)
        ax1.set_yticklabels(ylabels[0:max_range+1], fontsize=8)
        ax1.set_title(title_mid_str, fontsize=10)
            #f'SGP {site_id}, RWP, Short Pulse, Mean Velocity [m/s] (+ downward)', fontsize=10)
        ax1.grid(True)
        ax1.set_axisbelow(True)
        ax1.set_ylabel('Height [km]', fontsize=8)
        ax1.text(.05, .8, '(b)', fontsize=8, horizontalalignment='center',
                 verticalalignment='center', transform=ax1.transAxes, 
                 bbox=dict(facecolor='white', edgecolor='black'))

        ## Bottom Panel ##
        Vsig = ax2.pcolormesh(plot_time, (plot_range - plot_drange/2) /
                          1000, np.transpose(plot_bot[:-1, :-1]), shading='flat', 
                          vmin=bot_color_values[0], vmax=bot_color_values[2], cmap=cmap)

        cbar = fig.colorbar(Vsig, ax=ax2, ticks=np.arange(bot_color_values[0], 
                            bot_color_values[2]+bot_color_values[1], bot_color_values[1]), aspect=4)
        cbar.ax.set_yticklabels(bot_color_ticks)
        ax2.axis([start_hour, end_hour, 0, max_range])
        if(plot_xlabel_flag):
            ax2.set_xticks(x)
            ax2.set_xticklabels(xlabels, fontsize=8)
        # end if(plot_xlabel_flag):
        ax2.set_yticks(y)
        ax2.set_yticklabels(ylabels[0:max_range+1], fontsize=8)
        ax2.set_title(title_bot_str, fontsize=10)
                #f'SGP {site_id}, RWP, Short Pulse, Spectrum Breadth (Vsig) [m/s]', fontsize=10)
        ax2.grid(True)
        ax2.set_axisbelow(True)
        ax2.set_ylabel('Height [km]', fontsize=8)
        ax2.set_xlabel('Hour of Day [UTC]', fontsize=8)
        ax2.text(.05, .8, '(c)', fontsize=8, horizontalalignment='center',
                 verticalalignment='center', transform=ax2.transAxes, 
                 bbox=dict(facecolor='white', edgecolor='black'))

        fig.tight_layout()

        # Save destination here
        plt.savefig(plot_filename, dpi=800)

        # close all of the figures
        plt.close('all')

    # end if(np.sum(f) > 3):

    return
