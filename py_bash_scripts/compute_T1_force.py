from turtle import pos
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os
from scipy.ndimage import gaussian_filter
from scipy.interpolate import UnivariateSpline
import pickle
import time as timelib
import sys
import pandas as pd
import re
import glob

class Quantity:
    def __init__(self, directory, delt=0.125):
        #self.exp = exp #name of the experiment
        self.input_dir = directory
        self.result_dir = directory + '/results/'
        if not os.path.exists(self.result_dir):
            os.makedirs(self.result_dir)
        self.domain = load_domain(self.input_dir)
        self.delt = delt
        
    def force_T1_a(self, normalise=True, plot_std=False, 
                   gap=10, prefilter=0, ignore_timesteps=20,
                   annular=False, beyond_gap=20):
        """
        Average force around a T1 transition for a circle of radius = 'gap'
        If normalise is False, the average is calculated by summing all forces
        and dividing by number of transitions. However, all the transitions are
        rescaled to the same length before taking the average. The missing 
        values are interpolated using spline interpolation.
        
        If 'normalise' is True, then each transition's force is first divided 
        by the average force of that T1 transition. This way, all Transitions 
        have an equal contribution.

        Parameters
        ----------
        normalise : bool, optional
            DESCRIPTION. The default is True. if True, then the force is 
            divided by the mean force for that particular T1 transition
        plot_std : bool, optional
            DESCRIPTION. The default is False. if True, then the standard 
            deviation bars are plotted
        gap : int, optional
            DESCRIPTION. The default is 10. Radius of the averaging circle 
            in terms of number of grid points
        prefilter : float, optional
            DESCRIPTION. The default is 0. Gaussian filter to be applied to
            stress before calculating norm of the force. 
        ignore_timesteps : int, optional
            DESCRIPTION. The default is 20. T1 transitions that start in those
            initial timesteps would be ignored
        annular : bool, optional
            DESCRIPTION. The default is False. Averaging is done in a annular 
            disk instead of a disk
        beyond_gap : int, optional
            DESCRIPTION. The default is 20. Outer radius of annular disk in 
            terms of number of grid points.

        Returns
        -------
        None.
        
        This program generates two plots
        1. Force as function of timesteps
        2. Multiple subplots where each plot corresponds to T1 transition of 
        lengths in the following ranges: [1, 0.25max], [0.25max, 0.5max], 
        [0.5max, 0.75max], [0.75max, max] where max is the length of largest T1
        transition  

        """
        force_status = False

        print('**********-------------------**********')
        print('2-norm of force during a T1 transition')
        if annular:
            print('annular mode: ',gap , beyond_gap)
        else:
            print('normal mode: ',gap , beyond_gap)
        
        st = timelib.time()
        grid_x = np.load(self.input_dir + '/stress_fields_500/grid_x.npy')
        grid_y = np.load(self.input_dir + '/stress_fields_500/grid_y.npy')
        steps_x = grid_x.shape[0]
        steps_y = grid_y.shape[1]
        dx = self.domain[0]/grid_x.shape[0]
        print('dx: ', dx)
        
        times = np.load(self.input_dir + '/stress_fields_500/timesteps.npy')
        print('max time: ', times[-1])

        #load T1 statistics
        with open(self.input_dir + '/T1_time_info.pkl', 'rb') as f:
            T1_time_info = pickle.load(f)
        T1_transitions = np.load(self.input_dir + '/T1_transitions.npy')
        T1_start = np.load(self.input_dir + '/T1_start.npy')
        T1_end = np.load(self.input_dir + '/T1_end.npy')
        

        #removing initial few T1 transitions
        ignore_ind = np.argwhere(T1_start > self.delt*ignore_timesteps)
        T1_end = T1_end[ignore_ind].flatten()
        T1_transitions = T1_transitions[ignore_ind].reshape(len(T1_end), 4)
        T1_start = T1_start[ignore_ind].flatten()
        T1_locations = T1_loc(self.input_dir, T1_transitions, T1_start, T1_end)
        T1_diff = T1_end - T1_start
        print('Total T1 transitions: ', len(T1_start))
        
        max_forces_T1 = np.zeros(len(T1_start)) #store maximum average force
        mean_forces_T1 = np.zeros(len(T1_start)) #store 'mean(in time)' average force
        std_forces_T1 = np.zeros(len(T1_start))

        fig, ax = plt.subplots(2,3, figsize=(14,12), facecolor='w')
        
        #T1 transition that lasts the longest, should be 21 time points long or
        #less because the algorithm used to detect T1 transitions can only 
        #handle T1 transitions of that length.
        
        max_length = int(np.max(T1_diff)/self.delt+1) 
        t_a = int(max_length/4)
        t_b = int(max_length/2)
        t_c = int(3*max_length/4)
        
        print('timestep divisions: ', t_a, t_b, t_c, max_length)
        
        forces_T1 = np.zeros([len(T1_transitions), max_length])
        forces_a = []
        forces_b = []
        forces_c = []
        forces_d = []
        
        for indT1, T1 in enumerate(T1_transitions):
            T1_pos = T1_locations[indT1]
            start_time = T1_start[indT1]
            end_time = T1_end[indT1]
            T1_times = np.arange(start_time, end_time+self.delt, self.delt)

            if not annular:
                index_grid = circle_coord(T1_pos, gap*dx, grid_x, grid_y)
            else:
                index_grid = annular_coord(T1_pos, gap*dx, beyond_gap*dx, grid_x, grid_y)
            
            forces = np.zeros(len(T1_times))
                
            for indT, time in enumerate(T1_times):
                if prefilter == 0:
                    if force_status:
                        try:
                            force_magnitude = np.load(
                                self.input_dir + '/stress_fields_500/force_norm' + '_{:06.3f}'.format(time) + '.npy')
                            
                        except:
                            sigma_00 = np.load(
                                self.input_dir + '/stress_fields_500/sigma_00' + '_{:06.3f}'.format(time) + '.npy')
                            sigma_11 = np.load(
                                self.input_dir + '/stress_fields_500/sigma_11' + '_{:06.3f}'.format(time) + '.npy')
                            sigma_01 = np.load(
                                self.input_dir + '/stress_fields_500/sigma_01' + '_{:06.3f}'.format(time) + '.npy')

                            force_magnitude = force_norm(
                                sigma_00, sigma_11, sigma_01, dx=dx, filter_width=prefilter)
                    else:
                        force_magnitude = np.load(
                                self.input_dir + '/stress_fields_500/free_energy' + '_{:06.3f}'.format(time) + '.npy')
                else:
                    sigma_00 = np.load(
                        self.input_dir + '/stress_fields_500/sigma_00' + '_{:06.3f}'.format(time) + '.npy')
                    sigma_11 = np.load(
                        self.input_dir + '/stress_fields_500/sigma_11' + '_{:06.3f}'.format(time) + '.npy')
                    sigma_01 = np.load(
                        self.input_dir + '/stress_fields_500/sigma_01' + '_{:06.3f}'.format(time) + '.npy')

                    force_magnitude = force_norm(
                        sigma_00, sigma_11, sigma_01, dx=dx, filter_width=prefilter)

                #for indf, filterwidth in enumerate(filters):
                #    force_magnitude = np.copy(rough_force_magnitude)
                #    if post_smoothen:
                #        force_magnitude = gaussian_filter(force_magnitude, filterwidth, mode = "wrap")

                force_in_neighbourhood = force_magnitude[index_grid[0], index_grid[1]] 
                forces[indT] = np.mean(force_in_neighbourhood)
                
            max_forces_T1[indT1] = np.max(forces)
            mean_forces_T1[indT1] = np.mean(forces)
            std_forces_T1[indT1] = np.std(forces)
            
            if normalise:
                forces /= np.mean(forces)
                
            old_indices = np.arange(0, len(forces))
            new_indices = np.linspace(0, len(forces)-1, max_length)
            new_forces = np.zeros(max_length)
            if(len(old_indices) > 3):
                spl = UnivariateSpline(old_indices, forces, k=3, s=0)
                # extending length of forces to match the length of the longest T1 transition
                new_forces = spl(new_indices)
            elif(len(old_indices) == 3):
                spl = UnivariateSpline(old_indices, forces, k=2, s=0)
                # extending length of forces to match the length of the longest T1 transition
                new_forces = spl(new_indices)
            elif(len(old_indices) == 2):
                new_forces = np.linspace(forces[0], forces[1], max_length)
            elif(len(old_indices) == 1):
                new_forces = np.ones(max_length)*forces[0]

            ax[0, 1].plot(forces)
            if(len(forces) <= t_a):
                old_indices = np.arange(0, len(forces))
                new_indices = np.linspace(0, len(forces)-1, t_a)
                new_forces_a = np.zeros(5)
                if(len(old_indices) > 3):
                    #ax[0,2].plot(forces)
                    spl = UnivariateSpline(old_indices, forces, k=3, s=0)
                    # extending length of forces to match the length of the longest T1 transition
                    new_forces_a = spl(new_indices)
                elif(len(old_indices) == 3):
                    spl = UnivariateSpline(old_indices, forces, k=2, s=0)
                    # extending length of forces to match the length of the longest T1 transition
                    new_forces_a = spl(new_indices)
                elif(len(old_indices) == 2):
                    new_forces_a = np.linspace(forces[0], forces[1], t_a)
                elif(len(old_indices) == 1):
                    new_forces_a = np.ones(5)*forces[0]
                forces_a.append(new_forces_a)
                #ax[0,2].plot(forces)
            elif(len(forces) <= t_b):
                old_indices = np.arange(0, len(forces))
                new_indices = np.linspace(0, len(forces)-1, t_b)
                spl = UnivariateSpline(old_indices, forces, k=3, s=0)  # 10
                # extending length of forces to match the length of the longest T1 transition
                new_forces_b = spl(new_indices)
                forces_b.append(new_forces_b)
                #ax[1,0].plot(forces)
            elif(len(forces) <= t_c):
                old_indices = np.arange(0, len(forces))
                new_indices = np.linspace(0, len(forces)-1, t_c)
                spl = UnivariateSpline(old_indices, forces, k=3, s=0)  # 15
                # extending length of forces to match the length of the longest T1 transition
                new_forces_c = spl(new_indices)
                forces_c.append(new_forces_c)
                #ax[1,1].plot(forces)
            elif(len(forces) < max_length+1):
                old_indices = np.arange(0, len(forces))
                new_indices = np.linspace(0, len(forces)-1, max_length)
                spl = UnivariateSpline(old_indices, forces, k=3, s=0)  # 21
                # extending length of forces to match the length of the longest T1 transition
                new_forces_d = spl(new_indices)
                forces_d.append(new_forces_d)
                #ax[1,2].plot(forces)

            forces_T1[indT1] = new_forces

        average_forces = np.mean(forces_T1, axis=0)
        std_forces = np.std(forces_T1, axis=0)

        avg_forces_a = np.mean(np.array(forces_a), axis=0)
        avg_forces_b = np.mean(np.array(forces_b), axis=0)
        avg_forces_c = np.mean(np.array(forces_c), axis=0)
        avg_forces_d = np.mean(np.array(forces_d), axis=0)

        std_forces_a = np.std(np.array(forces_a), axis=0)
        std_forces_b = np.std(np.array(forces_b), axis=0)
        std_forces_c = np.std(np.array(forces_c), axis=0)
        std_forces_d = np.std(np.array(forces_d), axis=0)

        if plot_std:
            ax[0, 0].errorbar(np.arange(0, np.max(T1_diff)+self.delt, self.delt),
                              average_forces, yerr=std_forces,
                              fmt='-', color='black',
                              ecolor='lightgray', elinewidth=3, capsize=0)
            ax[0, 2].errorbar(np.arange(0, len(avg_forces_a)),
                              avg_forces_a, yerr=std_forces_a,
                              fmt='-', color='black',
                              ecolor='lightgray', elinewidth=3, capsize=0)
            ax[1, 0].errorbar(np.arange(0, len(avg_forces_b)),
                              avg_forces_b, yerr=std_forces_b,
                              fmt='-', color='black',
                              ecolor='lightgray', elinewidth=3, capsize=0)
            ax[1, 1].errorbar(np.arange(0, len(avg_forces_c)),
                              avg_forces_c, yerr=std_forces_c,
                              fmt='-', color='black',
                              ecolor='lightgray', elinewidth=3, capsize=0)
            ax[1, 2].errorbar(np.arange(0, len(avg_forces_d)),
                              avg_forces_d, yerr=std_forces_d,
                              fmt='-', color='black',
                              ecolor='lightgray', elinewidth=3, capsize=0)

        else:
            ax[0, 0].plot(np.arange(0, np.max(T1_diff) +
                          self.delt, self.delt), average_forces, 'k')
            ax[0, 2].plot(avg_forces_a, 'k', linewidth=2)
            ax[1, 0].plot(avg_forces_b, 'k', linewidth=2)
            ax[1, 1].plot(avg_forces_c, 'k', linewidth=2)
            ax[1, 2].plot(avg_forces_d, 'k', linewidth=2)

        ax[0, 0].set_xlabel('Time (rescaled)')  # , fontsize=22)
        if force_status:
            ax[0, 0].set_ylabel('Forces during \n T1 transition')  # , fontsize=22)
        else:
            ax[0, 0].set_ylabel('Energy during \n T1 transition')  # , fontsize=22)
        ax[0, 0].grid()
        #ax[0,1].set_xlabel('Time (rescaled)', fontsize=22)
        #ax[0,1].set_ylabel('Forces during \n T1 transition', fontsize=22)
        ax[0, 1].grid()
        ax[0, 2].grid()
        ax[1, 0].grid()
        ax[1, 1].grid()
        ax[1, 2].grid()
        ax[0, 2].set_title('#T1 = ' + str(len(forces_a)))
        ax[1, 0].set_title('#T1 = ' + str(len(forces_b)))
        ax[1, 1].set_title('#T1 = ' + str(len(forces_c)))
        ax[1, 2].set_title('#T1 = ' + str(len(forces_d)))

        plt.tight_layout()

        if force_status:
            filename = 'force_T1'
        else:
            filename = 'energy_T1'
        if normalise:
            filename += '_n'
        if plot_std:
            filename += '_std'
        if annular:
            filename += '_ann'

        plt.savefig(self.result_dir + '/' + filename + '.png', dpi=200)
        et = timelib.time()
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')

        #single plot
        fig, ax = plt.subplots(facecolor='w')
        if plot_std:
            ax.errorbar(np.arange(0, np.max(T1_diff)+self.delt, self.delt),
                        average_forces, yerr=std_forces,
                        fmt='-', color='black',
                        ecolor='lightgray', elinewidth=3, capsize=0)
        else:
            ax.plot(np.arange(0, np.max(T1_diff)+self.delt, self.delt),
                    average_forces, 'k')
        ax.set_xlabel('Time (rescaled)', fontsize=22)
        if normalise == False:
            if force_status:
                ax.set_ylabel('forces during \n T1 transition', fontsize=22)
            else:
                ax.set_ylabel('energy during \n T1 transition', fontsize=22)
            
        else:
            if force_status:
                ax.set_ylabel('rescaled forces during \n T1 transition', fontsize=22)
            else:
                ax.set_ylabel('rescaled energy during \n T1 transition', fontsize=22)
            
        ax.grid()
        if annular:
            ax.set_title('range: ' + str(gap) + ', ' + str(beyond_gap))
        else:
            ax.set_title('gap: ' + str(gap))
        plt.tight_layout()
        plt.savefig(self.result_dir + '/' + filename + '_s.png', dpi=200)
        
        #saving the forces
        np.save(self.result_dir + '/' + filename + '.npy', average_forces)
        
        #plotting max forces
        fig, ax = plt.subplots()
        #Freedman Diaconis rule to calculate the 'appropriate' bin width
        q25, q75 = np.percentile(max_forces_T1, [25, 75])
        bin_width = 2 * (q75 - q25) * len(max_forces_T1) ** (-1/3)
        bins = int(round((max_forces_T1.max() - max_forces_T1.min()) / bin_width))
        
        N, bins, patches = ax.hist(max_forces_T1, bins=bins, rwidth=0.9)
        fracs = N / N.max()
        colornorm = colors.Normalize(fracs.min(), fracs.max())
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.viridis(colornorm(thisfrac))
            thispatch.set_facecolor(color)
        if force_status:
            ax.set_xlabel('Max force')
        else:
            ax.set_xlabel('Max energy')
        plt.savefig(self.result_dir + '/' + filename + '_maxdist.png', dpi=200)
            
        #plotting mean forces
        fig, ax = plt.subplots()
        #Freedman Diaconis rule to calculate the 'appropriate' bin width
        q25, q75 = np.percentile(mean_forces_T1, [25, 75])
        bin_width = 2 * (q75 - q25) * len(mean_forces_T1) ** (-1/3)
        bins =int(round((mean_forces_T1.max() - mean_forces_T1.min()) / bin_width))
        
        N, bins, patches = ax.hist(mean_forces_T1, bins=bins, rwidth=0.9)
        fracs = N / N.max()
        colornorm = colors.Normalize(fracs.min(), fracs.max())
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.viridis(colornorm(thisfrac))
            thispatch.set_facecolor(color)
        plt.savefig(self.result_dir + '/' + filename + '_meandist.png', dpi=200)
        if force_status:
            ax.set_xlabel('Mean force')
        else:
            ax.set_xlabel('Mean energy')
        #plotting std forces
        fig, ax = plt.subplots()
        #Freedman Diaconis rule to calculate the 'appropriate' bin width
        q25, q75 = np.percentile(std_forces_T1, [25, 75])
        bin_width = 2 * (q75 - q25) * len(std_forces_T1) ** (-1/3)
        bins = int(round((std_forces_T1.max() - std_forces_T1.min()) / bin_width))
        
        N, bins, patches = ax.hist(std_forces_T1, bins=bins, rwidth=0.9)
        fracs = N / N.max()
        colornorm = colors.Normalize(fracs.min(), fracs.max())
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.viridis(colornorm(thisfrac))
            thispatch.set_facecolor(color)
        plt.savefig(self.result_dir + '/' + filename + '_stddist.png', dpi=200)
        if force_status:
            ax.set_xlabel('std force')
        else:
            ax.set_xlabel('std energy')
        #plotting t1 Time
        fig, ax = plt.subplots()
        #Freedman Diaconis rule to calculate the 'appropriate' bin width
        q25, q75 = np.percentile(T1_diff, [25, 75])
        bin_width = 2 * (q75 - q25) * len(T1_diff) ** (-1/3)
        bins = int(round((T1_diff.max() - T1_diff.min()) / bin_width))
        
        N, bins, patches = ax.hist(T1_diff, bins=bins, rwidth=0.9)
        fracs = N / N.max()
        colornorm = colors.Normalize(fracs.min(), fracs.max())
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.viridis(colornorm(thisfrac))
            thispatch.set_facecolor(color)
        print('average T1 transition duration is ', np.mean(T1_diff))
        plt.savefig(self.result_dir + '/' + filename + '_timedist.png', dpi=200)
        ax.set_xlabel('T1 duration')

            
        #scatter plot of max force vs duration of T1
        fig, ax = plt.subplots()
        ax.scatter(max_forces_T1, T1_diff)
        if force_status:
            ax.set_xlabel('Max force')
        else:
            ax.set_xlabel('Max energy')
        ax.set_ylabel('T1 duration')
        plt.savefig(self.result_dir + '/' + filename + '_force_time.png', dpi=200)

        

    def force_prepost_T1(self, normalise=False, plot_std=False, equalise_T1num=False,
                               gaps=np.array([10, 20, 30]), beyond_index=70,
                               time_gaps = np.arange(0.0, 30.5, 0.5), 
                               prefilter=0, postfilter=0, 
                               start_index=20, end_index=-1):

        """
        Plots average force around a T1 transition, before and after a T1 transition
        starts and ends. The average force is calculated for disk and annular disk

        Parameters
        ----------
        normalise : bool, optional
            DESCRIPTION. The default is False. If True then force for each T1
            transtion is divided by hte mean force. Not yet implemented
        plot_std : bool, optional
            DESCRIPTION. The default is False. If True then standard deviation 
            is also plotted. Not yet implemented
        gaps : int, optional
            DESCRIPTION. The default is np.array([10, 20, 30]). The radius of 
            disk and inner radisu of annular disk in terms of number of grid 
            points
        beyond_index : int, optional
            DESCRIPTION. The default is 70. The outer radius of annular disk
        time_gaps : 1D array float, optional
            DESCRIPTION. The default is np.arange(0.0, 70.5, 0.5). The time
            array before and after a t1 transition for which average is 
            calculated. So if array is np.arange(0.0, 70.5, 0.5) then force is
            calculated for time upto 70 units before and after T1 transition
            and every 0.5 units is taken into account
        prefilter : float, optional
            DESCRIPTION. The default is 0. Gaussian filter size to apply to 
            stress before finding norm of force
        postfilter : float, optional. 
            DESCRIPTION. The default is 0. Gaussian filter to apply to force 
            norm before finding average force
        start_index : int, optional
            DESCRIPTION. The default is 100. data before this time index is ignored
        end_index : int, optional
            DESCRIPTION. The default is -1. data after this time index is ignored

        Returns
        -------
        None.

        """
        force_status = False
        print('**********-------------------**********')
        print('2-norm of force before and after a T1 transition')
        st = timelib.time()
        
        eps = 0.1
        beyond_gaps = [(beyond_index - g) for g in gaps] 
        
        grid_x = np.load(self.input_dir + '/stress_fields_500/grid_x.npy')
        grid_y = np.load(self.input_dir + '/stress_fields_500/grid_y.npy')
        steps_x = grid_x.shape[0]
        steps_y = grid_y.shape[1]
        dx = self.domain[0]/grid_x.shape[0]
        print('dx: ', dx)
        print('gaps: ', gaps)
        print('beyond gaps:', beyond_gaps)

        
        times = np.load(self.input_dir + '/stress_fields_500/timesteps.npy')
        selected_times = times[start_index:end_index]
        print('start time:', selected_times[0])
        print('end time:', selected_times[-1])
        
        #load T1 statistics
        with open(self.input_dir + '/T1_time_info.pkl', 'rb') as f:
            T1_time_info = pickle.load(f)
        T1_transitions = np.load(self.input_dir + '/T1_transitions.npy')
        T1_start = np.load(self.input_dir + '/T1_start.npy')
        T1_end = np.load(self.input_dir + '/T1_end.npy')
        T1_locations = T1_loc(self.input_dir, T1_transitions, T1_start, T1_end)
        print('Total T1 transitions: ', len(T1_start))
        
        if equalise_T1num:
            ignore_ind = np.argwhere(T1_start-time_gaps[-1] > selected_times[0])
            T1_end = T1_end[ignore_ind].flatten()
            T1_transitions = T1_transitions[ignore_ind].reshape(len(T1_end), 4)
            T1_start = T1_start[ignore_ind].flatten()
            
            ignore_ind = np.argwhere(T1_end+time_gaps[-1] < selected_times[-1])
            T1_end = T1_end[ignore_ind].flatten()
            T1_transitions = T1_transitions[ignore_ind].reshape(len(T1_end), 4)
            T1_start = T1_start[ignore_ind].flatten()
            
            T1_locations = T1_loc(self.input_dir, T1_transitions, T1_start, T1_end)
            T1_diff = T1_end - T1_start
            print('Total T1 transitions: ', len(T1_start))
        
        S0_arr, S1_arr = load_property(self.input_dir, 'nematic')
        S_mag = np.linalg.norm([S0_arr, S1_arr], axis=0)
        mean_S = np.mean(S_mag, axis = 1)
        std_S = np.std(S_mag, axis = 1)
        
        ranks = load_ranks(self.input_dir)
        perimeter = np.zeros([len(ranks), len(times)])
        radius = load_property(self.input_dir, 'radius')
        area=97.8#change this based on simulations
        area = np.pi*(radius**2)
        for rank in ranks:
            perimeter[rank] = np.load(self.input_dir + '/shape/perimeter_' + str(rank) + '.npy')[1:]
        shape_index = (perimeter/np.sqrt(area.T)).T
        
        avg_vel_mag = np.linalg.norm(load_property(self.input_dir, 'velocity'), axis=0)
        avg_acl_mag = np.zeros([len(times), len(ranks)])
        avg_acl_mag[1:, :] = np.linalg.norm(load_property(self.input_dir, 'acceleration'), axis=0)
        
        #Global average of cell elongation
        fig, ax = plt.subplots(facecolor='w')
        ax.errorbar(times[3:], mean_S[3:], yerr=std_S[3:], fmt='-', color = 'black', ecolor = 'lightgray', elinewidth=3, capsize=0)
        ax.set_xlabel('Time')
        ax.set_ylabel('$|S|$')
        ax.grid()
        ax.set_ylim(0, 14)
        plt.tight_layout()
        
        start_forces = np.zeros([len(T1_start), len(gaps), len(time_gaps)])
        end_forces = np.zeros([len(T1_start), len(gaps), len(time_gaps)])
        
        #below two arrays will be used as a mask
        #this will be 0 if at that time, the time overlaps with start time - time gap
        start_mask = np.ones([len(T1_start), len(gaps), len(time_gaps)], dtype=bool) 
        #this will be 0 if at that time, the time overlaps with end time + time gap
        end_mask = np.ones([len(T1_start), len(gaps), len(time_gaps)], dtype=bool) 
        
        #below arrays will store values of forces not around a T1 but a bit beyond it.
        start_forces_beyond = np.zeros([len(T1_start), len(beyond_gaps), len(time_gaps)])
        end_forces_beyond = np.zeros([len(T1_start), len(beyond_gaps), len(time_gaps)])
        start_mask_beyond = np.ones([len(T1_start), len(beyond_gaps), len(time_gaps)], dtype=bool)
        end_mask_beyond = np.ones([len(T1_start), len(beyond_gaps), len(time_gaps)], dtype=bool)
        
        #force mean across the whole domain before T1
        mean_force_total_s = np.zeros([len(T1_start), len(time_gaps)]) 
        mean_force_total_mask_s = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        #force mean across the whole domain after T1
        mean_force_total_e = np.zeros([len(T1_start), len(time_gaps)]) 
        mean_force_total_mask_e = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        
        #stores time data before and after T1 transitions
        T1_start_gaps = np.zeros([len(T1_start), len(time_gaps)])
        T1_end_gaps = np.zeros([len(T1_end), len(time_gaps)])
        for indT1, start in enumerate(T1_start):
            T1_start_gaps[indT1] = T1_start[indT1] - np.flip(time_gaps)
        for indT1, end in enumerate(T1_end):
            T1_end_gaps[indT1] = T1_end[indT1] + time_gaps
        
        index_grid_list = [[[] for T1 in T1_start] for gap in gaps]
        stencil_grid_list = [[[] for T1 in T1_start] for gap in gaps]
        for indT1 in range(len(T1_start)):
            for indg, gap in enumerate(gaps):
                loc = T1_locations[indT1]
                index_grid = circle_coord(loc, gap*dx, grid_x, grid_y)
                stencil_grid = annular_coord(loc, gap*dx, beyond_index*dx, grid_x, grid_y)
                index_grid_list[indg][indT1] = index_grid
                stencil_grid_list[indg][indT1] = stencil_grid
                
        Smag_start = np.zeros([len(T1_start), len(time_gaps)])
        Smag_start_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        Smag_end = np.zeros([len(T1_start), len(time_gaps)])
        Smag_end_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        
        #shape index
        p_start = np.zeros([len(T1_start), len(time_gaps)])
        p_start_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        p_end = np.zeros([len(T1_start), len(time_gaps)])
        p_end_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        
        #avg vel
        avg_vel_start = np.zeros([len(T1_start), len(time_gaps)])
        avg_vel_start_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        avg_vel_end = np.zeros([len(T1_start), len(time_gaps)])
        avg_vel_end_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        
        #avg acl
        avg_acl_start = np.zeros([len(T1_start), len(time_gaps)])
        avg_acl_start_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
        avg_acl_end = np.zeros([len(T1_start), len(time_gaps)])
        avg_acl_end_mask = np.ones([len(T1_start), len(time_gaps)], dtype=bool)
                
        for indt, time in enumerate(selected_times):
            if prefilter < eps:
                if force_status:
                    try:
                        force_magnitude = np.load(self.input_dir + '/stress_fields_500/force_norm' + '_{:06.3f}'.format(time) +'.npy')
                        
                    except:        
                        sigma_00 = np.load(self.input_dir + '/stress_fields_500/sigma_00' + '_{:06.3f}'.format(time) +'.npy')
                        sigma_11 = np.load(self.input_dir + '/stress_fields_500/sigma_11' + '_{:06.3f}'.format(time) +'.npy')
                        sigma_01 = np.load(self.input_dir + '/stress_fields_500/sigma_01' + '_{:06.3f}'.format(time) +'.npy')
        
                        force_magnitude = force_norm(sigma_00, sigma_11, sigma_01, dx=dx, filter_width = prefilter)
                else:
                    force_magnitude = np.load(
                        self.input_dir + '/stress_fields_500/free_energy' + '_{:06.3f}'.format(time) + '.npy')
            else:
                sigma_00 = np.load(self.input_dir + '/stress_fields_500/sigma_00' + '_{:06.3f}'.format(time) +'.npy')
                sigma_11 = np.load(self.input_dir + '/stress_fields_500/sigma_11' + '_{:06.3f}'.format(time) +'.npy')
                sigma_01 = np.load(self.input_dir + '/stress_fields_500/sigma_01' + '_{:06.3f}'.format(time) +'.npy')
    
                force_magnitude = force_norm(sigma_00, sigma_11, sigma_01, dx=dx, filter_width = prefilter)
            
            if postfilter > eps:
                force_magnitude = gaussian_filter(force_magnitude, postfilter, mode = "wrap")
            
            mean_force = np.mean(force_magnitude)
            for indT1 in range(len(T1_start)):
                
                for indtg in range(len(T1_start_gaps[indT1])):
                    if time == T1_start_gaps[indT1][indtg]:
                        mean_force_total_s[indT1][indtg] = mean_force
                        mean_force_total_mask_s[indT1][indtg] = 0
        
                        for indg, gap in enumerate(gaps):
                            index_grid = index_grid_list[indg][indT1]
                            force_in_neighbourhood = force_magnitude[index_grid[0], index_grid[1]]
                            start_forces[indT1][indg][indtg] = np.mean(force_in_neighbourhood)
                            start_mask[indT1][indg][indtg] = 0
                
                            #beyond forces
                            stencil_grid = stencil_grid_list[indg][indT1]
                            force_beyond = force_magnitude[stencil_grid[0], stencil_grid[1]]
                            start_forces_beyond[indT1][indg][indtg] = np.mean(force_beyond)
                            start_mask_beyond[indT1][indg][indtg] = 0
                            
                        Smag_start[indT1, indtg] = np.mean(S_mag[start_index+indt][T1_transitions[indT1]])
                        Smag_start_mask[indT1, indtg] = 0
                        p_start[indT1, indtg] = np.mean(shape_index[start_index+indt][T1_transitions[indT1]])
                        p_start_mask[indT1, indtg] = 0
                        avg_vel_start[indT1, indtg] = np.mean(avg_vel_mag[start_index+indt][T1_transitions[indT1]])
                        avg_vel_start_mask[indT1, indtg] = 0
                        avg_acl_start[indT1, indtg] = np.mean(avg_acl_mag[start_index+indt][T1_transitions[indT1]])
                        avg_acl_start_mask[indT1, indtg] = 0
                
                for indtg in range(len(T1_end_gaps[indT1])):
                    if time == T1_end_gaps[indT1][indtg]:
                        mean_force_total_e[indT1][indtg] = mean_force
                        mean_force_total_mask_e[indT1][indtg] = 0
                        for indg, gap in enumerate(gaps):
                            index_grid = index_grid_list[indg][indT1]
                            force_in_neighbourhood = force_magnitude[index_grid[0], index_grid[1]]
                            end_forces[indT1][indg][indtg] = np.mean(force_in_neighbourhood)
                            end_mask[indT1][indg][indtg] = 0
                
                            #beyond forces
                            stencil_grid = stencil_grid_list[indg][indT1]
                            force_beyond = force_magnitude[stencil_grid[0], stencil_grid[1]]
                            end_forces_beyond[indT1][indg][indtg] = np.mean(force_beyond)
                            end_mask_beyond[indT1][indg][indtg] = 0
                            
                        Smag_end[indT1, indtg] = np.mean(S_mag[start_index+indt][T1_transitions[indT1]])
                        Smag_end_mask[indT1, indtg] = 0
                        p_end[indT1, indtg] = np.mean(shape_index[start_index+indt][T1_transitions[indT1]])
                        p_end_mask[indT1, indtg] = 0
                        avg_vel_end[indT1, indtg] = np.mean(avg_vel_mag[start_index+indt][T1_transitions[indT1]])
                        avg_vel_end_mask[indT1, indtg] = 0
                        avg_acl_end[indT1, indtg] = np.mean(avg_acl_mag[start_index+indt][T1_transitions[indT1]])
                        avg_acl_end_mask[indT1, indtg] = 0
         
        """
        #The following optional part makes sure that all time instances have 
        #equal number of T1 transitions contributing to them
        if equalise_T1num:
            #assuming that the first time instant has the lowest number of T1 transtion 
            #contributing to it
            for indL in range(start_mask.shape[-1]):
                start_mask[:, :, indL] = start_mask[:, :, 0] 
            for indL in range(start_mask.shape[-1]):
                end_mask[:, :, indL] = start_mask[:, :, 0] 
        """
        #counting number of T1 transitions that contribute to the average
        start_count = np.sum(1-start_mask, axis=0)
        end_count = np.sum(1-end_mask, axis=0)
        
        #finding mean and standard deviation across T1 transitions
        start_forces_m = np.ma.array(data=start_forces, mask=start_mask)
        end_forces_m = np.ma.array(data=end_forces, mask=end_mask)
        start_forces_mean = np.mean(start_forces_m, axis=0)
        end_forces_mean = np.mean(end_forces_m, axis=0)
        start_forces_std = np.std(start_forces_m, axis=0)
        end_forces_std = np.std(end_forces_m, axis=0)
        
        #finding mean and standard deviation across T1 transitions
        start_forces_beyond_m = np.ma.array(data=start_forces_beyond, mask=start_mask_beyond)
        end_forces_beyond_m = np.ma.array(data=end_forces_beyond, mask=end_mask_beyond)
        start_forces_beyond_mean = np.mean(start_forces_beyond_m, axis=0)
        end_forces_beyond_mean = np.mean(end_forces_beyond_m, axis=0)
        start_forces_beyond_std = np.std(start_forces_beyond_m, axis=0)
        end_forces_beyond_std = np.std(end_forces_beyond_m, axis=0)
        
        if normalise:
            #finding mean across time which will be used to normalise
            start_forces_timemean = np.mean(start_forces_m, axis=2)
            end_forces_timemean = np.mean(end_forces_m, axis=2)
            #getting a shared mean for before and after T1
            forces_timemean = (start_forces_timemean + end_forces_timemean)/2 
            
            start_forces_beyond_timemean = np.mean(start_forces_beyond_m, axis=2)
            end_forces_beyond_timemean = np.mean(end_forces_beyond_m, axis=2)
            #getting a shared mean for before and after T1
            forces_beyond_timemean = (start_forces_beyond_timemean + 
                                      end_forces_beyond_timemean)/2
            
            #expanding dimensions to help with division
            forces_timemean = np.expand_dims(forces_timemean, -1)
            forces_beyond_timemean = np.expand_dims(forces_beyond_timemean, -1)

            #dividing by time mean for all T1 transitions to get time normalised versions
            start_forces_m/= forces_timemean
            end_forces_m /= forces_timemean
            start_forces_beyond_m /= forces_beyond_timemean
            end_forces_beyond_m /= forces_beyond_timemean
            
            start_forces_mean = np.mean(start_forces_m, axis=0)
            end_forces_mean = np.mean(end_forces_m, axis=0)
            start_forces_std = np.std(start_forces_m, axis=0)
            end_forces_std = np.std(end_forces_m, axis=0)
            
            start_forces_beyond_mean = np.mean(start_forces_beyond_m, axis=0)
            end_forces_beyond_mean = np.mean(end_forces_beyond_m, axis=0)
            start_forces_beyond_std = np.std(start_forces_beyond_m, axis=0)
            end_forces_beyond_std = np.std(end_forces_beyond_m, axis=0)
        
        #mean domain force
        mean_force_total_sm = np.ma.array(data=mean_force_total_s, mask=mean_force_total_mask_s)
        mean_force_total_em = np.ma.array(data=mean_force_total_e, mask=mean_force_total_mask_e)
        mean_forces_tsm = np.mean(mean_force_total_sm, axis=0) #mean across all times if mask is 0
        mean_forces_tem = np.mean(mean_force_total_em, axis=0) #mean across all times if mask is 0
        std_forces_tsm = np.std(mean_force_total_sm, axis=0) #std across all times if mask is 0
        std_forces_tem = np.std(mean_force_total_em, axis=0) #std across all times if mask is 0
        
        
        plot_times = np.concatenate((-np.flip(np.array(time_gaps)), np.array(time_gaps)))
            
        et = timelib.time()
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')
        mfc = ['magenta', 'cyan', 'yellow', 'lime', 'gray']
        dfc = ['red', 'blue', 'orange', 'green', 'black']
        mfc2 = ['yellow', 'magenta', 'cyan', 'lime', 'grey']
        
        mean_total_forces = np.concatenate((mean_forces_tsm, mean_forces_tem))
        std_total_forces = np.concatenate((std_forces_tsm, std_forces_tem))
        
        fig, ax = plt.subplots(figsize=(8,8), facecolor='w')
        for indg, g in enumerate(gaps):
            forces = np.concatenate((start_forces_mean[indg,:], end_forces_mean[indg,:])) #method2
            ax.plot(plot_times[::2], forces[::2], marker='.', 
                    markerfacecolor = mfc[indg], label = str(g), 
                    color='k', linewidth=1, markersize=10)
            if indg==0:
                save_forces = np.ma.filled(forces, np.nan)
        if not normalise:
            ax.plot(plot_times, mean_total_forces, linestyle = 'dashed', linewidth=3)
         
        ax.legend(fontsize=22)
        ax.grid()
        if force_status:
            ax.set_ylabel('Force', fontsize=22)
        else:
            ax.set_ylabel('Energy', fontsize=22)
        
        plt.tight_layout()
        
        postfix = '_'
        if equalise_T1num:
            postfix += 'eq_'
        if normalise:
            postfix += 'n_'
        if plot_std:
            postfix += 'std_'
        
        if force_status:
            filename = 'force_prepost_T1'
        else:
            filename = 'energy_prepost_T1'
        plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
        np.save(self.result_dir + '/' + filename + postfix + '.npy', save_forces)
        np.save(self.result_dir + '/' + filename + postfix + 'times.npy', plot_times)
        
        fig, ax = plt.subplots(figsize=(8,8), facecolor='w')
        for indg, g in enumerate(beyond_gaps):
            beyond_forces = np.concatenate((start_forces_beyond_mean[indg,:], 
                                            end_forces_beyond_mean[indg,:])) 
            ax.plot(plot_times[::2], beyond_forces[::2], marker='.', 
                    markerfacecolor = mfc[indg], label = str(g), 
                    color='k', linewidth=1, markersize=10)
            if indg == 0:
                save_beyond_forces = np.ma.filled(beyond_forces, np.nan)
        if not normalise:
            ax.plot(plot_times, mean_total_forces, linestyle = 'dashed', linewidth=3)
        ax.legend(fontsize=22)
        ax.grid()
        if force_status:
            ax.set_ylabel('Beyond Force', fontsize=22)
        else:
            ax.set_ylabel('Beyond Energy', fontsize=22)
        plt.tight_layout()
        
        if force_status:
            filename = 'beyondforce_prepost_T1'
        else:
            filename = 'beyondenergy_prepost_T1'
        plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
        np.save(self.result_dir + '/' + filename + postfix + '.npy', save_beyond_forces)
        
        fig, ax = plt.subplots(1,2,figsize=(16,8), facecolor='w')
        for indg, g in enumerate(gaps):
            forces = np.concatenate((start_forces_mean[indg,:], end_forces_mean[indg,:])) 
            ax[0].plot(plot_times[::2], gaussian_filter(forces[::2], 1), marker='.', 
                    markerfacecolor = mfc[indg], label = str(g), 
                    color='k', linewidth=1, markersize=10)
            
            beyond_forces = np.concatenate((start_forces_beyond_mean[indg,:], 
                                            end_forces_beyond_mean[indg,:])) 
            ax[1].plot(plot_times[::2], gaussian_filter(beyond_forces[::2], 1), marker='.', 
                    markerfacecolor = mfc[indg], label = str(beyond_gaps[indg]), 
                    color='k', linewidth=1, markersize=10)
            
        for axes in ax:
            axes.plot(plot_times, mean_total_forces, linestyle = 'dashed', 
                      linewidth=3, color='indigo')
            axes.grid()
            axes.legend(fontsize=22)
        plt.tight_layout()
        
        if force_status:
            ax[0].set_ylabel('Force (filtered)')
            ax[1].set_ylabel('Beyond force (filtered)')
        else:
            ax[0].set_ylabel('Energy (filtered)')
            ax[1].set_ylabel('Beyond energy (filtered)')
        
        if force_status:
            filename = 'force_prepost_T1'
        else:
            filename = 'energy_prepost_T1'
        
        postfix += 'sm_'
        plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
        
        #plotting std
        fig, ax = plt.subplots(figsize=(8,8), facecolor='w')
        for indg, g in enumerate(gaps):
            forces = np.concatenate((start_forces_mean[indg,:], end_forces_mean[indg,:]))
            forces_std = np.concatenate((start_forces_std[indg,:], end_forces_std[indg,:]))
            ax.plot(plot_times[::2], forces_std[::2], marker='.', 
                    markerfacecolor = mfc[indg], label = str(g), 
                    color='k', linewidth=1, markersize=10)
            if indg == 0:
                save_forces_std = np.ma.filled(forces_std, np.nan)
        #ax.plot(plot_times, mean_total_forces, linestyle = 'dashed', linewidth=3)
        #ax.legend(fontsize=22)
        ax.legend(bbox_to_anchor=(1.05,1.025), loc="best", fontsize=16)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.75)
        ax.grid()
        ax.set_ylim(bottom=0.0)
        
        if force_status:
            ax.set_ylabel('Force std', fontsize=22)
        else:
            ax.set_ylabel('Energy std', fontsize=22)
        
        plt.tight_layout()
        if force_status:
            filename = 'force_prepost_T1_std'
        else:
            filename = 'energy_prepost_T1_std'
        plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
        np.save(self.result_dir + '/' + filename + postfix + '.npy', save_forces_std)
        
        fig, ax = plt.subplots(figsize=(8,8), facecolor='w')
        for indg, g in enumerate(gaps[0:3]):
            forces = np.concatenate((start_forces_mean[indg,:], end_forces_mean[indg,:]))
            forces_std = np.concatenate((start_forces_std[indg,:], end_forces_std[indg,:]))
            ax.plot(plot_times[::2], forces[::2], 
                    color = dfc[indg], label = str(g), 
                    linewidth=1)
            ax.fill_between(plot_times, forces-forces_std, forces+forces_std,
                            alpha = 0.3, facecolor=mfc[indg], edgecolor=dfc[indg])
        #ax.plot(plot_times, mean_total_forces, linestyle = 'dashed', linewidth=3)
        ax.legend(fontsize=22)
        ax.grid()
        
        if force_status:
            ax.set_ylabel('Force std', fontsize=22)
        else:
            ax.set_ylabel('Energy std', fontsize=22)

        plt.tight_layout()
        
        fig, ax = plt.subplots(figsize=(8,8), facecolor='w')
        for indg, g in enumerate(beyond_gaps):
            if indg>=3:
                forces_beyond = np.concatenate((start_forces_beyond_mean[indg,:], end_forces_beyond_mean[indg,:]))
                forces_std_beyond = np.concatenate((start_forces_beyond_std[indg,:], end_forces_beyond_std[indg,:]))
                ax.plot(plot_times[::2], forces_beyond[::2], 
                        color = dfc[indg], label = str(beyond_gaps[indg]), 
                        linewidth=1)
                ax.fill_between(plot_times, forces_beyond-forces_std_beyond, forces_beyond+forces_std_beyond,
                                alpha = 0.3, facecolor=mfc[indg], edgecolor=dfc[indg])
        #ax.plot(plot_times, mean_total_forces, linestyle = 'dashed', linewidth=3)
        ax.legend(fontsize=22)
        ax.grid()
        if force_status:
            ax.set_ylabel('Annular Force std', fontsize=22)
        else:
            ax.set_ylabel('Annular Energy std', fontsize=22)
        
        plt.tight_layout()
        
        #plotting counts of T1 used for statistics
        fig, ax = plt.subplots(figsize=(8,8), facecolor='w')
        for indg, g in enumerate(gaps):
            count_T1 = np.concatenate((start_count[indg,:], end_count[indg,:]))
            ax.plot(plot_times, count_T1, color=mfc[indg], label=str(g))
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('T1 statistics count', fontsize=22)
        plt.tight_layout()
        filename = 'T1_statistics_count'
        plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
        np.save(self.result_dir + '/' + filename + postfix + '.npy', count_T1)
        
        #force jump statistics: jump in force before and after T1 transtion
        fig, ax = plt.subplots(facecolor='w')
        indg=0
        force_jumps = start_forces_m[:, indg,-1] -  end_forces_m[:, indg,0]
        print('count', len(force_jumps))
        #Freedman Diaconis rule to calculate the 'appropriate' bin width
        q25, q75 = np.percentile(force_jumps, [25, 75])
        bin_width = 2 * (q75 - q25) * len(force_jumps) ** (-1/3)
        bins = int(round((force_jumps.max() - force_jumps.min()) / bin_width))
        
        N, bins, patches = ax.hist(force_jumps, bins=bins, rwidth=0.9)
        fracs = N / N.max()
        colornorm = colors.Normalize(fracs.min(), fracs.max())
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.viridis(colornorm(thisfrac))
            thispatch.set_facecolor(color)
        
        ax.set_ylabel('Count')
        
        if force_status:
            ax.set_xlabel('Force jump during a T1 transition')
            filename = 'T1_force_jump'
        else:
            ax.set_xlabel('Energy jump during a T1 transition')
            filename = 'T1_energy_jump'
        
        plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
            
        if equalise_T1num:
            #pre force jump statistics: jump in force from -10 to 0
            fig, ax = plt.subplots(facecolor='w')
            indg=0
            delta_t = time_gaps[1]-time_gaps[0]
            jump_size = int(10/delta_t)
            force_jumps = start_forces_m[:, indg,-1] -  start_forces_m[:, indg, -jump_size]
            
            print('count', len(force_jumps))
            #Freedman Diaconis rule to calculate the 'appropriate' bin width
            q25, q75 = np.percentile(force_jumps, [25, 75])
            bin_width = 2 * (q75 - q25) * len(force_jumps) ** (-1/3)
            bins = int(round((force_jumps.max() - force_jumps.min()) / bin_width))
            
            N, bins, patches = ax.hist(force_jumps, bins=bins, rwidth=0.9)
            fracs = N / N.max()
            colornorm = colors.Normalize(fracs.min(), fracs.max())
            for thisfrac, thispatch in zip(fracs, patches):
                color = plt.cm.viridis(colornorm(thisfrac))
                thispatch.set_facecolor(color)
                
            ax.set_ylabel('Count')
            if force_status:
                ax.set_xlabel('Force jump before T1 transition')
                filename = 'T1_force_prejump'
            else:
                ax.set_xlabel('Energy jump before T1 transition')
                filename = 'T1_energy_prejump'
            plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
                
            #post force jump statistics: jump in force from +0 to 10
            fig, ax = plt.subplots(facecolor='w')
            indg=0
            delta_t = time_gaps[1]-time_gaps[0]
            jump_size = int(10/delta_t)
            force_jumps = end_forces_m[:, indg,0] -  end_forces_m[:, indg, jump_size]
            
            print('count', len(force_jumps))
            #Freedman Diaconis rule to calculate the 'appropriate' bin width
            q25, q75 = np.percentile(force_jumps, [25, 75])
            bin_width = 2 * (q75 - q25) * len(force_jumps) ** (-1/3)
            bins = int(round((force_jumps.max() - force_jumps.min()) / bin_width))
            
            N, bins, patches = ax.hist(force_jumps, bins=bins, rwidth=0.9)
            fracs = N / N.max()
            colornorm = colors.Normalize(fracs.min(), fracs.max())
            for thisfrac, thispatch in zip(fracs, patches):
                color = plt.cm.viridis(colornorm(thisfrac))
                thispatch.set_facecolor(color)
                
            ax.set_ylabel('Count')
            if force_status:
                ax.set_xlabel('Force jump after T1 transition')
                filename = 'T1_force_postjump'
            else:
                ax.set_xlabel('Energy jump after T1 transition')
                filename = 'T1_energy_postjump'
            plt.savefig(self.result_dir + '/' + filename + postfix + '.png', dpi=200)
        
        #finding mean and standard deviation of elongation across T1 transitions
        Smag_start_m = np.ma.array(data=Smag_start, mask=Smag_start_mask)
        Smag_end_m = np.ma.array(data=Smag_end, mask=Smag_end_mask)
        Smag_start_mean = np.mean(Smag_start_m, axis=0)
        Smag_end_mean = np.mean(Smag_end_m, axis=0)
        Smag_start_std = np.std(Smag_start_m, axis=0)
        Smag_end_std = np.std(Smag_end_m, axis=0)
        
        #finding mean and standard deviation of elongation across T1 transitions
        p_start_m = np.ma.array(data=p_start, mask=p_start_mask)
        p_end_m = np.ma.array(data=p_end, mask=p_end_mask)
        p_start_mean = np.mean(p_start_m, axis=0)
        p_end_mean = np.mean(p_end_m, axis=0)
        p_start_std = np.std(p_start_m, axis=0)
        p_end_std = np.std(p_end_m, axis=0)
        
        #finding mean and standard deviation of velocity across T1 transitions
        avg_vel_start_m = np.ma.array(data=avg_vel_start, mask=avg_vel_start_mask)
        avg_vel_end_m = np.ma.array(data=avg_vel_end, mask=avg_vel_end_mask)
        avg_vel_start_mean = np.mean(avg_vel_start_m, axis=0)
        avg_vel_end_mean = np.mean(avg_vel_end_m, axis=0)
        avg_vel_start_std = np.std(avg_vel_start_m, axis=0)
        avg_vel_end_std = np.std(avg_vel_end_m, axis=0)
        
        #finding mean and standard deviation of acceleration across T1 transitions
        avg_acl_start_m = np.ma.array(data=avg_acl_start, mask=avg_acl_start_mask)
        avg_acl_end_m = np.ma.array(data=avg_acl_end, mask=avg_acl_end_mask)
        avg_acl_start_mean = np.mean(avg_acl_start_m, axis=0)
        avg_acl_end_mean = np.mean(avg_acl_end_m, axis=0)
        avg_acl_start_std = np.std(avg_acl_start_m, axis=0)
        avg_acl_end_std = np.std(avg_acl_end_m, axis=0)
        
        fig, ax = plt.subplots(facecolor='w')
        S_mag_meanT1 = np.concatenate((Smag_start_mean, Smag_end_mean)) 
        S_mag_stdT1 = np.concatenate((Smag_start_std, Smag_end_std))
        ax.plot(plot_times[::2], S_mag_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        ax.fill_between(plot_times[::2], S_mag_meanT1[::2]-S_mag_stdT1[::2], S_mag_meanT1[::2]+S_mag_stdT1[::2],
                        alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('Elongation of T1 cells', fontsize=22)
        plt.tight_layout()
        plt.savefig(self.result_dir + '/T1_Smag' + postfix + '.png', dpi=200)
        np.save(self.result_dir + '/T1_Smag' + postfix + '.npy', np.ma.filled(S_mag_meanT1, np.nan))
        
        fig, ax = plt.subplots( facecolor='w')
        S_mag_meanT1 = np.concatenate((Smag_start_mean, Smag_end_mean)) 
        S_mag_stdT1 = np.concatenate((Smag_start_std, Smag_end_std))
        ax.plot(plot_times[::2], S_mag_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        #ax.fill_between(plot_times[::2], S_mag_meanT1[::2]-S_mag_stdT1[::2], S_mag_meanT1[::2]+S_mag_stdT1[::2],
        #                alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('Elongation of T1 cells', fontsize=22)
        plt.tight_layout()
        
        fig, ax = plt.subplots(facecolor='w')
        p_meanT1 = np.concatenate((p_start_mean, p_end_mean)) 
        p_stdT1 = np.concatenate((p_start_std, p_end_std))
        ax.plot(plot_times[::2], p_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        ax.fill_between(plot_times[::2], p_meanT1[::2]-p_stdT1[::2], p_meanT1[::2]+p_stdT1[::2],
                        alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('shape index of T1 cells', fontsize=22)
        plt.tight_layout()
        plt.savefig(self.result_dir + '/T1_shapeindex' + postfix + '.png', dpi=200)
        np.save(self.result_dir+ '/T1_shapeindex' + postfix + '.npy', np.ma.filled(p_meanT1, np.nan))
        
        fig, ax = plt.subplots(facecolor='w')
        p_meanT1 = np.concatenate((p_start_mean, p_end_mean)) 
        p_stdT1 = np.concatenate((p_start_std, p_end_std))
        ax.plot(plot_times[::2], p_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        #ax.fill_between(plot_times[::2], p_meanT1[::2]-p_stdT1[::2], p_meanT1[::2]+p_stdT1[::2],
        #                alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('shape index of T1 cells', fontsize=22)
        plt.tight_layout()
        
        #plotting individual shape index plots
        """
        for j in range(20):
            print(50+j)
            print(np.concatenate((p_start_m[50+j], p_end_m[50+j]))[::7])
        (p_start_m)
        print(p_start_m.shape)
        """
        fig, ax = plt.subplots(facecolor='w')
        for indp, row in enumerate(p_start_m[:10]):
            p_T1 = np.concatenate((p_start_m[indp], p_end_m[indp])) 
            ax.plot(plot_times, p_T1)
        ax.legend(fontsize=22)
        ax.grid()
        #ax.set_ylim(bottom=3.0)
        ax.set_ylabel('shape index of T1 cells', fontsize=22)
        plt.tight_layout()
        
        
        #velocity of T1 cells
        fig, ax = plt.subplots(facecolor='w')
        avg_vel_meanT1 = np.concatenate((avg_vel_start_mean, avg_vel_end_mean)) 
        avg_vel_stdT1 = np.concatenate((avg_vel_start_std, avg_vel_end_std))
        ax.plot(plot_times[::2], avg_vel_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        ax.fill_between(plot_times[::2], avg_vel_meanT1[::2]-avg_vel_stdT1[::2], avg_vel_meanT1[::2]+avg_vel_stdT1[::2],
                        alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('velocity of T1 cells', fontsize=22)
        plt.tight_layout()
        plt.savefig(self.result_dir + '/T1_vel' + postfix + '.png', dpi=200)
        np.save(self.result_dir+ '/T1_vel' + postfix + '.npy', np.ma.filled(avg_vel_meanT1, np.nan))
        
        fig, ax = plt.subplots(facecolor='w')
        avg_vel_meanT1 = np.concatenate((avg_vel_start_mean, avg_vel_end_mean)) 
        avg_vel_stdT1 = np.concatenate((avg_vel_start_std, avg_vel_end_std))
        ax.plot(plot_times[::2], avg_vel_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        #ax.fill_between(plot_times[::2], avg_vel_meanT1[::2]-avg_vel_stdT1[::2], avg_vel_meanT1[::2]+avg_vel_stdT1[::2],
        #                alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('velocity of T1 cells', fontsize=22)
        plt.tight_layout()
        plt.savefig(self.result_dir + '/T1_vel2' + postfix + '.png', dpi=200)
        
        #acceleration of T1 cells
        fig, ax = plt.subplots(facecolor='w')
        avg_acl_meanT1 = np.concatenate((avg_acl_start_mean, avg_acl_end_mean)) 
        avg_acl_stdT1 = np.concatenate((avg_acl_start_std, avg_acl_end_std))
        ax.plot(plot_times[::2], avg_acl_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        ax.fill_between(plot_times[::2], avg_acl_meanT1[::2]-avg_acl_stdT1[::2], avg_acl_meanT1[::2]+avg_acl_stdT1[::2],
                        alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('acceleration of T1 cells', fontsize=22)
        plt.tight_layout()
        plt.savefig(self.result_dir + '/T1_acl' + postfix + '.png', dpi=200)
        np.save(self.result_dir+ '/T1_acl' + postfix + '.npy', np.ma.filled(avg_acl_meanT1, np.nan))
        
        fig, ax = plt.subplots(facecolor='w')
        avg_acl_meanT1 = np.concatenate((avg_acl_start_mean, avg_acl_end_mean)) 
        avg_acl_stdT1 = np.concatenate((avg_acl_start_std, avg_acl_end_std))
        ax.plot(plot_times[::2], avg_acl_meanT1[::2], 
                    color='k', linewidth=3, markersize=10)
        #ax.fill_between(plot_times[::2], avg_acl_meanT1[::2]-avg_acl_stdT1[::2], avg_acl_meanT1[::2]+avg_acl_stdT1[::2],
        #                alpha = 0.3)
        ax.legend(fontsize=22)
        ax.grid()
        ax.set_ylabel('acceleration of T1 cells', fontsize=22)
        plt.tight_layout()
        plt.savefig(self.result_dir + '/T1_acl2' + postfix + '.png', dpi=200)
                
def circle_coord(point, r, grid_X, grid_Y):
    """
    inputs:
        point: 2D array of x,y cordinates
        r:  radius of circle 
        grid_x, grid_y: 2D arrays with grid points
    outputs:
        two 1D arrays of indices of grid_x and grid_y such that the grid points
        are within a circle of radius r from the 'point'
    Note: Works only if r is less than half the domain size.
    """
    domain_dimension = np.array([np.max(grid_X), np.max(grid_Y)])
    grid_Xc = grid_X.copy()
    grid_Yc = grid_Y.copy()
    grid_Xc[(grid_X-point[0])>0.5*domain_dimension[0]] -= domain_dimension[0]
    grid_Xc[(point[0]-grid_X)>0.5*domain_dimension[0]] += domain_dimension[0]
    grid_Yc[(grid_Y-point[1])>0.5*domain_dimension[1]] -= domain_dimension[1]
    grid_Yc[(point[1]-grid_Y)>0.5*domain_dimension[1]] += domain_dimension[1]
    diff_x = np.abs(grid_Xc-point[0])
    diff_y = np.abs(grid_Yc-point[1])
    dist = np.sqrt(diff_x**2 + diff_y**2)
    
    return np.where(dist <= r)

def annular_coord(point, r_i, r_o, grid_X, grid_Y):
    """
    inputs:
        point: 2D array of x,y cordinates
        r_i: radius of inner circle 
        r_o: radius of outer circle 
        grid_x, grid_y: 2D arrays with grid points
    outputs:
        two 1D arrays of indices of grid_x and grid_y such that the grid points
        are within a circle of radius r from the 'point'
    Note: Works only if r_o is less than half the domain size.
    """
    domain_dimension = np.array([np.max(grid_X), np.max(grid_Y)])
    grid_Xc = grid_X.copy()
    grid_Yc = grid_Y.copy()
    grid_Xc[(grid_X-point[0])>0.5*domain_dimension[0]] -= domain_dimension[0]
    grid_Xc[(point[0]-grid_X)>0.5*domain_dimension[0]] += domain_dimension[0]
    grid_Yc[(grid_Y-point[1])>0.5*domain_dimension[1]] -= domain_dimension[1]
    grid_Yc[(point[1]-grid_Y)>0.5*domain_dimension[1]] += domain_dimension[1]
    diff_x = np.abs(grid_Xc-point[0])
    diff_y = np.abs(grid_Yc-point[1])
    dist = np.sqrt(diff_x**2 + diff_y**2)
    
    return np.where((dist <= r_o) & (dist>=r_i))

def force_norm(sigma_00, sigma_11, sigma_01, 
               dx=0.2, 
               pbc = True, 
               smoothen=True, 
               filter_width=3): 
    """
    #assuming dx = dy
    #It is also better to calculate force from the smoothened stress field to 
    ignore irregularities especially inside the cell
    
    inputs: 1.  stress fields: sigma_00, sigma_11, sigma_01
                2D numpy arrays of same dimension as the grid
            2.  dx: distance between two closest grid points
            3. pbc: bool, true if periodic BC
            4. smoothen: bool, if true, then gaussian filter is applied to 
                         stress
    outputs: 2d array of norm of the force as function of grid coordinates
    
    """
    
    if pbc:
        sigma_00 = gaussian_filter(sigma_00, filter_width, mode = 'wrap')
        sigma_11 = gaussian_filter(sigma_11, filter_width, mode = 'wrap')
        sigma_01 = gaussian_filter(sigma_01, filter_width, mode = 'wrap')
        
        sigma_00d0 = np.gradient(np.pad(sigma_00, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        sigma_11d1 = np.gradient(np.pad(sigma_11, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
        sigma_01d1 = np.gradient(np.pad(sigma_01, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
        sigma_10d0 = np.gradient(np.pad(sigma_01, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        
    Fx = sigma_00d0 + sigma_01d1
    Fy = sigma_10d0 + sigma_11d1
    
    return np.linalg.norm(np.array([Fx, Fy]), axis = 0)

def load_domain(input_dir):
    grid = load_grid(input_dir)
    return np.max(grid[0]), np.max(grid[1])

def load_grid(input_dir):
    """
    Parameters
    ----------
    input_dir : string
        The directory with the phasefield directory that contains numpy grid 
        files

    Returns
    -------
    grid_X : 2D numpy float array
        Grid of x positions
    grid_Y : 2D numpy float array
        Grid of y positions

    """
    grid_X = np.load(input_dir + '/global_fields_200/grid_x.npy')
    grid_Y = np.load(input_dir + '/global_fields_200/grid_y.npy')
    
    return grid_X, grid_Y

def T1_loc(input_dir, T1_transitions, start, end):
    """
    Same position returned for all time steps for a given T1 transition

    Parameters
    ----------
    input_dir : TYPE
        DESCRIPTION.
    T1_transitions : TYPE
        DESCRIPTION.
    start : TYPE
        DESCRIPTION.
    end : TYPE
        DESCRIPTION.

    Returns
    -------
    locations : TYPE
        DESCRIPTION.

    """
    times = np.load(input_dir + '/stress_fields_500/timesteps.npy')
    
    locations = np.zeros([len(T1_transitions), 2])
    
    ranks = []
    positions_raw = []
    file_pattern = input_dir + '/neo_positions_p*.csv'
    for filename in glob.glob(file_pattern):
        tmp = pd.read_csv(filename)
        positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', 
                                  'total_interaction', 'neighbours', 
                                  'confine_interaction', 'growth_rate', 
                                  'S0full', 'S1full']])
        # we also need to extract the rank to build the bridge to vtk files
        ranks.append(int(re.findall(r'\d+', filename)[-1]))
    ranks = np.array(ranks)
    
    sorted_index = np.argsort(ranks)
    positions_raw.sort(key= lambda elem: elem.iloc[0]['rank'])
    ranks = ranks[sorted_index]
    print('times is of length', len(times))
    x0 = np.zeros([len(times), len(ranks)])
    x1 = np.zeros([len(times), len(ranks)])
    for rank in ranks:
        x0[:, rank] = positions_raw[rank]['x0']
        x1[:, rank] = positions_raw[rank]['x1']
    
    domain = load_domain(input_dir)
    for ind_t, t1 in enumerate(T1_transitions):
        mid_time = (end[ind_t] + start[ind_t])/2
        nearest_time_ind = np.argmin(np.abs(times-mid_time))

        cell_x = np.zeros(len(t1))
        cell_y = np.zeros(len(t1))
        
        for ind_c, cell in enumerate(t1):
            cell_x[ind_c] = x0[nearest_time_ind, cell]
            cell_y[ind_c] = x1[nearest_time_ind, cell]
        
        #Correcting for periodic domain
        if ((np.max(cell_x) - np.min(cell_x)) > (domain[0]/2)):
            cell_x[cell_x<domain[0]/2] += domain[0]
        if ((np.max(cell_y) - np.min(cell_y)) > (domain[1]/2)):
            cell_y[cell_y<domain[1]/2] += domain[1]
            
        avg_x = np.mean(cell_x)
        avg_y = np.mean(cell_y)
        
        if avg_x > domain[0]:
            avg_x -= domain[0]
        if avg_y > domain[1]:
            avg_y -= domain[1]
        
        locations[ind_t] = avg_x, avg_y
    return locations



def load_pos(input_dir, smoothen=False):
    """
    Parameters
    ----------
    input_dir : string
        The directory which contains the position directory with neopositions
        csv files
    smoothen : bool, optional
        If true, then each cells velocity is further smoothened so that cells have some velocity
        influence based on what their velocity was in the previous timestep. This can be helpful 
        when a large number of cells appear to be stationary during a given timestep.

    Returns
    -------
    positions_raw : list of panda databases
        Each database corresponds to a given cell. It stores information for 
        every *th timestep about different cell properties like position of 
        center of mass, velocity of center of mass, ..
    
    ranks : int array
        An array of ranks or cell indices
    """
    ranks = []
    positions_raw = []
    file_pattern = input_dir + '/neo_positions_p*.csv'
    for filename in glob.glob(file_pattern):
        tmp = pd.read_csv(filename)
        positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', 
                                  'total_interaction', 'neighbours', 
                                  'confine_interaction', 'growth_rate', 
                                  'S0full', 'S1full']])
        # we also need to extract the rank to build the bridge to vtk files
        ranks.append(int(re.findall(r'\d+', filename)[-1]))
    ranks = np.array(ranks)
    
    sorted_index = np.argsort(ranks)
    positions_raw.sort(key= lambda elem: elem.iloc[0]['rank'])
    ranks = ranks[sorted_index]
    """
    if smoothen:
        for rank in ranks:
            np.nan_to_num(positions_raw[rank]['v0'], nan=0.0)
            np.nan_to_num(positions_raw[rank]['v1'], nan=0.0)
            positions_raw[rank] = uniform_filter1d(positions_raw[rank], size=5)
            positions_raw[rank] = uniform_filter1d(positions_raw[rank], size=5)
    """
    return positions_raw, ranks


def load_property(input_dir, prop='velocity', smoothen=False):
    """

    Parameters
    ----------
    input_dir : string
        input directory
    prop : string, optional
        This is the name of the property for which you would like the field of.
        The possible options are 
        1. 'velocity'
        2. 'normalised nematic' : The two independeny shape tensor components 
            are normalised such that the magnitude of S is 1
        3. 'nematic' : Returns the two independent shape tensor components
        4. 'velocity angle' : Returns the velocity angle field in [pi, pi]
        5. 'nematic angle' : Returns the nematic orientation field in 
            [pi/2, pi/2]
        The default is 'velocity'.
    smoothen: bool, optional
        If true, then each cells velocity is further smoothened so that cells have some velocity
        influence based on what their velocity was in the previous timestep. This can be helpful 
        when a large number of cells appear to be stationary during a given timestep.

    Returns
    -------
    one or more 2D array
        First axis is time and the second axis is the property of the cell and 
        is indexed by rank

    """
    positions, ranks = load_pos(input_dir, smoothen)
    times = np.load(input_dir + '/stress_fields_500/timesteps.npy')
    if prop == 'position':
        x0 = np.zeros([len(times), len(ranks)])
        x1 = np.zeros([len(times), len(ranks)])
        for rank in ranks:
            x0[:, rank] = positions[rank]['x0']
            x1[:, rank] = positions[rank]['x1']
        return x0, x1
    
    if prop == 'velocity':
        v0 = np.zeros([len(times), len(ranks)])
        v1 = np.zeros([len(times), len(ranks)])
        for rank in ranks:
            v0[:, rank] = positions[rank]['v0']
            v1[:, rank] = positions[rank]['v1']
        return v0, v1
    
    if prop == 'acceleration':
        v0, v1 = load_property(input_dir, prop='velocity')
        a0 = np.diff(v0, axis=0)/(times[5]-times[4])
        a1 = np.diff(v1, axis=0)/(times[5]-times[4])
        return a0, a1
    
    if prop == 'normalised nematic':
        S0 = np.zeros([len(times), len(ranks)])
        S1 = np.zeros([len(times), len(ranks)])
        for rank in ranks:
            S0[:, rank] = positions[rank]['S0']
            S1[:, rank] = positions[rank]['S1']
        return S0, S1
        
    if prop == 'nematic':
        S0 = np.zeros([len(times), len(ranks)])
        S1 = np.zeros([len(times), len(ranks)])
        for rank in ranks:
            S0[:, rank] = positions[rank]['S0full']
            S1[:, rank] = positions[rank]['S1full']
        return S0, S1
    
    if prop == 'velocity angle':
        v0 = np.zeros([len(times), len(ranks)])
        v1 = np.zeros([len(times), len(ranks)])
        for rank in ranks:
            v0[:, rank] = positions[rank]['v0']
            v1[:, rank] = positions[rank]['v1']
        vel_orient = np.arctan2(v1, v0)
        return vel_orient
    
    if prop == 'nematic angle':
        S0 = np.zeros([len(times), len(ranks)])
        S1 = np.zeros([len(times), len(ranks)])
        for rank in ranks:
            S0[:, rank] = positions[rank]['S0']
            S1[:, rank] = positions[rank]['S1']
        nem_orient = np.multiply(np.sign(S1),
                                 (np.arctan(S0/np.abs(S1))
                                  /2.0 + np.pi/4.0))
        return nem_orient
    
    if prop == 'radius':
        rad = np.zeros([len(times), len(ranks)])
        for rank in ranks:
            rad[:, rank] = positions[rank]['r']
        return rad
    
    if prop == 'area':
        return np.pi*(load_property(input_dir, prop='radius')**2)

    
def load_ranks(input_dir):
    ranks = []
    positions_raw = []
    file_pattern = input_dir + '/neo_positions_p*.csv'
    for filename in glob.glob(file_pattern):
        tmp = pd.read_csv(filename)
        positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', 
                                  'total_interaction', 'neighbours', 
                                  'confine_interaction', 'growth_rate', 
                                  'S0full', 'S1full']])
        # we also need to extract the rank to build the bridge to vtk files
        ranks.append(int(re.findall(r'\d+', filename)[-1]))
    ranks = np.array(ranks)
    return ranks

    
a = Quantity(sys.argv[1], delt=0.125)

a.force_T1_a(normalise=False, plot_std=False)
a.force_prepost_T1(gaps=np.array([10, 20, 30, 40, 50]), normalise=False, equalise_T1num=False)
a.force_prepost_T1(gaps=np.array([10, 20, 30, 40, 50]), normalise=False, equalise_T1num=True)