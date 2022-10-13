import numpy as np
import csv
import glob
import re
import pandas as pd

def get_neighbour_relations(input_dir, start_ind, end_ind, stride):
    """
    Parameters
    ----------
    input_dir: string
        input directory which has a neighbours folder which contains all the neighbour*.dat files
    ranks: list of cells or ranks
    start_ind:
        Time index from which we start considering neighbour relations. Should be atleast as 
        big as the number of lines that are skipped
    end_ind:
        The last time index upto which neighbour relations are computed
    stride:
        The gap in the number of indices which are skipped. Typically the same as the stride 
        where the neopositions files are saved.
    Output
    ------
    neighbours_arr: 3d numpy array
        The first index is for the time. Then for each cell, there is a second dimension and 
        if that cell is a neighbour of another cell then the corresponding index of the third
        dimension is 1 else 0.
    """
    file_pattern = input_dir + '/neighbours_p*.dat'
    ranks = load_ranks(input_dir)
    neighbours_arr = np.zeros([len(np.arange(start_ind, end_ind, stride)), len(ranks), len(ranks)], dtype=int)
    for indc, cell in enumerate(ranks):
        print(indc)
        filename = input_dir + '/neighbours_p' + str(cell) + '.dat'
        with open(filename) as file:
            lines = file.readlines()
            for indl, line in enumerate(lines):
                lines[indl] = list(map(int,lines[indl].split(" ")[:-1]))
            lines = lines[start_ind:end_ind:stride]
            for indl, line in enumerate(lines):
                neighbours_arr[indl, cell][line] = 1
    return neighbours_arr

def load_ranks(input_dir):
    ranks = []
    file_pattern = input_dir + '/neo_positions_p*.csv'
    for filename in glob.glob(file_pattern):
        ranks.append(int(re.findall(r'\d+', filename)[-1]))
    ranks = np.array(ranks)
    sorted_index = np.argsort(ranks)
    ranks = ranks[sorted_index]
    return ranks


def stress(phi_all, neighbours, ranks, mu_all=np.zeros(10), Ca=20e-2, In=10e-2, epsilon=0.15, pbc_=True, dx_ = 0.5, a = 1.5, method=1):
    #Considers mu in terms of rescaled phi. Just add a factor of 2
    #returns the normal and shear stress fields calculated 
    #across the planes of the coordinate axis
    sigma_00 = np.zeros(phi_all[0].shape)
    sigma_11 = np.zeros(phi_all[0].shape)
    sigma_01 = np.zeros(phi_all[0].shape)
    
    def grad_phi(phi_n, dx=0.5, pbc=True):
        if pbc:
            grad_phix = np.gradient(np.pad(phi_n, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
            grad_phiy = np.gradient(np.pad(phi_n, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
        else:
            grad_phix = np.gradient(phi_n, dx, axis=0)
            grad_phiy = np.gradient(phi_n, dx, axis=1)
        
        return np.array([grad_phix, grad_phiy])

    def g(phi_n):
        #return 3/(2*(1-phi_n)*(1+phi_n))
        return 1

    def W(phi_n):
        return ((phi_n**2 - 1)**2)/4 

    def grad_W(phi_n):
        return phi_n**3 - phi_n
    
    def w_int(phi_k, a):
        return 1 - ((a+1)*(phi_k-1)**2/4) + (a*(phi_k-1)**4/16)

    def w_int_der(phi_k, a):
        return -(2*(a+1)*(phi_k-1)/4) + (4*a*(phi_k-1)**3/16)
    
    def laplace_phi2(gradphi_n, dx=dx_, pbc = pbc_):
        if pbc:
            grad_phixx = np.gradient(np.pad(gradphi_n[0], 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
            grad_phiyy = np.gradient(np.pad(gradphi_n[1], 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
        else:
            grad_phixx = np.gradient(gradphi_n[0], dx, axis=0)
            grad_phiyy = np.gradient(gradphi_n[1], dx, axis=1)
        return grad_phixx + grad_phiyy
    
    def f_mu_IN_nk(phi_n, phi_k, In=10e-2, a=1.5):
        #returns free energy and chemical potential of interaction energy
        w = w_int(phi_k, a)
        w_der = w_int_der(phi_n, a)
        f = (1/In)*((phi_n+1)/2)*(w)
        mu = (1/In)*((phi_n/2*w) + w_der*(phi_k+1)/2)
        return f, mu

    def f_mu_CH(phi_n, gradphi_n, epsilon, Ca):
        #returns free energy and chemical potential of Cahn Hilliard energy
        f = (1/Ca)*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))
        mu = (1/Ca)*((grad_W(phi_n)/epsilon) - (epsilon*(laplace_phi2(gradphi_n))))
        return f, mu
    
    def f_CH(phi_n, gradphi_n, epsilon, Ca):
        return (1/Ca)*(g(phi_n))*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))
    
    def f_IN_nk(phi_n, phi_k, In=10e-2, a=1.5):
        #returns free energy and chemical potential of interaction energy
        w = w_int(phi_k, a)
        f = (1/In)*((phi_n+1)/2)*(w)
        return f
    
    for cell_n in ranks:
        grad_phi_ = grad_phi(phi_all[cell_n], dx=dx_, pbc=pbc_)
        
        isotropic = np.zeros(phi_all[0].shape)
        
        if method == 1:
            f, mu = f_mu_CH(phi_all[cell_n], grad_phi_, epsilon, Ca)
            isotropic = f
            isotropic -= 2*mu*(phi_all[cell_n]+1)/2
        
        #for cell_k in [cell_k for cell_k in rank if cell_k!=cell_n]:#considers all cells 
        for cell_k in np.where(neighbours[cell_n]==1)[0]: #considers only neighbours
            f, mu = f_mu_IN_nk(phi_all[cell_n], phi_all[cell_k],epsilon, Ca)
            isotropic += f
            if method == 1:
                isotropic -= 2*mu*(phi_all[cell_n]+1)/2
        
        if method == 2: #doesn't calculate mu but uses mu_cell
            isotropic -= 2*mu_all[cell_n]*(phi_all[cell_n]+1)/2
            
        sigma_00 += isotropic - grad_phi_[0]**2
        sigma_11 += isotropic - grad_phi_[1]**2
        sigma_01 += -grad_phi_[0]*grad_phi_[1]
    return sigma_00, sigma_11, sigma_01