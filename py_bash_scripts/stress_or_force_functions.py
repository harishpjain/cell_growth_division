import numpy as np
dx=dx

def grad_phi(phi_n, dx=dx, pbc=True):
    if pbc:
        grad_phix = np.gradient(np.pad(phi_n, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        grad_phiy = np.gradient(np.pad(phi_n, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
    else:
        grad_phix = np.gradient(phi_n, dx, axis=0)
        grad_phiy = np.gradient(phi_n, dx, axis=1)
    
    return np.array([grad_phix, grad_phiy])

def curvature(phi_n, dx=dx, pbc=True):
    if pbc:
        grad_phix = np.gradient(np.pad((phi_n+1.0)/2.0, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        grad_phiy = np.gradient(np.pad((phi_n+1.0)/2.0, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]

        grad_norm = np.sqrt(grad_phix**2 + grad_phiy**2)

        grad_phixx = np.gradient(np.pad(grad_phix, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        grad_phiyy = np.gradient(np.pad(grad_phiy, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]

        laplace_phi = grad_phixx+grad_phiyy

        grad_of_grad_norm_x = np.gradient(np.pad(grad_norm, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        grad_of_grad_norm_y = np.gradient(np.pad(grad_norm, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]

        grad_phi_prod_grad_of_grad_norm = grad_phix*grad_of_grad_norm_x + grad_phiy*grad_of_grad_norm_y
        curv_ = (laplace_phi - grad_phi_prod_grad_of_grad_norm/grad_norm)/grad_norm
        curv_[grad_norm < 0.001] = 0.0
        return curv_
        

def g(phi_n):
    #return 3/(2*(1-phi_n)*(1+phi_n))
    return 1

def W(phi_n):
    return ((phi_n**2 - 1)**2)/4 

def grad_W(phi_n):
    return phi_n**3 - phi_n

def w_int(phi_k, a=1.5):
    return 1 - ((a+1)*(phi_k-1)**2/4) + (a*(phi_k-1)**4/16)

def w_int_der(phi_k, a=1.5):
    return -(2*(a+1)*(phi_k-1)/4) + (4*a*(phi_k-1)**3/16)

def laplace_phi2(gradphi_n, dx=dx, pbc = True):
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
    mu = (1/In)*((0.5*w) + w_der*(phi_k+1)/2)
    return f, mu

def f_mu_CH(phi_n, gradphi_n, epsilon=0.15, Ca=20e-2, dx=dx):
    #returns free energy and chemical potential of Cahn Hilliard energy
    f = (1/Ca)*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))
    mu = (1/Ca)*((grad_W(phi_n)/epsilon) - (epsilon*(laplace_phi2(gradphi_n, dx=dx))))
    return f, mu

def f_CH(phi_n, gradphi_n, epsilon=0.15, Ca=20e-2):
    #return (1/Ca)*(g(phi_n))*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))
    return (1/Ca)*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))

def f_IN_nk(phi_n, phi_k, In=10e-2, a=1.5):
    #returns free energy and chemical potential of interaction energy
    w = w_int(phi_k, a)
    f = (1/In)*((phi_n+1)/2)*(w)
    return f

def f_mu_IN_nk_neo(phi_n, phi_k, In=10e-2, a_rep=1.0, a_adh=1.5):
    #returns interaction free energy density 
    f_rep = 0.5*a_rep*((phi_n+1)**2)*((phi_k+1)**2)
    f_adh = 0.5*a_adh*((phi_k**2 - 1)**2)*((phi_n**2 - 1)**2)
    mu_rep = a_rep*(phi_n+1)*((phi_k+1)**2)
    mu_adh = 2.0*a_adh*(phi_n**3 - phi_n)*((phi_k**2 - 1)**2)
    return (f_adh + f_rep)/In, (mu_adh + mu_rep)/In