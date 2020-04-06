#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main function that is called to initialize and run phase-field dynamics
"""

from __future__ import print_function
import fipy as fp
import os
#from fipy.solvers.pysparse import LinearLUSolver as Solver
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from utils.input_parse import *
from utils.graphics import *
import utils.free_energy as f_en
import utils.sim_tools as sim_tools
import timeit
import subprocess

def return_noise(mesh,time_step,noise_strength):
    """
        return_noise returns Gaussian noise of a given variance = noise_strength
        for a particular mesh & time_step
    """
    p = fp.GaussianNoiseVariable(mesh=mesh,variance=noise_strength/(time_step*mesh.cellVolumes))    
    p = p - np.mean(p.value);
    return(p);
    
def run_CH(args):
    """
    Function takes in path to input params, output_folder, and optionally params file that
    are passed while running code. With these parameters, this function initializes and runs
    phase-field simulations while writing output to files.

    **Input variables**

    -   args.i = Path to input params files (required)
    -   args.o = path to output folder (required)
    -   args.p = path to parameter file (optional)
    -   args.pN     =   Nth parameter to use from input (optional)

    """

    # In[4]:
    seed = np.random.randint(1e8);
    input_parameters = input_parse(args.i);

    if args.p is not None:
        params = input_parse(args.p,params_flag=True)
        par_name = str(list(params.keys())[0])
        par_values = params[par_name];
        if args.pN is not None:
            par_values = [par_values[int(args.pN)-1]];
    else:
        par_name = 'nx';
        par_values = [(int(input_parameters['nx']))];

    for par in par_values:
        start = timeit.default_timer()

        input_parameters[par_name] = par;
        nx = int(input_parameters['nx'])
        dx = input_parameters['dx']
        c_alpha = input_parameters['c_alpha'];
        c_beta = input_parameters['c_beta'];
        kappa = input_parameters['kappa']
        M_protein = input_parameters['M_protein']
        M_rna = input_parameters['M_rna']
        dimension = int(input_parameters['dimension'])
        plot_flag = bool(input_parameters['plot_flag'])
        rho_s = input_parameters['rho_s']
        rho_r = input_parameters['rho_r']
        chi = input_parameters['chi']
        changing_chi = int(input_parameters['changing_chi']);
        fh = int(input_parameters['fh']);

        """ Unpack rates from the parameters """
        k_production = input_parameters['k_production']
        k_degradation = input_parameters['k_degradation'];

        """ Define size of initial nucleus """
        nucleus_size = int(input_parameters['nucleus_size']);

        """ Initial RNA & protein concentrations """
        phi_p_0 = input_parameters['phi_p_0'];
        phi_r_0 = input_parameters['phi_r_0'];



        if 'a' in input_parameters.keys():
            a = float(input_parameters['a']);
        else:
            a=0.0;
        if 'b' in input_parameters.keys():
            b = float(input_parameters['b']);
        else:
            b=0.0;
        if 'c' in input_parameters.keys():
            c = float(input_parameters['c']);
        else:
            c=0.0;
            
        if 'noise_strength' in input_parameters.keys():
            noise_strength = float(input_parameters['noise_strength']);
        else:
            noise_strength=0.0;

        if 'm1' in input_parameters.keys():
            m1 = float(input_parameters['m1']);
        else:
            m1=1.0;
 
        if 'kp_noise' in input_parameters.keys():
            kp_noise = float(input_parameters['kp_noise']);
        else:
            kp_noise=0.0;

        if 'chi_ps' in input_parameters.keys():
            chi_ps = float(input_parameters['chi_ps']);
        else:
            chi_ps=None;
            
        if 'chi_rs' in input_parameters.keys():
            chi_rs = float(input_parameters['chi_rs']);
        else:
            chi_rs=0.0;
            
        if 'r' in input_parameters.keys():
            r = float(input_parameters['r']);
        else:
            r=1.0;


            
        if 'seed' in input_parameters.keys():
            seed=int(input_parameters['seed'])

        fp.numerix.random.seed(seed);


        """
        Set-up the appropriate choice of free-energy
            fh is a flag for employing Flory-Huggins instead of double-well
            changing_chi ==2 uses the gaussian form & 1 == uses double-well LG expression
            changing_chi ==0 is not changing_chi and there for backwards compatability
            rho_s/rho_r is height of double-well potential for protein/RNA respectively
            kappa is surface tension parameter for protein
            chi is value of pairwise interaction
            Y is value of landau-ginzburg like three-way interaction
            mu_r chooses whether you use D-R (mu_r=0) or chemical potential fo RNA (mu_r=1)
            a,ratio, and p characterize the gaussian form of chi
        """

        if not fh:

            if changing_chi==2:
                FE = f_en.free_energy_changing_chi(c_alpha=c_alpha,c_beta=c_beta,rho_s=rho_s,rho_r=rho_r,chi=chi,kappa=kappa,a=input_parameters['a'],ratio=input_parameters['ratio'],p=input_parameters['p'])
            elif changing_chi==1:
                FE = f_en.free_energy_changing_chi_LG(c_alpha=c_alpha,c_beta=c_beta,rho_s=rho_s,rho_r=rho_r ,chi=chi,kappa=kappa,a=a,b=b,c=c)
            else:
                FE = f_en.free_energy(c_alpha=c_alpha,c_beta=c_beta,rho_s=rho_s,chi=chi,kappa=kappa)

        else:

            if changing_chi==2:
                FE = f_en.free_energy_FH_changing_chi(c_alpha=c_alpha,c_beta=c_beta,rho_s=rho_s,chi=chi,kappa=kappa,a=input_parameters['a'],ratio=input_parameters['ratio'],p=input_parameters['p']);
            elif changing_chi==1:
                FE = f_en.free_energy_FH_changing_chi_LG(c_alpha=c_alpha,c_beta=c_beta,rho_s=rho_s,rho_r=rho_r,chi=chi,kappa=kappa,a=a,b=b,c=c,chi_ps=chi_ps,chi_rs=chi_rs,r=r)
            else:
                FE = f_en.free_energy_FH(c_alpha=c_alpha,c_beta=c_beta,rho_s=rho_s,chi=chi,kappa=kappa);

        """
        Define the parameters that dictate reaction kinetics
            if multiplier is specified, so must  t_change.
            Then after t_change has passed, the simulation will multiply k_production by multiplier
            threshold will ensure production only at phi_p>=threshold values
        """


        if 'multiplier' in input_parameters.keys():
            rates = f_en.RNA_reactions(k_production=k_production,k_degradation=k_degradation,threshold=input_parameters['threshold'],t_change=input_parameters['t_change'],multiplier=input_parameters['multiplier'],m1=m1,kp_noise=kp_noise);
        else:
            rates = f_en.RNA_reactions(k_production=k_production,k_degradation=k_degradation,threshold=input_parameters['threshold'],m1=m1,kp_noise=kp_noise);



        if dimension==2:
            if int(input_parameters['circ_flag']):
                mesh = sim_tools.create_circular_mesh(radius=float(nx)*dx/2,cellSize=dx*1.5)
            else:
                mesh = fp.Grid2D(nx=nx, ny=nx, dx=dx, dy=dx)
                mesh = mesh-float(nx)*dx*0.5
        elif dimension==3:
            mesh = fp.Grid3D(nx=nx, ny=nx,nz=nx, dx=dx, dy=dx,dz=dx)
            mesh = mesh-float(nx)*dx*0.5

        phi_p = fp.CellVariable(mesh=mesh, name=r'$\phi_{prot}$', hasOld=True,value = phi_p_0)
        phi_r = fp.CellVariable(mesh=mesh, name=r'$\phi_{RNA}$', hasOld=True,value = phi_r_0)
        phi_p[:] =fp.GaussianNoiseVariable(mesh=mesh,mean=phi_p_0,variance=0.1*phi_p_0).value
        phi_p[phi_p<phi_p_0*0.9] = phi_p_0*0.9;
        phi_p[phi_p>phi_p_0*1.1] = phi_p_0*1.1;
        print(min(phi_p),max(phi_p),np.mean(phi_p))


        phi_r[:] =fp.GaussianNoiseVariable(mesh=mesh,mean=phi_r_0,variance=0.1*phi_p_0).value
        phi_r[phi_r<phi_r_0*0.9] = phi_r_0*0.9;
        phi_r[phi_r>phi_r_0*1.1] = phi_r_0*1.1;

        print(min(phi_r),max(phi_r),np.mean(phi_r))

        # # We nucleate a high dense region at the center of the grid
        # # array of sample $\phi_{a}$-values:

        # In[5]:


        sim_tools.nucleate_seed(mesh,phi_p,phia_value=0.9*(c_beta),nucleus_size=nucleus_size,dimension=dimension)



        # ## Define relevant equations for this system

        # In[6]:
        t = fp.Variable(0.0)
        dt = input_parameters['dt'];
        dt_max = input_parameters['dt_max'];
        dt_min = input_parameters['dt_min'];
        tolerance = input_parameters['tolerance'];
        total_steps = int(input_parameters['total_steps']);
        checkpoint = int(input_parameters['checkpoint']);
        if 'text_log' in input_parameters.keys():
            text_log = int(input_parameters['text_log'])
        else:
            text_log = checkpoint;
        duration = input_parameters['duration'];
        time_step = fp.Variable(dt)
        print(time_step)
        
        eqn0 = fp.TransientTerm(coeff=1.,var=phi_p) == fp.DiffusionTerm(coeff=M_protein*FE.dmu_p_dphi_p(phi_p,phi_r),var=phi_p) + fp.DiffusionTerm(coeff=M_protein*FE.dmu_p_dphi_r(phi_p,phi_r),var=phi_r) - fp.DiffusionTerm(coeff=(M_protein,FE.kappa),var=phi_p) + M_protein*return_noise(mesh,time_step,noise_strength);

        if not int(input_parameters['mu_r']):
            eqn1 = fp.TransientTerm(coeff=1.,var=phi_r) == fp.DiffusionTerm(coeff=M_rna,var=phi_r) + rates.production(phi_p,phi_r,t) - rates.degradation(phi_p,phi_r) + M_rna*return_noise(mesh,time_step,noise_strength);
        else:
            eqn1 = fp.TransientTerm(coeff=1.,var=phi_r) == fp.DiffusionTerm(coeff=M_rna*FE.dmu_p_dphi_r(phi_p,phi_r),var=phi_p) + fp.DiffusionTerm(coeff=M_rna*FE.dmu_r_dphi_r(phi_p,phi_r),var=phi_r) + rates.production(phi_p,phi_r,t) - rates.degradation(phi_p,phi_r) +M_rna* return_noise(mesh,time_step,noise_strength);




        # ## Generate output directory strcuctures

        # In[9]:

        """
        Generates overall output_structure
                Folder = 'Output/' +  output_folder_specified + 'simulation_params' + 'seed'

            stats file will contain
                step number, time, dt, Radius, Pmin, Pmax, Rmin,Rmax,Pavg,Ravg, f,t_sim
        """
        output_directory = 'Output/' + args.o +'/'
        traj_dir = 'L_' +  str(round(nx*dx,2)) + '_phi_p_0_'+str(phi_p_0) + '_phi_r_0_'+ str(phi_r_0) + '_chiPR_'+ str(chi) + '_k_production_'+ str(k_production) +'_k_degradation_'+ str(k_degradation) + '_d_' + str(dimension);
        traj_dir = traj_dir + '_a_' + str(a) + '_b_' + str(b)+ '_c_' + str(c);
        traj_dir = traj_dir + '_rhos_' + str(rho_s) + '_rhor_' + str(rho_r)+ '_kappa_' + str(kappa);
        traj_dir = traj_dir + '_ca_' + str(c_alpha) + '_cb_' + str(c_beta);
        traj_dir = traj_dir + '_param_' + str(par_name) + '_' + str(par);
        rand_dir_id = '/' + str(seed) + '/'
        output_dir = output_directory + traj_dir + rand_dir_id;
        os.makedirs(output_dir);
        os.makedirs(output_dir+ 'Images/');
        os.makedirs(output_dir+ 'Mesh/');
        print(output_dir)

        with open(output_dir+ "/stats.txt", 'w+') as stats:
            stats.write("\t".join(["step", "t", "dt",'r', "Pmin", "Pmax", 'Rmin','Rmax',"Pavg","Ravg", "f","t_sim"]) + "\n")

        write_input_params(output_dir + '/input_params.txt',input_parameters)

        # In[10]:


        # ## Solve the Equation

        # To solve the equation a simple time stepping scheme is used which is decreased or increased based on whether the residual decreases or increases. A time step is recalculated if the required tolerance is not reached. In addition, the time step is kept under 1 unit. The data is saved out every 10 steps.

        elapsed = 0.0
        steps = 0


        phi_p.updateOld()
        phi_r.updateOld()



        while (elapsed <= duration) and (steps <= total_steps) and (dt>dt_min):

            res1 = eqn1.sweep(dt=dt)
            res0 = eqn0.sweep(dt=dt)

            if max(res0,res1) > tolerance:

                # anything in this loop will only be executed every $checkpoint steps
                if (steps % checkpoint == 0):
                    if (changing_chi==1):
                        fp.TSVViewer(vars=[phi_p,phi_r,FE.chi_eff(phi_r,phi_p)]).plot(filename=output_dir +"Mesh/mesh_{step}.txt".format(step=steps))
                    elif (changing_chi==2):
                        fp.TSVViewer(vars=[phi_p,phi_r,FE.chi_eff(phi_r)]).plot(filename=output_dir +"Mesh/mesh_{step}.txt".format(step=steps))

                    else:
                        fp.TSVViewer(vars=[phi_p,phi_r]).plot(filename=output_dir +"Mesh/mesh_{step}.txt".format(step=steps))
                        
                    if (dimension==2) and (plot_flag):

                        fig, ax  =plt.subplots()
                        cs = ax.tricontourf(mesh.x.value,mesh.y.value,phi_p.value,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(0,1.15*c_beta,256))
                        fig.colorbar(cs)
                        ax.set_title(phi_p.name)
                        fig.savefig(fname=output_dir +'Images/P_step_{step}.png'.format(step=steps),dpi=300,format='png')
    #                    fp.MatplotlibViewer(vars=phi_p,levels=np.linspace(0,1.0,10),cmap=plt.cm.get_cmap('Blues')).plot(filename=output_dir +'Images/A_step_{step}.png'.format(step=steps))
                        if input_parameters['svg_flag']:
                            for c in cs.collections:
                                c.set_edgecolor("face")
                            fig.savefig(fname=output_dir +'Images/P_step_{step}.svg'.format(step=steps),dpi=600,format='svg')
                        plt.close()


    #                    fp.MatplotlibViewer(vars=phi_r,datamin=0,datamax=0.35,cmap=plt.cm.get_cmap('PuRd')).plot(filename=output_dir +'Images/B_step_{step}.png'.format(step=steps))

                        fig, ax  =plt.subplots()
                        cs = ax.tricontourf(mesh.x.value,mesh.y.value,phi_r.value,cmap=plt.cm.get_cmap("PuRd"),levels=np.linspace(0,2.5e-1+1.15*k_production*c_beta/(k_degradation+1e-9),256))
                        fig.colorbar(cs)
                        ax.set_title(phi_r.name)
                        fig.savefig(fname=output_dir +'Images/R_step_{step}.png'.format(step=steps),dpi=300,format='png')
                        if input_parameters['svg_flag']:
                            for c in cs.collections:
                                c.set_edgecolor("face")
                            fig.savefig(fname=output_dir +'Images/R_step_{step}.svg'.format(step=steps),dpi=600,format='svg')
                        plt.close()


                        if (changing_chi):
                            fig, ax  =plt.subplots()
                            cs = ax.tricontourf(mesh.x,mesh.y,FE.chi_eff(phi_r,phi_p).value,cmap=plt.cm.get_cmap("RdYlGn"),levels=np.linspace(-FE.chi-1e-3,FE.chi+1e-3,256))
                            fig.colorbar(cs)
                            ax.set_title('$ \chi $')
                            fig.savefig(fname=output_dir +'Images/chi_step_{step}.png'.format(step=steps),dpi=300,format='png')
                            if input_parameters['svg_flag']:
                                for c in cs.collections:
                                    c.set_edgecolor("face")
                                fig.savefig(fname=output_dir +'Images/chi_step_{step}.svg'.format(step=steps),dpi=600,format='svg')
                            plt.close()
                
                if (steps % text_log ==0):                        
                    with open(output_dir+ "/stats.txt", 'a') as stats:
                        stats.write("\t".join([str(it) for it in [steps, t.value, dt, sim_tools.get_radius(phi_p,mesh,dimension=dimension,threshold=0.5*(c_alpha+c_beta)), min(phi_p), max(phi_p), min(phi_r), max(phi_r),np.mean(phi_p*mesh.cellVolumes),np.mean(phi_r*mesh.cellVolumes),
                                                              np.sum((FE.f(phi_p,phi_r)*mesh.cellVolumes).value), str(round((timeit.default_timer()-start),2))]]) + "\n")

                steps += 1
                elapsed += dt
                t.value = t.value +dt
                
                dt *= 1.1
                dt = min(dt, dt_max)
                time_step.value = dt;
                phi_p.updateOld()
                phi_r.updateOld()

            else:
                dt *= 0.8
                time_step.value = dt;
                phi_p[:] = phi_p.old
                phi_r[:] = phi_r.old


        if dimension==2 and (plot_flag):
            save_movie(output_dir +'Images/',duration=0.25)
            if input_parameters['svg_flag']:
                bash_cmd = 'rm '+ output_dir +'Images/*.png'
                res = subprocess.check_output(['bash','-c',bash_cmd])

        elif dimension==3 and (plot_flag):
            generate_images_3D(output_dir,label_idx=3,N=nx,colormap="Blues",vmin=0.0,vmax=1.0,opacity=0.2)
            generate_images_3D(output_dir,label_idx=4,N=nx,colormap="PuRd",vmin=0.0,vmax=0.35,opacity=0.2)
            save_movie(output_dir +'Images/',duration=0.25)
    # In[11]:


#    print(res0,res1)


if __name__ == "__main__":
    """
        Function is called when python code is run on command line and calls run_CH
        to initialize the simulation
    """
    parser = argparse.ArgumentParser(description='Take output filename to run CH simulations')
    parser.add_argument('--i',help="Name of input params", required = True);
    parser.add_argument('--p',help="Name of parameter file", required = False);
    parser.add_argument('--pN',help="Parameter number from file (indexed from 1)", required = False);

    parser.add_argument('--o',help="Name of output folder", required = True);
    args = parser.parse_args();

    run_CH(args);
