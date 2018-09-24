#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# adapted by An Ghysels (August 2012)
#


import numpy as np

from MCState import MCState
from transitions import Transitions, RadTransitions
from log import Logger, write_average_from_pic


def do_mc_cycles(MC,logger):
    # MONTE CARLO OPTIMIZATION    # TODO this function can become function of MCState object
    print "\n MC-move log-like acc(v) acc(w)"

    MC.init_log_like()
    logger.log(0,MC)


    for imc in range(MC.nmc):  # loop over Monte Carlo moves

        if MC.lmax > 0:
            MC.mcmove_diffusion_radial()

        else:
            choice = np.random.rand()
            if choice < 0.5 and MC.model.ncosF != 1 and MC.dv > 0.:  # do not update if only a flat basis function
                # potential move
                MC.mcmove_potential()
            elif MC.dw > 0.:
                # diffusion move
                MC.mcmove_diffusion()
            if MC.move_timezero and MC.dtimezero > 0:
                # time offset move
                MC.mcmove_timezero()

        # print
        printfreq = 100  # TODO make this optional
        MC.print_intermediate(imc,printfreq)
        logger.log(imc+1,MC)
 
        # update
        MC.update_movewidth(imc)
        MC.update_temp(imc)


def find_parameters(filenames,pbc,model,
      dv,dw,dwrad,D0,dtimezero,temp,temp_end,nmc,nmc_update,seed,outfile, ncosF,ncosD,ncosDrad,
      move_timezero,initfile,k,
      lmax,reduction,pull): 
    print "python program to extract diffusion coefficient and free energy from transition counts"
    print "copyright: Gerhard Hummer (NIH, July 2012)"
    print "adapted by An Ghysels (August 2012)\n"

    if seed is not None:
        np.random.seed(seed)

    # start Monte Carlo object
    MC = MCState(pbc,lmax)
    # settings
    MC.set_MC_params(dv,dw,dwrad,D0,dtimezero,temp,nmc,nmc_update,move_timezero,k,temp_end=temp_end,)
    #MC.print_MC_params()

    # INPUT and INITIALIZATION model/MC
    if MC.do_radial:
        data = RadTransitions(filenames)
    else:
        data = Transitions(filenames,reduction=reduction)
    MC.set_model(model,data,ncosF,ncosD,ncosDrad,pull)

    # USE INFO from INITFILE
    if initfile is not None:
        import sys
        f = sys.stdout
        MC.use_initfile(initfile)
        MC.print_MC_params(f)
        MC.print_coeffs_laststate(f)

    logger = Logger(MC)

    # MONTE CARLO OPTIMIZATION
    do_mc_cycles(MC,logger)

    # print final results (potential and diffusion coefficient)
    #----------------------------------------------------------
    # choose filename for pickle object
    if outfile is None:
        import sys
        f = sys.stdout
        picfile = "mc.pic"
    else:
        f = file(outfile,"w+")  # print final model to a file
        picfile = outfile+".pic"

    # print to screen
    #MC.print_log_like()
    MC.print_statistics()

    MC.print_laststate(f,final=True)  # print model, coeffs
    if outfile is not None:
        f.close()

    logger.model = MC.model   # this is not a hard copy
    logger.dump(picfile)
    logger.statistics(MC)  #st=1000)

    # write average of pic file
    ave_file = picfile+".ave.dat"
    write_average_from_pic(picfile,ave_file)

    return()

