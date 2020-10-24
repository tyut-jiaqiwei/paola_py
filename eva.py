
from numpy import pi,sin,cos,exp
import numpy as np
from function import *
from PIXMATSIZE2 import PIXMATSIZE
from PSFOTFTSC2 import PSFOTFTSC
from segpos2 import SEGPOS
from PSD_ANISO_SERVO_GLAO_STAR2  import PSD_ANISO_SERVO_GLAO_STAR
import random
import readallnew
import time
from paola import paola


def eva(lgslist,fov, atmlist, aperture):




    #height = profdata[indx].reshape((100,2))[:,0]
    height = atmlist[:, 0]
    # dcn2 = profdata[indx].reshape((100,2))[:,1]
    dcn2 = atmlist[:, 1]
    windir = atmlist[:, 3]* pi / 180
    #print(atmlist[:, 3])
    wind = np.zeros([2,100])
    v = atmlist[:, 2]
    ang06 = lgslist
    #print(atmlist[:, 2])

    paola_start = time.time()

    # height = np.array([48.0, 162, 324, 649, 1299, 2598, 5196, 10392, 20785]) # 0 = ground level
    # dcn2 = np.array([53.28, 1.45, 3.5, 9.57, 10.83, 4.37, 6.58, 3.71, 6.71]) # Cn2 distribution
    dcn2 = dcn2 / sum(dcn2)
    # v = np.array([15.0, 13, 13, 9, 9, 15, 25, 40, 21])                      # wind velocity profile
    # windir = np.array([0, 10, -20, 30, -60, -70, -80, 90, 0]) * pi / 180
    vx = v * np.cos(windir)
    vy = v * np.sin(windir)
    # wind = np.zeros([2,9])
    wind[0:] = vx
    wind[1:] = vy
    w0 = 0.7
    L0 = 27.0
    ZA = 0

    '''GLAO MODES'''
    '''we will use here a simple telescope architecture - just a monolithic mirror'''
    mir = SEGPOS(aperture,1.8) #SEGPOS函数
    dxf = -1
    n_psf = -1 # default value = such that FoV = 8 times seeing limited PSF FWHM
    lam =1.25
    wfspitch = -1 # means that WFS pitch = r0 @ lambda
    dmh = 0 # conjugation height of the DM, here pupil level
    dm_params={'dmtf':-1,'actpitch':-1,'dm_height':dmh}                                                     ####字典的键要不要用字符，IDL中不用，py中用
    wfs_params={'wfs_pitch':wfspitch}
    ang = 20 #science object off-axis angle [asec]
    ori = 60 #science object position angle in deg/x-axis
    wfs_int=10
    lag=5 #loop time lag (WFS reading + DM commands calculation) in msec
    '''Here we have 6 stars on a circle of radius 60 arcsec'''
    #FoVrad=420 # [asec]
    FoVrad = fov
    #print('ang06='+str(2*pi/6*pi*np.array([0,1,2,3,4,5])))
    glao_wfs={'type':'star','ang':ang06}
    info=1

    psg2=PIXMATSIZE(mir,dxf,n_psf,lam,w0,L0,ZA,height,dcn2,wind,wfspitch,dmh,ang,ori,wfs_int,lag,glao_wfs,info)
    #now we need the telescope OTF
    tsc=PSFOTFTSC(mir,psg2)
    # GLAO modeling, giving a NEA for the WFS noise error
    gs_weight=-1 # all NGS of the constellation are given the same weight
    wfs_nea=np.zeros(6) +0.02 # WFS Noise Equivalent Angle / NGS, in asec
    # print("wfs_nea="+str(wfs_nea))

    glao_star1=paola('glao',psg2,tsc,w0,L0,ZA,height,dcn2,wind,dm_params,wfs_params,ang,ori,wfs_int,lag,'open',1,glao_wfs,gs_weight,wfs_nea,\
                  INFO=1,OTF=1,PSF=1,SF=1,ONLY_PSD=1,LOGCODE='star1')

    paola_end = time.time()
    # print('评价函数（PSD）程序总用时' + str(paola_end - paola_start))

    return glao_star1
