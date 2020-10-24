from tkinter import messagebox
from function import *
from VALID_INPUT1 import VALID_INPUT
from math import pi,sin,cos,exp


def err_exit(FUN=0,ERR=0):
    if  not(isinstance(FUN,str)) or not(isinstance(ERR,str)):
        return
    else:
        print('#############################')
        print('ERROR IN THE FUNCTION/PROCEDURE : '+FUN)
        print('ERROR MESSAGE : '+ERR)
        print('#############################')



def SEGPOS(FIRST_INPUT=None,DINT=None):

    # n_params = n_param([FIRST_INPUT, DINT])  当输入两个参数时，为整体镜的标准输入，在该程序中只考虑这种情况
    DEXT = FIRST_INPUT
    NSD = 1
    GAP = 0
    TYPE = 'DISC'
    # VALID_INPUT('SEGPOS.PRO', 'DEXT', DEXT, 'real', 0, 'no', '++', 'free')
    # VALID_INPUT('SEGPOS.PRO', 'DINT', DINT, 'real', 0, 'no', '0+', 'free')
    # VALID_INPUT('SEGPOS.PRO', 'NSD', NSD, 'integer', 0, 'no', '++', 'free')
    # VALID_INPUT('SEGPOS.PRO', 'GAP', GAP, 'real', 0, 'no', '0+', 'free')
    # VALID_INPUT('SEGPOS.PRO', 'TYPE', TYPE, 'string', 0, 'no', ['hexa', 'disc', 'squa', 'HEXA', 'DISC', 'SQUA'], 'free')
    if DINT >= DEXT :
        messagebox.showinfo('','WRONG OR INCOMPATIBLE DEXT & DINT')

    # 在此模型下 以下两句无用
    if NSD > 1 :
        if GAP >= float(DEXT)/(NSD-1) :
            messagebox.showinfo('','GAP TOO LARGE, SEGMENT SIZE BECOMES NULL !')
    if TYPE != 'SQUA' and NSD % 2 == 0 :
        messagebox.showinfo('','NSD MUST BE ODD WITH THIS SEGMENT SHAPE')

    pos=np.zeros(shape=(1,2))

    tns=1

    if TYPE.upper()  == 'DISC':
        swd = float(DEXT)
        ssz = 0.5*float(DEXT)
        dbs = 0.0

    if TYPE.upper()  == 'DISC' :
        surf=tns*pi*ssz**2
    if tns == 1 :
        surf=surf-pi*(0.5*DINT)**2


    if TYPE.upper()  == 'DISC' :
        dextmax=2*ssz

    dxy = np.zeros(shape=(2, int(tns)))

    result={'TNS':tns,'SWD':swd,'SSZ':ssz*2,'DBS':dbs,'POS':pos,'DEXTMAX':dextmax,'DEXT':DEXT,'DINT':DINT,'TYPE':TYPE,'GAP':GAP,'NSD':NSD,'SURF':surf,'REMSEG':'n'}

    return(result)

