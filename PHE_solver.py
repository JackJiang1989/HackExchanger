from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp


def cal_N(flowrate,Cp,K, passes, A_oneplate):
    return passes*K*A_oneplate/(flowrate*Cp)

def passflag(number):
    if number%2==0:
        return 1
    if number%2==1:
        return -1

def hcflag(number):
    if number==0:
        return 1
    if number==1:
        return -1

def cal_Nmatrix(passchannlist, N, posflag=[1,1]):
#starter=1 means flow from left to right
#passchannlist=[[0,2],[4,6],[8,10],[12,14]]
#hcflag,hot=1, cold=-1
#below + or - is based on pass number is odd and hot
#inletflag, lefttoright=-1, righttoleft=1
    total_chan=len(passchannlist[0])*len(passchannlist[0][0])*2
    return_list=np.zeros((total_chan, total_chan))
    #if passchannlist[0][0] == 0: #if there is 0, means the first channel is there, so we set it
    #    return_list[0][0] = 1 * hcflag * posflag * N
    #    return_list[0][1] = -1 * hcflag * posflag * N
    #    del passchannlist[0][0] #delete otherwise it will appear in the loop below
    #elif passchannlist[0][0] == 1: #if no 0, means there is the last channel, then we just set it
    #    return_list[-1][-1] = -1 * hcflag * passflag(len(passchannlist)) * inletflag * N #depending on the total passes number, the flow direction changed
    #    return_list[-1][-2] = 1 * hcflag * passflag(len(passchannlist)) * inletflag * N
    #    del passchannlist[-1][-1] #delete otherwise it will appear in the loop below
    #elif passchannlist[-1][-1] == 1:
    #    return_list[-1][-1] = -1 * hcflag * inletflag * N
    #    return_list[-1][-2] = 1 * hcflag * inletflag * N
    #    del passchannlist[0][0]
    #elif passchannlist[-1][-1] == 0:
    #    return_list[0][0] = 1 * hcflag * posflag * passflag(len(passchannlist)) * N
    #    return_list[0][1] = -1 * hcflag * posflag * passflag(len(passchannlist)) * N
    #    del passchannlist[-1][-1]

#posflag==1入口位置从左到右，入口位置posflag==-1从右到左
    for sn,side in enumerate(passchannlist):
        for pn, pl in enumerate(side):
            for c in pl:
                coff_flag=hcflag(sn)*posflag[sn]*passflag(pn)*N[sn]

                if c!=0:return_list[c][c-1]=-1*coff_flag

                if c!=0 and c!=(total_chan-1):return_list[c][c]=2*coff_flag
                elif c==0:return_list[c][c]=1*coff_flag
                elif c==(total_chan-1):return_list[c][c]=1*coff_flag

                if c!=(total_chan-1):return_list[c][c+1]=-1*coff_flag
    return return_list

def constr_fun(total_channels, N):
    def fun(x, y):
        dydx=np.dot(N,y)
        t = tuple(dydx[i] for i in range(total_channels))
        return np.vstack(t)
    return fun

def ya_or_yb(p,flag, y):
#0,ya,1=>ya
#0,ya,0=>yb
#1,ya,1=>yb
#1,ya,0=>ya

    if (p%2)^flag:
        return y
    else:
        if y == 'ya':
            return 'yb'
        if y== 'yb':
            return 'ya'

# print ya_or_yb(0,1,'yb')
# print ya_or_yb(1,1,'yb')
# print ya_or_yb(0,0,'yb')
# print ya_or_yb(1,0,'yb')

def constr_bc(passchannlist, st_y,inlet_t, outlet_t):
# passchannlist=[ [[0,2],[4,6]],[[1,3],[5,7]] ]
# st_y=['ya','ya'] st_y[0]=>first side, st_y[1]=>second side
# inlet_t=[60,None]
# outlet_t=[None,30]
    def bc_test(ya,yb):
        bc_list=[]
        ydict={'ya':ya, 'yb':yb}
        channel_nr=[len(passchannlist[0][0]),len(passchannlist[1][0])]
        for sn,side in enumerate(passchannlist):   #[[ [0,2],[4,6] ],[ [1,3],[5,7] ]]
            for pn, pl in enumerate(side):
                temp4sum_mid=[]
                for c in pl:
                    temp4sum_mid.append(ydict[ya_or_yb(pn,0,st_y[sn])][c])
                    if pn != 0:
                        bc_list.append(current_sum/channel_nr[sn] - ydict[ya_or_yb(pn,1,st_y[sn])][c])
#                        print 1
                    elif inlet_t[sn] and pn==0:
                        bc_list.append(ydict[st_y[sn]][c]-inlet_t[sn])
#                        print 2
                    elif outlet_t[sn] and pn == (len(passchannlist)-1):
                        bc_list.append(ydict[ya_or_yb(pn,0,st_y[sn])][c]-outlet_t[sn])
#                        print 3
#                print bc_list
                current_sum=sum(temp4sum_mid)
        return np.array(bc_list)
    return bc_test


def contr_passchannlist(hp, hc,cp,cc, flowdir):
    if flowdir=='cocurrent':
        h_lst = [range(0,2*h_passes*h_channels,2)[i:i+h_channels] for i in range (0,h_passes*h_channels,h_channels)]
        c_lst = [range(1,2*c_passes*c_channels,2)[i:i+c_channels] for i in range(0,c_passes*c_channels,c_channels)]
    if flowdir=='countercurrent':
        h_lst = [range(0,2*h_passes*h_channels,2)[i:i+h_channels] for i in range (0,h_passes*h_channels,h_channels)]
        c_lst = [range(1,2*c_passes*c_channels,2)[::-1][i:i+c_channels] for i in range(0,c_passes*c_channels,c_channels)]
    return [h_lst, c_lst]
     

h_flow=1 #m/s
c_flow=1 #m/s
K=7.675 #kw/m2k
h_Cp=4.175 #kj/kgC
c_Cp=4.175 #kj/kgC
#flowdir="countercurrent"
h_passes=2
h_channels=2
c_passes=2
c_channels=2
A_oneplate=0.18

position_pc=[1,1]
postion_bc=['yb','yb']


h_N=cal_N(h_flow,h_Cp,K, h_passes, A_oneplate)
c_N=cal_N(c_flow,c_Cp,K, c_passes, A_oneplate)

pclist = contr_passchannlist(h_passes,h_channels,c_passes,c_channels,'countercurrent',position_pc)
N1 = cal_Nmatrix
tc1=h_passes*h_channels+c_passes*c_channels
fun1 = constr_fun(tc1,N1)
bc1 = constr_bc(pclist, postion_bc,inlet_t, outlet_t)

x1=np.linspace(0,1,20)
y_starter = np.zeros((tc1, x1.size))

r1 = solve_bvp(fun1, bc1, x1, y_starter)