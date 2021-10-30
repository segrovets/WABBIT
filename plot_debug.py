#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import csv

ddebug_path = './ddir/forces_rk.t'
insect_state_path = './ddir/insect_state_vector.t'
tsubame_path = '/home/ii/titech/wabbit/falling-sphere-tsub/forces_rk.t'
python_path = '/home/ii/titech/pendulum-tree/python_pendulum.t'

ddebug = []
tsub   = []
insect_state = []
python_pend = []

def _getForces(path, xf):
    with open(path, newline = '\n') as f:
        for row in f:
            columns = row.split(' ')
            columns = [float(dat.strip()) for dat in columns if (dat != '')]
            xf.append(columns)

                

_getForces(ddebug_path, ddebug)
_getForces(insect_state_path, insect_state)
_getForces(python_path, python_pend)
#_getForces(tsubame_path, tsub)

fig = plt.figure(figsize=(10,20))
ax = fig.add_subplot(311)
ax.grid()
#plt.plot(tsub[1],tsub[0],ddebug[1],ddebug[0],'*')
_time =[d[0] for d in ddebug]
_timep = [p[0] for p in python_pend]
rho = 1
mu = 0.01
scaling_factor = rho/mu**2
rhop = 1.2
mup = (1.458e-6)*297**(3./2)/(297+110.4)
scaling_factor_python = rhop/mup**2
k=0

for i in range(1,4):
    ax.plot(_time,[d[i]*scaling_factor for d in ddebug],'*')
ax.plot(_time,[np.sqrt(d[1]**2 + d[3]**2)*scaling_factor for d in ddebug],".")
ax.plot([_timep[i]-_timep[5000] for i in range(5000,len(_timep))], [python_pend[p][-1]*scaling_factor_python for p in range(5000,len(python_pend))], "+")
# t
#plt.legend(['tsubame-oldwabbit','currently-debugging'])
ax.set(ylim=[-1000,3000])
ax.legend(['fx','fy','fz','ft',"f-python"])

ar = [max(d[1:4])*scaling_factor for d in ddebug ]
arp = [p[-1]*scaling_factor_python for p in python_pend[0:len(_time)*10] ]
ar2 =[min(d[1:3])*scaling_factor for d in ddebug ]
arp2 = [p[-1]*scaling_factor_python for p in python_pend[0:len(_time)*10] ]
ax.set(title='debug verification inverse pendulum tree',ylabel='forces * rho / mu ^ 2',xlim=[0,_time[-1]],ylim=[-10000,15000])

# This is in INSECT%STATE:
# STATE(1) : x-position of body
# STATE(2) : y-position of body
# STATE(3) : z-position of body
# STATE(4) : x-velocity of body
# STATE(5) : y-velocity of body
# STATE(6) : z-velocity of body
# STATE(7) : 1st component of body quaternion
# STATE(8) : 2nd component of body quaternion
# STATE(9) : 3rd component of body quaternion
# STATE(10) : 4th component of body quaternion
# STATE(11) : x-angular velocity of body *not sure if around the body mass center?
# STATE(12) : y-angular velocity of body
# STATE(13) : z-angular velocity of body
# STATE(14) : angle x -- around origin ""my additions
# STATE(15) : angle y -- around origin
# STATE(16) : angle z -- around origin
# STATE(17) : angular velocity -- around origin
# STATE(18) : angular velocity -- around origin
# STATE(19) : angular velocity -- around origin
# STATE(20-25) : 0
x = [i[0] for i in insect_state]
y = [i[1] for i in insect_state]

_time = np.linspace(_time[0],_time[-1],len(x)) 

markers = ['*','.','-']

ax2 = fig.add_subplot(312)
#for x in insect_state:
#    print(np.arcsin(x[1]/8.), np.arcsin(x[1]/8-np.pi/2)*180/np.pi)
    
ang_pos = [np.arccos(x[3]/8)*180/np.pi for x in insect_state] 
#ang_vel = [i[j+3]*180.0/np.pi for i in insect_state]   

pend_time =[_timep[i]-_timep[5000] for i in range(5000,len(_timep))]
pend_y = [python_pend[p][1] for p in range(5000,len(python_pend))]
ax2.plot(_time, ang_pos, '*')#, _time, ang_vel)
ax2.plot(pend_time, pend_y, '+')
ax2.legend(['angular position x', 'python angular position'])
ax2.set(xlim=[0,_time[-1]])
ax2.set(ylabel='deg')

ax3 = fig.add_subplot(313)
ang_vel = [np.sqrt(i[4]**2. + i[6]**2)/8. for i in insect_state] # sqrt( Vx^2 + Vz^2)/ r   
ax3.plot(_time, ang_pos, markers[k])#, _time, ang_vel)
pend_y = [python_pend[p][2] for p in range(5000,len(python_pend))]
ax3.plot(pend_time, pend_y, '+')
ax3.legend(['angular velocity x', "python ang vel"])
#ax2.title('ang pos of inverse pendulum tree')
ax3.set(xlabel='time',ylabel='deg',xlim=[0,_time[-1]])

plt.savefig('/home/ii/opt/WABBIT/comparison_debug2.png')
plt.show()
