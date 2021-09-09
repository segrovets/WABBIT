#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import csv

ddebug_path = './ddir/forces_rk.t'
tsubame_path = '/home/ii/titech/wabbit/falling-sphere-tsub/forces_rk.t'

ddebug = [[],[]]
tsub   = [[],[]]

def _getForces(path, xf):
    with open(path, newline = '\n') as forces:
        reader = csv.reader(forces, delimiter=' ')
        for row in reader:
                xf[1].append(float(row[1]))
                xf[0].append(float(row[-1]))
                

_getForces(ddebug_path, ddebug)
_getForces(tsubame_path, tsub)

plt.plot(tsub[1],tsub[0],ddebug[1],ddebug[0],'*')
plt.legend(['tsubame-oldwabbit','currently-debugging'])
plt.title('debug verification falling sphere')
plt.xlabel('time "forces_rk.t[0]"')
plt.ylabel('forces "force_rk.t[-1]"')
plt.xscale('log')
plt.savefig('/home/ii/opt/WABBIT/comparison_debug.png')
plt.show()



