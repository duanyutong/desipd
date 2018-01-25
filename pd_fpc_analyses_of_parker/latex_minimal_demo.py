# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 16:12:53 2016

@author: Duan Yutong (dyt@physics.bu.edu, dyt@lbl.gov)

"""

#latex minimal example
import numpy as np
import matplotlib
import os
#matplotlib.use('PS')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.unicode']=True
import matplotlib.pyplot as plt
# plt.switch_backend('PS')

string1 = r'z=$\mathrm{value}^{upper}_{lower}$'.format(
                value='{adu}',
                upper='{+' + str(0.01) + '}',
                lower='{-' + str(0.01) + '}')

dir_save = r'C:\Users\givoltage\Downloads\linear'

fig = plt.figure(figsize=(3,1))
fig.text(0.1,0.5, string1, size=24,va='center')
fig.savefig(os.path.join(dir_save, 'source-exptime.png'))

t = np.arange(0.0, 1.0 + 0.01, 0.01)
s = np.cos(4 * np.pi * t) + 2

# === another demo ===

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(t, s)

plt.xlabel(r'\textbf{time} /s \cdot)')
# plt.ylabel(r'$\text{} \cdot $'.format('{ADU}'),fontsize=16)
# plt.ylabel(r'\textit{voltage} (mV)',fontsize=16)
plt.title(r'\TeX\ is Number ' +
          r'$\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}!\cdot$',
          fontsize=16, color='gray')
# Make room for the ridiculously large title.
plt.subplots_adjust(top=0.8)

os.chdir(r'C:\Users\givoltage\Downloads\linear')

plt.savefig('tex_demo')

# from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('issue5076.pdf')
#pp.savefig(fig)
#pp.close()