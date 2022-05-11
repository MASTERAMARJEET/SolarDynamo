import math
import numpy as np
import matplotlib.pyplot as plt

# Data obtained from http://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/562/L10#/browse

input_file_name = "GNhs_y.txt"

f = open(input_file_name, "r")
f1 = f.readlines()

grp = []
yrs = []
for x in f1:
    val = float(x[7:13].strip())
    yrs.append(float(x[:6]))
    grp.append(val)
    # print(val)

grp = np.array(grp)
f.close()

decadal_yrs = []
decadal_grp = []
count = 0
for ii in range(len(grp)):
    count += 1
    if count == 10:
        decadal_yrs.append(np.mean(yrs[ii-9:ii+1]))
        decadal_grp.append(np.mean(grp[ii-9:ii+1]))
        count = 0

decadal_yrs = np.array(decadal_yrs)
decadal_grp = np.array(decadal_grp)

aa, bb, cc, dd = np.loadtxt("sn.dat", unpack = True)
xx = np.linspace(-1150, 1950, 310, endpoint = True)

sigmas = (dd - cc)/(2*1.96)
n_draws = int(1e3)
new_data = np.ndarray.flatten(np.outer(sigmas,np.random.randn(n_draws)).T + bb)

# plt.figure(figsize = (8, 3.6))
# plt.plot(xx, bb, 'r', linewidth = 1.1)
# plt.fill_between(xx, cc, dd, color='black', alpha=0.3)
# plt.plot(decadal_yrs, 12*decadal_grp, 'b.-')
# # plt.errorbar(xx, bb, yerr = cc)
# # plt.errorbar(xx, bb, yerr = dd)
# # , 'k', alpha = 0.5, linewidth = 1.2
# plt.xlim(-1150, 2000)
# plt.ylim(0, 90)
# plt.axhline(y = 21, color = 'k', linestyle = '--')
# # plt.axhline(y = 67)
# plt.ylabel('Sunspot Number', fontsize = 14)
# plt.xlabel('Years (-BC/AD)', fontsize = 14)
# plt.savefig('RconsSSN.png', dpi = 300)
# plt.show()

# max_range = math.ceil((np.max(bb)-np.min(bb)))
bin_size = int(1e2)
plt.hist(new_data, bins = bin_size, density = True, color = 'r', alpha = 0.5, label='Reconstructed SSN')
plt.hist(12*decadal_grp, bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90], histtype = 'step', density = True, color = 'blue', label='Observed SSN')
plt.ylabel('PDF', fontsize = 14)
plt.xlabel('Sunspot number', fontsize = 14)
plt.axvline(x = 21, color = 'k', linestyle = '--')
plt.legend()
plt.show()
