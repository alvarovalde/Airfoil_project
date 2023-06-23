import numpy as np
import matplotlib.pyplot as plt
import scipy


def chev_grid_creator(n_points):
    xi = np.cos(np.pi / (2 * (n_points + 1)) * (2 * np.linspace(0, n_points + 1, n_points + 1) - 1))
    return ((-xi + 1) / 2) * (1 - 0) + 0



xf = np.array([0.05, 0.1, 0.15, 0.2, 0.25])
m = np.array([0.0580, 0.1260, 0.2025, 0.2900, 0.3910])

z = np.polyfit(xf, m, 3)
f = np.poly1d(z)

x_new = chev_grid_creator(10)
y_new = f(x_new)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(xf,m,color='red')
ax.scatter(x_new,y_new,color='blue')
plt.show()



Cl = 0.3
xf = np.array([0.05,0.1,0.15,0.2,0.25])
K1 = np.array([361.4,51.64,15.957,6.643,3.230])
K1 = K1*Cl/0.3

f = scipy.interpolate.interp1d(xf,K1,bounds_error=False,fill_value='extrapolate',kind='quadratic',axis=0)


Xf_new = np.linspace(0,1,50)
new_K1 = f(Xf_new)

plt.scatter(xf,K1,color='blue')
plt.scatter(Xf_new,new_K1,color='red')
plt.show()

#===========================================================================================

Cl = 0.3

xf = np.array([0.05,0.1,0.15,0.2,0.25])
K1 = np.array([361.4,51.64,15.957,6.643,3.230])
K1 = K1*Cl/0.3

# Xf_new = chev_grid_creator(100)
# f = scipy.interpolate.CubicSpline(xf,K1)
# y_new = f(Xf_new)
# plt.scatter(xf,K1,color='blue')
# plt.scatter(Xf_new,y_new,color='red')
# plt.show()

#-------------------------------------------------------------------------------------------------------------------



# Define the function that approaches zero as x tends to infinity



plt.scatter(xf, K1, label='Original Data')
plt.scatter(x_interpolation, y_interpolation, label='Interpolation',color='red')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.show()
