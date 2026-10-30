import numpy as np
import matplotlib.pyplot as plt

E  = 70e9      # Pa
A  = 12.5e-4   # m^2 (this is A0)
L  = 0.5
P  = -5e3

meshes = [1, 2, 4, 8, 16, 32, 64]
results = {}                       # results[n] = (x, u, stress)
freeEndError = []

# exact solution
uExact = lambda x: -(abs(P)*L/(E*A))*np.log(2/(1 + x/L))
sExact = lambda x: abs(P)/(A*(1 + x/L))

for numberElements in meshes:
    nodeCoordinates = np.linspace(0, L, numberElements + 1)
    numberNodes = numberElements + 1
    elementNodes = [[i, i + 1] for i in range(numberElements)]

    GDOF = numberNodes
    displacements = np.zeros(GDOF)

    F = np.zeros(GDOF)
    F[0] = P

    stiffness = np.zeros((GDOF, GDOF))
    Lelem = np.zeros(numberElements)
    for element, (n1, n2) in enumerate(elementNodes):
        x1, x2 = nodeCoordinates[n1], nodeCoordinates[n2]
        Lelem[element] = x2 - x1
        Ae = A*(1 + ((x1 + x2)/2)/L)
        stiffness[np.ix_([n1, n2], [n1, n2])] += E*Ae/Lelem[element]*np.array([[1, -1],
                                                                             [-1, 1]])

    bc = [numberNodes - 1]
    free = [i for i in range(GDOF) if i not in bc]
    displacements[free] = np.linalg.solve(stiffness[np.ix_(free, free)], F[free])

    # element stress (constant per element)
    stress = np.array([E*(displacements[n2] - displacements[n1])/Lelem[e]
                       for e, (n1, n2) in enumerate(elementNodes)])

    results[numberElements] = (nodeCoordinates, displacements, stress)
    freeEndError.append(abs(displacements[0] - uExact(0)))

# ---------- plots ----------
xf = np.linspace(0, L, 500)
plotMeshes = [1, 2, 4, 8, 16]

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
for n in plotMeshes:
    x, u, s = results[n]

    # displacement error (FE is linear inside each element)
    ax[0].plot(xf, (np.interp(xf, x, u) - uExact(xf))*1e6, label=f'{n} el')

    # stress, plotted as a staircase
    ax[1].plot(np.repeat(x, 2)[1:-1], np.repeat(s, 2)/1e6, label=f'{n} el')

    # stress error on the fine grid
    idx = np.clip(np.searchsorted(x, xf, side='right') - 1, 0, n - 1)
    ax[2].plot(xf, (s[idx] - sExact(xf))/1e6, label=f'{n} el')

ax[1].plot(xf, sExact(xf)/1e6, 'k--', label='exact')

ax[0].set(title='Displacement error (FE - exact)', xlabel='x (m)', ylabel='error (µm)')
ax[1].set(title='Stress', xlabel='x (m)', ylabel='σ (MPa)')
ax[2].set(title='Stress error (FE - exact)', xlabel='x (m)', ylabel='error (MPa)')
for a in ax:
    a.grid(True); a.legend()
plt.tight_layout()

# convergence of the free-end displacement
plt.figure()
plt.loglog(meshes, freeEndError, 'o-')
plt.xlabel('number of elements')
plt.ylabel('|u_FE(0) - u_exact(0)| (m)')
plt.title('Convergence at the free end')
plt.grid(True, which='both')
plt.show()