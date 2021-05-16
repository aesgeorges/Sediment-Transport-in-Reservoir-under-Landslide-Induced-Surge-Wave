# Time advancement functions
# Each function performs a single time-step for a corresponding component
# adv function calls all function to perform a time-step for all components
# Components are:
# particles, u, v, temp, rho, q2, q21, l, kz, nu_t, kq

# Velocity (U,V)
def velocity(c): 
    for i in range(c.N-2):
        return
