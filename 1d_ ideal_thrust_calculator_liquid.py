import math 

# Constants
g0 = 9.81  # m/s^2

# Inputs 
Pc = 3e6        # chamber pressure (Pa)
Tc = 3500       # chamber temperature (K)
gamma = 1.22
R = 355         # J/(kg*K) typical combustion gas
At = 0.001      # throat area (m^2)
Pa = 101325     # ambient pressure (Pa)

# Exhaust velocity
ve = math.sqrt(
    (2 * gamma / (gamma - 1)) * R * Tc *
    (1 - (Pa / Pc) ** ((gamma - 1) / gamma))
)

# Mass flow rate
mdot = (
    At * Pc *
    math.sqrt(gamma / (R * Tc)) *
    (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
)

# Thrust
thrust = mdot * ve

# Specific impulse
Isp = ve / g0

# Output
print(f"Exhaust velocity: {ve:.1f} m/s")
print(f"Mass flow rate: {mdot:.2f} kg/s")
print(f"Thrust: {thrust:.1f} N")
print(f"Specific impulse: {Isp:.1f} s")
