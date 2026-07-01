import math

# Constants
g0 = 9.81  # m/s^2

# Inputs
Tc = 3200       # chamber temperature (K)
gamma = 1.16
R = 291         # J/(kg*K) typical combustion gas
At = 1.256e-5      # throat area (m^2)
Pa = 101325     # ambient pressure (Pa)
Ae = 7.23e-5   # exit area (m^2)
exit_deg = 15   # exit cone angle in degrees 

# Propellant properties
a = .000155       # burn rate coefficient (kg/m^2/s) (rough estimate based on typical hybrid propellant)
n = .5       # burn rate exponent (rough estimate based on typical hybrid propellant)
Pp = 930      # propellant density (kg/m^3)
# Grain geometry
L = .12 # grain length in meters
r_core = .0064 # radius of hollow core in meters 
Ro = .0254 # radius of grain edge touching the case in meters

#define time steps 
t = 0
dt = .01

#store data for graph
time_data = []
thrust_data = []
pressure_data = []
of_ratio_data = []

total_impulse = 0

#define mass flow rate of oxidizer (will be made into a function later)
mdot_ox = .035 # mass flow rate of oxidizer (kg/s) - this is a guess, will need to be refined based on desired O/F ratio

#start sim loop
while r_core < Ro:
    Ab = (2* math.pi * r_core *L)     # burning surface area (m^2)
#check for ratio issues
    if Ae<At:
        print("Error: Exit area must be greater than throat area for supersonic flow.")
        break # stop the simulation if the area ratio is invalid
# refine Me
    else:# solve for mach exit number (roughly)
# rough starting point
        Me = 3.0
        area_ratio = Ae/At

        for _ in range(10):
            f = (1/Me) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * Me**2))**((gamma+1)/(2*(gamma-1))) - area_ratio
            df = (f + area_ratio) *(Me**2 - 1) / (Me * (1 + (gamma -1) / 2 * Me**2))
            Me = abs(Me - f/df)
#check flow type
            if Me >1.0:
                flow_type = "supersonic"
            elif Me ==  1.0:
                flow_type = "sonic (choked)"
            else:
                flow_type = "subsonic"

# solve for C
        C = math.sqrt(gamma/(R*Tc))*(2/(gamma + 1))**((gamma + 1)/(2*(gamma - 1)))

#burn rate
        Aport = math.pi * (r_core**2) # cross sectional area of the core (m^2)
        
        Gox = mdot_ox/Aport # mass flux of oxidizer (kg/m^2/s)

        Rb = a * Gox**n #burn rate (m/s)

# mass flow rate of fuel 
        mdot_fuel = Pp * Ab * Rb  # mass flow rate of fuel (kg/s)

        mdot = mdot_ox + mdot_fuel # total mass flow rate (kg/s)

        O_F_ratio = mdot_ox/mdot_fuel # oxidizer to fuel ratio

# solve for flow rate coefficient

        C = math.sqrt(gamma/(R*Tc))*(2/(gamma + 1))**((gamma + 1)/(2*(gamma - 1))) 

# solve for chamber pressure
        Pc = (mdot / (At * C))

#solve for exit pressure
        Pe = Pc*((1 + ((gamma - 1)/2)*(Me**2))**(-gamma/(gamma - 1)))

# Exhaust velocity
        ve = math.sqrt(
            (2 * gamma / (gamma - 1)) * R * Tc *
            (1 - (Pe / Pc) ** ((gamma - 1) / gamma))
        )

# Mass flow rate (throat)/(choked flow)
        Mtdot = ((At * Pc) * math.sqrt(gamma / (R * Tc))*(2 / (gamma + 1))** ((gamma +1)/(2*(gamma -1))) )

# Thrust
        alpha_rad = math.radians(exit_deg)
        lambda_coeff = (1 + math.cos(alpha_rad)) / 2
        thrust = (lambda_coeff * mdot * ve) + (Pe - Pa)*Ae

#thrust in pounds of force
        thrust_lbf = thrust/4.44822
# Specific impulse
        Isp = thrust / (mdot*g0)

#store data for graph
        time_data.append(t)
        thrust_data.append(thrust)
        pressure_data.append(Pc/1e6) #keep in MPa for graph 
        total_impulse += (thrust*dt)
        of_ratio_data.append(O_F_ratio)

#advance time and grain geometry
        r_core = r_core + (Rb*dt)
        t = t + dt

# Output
print(f"Final Exhaust velocity: {ve:.1f} m/s")
print(f"Final Thrust: {thrust:.1f} N")
print(f"Final Chamber Pressure: {Pc/1e6:.2f} MPa")
print(f"Total Impulse: {total_impulse:.2f} N-s")
print(f"Burn complete. Total burn time: {t:.2f} seconds")

#plot thrust curve and of ratio curve
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(time_data, thrust_data, label="Thrust (N)", color="red", linewidth=2)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Thrust (N)", color="red")
ax1.tick_params(axis='y', labelcolor="red")
ax1.set_ylim(bottom=0)

# Create a second y-axis to show the O/F Ratio changing!
ax2 = ax1.twinx()
ax2.plot(time_data, of_ratio_data, label="O/F Ratio", color="blue", linestyle="--")
ax2.set_ylabel("O/F Ratio", color="blue")
ax2.tick_params(axis='y', labelcolor="blue")

plt.title("Hybrid Motor Thrust Curve & O/F Shift", fontweight="bold")
plt.grid(True, linestyle="--", alpha=0.7)
fig.tight_layout()
plt.show()