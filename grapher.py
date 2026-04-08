import numpy as np
import matplotlib.pyplot as plt

# =========================
# 1. FREE FALL
# =========================
try:
    data = np.loadtxt("freefall.csv")
    t = data[:,0]
    z = data[:,1]

    g = 9.81
    z0 = z[0]
    z_analytical = z0 - 0.5 * g * t**2

    plt.figure()
    plt.plot(t, z, label="Numerical")
    plt.plot(t, z_analytical, '--', label="Analytical")
    plt.xlabel("Time")
    plt.ylabel("Height (z)")
    plt.title("Free Fall: Numerical vs Analytical")
    plt.legend()
    plt.grid()
    plt.savefig("freefall_comparison.png")
    plt.show()

except:
    print("freefall.csv not found")

# =========================
# 2. CONSTANT VELOCITY
# =========================
try:
    data = np.loadtxt("constvel.csv")
    t = data[:,0]
    x = data[:,1]
    v = data[:,2]

    plt.figure()
    plt.plot(t, x)
    plt.xlabel("Time")
    plt.ylabel("Position (x)")
    plt.title("Constant Velocity: Position vs Time")
    plt.grid()
    plt.savefig("constvel_position.png")
    plt.show()

    plt.figure()
    plt.plot(t, v)
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.title("Constant Velocity: Velocity vs Time")
    plt.grid()
    plt.savefig("constvel_velocity.png")
    plt.show()

except:
    print("constvel.csv not found")

# =========================
# 3. BOUNCE TEST
# =========================
try:
    data = np.loadtxt("bounce.csv")
    t = data[:,0]
    z = data[:,1]
    vz = data[:,2]

    plt.figure()
    plt.plot(t, z)
    plt.xlabel("Time")
    plt.ylabel("Height (z)")
    plt.title("Bounce: Height vs Time")
    plt.grid()
    plt.savefig("bounce_height.png")
    plt.show()

    plt.figure()
    plt.plot(t, vz)
    plt.xlabel("Time")
    plt.ylabel("Velocity (z)")
    plt.title("Bounce: Velocity vs Time")
    plt.grid()
    plt.savefig("bounce_velocity.png")
    plt.show()

    # Kinetic Energy
    m = 1.0
    KE = 0.5 * m * vz**2

    plt.figure()
    plt.plot(t, KE)
    plt.xlabel("Time")
    plt.ylabel("Kinetic Energy")
    plt.title("Bounce: Kinetic Energy vs Time")
    plt.grid()
    plt.savefig("bounce_energy.png")
    plt.show()

except:
    print("bounce.csv not found")

# =========================
# 4. ERROR vs TIMESTEP (optional)
# =========================
try:
    data = np.loadtxt("error.csv")
    dt = data[:,0]
    err = data[:,1]

    plt.figure()
    plt.loglog(dt, err, 'o-')
    plt.xlabel("Timestep")
    plt.ylabel("Error")
    plt.title("Error vs Timestep")
    plt.grid()
    plt.savefig("error_vs_dt.png")
    plt.show()

except:
    print("error.csv not found")