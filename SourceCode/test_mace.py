from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary, ZeroRotation
import warnings
from potential import MACEPotential
from ase.md.verlet import VelocityVerlet
from ase.io.trajectory import Trajectory



#pip install cuequivariance cuequivariance-torch to remove one of the warnings
warnings.filterwarnings("ignore", module="mace")
warnings.filterwarnings("ignore", module="torch")
warnings.filterwarnings("ignore", module="e3nn")

a=FaceCenteredCubic(size = (5,5,5), symbol= "Cu",pbc = True)
pot = MACEPotential("MACEModels/mace-mpa-0-medium.model")
a.calc = pot.getPotentialCalculator()

E = a.get_potential_energy()
F = a.get_forces()
print(f"Static check: E={E} eV, max|F|={abs(F).max():} eV/Å")

# 2) initiera hastigheter och kör 3 MD-steg så KE ≠ 0
MaxwellBoltzmannDistribution(a, temperature_K=300, force_temp=True)
Stationary(a); ZeroRotation(a)

dyn = VelocityVerlet(a, timestep=1.0)  # 1 fs
traj = Trajectory("mace_test.traj", "w", atoms=a)

def sample():
    Ep = a.get_potential_energy()
    Ek = a.get_kinetic_energy()
    Et = Ep + Ek
    a.info["Ep"] = float(Ep); a.info["Ek"] = float(Ek); a.info["Et"] = float(Et)
    traj.write()
    print(f"MD: Ep={Ep:} Ek={Ek:} Et={Et:}")
    try:
        print(a.get_stress())
    except Exception as e:
        print("Stress unavailable:", e)

dyn.attach(sample, interval=1)
dyn.run(3)
traj.close()



print(a.get_kinetic_energy())