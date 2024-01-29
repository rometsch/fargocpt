from dataclasses import dataclass
import numpy as np
from types import SimpleNamespace
import argparse


@dataclass
class FargoCPTField:
    outputdir: str
    snapshotid: str
    name: str

    def __post_init__(self):
        self.filename = f"{self.outputdir}/snapshots/{self.snapshotid}/{self.name}.dat"
        self.grid = get_fargo_grid(self.outputdir)
        if self.name == "vrad":
            self.R, self.Phi = np.meshgrid(self.grid.ri, self.grid.phic, indexing="ij")
        elif self.name == "vtheta":
            self.R, self.Phi = np.meshgrid(self.grid.rc, self.grid.phii[:-1], indexing="ij")
        else:
            self.R, self.Phi = np.meshgrid(self.grid.rc, self.grid.phic, indexing="ij")
        self._load()

    def _load(self):
        self._data = np.fromfile(self.filename, dtype=np.float64)
        self._data = self._data.reshape(self.R.shape[0], self.R.shape[1])

    def save(self, altid=None):
        if altid is not None:
            filename = f"{self.outputdir}/snapshots/{altid}/{self.name}.dat"
        self._data.tofile(self.filename)

    @property
    def array(self):
        return self._data

    @array.setter
    def array(self, data):
        if not self._data.shape == data.shape:
            raise ValueError("Shape of data does not match shape of field")
        self._data = data
        

def get_fargo_grid(outputdir):

    Nrad, Naz = np.genfromtxt(f"{outputdir}/dimensions.dat", usecols=(4,5), dtype=int, unpack=True)

    ri = np.genfromtxt(f"{outputdir}/used_rad.dat")
    phii = np.linspace(0, 2*np.pi, Naz+1)
    Ri, Phii = np.meshgrid(ri, phii, indexing="ij")
    Xi = Ri*np.cos(Phii)
    Yi = Ri*np.sin(Phii)

    rc = 2/3*(ri[1:]**2/(ri[1:]+ri[:-1]) + ri[:-1]) # approx center in polar coords
    phic = 0.5*(phii[1:]+phii[:-1])
    Rc, Phic = np.meshgrid(rc, phic, indexing="ij")
    Xc = Rc*np.cos(Phic)
    Yc = Rc*np.sin(Phic)

    dphi = phii[1] - phii[0]
    dr = ri[1:] - ri[:-1]
    A = 0.5*(Ri[1:,1:]**2 - Ri[:-1,1:]**2)*dphi

    return SimpleNamespace(
        Nrad=Nrad, Naz=Naz,
        ri=ri, phii=phii, Ri=Ri, Phii=Phii, Xi=Xi, Yi=Yi,
        rc=rc, phic=phic, Rc=Rc, Phic=Phic, Xc=Xc, Yc=Yc, 
        dphi=dphi, dr=dr, A=A
    )
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outputdir")
    opts = parser.parse_args()
    
    field = FargoCPTField(opts.outputdir, 0, "Sigma")
    R = field.R
    print(R)
    Phi = field.Phi
    print(Phi)
    R0 = 4
    phi1 = np.pi
    phi2 = np.pi/2
    Deltar = 1
    Deltaphi = 0.3
    
    gauss1 = np.exp(-0.5*(R-R0)**2/Deltar**2) * np.exp(-0.5*((Phi-phi1))**2/Deltaphi**2)
    gauss2 = np.exp(-0.5*(R-R0)**2/Deltar**2) * np.exp(-0.5*((Phi-phi2))**2/Deltaphi**2)
    
    krad = np.argmin(np.abs(field.grid.ri - R0))
    field.array = field.array[krad,0] * (gauss1 + gauss2)
    
    field.save()
    field.save("reference")

if __name__ == "__main__":
    main()