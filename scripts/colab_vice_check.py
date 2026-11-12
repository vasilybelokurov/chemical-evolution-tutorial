# Paste into ONE Colab cell (fresh runtime) and run. Prints timing, versions and a smoke result.
import subprocess, sys, time, platform
t0 = time.time()
# VICE 1.3.1 is sdist-only and fails to Cythonize with Cython >= 3.1 ('long' builtin removed),
# so pin Cython < 3.1 and build without pip's isolated build environment.
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "cython>=3.0,<3.1", "setuptools"], check=True)
r = subprocess.run([sys.executable, "-m", "pip", "install", "--no-build-isolation", "vice==1.3.1"],
                   capture_output=True, text=True)
print("pip exit", r.returncode, f"| install {time.time()-t0:.0f} s | python {platform.python_version()} | {platform.platform()}")
if r.returncode: print(r.stdout[-3000:], r.stderr[-3000:])

import os, tempfile, numpy as np, vice
os.chdir(tempfile.mkdtemp())
sz = vice.singlezone(name="smoke", elements=("fe", "mg"), mode="ifr",
                     func=lambda t: 9.1, tau_star=2.0, eta=2.5, dt=0.01)
sz.run(np.linspace(0, 10, 201), overwrite=True)
h = vice.history("smoke")
feh, mgfe = np.array(h["[fe/h]"]), np.array(h["[mg/fe]"])
print(f"vice {vice.__version__}  [Fe/H]_end={feh[-1]:.3f}  [Mg/Fe] {mgfe[1]:.3f}->{mgfe[-1]:.3f}")
# macOS reference (Py 3.12 and 3.13): [Fe/H]_end=-0.184, [Mg/Fe] 0.567->-0.448
assert abs(feh[-1] + 0.184) < 0.002 and abs(mgfe[-1] + 0.448) < 0.002
print("SMOKE OK")
