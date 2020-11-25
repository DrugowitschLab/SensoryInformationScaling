# Information Scaling in Large Neural Populations

All the simulation data generated uing `simPopActivity` function will be stored in this folder with a name in this format: `sim[simid]` where `simid` is the integer representing the simulation id and ranges from 1 to 4. Simulation details for different `simid` numbers are expressed as follows:
- `simid = 1` - Gaussian population with limited information, $$N=1000$$, $$T=1000$$, $$I_\infty=20$$
- `simid = 2` - Gaussian population with unlimited information,  $$N=1000$$, $$T=1000$$, $$I_\infty=\infty$$
- `simid = 3` - LNP population with limited information, $$N=2500$$, $$T=2500$$
- `simid = 4` - LNP population with unlimited information, $$N=2500$$, $$T=2500$$

See Sec. 4.1 and Sec 4.2 of the SI for more details on the multivariate Gaussian population model and LNP model. 