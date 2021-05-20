
# FOSLCs 
Author: Takahiro Suzuki

lastpudate: 2021-05-20

---
Analysis scripts used for Yoshino et al. "Embryo-free reconstitution of ovarian follicles from mouse pluripotent stem cells# ### (2021)

### CellRanger
- cell_ranger_job_script.sh
script for cell rangerbanalysis
- cell_ranger_job_script_queuing.sh
queuing code for the cell_ranger_job_script.sh

### R
- Functions.r
functions for the analysis. This should firstly be loaded in advance of the R main analysis.
- all_ells.r
A script for the all sample (E10, E11, E12, E13, E14, induced D6) analysis.
- E12_D6.r
A script for the comaprison between E12 and induced D6.
