# The probabilities of Clovis atlatl use

This repository contains the manuscript, R code, source data, Stan models, saved posterior draws, to reproduce results of **The probabilities of Clovis atlatl use**
by *Erik Otárola-Castillo* and *Benjamin J. Schoville*. It follows the reproducible R Markdown and Stan framework used in Otárola-Castillo and Torquato (2018) and Otárola-Castillo et al. (2020, 2022, 2023).

The archaeological approach utilized here builds on earlier Bayesian examples for assigning hafting probabilities to observed projectile points. These examples were outlined by Otárola-Castillo and Torquato (2018) and Otárola-Castillo et al. (2023). Additional routines for working with posterior draws from models fitted in Stan through the R brms package, from scratch - so to speak, appear in Perri et al. (2019).

Like in the references above, the code presented here is not optimized to run as efficient software, rather to present the numerical work in a progressive and logical manner.

## Run in RStudio

Open `Clovis_atlatl_probabilities.Rproj` in RStudio.

To create the manuscript, open `manuscript/letter.Rmd` and select Knit to the desired format
from the Knit menu.

To create a document out of analysis code, open `code.Rmd` and click Knit.
The default `report` mode reads `results/frozen/run1/` and rebuilds
the tables and figures from the frozen results, without refitting models or recalculating results.

To fit the models again, change `mode` at the top of `code.Rmd` to `refit`
and click Knit. To only recalculate from saved draws, use `recalculate`.

New results go to `results/refit/` or `results/recalculate/`, with `run_id`
as the subfolder when supplied. Reusing an output folder replaces its files.

## Repository files

* `data/`: Clovis and reference measurements, radiocarbon records, calibration
data, priors, and model settings.
* `models/`: Stan models for chronology, Clovis morphology, and reference morphology.
* `scripts/`: Includes the main controller, `code.Rmd`, code for analysis and report, and helper functions for manuscript tables and figures.
* `results/`: frozen results, posterior draws, diagnostics, and folders for refits and recalculations.
* `manuscript/`: letter source, Word template, rendered drafts, and publication figures and tables in its `figures/` and `tables/` folders, document format reference and bibliography, and citation style in its `formatting/` folder.

## References Cited

* **Otárola‐Castillo, E.R.**, Torquato, M.G. and Buck, C.E., 2023. The Bayesian inferential paradigm in archaeology. Handbook of Archaeological Sciences, 2, pp.1193-1209. https://github.com/dreoc/the-bayesian-inferential-paradigm-in-archaeology
* **Otárola‐Castillo, E.R.**, Torquato, M.G., Wolfhagen, J., Hill Jr, M.E. and Buck, C.E., 2022. Beyond chronology, using Bayesian inference to evaluate hypotheses in archaeology. Advances in Archaeological Practice, 10(4), pp.397-413. https://github.com/dreoc/Bayesian-Archaeology-Beyond-Chronology
* **Otárola‐Castillo, E.R.**, Torquato, M.G. and Hill, M.E., 2020. Intensification mechanisms driving dietary change among the Great Plains big game hunters of North America. Defining and Measuring Diversity in Archaeology, pp.125-158. https://github.com/dreoc/diet-diversity-mechanisms-paper
* Perri, A.R., Koster, J.M., **Otárola‐Castillo, E.R.**, Burns, J.L. and Cooper, C.G., 2019. Dietary variation among indigenous Nicaraguan horticulturalists and their dogs: An ethnoarchaeological application of the Canine Surrogacy Approach. Journal of Anthropological Archaeology, 55, p.101066. https://osf.io/dbxe9/overview
* **Otárola‐Castillo, E.R.** and Torquato, M.G., 2018. Bayesian statistics in archaeology. Annual Review of Anthropology, 47(1), pp.435-453. https://doi.org/10.1146/annurev-anthro-102317-045834

