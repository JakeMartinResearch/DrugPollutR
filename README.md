# PollutionScopeR

PollutionScopeR is an R Shiny application for exploring aquatic environmental occurrence data and using those data to build environmentally grounded dose series for ecotoxicology experiments.

Hosted app: [https://jakemartinresearch.shinyapps.io/PollutionScopeR/](https://jakemartinresearch.shinyapps.io/PollutionScopeR/)

## Background

Pharmaceutical pollution is a growing global environmental concern, and experimental studies often test concentrations that are much higher than those typically reported in real-world aquatic systems. A major practical reason for this mismatch is that environmental occurrence data are large, fragmented, inconsistently indexed, and time-consuming to clean before they can be used for study design.

PollutionScopeR was built to reduce that friction. It provides a simple interface for searching a harmonized occurrence dataset, summarizing observed concentrations, visualizing their distribution, and translating those observations into candidate dose series for laboratory studies. The app is intended as a decision-support tool that helps researchers justify exposure concentrations with reference to measured environmental data.

## What The App Does

PollutionScopeR is organized into two linked modules.

### ChemicalExploreR

ChemicalExploreR helps users find and summarize occurrence records for a compound of interest.

Features:

- Search by compound name, CAS number, or CID
- Select one or more matching search terms
- Group results by `parent_compound` when related forms should be treated together
- Filter records by matrix type: `surfacewater`, `effluent`, or `Combined`
- View matched compounds, summary statistics, concentration density plots, and filtered records
- Download filtered datasets and summaries with an embedded citation note

Summary statistics and plots are calculated from positive environmental detections only.

### DoseSelectR

DoseSelectR uses the environmental concentration data returned from the search workflow to propose a dose series for ecotoxicology studies.

Features:

- Choose the central tendency used to anchor the series: `median`, `mean`, or `mode`
- Set the target number of doses
- Adjust the dose spacing factor
- Generate a proposed dose table and visualize it against the observed distribution on raw and log scales
- Download the resulting dose table

The current workflow uses only positive concentration values, requires enough data to define a usable distribution, and applies a lower quantification threshold to avoid implausibly low doses.

## Data Foundation

This database was compiled for Martin et al (2025) Environ. Sci. Technol. Lett. 2025, 12, 10, 1308-1313 (https://doi.org/10.1021/acs.estlett.5c00665). It is based on a filtered synthesis of three publicly available datasets: (1) the NORMAN EMPODAT database for chemical occurrence (accessed 18/03/2025), (2) the Umweltbundesamt Pharmaceuticals in the Environment database (PHARMS-UBA; accessed 19/12/2024), and (3) Wilkinson et al. (2022) Pharmaceutical Pollution of the World's Rivers database. Data were restricted to entries reported in mass per volume of water (e.g., µg/L) and relevant to surface water and wastewater matrices (for details of the filtering process, please refer to Martin et al. (2025)

In the app, the processed datasets are stored locally as:

- `env_data.rds`
- `compound_key.rds`

## Repository Structure

- `app.R`: main Shiny application
- `DoseSelectR.R`: helper function for dose selection and plotting
- `env_data.rds`: processed environmental occurrence data used by the app
- `compound_key.rds`: search key used for compound lookup

## Running Locally

### Requirements

Install a recent version of R and the packages used by the app:

- `shiny`
- `bslib`
- `shinyWidgets`
- `dplyr`
- `ggplot2`
- `data.table`
- `DT`
- `scales`
- `patchwork`

### Launch

From the project directory in R:

```r
shiny::runApp()
```

Or open `app.R` in RStudio and run the app.

## Typical Workflow

1. Search for a compound in ChemicalExploreR by name, CAS, or CID.
2. Review matched compounds and choose whether to group by parent compound.
3. Filter by matrix type and inspect summary statistics, plots, and underlying records.
4. Move to DoseSelectR and choose a central tendency, number of doses, and spacing factor.
5. Generate a candidate dose series and download the outputs if needed.

## Intended Use And Scope

PollutionScopeR is designed to support experimental planning, not to replace researcher judgment. Users may still decide to test concentrations above currently observed environmental levels, include hotspot scenarios, or add positive controls for specific experimental aims. The app is best used as a transparent starting point for environmentally informed study design.

## Citation

If you use the app or its compiled data in research outputs, please cite the associated paper and the original occurrence data sources:

- Martin, J. M.; Brand, J. A.; McCallum, E. S. (2025). *Aligning Behavioral Ecotoxicology with Real-World Water Concentrations: Current Minimum Tested Levels for Pharmaceuticals Far Exceed Environmental Reality*. Environmental Science & Technology Letters. [https://doi.org/10.1021/acs.estlett.5c00665](https://doi.org/10.1021/acs.estlett.5c00665)
- NORMAN EMPODAT
- PHARMS-UBA
- Wilkinson et al. (2022). *Pharmaceutical Pollution of the World's Rivers*. Proceedings of the National Academy of Sciences. [https://doi.org/10.1073/pnas.2113947119](https://doi.org/10.1073/pnas.2113947119)

## Acknowledgement

This README was drafted from the app source code and the manuscript background in `sandpit/PollutionScopeR_submit.docx`, then condensed for repository use.
