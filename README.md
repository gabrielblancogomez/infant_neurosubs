# Neurosubtyping infants using EEG

This study examines EEG measures in infants (n=144) from the EEG-IP dataset, using both latent profile analysis (LPA) and hierarchical clustering (HC) to identify distinct neurophysiological subtypes. Three distinct classes were identified through each method, showing different language development trajectories over time. The following code is used to extract EEG features, defines clusters and does statistical analyses to test differences in language directories. 




## Data summary visualization
A visual summary of all graphs and analyzes can be found here:[Infant EEG Dashboard](https://https://gabrielblancogomez.github.io/infant_neurosubs_dashboard): 

- **Study Overview**: Summary of sample demographics and analytical methods.
- **Class Identification**: EEG-based subtypes identified using LPA and HC.
- **Feature Distributions**: How EEG features differ across classes.
- **UMAP Projections**: Low-dimensional visualization of EEG feature clustering.
- **Language Trajectories**: Longitudinal language development outcomes by class.
- **36-Month Language Outcomes**: Final language scores by EEG subtype.

**Key Finding:**  
Early EEG measures can predict distinct language development trajectories, highlighting potential neurophysiological markers of developmental outcomes.

---
## Project structure

| File/Script | Description |
|:------------|:------------|
| `Feature_extraction.ipynb` | Extracts EEG features (gamma power, frontal power, Connectivity Auditory Network and Connectivity Speech Network. |
| `Hierarchical_clustering_SKklearn.ipynb` | Performs hierarchical clustering on EEG features to identify distinct neurophysiological subgroups. |
| `NBClust_selection.R` | Determines the optimal number of clusters for EEG data using various statistical indices. |
| `Linear_Mixed_effects_modelling_language.R` | Analyzes language development trajectories over time using linear mixed-effects models. |
| `Visualization_brain_regions.R` | Generates brain region visualizations for extracted EEG measures. |
| `Compare_HC_LPA.ipynb` | Compares LPA and HC clustering results and examines agreement between the two methods. |
| `Tables_and_Figures.ipynb` | Compiles final tables and figures used in the dashboard and manuscript. |

## Roadmap


  1. Clone this repository.
  2. Run `Feature_extraction.ipynb` to generate EEG feature datasets.
  3. Apply clustering with `Hierarchical_clustering_SKklearn.ipynb` and `NBClust_selection.R`.
  4. Analyze language outcomes using `Linear_Mixed_effects_modelling_language.R`.
  5. Visualize outputs with `Visualization_brain_regions.R` and `Tables_and_Figures.ipynb`.

## Authors

- Gabriel Blanco Gomez [@gabrielblancogomez](https://www.github.com/gabrielblancogomez)


## EEG Data access

Data used in all these analyses comes from this study and its available upon reasonable request: 
 
van Noordt, S., Desjardins, J.A., Huberty, S. et al. EEG-IP: an international infant EEG data integration platform for the study of risk and resilience in autism and related conditions. Mol Med 26, 40 (2020). https://doi.org/10.1186/s10020-020-00149-3


## Acknowledgements

The authors would like to thank all families who took part in this study, as well as the mentors who helped with this project. Special thanks to Myriam Beauchamps, James Desjardins, Diksha Srishyla, Scott Huberty, Julie Scorah and Shoi Shi.

## Badges



[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)


## Requirements
- Python 3.8+
- R 4.0+
- Packages:
  - Python: `numpy`, `pandas`, `scikit-learn`, `plotly`, `umap-learn`
  - R: `NbClust`, `lme4`, `ggplot2`

## Support

For support, email gabriel.blancogomez@mail.mcgill.ca


## License

[MIT](https://choosealicense.com/licenses/mit/)

