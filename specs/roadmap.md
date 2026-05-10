# Project Roadmap

## Milestone 1: Robustness and Stability (Current Focus)
The immediate priority is to ensure the existing feature set is production-ready and accessible to a broad audience.
1. **Ease of Installation:** Streamline the PyPI and source installation processes across different OS and hardware (e.g., Apple Silicon).
2. **Automated Testing:** Expand the test suite (unit and integration tests) to cover edge cases in low-coverage deconvolution.
3. **Stress Testing:** Benchmark the performance of LLSE and MMSE on large-scale datasets and diverse atlas configurations.
4. **CLI Polishing:** Improve error messaging, help documentation, and user feedback during long-running deconvolution tasks.

## Milestone 2: Handling cfDNA Biases
Refine the models to better handle the biological and technical nuances of cell-free DNA.
1. **End-Repair Bias Mitigation:** Develop methods to account for the loss of methylation signal near read ends caused by library preparation.
2. **MNase Digestion Modeling:** Incorporate biases introduced by nucleosome-aware digestion into the deconvolution framework.
3. **Read-Length Awareness:** Adjust error parameters dynamically based on read length and mapping quality.

## Milestone 3: Support for Unknown Cell Types
Move beyond closed-reference deconvolution to support more realistic clinical scenarios.
1. **Latent Cell Type Detection:** Allow the model to account for signal that does not match any reference in the atlas.
2. **Atlas Refinement:** Implement iterative updates to reference frequencies based on observed sample data (similar to `CelFiE`).

## Milestone 4: Multi-Omic Integration
Expand the utility of Nanomix by integrating other data modalities.
1. scRNA-seq Pairing: Use single-cell RNA sequencing to derive more accurate and context-specific methylation signatures.
2. Dynamic Atlas Construction: Automate the selection of differentially methylated regions (DMRs) from diverse public databases.
3. Clinical Validation: Large-scale evaluation on clinical liquid biopsy cohorts for cancer monitoring and early detection.

## Milestone 5: Web Interface and Visualization
Translate the power of the Nanomix CLI into a user-friendly web application for broader accessibility.
1. **Interactive Dashboard:** Develop a web-based interface for uploading methylome files and selecting deconvolution parameters.
2. **Visual Analytics:** Implement interactive charts (e.g., stacked bar plots, heatmaps) for exploring mixture proportions and fragment assignments.
3. **Cloud Deployment:** Provide containerized deployment options (e.g., Docker) and hosted versions for institutional use.
4. **API Layer:** Create a RESTful API to allow external services to programmatically invoke Nanomix deconvolution models.

