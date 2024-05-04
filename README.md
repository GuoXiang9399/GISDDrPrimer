# GISDDrPrimer

GISDDrPrimer is a Shiny R based application that allows searching for published/novel primer pairs targeting the DENV fragment and visualizing their alignment to the reference genome.


GISDDrPrimer shows the amplicon and any variation found in the amplicon as well as in the primer binding regions. GISDDrPrimer also provides a list of in-house designed conserved and degenerate primer pairs across the viral genome and presents information on occurrence and lineage of mutations with the degenerate codes, alignment to the genome.

## Primer Bank

<img src="Figure/PrimerBank.png" width="100%" style="display: block; margin: auto;" />

### Features  
  
- **Method Selection**: Users can choose from "All", "Detection", "Sanger sequence", or "High through sequence".  
- **Serotype Selection**: Users can choose from "All", "Universal", "D1", "D2", "D3", or "D4".  
- **Target Selection**: Users can choose from "All", "C-prM", "E", "whole genome", or "Others".  
- **Data Table Output**: Displays the data that matches the selected criteria.  
- **Download Functionality**: Allows users to download the filtered data.  

## Primer Evaluation

<img src="Figure/PrimerEvaluation.png" width="100%" style="display: block; margin: auto;" />

### Features  
  
- **Primer Input**: Users can enter forward and reverse primer sequences in the text boxes.  
- **Serotype Selection**: Users can choose from "All", "Universal", "D1", "D2", "D3", or "D4".  
- **Primer Evaluation**: Displays the primer evaluation results associated with the selected primers and serotype.  

## Usage  
  
To run this application, you need to have R and Shiny installed on your system. Save the code as an R script file (e.g., `app.R`) and run it in an R environment. This will launch the Shiny application locally, and you can access it through a web browser.  
  
## Contributions and Feedback  
  
If you have any suggestions, feedback, or want to contribute code, please contact us through the [GitHub repository link] (replace with the actual GitHub repository link).  