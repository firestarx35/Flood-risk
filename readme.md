# Flood Risk Prediction Tool

## Overview

Flood risk is a significant concern, affecting both global and UK-based human settlements. Predicted climate changes suggest an increase in the frequency and severity of floods. Floods can occur from various sources such as surface water, rivers, and seas. To better assess and manage these risks, we have developed an efficient tool leveraging Machine Learning (ML) techniques.

![Flooding in York, UK](images/York_Floods_2015.jpg)
<p align="center"><i> Picture taken by <a href="https://flickr.com/photos/60539443@N05/23675629939">J. Scott</a>, used under <a href="https://creativecommons.org/licenses/by/2.0/deed.en">CC-BY</a></i></p>

For a detailed synopsis and definitions, please refer to the `ProjectDescription.ipynb` file.

## Features of the Tool

- Based on an arbitrary postcode in London, the tool can provide:
  - Flood risk on a 10 point scale.
  - Estimate of the median house price.
- The tool estimates the Local Authority & flood risk for arbitrary locations.
- It can find the rainfall and water level near a given postcode using provided rainfall, river and tide level data, or by fetching the data from the [Met Office API](http://environment.data.gov.uk/flood-monitoring/data/readings.csv?).

## Installation

This application requires Python 3.9 and above. Create a new Conda environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
```

Note: PyQtWebEngine package does not work properly with Apple M1 chips. If possible, use a Mac with an Intel CPU.

How it Works
The tool employs Pandas for preprocessing, Scikit-learn for prediction, Matplotlib and Folium for visualisation, and PyQt5 for the user interface.

Below are some links to visualisation examples:

Flood Risk Label and Tide
Flood risk region
England Flood Risk

Apologies for the confusion. Below is the complete markdown for the readme.md file:

markdown
Copy code
# Flood Risk Prediction Tool

## Overview

Flood risk is a significant concern, affecting both global and UK-based human settlements. Predicted climate changes suggest an increase in the frequency and severity of floods. Floods can occur from various sources such as surface water, rivers, and seas. To better assess and manage these risks, we have developed an efficient tool leveraging Machine Learning (ML) techniques.

![Flooding in York, UK](images/York_Floods_2015.jpg)
<p align="center"><i> Picture taken by <a href="https://flickr.com/photos/60539443@N05/23675629939">J. Scott</a>, used under <a href="https://creativecommons.org/licenses/by/2.0/deed.en">CC-BY</a></i></p>

For a detailed synopsis and definitions, please refer to the `ProjectDescription.ipynb` file.

## Features of the Tool

- Based on an arbitrary postcode in London, the tool can provide:
  - Flood risk on a 10 point scale.
  - Estimate of the median house price.
- The tool estimates the Local Authority & flood risk for arbitrary locations.
- It can find the rainfall and water level near a given postcode using provided rainfall, river and tide level data, or by fetching the data from the [Met Office API](http://environment.data.gov.uk/flood-monitoring/data).

## Installation

This application requires Python 3.6 and above. Create a new Conda environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
```

#### Note: PyQtWebEngine package does not work properly with Apple M1 chips. If possible, use a Mac with an Intel CPU.

## How it Works
The tool employs Pandas for preprocessing, Scikit-learn for prediction, Matplotlib and Folium for visualisation, and PyQt5 for the user interface.

Below are some links to visualisation examples:

![Flood Risk Label and Tide](images/Flood_labels.png)
![Flood risk region](images/flood_prediction.png)
![England Flood Risk](images/EnglandFloodRisk.png)
![UK soil types](images/UKSoilTypes.png)

## File/Folder Structure

preprocessing_tool: Contains functions used for preprocessing.
Extra (folder): Contains unused modules. See docstrings in files for more information.

## User Guide
Modelling (predictions)
Models are available for risk label, house price, and local authority prediction. Each model offers a few methods. Please note that method = 1 is recommended for all models. For more examples of running the models, refer to `trainModels.ipynb`.


## Data Visualisation (GUI Design)
The package used to create the GUI is PyQt5 (requires additional external packages). Remember, PyQtWebEngine package does not work properly with Apple M1 chips. Use a Mac with an Intel CPU instead.

There are two GUI tools available in the main folder: `CoordinateTool.py` and `GraphView.py`, and two folders used to provide input data: coordinate_file and result_file.

#### Detailed usage for CoordinateTool.py and GraphView.py is provided in their respective sections.

#### To run GUIs, use the following commands:

```bash
python GraphView.py
python CoordinateTool.py
```

## Documentation and Report
Documentation can be found in the docs folder (see index.html). A project report is also located in the docs folder.

## Testing
The tool includes several tests to verify its functionality. Ensure you have pytest installed, and run the tests with the following command:

```bash
python -m pytest --doctest-modules flood_tool
```

## Reading List
A Guide to Coordinate Systems in Great Britain
Information on Postcode Validity