# Flood Risk Prediction tool

## Deadlines
-  *code 12pm GMT Friday 25th November 2022*
-  *presentation/ one page report 4pm GMT Friday 25th November 2022*

### Key Requirements

Your project must provide the following:

 - at least one analysis method to estimate a number of properties for unlabelled postcodes extrapolated from sample data which is provided to you:
    - Flood risk (on a 10 point scale).
    - Median house price.
 - at least one analysis method to estimate the Local Authority & flood risk of arbitrary locations.
 - a method to find the rainfall and water level near a given postcode from provided rainfall, river and tide level data, or by looking online.

 You should also provide visualization and analysis tools for the postcode, rainfall, river & tide data provided to you, ideally in a way which will identify potential areas at immediate risk of flooding.
 
 Your code should have installation instructions and basic documentation, either as docstrings for functions & class methods, a full manual or both.

![London postcode density](images/LondonPostcodeDensity.png)
![England Flood Risk](images/EnglandFloodRisk.png)
![UK soil types](images/UKSoilTypes.png)

This README file *should be updated* over the course of your group's work to represent the scope and abilities of your project.

### Assessment

 - your code will be assessed for its speed (both at training and prediction) & accuracy.
 - Your code should include tests of its functionality.
 - Additional marks will be awarded for high code quality and a clean, well organised repository.

 ### Installation Guide

Packages required to install that are not in environment.yml: Pyqt5
   Note: PyQtWebEngine package does not work properly with Apple M1 chips. Instead, use Mac with an Intel CPU.

Added python files/ folders:

   - preprocessing_tool: contains functions used for preprocessing.
   - Extra (folder): contains unused modules (see docstrings in files for more information).

### User instructions

**Modelling (predictions):**

   risk label prediction method: {0: 'all_zero_risk', 1: knn_model}

   house price prediction method: { 0: 'all_england_median',
                                    1: knn_model,
                                    2: ridge_model,
                                    3: lasso_model}

   local authority prediction method: { 0:'Do Nothing',
                                       1: knn_model,
                                       2: decisiontree_model,
                                       3: random_forest_model}

   *Note: method = 1 is the best one in all methods*
   More examples for running the models can be found in trainModels.ipynb.

**Data Visualisation(GUI Design):**

   The package used to create the GUI is PyQt5 (requires additional external packages). However, note that the PyQtWebEngine package does not work properly with Apple M1 chips. Instead, use Mac with an Intel CPU.

   There are two graphic user interfaces located in the main folder: ‘CoordinateTool.py’ and ‘GraphView.py’.
   Moreover, there are two folders located in the main folder, used to provide input data: ‘coordinate_file’ and ‘result_file’.

   ‘CoordinateTool.py’: the GUI converts latitude and longitude to easting and northing (and vice versa). It takes as input:

   - a single pair of values, or

   - a list of pairs of values.
         Note: this is done by uploading a csv file. For example, if ‘example.csv’ contains all the pairs of values, type ‘example’ in the input box. The user has to make sure the csv file only contains two columns (‘Latitude’ and ‘Longitude’, or ‘Easting’ and ’Northing’). The csv file should be stored in the ‘coordinate_file’ folder. There is an example in the folder to play with.

   ‘GraphView.py’: The GUI displays a data visualisation graph. The user should input a csv file (stored in the ‘result_file’ folder), with columns: ‘Latitude’, ‘Longitude’, ‘Risk Label’, ‘Property Median Price’, and ‘Risk Level’. The user can then choose different options to plot the graph they want.
         Note1: Due to the GUI power limitations, the GUI cannot deal with high-density data (it displays a white image). Therefore, the GUI displays at most 1000 data points in each plot.
         Note2: The plotting API has more functionalities that the GUI, see details in jupyter notebook (DataVisualization.ipynb).

   To run GUIs, use commands:

      python GraphView.py
      python CoordinateTool.py

### Documentation

Documentation can be found in the docs folder (index.html).

### Testing

The tool includes several tests, which you can use to check its operation on your system. With [pytest](https://doc.pytest.org/en/latest) installed, these can be run with

```bash
python -m pytest --doctest-modules flood_tool
```

### Reading list

 - (A guide to coordinate systems in Great Britain)[https://webarchive.nationalarchives.gov.uk/20081023180830/http://www.ordnancesurvey.co.uk/oswebsite/gps/information/coordinatesystemsinfo/guidecontents/index.html]

 - (Information on postcode validity)[https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/283357/ILRSpecification2013_14Appendix_C_Dec2012_v1.pdf]
