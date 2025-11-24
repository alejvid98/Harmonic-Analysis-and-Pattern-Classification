# Harmonic Analysis and Pattern Classification of Electrocardiograms for Heart Disease Diagnosis


## Manuscript ID: 10091  
**Authors:**


- Alejandro Vidales Esquivel: He received the B.Sc. degree in Electronic Engineering from Universidad Michoacana de San Nicolas de Hidalgo (UMSNH), Morelia, Mexico, in 2021, 
and the M.Sc. degree in Electrical Engineering from UMSNH in 2025. His research focuses on harmonic analysis of ECG signals, development of state observers for system estimation,
biomedical signal modeling, and the application of machine learning techniques for computer-based medical diagnosis.
- Fernando Ornelas Tellez: He received the B.Sc. degree from the Instituto Tecnologico de Morelia (ITM), Morelia, Mexico, in 2005 and the M.Sc. and D.Sc. degrees in electrical engineering from the Advanced Studies and Research Center, National Polytechnic Institute (CINVESTAV-IPN),Guadalajara, Mexico, in 2008 and 2011, respectively.
Since 2012, he has been with the Universidad Michoacana de San Nicolas de Hidalgo, where he is currently a professor of Electrical Engineering graduate programs. 
His research interest centers on optimal control, neural control, sliding modes control, and passivity, and their applications to smart grids, power electronics, mechanical systems, and electrical machines.
- Jose Ortiz Bejar: He received the Ph.D. degree in data science from INFOTEC, Unidad Aguascalientes, in 2020. He is currently an associate professor with the Michoacan University of San Nicolás de Hidalgo. He has served as an Organizer of ROPEC and IEEE T\&DLA conferences. 
His research interest centers on machine learning applied to power systems analysis, text categorization, and clustering.

--

## 📁 Files

The following files are used to obtain the research results.

| File |  
|------|
| example estimator.nb | 
| base_datos.csv |
| main.m| 
| armonicosv2.slx |
| clasificador_armonicos.slx |
| knnSFunction.m|
| predictKNN.m |
| knnModel.mat |

### 📝 Description of Files

- **example estimator.nb:** Mathematica file that includes an example of the harmonic estimator used in the methodology.
- **base_datos.csv:** Contains the database used in the study, specifically the amplitudes of the harmonics extracted from ECG signals.
- **main.m:** Main execution script in MATLAB where the system parameters are defined, and the observer gains for the estimator are calculated. This script is responsible for controlling the workflow, calling the required Simulink models, interpreting the results, and finally determining the pathology.
- **armonicosv2.slx:** Simulink model used for harmonic estimation. This model must be opened before executing `main.m`, but no manual intervention is required inside the model, as its execution is fully controlled from `main.m`.
- **clasificador_armonicos.slx:** Simulink model responsible for classifying the ECG signal based on the estimated harmonic values. This model also needs to be open before execution but does not require manual actions. It receives the estimation results from the previous model and outputs a classification label.
- **knnSFunction.m:** This file is required by the `clasificador_armonicos.slx` model. It acts as a custom Simulink function that connects the harmonic data with the trained KNN classifier. No manual configuration is needed — the function works automatically during model execution.
- **predictKNN.m:** This file is used by the `knnSFunction.m` block inside the `clasificador_armonicos.slx` model. It loads the trained classifier stored in `knnModel.mat` and receives the harmonic values as input. Then, it processes the data and returns the final prediction label based on the trained KNN model. No modification or manual interaction is required — the function runs automatically as part of the Simulink classification process.
- **knnModel.mat:** This file contains the pre-trained K-Nearest Neighbors model.  



During execution, `main.m` first runs `armonicosv2.slx` to obtain the harmonic estimation. Once completed, the results are automatically passed to `clasificador_armonicos.slx`, which outputs a numeric label corresponding to a pathology:

- **0 → Healthy**
- **1 → Arrhythmia**
- **2 → Myocardial Infarction**
- **3 → Heart Failure**

Finally, `main.m` interprets the label and outputs the corresponding diagnosis in text format.

To run the diagnosis process, the ECG signals to be analyzed must be placed inside a folder named **`data`**.  
The system reads the signals automatically from this folder, so maintaining the folder name and structure is required.  
The dataset (`base_datos.csv`) is **not needed for normal execution** — it is only included in case the user wishes to retrain or test new classifiers.


