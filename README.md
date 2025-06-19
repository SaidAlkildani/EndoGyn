# EndoGyn

![Endogyn project image](images/Endogyn_banner_bg_1200px.jpg)

An explorative Data Science and AI capstone project, where ovarian endometriosis and ovarian cancer are explored using machine learning. 

## Background

**Endometriosis** is an understudied and painful condition impacting 10% of women globally [1]. 
* Research suggests an average of 8 years until diagnosis due to lack of disease knowledge [2]. It is defined by the uterine tissue growth outside of the uterus. Additionally, those affected by endometriosis are at a higher risk of gynecologic cancers as endometriosis shares many of the cancer's characteristics. 

**Ovarian cancer** is the most common endometriosis-associated cancers [3], as well as the most lethal across gynecologic cancers. Many endometriosis hallmarks are shared with cancer hallmarks: angiogenesis, proliferation, evading the immune system, adhering, metastsis (according to one of the two pathogensis hypotheses) [4]. This raises the question of whether there are conserved, as well as diverging, molecular programs at play. 
## Data
| Condition| Dataset | Size | Patients | Cell Count| Cell Count post ETL |
| -------- | ------- | ------- | ------- | ------- | ------- |
| Ovarian Endometriosis  |  GEO214411 | 1.5 gb | 6 endometriosis, 7 controls | 155361| 61031|
| Ovarian Cancer | GSE184880  | 0.5 gb | 7 cancer, 5 controls | 65820 | 32520|

## Methodologies
### Data Allocation and Preprocessing
### Data Integration
### Machine Learning 

## Requirements
* pyenv with Python: 3.11.3
* Further package requirements are included as .txt files. To set up the virtual environment for the ETL stage:
### **`macOS`** type the following commands : 
- Install the virtual environment and the required packages by following commands:

    ```BASH
    pyenv local 3.11.3
    python3 -m venv .venv_ETL
    source .venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements_ETL.txt #choose the correct .txt file
    ```
### **`Windows`** type the following commands :

- Install the virtual environment and the required packages by following commands.

   For `PowerShell` CLI :

    ```PowerShell
    pyenv local 3.11.3
    python -m venv .venv_ETL
    .venv\Scripts\Activate.ps1
    python -m pip install --upgrade pip
    pip install -r requirements_ETL.txt #choose the correct .txt file
    ```

    For `Git-bash` CLI :
  
    ```BASH
    pyenv local 3.11.3
    python -m venv .venv_ETL
    source .venv/Scripts/activate
    python -m pip install --upgrade pip
    pip install -r requirements_ETL.txt #choose the correct .txt file
    ```

## References
[1] Zondervan KT, Becker CM, Missmer SA. Endometriosis. N Engl J Med. 2020 Mar 26;382(13):1244-1256. doi: 10.1056/NEJMra1810764. PMID: 32212520.

[2] Vincent K, Horne AW. Current Pathways of Care Continue to Fail Those With a Diagnosis of Endometriosis. BJOG. 2025 Jun 16. doi: 10.1111/1471-0528.18245. Epub ahead of print. PMID: 40524487.

[3] Krawczyk N, Banys-Paluchowski M, Schmidt D, Ulrich U, Fehm T. Endometriosis-associated Malignancy. Geburtshilfe Frauenheilkd. 2016 Feb;76(2):176-181. doi: 10.1055/s-0035-1558239. PMID: 26941451; PMCID: PMC4771509.

[4] Siufi Neto J, Kho RM, Siufi DF, Baracat EC, Anderson KS, Abrão MS. Cellular, histologic, and molecular changes associated with endometriosis and ovarian cancer. J Minim Invasive Gynecol. 2014 Jan-Feb;21(1):55-63. doi: 10.1016/j.jmig.2013.07.021. Epub 2013 Aug 17. PMID: 23962574.
