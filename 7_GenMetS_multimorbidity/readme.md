
## Traits available in UKBiobank and used in this study.
| FieldID | Field                                      | ValueType           | Units         | Trait In This Study            |
|---------|--------------------------------------------|---------------------|---------------|---------------------------------|
| 31      | Sex                                        | Categorical single  |               | Sex                             |
| 34      | Year of birth                              | Integer             | years         | Year of birth                   |
| 48      | Waist circumference                        | Continuous          | cm            | WC                              |
| 93      | Systolic blood pressure, manual reading    | Integer             | mmHg          | SBP manual                      |
| 94      | Diastolic blood pressure, manual reading   | Integer             | mmHg          | DBP manual                      |
| 1070    | Time spent watching television (TV)        | Integer             | hours/day     | TV hours                        |
| 1160    | Sleep duration                             | Integer             | hours/day     | Sleep hours                     |
| 4079    | Diastolic blood pressure, automated reading| Integer             | mmHg          | DBP*                            |
| 4080    | Systolic blood pressure, automated reading | Integer             | mmHg          | SBP*                            |
| 6138    | Qualifications                             | Categorical multiple|               | Education                       |
| 20003   | Treatment/medication code                  | Categorical multiple|               | Medication                      |
| 21000   | Ethnic background                          | Categorical single  |               | Self-reported Ethnicity         |
| 21001   | Body mass index (BMI)                      | Continuous          | kg/m2         | BMI                             |
| 21003   | Age when attended assessment centre        | Integer             | years         | Age                             |
| 22001   | Genetic sex                                | Categorical single  |               | Genetic sex                     |
| 22189   | Townsend deprivation index at recruitment  | Continuous          |               | Socioeconomic status            |
| 22040   | Summed MET minutes per week for all activity| Continuous          | minutes/week | Physical Activity               |
| 30750   | Glycated haemoglobin (HbA1c)              | Continuous          | mmol/mol      | HbA1c                           |
| 30760   | High Density Lipoprotein (HDL)            | Continuous          | mmol/L        | HDL                             |
| 30780   | Low Density Lipoprotein (LDL)             | Continuous          | mmol/L        | LDL                             |
| 30870   | Triglycerides (TG)                        | Continuous          | mmol/L        | TG                              |
| 41270   | Diagnoses - ICD10                         | Categorical multiple|               | ICD10                           |
| 41280   | Date of first in-patient diagnosis - ICD10| Date                |               | ICD10 disease onset             |

*Automated DBP/SBP readings were the preferred source of data. However, if unavailable, manual readings were used.


## ICD10 codes for disease diagnoses from electronic record in UK Biobank used in the multi-morbidity analysis.

| Diseases                            | Case_ICD10_codes                                                               | non-control_ICD_codes                                                                                                                                     |
|-------------------------------------|--------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| Type 2 Diabetes Mellitus            | E11                                                                            | E10, E12, E13, and E14                                                                                                                                    |
| Heart Failure                       | I11.0, I13.0, I13.2, I25.5, I42.0, I42.5, I42.8, I42.9, I50.0, I50.1, and I50.9                                       |                                                                                                                                                           |
| Hypertension                        | I10, I15                                                                       |                                                                                                                                                           |
| Non-Alcoholic Fatty Liver Disease (NAFLD) | K76.0, K75.8                                                                 | K70, B16, B17, B18, B19, K83.0, K74.3, K75.4, E83.1, E83.0B, E88.0, I82.0, K76.5, K73.9, K73.2, K74.4, K74.5                                               |
| Stroke                              | I60, I61, I63, I64, H34.1, and G45                                             |                                                                                                                                                           |
| Coronary Artery Disease (CAD)       | I20-I25                                                                        |                                                                                                                                                           |
| Myocardial Infarction               | I20, I21                                                                       |                                                                                                                                                           |
