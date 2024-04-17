# Files:
XGboost_three_models.py
# SHAP value
When using SHAP (SHapley Additive exPlanations) to interpret machine learning models, especially for models like XGBoost, normalization is not necessary.
# Covariates selection
In our model, we select 6 covariates that potentially contributes to the Multimorbidity: Age, social economical status, TV hour, education, physical activity, and sleep hour. These covariates are selected based on the previous literatures.
## 1. Age: 
### What's the variable:
Age is a continuous variable, representing the number of years a person has lived.
### Reason why it's important to Multimorbidity: 
Age is a known risk factor for metabolic syndrome. As individuals get older, they are more likely to experience hormonal changes, decreased muscle mass, and other physiological changes that can contribute to the conditions associated with metabolic syndrome.
### Citation: 
[Barnett, K., Mercer, S. W., Norbury, M., Watt, G., Wyke, S., & Guthrie, B. (2012). Epidemiology of multimorbidity and implications for health care, research, and medical education: a cross-sectional study. The Lancet, 380(9836), 37-43.](https://doi.org/10.1016/S0140-6736(12)60240-2)
## 2. Social economical status
### What's the variable: 
Social economic status (SES) encompasses an individual's economic and sociological standing, usually combined factors like income, education, and occupation.
### Reason why it's important to Multimorbidity: 
Individuals with lower socio-economic status tend to have higher rates of multi-morbidity due to various factors like limited access to health care, poorer living conditions, and increased health risk behaviors.
### Citation: 
[Katikireddi, S. V., Skivington, K., Leyland, A. H., Hunt, K., & Mercer, S. W. (2017). The contribution of risk factors to socioeconomic inequalities in multimorbidity across the lifecourse: a longitudinal analysis of the Twenty-07 cohort. BMC medicine, 15(1), 152.](https://doi.org/10.1186/s12916-017-0913-6)
## 3. TV hour
### What's the variable: 
TV hour represents the average number of hours an individual spends watching television daily.

### Reason why it's important to Multimorbidity: 
Extended TV viewing hours are associated with sedentary behavior, which can contribute to the development of multiple chronic diseases and, therefore, multi-morbidity.
### Citation: 
[Stamatakis, E., Hamer, M., & Dunstan, D. W. (2011). Screen-based entertainment time, all-cause mortality, and cardiovascular events: Population-based study with ongoing mortality and hospital events follow-up. Journal of the American College of Cardiology, 57(3), 292-299.](https://doi.org/10.1016/j.jacc.2010.05.065)
## 4. Education
### What's the variable: 
Education is typically categorized by the highest level of schooling or degree an individual has achieved.

### Reason why it's important to Multimorbidity: 
Lower levels of education are linked to higher rates of multi-morbidity, possibly due to reduced health literacy and lower health-promoting behaviors.
### Citation: 
[Schiøtz, M. L., Stockmarr, A., Høst, D., Glümer, C., & Frølich, A. (2017). Social disparities in the prevalence of multimorbidity – A register-based population study. BMC public health, 17(1), 422.](https://doi.org/10.1186/s12889-017-4314-8)
## 5. Physical activity
### What's the variable: 
Regular physical activity can prevent or manage multiple chronic diseases, reducing the risk of multi-morbidity.

### Reason why it's important to Multimorbidity: 
Regular physical activity can help regulate blood sugar, improve cholesterol levels, and maintain a healthy weight, thus decreasing the risk of metabolic syndrome.
### Citation: 
[Pedersen, B. K., & Saltin, B. (2015). Exercise as medicine - evidence for prescribing exercise as therapy in 26 different chronic diseases. Scandinavian journal of medicine & science in sports, 25(S3), 1-72.](https://doi.org/10.1111/sms.12581)
## 6. Sleep hour
### What's the variable: 
Sleep hour refers to the average number of hours an individual sleeps per night.

### Reason why it's important to Multimorbidity: 
Both short and long sleep durations are associated with an increased risk for several chronic conditions, which can collectively contribute to multi-morbidity.
### Citation: 
[Grandner, M. A., Hale, L., Moore, M., & Patel, N. P. (2010). Mortality associated with short sleep duration: The evidence, the possible mechanisms, and the future. Sleep medicine reviews, 14(3), 191-203.](https://doi.org/10.1016/j.smrv.2009.07.006)

## General covariates paper for Multimorbidity
[Barnett, K., Mercer, S. W., Norbury, M., Watt, G., Wyke, S., & Guthrie, B. (2012). Epidemiology of multimorbidity and implications for health care, research, and medical education: a cross-sectional study. The Lancet, 380(9836), 37-43.](https://doi.org/10.1016/S0140-6736(12)60240-2)

[Marmot, M., & Bell, R. (2012). Fair society, healthy lives. Public Health, 126, S4-S10](https://doi.org/10.1016/j.puhe.2012.05.014)

